"""
PDB Filter Worker

Core functionality for filtering PDB structures by component type with interface analysis.
"""

import logging
import os
import tempfile
import numpy as np
from collections import deque, defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

logger = logging.getLogger(__name__)


try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    logger.warning("NetworkX not available. Chain topology analysis will use simplified algorithms.")

try:
    import freesasa
    FREESASA_AVAILABLE = True
except ImportError:
    FREESASA_AVAILABLE = False
    logger.warning("FreeSASA not available. SASA calculations will be disabled.")

from Bio.PDB import PDBIO, PDBParser, Selection
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.NeighborSearch import NeighborSearch

import proprep.structure_prep.chem_comp_dict_fetcher as ccd

        

class WaterAnalyzer:
    """Water analysis for structural assessment."""
    
    def __init__(self, structure: Structure, model_idx: int = 0):
        """Initialize water analyzer."""
        self.structure = structure
        self.model = structure[model_idx]
        
        self.parameters = {
            'metal_distance_cutoff': 2.5,
            'hbond_distance_cutoff': 3.5,
            'hbond_atoms': 'N,O,S',
            'max_hbonds_per_water': 4,
            'interface_distance_cutoff': 5.0,
            'interface_bsa_threshold': 50.0,
            'sasa_probe_radius': 1.4,
            
            # Burial analysis parameters
            'burial_radius': 5.0,                    
            'burial_atom_types': 'protein,hetero,metal',   # metal ions are their own class; a bound ion must shield its water   
            'burial_weighting': 'count',             
        
            # NEW: Network analysis parameters
            'network_type': 'water_only',           # water_only/water_protein_water/all_hbonds
            'network_min_cluster_size': 2,          # Minimum cluster size to report
            'network_show_isolated': True,          # Show single waters
            'network_max_display_size': 8,          # Max network size for detailed ASCII
        } 

        # SASA calculation cache
        self._burial_cache: Dict[str, Any] = {}   # occluder set, KD-tree, enclosure grid (see _burial_context)

    def get_hbond_atoms(self):
        """Get list of atoms to consider for H-bonding."""
        atoms_str = self.parameters.get('hbond_atoms', 'N,O,S')
        return [atom.strip() for atom in atoms_str.split(',')]        

    def set_parameters(self, **kwargs):
        """Update analysis parameters."""
        self.parameters.update(kwargs)
        if 'burial_atom_types' in kwargs or 'sasa_probe_radius' in kwargs:
            self._burial_cache.clear()
        
    def get_metal_atoms(self) -> List[Tuple[Any, str, str]]:
        """Find all metal atoms in the structure."""
        metal_ions = {
            'FE', 'ZN', 'MN', 'CU', 'CO', 'NI', 'MG', 'CA', 
            'NA', 'K', 'LI', 'HG', 'CD', 'PT', 'AU', 'AG'
        }
        
        metals = []
        for chain in self.model:
            for residue in chain:
                if residue.resname.strip() in metal_ions:
                    for atom in residue:
                        metals.append((atom, chain.id, f"{residue.resname}{residue.id[1]}"))
        return metals
        
    def calculate_hydrogen_bonds(self, water_oxygen, chain: Chain) -> Dict[str, Any]:
        """Hydrogen-bond partners of a water from heavy-atom distances alone.

        Candidates are the N/O/S atoms within ``hbond_distance_cutoff`` of the
        oxygen (3.5 Å by default, the Baker-Hubbard donor-acceptor criterion).
        No angular test: hydrogens are usually absent at this point and the
        heavy-atom positions do not fix where they point, so any "angle score"
        would be invented. Candidates are ranked by distance and truncated to
        ``max_hbonds_per_water``, so ``total`` is a count of plausible partners
        within the cutoff, not a hydrogen-bond census.
        """
        hbond_counts = {'protein': 0, 'water': 0, 'hetero': 0, 'total': 0}
        hbond_details = []  # Store details of each H-bond
        
        hbond_atoms = self.get_hbond_atoms()
        protein_residues = {
            'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP',
            'SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN', 'ASP', 'GLU', 
            'LYS', 'ARG', 'HIS', 'MSE'
        }
        
        from Bio.PDB.vectors import Vector
        water_pos = Vector(water_oxygen.coord)
        
        potential_partners = []
        
        # Find potential H-bond partners within distance cutoff
        for search_chain in self.model:
            for residue in search_chain:
                for atom in residue:
                    # Identity, not ==: Biopython's Atom.__eq__ compares by name
                    # (and parent), so == would also skip every backbone O and
                    # every other water's O.
                    if atom is water_oxygen:
                        continue
                        
                    if atom.element not in hbond_atoms:
                        continue
                        
                    distance = (Vector(atom.coord) - water_pos).norm()
                    if distance <= self.parameters['hbond_distance_cutoff']:
                        potential_partners.append({
                            'atom': atom,
                            'residue': residue,
                            'distance': distance,
                            'atom_name': atom.name,
                            'res_name': residue.resname.strip(),
                            'res_num': residue.id[1],
                            'chain_id': residue.parent.id
                        })
        
        # Closest first
        potential_partners.sort(key=lambda x: x['distance'])
        
        # Accept up to max_hbonds_per_water with best scores
        max_hbonds = self.parameters.get('max_hbonds_per_water', 4)
        accepted_partners = potential_partners[:max_hbonds]
        
        # Categorize accepted H-bonds
        for partner in accepted_partners:
            residue = partner['residue']
            
            # Store detailed information
            hbond_details.append({
                'partner_atom': partner['atom_name'],
                'partner_residue': f"{partner['res_name']}{partner['res_num']}",
                'partner_chain': partner['chain_id'],
                'distance': partner['distance'],
            })
            
            if partner['res_name'] in protein_residues:
                hbond_counts['protein'] += 1
            elif partner['res_name'] in {'HOH', 'WAT'}:
                hbond_counts['water'] += 1
            else:
                hbond_counts['hetero'] += 1
                
        hbond_counts['total'] = len(accepted_partners)
        hbond_counts['details'] = hbond_details
        
        return hbond_counts

    def calculate_interfaces(self) -> Dict[str, any]:
        """
        Calculate interface information for the structure.

        Returns a dictionary with chain IDs for quick multi-chain detection.
        This simplified approach identifies chains that could form interfaces.
        """
        interface_data = {
            'chain_ids': set(),
            'has_multiple_chains': False
        }

        # Collect all chain IDs with protein/non-water residues
        for chain in self.model:
            for residue in chain:
                # Skip water molecules
                if residue.resname.strip() not in {'HOH', 'WAT'}:
                    interface_data['chain_ids'].add(chain.id)
                    break  # Found at least one protein residue in this chain

        interface_data['has_multiple_chains'] = len(interface_data['chain_ids']) > 1

        return interface_data

    def is_at_interface(self, water_residue: Residue, interface_data: Dict = None) -> tuple:
        """
        Check if water is at a protein-protein interface.

        A water is considered at an interface if it is simultaneously close to
        atoms from two or more different chains within the cutoff distance.

        Returns:
            tuple: (is_at_interface: bool, bridged_chains: set) where bridged_chains
                   contains the chain IDs that the water is bridging
        """
        if not interface_data or not interface_data.get('has_multiple_chains', False):
            return (False, set())

        from Bio.PDB.vectors import Vector
        water_pos = Vector(water_residue['O'].coord)
        cutoff = self.parameters['interface_distance_cutoff']

        # Track which chains have atoms near this water
        nearby_chains = set()

        for chain in self.model:
            for residue in chain:
                # Skip water molecules
                if residue.resname.strip() in {'HOH', 'WAT'}:
                    continue

                try:
                    # Check if any atom in this residue is close to the water
                    for atom in residue:
                        distance = (Vector(atom.coord) - water_pos).norm()
                        if distance <= cutoff:
                            nearby_chains.add(chain.id)
                            break  # Found close atom, move to next residue
                except:
                    continue

        # Water is at interface only if it's near 2+ different chains
        is_at_interface = len(nearby_chains) >= 2
        return (is_at_interface, nearby_chains if is_at_interface else set())
        
    # Radii for solvent accessibility (Bondi, J. Phys. Chem. 1964), heavy atoms only.
    # Hydrogens are never occluders: the united-atom convention every SASA program uses.
    # Elements not listed (metals) take SASA_DEFAULT_RADIUS.
    SASA_RADII = {
        'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
        'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98, 'SE': 1.90,
    }
    SASA_DEFAULT_RADIUS = 2.00
    # A water is one sphere of radius 1.4 Å whether it is the probe or a
    # crystallographic water (the query oxygen, and other waters when they are
    # occluders): half the O···O distance in liquid water, the Lee-Richards
    # probe convention. Bondi's 1.52 Å is the radius of an oxygen ATOM inside
    # a molecule, a different quantity; using it here gave the same species two
    # sizes and a 2.92 Å "contact" where two real waters touch at 2.80 Å.
    WATER_OXYGEN_RADIUS = 1.40
    # wwPDB validation lists heavy-atom pairs closer than this as close contacts
    # (REMARK 500). A water this close to a C/N/O/S/P atom overlaps it; metals are
    # excluded because 2.0-2.2 Å is ordinary coordination.
    CLASH_DISTANCE = 2.2
    CLASH_ELEMENTS = {'C', 'N', 'O', 'S', 'P'}
    # Grid resolution for the enclosure test. This is a resolution, not a threshold:
    # the bulk/enclosed answer converges as it shrinks. 0.5 Å resolves any channel a
    # 1.4 Å probe can pass through. Grids above ENCLOSURE_MAX_CELLS are coarsened.
    ENCLOSURE_GRID_SPACING = 0.5
    ENCLOSURE_MAX_CELLS = 30_000_000
    # Lee-Richards integrates exactly within each z-slice; the slice count sets how
    # thin an exposed sliver can be and still register. 100 slices over a water's
    # 5.8 Å accessible sphere is 0.06 Å.
    SASA_SLICES = 100

    def _burial_occluders(self):
        """Heavy atoms of the residue classes in ``burial_atom_types``: the atoms that
        can shield a water from solvent. Waters are never occluders unless 'water' is
        listed explicitly: computed among its neighbours, a surface water in a full
        hydration shell would look buried."""
        atom_types = [t.strip().lower() for t in self.parameters['burial_atom_types'].split(',')]
        atoms = []
        for chain in self.model:
            for residue in chain:
                if self._classify_residue_for_burial(residue) not in atom_types:
                    continue
                for atom in residue:
                    element = (atom.element or '').strip().upper()
                    if element == 'H' or element == 'D':
                        continue
                    atoms.append(atom)
        coords = np.array([a.coord for a in atoms], dtype=float).reshape(-1, 3)
        is_water = np.array([self._classify_residue_for_burial(a.get_parent()) == 'water' for a in atoms], dtype=bool)
        radii = np.array([self.SASA_RADII.get((a.element or '').strip().upper(), self.SASA_DEFAULT_RADIUS)
                          for a in atoms], dtype=float)
        radii[is_water] = self.WATER_OXYGEN_RADIUS
        return atoms, coords, radii, is_water

    def _burial_context(self) -> Dict[str, Any]:
        """Build (once per occluder set) everything the burial metrics share."""
        if self._burial_cache:
            return self._burial_cache
        from scipy.spatial import cKDTree
        atoms, coords, radii, is_water = self._burial_occluders()
        probe = float(self.parameters['sasa_probe_radius'])
        ctx = {'atoms': atoms, 'coords': coords, 'radii': radii, 'is_water': is_water, 'probe': probe,
               'tree': cKDTree(coords) if len(coords) else None}
        # Two accessible spheres (r_i + probe, r_j + probe) intersect only if their centres
        # are closer than r_i + r_j + 2 * probe. No atom beyond that distance can remove any
        # of the water's accessible area, so it need not be given to the SASA calculation.
        r_max = float(radii.max()) if len(radii) else self.SASA_DEFAULT_RADIUS
        ctx['sasa_cutoff'] = self.WATER_OXYGEN_RADIUS + r_max + 2.0 * probe
        ctx['enclosure'] = None   # built lazily by _enclosure_grid()
        self._burial_cache = ctx
        return ctx

    def calculate_water_sasa(self, water_oxygen_coord) -> float:
        """Solvent-accessible surface area (Å²) of a water oxygen against the occluders
        alone, Lee-Richards with the configured probe. Only atoms within
        ``sasa_cutoff`` are passed to FreeSASA; the result is identical to using the
        whole structure (see _burial_context)."""
        import math
        ctx = self._burial_context()
        probe, r_w = ctx['probe'], self.WATER_OXYGEN_RADIUS
        isolated = 4.0 * math.pi * (r_w + probe) ** 2
        if ctx['tree'] is None:
            return isolated
        o = np.asarray(water_oxygen_coord, dtype=float)
        idx = ctx['tree'].query_ball_point(o, ctx['sasa_cutoff'])
        # When waters are occluders the query water is among them; a sphere coincident
        # with itself is a degenerate case Lee-Richards gets wrong, so drop it.
        idx = [i for i in idx if np.linalg.norm(ctx['coords'][i] - o) > 1e-3]
        if not idx:
            return isolated
        if not FREESASA_AVAILABLE:
            raise RuntimeError("FreeSASA is required for water burial analysis")
        coords = np.vstack([ctx['coords'][idx], o[None, :]]).ravel().tolist()
        radii = ctx['radii'][idx].tolist() + [r_w]
        result = freesasa.calcCoord(coords, radii, freesasa.Parameters({
            'probe-radius': probe, 'algorithm': freesasa.LeeRichards, 'n-slices': self.SASA_SLICES}))
        return float(result.atomArea(len(radii) - 1))

    def _enclosure_grid(self) -> Dict[str, Any]:
        """Flood fill of probe-centre-accessible space.

        A grid point is accessible if a probe centred there overlaps no occluder, i.e.
        its distance to every atom exceeds r_atom + probe. Connected components of the
        accessible points are labelled; the component touching the box boundary is
        bulk solvent. A water whose surroundings are not part of that component sits
        in a cavity the solvent cannot reach without the protein moving.

        Waters are never walls here, whatever ``burial_atom_types`` says: a water
        cannot be trapped by other, equally mobile, waters.
        """
        from scipy import ndimage
        ctx = self._burial_context()
        if ctx['enclosure'] is not None:
            return ctx['enclosure']
        keep = ~ctx['is_water']
        coords, radii, probe = ctx['coords'][keep], ctx['radii'][keep], ctx['probe']
        if len(coords) == 0:
            ctx['enclosure'] = {'empty': True}
            return ctx['enclosure']
        reach = radii + probe
        pad = float(reach.max()) + 1.0            # boundary layer is guaranteed accessible
        lo = coords.min(axis=0) - pad
        hi = coords.max(axis=0) + pad
        h = self.ENCLOSURE_GRID_SPACING
        n_cells = np.prod(np.ceil((hi - lo) / h) + 1)
        if n_cells > self.ENCLOSURE_MAX_CELLS:
            h = float(h * (n_cells / self.ENCLOSURE_MAX_CELLS) ** (1.0 / 3.0))
            logger.info(f"Enclosure grid coarsened to {h:.2f} Å to stay within {self.ENCLOSURE_MAX_CELLS:,} cells")
        shape = tuple(int(x) for x in np.ceil((hi - lo) / h) + 1)
        axes = [lo[i] + h * np.arange(shape[i]) for i in range(3)]
        blocked = np.zeros(shape, dtype=bool)
        for c, r in zip(coords, reach):
            i0 = np.maximum(np.floor((c - r - lo) / h).astype(int), 0)
            i1 = np.minimum(np.ceil((c + r - lo) / h).astype(int) + 1, shape)
            if np.any(i1 <= i0):
                continue
            gx, gy, gz = np.meshgrid(axes[0][i0[0]:i1[0]] - c[0], axes[1][i0[1]:i1[1]] - c[1],
                                     axes[2][i0[2]:i1[2]] - c[2], indexing='ij')
            blocked[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] |= (gx * gx + gy * gy + gz * gz) < r * r
        accessible = ~blocked
        labels, n_components = ndimage.label(accessible)
        bulk_label = int(labels[0, 0, 0])            # the padded boundary is accessible everywhere
        n_pockets = int(n_components) - (1 if bulk_label else 0)
        ctx['enclosure'] = {'empty': False, 'lo': lo, 'h': h, 'shape': shape,
                            'labels': labels, 'bulk_label': bulk_label, 'n_pockets': n_pockets}
        logger.debug(f"Enclosure grid {shape} at {h:.2f} Å: {n_components} accessible components, "
                     f"{n_pockets} not connected to bulk")
        return ctx['enclosure']

    def classify_water_enclosure(self, water_oxygen_coord) -> str:
        """'bulk' if a bulk-connected probe position lies within touching distance of
        the water (r_water + probe, plus one grid cell), otherwise 'enclosed'."""
        grid = self._enclosure_grid()
        if grid.get('empty'):
            return 'bulk'
        ctx = self._burial_context()
        o = np.asarray(water_oxygen_coord, dtype=float)
        reach = self.WATER_OXYGEN_RADIUS + ctx['probe'] + grid['h']
        lo, h, shape, labels = grid['lo'], grid['h'], grid['shape'], grid['labels']
        i0 = np.maximum(np.floor((o - reach - lo) / h).astype(int), 0)
        i1 = np.minimum(np.ceil((o + reach - lo) / h).astype(int) + 1, shape)
        if np.any(i1 <= i0):
            return 'bulk'          # beyond the padded box: nothing there to enclose it
        sub = labels[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]]
        gx, gy, gz = np.meshgrid(lo[0] + h * np.arange(i0[0], i1[0]) - o[0],
                                 lo[1] + h * np.arange(i0[1], i1[1]) - o[1],
                                 lo[2] + h * np.arange(i0[2], i1[2]) - o[2], indexing='ij')
        within = (gx * gx + gy * gy + gz * gz) <= reach * reach
        return 'bulk' if np.any(sub[within] == grid['bulk_label']) else 'enclosed'

    def analyze_water(self, water_residue: Residue, interface_data: Dict = None) -> Dict[str, Any]:
        """Comprehensive analysis of a single water molecule."""
        try:
            water_oxygen = water_residue['O']
        except KeyError:
            logger.warning(f"Water residue {water_residue.id} missing oxygen atom")
            return {}
            
        results = {
            'residue_number': water_residue.id[1],
            'residue_name': water_residue.resname,
            'chain_id': water_residue.parent.id,
            'b_factor': water_oxygen.bfactor,
            'protein_median_b': self.protein_median_bfactor(),
        }
        if results['protein_median_b']:
            results['b_factor_ratio'] = water_oxygen.bfactor / results['protein_median_b']
        
        # Metal distance analysis
        metals = self.get_metal_atoms()
        min_metal_distance = float('inf')
        closest_metal = None
        
        from Bio.PDB.vectors import Vector
        water_pos = Vector(water_oxygen.coord)
        for metal_atom, chain_id, residue_info in metals:
            distance = (Vector(metal_atom.coord) - water_pos).norm()
            if distance < min_metal_distance:
                min_metal_distance = distance
                closest_metal = f"{residue_info} (Chain {chain_id})"
                
        results['metal_distance'] = min_metal_distance if min_metal_distance != float('inf') else None
        results['closest_metal'] = closest_metal
        results['coordinating_metal'] = min_metal_distance <= self.parameters['metal_distance_cutoff'] if min_metal_distance != float('inf') else False
        
        # Hydrogen bond analysis
        hbond_counts = self.calculate_hydrogen_bonds(water_oxygen, water_residue.parent)
        results.update(hbond_counts)
        
        # Burial: solvent accessibility of the oxygen and bulk connectivity
        burial_data = self.calculate_burial_analysis(water_residue)
        results.update(burial_data)

        # Interface analysis
        at_interface, bridged_chains = self.is_at_interface(water_residue, interface_data)
        results['at_interface'] = at_interface
        results['bridged_chains'] = sorted(list(bridged_chains))  # Store as sorted list for display

        return results
        
    def protein_median_bfactor(self) -> Optional[float]:
        """Median B-factor of the protein heavy atoms, the structure's own yardstick
        for "ordered". B-factors scale with resolution and refinement, so an
        absolute cutoff means different things in different entries; a water
        below the protein median is at least as ordered as a typical protein atom
        in *this* structure. None when the model has no protein atoms."""
        if hasattr(self, '_protein_median_b'):
            return self._protein_median_b
        bs = [a.bfactor for chain in self.model for res in chain if res.id[0] == ' '
              for a in res if (a.element or '').strip().upper() not in ('H', 'D')]
        self._protein_median_b = float(np.median(bs)) if bs else None
        return self._protein_median_b

    def _classify_residue_for_burial(self, residue) -> str:
        """Classify residue type for burial analysis."""
        resname = residue.resname.strip()

        # Standard amino acids
        standard_aa = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }

        # Water molecules
        if resname in ['HOH', 'WAT', 'TIP', 'SOL']:
            return 'water'

        # Metal ions - comprehensive periodic table list
        # (matches the METALS set from comprehensive_redox_detector)
        metal_ions = {
            # Group 1: Alkali metals
            'LI', 'NA', 'K', 'RB', 'CS', 'FR',
            # Group 2: Alkaline earth metals
            'BE', 'MG', 'CA', 'SR', 'BA', 'RA',
            # Group 3-12: Transition metals
            'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
            'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
            'LA', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
            'AC', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN',
            # Lanthanides
            'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU',
            # Actinides
            'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR',
            # Post-transition metals
            'AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'PO',
            # Metalloids that often behave as metals in coordination
            'SB', 'TE'
        }
        if resname in metal_ions:
            return 'metal'

        # Standard amino acids
        if resname in standard_aa:
            return 'protein'

        # Everything else is hetero
        return 'hetero'

    # Van der Waals radii from Charry & Tkatchenko (2024), Table S1
    # Reference: Charry, J. & Tkatchenko, A. "van der Waals Radii of Free and Bonded
    #            Atoms From Hydrogen (Z=1) to Oganesson (Z=118)"
    #            J. Chem. Theory Comput. 2024, DOI: 10.1021/acs.jctc.4c00784
    # Values in Angstroms (converted from bohr: 1 bohr = 0.529177 Å)
    # Complete table for all elements Z=1 to Z=118
    VDW_RADII = {
        # Row 1
        'H': 1.674, 'HE': 1.415,
        # Row 2
        'LI': 2.799, 'BE': 2.270, 'B': 2.079, 'C': 1.910, 'N': 1.798, 'O': 1.714, 'F': 1.631, 'NE': 1.554,
        # Row 3
        'NA': 2.796, 'MG': 2.484, 'AL': 2.412, 'SI': 2.267, 'P': 2.140, 'S': 2.063, 'CL': 1.981, 'AR': 1.906,
        # Row 4
        'K': 3.038, 'CA': 2.792, 'SC': 2.597, 'TI': 2.608, 'V': 2.557, 'CR': 2.540, 'MN': 2.471, 'FE': 2.440,
        'CO': 2.394, 'NI': 2.356, 'CU': 2.343, 'ZN': 2.278, 'GA': 2.363, 'GE': 2.288, 'AS': 2.196, 'SE': 2.185,
        'BR': 2.087, 'KR': 2.022,
        # Row 5
        'RB': 3.078, 'SR': 2.872, 'Y': 2.794, 'ZR': 2.651, 'NB': 2.599, 'MO': 2.557, 'TC': 2.522, 'RU': 2.490,
        'RH': 2.456, 'PD': 2.153, 'AG': 2.395, 'CD': 2.335, 'IN': 2.453, 'SN': 2.382, 'SB': 2.313, 'TE': 2.272,
        'I': 2.227, 'XE': 2.167,
        # Row 6
        'CS': 3.179, 'BA': 3.009, 'LA': 2.909, 'CE': 2.890, 'PR': 2.912, 'ND': 2.896, 'PM': 2.881, 'SM': 2.864,
        'EU': 2.845, 'GD': 2.785, 'TB': 2.813, 'DY': 2.802, 'HO': 2.779, 'ER': 2.764, 'TM': 2.748, 'YB': 2.734,
        'LU': 2.728, 'HF': 2.619, 'TA': 2.499, 'W': 2.466, 'RE': 2.437, 'OS': 2.407, 'IR': 2.389, 'PT': 2.349,
        'AU': 2.254, 'HG': 2.235, 'TL': 2.363, 'PB': 2.342, 'BI': 2.349, 'PO': 2.320, 'AT': 2.304, 'RN': 2.245,
        # Row 7
        'FR': 3.076, 'RA': 2.967, 'AC': 2.886, 'TH': 2.916, 'PA': 2.775, 'U': 2.705, 'NP': 2.768, 'PU': 2.716,
        'AM': 2.710, 'CM': 2.748, 'BK': 2.693, 'CF': 2.683, 'ES': 2.673, 'FM': 2.657, 'MD': 2.641, 'NO': 2.644,
        'LR': 3.080, 'RF': 2.651, 'DB': 2.304, 'SG': 2.288, 'BH': 2.272, 'HS': 2.254, 'MT': 2.236, 'DS': 2.216,
        'RG': 2.218, 'CN': 2.175, 'NH': 2.185, 'FL': 2.207, 'MC': 2.482, 'TS': 2.508, 'OG': 2.414,
    }

    def _calculate_burial_weight(self, distance: float, weighting_scheme: str, atom=None) -> float:
        """
        Calculate burial weight based on distance and scheme.

        Args:
            distance: Distance to atom in Angstroms
            weighting_scheme: 'count', 'distance', or 'vdw'
            atom: BioPython atom object (needed for vdw scheme)

        Returns:
            Weight contribution from this atom
        """
        if weighting_scheme == 'count':
            return 1.0
        elif weighting_scheme == 'distance':
            return 1.0 / distance if distance > 0 else 1.0
        elif weighting_scheme == 'vdw':
            # vdW-scaled weighting: weight based on vdW radius overlap
            if atom is None:
                return 1.0 / distance if distance > 0 else 1.0

            element = atom.element.strip().upper()
            vdw_radius = self.VDW_RADII.get(element, 1.7)  # Default to ~carbon radius

            # Weight by (sum of vdW radii) / distance
            # This gives higher weight when atoms are within vdW contact
            water_vdw = self.VDW_RADII['O']  # Water oxygen radius
            combined_vdw = water_vdw + vdw_radius

            if distance > 0:
                return combined_vdw / distance
            else:
                return 1.0
        else:
            return 1.0

    def calculate_burial_analysis(self, water_residue: Residue) -> Dict[str, Any]:
        """Burial of one water from geometry alone.

        * ``burial_sasa``: accessible area (Å²) of the oxygen against the occluders,
          1.4 Å probe. Zero means no water-sized probe can touch it.
        * ``burial_covered_pct``: 100 × (1 − sasa / area of an isolated water oxygen),
          so the same number reads as "how much of this water is covered".
        * ``burial_access``: 'bulk' or 'enclosed' from the flood fill.
        * ``burial_closest_distance`` / ``burial_closest_atom``: nearest occluder.
        * ``burial_category``: 'Clash' (overlaps a protein atom by the wwPDB
          close-contact criterion), 'Enclosed' (no path to bulk solvent), 'Buried'
          (touchable by nothing, but bulk-connected: bottom of a cleft), 'Exposed'.
        No calibration constants: the only inputs are the probe radius and the vdW radii.
        """
        try:
            water_oxygen = water_residue['O']
        except KeyError:
            logger.warning(f"Water residue {water_residue.id} missing oxygen atom")
            return {}

        import math
        ctx = self._burial_context()
        o = water_oxygen.coord.astype(float)
        sasa = self.calculate_water_sasa(o)
        access = self.classify_water_enclosure(o)

        closest_distance, closest_atom, clash = None, None, False
        if ctx['tree'] is not None:
            d, i = ctx['tree'].query(o)
            atom = ctx['atoms'][int(i)]
            res = atom.get_parent()
            closest_distance = float(d)
            closest_atom = f"{res.resname}{res.id[1]} {atom.get_id()}"
            element = (atom.element or '').strip().upper()
            clash = element in self.CLASH_ELEMENTS and closest_distance < self.CLASH_DISTANCE

        if clash:
            category = 'Clash'
        elif access == 'enclosed':
            category = 'Enclosed'
        elif sasa == 0.0:
            category = 'Buried'
        else:
            category = 'Exposed'

        isolated = 4.0 * math.pi * (self.WATER_OXYGEN_RADIUS + ctx['probe']) ** 2
        return {
            'burial_sasa': sasa,
            'burial_sasa_isolated': isolated,
            'burial_covered_pct': min(100.0, max(0.0, 100.0 * (1.0 - sasa / isolated))),
            'burial_access': access,
            'burial_closest_distance': closest_distance,
            'burial_closest_atom': closest_atom,
            'burial_category': category,
        }

    def calculate_burial_profile(self, water_residue: Residue, 
                            min_radius: float = 2.0, 
                            max_radius: float = 8.0, 
                            step: float = 0.5) -> Dict[str, Any]:
        """Calculate burial profile across multiple radii."""
        try:
            water_oxygen = water_residue['O']
        except KeyError:
            logger.warning(f"Water residue {water_residue.id} missing oxygen atom")
            return {}
        
        from Bio.PDB.vectors import Vector
        import numpy as np
        
        # Get atom types and weighting from current parameters
        atom_types_str = self.parameters['burial_atom_types']
        weighting = self.parameters['burial_weighting']
        atom_types = [t.strip().lower() for t in atom_types_str.split(',')]
        
        # Build list of all relevant atoms
        all_atoms = []
        for chain in self.model:
            for residue in chain:
                if residue is water_residue:
                    continue
                res_type = self._classify_residue_for_burial(residue)
                if res_type in atom_types:
                    for atom in residue:
                        all_atoms.append(atom)
        
        if not all_atoms:
            return {}
        
        # Calculate burial at each radius
        radii = np.arange(min_radius, max_radius + step, step)
        burial_counts = []
        water_pos = Vector(water_oxygen.coord)
        
        # Pre-calculate all atom distances for efficiency
        atom_distances = []
        for atom in all_atoms:
            if atom is water_oxygen:
                continue
            distance = (Vector(atom.coord) - water_pos).norm()
            atom_distances.append((atom, distance))
        
        # Sort by distance for efficient radius searching
        atom_distances.sort(key=lambda x: x[1])
        
        for radius in radii:
            total_weight = 0.0
            for atom, distance in atom_distances:
                if distance > radius:
                    break  # All remaining atoms are farther
                weight = self._calculate_burial_weight(distance, weighting)
                total_weight += weight
            burial_counts.append(total_weight)
        
        # Analyze the profile
        burial_counts = np.array(burial_counts)
        
        # Find saturation point (where burial levels off)
        saturation_radius = None
        saturation_count = None
        if len(burial_counts) > 3:
            for i in range(2, len(burial_counts)):
                if burial_counts[i-1] > 0:
                    increase_rate = (burial_counts[i] - burial_counts[i-1]) / burial_counts[i-1]
                    if increase_rate < 0.1:  # Less than 10% increase
                        saturation_radius = radii[i]
                        saturation_count = burial_counts[i]
                        break
        
        # Find steepest rise region
        steepest_start = None
        steepest_end = None
        max_slope = 0
        if len(burial_counts) > 1:
            slopes = np.diff(burial_counts) / step
            max_slope_idx = np.argmax(slopes)
            max_slope = slopes[max_slope_idx]
            steepest_start = radii[max_slope_idx]
            steepest_end = radii[max_slope_idx + 1]
        
        return {
            'radii': radii.tolist(),
            'burial_counts': burial_counts.tolist(),
            'saturation_radius': saturation_radius,
            'saturation_count': saturation_count,
            'steepest_start': steepest_start,
            'steepest_end': steepest_end,
            'max_slope': max_slope,
            'final_count': burial_counts[-1] if len(burial_counts) > 0 else 0
        }

    def generate_burial_profile_ascii(self, profile_data: Dict[str, Any], 
                                    width: int = 20, height: int = 12) -> str:
        """Generate ASCII art visualization of burial profile."""
        if not profile_data or 'radii' not in profile_data:
            return "No profile data available"
        
        radii = profile_data['radii']
        counts = profile_data['burial_counts']
        
        if not radii or not counts:
            return "No data to plot"
        
        import numpy as np
        
        # Normalize counts to fit in height
        max_count = max(counts) if counts else 1
        if max_count == 0:
            max_count = 1
        
        # Create the grid
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        # Fill the grid with data
        for i, count in enumerate(counts):
            if i >= width:
                break
            
            # Calculate bar height
            bar_height = int((count / max_count) * (height - 1))
            
            # Fill column from bottom up
            for j in range(bar_height):
                row_idx = height - 1 - j
                if 0 <= row_idx < height:
                    grid[row_idx][i] = '█'
        
        # Convert grid to string
        chart_lines = []
        
        # Add count scale on left
        for i, row in enumerate(grid):
            row_count = int(max_count * (height - 1 - i) / (height - 1))
            scale_str = f"{row_count:3d} |"
            chart_lines.append(scale_str + ''.join(row))
        
        # Create x-axis with tick marks
        tick_line = "    "
        label_line = "    "
        
        # Determine which positions get tick marks and labels
        tick_positions = []
        if len(radii) <= width:
            # Show ticks at regular intervals
            interval = max(2, len(radii) // 6)  # Every 2-4 positions
            for i in range(0, len(radii), interval):
                if i < width:
                    tick_positions.append(i)
            # Always include the last position
            if len(radii) - 1 < width and (len(radii) - 1) not in tick_positions:
                tick_positions.append(len(radii) - 1)
        
        # Build tick and label lines
        for i in range(width):
            if i in tick_positions and i < len(radii):
                tick_line += "│"
                # Format radius value
                radius_val = radii[i]
                if radius_val == int(radius_val):
                    label = f"{int(radius_val)}"
                else:
                    label = f"{radius_val:.1f}"
                
                # Center the label under the tick
                label_line += label[0] if len(label) > 0 else " "
            else:
                tick_line += "─"
                label_line += " "
        
        chart_lines.append(tick_line)
        chart_lines.append(label_line)
        chart_lines.append("    Radius (Å)")
        
        # Convert to string
        chart_str = "\n".join(chart_lines)
        
        # Add analysis annotations
        annotations = []
        if profile_data.get('saturation_radius'):
            annotations.append(f"Saturation: ~{profile_data['saturation_radius']:.1f}Å")
        if profile_data.get('steepest_start') and profile_data.get('steepest_end'):
            annotations.append(f"Steep rise: {profile_data['steepest_start']:.1f}-{profile_data['steepest_end']:.1f}Å")
        
        if annotations:
            chart_str += "\n\nKey features: " + ", ".join(annotations)
        
        return chart_str

    def calculate_directional_burial(self, water_residue: Residue, 
                                radius: float = None) -> Dict[str, Any]:
        """Calculate burial in 8 directional sectors around a water molecule."""
        try:
            water_oxygen = water_residue['O']
        except KeyError:
            logger.warning(f"Water residue {water_residue.id} missing oxygen atom")
            return {}
        
        from Bio.PDB.vectors import Vector
        from Bio.PDB.NeighborSearch import NeighborSearch
        import numpy as np
        
        # Use provided radius or parameter default
        if radius is None:
            radius = self.parameters['burial_radius']
        
        # Get atom types and weighting from current parameters
        atom_types_str = self.parameters['burial_atom_types']
        weighting = self.parameters['burial_weighting']
        atom_types = [t.strip().lower() for t in atom_types_str.split(',')]
        
        # Build neighbor search for all relevant atoms
        all_atoms = []
        for chain in self.model:
            for residue in chain:
                if residue is water_residue:
                    continue
                res_type = self._classify_residue_for_burial(residue)
                if res_type in atom_types:
                    for atom in residue:
                        all_atoms.append(atom)
        
        if not all_atoms:
            return {}
        
        # Find atoms within radius
        ns = NeighborSearch(all_atoms)
        water_pos = Vector(water_oxygen.coord)
        nearby_atoms = ns.search(water_oxygen.coord, radius)
        
        # Define 8 directional sectors (45° each)
        sectors = {
            'N': {'min': 337.5, 'max': 22.5, 'count': 0, 'weight': 0.0},
            'NE': {'min': 22.5, 'max': 67.5, 'count': 0, 'weight': 0.0},
            'E': {'min': 67.5, 'max': 112.5, 'count': 0, 'weight': 0.0},
            'SE': {'min': 112.5, 'max': 157.5, 'count': 0, 'weight': 0.0},
            'S': {'min': 157.5, 'max': 202.5, 'count': 0, 'weight': 0.0},
            'SW': {'min': 202.5, 'max': 247.5, 'count': 0, 'weight': 0.0},
            'W': {'min': 247.5, 'max': 292.5, 'count': 0, 'weight': 0.0},
            'NW': {'min': 292.5, 'max': 337.5, 'count': 0, 'weight': 0.0}
        }
        
        # Analyze each nearby atom
        for atom in nearby_atoms:
            if atom is water_oxygen:
                continue
            
            # Calculate direction vector
            atom_pos = Vector(atom.coord)
            direction = atom_pos - water_pos
            distance = direction.norm()
            
            if distance == 0:
                continue
            
            # Calculate azimuthal angle (0-360°)
            angle_rad = np.arctan2(direction[1], direction[0])  # y, x for proper orientation
            angle_deg = np.degrees(angle_rad)
            if angle_deg < 0:
                angle_deg += 360
            
            # Determine sector
            sector = None
            for sector_name, sector_data in sectors.items():
                min_angle = sector_data['min']
                max_angle = sector_data['max']
                
                # Handle wrap-around for North sector
                if min_angle > max_angle:  # North sector case
                    if angle_deg >= min_angle or angle_deg <= max_angle:
                        sector = sector_name
                        break
                else:
                    if min_angle <= angle_deg < max_angle:
                        sector = sector_name
                        break
            
            if sector:
                sectors[sector]['count'] += 1
                weight = self._calculate_burial_weight(distance, weighting)
                sectors[sector]['weight'] += weight
        
        # Identify primary and secondary burial directions
        sorted_sectors = sorted(sectors.items(), key=lambda x: x[1]['weight'], reverse=True)
        primary_direction = sorted_sectors[0][0] if sorted_sectors[0][1]['weight'] > 0 else None
        secondary_direction = sorted_sectors[1][0] if len(sorted_sectors) > 1 and sorted_sectors[1][1]['weight'] > 0 else None

        # Find the ACTUAL least buried direction (not just opposite of primary)
        least_buried_sectors = sorted(sectors.items(), key=lambda x: x[1]['weight'])
        pocket_opening = least_buried_sectors[0][0] if least_buried_sectors[0][1]['weight'] >= 0 else None

        # Calculate pattern metrics for better classification
        sector_counts = [s['count'] for s in sectors.values()]
        non_zero_sectors = sum(1 for count in sector_counts if count > 0)
        min_count = min(sector_counts) if sector_counts else 0
        max_count = max(sector_counts) if sector_counts else 0
        count_range = max_count - min_count
        avg_count = sum(sector_counts) / len(sector_counts) if sector_counts else 0

        # Determine burial pattern type
        # The labels are conventions, not findings: the sector-count range (max − min)
        # is compared with the mean sector count at 0.5× and 1.5×. Each label states
        # its rule so a reader can weigh it against the counts themselves.
        if count_range == 0:
            pattern_type = "No burial (all sectors empty)"
        elif count_range <= avg_count * 0.5:
            pattern_type = "Uniform (sector range ≤ 0.5 × mean)"
        elif count_range <= avg_count * 1.5:
            pattern_type = "Moderately directional (range ≤ 1.5 × mean)"
        else:
            pattern_type = "Highly directional (range > 1.5 × mean)"

        return {
            'sectors': sectors,
            'primary_direction': primary_direction,
            'secondary_direction': secondary_direction,
            'pocket_opening': pocket_opening,
            'total_weight': sum(s['weight'] for s in sectors.values()),
            'max_sector_weight': max(s['weight'] for s in sectors.values()) if sectors else 0,
            'min_sector_count': min_count,
            'max_sector_count': max_count,
            'count_range': count_range,
            'pattern_type': pattern_type,
            'analysis_radius': radius
        }

    def generate_directional_compass_ascii(self, directional_data: Dict[str, Any]) -> str:
        """Generate ASCII compass visualization of directional burial."""
        # Add error checking at the start
        if not directional_data:
            return "No directional data available - analysis may have failed"
        
        if 'sectors' not in directional_data:
            return f"Invalid directional data - missing sectors. Available keys: {list(directional_data.keys())}"
        
        sectors = directional_data['sectors']
        
        if not sectors:
            return "No sector data available - no atoms found within analysis radius"
        
        # Verify all expected sectors exist
        expected_sectors = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
        missing_sectors = [s for s in expected_sectors if s not in sectors]
        if missing_sectors:
            return f"Missing sector data: {missing_sectors}. Available sectors: {list(sectors.keys())}"
        
        # Get max count for scaling (use count, not weight for symbol scaling)
        try:
            max_count = max(s['count'] for s in sectors.values()) if sectors else 1
            if max_count == 0:
                max_count = 1
        except (KeyError, TypeError) as e:
            return f"Error accessing sector count data: {e}. Sector structure: {sectors}"
        
        # Create visual representation based on ATOM COUNT not weight
        def get_symbol(count):
            try:
                if count == 0:
                    return '○'
                elif count <= max_count * 0.25:
                    return '●'
                elif count <= max_count * 0.5:
                    return '●●'
                elif count <= max_count * 0.75:
                    return '●●●'
                else:
                    return '████'
            except:
                return '?'
        
        # Build compass using atom counts with error handling
        try:
            n_sym = get_symbol(sectors['N']['count'])
            ne_sym = get_symbol(sectors['NE']['count'])
            e_sym = get_symbol(sectors['E']['count'])
            se_sym = get_symbol(sectors['SE']['count'])
            s_sym = get_symbol(sectors['S']['count'])
            sw_sym = get_symbol(sectors['SW']['count'])
            w_sym = get_symbol(sectors['W']['count'])
            nw_sym = get_symbol(sectors['NW']['count'])
        except KeyError as e:
            return f"Missing sector count data: {e}. Check sector structure: {sectors}"
        
        # Format compass with proper spacing and consistent width
        compass_lines = []
        
        try:
            # Build compass with centered alignment
            compass_lines.append(f"      {n_sym}({sectors['N']['count']})      ")
            compass_lines.append(f"   {nw_sym}({sectors['NW']['count']})   ⚬   {ne_sym}({sectors['NE']['count']})   ")
            compass_lines.append(f"{w_sym}({sectors['W']['count']})  ●●●  {e_sym}({sectors['E']['count']})")
            compass_lines.append(f"   {sw_sym}({sectors['SW']['count']})  ●●●  {se_sym}({sectors['SE']['count']})   ")
            compass_lines.append(f"      {s_sym}({sectors['S']['count']})      ")
            
            compass_str = "\n".join(compass_lines)
            
            # Add legend and analysis
            compass_str += "\n\n⚬ = Water molecule, ○ = No atoms, ● = Protein atoms"
            compass_str += "\nNumbers in () show actual atom count per sector"
            compass_str += "\nSymbol intensity shows relative burial density"
            
            # Add directional analysis with error handling
            if directional_data.get('primary_direction'):
                primary = directional_data['primary_direction']
                primary_count = sectors[primary]['count']
                compass_str += f"\nPrimary burial: {primary} ({primary_count} atoms - highest density)"
                
            if directional_data.get('secondary_direction'):
                secondary = directional_data['secondary_direction']
                secondary_count = sectors[secondary]['count']
                compass_str += f"\nSecondary burial: {secondary} ({secondary_count} atoms)"

            if directional_data.get('pocket_opening'):
                opening = directional_data['pocket_opening']
                opening_count = sectors[opening]['count']
                compass_str += f"\nLeast buried direction: {opening} ({opening_count} atoms - potential access route)"

            # Add interpretation with error handling
            total_atoms = sum(s['count'] for s in sectors.values())
            compass_str += f"\nTotal surrounding atoms: {total_atoms}"

            # Use the improved pattern classification
            pattern_type = directional_data.get('pattern_type', 'Unknown pattern')
            compass_str += f"\nBurial pattern: {pattern_type}"

            # Add range information for context
            min_count = directional_data.get('min_sector_count', 0)
            max_count = directional_data.get('max_sector_count', 0)
            count_range = directional_data.get('count_range', 0)
            compass_str += f"\nCount range: {min_count}-{max_count} atoms (variation: {count_range})"

            # Add biological interpretation
            if count_range == 0:
                compass_str += "\nInterpretation: No surrounding atoms detected"
            elif count_range <= 2:
                compass_str += "\nInterpretation: Relatively uniform environment (bulk solvent)"
            elif count_range <= 4:
                compass_str += "\nInterpretation: Moderate asymmetry (surface water or partial cavity)"
            else:
                compass_str += "\nInterpretation: High asymmetry (binding site, crevice, or active site)"

            # Noe about analysis interpretation 
            compass_str += "\n\nNote: Compass directions depend on PDB coordinate frame."
            compass_str += "\nFocus on burial patterns and asymmetry, not absolute directions."

            return compass_str
            
        except Exception as e:
            return f"Error generating compass visualization: {e}\nDirectional data: {directional_data}"
        
    def analyze_water_comprehensive(self, water_residue: Residue) -> Dict[str, Any]:
        """Most comprehensive water analysis including all burial methods."""
        # Start with basic analysis (includes Phase 1 burial)
        analysis = self.analyze_water(water_residue)
        
        # Add burial profile (Phase 2)
        profile_data = self.calculate_burial_profile(water_residue)
        analysis['burial_profile'] = profile_data
        
        if profile_data:
            ascii_chart = self.generate_burial_profile_ascii(profile_data)
            analysis['burial_profile_ascii'] = ascii_chart
        
        # Add directional analysis (Phase 3)
        directional_data = self.calculate_directional_burial(water_residue)
        analysis['directional_burial'] = directional_data
                    
        if directional_data:
            compass_ascii = self.generate_directional_compass_ascii(directional_data)
            analysis['directional_compass'] = compass_ascii
        
        return analysis

    def build_water_network(self, water_residues: List[Residue]):
        """Build water-water H-bond network."""
        # Check if NetworkX is available
        if not NETWORKX_AVAILABLE:
            logger.warning("NetworkX not available - cannot perform network analysis")
            return None
        
        import networkx as nx
        
        network = nx.Graph()
        
        # Add all waters as nodes with their data
        for water in water_residues:
            network.add_node(water.id[1], 
                            residue_name=water.resname,
                            chain_id=water.parent.id,
                            residue_obj=water)
        
        # Build water-water connections based on H-bond analysis
        water_connections = {}
        
        for water in water_residues:
            try:
                water_oxygen = water['O']
                hbond_data = self.calculate_hydrogen_bonds(water_oxygen, water.parent)
                
                # Extract water partners from H-bond analysis
                water_partners = []
                if 'water_partners' in hbond_data:
                    water_partners = hbond_data['water_partners']
                elif 'water' in hbond_data and hbond_data['water'] > 0:
                    # If water partners not explicitly stored, find them
                    water_partners = self._find_water_hbond_partners(water_oxygen, water_residues)
                
                water_connections[water.id[1]] = water_partners
                
            except Exception as e:
                logger.debug(f"Error processing water {water.id} for network: {e}")
                water_connections[water.id[1]] = []
        
        # Add edges for water-water H-bonds
        for water_id, partners in water_connections.items():
            for partner_id in partners:
                if partner_id in [w.id[1] for w in water_residues]:
                    # Only add if both waters exist and avoid self-loops
                    if water_id != partner_id and network.has_node(partner_id):
                        network.add_edge(water_id, partner_id)
        
        return network

    def _find_water_hbond_partners(self, water_oxygen, water_residues: List[Residue]) -> List[int]:
        """Find water molecules H-bonded to this water."""
        from Bio.PDB.vectors import Vector
        from Bio.PDB.NeighborSearch import NeighborSearch
        
        partners = []
        hbond_cutoff = self.parameters['hbond_distance_cutoff']
        
        # Get all other water oxygens
        other_waters = []
        for water in water_residues:
            try:
                other_oxygen = water['O']
                if other_oxygen != water_oxygen:
                    other_waters.append((other_oxygen, water.id[1]))
            except KeyError:
                continue
        
        if not other_waters:
            return partners
        
        # Find nearby water oxygens
        water_pos = Vector(water_oxygen.coord)
        for other_oxygen, other_id in other_waters:
            distance = (Vector(other_oxygen.coord) - water_pos).norm()
            if distance <= hbond_cutoff:
                partners.append(other_id)
        
        return partners

    def analyze_water_networks(self, water_residues: List[Residue]) -> Dict[str, Any]:
        """Analyze water network topology and connectivity."""
        if not NETWORKX_AVAILABLE:
            return {
                'error': 'NetworkX not available',
                'total_waters': len(water_residues),
                'network_type': 'unavailable'
            }
        
        import networkx as nx
        
        # Build the network
        network = self.build_water_network(water_residues)
        if network is None:
            return {'error': 'Network construction failed'}
        
        # Basic network statistics
        components = list(nx.connected_components(network))
        
        # Filter components by minimum size
        min_size = self.parameters.get('network_min_cluster_size', 2)
        significant_components = [comp for comp in components if len(comp) >= min_size]
        isolated_waters = [comp for comp in components if len(comp) == 1]
        
        # Find hub waters (highly connected)
        hub_waters = []
        for node in network.nodes():
            degree = network.degree(node)
            if degree >= 3:  # 3 or more connections
                hub_waters.append((node, degree))
        
        # Sort hubs by degree
        hub_waters.sort(key=lambda x: x[1], reverse=True)
        
        # Find longest chains
        longest_chains = []
        for component in significant_components:
            subgraph = network.subgraph(component)
            if nx.is_tree(subgraph):  # Chain-like structure
                # Find diameter (longest path)
                try:
                    diameter = nx.diameter(subgraph)
                    # Find the actual path
                    for node1 in component:
                        for node2 in component:
                            if node1 != node2:
                                try:
                                    path = nx.shortest_path(subgraph, node1, node2)
                                    if len(path) - 1 == diameter:  # -1 because path length = nodes - 1
                                        longest_chains.append(path)
                                        break
                                except:
                                    continue
                        if longest_chains:
                            break
                except:
                    pass
        
        return {
            'network': network,
            'total_waters': network.number_of_nodes(),
            'total_connections': network.number_of_edges(),
            'total_components': len(components),
            'significant_components': significant_components,
            'isolated_waters': [list(comp)[0] for comp in isolated_waters],
            'hub_waters': hub_waters,
            'longest_chains': longest_chains,
            'largest_cluster_size': len(max(components, key=len)) if components else 0,
            'network_type': self.parameters.get('network_type', 'water_only'),
            'show_isolated': self.parameters.get('network_show_isolated', True)
        }

    def generate_network_ascii(self, network_analysis: Dict[str, Any]) -> str:
        """Generate ASCII representation of water networks."""
        if 'error' in network_analysis:
            return f"Network Analysis Error: {network_analysis['error']}"
        
        if not NETWORKX_AVAILABLE:
            return "NetworkX not available for network visualization"
        
        import networkx as nx
        
        network = network_analysis.get('network')
        if network is None:
            return "No network data available"
        
        significant_components = network_analysis.get('significant_components', [])
        isolated_waters = network_analysis.get('isolated_waters', [])
        hub_waters = network_analysis.get('hub_waters', [])
        longest_chains = network_analysis.get('longest_chains', [])
        
        # Start building ASCII output
        ascii_lines = []
        ascii_lines.append(f"Water Network Analysis ({network_analysis['network_type']})")
        ascii_lines.append("=" * 60)
        ascii_lines.append("")
        
        # Summary statistics
        ascii_lines.append("Network Summary:")
        ascii_lines.append(f"  Total waters: {network_analysis['total_waters']}")
        ascii_lines.append(f"  Total H-bond connections: {network_analysis['total_connections']}")
        ascii_lines.append(f"  Connected clusters: {len(significant_components)}")
        ascii_lines.append(f"  Isolated waters: {len(isolated_waters)}")
        ascii_lines.append(f"  Largest cluster: {network_analysis['largest_cluster_size']} waters")
        ascii_lines.append("")
        
        # Hub waters
        if hub_waters:
            ascii_lines.append("Hub Waters (≥3 connections):")
            for water_id, degree in hub_waters[:5]:  # Show top 5
                ascii_lines.append(f"  HOH{water_id}: {degree} connections")
            if len(hub_waters) > 5:
                ascii_lines.append(f"  ... and {len(hub_waters)-5} more hubs")
            ascii_lines.append("")
        
        # Individual network components
        max_display_size = self.parameters.get('network_max_display_size', 8)
        
        for i, component in enumerate(significant_components[:10]):  # Show first 10 components
            size = len(component)
            ascii_lines.append(f"Cluster {i+1} ({size} waters):")
            
            if size <= max_display_size:
                # Small networks: show detailed connectivity
                ascii_lines.extend(self._draw_small_network(network, component))
            else:
                # Large networks: show summary
                ascii_lines.extend(self._draw_large_network_summary(network, component))
            
            ascii_lines.append("")
        
        # Show longest chains
        if longest_chains:
            ascii_lines.append("Longest Water Chains:")
            for i, chain in enumerate(longest_chains[:3]):  # Show top 3
                chain_str = " → ".join([f"HOH{water_id}" for water_id in chain])
                ascii_lines.append(f"  Chain {i+1}: {chain_str} ({len(chain)} waters)")
            ascii_lines.append("")
        
        # Isolated waters
        if isolated_waters and network_analysis.get('show_isolated', True):
            ascii_lines.append("Isolated Waters (no H-bond connections):")
            isolated_str = ", ".join([f"HOH{water_id}" for water_id in isolated_waters[:20]])
            if len(isolated_waters) > 20:
                isolated_str += f" ... and {len(isolated_waters)-20} more"
            ascii_lines.append(f"  {isolated_str}")
            ascii_lines.append("")
        
        return "\n".join(ascii_lines)

    def _draw_small_network(self, network, component) -> List[str]:
        """Draw ASCII diagram for small networks (≤8 waters)."""
        import networkx as nx
        
        lines = []
        subgraph = network.subgraph(component)
        
        # For very small networks, show simple connectivity
        if len(component) <= 4:
            # Linear or simple star layout
            nodes = list(component)
            
            if len(nodes) == 2:
                lines.append(f"    HOH{nodes[0]} ═══ HOH{nodes[1]}")
            elif len(nodes) == 3:
                # Check if it's a line or triangle
                edges = list(subgraph.edges())
                if len(edges) == 2:
                    # Linear: A-B-C
                    center = None
                    for node in nodes:
                        if subgraph.degree(node) == 2:
                            center = node
                            break
                    if center:
                        neighbors = list(subgraph.neighbors(center))
                        lines.append(f"    HOH{neighbors[0]} ═══ HOH{center} ═══ HOH{neighbors[1]}")
                    else:
                        # Fallback
                        lines.append(f"    HOH{edges[0][0]} ═══ HOH{edges[0][1]} ═══ HOH{edges[1][1]}")
                else:
                    # Triangle
                    lines.append(f"    HOH{nodes[0]} ═══ HOH{nodes[1]}")
                    lines.append(f"      ║           ║")
                    lines.append(f"    HOH{nodes[2]} ═══════════")
            else:
                # For 4+ nodes, show edge list
                edges = list(subgraph.edges())
                lines.append("    Connections:")
                for edge in edges:
                    lines.append(f"      HOH{edge[0]} ═══ HOH{edge[1]}")
        else:
            # For larger small networks (5-8), show edge list format
            edges = list(subgraph.edges())
            lines.append("    Water connections:")
            for edge in edges:
                lines.append(f"      HOH{edge[0]} ═══ HOH{edge[1]}")
        
        return lines

    def _draw_large_network_summary(self, network, component) -> List[str]:
        """Summarize large networks with statistics."""
        import networkx as nx
        
        lines = []
        subgraph = network.subgraph(component)
        
        # Basic stats
        lines.append(f"    Waters: {len(component)}")
        lines.append(f"    Connections: {subgraph.number_of_edges()}")
        
        # Find hubs in this component
        local_hubs = [(node, subgraph.degree(node)) for node in component if subgraph.degree(node) >= 3]
        local_hubs.sort(key=lambda x: x[1], reverse=True)
        
        if local_hubs:
            lines.append(f"    Hub waters: {', '.join([f'HOH{node}({deg})' for node, deg in local_hubs[:3]])}")
        
        # Check if it's tree-like or has cycles
        if nx.is_tree(subgraph):
            lines.append("    Structure: Chain-like (no cycles)")
            try:
                diameter = nx.diameter(subgraph)
                lines.append(f"    Longest path: {diameter} connections")
            except:
                pass
        else:
            # Count cycles
            try:
                cycles = nx.minimum_cycle_basis(subgraph)
                if cycles:
                    lines.append(f"    Structure: Network with {len(cycles)} cycle(s)")
                else:
                    lines.append("    Structure: Complex network")
            except:
                lines.append("    Structure: Complex network")
        
        return lines

    def analyze_water_with_networks(self, water_residues: List[Residue]) -> Dict[str, Any]:
        """Comprehensive water analysis including network topology."""
        # Start with individual water analysis
        water_analyses = []
        for water_residue in water_residues:
            analysis = self.analyze_water_comprehensive(water_residue)
            water_analyses.append(analysis)
        
        # Add network analysis
        network_analysis = self.analyze_water_networks(water_residues)
        
        # Generate network ASCII
        network_ascii = self.generate_network_ascii(network_analysis)
        
        return {
            'individual_analyses': water_analyses,
            'network_analysis': network_analysis,
            'network_ascii': network_ascii,
            'total_waters': len(water_residues)
        }

class ComponentClassifier:
    """Static class to classify standard components with improved special residue handling"""

    PROTEIN_RESIDUES = {
        "GLY",
        "ALA",
        "VAL",
        "LEU",
        "ILE",
        "PRO",
        "PHE",
        "TYR",
        "TRP",
        "SER",
        "THR",
        "CYS",
        "MET",
        "ASN",
        "GLN",
        "ASP",
        "GLU",
        "LYS",
        "ARG",
        "HIS",
        "MSE",
    }

    DNA_RESIDUES = {"DA", "DC", "DG", "DT"}
    RNA_RESIDUES = {"A", "C", "G", "U"}
    WATER_RESIDUES = {"HOH", "WAT"}

    # Component-type keys produced by classify_residue for the four structural
    # classes. Everything else is a "special" type whose key is the CCD chemical
    # name (or "unknown (RESNAME)") and is displayed verbatim.
    STANDARD_TYPES = ["amino_acid", "dna_base", "rna_base", "water"]

    # Keys whose human-readable label is not just the title-cased key.
    DISPLAY_OVERRIDES = {"amino_acid": "Protein"}

    _ccd_cache: Dict[str, Union[str, Dict[str, str]]] = {}

    @classmethod
    def display_name(cls, comp_type: str) -> str:
        """
        Human-readable label for a component type key.

        Single source of truth for how component types are shown to the user,
        so the composition table, the per-component prompts and the saved
        filter-selection summary all agree.

        Args:
            comp_type: Component type key from classify_residue

        Returns:
            Display label, e.g. "Protein", "Dna Base", "FE (III) ION"
        """
        if comp_type in cls.DISPLAY_OVERRIDES:
            return cls.DISPLAY_OVERRIDES[comp_type]
        if comp_type in cls.STANDARD_TYPES:
            return comp_type.replace("_", " ").title()
        # Special types are already human-readable chemical names; recasing
        # them mangles CCD names like "PROTOPORPHYRIN IX CONTAINING FE".
        return comp_type

    @classmethod
    def classify_residue(cls, residue: Residue, ccd_parser: Any = None) -> str:
        """
        Classify a residue into component type with caching.

        Args:
            residue: Biopython Residue object
            ccd_parser: CCD parser for special residue lookup

        Returns:
            Component type string
        """
        res_name = residue.resname.strip().upper()

        if res_name in cls._ccd_cache:
            return cls._ccd_cache[res_name]

        if res_name in cls.PROTEIN_RESIDUES:
            cls._ccd_cache[res_name] = "amino_acid"
            return "amino_acid"
        elif res_name in cls.DNA_RESIDUES:
            cls._ccd_cache[res_name] = "dna_base"
            return "dna_base"
        elif res_name in cls.RNA_RESIDUES:
            cls._ccd_cache[res_name] = "rna_base"
            return "rna_base"
        elif res_name in cls.WATER_RESIDUES:
            cls._ccd_cache[res_name] = "water"
            return "water"

        if ccd_parser:
            try:
                residue_data = ccd_parser.get_residue_data(res_name)

                if "error" in residue_data:
                    result = f"unknown ({res_name})"
                    cls._ccd_cache[res_name] = result
                    return result

                ccd_name = residue_data.get("name", f"unknown ({res_name})")
                cls._ccd_cache[res_name] = ccd_name
                return ccd_name

            except Exception:
                result = f"unknown ({res_name})"
                cls._ccd_cache[res_name] = result
                return result

        result = f"unknown ({res_name})"
        cls._ccd_cache[res_name] = result
        return result

    @classmethod
    def classify_residue_hplusplus(cls, residue: Residue) -> str:
        """
        Classify a residue from H++ structures using backbone atom detection.

        H++ structures use AMBER forcefield naming conventions (HID, HIE, HIP, CYX, etc.)
        which may not be recognized by CCD. This method uses backbone atoms (N, CA, C, O)
        to identify amino acids, and falls back to standard residue lists for others.

        Args:
            residue: Biopython Residue object

        Returns:
            Component type string
        """
        res_name = residue.resname.strip().upper()

        # Check cache first
        if res_name in cls._ccd_cache:
            return cls._ccd_cache[res_name]

        # Standard residue checks (water, DNA, RNA work with AMBER naming too)
        if res_name in cls.WATER_RESIDUES:
            cls._ccd_cache[res_name] = "water"
            return "water"
        elif res_name in cls.DNA_RESIDUES:
            cls._ccd_cache[res_name] = "dna_base"
            return "dna_base"
        elif res_name in cls.RNA_RESIDUES:
            cls._ccd_cache[res_name] = "rna_base"
            return "rna_base"

        # Backbone-based amino acid detection (works for AMBER-named residues)
        atom_names = {atom.name.strip().upper() for atom in residue.get_atoms()}
        if {'N', 'CA', 'C', 'O'}.issubset(atom_names):
            cls._ccd_cache[res_name] = "amino_acid"
            return "amino_acid"

        # Fallback: use residue name as component type
        result = res_name
        cls._ccd_cache[res_name] = result
        return result


class ChainTopologyAnalyzer:
    """Analyzes protein chain connectivity to determine topology and linear sequences."""
    
    def __init__(self, interface_data: Dict[str, List[str]], interface_areas: Dict[str, Dict[str, float]] = None):
        """
        Initialize with chain interface data.
        
        Args:
            interface_data: Dictionary mapping chain IDs to lists of interfacing chains
            interface_areas: Optional buried surface area data for weighting
        """
        self.interface_data = interface_data
        self.interface_areas = interface_areas or {}
        
        if NETWORKX_AVAILABLE:
            self.graph = self._build_networkx_graph()
        else:
            self.graph = self._build_simple_graph()
    
    def _build_networkx_graph(self):
        """Build a NetworkX graph from interface data."""
        G = nx.Graph()
        
        # Add all chains as nodes
        for chain_id in self.interface_data:
            G.add_node(chain_id)
        
        # Add edges for interfaces with weights (BSA)
        for chain_id, neighbors in self.interface_data.items():
            for neighbor in neighbors:
                weight = self.interface_areas.get(chain_id, {}).get(neighbor, 1.0)
                G.add_edge(chain_id, neighbor, weight=weight)
        
        return G
    
    def _build_simple_graph(self):
        """Build a simple adjacency list representation."""
        return {chain_id: set(neighbors) for chain_id, neighbors in self.interface_data.items()}
    
    def determine_linear_sequence(self) -> Optional[List[str]]:
        """
        Determine the linear sequence if the structure is linear or has a main chain.
        
        Returns:
            List of chain IDs in linear order, or None if not applicable
        """
        if NETWORKX_AVAILABLE:
            return self._determine_linear_sequence_nx()
        else:
            return self._determine_linear_sequence_simple()
    
    def _determine_linear_sequence_nx(self) -> Optional[List[str]]:
        """NetworkX-based linear sequence determination."""
        if not nx.is_connected(self.graph):
            return None
        
        # Find terminal nodes (degree 1)
        terminal_nodes = [node for node, degree in self.graph.degree() if degree == 1]
        
        if len(terminal_nodes) < 2:
            return None
        
        # Find longest path between any two terminals
        best_path = []
        for i in range(len(terminal_nodes)):
            for j in range(i + 1, len(terminal_nodes)):
                try:
                    path = nx.shortest_path(self.graph, terminal_nodes[i], terminal_nodes[j])
                    if len(path) > len(best_path):
                        best_path = path
                except nx.NetworkXNoPath:
                    continue
        
        return best_path if best_path else None
    
    def _determine_linear_sequence_simple(self) -> Optional[List[str]]:
        """Simple algorithm without NetworkX."""
        # Find terminal nodes (degree 1)
        terminal_nodes = [chain for chain, neighbors in self.graph.items() if len(neighbors) == 1]
        
        if len(terminal_nodes) < 2:
            return None
        
        # Use BFS to find path between terminals
        def find_path(start, end, graph):
            queue = deque([(start, [start])])
            visited = {start}
            
            while queue:
                current, path = queue.popleft()
                
                if current == end:
                    return path
                
                for neighbor in graph[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append((neighbor, path + [neighbor]))
            
            return None
        
        # Find longest path between any two terminals
        best_path = []
        for i in range(len(terminal_nodes)):
            for j in range(i + 1, len(terminal_nodes)):
                path = find_path(terminal_nodes[i], terminal_nodes[j], self.graph)
                if path and len(path) > len(best_path):
                    best_path = path
        
        return best_path if best_path else None
    
    def get_topology_info(self) -> Dict[str, any]:
        """Get basic topology information."""
        if NETWORKX_AVAILABLE:
            return self._get_topology_info_nx()
        else:
            return self._get_topology_info_simple()
    
    def _get_topology_info_nx(self) -> Dict[str, any]:
        """NetworkX-based topology analysis."""
        terminal_nodes = [node for node, degree in self.graph.degree() if degree == 1]
        hub_nodes = [node for node, degree in self.graph.degree() if degree > 2]
        num_nodes = self.graph.number_of_nodes()
        num_edges = self.graph.number_of_edges()

        # For monomeric structures (single chain, no interfaces), is_connected should be False
        is_connected = nx.is_connected(self.graph) if num_nodes > 1 or num_edges > 0 else False

        return {
            'is_connected': is_connected,
            'num_chains': num_nodes,
            'num_interfaces': num_edges,
            'terminal_chains': terminal_nodes,
            'branching_points': hub_nodes,
            'topology_type': self._classify_topology_nx()
        }
    
    def _get_topology_info_simple(self) -> Dict[str, any]:
        """Simple topology analysis."""
        terminal_nodes = [chain for chain, neighbors in self.graph.items() if len(neighbors) == 1]
        hub_nodes = [chain for chain, neighbors in self.graph.items() if len(neighbors) > 2]
        num_nodes = len(self.graph)
        num_edges = sum(len(neighbors) for neighbors in self.graph.values()) // 2

        # For monomeric structures (single chain, no interfaces), is_connected should be False
        is_connected = self._is_connected_simple() if num_nodes > 1 or num_edges > 0 else False

        return {
            'is_connected': is_connected,
            'num_chains': num_nodes,
            'num_interfaces': num_edges,
            'terminal_chains': terminal_nodes,
            'branching_points': hub_nodes,
            'topology_type': self._classify_topology_simple()
        }
    
    def _is_connected_simple(self) -> bool:
        """Check if graph is connected using BFS."""
        if not self.graph:
            return True
        
        start = next(iter(self.graph))
        visited = set()
        queue = deque([start])
        
        while queue:
            current = queue.popleft()
            if current in visited:
                continue
            visited.add(current)
            queue.extend(self.graph[current] - visited)
        
        return len(visited) == len(self.graph)
    
    def _classify_topology_nx(self) -> str:
        """Classify topology using NetworkX."""
        num_nodes = self.graph.number_of_nodes()
        num_edges = self.graph.number_of_edges()

        # Single chain with no interfaces = monomeric
        if num_nodes == 1 and num_edges == 0:
            return 'monomeric'

        if not nx.is_connected(self.graph):
            return 'disconnected'

        # Check for actual cycles in the graph
        has_cycle = False
        try:
            # For undirected graphs, find_cycle will detect any cycle
            nx.find_cycle(self.graph)
            has_cycle = True
        except nx.NetworkXNoCycle:
            has_cycle = False

        # Classify based on branching and cycles
        hub_nodes = [node for node, degree in self.graph.degree() if degree > 2]

        if has_cycle:
            # There's a cycle - could be simple cycle or complex with branches
            if len(hub_nodes) == 0:
                # Simple cycle with no branching (all nodes degree 2)
                return 'cyclic'
            else:
                # Cycle with branches
                return 'complex_cyclic'
        else:
            # No cycle - linear or branched tree
            if len(hub_nodes) == 0:
                return 'linear'
            elif len(hub_nodes) == 1:
                return 'star_branched'
            else:
                return 'multi_branched'
    
    def _classify_topology_simple(self) -> str:
        """Classify topology using simple algorithms."""
        num_nodes = len(self.graph)
        num_edges = sum(len(neighbors) for neighbors in self.graph.values()) // 2

        # Single chain with no interfaces = monomeric
        if num_nodes == 1 and num_edges == 0:
            return 'monomeric'

        if not self._is_connected_simple():
            return 'disconnected'

        # Simple cycle detection using DFS
        has_cycle = self._has_cycle_simple()

        hub_nodes = [chain for chain, neighbors in self.graph.items() if len(neighbors) > 2]

        if has_cycle:
            # There's a cycle
            if len(hub_nodes) == 0:
                return 'cyclic'
            else:
                return 'complex_cyclic'
        else:
            # No cycle - linear or branched tree
            if len(hub_nodes) == 0:
                return 'linear'
            elif len(hub_nodes) == 1:
                return 'star_branched'
            else:
                return 'multi_branched'

    def _has_cycle_simple(self) -> bool:
        """Detect if graph has a cycle using DFS."""
        if not self.graph:
            return False

        visited = set()

        def dfs(node, parent):
            visited.add(node)
            for neighbor in self.graph[node]:
                if neighbor not in visited:
                    if dfs(neighbor, node):
                        return True
                elif neighbor != parent:
                    # Visited neighbor that's not the parent = cycle found
                    return True
            return False

        # Check from first node
        start = next(iter(self.graph))
        return dfs(start, None)

    def draw_topology_ascii(self) -> str:
        """
        Generate ASCII art representation of the chain topology.
        
        Returns:
            String containing ASCII art representation
        """
        if NETWORKX_AVAILABLE:
            return self._draw_topology_ascii_nx()
        else:
            return self._draw_topology_ascii_simple()

    def _draw_topology_ascii_nx(self) -> str:
        """NetworkX-based ASCII topology drawing."""
        import io
        
        output = io.StringIO()
        
        # Get topology info
        topology_info = self.get_topology_info()
        topology_type = topology_info['topology_type']
        
        if topology_type == 'linear':
            return self._draw_linear_ascii()
        elif topology_type == 'star_branched':
            return self._draw_star_ascii()
        elif topology_type == 'multi_branched':
            return self._draw_tree_ascii()
        elif topology_type == 'cyclic':
            return self._draw_cyclic_ascii()
        elif topology_type == 'complex_cyclic':
            return self._draw_complex_cyclic_ascii()
        else:
            return self._draw_generic_ascii()

    def _draw_linear_ascii(self) -> str:
        """Draw linear chain arrangement."""
        linear_seq = self.determine_linear_sequence()
        if not linear_seq:
            return "Unable to determine linear sequence"
        
        # Create horizontal linear representation
        lines = []
        lines.append("Linear Chain Topology:")
        lines.append("")
        
        # Top line with chain IDs
        chain_line = ""
        connector_line = ""
        
        for i, chain in enumerate(linear_seq):
            if i == 0:
                chain_line += f"[{chain}]"
                connector_line += " │ "
            else:
                chain_line += f"───[{chain}]"
                connector_line += "    │ "
        
        lines.append(chain_line)
        lines.append("")
        
        # Add interface strengths if available
        if self.interface_areas:
            strength_line = ""
            for i in range(len(linear_seq) - 1):
                chain1, chain2 = linear_seq[i], linear_seq[i + 1]
                area = self.interface_areas.get(chain1, {}).get(chain2, 0)
                strength = self._get_strength_symbol(area)
                if i == 0:
                    strength_line += f"  {strength}  "
                else:
                    strength_line += f"  {strength}  "
            lines.append(f"Interfaces: {strength_line}")
            lines.append("")
        
        return "\n".join(lines)

    def _draw_star_ascii(self) -> str:
        """Draw star/hub topology."""
        lines = []
        lines.append("Star Topology:")
        lines.append("")
        
        # Find the hub (node with highest degree)
        hub = None
        max_degree = 0
        
        if NETWORKX_AVAILABLE:
            for node, degree in self.graph.degree():
                if degree > max_degree:
                    max_degree = degree
                    hub = node
        else:
            for node, neighbors in self.graph.items():
                if len(neighbors) > max_degree:
                    max_degree = len(neighbors)
                    hub = node
        
        if not hub:
            return "Unable to identify hub"
        
        # Get connected chains
        if NETWORKX_AVAILABLE:
            connected = list(self.graph.neighbors(hub))
        else:
            connected = list(self.graph[hub])
        
        # Create star pattern
        if len(connected) <= 4:
            # Simple cross pattern
            lines.append("    [{}]".format(connected[0] if len(connected) > 0 else " "))
            lines.append("     │")
            hub_line = ""
            if len(connected) > 3:
                hub_line += f"[{connected[3]}]───"
            else:
                hub_line += "    "
            hub_line += f"[{hub}]"
            if len(connected) > 1:
                hub_line += f"───[{connected[1]}]"
            lines.append(hub_line)
            if len(connected) > 2:
                lines.append("     │")
                lines.append("    [{}]".format(connected[2]))
        else:
            # Circular arrangement for many connections
            lines.append(f"Hub [{hub}] connected to:")
            for i, chain in enumerate(connected):
                strength = ""
                if self.interface_areas:
                    area = self.interface_areas.get(hub, {}).get(chain, 0)
                    strength = f" {self._get_strength_symbol(area)}"
                lines.append(f"  ├─{strength}─ [{chain}]")
        
        lines.append("")
        return "\n".join(lines)

    def _draw_tree_ascii(self) -> str:
        """Draw tree/branched topology."""
        lines = []
        lines.append("Branched Tree Topology:")
        lines.append("")
        
        # Find longest path as main backbone
        linear_seq = self.determine_linear_sequence()
        if not linear_seq:
            return self._draw_generic_ascii()
        
        # Draw main backbone horizontally
        backbone_line = ""
        for i, chain in enumerate(linear_seq):
            if i == 0:
                backbone_line += f"[{chain}]"
            else:
                backbone_line += f"───[{chain}]"
        
        lines.append("Main backbone:")
        lines.append(backbone_line)
        lines.append("")
        
        # Find and draw branches
        branches_found = False
        
        for i, chain in enumerate(linear_seq):
            # Find what else is connected to this chain
            if NETWORKX_AVAILABLE:
                neighbors = set(self.graph.neighbors(chain))
            else:
                neighbors = set(self.graph[chain])
            
            # Remove backbone neighbors
            backbone_neighbors = set()
            if i > 0:
                backbone_neighbors.add(linear_seq[i-1])
            if i < len(linear_seq) - 1:
                backbone_neighbors.add(linear_seq[i+1])
            
            branch_chains = neighbors - backbone_neighbors
            
            if branch_chains:
                branches_found = True
                lines.append(f"Branches from [{chain}]:")
                for branch in sorted(branch_chains):
                    strength = ""
                    if self.interface_areas:
                        area = self.interface_areas.get(chain, {}).get(branch, 0)
                        strength = f" {self._get_strength_symbol(area)}"
                    lines.append(f"  │")
                    lines.append(f"  └─{strength}─ [{branch}]")
                lines.append("")
        
        if not branches_found:
            lines.append("No branches detected")
        
        return "\n".join(lines)

    def _draw_cyclic_ascii(self) -> str:
        """Draw cyclic topology."""
        lines = []
        lines.append("Cyclic Topology:")
        lines.append("")
        
        # For simple cycles, try to arrange in a circle
        chains = list(self.interface_data.keys())
        
        if len(chains) <= 6:
            # Small cycle - draw as polygon
            if len(chains) == 3:
                lines.extend([
                    "    [{}]".format(chains[0]),
                    "   ╱   ╲",
                    "[{}]─────[{}]".format(chains[1], chains[2])
                ])
            elif len(chains) == 4:
                lines.extend([
                    "[{}]─────[{}]".format(chains[0], chains[1]),
                    " │       │",
                    " │       │", 
                    "[{}]─────[{}]".format(chains[3], chains[2])
                ])
            else:
                # Larger cycles - list format
                lines.append("Cycle order:")
                cycle_line = " → ".join(chains) + " → " + chains[0]
                lines.append(cycle_line)
        else:
            lines.append(f"Large cycle with {len(chains)} chains")
            lines.append("Partial representation:")
            for i, chain in enumerate(chains[:6]):
                lines.append(f"[{chain}] → [{chains[(i+1) % len(chains)]}]")
            if len(chains) > 6:
                lines.append("...")
        
        lines.append("")
        return "\n".join(lines)

    def _draw_complex_cyclic_ascii(self) -> str:
        """Draw complex cyclic topology (cycle with branches)."""
        lines = []
        lines.append("Complex Cyclic Topology (contains cycles with branching):")
        lines.append("")

        # Show basic statistics
        chains = list(self.interface_data.keys())
        hub_nodes = [chain for chain, neighbors in self.interface_data.items() if len(neighbors) > 2]

        lines.append(f"Total chains: {len(chains)}")
        lines.append(f"Branching points: {len(hub_nodes)} ({', '.join(hub_nodes) if hub_nodes else 'none'})")
        lines.append("")

        # Show connectivity for all chains
        lines.append("Chain connectivity:")
        for chain, neighbors in self.interface_data.items():
            if neighbors:
                neighbor_str = ", ".join(sorted(neighbors))
                lines.append(f"  [{chain}] ↔ {neighbor_str}")

        lines.append("")
        return "\n".join(lines)

    def _draw_generic_ascii(self) -> str:
        """Draw generic/complex topology."""
        lines = []
        lines.append("Complex Topology:")
        lines.append("")
        
        # Just list all connections
        for chain, neighbors in self.interface_data.items():
            if neighbors:
                lines.append(f"[{chain}] connected to:")
                for neighbor in sorted(neighbors):
                    strength = ""
                    if self.interface_areas:
                        area = self.interface_areas.get(chain, {}).get(neighbor, 0)
                        strength = f" {self._get_strength_symbol(area)}"
                    lines.append(f"  ├─{strength}─ [{neighbor}]")
            else:
                lines.append(f"[{chain}] (isolated)")
            lines.append("")
        
        return "\n".join(lines)

    def _draw_topology_ascii_simple(self) -> str:
        """Simple ASCII drawing without NetworkX."""
        # Use the same logic but without NetworkX-specific functions
        topology_info = self.get_topology_info()
        
        if topology_info['topology_type'] == 'linear':
            return self._draw_linear_ascii()
        else:
            return self._draw_generic_ascii()

    def _get_strength_symbol(self, area: float) -> str:
        """Convert interface area to visual strength symbol."""
        if area > 2000:
            return "═══"  # Very strong
        elif area > 1000:
            return "━━━"  # Strong
        elif area > 500:
            return "───"  # Medium
        elif area > 200:
            return "┄┄┄"  # Weak
        else:
            return "···"  # Very weak


# Viewer annotation keys used for water halos. Every key is cleared before a
# new set is drawn, so a category that disappears after a parameter change
# disappears from the viewer too.
WATER_HALO_KEYS = (
    "water_metal", "water_interface", "water_b_below_median",
    "water_hbond_1", "water_hbond_2", "water_hbond_3", "water_hbond_4",
    "water_burial_clash", "water_burial_enclosed", "water_burial_zero",
    "water_burial_covered90", "water_burial_covered50",
)


class PDBFilterWorker:
    """Handles PDB structure filtering operations."""

    def __init__(self, filename: str, existing_structure: Optional[Structure] = None, processor=None):
        """
        Initialize PDB filter worker.

        Args:
            filename: Path to the PDB file
            existing_structure: Optional pre-loaded structure object
            processor: Optional ProPrep processor for session recording context
        """
        self.filename = filename
        self.processor = processor
        self.parser = PDBParser(QUIET=True)
        self.ccd_parser = ccd.CCDParser(use_cache=True)
        self.filtered_structure = None
        self.filter_selections = {}
        self.interface_areas = {}

        # Detect H++ structures (AMBER forcefield naming, hydrogens already added)
        self.is_hplusplus_structure = self._detect_hplusplus_structure(filename)

        try:
            if existing_structure:
                self.structure = existing_structure
            else:
                self.structure = self.parser.get_structure("protein", filename)
        except Exception as e:
            logger.error(f"Error loading PDB file: {e}")
            raise

    def _parse_selection_string(self, selection: str) -> List[int]:
        """Parse selection string like '1,3,5-8,10' into list of indices."""
        indices = []
        parts = selection.split(',')
        
        for part in parts:
            part = part.strip()
            if '-' in part:
                start, end = map(int, part.split('-'))
                indices.extend(range(start, end + 1))
            else:
                indices.append(int(part))
                
        return sorted(set(indices))

    def filter_water_with_analysis(self, chain: Chain, water_residues: List[Residue], console=None) -> Set[int]:
        """Water filtering with structural analysis and improved navigation."""
        if not console:
            # Fallback to simple filtering if no console
            return {residue.id[1] for residue in water_residues}
            
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        from rich.panel import Panel
        
        if not water_residues:
            console.print("[yellow]No water molecules found in this chain.[/yellow]")
            return set()
            
        console.print(f"\n[bold]Water Analysis for Chain {chain.id}[/bold]")
        console.print(f"Found {len(water_residues)} water molecules")
        
        # Configure analysis parameters
        analyzer = WaterAnalyzer(self.structure, model_idx=0)
        
        # Show initial overview
        self._display_analysis_overview(analyzer, console)
        
        # MAIN NAVIGATION LOOP - allows returning to metrics selection
        while True:
            # Select metrics to display
            metrics_available = {
                '1': 'Distance to nearest metal ion',
                '2': 'Hydrogen bond analysis',
                '3': 'B-factor/thermal information',
                '4': 'Water burial: solvent-accessible area and enclosure',
                '5': 'Interface proximity',
                '6': 'Multi-radius atom-count profiling',
                '7': 'Directional atom-count analysis (8 sectors)',
                '8': 'Water network analysis (clusters and chains)',
                '9': 'Return with no waters selected and proceed with component filtering'
            }

            console.print("\n[bold]Select Analysis Metrics to Display:[/bold]")
            for key, desc in metrics_available.items():
                console.print(f"  {key:>2}. {desc}")

            metric_choice = prompt_with_context(
                self.processor,
                "Select metrics to display (comma-separated or single choice)",
                choices=list(metrics_available.keys()),  # Add choices validation
                default="6",
                module="PDB Filter - Water Analysis",
                description="Select water analysis metrics to display",
                options_map=metrics_available,
            )

            # Handle return to component menu
            if metric_choice == '9':
                return set()

            # Parse selected metrics (comma-separated list)
            selected_metrics = metric_choice.split(',')
                            
            # Show parameters relevant to selected metrics
            has_adjustable_params = self._display_selected_parameters(analyzer, selected_metrics, console)

            # Analysis loop - allow parameter adjustment and re-analysis
            while True:
                if has_adjustable_params and confirm_with_context(
                    self.processor,
                    "Would you like to modify any parameters?",
                    default=False,
                    module="PDB Filter - Water Analysis",
                    description="Modify water analysis parameters",
                ):
                    self._modify_analysis_parameters(analyzer, selected_metrics, console)
                    # Show updated parameters after modification
                    has_adjustable_params = self._display_selected_parameters(analyzer, selected_metrics, console)
                else:
                    break

            # Run analysis and display results
            analysis_loop_active = True
            while analysis_loop_active:
                console.print("\n[bold blue]Analyzing water molecules...[/bold blue]")
                
                # Show explanation for advanced analysis if selected
                if '6' in selected_metrics and len(selected_metrics) == 1:
                    self._display_multiradius_explanation(analyzer, console)
                    if not confirm_with_context(
                        self.processor,
                        "Proceed with multi-radius profiling?",
                        default=True,
                        module="PDB Filter - Water Analysis",
                        description="Proceed with multi-radius profiling",
                    ):
                        analysis_loop_active = False
                        break  # Exit analysis loop and return to metrics selection

                if '7' in selected_metrics and len(selected_metrics) == 1:
                    self._display_directional_explanation(analyzer, console)
                    if not confirm_with_context(
                        self.processor,
                        "Proceed with directional analysis?",
                        default=True,
                        module="PDB Filter - Water Analysis",
                        description="Proceed with directional analysis",
                    ):
                        analysis_loop_active = False
                        break  # Exit analysis loop and return to metrics selection

                if '8' in selected_metrics and len(selected_metrics) == 1:
                    self._display_network_explanation(analyzer, console)
                    if not confirm_with_context(
                        self.processor,
                        "Proceed with water network analysis?",
                        default=True,
                        module="PDB Filter - Water Analysis",
                        description="Proceed with water network analysis",
                    ):
                        analysis_loop_active = False
                        break  # Exit analysis loop and return to metrics selection

                # Get interface data if needed
                interface_data = {}
                if '5' in selected_metrics:
                    interface_data = analyzer.calculate_interfaces()
                    
                # Determine which analysis method to use based on selected metrics
                if '6' in selected_metrics or '7' in selected_metrics or '8' in selected_metrics or '9' in selected_metrics or '10' in selected_metrics:
                    # Comprehensive analysis needed for advanced features
                    console.print("[blue]Running comprehensive analysis (includes profiling, directional, and network data)...[/blue]")
                    
                    if '8' in selected_metrics:
                        # Network analysis requires special handling
                        comprehensive_analysis = analyzer.analyze_water_with_networks(water_residues)
                        water_analyses = comprehensive_analysis['individual_analyses']
                        network_data = comprehensive_analysis
                    else:
                        water_analyses = []
                        for water in water_residues:
                            analysis = analyzer.analyze_water_comprehensive(water)
                            if analysis:
                                water_analyses.append(analysis)
                        network_data = None
                else:
                    # Standard analysis
                    water_analyses = []
                    network_data = None
                    for water in water_residues:
                        analysis = analyzer.analyze_water(water, interface_data)
                        if analysis:
                            water_analyses.append(analysis)
                        
                # Display results table
                self._display_water_analysis_table(water_analyses, selected_metrics, console, analyzer)

                # Halo the waters so the user can read the tables and see the
                # spatial distribution at the same time: one group per fact a
                # displayed metric establishes, each labelled with its rule and
                # cutoff, overlaps allowed (no combined category).
                try:
                    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                    for lbl in WATER_HALO_KEYS:
                        _viewer.unhighlight(lbl)
                    for key, display, color, items in self._water_halo_groups(
                            water_analyses, selected_metrics, analyzer):
                        clauses = [
                            f"(:{a.get('chain_id', chain.id)} and {a['residue_number']})"
                            for a in items
                        ]
                        if clauses:
                            _viewer.highlight(
                                " or ".join(clauses),
                                style="ball+stick",
                                color=color,
                                label=key,
                                display_label=display,
                            )
                except Exception:
                    pass

                # NEW: Add display logic for advanced metrics
                if '6' in selected_metrics:
                    self._display_burial_profile_table(water_analyses, console)

                if '7' in selected_metrics:
                    self._display_directional_analysis_table(water_analyses, console)
                    
                if '8' in selected_metrics:
                    self._display_water_network_analysis(network_data, console)

                # UPDATED OPTIONS MENU with better navigation
                console.print("\n[bold]Analysis Options:[/bold]")
                console.print("1. Adjust analysis parameters and re-run", highlight=False)
                console.print("2. Return to metrics selection menu", highlight=False)  # NEW OPTION
                console.print("3. Proceed with water selection", highlight=False)       # RENUMBERED
                console.print("4. Cancel and return to component menu", highlight=False) # RENUMBERED
                
                choice = prompt_with_context(
                    self.processor,
                    "Enter choice",
                    choices=['1', '2', '3', '4'],
                    default='4',
                    module="PDB Filter - Water Analysis",
                    description="Water analysis options",
                    options_map={
                        "1": "Adjust analysis parameters and re-run",
                        "2": "Return to metrics selection menu",
                        "3": "Proceed with water selection",
                        "4": "Cancel and return to component menu",
                    },
                )

                if choice == '1':
                    self._modify_analysis_parameters(analyzer, selected_metrics, console)
                    continue
                elif choice == '2':
                    # Return to metrics selection - break out of analysis loop
                    analysis_loop_active = False
                    break
                elif choice == '3':
                    # Proceed with water selection
                    return self._handle_final_water_selection(water_residues, water_analyses, console)
                elif choice == '4':
                    return set()  # Cancel and return to component menu
            
            # If we get here, user chose option 3 (return to metrics selection)
            # The main while loop will restart metrics selection

    def _handle_final_water_selection(self, water_residues, water_analyses, console):
        """Handle the final water selection step using residue numbers."""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        from rich.table import Table
        
        console.print("\n[bold]Water Selection[/bold]")
        
        # Show a summary table of all analyzed waters with their residue numbers
        summary_table = Table(title="Summary of Analyzed Waters", show_lines=True)
        summary_table.add_column("Residue #", style="bold yellow", width=10)
        summary_table.add_column("Name", style="yellow", width=8)
        summary_table.add_column("Chain", style="cyan", width=8)
        summary_table.add_column("Quick Info", style="grey50", width=30)
        
        # Sort by residue number for easy reference
        sorted_analyses = sorted(water_analyses, key=lambda x: x['residue_number'])
        
        for analysis in sorted_analyses:
            # Create quick info summary
            info_parts = []
            
            if analysis.get('metal_distance') and analysis['metal_distance'] <= 5.0:
                info_parts.append(f"Metal: {analysis['metal_distance']:.1f}Å")
            
            if analysis.get('total', 0) > 0:
                info_parts.append(f"H-bonds: {analysis['total']}")
            
            if analysis.get('b_factor') is not None:
                ratio = analysis.get('b_factor_ratio')
                if ratio is not None:
                    info_parts.append(f"B {analysis['b_factor']:.0f} ({ratio:.1f}× protein median)")
                else:
                    info_parts.append(f"B {analysis['b_factor']:.0f}")
            
            if analysis.get('at_interface'):
                info_parts.append("Interface")
            
            quick_info = "; ".join(info_parts) if info_parts else "Standard water"
            
            summary_table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                analysis['residue_name'],
                analysis['chain_id'],
                quick_info
            )
        
#       console.print(summary_table)
        
        # Instructions using residue numbers
        console.print("\n[bold cyan]Selection Instructions:[/bold cyan]")
        console.print("• Use [bold]residue numbers[/bold] from the Residue # column above")
        console.print("• Examples: '215,847,102' (select specific waters)")
        console.print("• Examples: '100-150' (select range of residue numbers)")
        console.print("• Examples: '215,300-310,847' (mix specific numbers and ranges)")
        console.print("• Use 'all' to select all analyzed waters")
        console.print("• Use 'none' to select no waters")
        
        # Show available residue number range
        if sorted_analyses:
            all_res_nums = sorted([analysis['residue_number'] for analysis in sorted_analyses])
            
            if len(all_res_nums) <= 15:
                console.print(f"\n[grey50]Available: {', '.join(map(str, all_res_nums))}[/grey50]")
            else:
                min_res = min(all_res_nums)
                max_res = max(all_res_nums)
                console.print(f"\n[grey50]Available: {min_res} to {max_res} ({len(all_res_nums)} waters total)[/grey50]")
        
        while True:
            choice = prompt_with_context(
                self.processor,
                "\nSelect waters by residue numbers",
                default="all",
                module="PDB Filter - Water Selection",
                description="Select waters by residue numbers, ranges, 'all', or 'none'",
            )
            
            if choice.lower() == "all":
                selected_residues = {analysis['residue_number'] for analysis in water_analyses}
                break  # Valid selection, exit loop
            elif choice.lower() == "none":
                selected_residues = set()
                break  # Valid selection, exit loop
            else:
                try:
                    # Parse the selection (numbers and ranges)
                    selected_residue_numbers = self._parse_selection_string(choice)
                    
                    # Validate that these residue numbers exist in our analysis
                    available_residue_numbers = {analysis['residue_number'] for analysis in water_analyses}
                    invalid_numbers = [num for num in selected_residue_numbers if num not in available_residue_numbers]
                    
                    if invalid_numbers:
                        console.print(f"[red]Invalid residue numbers: {invalid_numbers}[/red]")
                        valid_numbers = sorted(available_residue_numbers)
                        if len(valid_numbers) <= 20:
                            console.print(f"[yellow]Available: {', '.join(map(str, valid_numbers))}[/yellow]")
                        else:
                            console.print(f"[yellow]Available: {min(valid_numbers)} to {max(valid_numbers)}[/yellow]")
                        console.print("[cyan]Please try again with valid residue numbers.[/cyan]")
                        continue  # Invalid selection, retry
                    
                    selected_residues = set(selected_residue_numbers)
                    break  # Valid selection, exit loop
                    
                except (ValueError, IndexError) as e:
                    console.print(f"[red]Invalid selection format: {e}[/red]")
                    console.print("[yellow]Please use residue numbers, ranges (100-150), 'all', or 'none'[/yellow]")
                    console.print("[cyan]Please try again.[/cyan]")
                    continue  # Invalid format, retry
                        
        # Confirm selection with detailed info
        kept_count = len(selected_residues)
        total_water_count = len(water_residues)
        discarded_count = total_water_count - kept_count
        
        console.print(f"\n[green]Keep {kept_count} waters, discard {discarded_count} waters.[/green]")
        
        if kept_count > 0:
            # Show which specific waters are being kept
            kept_waters = [analysis for analysis in sorted_analyses if analysis['residue_number'] in selected_residues]
            if len(kept_waters) <= 10:
                kept_list = [f"HOH{w['residue_number']}" for w in kept_waters]
                console.print(f"[grey50]Keeping: {', '.join(kept_list)}[/grey50]")
            else:
                kept_list = [f"HOH{w['residue_number']}" for w in kept_waters[:10]]
                console.print(f"[grey50]Keeping: {', '.join(kept_list)} and {len(kept_waters)-10} more...[/grey50]")
        
        # Replace the per-category water halos with one green ball+stick
        # rep showing exactly the kept set, so the user is confirming
        # against the visible structural answer rather than just a count.
        # Clear all category labels (they're stale once a final pick
        # exists), then halo just the kept residues. Pick chain id from
        # the first analysis since all waters in this call share a chain.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            for lbl in WATER_HALO_KEYS:
                _viewer.unhighlight(lbl)
            if selected_residues and water_analyses:
                chain_id = water_analyses[0].get('chain_id', '')
                clauses = [f"(:{chain_id} and {r})" for r in sorted(selected_residues)]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color="#33a02c",
                    label="water_final_keep",
                )
            else:
                _viewer.unhighlight("water_final_keep")
        except Exception:
            pass

        if confirm_with_context(
            self.processor,
            "Confirm selection?",
            default=True,
            module="PDB Filter - Water Selection",
            description="Confirm water residue selection",
        ):
            return selected_residues
        else:
            return set()

    def _display_analysis_overview(self, analyzer, console):
        """Display initial overview of all available analysis methods."""
        from rich.panel import Panel

        overview_text = "Available Analysis Methods:\n\n"
        overview_text += f"  • Metal Coordination: Distance to metal ions\n"
        overview_text += f"  • Hydrogen Bonding: N/O/S partners within a heavy-atom distance cutoff, closest first\n"
        overview_text += f"  • Thermal Factors: B-factor assessment for structural ordering\n"
        overview_text += f"  • Water Burial: Accessible area of each water oxygen (1.4 Å probe) and whether it is enclosed from bulk solvent\n"
        overview_text += f"  • Interface Proximity: Waters at protein-protein interfaces\n"
        overview_text += f"  • Multi-Radius Profiling: Atom count around each water vs radius, with ASCII charts\n"
        overview_text += f"  • Directional Analysis: Atom count in 8 compass sectors around each water\n"
        overview_text += f"  • Network Analysis: Water connectivity and clustering topology"  

        panel = Panel(overview_text, title="Water Analysis Overview", border_style="blue", expand=False)
        console.print(panel)
        
    def _display_selected_parameters(self, analyzer, selected_metrics, console):
        """Display parameters relevant only to selected metrics."""
        from rich.panel import Panel
        
        params_text = "Parameters for Selected Analyses:\n\n"
        has_parameters = False
        
        if '1' in selected_metrics:  # Metal coordination
            params_text += f"Metal Coordination Analysis:\n"
            params_text += f"  • Default coordination distance: {analyzer.parameters['metal_distance_cutoff']} Å (adjustable during analysis)\n"
            params_text += f"  • Dynamic cutoff adjustment changes which waters are analyzed\n\n"
            has_parameters = True
                        
        if '2' in selected_metrics:  # Hydrogen bonds
            params_text += f"Hydrogen Bond Analysis:\n"
            params_text += f"  • H-bond distance cutoff: {analyzer.parameters['hbond_distance_cutoff']} Å\n"
            params_text += f"  • H-bond atoms: {analyzer.parameters['hbond_atoms']}\n"
            params_text += f"  • Max H-bonds per water: {analyzer.parameters['max_hbonds_per_water']}\n"
            params_text += f"  • Ranking: by distance, closest first; no angular test (hydrogens are usually absent, so heavy atoms do not fix where they point)\n\n"
            has_parameters = True
            
        if '3' in selected_metrics:  # B-factor
            params_text += f"B-factor Analysis:\n"
            params_text += f"  • Uses atomic B-factors from PDB file (no adjustable parameters)\n\n"
            
        if '4' in selected_metrics:  # SASA/Burial
            params_text += f"Water Burial (solvent accessibility):\n"
            params_text += f"  • Probe radius: {analyzer.parameters['sasa_probe_radius']} Å (a water molecule)\n"
            params_text += f"  • Occluding atoms: {analyzer.parameters['burial_atom_types']} (other waters never occlude)\n"
            params_text += f"  • Atomic radii: Bondi; enclosure grid: {analyzer.ENCLOSURE_GRID_SPACING} Å\n\n"
            has_parameters = True
            
        if '5' in selected_metrics:  # Interface
            params_text += f"Interface Proximity Analysis:\n"
            params_text += f"  • Interface proximity cutoff: {analyzer.parameters['interface_distance_cutoff']} Å\n\n"
            has_parameters = True
            
        if '6' in selected_metrics:  # Multi-radius profiling
            params_text += f"Multi-Radius Profiling (Phase 2):\n"
            params_text += f"  • Analysis method: Counts atoms at multiple radii (2.0-8.0 Å in 0.5 Å steps)\n"
            params_text += f"  • Atom types analyzed: {analyzer.parameters['burial_atom_types']}\n"
            params_text += f"  • Weighting scheme: {analyzer.parameters['burial_weighting']}\n"
            params_text += f"  • Saturation detection: <10% increase between consecutive radii\n"
            params_text += f"  • Steepest rise: Identifies radius range with maximum burial slope\n"
            params_text += f"  • Generates: ASCII profile charts showing burial vs radius\n\n"
            
        if '7' in selected_metrics:  # Directional analysis
            params_text += f"Directional Analysis (Phase 3):\n"
            params_text += f"  • Analysis method: 8-sector compass (N, NE, E, SE, S, SW, W, NW)\n"
            params_text += f"  • Uses burial radius: {analyzer.parameters['burial_radius']} Å\n"
            params_text += f"  • Atom types analyzed: {analyzer.parameters['burial_atom_types']}\n"
            params_text += f"  • Identifies: Primary burial direction and potential pocket openings\n"
            params_text += f"  • Generates: ASCII compass charts with directional burial data\n\n"
            
        if '8' in selected_metrics:  # Water network analysis
            params_text += f"Water Network Analysis:\n"
            params_text += f"  • Network type: {analyzer.parameters['network_type']}\n"
            params_text += f"    [grey50](Note: Only 'water_only' is currently implemented.\n"
            params_text += f"     Future: 'water_protein_water' and 'all_hbonds' network types)[/grey50]\n"
            params_text += f"  • H-bond distance cutoff: {analyzer.parameters['hbond_distance_cutoff']} Å\n"  # NEW
            params_text += f"  • H-bond atoms: {analyzer.parameters['hbond_atoms']}\n"  # NEW
            params_text += f"  • Minimum cluster size: {analyzer.parameters['network_min_cluster_size']}\n"
            params_text += f"  • Show isolated waters: {analyzer.parameters['network_show_isolated']}\n"
            params_text += f"  • Max detailed display size: {analyzer.parameters['network_max_display_size']}\n"
            params_text += f"  • Method: H-bond connectivity analysis\n\n"
            has_parameters = True   
                    
        # Handle case where no adjustable parameters exist
        if not has_parameters:
            params_text = "Selected analyses use fixed parameters or PDB file data.\n"
            params_text += "No user-adjustable parameters for current selection.\n\n"
            
            if '6' in selected_metrics or '7' in selected_metrics:
                params_text += "Advanced Analysis Details:\n\n"
                
            if '6' in selected_metrics:
                params_text += "Multi-Radius Profiling:\n"
                params_text += f"  • Analyzes burial at radii from 2.0 to 8.0 Å (0.5 Å steps)\n"
                params_text += f"  • Counts {analyzer.parameters['burial_atom_types']} atoms around each water\n"
                params_text += f"  • Saturation point: Where burial increase drops below 10%\n"
                params_text += f"  • Steepest rise: Radius range with maximum burial slope\n\n"
                
            if '7' in selected_metrics:
                params_text += "Directional Analysis:\n"
                params_text += f"  • Divides space into 8 compass sectors (45° each)\n"
                params_text += f"  • Uses {analyzer.parameters['burial_radius']} Å radius\n"
                params_text += f"  • Identifies burial directionality and pocket openings\n"
        
        # Remove trailing newlines
        params_text = params_text.rstrip()
        
        panel = Panel(params_text, title="Analysis Parameters", border_style="green", expand=False)
        console.print(panel)
        
        # Return whether there are adjustable parameters
        return has_parameters

    def _modify_analysis_parameters(self, analyzer, selected_metrics, console):
        """Allow user to modify parameters for selected analyses only."""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context, IntPrompt, FloatPrompt

        modified = False

        # Build parameter choices based on selected metrics
        param_choices = {}
        choice_num = 1
        
        if '1' in selected_metrics:  # Metal coordination
            param_choices[str(choice_num)] = ('metal_distance_cutoff', 'Metal coordination distance', 'Å')
            choice_num += 1
            
        if '2' in selected_metrics:  # Hydrogen bonds
            param_choices[str(choice_num)] = ('hbond_distance_cutoff', 'H-bond distance cutoff', 'Å')
            choice_num += 1
            param_choices[str(choice_num)] = ('hbond_atoms', 'H-bond atoms (comma-separated)', '')
            choice_num += 1
            param_choices[str(choice_num)] = ('max_hbonds_per_water', 'Max H-bonds per water', '')
            choice_num += 1
        
        if '4' in selected_metrics:  # Burial (solvent accessibility)
            param_choices[str(choice_num)] = ('sasa_probe_radius', 'Probe radius', 'Å')
            choice_num += 1
            param_choices[str(choice_num)] = ('burial_atom_types', 'Occluding atom types', '')
            choice_num += 1

        if '6' in selected_metrics or '7' in selected_metrics:  # Count-based profiles
            param_choices[str(choice_num)] = ('burial_radius', 'Burial radius (3.0-8.0)', 'Å')
            choice_num += 1
            if '4' not in selected_metrics:
                param_choices[str(choice_num)] = ('burial_atom_types', 'Atom types to count', '')
                choice_num += 1
            param_choices[str(choice_num)] = ('burial_weighting', 'Weighting scheme (count/distance/vdw)', '')
            choice_num += 1
            
        if '5' in selected_metrics:  # Interface
            param_choices[str(choice_num)] = ('interface_distance_cutoff', 'Interface proximity cutoff', 'Å')
            choice_num += 1
            
        if '8' in selected_metrics:  # Network analysis
            param_choices[str(choice_num)] = ('hbond_distance_cutoff', 'H-bond distance cutoff', 'Å')
            choice_num += 1
            param_choices[str(choice_num)] = ('hbond_atoms', 'H-bond atoms (comma-separated)', '')
            choice_num += 1
            param_choices[str(choice_num)] = ('network_min_cluster_size', 'Minimum cluster size to report', '')
            choice_num += 1
            param_choices[str(choice_num)] = ('network_show_isolated', 'Show isolated waters (True/False)', '')
            choice_num += 1
            param_choices[str(choice_num)] = ('network_max_display_size', 'Max network size for detailed ASCII', '')
            choice_num += 1
            
        param_choices[str(choice_num)] = ('done', 'Done with changes', '')
        
        if len(param_choices) == 1:  # Only "done" option
            console.print("[yellow]\nNo adjustable parameters for selected analyses.\n[/yellow]")
            return
        
        while True:
            console.print("\n[bold]Which parameter would you like to change?[/bold]")
            for key, (param_key, desc, unit) in param_choices.items():
                if param_key == 'done':
                    console.print(f"  {key}. {desc}")
                else:
                    current_val = analyzer.parameters[param_key]
                    console.print(f"  {key}. {desc} ({current_val} {unit})")
                    
            param_options_map = {k: v[1] for k, v in param_choices.items()}
            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=list(param_choices.keys()),
                default=str(choice_num),
                module="PDB Filter - Water Analysis",
                description="Select parameter to modify",
                options_map=param_options_map,
            )
            
            if param_choices[choice][0] == 'done':
                break
                
            param_key, desc, unit = param_choices[choice]
            current_val = analyzer.parameters[param_key]
            
            try:
                if param_key == 'hbond_atoms':
                    new_val_str = prompt_with_context(
                        self.processor,
                        f"Enter new {desc.lower()} (e.g., N,O,S)",
                        default=str(current_val),
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc}",
                    )
                    analyzer.set_parameters(**{param_key: new_val_str})
                    console.print(f"[green]Updated: {desc} = {new_val_str}[/green]")
                    modified = True
                elif param_key == 'burial_atom_types':
                    console.print("\n[cyan]Atom types to include in burial calculation:[/cyan]")
                    console.print("  • protein - Standard amino acids")
                    console.print("  • hetero - Non-standard residues (ligands, cofactors)")
                    console.print("  • water - Other water molecules")
                    console.print("  • metal - Metal ions")
                    console.print("  Classes: protein, hetero, metal, water. Example: 'protein,hetero,metal' (default) "
                                  "or 'protein,hetero,metal,water' to let other waters shield too\n")
                    new_val_str = prompt_with_context(
                        self.processor,
                        f"Enter atom types (comma-separated)",
                        default=str(current_val),
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc}",
                    )
                    analyzer.set_parameters(**{param_key: new_val_str})
                    console.print(f"[green]Updated: {desc} = {new_val_str}[/green]")
                    modified = True
                elif param_key == 'burial_weighting':
                    console.print("\n[cyan]Weighting schemes:[/cyan]")
                    console.print("  • count - Each atom = 1 point (simple count)")
                    console.print("  • distance - Weight = 1/distance (closer atoms contribute more)")
                    console.print("  • vdw - Weight = (r_water + r_atom)/distance (vdW radius-based weighting)\n")
                    new_val_str = prompt_with_context(
                        self.processor,
                        f"Enter weighting scheme",
                        default=str(current_val),
                        choices=['count', 'distance', 'vdw'],
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc}",
                        options_map={
                            "count": "Each atom = 1 point",
                            "distance": "Weight by 1/distance",
                            "vdw": "Weight by vdW radius / distance",
                        },
                    )
                    analyzer.set_parameters(**{param_key: new_val_str})
                    console.print(f"[green]Updated: {desc} = {new_val_str}[/green]")
                    modified = True
                elif param_key == 'network_show_isolated':
                    new_val = confirm_with_context(
                        self.processor,
                        f"Show isolated waters?",
                        default=current_val,
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc}",
                    )
                    analyzer.set_parameters(**{param_key: new_val})
                    console.print(f"[green]Updated: {desc} = {new_val}[/green]")
                    modified = True
                elif param_key in ['max_hbonds_per_water', 'network_min_cluster_size', 'network_max_display_size']:
                    new_val_str = prompt_with_context(
                        self.processor,
                        f"Enter new {desc.lower()}",
                        default=str(current_val),
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc}",
                    )
                    new_val = int(new_val_str)
                    analyzer.set_parameters(**{param_key: new_val})
                    console.print(f"[green]Updated: {desc} = {new_val}[/green]")
                    modified = True
                elif param_key == 'burial_radius':
                    new_val_str = prompt_with_context(
                        self.processor,
                        f"Enter new {desc.lower()} [{unit}] (3.0-8.0)",
                        default=str(current_val),
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc} (3.0-8.0 {unit})",
                    )
                    new_val = max(3.0, min(8.0, float(new_val_str)))
                    analyzer.set_parameters(**{param_key: new_val})
                    console.print(f"[green]Updated: {desc} = {new_val} {unit}[/green]")
                    modified = True
                else:
                    new_val_str = prompt_with_context(
                        self.processor,
                        f"Enter new {desc.lower()} [{unit}]",
                        default=str(current_val),
                        module="PDB Filter - Water Analysis",
                        description=f"Update {desc}",
                    )
                    new_val = float(new_val_str)
                    analyzer.set_parameters(**{param_key: new_val})
                    console.print(f"[green]Updated: {desc} = {new_val} {unit}[/green]")
                    modified = True
            except Exception as e:
                console.print(f"[red]Invalid input: {e}[/red]")
                
        return modified 

    def _water_halo_groups(self, analyses: List[Dict], selected_metrics: List[str], analyzer):
        """Group waters for viewer halos: ``[(key, display_label, color, [analysis, ...])]``.

        One group per fact a *displayed* metric establishes, each labelled with
        its rule and cutoff. A water can belong to several groups (a
        metal-bound water that is also enclosed appears in both); there is no
        precedence and no combined category, in keeping with the rest of the
        analysis, which reports each measurement on its own. Metrics 6 to 8 are
        display-only and draw nothing.
        """
        p = analyzer.parameters
        spec = []
        if '1' in selected_metrics:
            spec.append(('water_metal', f"Waters: metal within {p['metal_distance_cutoff']} Å", '#e31a1c',
                         lambda a: bool(a.get('coordinating_metal'))))
        if '2' in selected_metrics:
            # One group per partner count: the count is the datum, not a threshold.
            palette = {1: '#cab2d6', 2: '#9e7bb5', 3: '#6a3d9a', 4: '#3f1f66'}
            for n in (1, 2, 3, 4):
                spec.append((f'water_hbond_{n}',
                             f"Waters: {n} H-bond partner{'s' if n > 1 else ''} (≤ {p['hbond_distance_cutoff']} Å)",
                             palette[n], (lambda n: lambda a: a.get('total', 0) == n)(n)))
        if '3' in selected_metrics:
            median_b = analyzer.protein_median_bfactor()
            label = (f"Waters: B below protein median ({median_b:.0f} Å²)" if median_b
                     else "Waters: B below protein median")
            spec.append(('water_b_below_median', label, '#1f78b4',
                         lambda a: a.get('b_factor_ratio') is not None and a['b_factor_ratio'] < 1.0))
        if '4' in selected_metrics:
            touch = analyzer.WATER_OXYGEN_RADIUS + float(p['sasa_probe_radius'])
            spec += [
                ('water_burial_clash', f"Waters: clash, atom < {analyzer.CLASH_DISTANCE} Å", '#d6008f',
                 lambda a: a.get('burial_category') == 'Clash'),
                ('water_burial_enclosed', f"Waters: enclosed (no bulk path within {touch:.2f} Å)", '#b30000',
                 lambda a: a.get('burial_category') == 'Enclosed'),
                ('water_burial_zero', "Waters: SASA 0 Å², bulk-connected", '#ff7f00',
                 lambda a: a.get('burial_category') == 'Buried'),
                ('water_burial_covered90', "Waters: covered ≥ 90%", '#33a02c',
                 lambda a: a.get('burial_category') == 'Exposed' and a.get('burial_covered_pct', 0) >= 90),
                ('water_burial_covered50', "Waters: covered 50–90%", '#a6d96a',
                 lambda a: a.get('burial_category') == 'Exposed' and 50 <= a.get('burial_covered_pct', 0) < 90),
            ]
        if '5' in selected_metrics:
            spec.append(('water_interface', f"Waters: within {p['interface_distance_cutoff']} Å of two chains", '#ff7f00',
                         lambda a: bool(a.get('at_interface'))))
        groups: List[Tuple[str, str, str, list]] = []
        for key, display, color, pred in spec:
            items = [a for a in analyses if pred(a)]
            if items:
                groups.append((key, display, color, items))
        return groups

    def _display_water_analysis_table(self, analyses: List[Dict], selected_metrics: List[str], console, analyzer):
        """Display metric-specific analysis results with separate tables."""
        if not analyses:
            console.print("[yellow]No water analysis results to display.[/yellow]")
            return
            
        # Display separate table for each selected metric
        for metric in selected_metrics:
            if metric == '1':  # Metal distance
                self._display_metal_distance_table(analyses, console, analyzer)
            elif metric == '2':  # Hydrogen bonds
                self._display_hydrogen_bond_table(analyses, console)
            elif metric == '3':  # B-factor
                self._display_bfactor_table(analyses, console)
            elif metric == '4':  # SASA
                self._display_sasa_table(analyses, console, analyzer)
            elif metric == '5':  # Interface
                self._display_interface_table(analyses, console, analyzer)

    def _display_metal_distance_table(self, analyses: List[Dict], console, analyzer):
        """Display metal coordination distance analysis with dynamic filtering."""
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context
        
        # Get current cutoff
        current_cutoff = analyzer.parameters['metal_distance_cutoff']
        
        while True:
            # Find waters with metals detected (using current cutoff)
            metal_analyses = []
            for analysis in analyses:
                metal_dist = analysis.get('metal_distance')
                if metal_dist is not None and metal_dist <= current_cutoff:
                    metal_analyses.append(analysis)
            
            # Sort by distance (smallest first)
            metal_analyses.sort(key=lambda x: x.get('metal_distance', float('inf')))
            
            # Display filter options
            console.print(f"\n[bold]Metal Coordination Analysis - Current Filter: ≤{current_cutoff}Å[/bold]")
            
            if not metal_analyses:
                console.print(f"[yellow]No waters found within {current_cutoff}Å of metals[/yellow]")
                console.print("\nFilter Options:")
                console.print("1. Try larger distance cutoff (10.0Å)", highlight=False)
                console.print("2. Skip metal analysis", highlight=False)
                
                choice = prompt_with_context(
                    self.processor,
                    "Enter choice",
                    choices=['1', '2'],
                    default='1',
                    module="PDB Filter - Metal Analysis",
                    description="Metal-analysis fallback when no waters found",
                    options_map={
                        "1": "Try larger distance cutoff (10.0Å)",
                        "2": "Skip metal analysis",
                    },
                )
                if choice == '1':
                    current_cutoff = 10.0
                    analyzer.set_parameters(metal_distance_cutoff=current_cutoff)
                    continue
                return
            
            # Count waters in different ranges for filter options
            coordinating_count = len([a for a in metal_analyses if a.get('metal_distance', 999) <= 2.5])
            close_count = len([a for a in metal_analyses if a.get('metal_distance', 999) <= 5.0])
            
            console.print("\nFilter Options:")
            console.print(f"1. Show coordinating waters (≤2.5Å) ({coordinating_count} waters)")
            console.print(f"2. Show close proximity waters (≤5.0Å) ({close_count} waters)")
            console.print(f"3. Show current filter (≤{current_cutoff}Å) ({len(metal_analyses)} waters)")
            console.print("4. Custom distance cutoff", highlight=False)
            console.print("5. Proceed with current display", highlight=False)
            
            filter_choice = prompt_with_context(
                self.processor,
                "Select filter",
                choices=['1', '2', '3', '4', '5'],
                default='5',
                module="PDB Filter - Metal Analysis",
                description="Select metal coordination filter",
                options_map={
                    "1": "Show coordinating waters (≤2.5Å)",
                    "2": "Show close proximity waters (≤5.0Å)",
                    "3": f"Show current filter (≤{current_cutoff}Å)",
                    "4": "Custom distance cutoff",
                    "5": "Proceed with current display",
                },
            )
            
            if filter_choice == '1':
                display_cutoff = 2.5
            elif filter_choice == '2':
                display_cutoff = 5.0
            elif filter_choice == '3':
                display_cutoff = current_cutoff
            elif filter_choice == '4':
                new_cutoff = self._get_distance_cutoff(console, current_cutoff)
                if new_cutoff:
                    current_cutoff = new_cutoff
                    analyzer.set_parameters(metal_distance_cutoff=current_cutoff)
                    continue
                else:
                    display_cutoff = current_cutoff
            else:
                display_cutoff = current_cutoff
            
            # Filter the analyses based on the display cutoff
            filtered_analyses = []
            for analysis in analyses:
                metal_dist = analysis.get('metal_distance')
                if metal_dist is not None and metal_dist <= display_cutoff:
                    filtered_analyses.append(analysis)
            
            # Sort by distance (smallest first)
            filtered_analyses.sort(key=lambda x: x.get('metal_distance', float('inf')))
            
            # Display the table with filtered results
            table = Table(
                title=f"Metal Coordination Analysis ({len(filtered_analyses)} waters ≤{display_cutoff}Å)",
                show_lines=True  
            )
            table.add_column("Residue #", style="bold yellow", width=10)
            table.add_column("Name", style="yellow", width=8)
            table.add_column("Distance (Å)", style="red", width=12)
            table.add_column("Closest Metal", style="red", width=20)
            table.add_column("Coordinating", style="green", width=12)
            
            for analysis in filtered_analyses:
                metal_dist = analysis.get('metal_distance')
                
                # Color coding: red for coordinating distances (≤2.5Å), yellow for close (≤5Å)
                if metal_dist <= 2.5:
                    dist_str = f"[bold red]{metal_dist:.2f}[/bold red]"
                    coord_str = "[bold green]YES[/bold green]"
                elif metal_dist <= 5.0:
                    dist_str = f"[yellow]{metal_dist:.2f}[/yellow]"
                    coord_str = "Close"
                else:
                    dist_str = f"{metal_dist:.2f}"
                    coord_str = "No"
                    
                closest_metal = analysis.get('closest_metal', 'None')
                
                table.add_row(
                    f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                    analysis['residue_name'],
                    dist_str,
                    closest_metal[:20],
                    coord_str
                )
                
            console.print(table)
            
            # Add legend
            coord_cutoff = analyzer.parameters['metal_distance_cutoff']
            legend_text = "\n[bold]Legend:[/bold]\n"
            legend_text += f"• [bold red]Red[/bold red] distances: ≤[cyan]{coord_cutoff}Å[/cyan] (coordinating distance)\n"
            legend_text += f"• [yellow]Yellow[/yellow] distances: [cyan]{coord_cutoff}-5.0Å[/cyan] (close proximity)\n"
            legend_text += f"• Showing waters within [cyan]{display_cutoff}Å[/cyan] of nearest metal\n"
            console.print(legend_text)
            
            break
        
    def _get_distance_cutoff(self, console, current_cutoff):
        """Get new distance cutoff from user."""
        from proprep.utils.prompts import prompt_with_context
        
        try:
            new_cutoff_str = prompt_with_context(
                self.processor,
                f"Enter new distance cutoff [Å]",
                default=str(current_cutoff),
                module="PDB Filter - Metal Analysis",
                description="New metal-distance cutoff in Å",
            )
            new_cutoff = float(new_cutoff_str)
            if new_cutoff <= 0:
                console.print("[red]Distance must be positive[/red]")
                return None
            return new_cutoff
        except ValueError:
            console.print("[red]Invalid distance value[/red]")
            return None
        
    def _display_hydrogen_bond_table(self, analyses: List[Dict], console):
        """Display hydrogen bond analysis with dynamic filtering."""
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context
        
        # Sort by total H-bonds (highest first) 
        hbond_analyses = sorted(analyses, key=lambda x: x.get('total', 0), reverse=True)
        
        while True:
            # Calculate filter statistics
            all_count = len(hbond_analyses)
            highly_connected = [a for a in hbond_analyses if a.get('total', 0) >= 4]
            well_connected = [a for a in hbond_analyses if a.get('total', 0) >= 2]
            with_protein = [a for a in hbond_analyses if a.get('protein', 0) > 0]
            with_water = [a for a in hbond_analyses if a.get('water', 0) > 0]
            with_hetero = [a for a in hbond_analyses if a.get('hetero', 0) > 0]
            
            console.print(f"\n[bold]Hydrogen Bond Analysis - Display Filter Options[/bold]")
            console.print(f"1. Show all waters ({all_count} waters)")
            console.print(f"2. Show highly connected (4 H-bonds) ({len(highly_connected)} waters)")
            console.print(f"3. Show well connected (≥2 H-bonds) ({len(well_connected)} waters)")
            console.print(f"4. Show waters with protein H-bonds ({len(with_protein)} waters)")
            console.print(f"5. Show waters with water-water H-bonds ({len(with_water)} waters)")
            console.print(f"6. Show waters with hetero H-bonds ({len(with_hetero)} waters)")
            console.print("7. Custom minimum H-bond count", highlight=False)
            console.print("8. Proceed with water selection", highlight=False)
            
            filter_choice = prompt_with_context(
                self.processor,
                "Select display filter",
                choices=['1', '2', '3', '4', '5', '6', '7', '8'],
                default='8',
                module="PDB Filter - H-Bond Analysis",
                description="Select H-bond display filter",
                options_map={
                    "1": "Show all waters",
                    "2": "Show highly connected (4 H-bonds)",
                    "3": "Show well connected (≥2 H-bonds)",
                    "4": "Show waters with protein H-bonds",
                    "5": "Show waters with water-water H-bonds",
                    "6": "Show waters with hetero H-bonds",
                    "7": "Custom minimum H-bond count",
                    "8": "Proceed with water selection",
                },
            )
            
            if filter_choice == '1':
                filtered_analyses = hbond_analyses
                filter_desc = "All waters"
            elif filter_choice == '2':
                filtered_analyses = highly_connected
                filter_desc = "Highly connected (≥4 H-bonds)"
            elif filter_choice == '3':
                filtered_analyses = well_connected
                filter_desc = "Well connected (≥2 H-bonds)"
            elif filter_choice == '4':
                filtered_analyses = with_protein
                filter_desc = "With protein H-bonds"
            elif filter_choice == '5':
                filtered_analyses = with_water
                filter_desc = "With water-water H-bonds"
            elif filter_choice == '6':
                filtered_analyses = with_hetero
                filter_desc = "With hetero H-bonds"
            elif filter_choice == '7':
                min_hbonds = self._get_minimum_hbonds(console)
                if min_hbonds is not None:
                    filtered_analyses = [a for a in hbond_analyses if a.get('total', 0) >= min_hbonds]
                    filter_desc = f"≥{min_hbonds} H-bonds"
                else:
                    continue
            else:
                break
            
            # Display filtered table
            self._display_hbond_table_data(filtered_analyses, filter_desc, console)
            
            # Ask if user wants to change filter
            console.print("\nDisplay Options:")
            console.print("1. Change display filter", highlight=False)
            console.print("2. Return to the Analysis Options menu", highlight=False)
            
            proceed_choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=['1', '2'],
                default='2',
                module="PDB Filter - H-Bond Analysis",
                description="Continue filtering or return to analysis menu",
                options_map={"1": "Change display filter", "2": "Return to Analysis Options menu"},
            )
            if proceed_choice == '2':
                break

    def _get_minimum_hbonds(self, console):
        """Get minimum H-bond count from user."""
        from proprep.utils.prompts import prompt_with_context
        
        try:
            min_hbonds_str = prompt_with_context(
                self.processor,
                "Enter minimum H-bond count",
                default="2",
                module="PDB Filter - H-Bond Analysis",
                description="Minimum H-bond count for filter",
            )
            min_hbonds = int(min_hbonds_str)
            if min_hbonds < 0:
                console.print("[red]H-bond count must be non-negative[/red]")
                return None
            return min_hbonds
        except ValueError:
            console.print("[red]Invalid H-bond count[/red]")
            return None
        
    def _display_hbond_table_data(self, filtered_analyses, filter_desc, console):
        """Display the actual H-bond table data."""
        from rich.table import Table
        
        table = Table(
            title=f"Hydrogen Bond Analysis ({len(filtered_analyses)} waters) - {filter_desc}",
            show_lines=True
        )
        table.add_column("Residue #", style="bold yellow", width=10)
        table.add_column("Name", style="yellow", width=8)
        table.add_column("Total H-bonds", style="green", width=12)
        table.add_column("Protein", style="blue", width=8)
        table.add_column("Water", style="cyan", width=8)
        table.add_column("Hetero", style="magenta", width=8)
        table.add_column("H-bond Partners", style="green", width=30)
        table.add_column("Category", style="green", width=15)
        
        for analysis in filtered_analyses:
            total_hbonds = analysis.get('total', 0)
            prot = analysis.get('protein', 0)
            water = analysis.get('water', 0)
            hetero = analysis.get('hetero', 0)
            
            # Get H-bond details
            hbond_details = analysis.get('details', [])
            
            # Format partner information - show ALL partners (up to 4)
            if hbond_details:
                partners = []
                for detail in hbond_details:
                    partner_str = f"{detail['partner_atom']}({detail['partner_residue']})"
                    partners.append(partner_str)
                    
                partners_str = ", ".join(partners)
            else:
                partners_str = "None"
            
            # Color coding for H-bond count
            if total_hbonds >= 4:
                total_str = f"[bold green]{total_hbonds}[/bold green]"
                category = "Highly connected"
            elif total_hbonds >= 2:
                total_str = f"[green]{total_hbonds}[/green]"
                category = "Well connected"
            else:
                total_str = str(total_hbonds)
                category = "Poorly connected"
                
            table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                analysis['residue_name'],
                total_str,
                str(prot),
                str(water),
                str(hetero),
                partners_str,
                category
            )
            
        console.print(table)
        
        # Add legend
        legend_text = "\n[bold]Legend:[/bold]\n"
        legend_text += "• [bold green]Bold green[/bold green]: ≥4 H-bonds (highly connected)\n"
        legend_text += "• [green]Green[/green]: 2-3 H-bonds (well connected)\n"
        legend_text += "• Partners shown as: AtomName(ResidueName+Number)\n"
        legend_text += f"• Partners are N/O/S atoms within the heavy-atom cutoff, closest first; no angular test\n"
        console.print(legend_text)

    def _display_bfactor_table(self, analyses: List[Dict], console):
        """Display B-factor analysis."""
        from rich.table import Table
        
        # Sort by B-factor (lowest first)
        bfactor_analyses = sorted(analyses, key=lambda x: x.get('b_factor', float('inf')))
        
        median_b = next((a.get('protein_median_b') for a in analyses if a.get('protein_median_b')), None)
        title = f"B-factor Analysis ({len(analyses)} waters"
        title += f"; protein heavy-atom median B = {median_b:.1f} Å²)" if median_b else ")"
        table = Table(title=title)
        table.add_column("Residue #", style="bold yellow", width=10)
        table.add_column("Name", style="yellow", width=8)
        table.add_column("B (Å²)", style="blue", width=10, justify="right")
        table.add_column("× protein median", style="blue", width=17, justify="right")
        table.add_column("Category", style="blue", width=22)

        for analysis in bfactor_analyses:
            bfactor = analysis.get('b_factor', 0)
            ratio = (bfactor / median_b) if median_b else None
            if ratio is None:
                bfactor_str, ratio_str, category = f"{bfactor:.1f}", "-", "no protein reference"
            elif ratio < 1.0:
                bfactor_str, ratio_str, category = f"[blue]{bfactor:.1f}[/blue]", f"[blue]{ratio:.2f}[/blue]", "below protein median"
            else:
                bfactor_str, ratio_str, category = f"{bfactor:.1f}", f"{ratio:.2f}", "above protein median"
            table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                analysis['residue_name'],
                bfactor_str,
                ratio_str,
                category,
            )

        console.print(table)

        legend_text = "\n[bold]Legend:[/bold] B-factors depend on resolution and refinement, so waters are compared with "
        legend_text += "this structure's own protein heavy atoms rather than an absolute cutoff. "
        legend_text += "[blue]Blue[/blue] = B below the protein median (at least as ordered as a typical protein atom); "
        legend_text += "the ratio column is the number, sort order is ascending B."
        console.print(legend_text)

    def _display_sasa_table(self, analyses: List[Dict], console, analyzer):
        """Display water burial: accessible area of the oxygen and bulk connectivity."""
        from rich.table import Table
        from rich.panel import Panel

        probe = analyzer.parameters['sasa_probe_radius']
        isolated = analyses[0].get('burial_sasa_isolated', 0.0) if analyses else 0.0

        touch = analyzer.WATER_OXYGEN_RADIUS + probe
        grid_h = analyzer.ENCLOSURE_GRID_SPACING
        reach = touch + grid_h
        info_text = "[bold cyan]About Water Burial[/bold cyan]\n\n"
        info_text += ("Two geometric questions are asked about each water oxygen, at two different distances. "
                      f"The probe is a {probe} Å sphere, i.e. another water. A probe [bold]touches[/bold] the oxygen when its "
                      f"centre is exactly r_water + r_probe = {analyzer.WATER_OXYGEN_RADIUS} + {probe} = {touch:.2f} Å away.\n\n")
        info_text += f"[bold]1. Can anything touch it? (SASA, measured on the {touch:.2f} Å contact sphere)[/bold]\n"
        info_text += (f"The accessible area is the part of that {touch:.2f} Å sphere where a probe centre can sit without "
                      f"overlapping any {analyzer.parameters['burial_atom_types']} atom (Lee-Richards, Bondi radii; other "
                      f"waters never occlude). A free water has the whole sphere, {isolated:.0f} Å². 0 Å² means no point on "
                      "the contact sphere is free: nothing can be in contact with this water.\n\n")
        info_text += f"[bold]2. Does the free space around it reach the outside? (Access, looked for within {reach:.2f} Å)[/bold]\n"
        info_text += (f"Every position where a probe centre fits is marked on a {grid_h} Å grid and the marked positions "
                      "are joined into connected regions; the region touching the box edge is bulk solvent. The test then "
                      f"asks whether any bulk-connected position lies within {touch:.2f} + {grid_h} = {reach:.2f} Å of the oxygen "
                      "(contact distance plus one grid cell of tolerance). If none does, the water is [bold]enclosed[/bold]: "
                      "it cannot leave, or be replaced, without the protein moving. Note this looks slightly beyond the "
                      "contact sphere, which is why a water can have 0 Å² and still be bulk-connected: a probe can come down "
                      "an open cleft to just above it but not the last fraction of an Ångström to touch it.\n\n")
        info_text += "[bold]Covered %:[/bold]\n"
        info_text += (f"100 × (1 − SASA / {isolated:.0f} Å²): the share of the water's accessible surface that the "
                      "structure takes away. 0% is a free water, 100% is untouchable.\n\n")
        info_text += "[bold]Category rules (applied in this order; the first that matches wins):[/bold]\n"
        info_text += (f"1. [bold magenta]Clash[/bold magenta]    nearest C/N/O/S/P atom < {analyzer.CLASH_DISTANCE} Å "
                      "(wwPDB close-contact criterion; metals exempt, coordination is 2.0-2.2 Å)\n")
        info_text += (f"2. [bold red]Enclosed[/bold red] no bulk-connected probe position within {reach:.2f} Å of the oxygen: "
                      "the free space around it, if any, is a sealed pocket\n")
        info_text += ("3. [bold dark_orange3]Buried[/bold dark_orange3]   SASA = 0.0 Å² (nothing can touch it) but open space "
                      f"within {reach:.2f} Å reaches bulk: a water at the bottom of a cleft too narrow to enter. "
                      f"A thin class by construction (the {touch:.2f}-{reach:.2f} Å band); often empty\n")
        info_text += "4. [cyan]Exposed[/cyan]  everything else; read the SASA and Covered columns\n\n"
        info_text += ("The rules are yes/no statements, not bins on the continuum: a 99%-covered water is still "
                      "Exposed. The numbers are shown so you can draw your own line.")

        panel = Panel(info_text, title="[bold]Method Information[/bold]", border_style="cyan", width=min(100, console.width))
        console.print()
        console.print(panel)
        console.print()

        order = {'Clash': 0, 'Enclosed': 1, 'Buried': 2, 'Exposed': 3}
        burial_analyses = sorted(
            analyses, key=lambda x: (order.get(x.get('burial_category'), 3), x.get('burial_sasa', 0.0)))

        counts = {k: sum(1 for a in analyses if a.get('burial_category') == k) for k in order}
        table = Table(title=f"Water Burial ({len(analyses)} waters: {counts['Enclosed']} enclosed, "
                            f"{counts['Buried']} buried, {counts['Clash']} clashing)")
        table.add_column("Residue #", style="bold yellow", width=10)
        table.add_column("Name", style="yellow", width=8)
        table.add_column("SASA (Å²)", style="cyan", width=10, justify="right")
        table.add_column("Covered", style="cyan", width=8, justify="right")
        table.add_column("Access", style="blue", width=10)
        table.add_column("Nearest atom", style="grey50", width=18)
        table.add_column("Category", style="dark_orange3", width=12)

        style = {'Clash': 'bold magenta', 'Enclosed': 'bold red', 'Buried': 'bold dark_orange3', 'Exposed': 'cyan'}
        for analysis in burial_analyses:
            category = analysis.get('burial_category', 'Exposed')
            sasa = analysis.get('burial_sasa', 0.0)
            covered = analysis.get('burial_covered_pct')
            if covered is None:
                covered = 100.0 * (1.0 - sasa / isolated) if isolated else 0.0
            d = analysis.get('burial_closest_distance')
            nearest = f"{analysis.get('burial_closest_atom', '')} {d:.2f} Å" if d is not None else "-"
            table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                analysis['residue_name'],
                f"[{style[category]}]{sasa:.1f}[/{style[category]}]",
                f"{covered:.0f}%",
                analysis.get('burial_access', '-'),
                nearest,
                f"[{style[category]}]{category}[/{style[category]}]",
            )

        console.print(table)

        legend = "\n[bold]Legend:[/bold] "
        legend += f"SASA = area a {probe} Å probe can touch (isolated water {isolated:.0f} Å²); "
        legend += f"Covered = 100 × (1 − SASA/{isolated:.0f}); "
        legend += "Access = bulk / enclosed from the flood fill.\n"
        legend += (f"Category: [bold magenta]Clash[/bold magenta] nearest C/N/O/S/P < {analyzer.CLASH_DISTANCE} Å  →  "
                   f"[bold red]Enclosed[/bold red] no bulk-connected probe position within {reach:.2f} Å  →  "
                   f"[bold dark_orange3]Buried[/bold dark_orange3] SASA = 0.0 on the {touch:.2f} Å contact sphere, "
                   f"but bulk reachable within {reach:.2f} Å  →  [cyan]Exposed[/cyan] otherwise.")
        console.print(legend)

    def _display_interface_table(self, analyses: List[Dict], console, analyzer):
        """Display interface proximity analysis."""
        from rich.table import Table
        from rich.panel import Panel

        # Display information panel
        info_text = "[bold cyan]About Interface Proximity Analysis[/bold cyan]\n\n"
        info_text += "[bold]Purpose:[/bold]\n"
        info_text += "Identifies waters that bridge between different protein chains, which often play\n"
        info_text += "important structural or functional roles at protein-protein interfaces.\n\n"

        info_text += "[bold]Method:[/bold]\n"
        info_text += f"• Cutoff distance: {analyzer.parameters['interface_distance_cutoff']}Å\n"
        info_text += "• A water is considered at an interface if it is simultaneously close to\n"
        info_text += "  atoms from two or more different chains\n\n"

        info_text += "[bold]Interpretation:[/bold]\n"
        info_text += "• Interface waters often stabilize protein-protein interactions\n"
        info_text += "• These waters are typically structurally conserved and functionally important\n"
        info_text += "• May mediate inter-chain hydrogen bonding networks"

        panel = Panel(info_text, title="[bold]Method Information[/bold]", border_style="cyan", expand=False)
        console.print(panel)
        console.print()

        # Filter and sort interface waters
        interface_analyses = [a for a in analyses if a.get('at_interface', False)]
        interface_analyses.sort(key=lambda x: x.get('residue_number', 0))

        if len(interface_analyses) == 0:
            console.print("[yellow]No waters detected at protein-protein interfaces.[/yellow]")
            console.print("[grey50]Note: This analysis requires structures with multiple protein chains.[/grey50]")
            return

        table = Table(title=f"Interface Proximity Analysis ({len(interface_analyses)} waters at interfaces)")
        table.add_column("Residue #", style="bold yellow", width=10)
        table.add_column("Name", style="yellow", width=8)
        table.add_column("Bridged Chains", style="purple", width=15)
        table.add_column("Category", style="purple", width=20)

        for analysis in interface_analyses:
            bridged_chains = analysis.get('bridged_chains', [])
            chains_str = ", ".join(bridged_chains) if bridged_chains else "N/A"

            table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                analysis['residue_name'],
                f"[bold purple]{chains_str}[/bold purple]",
                "Interface water"
            )

        console.print(table)

        # Add legend
        legend_text = "\n[bold]Legend:[/bold]\n"
        legend_text += "• [bold purple]Bridged Chains[/bold purple]: Chain IDs that the water molecule is bridging\n"
        legend_text += f"• Interface waters are within {analyzer.parameters['interface_distance_cutoff']}Å of atoms from multiple chains"
        console.print(legend_text)

    def _display_burial_profile_table(self, analyses: List[Dict], console):
        """Display burial profile analysis with ASCII charts for selected waters."""
        from rich.table import Table
        from rich.panel import Panel
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        
        if not analyses:
            console.print("[yellow]No water analyses available for profiling.[/yellow]")
            return
        
        # Filter analyses with profile data
        profile_analyses = []
        for analysis in analyses:
            if 'burial_profile' in analysis and analysis['burial_profile']:
                profile_analyses.append(analysis)
        
        if not profile_analyses:
            console.print("[yellow]No burial profile data available.[/yellow]")
            return
        
        profile_analyses.sort(key=lambda x: x['burial_profile'].get('final_count', 0), reverse=True)
        
        # Summary table
        table = Table(title=f"Multi-Radius Burial Profiles ({len(profile_analyses)} waters)")
        table.add_column("Residue #", style="bold yellow", width=10)
        table.add_column("Final Count", style="dark_orange3", width=12)
        table.add_column("Saturation (<10%/step)", style="blue", width=22)
        table.add_column("Steep Rise", style="green", width=15)
        
        for analysis in profile_analyses:
            profile = analysis['burial_profile']
            
            final_count = f"{profile.get('final_count', 0):.1f}"
            
            saturation = "N/A"
            if profile.get('saturation_radius'):
                saturation = f"{profile['saturation_radius']:.1f}Å"
            
            steep_rise = "N/A"
            if profile.get('steepest_start') and profile.get('steepest_end'):
                steep_rise = f"{profile['steepest_start']:.1f}-{profile['steepest_end']:.1f}Å"
            
            table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                final_count,
                saturation,
                steep_rise
            )
        
        console.print(table)
        console.print(
            "\n[bold]Legend:[/bold] counts are weighted atom counts within each radius (2.0 to 8.0 Å in 0.5 Å "
            "steps). [bold]Saturation[/bold] = first radius at which the count grows by less than 10% over the "
            "previous step; [bold]Steep Rise[/bold] = the 0.5 Å step with the largest increase. The 10% is a "
            "convention for labelling the curve, not a property of the water; the profile itself is the result."
        )
        
        # Offer to show detailed profiles
        if confirm_with_context(
            self.processor,
            "\nShow detailed burial profiles for specific waters?",
            default=True,
            module="PDB Filter - Burial Profile",
            description="Show detailed burial profiles for specific waters",
        ):
            while True:
                available_nums = [str(a['residue_number']) for a in profile_analyses]
                
                console.print(f"\n[grey50]Available waters: {', '.join(available_nums)}[/grey50]")
                choice = prompt_with_context(
                    self.processor,
                    "Enter residue number for detailed profile (or 'done' to continue)",
                    default="done",
                    module="PDB Filter - Burial Profile",
                    description="Residue number for burial profile view",
                )
                
                if choice.lower() == 'done':
                    break
                
                try:
                    res_num = int(choice)
                    selected_analysis = None
                    for analysis in profile_analyses:
                        if analysis['residue_number'] == res_num:
                            selected_analysis = analysis
                            break
                    
                    if selected_analysis and 'burial_profile_ascii' in selected_analysis:
                        title = f"Burial Profile for HOH{res_num}"
                        ascii_chart = selected_analysis['burial_profile_ascii']
                        panel = Panel(ascii_chart, title=title, border_style="blue")
                        console.print(panel)
                    else:
                        console.print(f"[yellow]No ASCII profile available for HOH{res_num}[/yellow]")
                        
                except ValueError:
                    console.print(f"[red]Invalid residue number: {choice}[/red]")

    def _display_directional_analysis_table(self, analyses: List[Dict], console):
        """Display directional burial analysis with compass charts."""
        from rich.table import Table
        from rich.panel import Panel
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        
        if not analyses:
            console.print("[yellow]No water analyses available for directional analysis.[/yellow]")
            return
        
        # Filter analyses with directional data
        directional_analyses = []
        for analysis in analyses:
            if 'directional_burial' in analysis and analysis['directional_burial']:
                directional_analyses.append(analysis)
        
        if not directional_analyses:
            console.print("[yellow]No directional burial data available.[/yellow]")
            return
        
        directional_analyses.sort(
            key=lambda x: x['directional_burial'].get('total_weight', 0), 
            reverse=True
        )
        
        # Summary table
        table = Table(title=f"Directional Burial Analysis ({len(directional_analyses)} waters)")
        table.add_column("Residue #", style="bold yellow", width=10)
        table.add_column("Total Weight", style="dark_orange3", width=12)
        table.add_column("Primary Dir", style="blue", width=12)
        table.add_column("Least Buried", style="cyan", width=13)
        table.add_column("Pattern", style="grey50", width=44)

        for analysis in directional_analyses:
            directional = analysis['directional_burial']

            total_weight = f"{directional.get('total_weight', 0):.1f}"
            primary_dir = directional.get('primary_direction', 'N/A')
            pocket_opening = directional.get('pocket_opening', 'N/A')
            pattern = directional.get('pattern_type', 'N/A')

            table.add_row(
                f"[bold yellow]{analysis['residue_number']}[/bold yellow]",
                total_weight,
                primary_dir,
                pocket_opening,
                pattern,
            )
        
        console.print(table)
        console.print(
            "\n[bold]Legend:[/bold] sectors are 45° wedges of azimuth in the xy-plane of the coordinate frame, "
            "so directions are relative to the crystal frame, not the protein. [bold]Primary[/bold] = sector with "
            "the largest weighted count, [bold]Least Buried[/bold] = the smallest. [bold]Pattern[/bold] compares "
            "the sector-count range (max − min) with the mean sector count; the 0.5× and 1.5× in each label are "
            "conventions for naming the shape, not properties of the water. Compass glyphs (●, ●●, ●●●, ████) are "
            "quarters of that water's own largest sector."
        )
        
        # Offer to show detailed compass charts
        if confirm_with_context(
            self.processor,
            "\nShow detailed directional compass for specific waters?",
            default=True,
            module="PDB Filter - Directional Analysis",
            description="Show detailed compass for specific waters",
        ):
            while True:
                available_nums = [str(a['residue_number']) for a in directional_analyses]
                
                console.print(f"\n[grey50]Available waters: {', '.join(available_nums)}[/grey50]")
                choice = prompt_with_context(
                    self.processor,
                    "Enter residue number for compass view (or 'done' to continue)",
                    default="done",
                    module="PDB Filter - Directional Analysis",
                    description="Residue number for compass view",
                )
                
                if choice.lower() == 'done':
                    break
                
                try:
                    res_num = int(choice)
                    selected_analysis = None
                    for analysis in directional_analyses:
                        if analysis['residue_number'] == res_num:
                            selected_analysis = analysis
                            break
                    
                    if selected_analysis and 'directional_compass' in selected_analysis:
                        compass_chart = selected_analysis['directional_compass']
                        if compass_chart and isinstance(compass_chart, str):
                            title = f"Directional Burial Compass for HOH{res_num}"
                            panel = Panel(compass_chart, title=title, border_style="blue")
                            console.print(panel)
                        else:
                            console.print(f"[yellow]Compass chart data is invalid for HOH{res_num}[/yellow]")
                            console.print(f"[grey50]Debug: compass_chart = {repr(compass_chart)}[/grey50]")
                    else:
                        console.print(f"[yellow]No compass chart available for HOH{res_num}[/yellow]")
                        if selected_analysis:
                            available_keys = list(selected_analysis.keys())
                            console.print(f"[grey50]Available data keys: {available_keys}[/grey50]")
                                                    

                except ValueError:
                    console.print(f"[red]Invalid residue number: {choice}[/red]")

    def _display_multiradius_explanation(self, analyzer, console):
        """Display detailed explanation of multi-radius profiling."""
        from rich.panel import Panel
        
        explanation = f"""[bold cyan]Multi-Radius Burial Profiling Explanation[/bold cyan]

    [bold]What it does:[/bold]
    - Analyzes burial at multiple radii from 2.0 to 8.0 Å in 0.5 Å steps
    - Counts {analyzer.parameters['burial_atom_types']} atoms around each water at each radius
    - Creates a "burial profile" showing how burial changes with distance

    [bold]Atom types analyzed:[/bold]
    - protein: Standard amino acids (ALA, VAL, PHE, etc.)
    - hetero: Non-standard residues, ligands, cofactors
    - Water molecules and metals are excluded from the count

    [bold]Key metrics calculated:[/bold]
    - [yellow]Saturation point:[/yellow] Radius where burial levels off (<10% increase)
    - [yellow]Steepest rise:[/yellow] Radius range with maximum burial slope (first shell)
    - [yellow]Final count:[/yellow] Total burial at maximum radius (8.0 Å)

    [bold]Visual output:[/bold]
    - ASCII chart showing burial count vs radius
    - Annotations highlighting saturation and steep rise regions
    - Interactive display of individual water profiles

    This helps identify:
    - Waters in tight binding pockets (early saturation)
    - Waters in protein cavities (gradual saturation)
    - Surface waters (no saturation within 8 Å)"""
        
        panel = Panel(explanation, title="Multi-Radius Profiling Guide", border_style="blue", expand=False)
        console.print(panel)

    def _display_directional_explanation(self, analyzer, console):
        """Display detailed explanation of directional burial analysis."""
        from rich.panel import Panel

        explanation = f"""[bold cyan]Directional Burial Analysis Explanation[/bold cyan]

[bold]What it does:[/bold]
- Divides the space around each water into 8 compass sectors (45° each)
- Counts {analyzer.parameters['burial_atom_types']} atoms in each direction
- Identifies where protein/ligand atoms are concentrated around the water

[bold]How it works:[/bold]
- Uses {analyzer.parameters['burial_radius']} Å radius sphere around water
- Calculates angle from water to each nearby atom (0°=East, 90°=North, etc.)
- Assigns atoms to sectors: N, NE, E, SE, S, SW, W, NW
- Counts atoms and calculates burial weight per sector

[bold yellow]⚠️  Important Note About Directions:[/bold yellow]
- Compass directions (N/S/E/W) are relative to the PDB coordinate frame
- Rotating the protein would change which direction is labeled "North"
- The biological meaning lies in the PATTERN, not the absolute directions
- Compare burial asymmetry and accessibility, not specific compass labels

[bold]Compass interpretation:[/bold]
- ⚬ = Water molecule at center
- ○ = No atoms in that direction (0 count)
- ● / ●● / ●●● / ████ = up to 25% / 50% / 75% / above 75% of this water's largest sector count
- Numbers in () = Actual atom count in that sector

[bold]Pattern label (a naming convention, stated with each label):[/bold]
- Uniform: sector-count range (max − min) ≤ 0.5 × mean sector count
- Moderately directional: range ≤ 1.5 × mean
- Highly directional: range > 1.5 × mean

[bold]Key metrics:[/bold]
- [yellow]Primary direction:[/yellow] Sector with highest atom density
- [yellow]Secondary direction:[/yellow] Sector with second-highest density  
- [yellow]Least buried direction:[/yellow] Sector with lowest atom density (potential access)
- [yellow]Total weight:[/yellow] Sum of all directional burial weights

[bold]Biological meaning:[/bold]
- High burial directions → Protein contact surfaces
- Low burial directions → Solvent-accessible regions or pocket openings
- Asymmetric patterns → Waters in binding sites or cavities
- Uniform patterns → Waters in bulk solvent

[bold]Focus on patterns, not absolute directions![/bold]"""
        
        panel = Panel(explanation, title="Directional Burial Analysis Guide", border_style="blue", expand=False)
        console.print(panel)
        
    def _display_water_network_analysis(self, network_data: Dict, console):
        """Display water network analysis results."""
        from rich.panel import Panel
        from rich.table import Table
        from proprep.utils.prompts import prompt_with_context, confirm_with_context
        
        if not network_data:
            console.print("[yellow]No network data available.[/yellow]")
            return
        
        network_analysis = network_data.get('network_analysis', {})
        
        if 'error' in network_analysis:
            console.print(f"[red]Network Analysis Error: {network_analysis['error']}[/red]")
            return
        
        # Display summary table
        table = Table(title=f"Water Network Summary ({network_analysis.get('total_waters', 0)} waters)")
        table.add_column("Metric", style="cyan", width=25)
        table.add_column("Value", style="yellow", width=15)
        table.add_column("Description", style="grey50", width=40)
        
        table.add_row(
            "Total Waters",
            str(network_analysis.get('total_waters', 0)),
            "All water molecules analyzed"
        )
        table.add_row(
            "H-bond Connections", 
            str(network_analysis.get('total_connections', 0)),
            "Direct water-water H-bonds"
        )
        table.add_row(
            "Connected Clusters",
            str(len(network_analysis.get('significant_components', []))),
            f"Groups of ≥{network_analysis.get('network_analysis', {}).get('network_min_cluster_size', 2)} waters"
        )
        table.add_row(
            "Isolated Waters",
            str(len(network_analysis.get('isolated_waters', []))),
            "Waters with no H-bond connections"
        )
        table.add_row(
            "Largest Cluster",
            str(network_analysis.get('largest_cluster_size', 0)),
            "Size of biggest connected network"
        )
        table.add_row(
            "Hub Waters",
            str(len(network_analysis.get('hub_waters', []))),
            "Waters with ≥3 connections"
        )
        
        console.print(table)
        
        # Show detailed ASCII if requested
        if confirm_with_context(
            self.processor,
            "\nShow detailed network ASCII visualization?",
            default=True,
            module="PDB Filter - Network Analysis",
            description="Show detailed network ASCII visualization",
        ):
            network_ascii = network_data.get('network_ascii', 'No ASCII available')
            panel = Panel(network_ascii, title="Water Network Topology", border_style="blue", expand=False)
            console.print(panel)
        
        # Offer hub water details
        hub_waters = network_analysis.get('hub_waters', [])
        if hub_waters and confirm_with_context(
            self.processor,
            "\nShow hub water details?",
            default=False,
            module="PDB Filter - Network Analysis",
            description="Show hub water details",
        ):
            hub_table = Table(title="Hub Waters (Highly Connected)")
            hub_table.add_column("Water ID", style="bold yellow", width=12)
            hub_table.add_column("Connections", style="green", width=12)
            hub_table.add_column("Chain", style="cyan", width=8)
            
            for water_id, degree in hub_waters[:10]:  # Show top 10
                # Try to get chain info from the analysis
                chain_id = "?"
                for analysis in network_data.get('individual_analyses', []):
                    if analysis.get('residue_number') == water_id:
                        chain_id = analysis.get('chain_id', '?')
                        break
                
                hub_table.add_row(f"HOH{water_id}", str(degree), chain_id)
            
            console.print(hub_table)

    def _display_network_explanation(self, analyzer, console):
        """Display detailed explanation of water network analysis."""
        from rich.panel import Panel
        
        explanation = f"""[bold cyan]Water Network Analysis Explanation[/bold cyan]

    [bold]What it does:[/bold]
    • Builds a network where waters are nodes and H-bonds are connections
    • Identifies clusters of interconnected waters
    • Finds hub waters that coordinate multiple water molecules
    • Analyzes network topology (chains, cycles, isolated waters)

    [bold]How it works:[/bold]
    • Uses existing H-bond analysis to find water-water connections
    • Applies graph theory algorithms to analyze connectivity patterns
    • Creates ASCII visualizations of network structures
    • Identifies biologically relevant water organization

    [bold]Network types (Phase 1):[/bold]
    • [yellow]water_only:[/yellow] Direct water-water H-bonds only
    • Future phases will include water-protein-water bridges

    [bold]Key metrics:[/bold]
    • [yellow]Connected clusters:[/yellow] Groups of H-bonded waters
    • [yellow]Hub waters:[/yellow] Waters with ≥3 connections (network centers)
    • [yellow]Isolated waters:[/yellow] Waters with no H-bond partners
    • [yellow]Water chains:[/yellow] Linear sequences of connected waters

    [bold]Biological significance:[/bold]
    • Large networks → Extensive H-bond networks (active sites, channels)
    • Hub waters → Structurally important coordination centers
    • Isolated waters → Potential removal candidates
    • Long chains → Water wires for proton transfer

    [bold]ASCII visualization:[/bold]
    • Small networks: Detailed connectivity diagrams
    • Large networks: Statistical summaries
    • Connection symbols: ═══ (water-water H-bonds)

    [bold]Applications:[/bold]
    • Identify conserved water networks
    • Find structurally important water clusters
    • Locate potential water channels
    • Guide water molecule retention decisions"""
        
        panel = Panel(explanation, title="Water Network Analysis Guide", border_style="blue", expand=False)
        console.print(panel)

    def interactive_filter(self) -> Optional[Structure]:
        """
        Interactive filtering of PDB structure.

        Returns:
            Filtered structure or None if cancelled
        """
        selected_model_idx = self._get_model_selection()
        selected_model = self.structure[selected_model_idx]

        selected_chain_ids = self._get_chain_selection(selected_model)

        filter_selections: Dict[str, Dict[str, Set[int]]] = {}

        for chain_id in selected_chain_ids:
            chain = selected_model[chain_id]
            chain_selections: Dict[str, Set[int]] = {}
            chain_composition = self.analyze_chain_composition(chain)

            for comp_type, residue_counts in chain_composition.items():
                display_type = ComponentClassifier.display_name(comp_type)

                retention_choice = self._get_component_retention_choice(
                    chain_id, display_type
                )

                if retention_choice == "d":
                    continue

                if retention_choice == "s":
                    selected_residues = self._filter_component_type(chain, comp_type)
                    chain_selections[comp_type] = selected_residues
                else:
                    # Use H++-specific classification for AMBER-named residues
                    if self.is_hplusplus_structure:
                        chain_selections[comp_type] = {
                            residue.id[1]
                            for residue in chain
                            if ComponentClassifier.classify_residue_hplusplus(residue)
                            == comp_type
                        }
                    else:
                        chain_selections[comp_type] = {
                            residue.id[1]
                            for residue in chain
                            if ComponentClassifier.classify_residue(
                                residue, self.ccd_parser
                            )
                            == comp_type
                        }

            filter_selections[chain_id] = chain_selections

        self.filter_selections = filter_selections
        filtered_structure = self._review_selections(
            selected_model_idx, selected_chain_ids, filter_selections
        )
        return filtered_structure

    def _get_model_selection(self) -> int:
        """
        Get model selection (to be overridden by UI layer).

        Returns:
            Selected model index
        """
        models = list(range(len(self.structure)))
        if len(models) == 1:
            return models[0]
        return 0  # Default to first model

    def _get_chain_selection(self, model: Model) -> List[str]:
        """
        Get chain selection (to be overridden by UI layer).

        Args:
            model: Selected model

        Returns:
            List of selected chain IDs
        """
        return [chain.id for chain in model]

    def _get_component_retention_choice(self, chain_id: str, display_type: str) -> str:
        """
        Get component retention choice (to be overridden by UI layer).

        Args:
            chain_id: Chain identifier
            display_type: Display name for component type

        Returns:
            Choice string ('r', 's', or 'd')
        """
        return "r"  # Default to retain

    def _filter_component_type(self, chain: Chain, comp_type: str) -> Set[int]:
        """
        Filter specific component type within a chain (to be overridden by UI layer).

        Args:
            chain: Chain to filter
            comp_type: Component type to filter

        Returns:
            Set of selected residue numbers
        """
        residues = self.get_component_residues(chain, comp_type)
        return {residue.id[1] for residue in residues}

    def _review_selections(
        self,
        model_idx: int,
        chain_ids: List[str],
        filter_selections: Dict[str, Dict[str, Set[int]]],
    ) -> Optional[Structure]:
        """
        Review and confirm filter selections (to be overridden by UI layer).

        Args:
            model_idx: Selected model index
            chain_ids: Selected chain IDs
            filter_selections: Detailed filter selections

        Returns:
            Filtered structure or None if cancelled
        """
        filtered_structure = self.apply_filters(model_idx, chain_ids, filter_selections)
        self.filtered_structure = filtered_structure

        serializable_selections = {}
        for chain_id, chain_data in filter_selections.items():
            serializable_selections[chain_id] = {}
            for comp_type, residue_set in chain_data.items():
                serializable_selections[chain_id][comp_type] = list(residue_set)

        self.filter_selections = serializable_selections
        return filtered_structure

    def analyze_chain_composition(self, chain: Chain, redox_sites: List = None) -> Dict[str, Dict[str, Any]]:
        """
        Analyze the composition of a chain with detailed residue counts and optional redox site information.

        Args:
            chain: Chain to analyze
            redox_sites: Optional list of RedoxSite objects for integration

        Returns:
            Nested dictionary of component types and their residue/redox information
            Format: {component_type: {residue_name: {"count": int, "redox_residues": {resid: [site_ids]}}}}
        """
        composition: Dict[str, Dict[str, Any]] = {}
        
        # Build redox site lookup for this chain
        chain_redox_lookup = {}
        if redox_sites:
            for site in redox_sites:
                # Check centers
                if hasattr(site, 'centers') and site.centers:
                    for center in site.centers:
                        if center.chain == chain.id:
                            key = (center.resname, center.resid)
                            if key not in chain_redox_lookup:
                                chain_redox_lookup[key] = []
                            site_type = getattr(site, 'site_type', 'unknown')
                            chain_redox_lookup[key].append(f"{site.site_id} ({site_type})")
                
                # Check other atoms
                if hasattr(site, 'atoms') and site.atoms:
                    for atom in site.atoms:
                        if atom.chain == chain.id:
                            key = (atom.resname, atom.resid)
                            if key not in chain_redox_lookup:
                                chain_redox_lookup[key] = []
                            site_type = getattr(site, 'site_type', 'unknown')
                            chain_redox_lookup[key].append(f"{site.site_id} ({site_type})")

        for residue in chain:
            # Use H++-specific classification for AMBER-named residues
            if self.is_hplusplus_structure:
                comp_type = ComponentClassifier.classify_residue_hplusplus(residue)
            else:
                comp_type = ComponentClassifier.classify_residue(residue, self.ccd_parser)

            if comp_type not in composition:
                composition[comp_type] = {}

            res_name = residue.resname.strip().upper()
            resid = residue.id[1]
            
            if res_name not in composition[comp_type]:
                composition[comp_type][res_name] = {
                    "count": 0,
                    "redox_residues": {}
                }

            composition[comp_type][res_name]["count"] += 1
            
            # Check if this residue is involved in redox sites
            redox_key = (res_name, resid)
            if redox_key in chain_redox_lookup:
                composition[comp_type][res_name]["redox_residues"][resid] = chain_redox_lookup[redox_key]

        return composition

    def identify_chain_interfaces_by_bsa(
        self, model: Model, area_threshold: float = 200.0
    ) -> Dict[str, List[str]]:
        """
        Identify chain interfaces based on buried surface area using FreeSASA.

        Args:
            model: Model to analyze
            area_threshold: Minimum buried surface area (Å²) to consider chains as interfacing

        Returns:
            Dictionary mapping chain IDs to lists of interfacing chain IDs
        """
        # Skip freeSASA for H++ structures (AMBER residue names may cause errors)
        if self.is_hplusplus_structure:
            logger.debug(
                "H++ structure detected. Using distance-based interface detection "
                "instead of freeSASA (AMBER residue names may not be compatible)."
            )
            return self._identify_interfaces_by_distance(model)

        if not FREESASA_AVAILABLE:
            logger.warning(
                "FreeSASA library not found. Falling back to distance-based interface detection."
            )
            return self._identify_interfaces_by_distance(model)

        interfaces = {chain.id: [] for chain in model}
        interface_areas = {chain.id: {} for chain in model}

        # FreeSASA's Python binding ignores HETATM lines by default, so a chain
        # made up entirely of heteroatoms (cofactor groups, ions, waters)
        # extracts to a file with zero readable atoms and raises
        # "input had no valid ATOM or HETATM lines". Such chains are never
        # protein-protein interface participants anyway, so restrict SASA-based
        # interface detection to polymer chains (those with at least one
        # standard ATOM residue, i.e. residue.id[0] == " ").
        chain_ids = [
            chain.id for chain in model
            if any(residue.id[0] == " " for residue in chain)
        ]
        skipped = [chain.id for chain in model if chain.id not in chain_ids]
        if skipped:
            logger.debug(
                "Skipping non-polymer chain(s) %s for SASA interface detection "
                "(no standard ATOM residues; FreeSASA reads no atoms from them).",
                ", ".join(skipped),
            )

        param_dict = {
            "algorithm": freesasa.LeeRichards,
            "probe-radius": 1.4,
            "n-points": 100,
            "n-slices": 20,
        }

        parameters = freesasa.Parameters(param_dict)
        potential_interfaces = self._prefilter_chain_pairs(model)

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                single_chain_sasa = {}

                for chain_id in chain_ids:
                    chain_path = os.path.join(temp_dir, f"chain_{chain_id}.pdb")

                    class ChainSelect:
                        def __init__(self, chain_id):
                            self.chain_id = chain_id

                        def accept_model(self, model):
                            return True

                        def accept_chain(self, chain):
                            return chain.id == self.chain_id

                        def accept_residue(self, residue):
                            return True

                        def accept_atom(self, atom):
                            return True

                    io = PDBIO()
                    io.set_structure(model.get_parent())
                    io.save(chain_path, ChainSelect(chain_id))

                    try:
                        structure = freesasa.Structure(chain_path)
                        result = freesasa.calc(structure, parameters)
                        single_chain_sasa[chain_id] = result.totalArea()
                    except Exception as e:
                        # Degrade gracefully: skip just this chain rather than
                        # abandoning SASA for the whole structure. Pairs that
                        # reference a chain with no computed SASA are skipped
                        # below (guarded on single_chain_sasa membership).
                        logger.warning(
                            f"FreeSASA failed to process chain {chain_id}: {str(e)}. "
                            "Excluding it from SASA-based interface detection."
                        )
                        continue

                # If we had polymer chains to measure but FreeSASA produced
                # nothing for any of them, SASA is unusable here — fall back to
                # distance-based detection rather than returning empty interfaces.
                if chain_ids and not single_chain_sasa:
                    logger.warning(
                        "FreeSASA produced no single-chain areas for any polymer "
                        "chain. Falling back to distance-based interface detection."
                    )
                    return self._identify_interfaces_by_distance(model)

                for chain1, chain2 in potential_interfaces:
                    # Skip pairs that reference a chain we have no single-chain
                    # SASA for (non-polymer chains, or chains FreeSASA failed on).
                    # _prefilter_chain_pairs proposes pairs from raw atom
                    # proximity, so it can include such chains.
                    if chain1 not in single_chain_sasa or chain2 not in single_chain_sasa:
                        continue

                    pair_path = os.path.join(temp_dir, f"pair_{chain1}_{chain2}.pdb")

                    class PairSelect:
                        def __init__(self, chain1, chain2):
                            self.chain1 = chain1
                            self.chain2 = chain2

                        def accept_model(self, model):
                            return True

                        def accept_chain(self, chain):
                            return chain.id in (self.chain1, self.chain2)

                        def accept_residue(self, residue):
                            return True

                        def accept_atom(self, atom):
                            return True

                    io = PDBIO()
                    io.set_structure(model.get_parent())
                    io.save(pair_path, PairSelect(chain1, chain2))

                    try:
                        structure = freesasa.Structure(pair_path)
                        result = freesasa.calc(structure, parameters)
                        pair_sasa = result.totalArea()
                    except Exception as e:
                        logger.warning(
                            f"FreeSASA failed to process chain pair {chain1}-{chain2}: {str(e)}. "
                            "Skipping this pair in SASA-based interface detection."
                        )
                        continue

                    bsa = single_chain_sasa[chain1] + single_chain_sasa[chain2] - pair_sasa

                    if bsa > area_threshold:
                        interfaces[chain1].append(chain2)
                        interfaces[chain2].append(chain1)
                        interface_areas[chain1][chain2] = bsa
                        interface_areas[chain2][chain1] = bsa
        except Exception as e:
            logger.warning(
                f"FreeSASA BSA calculation failed: {str(e)}. "
                "Falling back to distance-based interface detection."
            )
            return self._identify_interfaces_by_distance(model)

        self.interface_areas = interface_areas
        return interfaces

    def _identify_interfaces_by_distance(
        self, model: Model, distance_cutoff: float = 4.5
    ) -> Dict[str, List[str]]:
        """
        Identify chain interfaces based on distance (fallback method).

        Args:
            model: Model to analyze
            distance_cutoff: Maximum distance (Å) between atoms to consider chains as interfacing

        Returns:
            Dictionary mapping chain IDs to lists of interfacing chain IDs
        """
        interfaces = {chain.id: [] for chain in model}
        potential_interfaces = self._prefilter_chain_pairs(model, distance_cutoff)

        for chain1, chain2 in potential_interfaces:
            interfaces[chain1].append(chain2)
            interfaces[chain2].append(chain1)

        self.interface_areas = {chain_id: {} for chain_id in interfaces}
        for chain1, chain2 in potential_interfaces:
            self.interface_areas[chain1][chain2] = 300
            self.interface_areas[chain2][chain1] = 300

        return interfaces

    def _prefilter_chain_pairs(self, model: Model, distance_cutoff: float = 4.5) -> set:
        """
        Pre-filter chain pairs based on distance to avoid unnecessary BSA calculations.

        Args:
            model: Model to analyze
            distance_cutoff: Maximum distance (Å) between atoms to consider chains as potentially interfacing

        Returns:
            Set of chain ID pairs that might form interfaces
        """
        atoms = []
        atom_chains = {}

        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom)
                    atom_chains[atom] = chain.id

        ns = NeighborSearch(atoms)
        potential_interfaces = set()
        pairs = ns.search_all(distance_cutoff, level="A")

        for atom1, atom2 in pairs:
            chain1 = atom_chains[atom1]
            chain2 = atom_chains[atom2]

            if chain1 != chain2:
                potential_interfaces.add(tuple(sorted([chain1, chain2])))

        return potential_interfaces

    def get_component_residues(self, chain: Chain, comp_type: str) -> List[Residue]:
        """
        Get residues of a specific component type within a chain.

        Args:
            chain: Chain to filter
            comp_type: Component type to filter

        Returns:
            List of residues of the specified type
        """
        residues = []
        for residue in chain:
            # Use H++-specific classification for AMBER-named residues
            if self.is_hplusplus_structure:
                residue_type = ComponentClassifier.classify_residue_hplusplus(residue)
            else:
                residue_type = ComponentClassifier.classify_residue(
                    residue, self.ccd_parser
                )
            if residue_type == comp_type:
                residues.append(residue)
        return residues

    def apply_filters(
        self,
        model_idx: int,
        chain_ids: List[str],
        filter_selections: Dict[str, Dict[str, Set[int]]],
    ) -> Structure:
        """
        Apply filters to the structure.

        Args:
            model_idx: Model index to filter
            chain_ids: Chain IDs to include
            filter_selections: Detailed filter selections

        Returns:
            Filtered structure
        """
        filtered_structure = self.parser.get_structure("filtered", self.filename)
        filtered_structure.detach_parent()

        model = filtered_structure[model_idx]

        for i in reversed(range(len(filtered_structure))):
            if i != model_idx:
                del filtered_structure[i]

        for chain in list(model):
            if chain.id not in chain_ids:
                model.detach_child(chain.id)
                continue

            chain_filters = filter_selections.get(chain.id, {})

            for residue in list(chain):
                # Use H++-specific classification for AMBER-named residues
                if self.is_hplusplus_structure:
                    comp_type = ComponentClassifier.classify_residue_hplusplus(residue)
                else:
                    comp_type = ComponentClassifier.classify_residue(
                        residue, self.ccd_parser
                    )

                if comp_type in chain_filters:
                    if residue.id[1] not in chain_filters[comp_type]:
                        chain.detach_child(residue.id)
                else:
                    chain.detach_child(residue.id)

        return filtered_structure

    def save_filtered_structure(
        self, filtered_structure: Structure, output_filename: str
    ):
        """
        Save filtered structure to a new PDB file.

        Args:
            filtered_structure: Filtered structure to save
            output_filename: Output PDB filename
        """
        io = PDBIO()
        io.set_structure(filtered_structure)
        io.save(output_filename)
        logger.debug(f"Filtered structure saved to {output_filename}")

    def get_structure_info(
        self, structure: Optional[Structure] = None
    ) -> Dict[str, Any]:
        """
        Get information about the structure.

        Args:
            structure: Structure to analyze (defaults to self.structure)

        Returns:
            Dictionary with structure information
        """
        if structure is None:
            structure = self.structure

        if not structure:
            return {}

        info = {"models": len(structure), "chains": {}, "hetero_groups": set()}

        if len(structure) > 0:
            model = structure[0]
            chains_info = {}

            for chain in model:
                chain_id = chain.id
                residue_count = len(chain)

                residue_types = {}
                for residue in chain:
                    # Use H++-specific classification for AMBER-named residues
                    if self.is_hplusplus_structure:
                        comp_type = ComponentClassifier.classify_residue_hplusplus(residue)
                    else:
                        comp_type = ComponentClassifier.classify_residue(
                            residue, self.ccd_parser
                        )

                    if comp_type not in residue_types:
                        residue_types[comp_type] = 0
                    residue_types[comp_type] += 1

                    if residue.id[0] != " ":
                        info["hetero_groups"].add(residue.resname)

                chains_info[chain_id] = {
                    "residue_count": residue_count,
                    "residue_types": residue_types,
                }

            info["chains"] = chains_info

        info["hetero_groups"] = list(info["hetero_groups"])
        return info

    def get_filter_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the filtering process.

        Returns:
            Dictionary with filter statistics
        """
        if not self.filtered_structure or not self.filter_selections:
            return {}

        stats = {
            "original": self.get_structure_info(),
            "filtered": {},
            "removed": {},
            "selections": self.filter_selections,
        }

        if self.filtered_structure:
            stats["filtered"] = self.get_structure_info(self.filtered_structure)

            if "chains" in stats["original"] and "chains" in stats["filtered"]:
                removed_chains = set(stats["original"]["chains"].keys()) - set(
                    stats["filtered"]["chains"].keys()
                )
                stats["removed"]["chains"] = list(removed_chains)

                removed_residues = {}
                for chain_id in stats["original"]["chains"]:
                    if chain_id in stats["filtered"]["chains"]:
                        orig_count = stats["original"]["chains"][chain_id][
                            "residue_count"
                        ]
                        filt_count = stats["filtered"]["chains"][chain_id][
                            "residue_count"
                        ]

                        if orig_count > filt_count:
                            removed_residues[chain_id] = orig_count - filt_count

                stats["removed"]["residues"] = removed_residues

                orig_hetero = set(stats["original"]["hetero_groups"])
                filt_hetero = set(stats["filtered"].get("hetero_groups", []))
                stats["removed"]["hetero_groups"] = list(orig_hetero - filt_hetero)

        return stats

    def get_available_models(self) -> List[int]:
        """Get list of available model indices."""
        return list(range(len(self.structure)))

    def get_available_chains(self, model_idx: int) -> List[str]:
        """Get list of available chain IDs for a model."""
        return [chain.id for chain in self.structure[model_idx]]

    def get_model_chain_info(self, model_idx: int) -> Dict[str, Dict[str, Any]]:
        """
        Get detailed information about chains in a model including topology analysis.

        Args:
            model_idx: Model index

        Returns:
            Dictionary with chain information including interfaces and topology
        """
        model = self.structure[model_idx]
        interfaces = self.identify_chain_interfaces_by_bsa(model)

        chain_info = {}
        for chain in model:
            composition = self.analyze_chain_composition(chain)
            chain_interfaces = interfaces.get(chain.id, [])

            chain_info[chain.id] = {
                "residue_count": len(chain),
                "composition": composition,
                "interfaces": chain_interfaces,
                "interface_areas": self.interface_areas.get(chain.id, {}),
            }

        # Add topology analysis
        try:
            topology_analyzer = ChainTopologyAnalyzer(interfaces, self.interface_areas)
            topology_info = topology_analyzer.get_topology_info()
            linear_sequence = topology_analyzer.determine_linear_sequence()

            # Only add ASCII art for non-monomeric structures
            if topology_info.get('topology_type') != 'monomeric':
                topology_info['ascii_art'] = topology_analyzer.draw_topology_ascii()

            # Add topology info to the result
            chain_info['_topology'] = topology_info
            if linear_sequence:
                chain_info['_topology']['linear_sequence'] = linear_sequence
                chain_info['_topology']['linear_sequence_str'] = ' — '.join(linear_sequence)  # ADD THIS LINE

        except Exception as e:
            logger.warning(f"Could not perform topology analysis: {e}")
            # Add empty topology info as fallback
            chain_info['_topology'] = {
                'is_connected': True,
                'num_chains': len(chain_info),
                'num_interfaces': sum(len(info['interfaces']) for info in chain_info.values()) // 2,
                'terminal_chains': [],
                'branching_points': [],
                'topology_type': 'unknown'
            }

        return chain_info

    def get_filter_results(self) -> Dict[str, Any]:
        """
        Get the results of the filtering process for integration with other modules.

        Returns:
            Dictionary containing filter information and structure
        """
        results = {
            "original_file": self.filename,
            "filtered_structure": self.filtered_structure,
            "filter_selections": self.filter_selections,
            "statistics": self.get_filter_statistics(),
        }

        return results

    def export_filter_statistics(self, output_file: str) -> bool:
        """
        Export filter statistics to a JSON file.

        Args:
            output_file: Output file path

        Returns:
            True if successful, False otherwise
        """
        import json

        stats = self.get_filter_statistics()
        if not stats:
            return False

        try:
            with open(output_file, "w") as f:
                json.dump(stats, f, indent=2)
            return True
        except Exception as e:
            logger.error(f"Error exporting filter statistics: {e}")
            return False

    def _detect_hplusplus_structure(self, filename: str) -> bool:
        """
        Detect if this PDB file was generated by the H++ server.

        H++ structures have REMARK containing "Created by http://biophysics.cs.vt.edu/H++"
        and use AMBER forcefield residue naming conventions (HID, HIE, HIP, CYX, etc.)

        Args:
            filename: Path to the PDB file

        Returns:
            True if H++ structure detected, False otherwise
        """
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith('REMARK'):
                        if 'biophysics.cs.vt.edu/H++' in line or 'H++ server' in line.lower():
                            logger.debug(f"H++ structure detected in {filename}")
                            return True
                    # Stop searching after ATOM/HETATM records start
                    elif line.startswith(('ATOM', 'HETATM')):
                        break
        except Exception as e:
            logger.warning(f"Error detecting H++ structure: {e}")

        return False

    def clear_cache(self):
        """Clear the component classifier cache."""
        ComponentClassifier._ccd_cache.clear()

    def analyze_chain_topology(self, model_idx: int = 0) -> Optional[ChainTopologyAnalyzer]:
        """
        Analyze the topology of chain interfaces in the structure.
        
        Args:
            model_idx: Model index to analyze
            
        Returns:
            ChainTopologyAnalyzer instance or None if no structure loaded
        """
        if not self.structure or model_idx >= len(self.structure):
            return None
        
        model = self.structure[model_idx]
        interfaces = self.identify_chain_interfaces_by_bsa(model)
        
        return ChainTopologyAnalyzer(interfaces, self.interface_areas)

    def get_linear_chain_sequence(self, model_idx: int = 0) -> Optional[List[str]]:
        """
        Get the linear sequence of chains if the structure is linear or has a main backbone.
        
        Args:
            model_idx: Model index to analyze
            
        Returns:
            List of chain IDs in linear order, or None if not applicable
        """
        analyzer = self.analyze_chain_topology(model_idx)
        if analyzer:
            return analyzer.determine_linear_sequence()
        return None

    def get_topology_summary(self, model_idx: int = 0) -> Optional[Dict[str, any]]:
        """
        Get a summary of the chain topology.
        
        Args:
            model_idx: Model index to analyze
            
        Returns:
            Dictionary containing topology information
        """
        analyzer = self.analyze_chain_topology(model_idx)
        if analyzer:
            topology_info = analyzer.get_topology_info()
            linear_sequence = analyzer.determine_linear_sequence()
            if linear_sequence:
                topology_info['linear_sequence'] = linear_sequence
                topology_info['linear_sequence_str'] = ' — '.join(linear_sequence)
            return topology_info
        return None

    def print_topology_analysis(self, model_idx: int = 0):
        """
        Print a detailed topology analysis to console.
        
        Args:
            model_idx: Model index to analyze
        """
        topology_info = self.get_topology_summary(model_idx)
        
        if not topology_info:
            if hasattr(self, 'console') and self.console:
                self.console.print("[yellow]No topology information available[/yellow]")
            else:
                print("No topology information available")
            return
        
        # Use rich console if available, otherwise regular print
        if hasattr(self, 'console') and self.console:
            self.console.print("\n[bold underline]Chain Topology Analysis[/bold underline]")
            # highlight=False so the integer values aren't auto-cyaned.
            self.console.print(f"[blue]Number of chains:[/blue] {topology_info['num_chains']}", highlight=False)
            self.console.print(f"[blue]Number of interfaces:[/blue] {topology_info['num_interfaces']}", highlight=False)
            self.console.print(f"[blue]Topology type:[/blue] {topology_info['topology_type']}", highlight=False)
            self.console.print(f"[blue]Connected:[/blue] {topology_info['is_connected']}", highlight=False)
            
            if topology_info.get('terminal_chains'):
                terminal_str = ', '.join(topology_info['terminal_chains'])
                self.console.print(f"[dark_orange3]Terminal chains:[/dark_orange3] {terminal_str}")
            
            if topology_info.get('branching_points'):
                branching_str = ', '.join(topology_info['branching_points'])
                self.console.print(f"[red]Branching points:[/red] {branching_str}")
            
            if topology_info.get('linear_sequence_str'):
                self.console.print(f"[green]Linear sequence:[/green] {topology_info['linear_sequence_str']}")
        else:
            print("\nChain Topology Analysis")
            print("=" * 30)
            print(f"Number of chains: {topology_info['num_chains']}")
            print(f"Number of interfaces: {topology_info['num_interfaces']}")
            print(f"Topology type: {topology_info['topology_type']}")
            print(f"Connected: {topology_info['is_connected']}")
            
            if topology_info.get('terminal_chains'):
                print(f"Terminal chains: {', '.join(topology_info['terminal_chains'])}")
            
            if topology_info.get('branching_points'):
                print(f"Branching points: {', '.join(topology_info['branching_points'])}")
            
            if topology_info.get('linear_sequence_str'):
                print(f"Linear sequence: {topology_info['linear_sequence_str']}")


# Keep the original PDBFilterTool class for backward compatibility
class PDBFilterTool(PDBFilterWorker):
    """
    Legacy PDBFilterTool class that extends PDBFilterWorker with console UI.
    Maintains backward compatibility while using the new worker architecture.
    """

    def __init__(self, filename: str, existing_structure: Optional[Structure] = None, processor=None):
        """Initialize with console support."""
        super().__init__(filename, existing_structure, processor=processor)
        try:
            from rich.console import Console

            self.console = Console()
        except ImportError:
            self.console = None

    def _get_model_selection(self) -> int:
        """Prompt user to select a model."""
        if not self.console:
            return super()._get_model_selection()

        models = list(range(len(self.structure)))

        if len(models) == 1:
            return models[0]

        self.console.print("\n[bold underline]Available Models[/bold underline]")
        for i, model_idx in enumerate(models):
            self.console.print(f"{i+1}. Model {model_idx}")

        from proprep.utils.prompts import prompt_with_context

        while True:
            model_options_map = {str(i + 1): f"Model {m}" for i, m in enumerate(models)}
            choice = prompt_with_context(
                self.processor,
                "Select model",
                choices=[str(i + 1) for i in range(len(models))],
                default="1",
                module="PDB Filter",
                description="Select PDB model to process",
                options_map=model_options_map,
            )
            return models[int(choice) - 1]

    def _get_chain_selection(self, model: Model) -> List[str]:
        """Prompt user to select chains in the given model."""
        if not self.console:
            return super()._get_chain_selection(model)

        chains = list(model)

        self.console.print(
            "\n[bold underline]Chain Interface Analysis[/bold underline]"
        )
        self.console.print(
            "Calculating buried surface area between chains. This may take a moment..."
        )

        interfaces = self.identify_chain_interfaces_by_bsa(model)

        self.console.print("\n[bold underline]Available Chains[/bold underline]")
        for i, chain in enumerate(chains):
            if interfaces[chain.id]:
                sorted_interfaces = sorted(
                    [
                        (c, self.interface_areas[chain.id].get(c, 0))
                        for c in interfaces[chain.id]
                    ],
                    key=lambda x: x[1],
                    reverse=True,
                )
                interface_text = ", ".join(
                    [f"{c} ({int(area)}Å²)" for c, area in sorted_interfaces]
                )
                interface_display = f"[green]Interfaces with: {interface_text}[/green]"
            else:
                interface_display = "[yellow]No interfaces detected[/yellow]"

            self.console.print(
                f"{i+1}. Chain {chain.id} - {len(chain)} residues - {interface_display}"
            )

        self._display_chain_interface_heatmap(model, interfaces)

        from proprep.utils.prompts import prompt_with_context

        while True:
            chain_options_map = {
                str(i + 1): f"Chain {chain.id} ({len(chain)} residues)"
                for i, chain in enumerate(chains)
            }
            chain_options_map["all"] = "All chains"
            choice = prompt_with_context(
                self.processor,
                "Select chains (comma-separated indices, or 'all')",
                default="all",
                module="PDB Filter",
                description="Select chains (comma-separated indices or 'all')",
                options_map=chain_options_map,
            )

            if choice.lower() == "all":
                return [chain.id for chain in chains]

            try:
                selected_indices = [int(x.strip()) for x in choice.split(",")]

                for idx in selected_indices:
                    if idx < 1 or idx > len(chains):
                        self.console.print(
                            f"[bold red]Invalid chain index: {idx}. Must be between 1 and {len(chains)}[/bold red]"
                        )
                        raise ValueError(f"Chain index out of range: {idx}")

                return [chains[idx - 1].id for idx in selected_indices]

            except ValueError:
                self.console.print(
                    f"[bold red]Invalid input: {choice}. Please enter comma-separated numbers or 'all'[/bold red]"
                )

    def _get_component_retention_choice(self, chain_id: str, display_type: str) -> str:
        """Get component retention choice from user."""
        if not self.console:
            return super()._get_component_retention_choice(chain_id, display_type)

        self.console.print(f"\n[bold]Chain {chain_id} - {display_type}[/bold]")
        self.console.print("\\[r] Retain entire component")
        self.console.print("\\[s] Select specific residues")
        self.console.print("\\[d] Discard entire component")

        from proprep.utils.prompts import prompt_with_context

        return prompt_with_context(
            self.processor,
            "Choose option",
            choices=["r", "s", "d"],
            default="r",
            module="PDB Filter",
            description=f"Retention choice for Chain {chain_id} {display_type}",
            options_map={
                "r": "Retain entire component",
                "s": "Select specific residues",
                "d": "Discard entire component",
            },
        )

    def _filter_component_type(self, chain: Chain, comp_type: str) -> Set[int]:
        """Filter specific component type within a chain with user interaction."""
        if not self.console:
            return super()._filter_component_type(chain, comp_type)

        residues = self.get_component_residues(chain, comp_type)

        display_type = ComponentClassifier.display_name(comp_type)

        from rich.table import Table

        table = Table(title=f"Residues in Chain {chain.id} - {display_type}")
        table.add_column("Select", style="cyan")
        table.add_column("Record Type", style="yellow")
        table.add_column("Residue Name", style="magenta")
        table.add_column("Residue Number", style="green")

        for residue in residues:
            record_type = "Standard" if residue.id[0] == " " else "HETATM"
            table.add_row(
                str(residues.index(residue) + 1),
                record_type,
                residue.resname,
                str(residue.id[1]),
            )

        self.console.print(table)

        from proprep.utils.prompts import prompt_with_context

        residue_options_map = {
            str(i + 1): f"{r.resname} {r.id[1]}" for i, r in enumerate(residues)
        }
        residue_options_map["all"] = "All residues"
        residue_options_map["none"] = "No residues"
        choice = prompt_with_context(
            self.processor,
            "Select residues (comma-separated indices, 'all', or 'none')",
            default="all",
            module="PDB Filter",
            description=f"Select residues in Chain {chain.id} {display_type}",
            options_map=residue_options_map,
        )

        if choice.lower() == "all":
            return {residue.id[1] for residue in residues}
        elif choice.lower() == "none":
            return set()

        selected_indices = [int(x.strip()) - 1 for x in choice.split(",")]
        return {residues[idx].id[1] for idx in selected_indices}

    def _review_selections(
        self,
        model_idx: int,
        chain_ids: List[str],
        filter_selections: Dict[str, Dict[str, Set[int]]],
    ) -> Optional[Structure]:
        """Review and confirm filter selections with user interaction."""
        if not self.console:
            return super()._review_selections(model_idx, chain_ids, filter_selections)

        selections_text = f"Model: {model_idx}\n\n"
        for chain_id, chain_filters in filter_selections.items():
            selections_text += f"Chain {chain_id}:\n"
            if not chain_filters:
                selections_text += "  - No filters applied\n"
            else:
                for comp_type, residues in chain_filters.items():
                    selections_text += (
                        f"  - {comp_type.capitalize()}: {len(residues)} residues\n"
                    )

        from rich.panel import Panel
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        self.console.print(Panel(selections_text, title="Filter Selections", expand=False))

        if confirm_with_context(
            self.processor,
            "Do you want to apply these filters?",
            module="PDB Filter",
            description="Apply chain/residue filter selections",
        ):
            filtered_structure = self.apply_filters(
                model_idx, chain_ids, filter_selections
            )
            self.filtered_structure = filtered_structure

            serializable_selections = {}
            for chain_id, chain_data in filter_selections.items():
                serializable_selections[chain_id] = {}
                for comp_type, residue_set in chain_data.items():
                    serializable_selections[chain_id][comp_type] = list(residue_set)

            self.filter_selections = serializable_selections

            if confirm_with_context(
                self.processor,
                "Do you want to save the filtered structure?",
                module="PDB Filter",
                description="Save filtered structure to disk",
            ):
                output_filename = prompt_with_context(
                    self.processor,
                    "Enter output filename",
                    default="filtered_structure.pdb",
                    module="PDB Filter",
                    description="Output filename for filtered structure",
                )
                self.save_filtered_structure(filtered_structure, output_filename)

            return filtered_structure

        return None

    def _display_topology_analysis(self, model_idx: int = 0):
        """Display topology analysis results."""
        if not hasattr(self, 'console') or not self.console:
            return
        
        try:
            topology_info = self.get_topology_summary(model_idx)
            if not topology_info:
                return
            
            from rich.panel import Panel
            
            info_lines = [
                f"[blue]Chains:[/blue] {topology_info['num_chains']}",
                f"[blue]Interfaces:[/blue] {topology_info['num_interfaces']}",
                f"[blue]Topology:[/blue] {topology_info['topology_type']}",
                f"[blue]Connected:[/blue] {topology_info['is_connected']}"
            ]

            if topology_info.get('terminal_chains'):
                info_lines.append(f"[dark_orange3]Terminals:[/dark_orange3] {', '.join(topology_info['terminal_chains'])}")
            
            if topology_info.get('branching_points'):
                info_lines.append(f"[red]Branches:[/red] {', '.join(topology_info['branching_points'])}")
            
            if topology_info.get('linear_sequence_str'):
                info_lines.append(f"[green]Linear Order:[/green] {topology_info['linear_sequence_str']}")
            
            # Text (not str) so Rich's number highlighter doesn't repaint
            # the panel values in cyan; label colors come from the markup.
            from rich.text import Text
            panel_content = Text.from_markup("\n".join(info_lines))
            self.console.print(Panel(panel_content, title="Chain Topology", border_style="blue", expand=False))
            
        except Exception as e:
            logger.warning(f"Could not display topology analysis: {e}")

    def _display_chain_interface_heatmap(
        self, model: Model, interfaces: Dict[str, List[str]]
    ):
        """Display a visual heatmap of chain interfaces."""
        if not self.console:
            return

        from rich.table import Table

        chains = sorted([chain.id for chain in model])

        # header_style="bold blue" so the chain-ID headers and row labels
        # don't ride the bold-default-foreground, which is invisible on a
        # white background (see pdb_filter.py for the full rationale).
        table = Table(
            title="Chain Interface Map (Buried Surface Area in Å²)",
            header_style="bold blue",
        )
        table.add_column("")
        for chain_id in chains:
            table.add_column(f"{chain_id}", justify="center")

        for chain1 in chains:
            row = [f"[bold blue]{chain1}[/bold blue]"]
            for chain2 in chains:
                if chain1 == chain2:
                    row.append("[grey50]■[/grey50]")
                elif chain2 in interfaces[chain1]:
                    bsa = int(self.interface_areas[chain1][chain2])

                    if bsa > 2000:
                        row.append(f"[bold red]{bsa}[/bold red]")
                    elif bsa > 1000:
                        row.append(f"[bold green]{bsa}[/bold green]")
                    elif bsa > 500:
                        row.append(f"[green]{bsa}[/green]")
                    elif bsa > 200:
                        row.append(f"[dark_orange3]{bsa}[/dark_orange3]")
                    else:
                        row.append(f"[grey50]{bsa}[/grey50]")
                else:
                    row.append("[grey50]0[/grey50]")

            table.add_row(*row)

        self.console.print(table)
        self.console.print("\n[bold]BSA Color Legend:[/bold]")
        self.console.print("[grey50]< 200 Å²[/grey50]: Minimal contact")
        self.console.print("[dark_orange3]200-500 Å²[/dark_orange3]: Small interface")
        self.console.print("[green]500-1000 Å²[/green]: Medium interface")
        self.console.print("[bold green]1000-2000 Å²[/bold green]: Large interface")
        self.console.print("[bold red]> 2000 Å²[/bold red]: Very large interface")

        self._display_topology_analysis(0)

    def _display_chain_composition(
        self, chain_composition: Dict[str, Dict[str, int]], chain_id: str
    ):
        """Display a hierarchical view of chain composition."""
        if not self.console:
            return

        from rich.table import Table

        table = Table(
            title=f"Chain {chain_id} Composition", show_lines=False,
            header_style="bold blue",
        )
        # Cap the width so long CCD chemical names wrap to several short lines
        # instead of stretching the column far wider than the others.
        table.add_column("Component Type", style="bold blue", max_width=18, overflow="fold")
        table.add_column("Residue\nName", style="default", justify="center")
        table.add_column("Count", style="default", justify="center")

        standard_types = ["amino_acid", "dna_base", "rna_base", "water"]

        standard_components = {
            t: chain_composition.get(t, {})
            for t in standard_types
            if t in chain_composition
        }

        special_components = {
            t: c for t, c in chain_composition.items() if t not in standard_types
        }

        for comp_type, residues in standard_components.items():
            sorted_residues = sorted(residues.items(), key=lambda x: x[1], reverse=True)
            display_type = ComponentClassifier.display_name(comp_type)
            type_total = sum(count for _, count in sorted_residues)
            table.add_row(
                display_type, "[Total]", str(type_total)
            )

            for residue, count in sorted_residues:
                table.add_row("", residue, str(count))

            if sorted_residues:
                table.add_row("", "", "")

        for comp_type, residues in sorted(special_components.items()):
            sorted_residues = sorted(residues.items(), key=lambda x: x[1], reverse=True)
            type_total = sum(count for _, count in sorted_residues)
            table.add_row(comp_type, "[Total]", str(type_total))

            for residue, count in sorted_residues:
                table.add_row("", residue, str(count))

            if sorted_residues:
                table.add_row("", "", "")

        self.console.print(table)

    def display_filter_summary(self):
        """Display a summary of the filtering process"""
        if not self.console:
            return

        if not self.filtered_structure or not self.filter_selections:
            self.console.print("[yellow]No filtering has been performed yet[/yellow]")
            return

        stats = self.get_filter_statistics()

        self.console.print("\n[bold]Original Structure Summary[/bold]")
        self.console.print(f"Models: {stats['original']['models']}")
        self.console.print(f"Chains: {len(stats['original']['chains'])}")

        total_residues = sum(
            chain["residue_count"] for chain in stats["original"]["chains"].values()
        )
        self.console.print(f"Total residues: {total_residues}")

        self.console.print("\n[bold]Filtered Structure Summary[/bold]")
        self.console.print(f"Models: {stats['filtered']['models']}")
        self.console.print(f"Chains: {len(stats['filtered']['chains'])}")

        total_residues = sum(
            chain["residue_count"] for chain in stats["filtered"]["chains"].values()
        )
        self.console.print(f"Total residues: {total_residues}")

        self.console.print("\n[bold]Removed Components[/bold]")
        if stats["removed"].get("chains"):
            self.console.print(f"Chains: {', '.join(stats['removed']['chains'])}")

        if stats["removed"].get("residues"):
            self.console.print("Residues:")
            for chain_id, count in stats["removed"]["residues"].items():
                self.console.print(f"  Chain {chain_id}: {count} residues removed")

        if stats["removed"].get("hetero_groups"):
            self.console.print(
                f"Hetero groups: {', '.join(stats['removed']['hetero_groups'])}"
            )

