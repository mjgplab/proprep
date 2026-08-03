"""
Trajectory Analysis for AMBER MD Simulations

Mirrors AMBERMonitor architecture for mdout files but operates on trajectory (.nc) files.
Provides structural analysis to complement energetic analysis.
"""

import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from collections import defaultdict

logger = logging.getLogger(__name__)

try:
    import pytraj as pt
except ImportError:
    pt = None
    logger.warning("pytraj not available - trajectory analysis will be disabled")


class TrajectoryAnalyzer:
    """
    Analyze AMBER trajectory files for structural properties.

    Mirrors AMBERMonitor pattern for consistency with mdout analysis.
    Stores data in self.data dictionary, provides statistics methods,
    and supports region-based analysis.
    """

    def __init__(self, topology: Union[str, Path], trajectory: Union[str, Path] = None,
                 traj_object=None):
        """
        Initialize trajectory analyzer.

        Args:
            topology: Path to .prmtop topology file
            trajectory: Path to .nc trajectory file (for single file)
            traj_object: Pre-loaded pytraj trajectory object (for multi-file concatenation)

        Raises:
            ImportError: If pytraj is not available
            FileNotFoundError: If topology or trajectory file not found
        """
        if pt is None:
            raise ImportError(
                "pytraj is required for trajectory analysis. "
                "Install with: conda install -c conda-forge pytraj ambertools"
            )

        self.topology = str(topology)

        # Load trajectory
        if traj_object is not None:
            # Pre-loaded trajectory (from concatenation)
            self.traj = traj_object
            logger.info(f"Initialized with pre-loaded trajectory: {len(traj_object)} frames")
        elif trajectory is not None:
            # Load single trajectory file
            trajectory = str(trajectory)
            logger.info(f"Loading trajectory from {trajectory}")
            self.traj = pt.load(trajectory, top=self.topology)
            logger.info(f"Loaded {len(self.traj)} frames")
        else:
            raise ValueError("Either trajectory path or traj_object must be provided")

        # Data storage (mirroring AMBERMonitor.data pattern)
        self.data = defaultdict(list)

        # System information
        self.system_info = {}
        self._extract_system_info()

        # Frame timing information
        self.frame_times = []  # ps
        self._extract_timing_info()

    def _extract_system_info(self):
        """Extract system composition from topology."""
        self.system_info = {
            'n_frames': len(self.traj),
            'n_atoms': self.traj.top.n_atoms,
            'n_residues': self.traj.top.n_residues,
            'box_type': 'periodic' if self.traj.top.box is not None else 'non-periodic'
        }

        # Try to determine protein atoms
        try:
            protein_mask = self.traj.top.select('!(:WAT,Na+,Cl-,K+)')
            self.system_info['n_protein_atoms'] = len(protein_mask)
        except:
            self.system_info['n_protein_atoms'] = self.system_info['n_atoms']

        logger.info(f"System: {self.system_info['n_atoms']} atoms, "
                   f"{self.system_info['n_residues']} residues, "
                   f"{self.system_info['n_frames']} frames")

    def _extract_timing_info(self):
        """Extract frame timing from trajectory metadata."""
        # Try to get times from trajectory
        try:
            if hasattr(self.traj, 'time'):
                self.frame_times = self.traj.time.tolist()
            else:
                # Fall back to frame indices * 2 ps (common default)
                self.frame_times = [i * 2.0 for i in range(len(self.traj))]
        except Exception as e:
            logger.warning(f"Could not extract timing info: {e}, using defaults")
            self.frame_times = [i * 2.0 for i in range(len(self.traj))]

    def calculate_rmsd(self, mask: str = '@CA', reference: int = 0,
                      label: str = None) -> np.ndarray:
        """
        Calculate RMSD over trajectory.

        Args:
            mask: AMBER atom selection mask (e.g., '@CA', '@C,CA,N,O', ':1-50@CA')
            reference: Reference frame index (0-based), or -1 for average structure
            label: Storage label (auto-generated if None)

        Returns:
            np.ndarray: RMSD values per frame (Angstroms)
        """
        try:
            logger.info(f"Calculating RMSD with mask '{mask}', reference frame {reference}")

            # Calculate RMSD using pytraj
            rmsd_values = pt.rmsd(self.traj, mask=mask, ref=reference)

            # Store with descriptive label
            key = label or f'rmsd_{mask.replace(",", "_").replace(":", "_").replace("@", "")}'
            self.data[key] = rmsd_values.tolist()
            self.data[f'{key}_mask'] = mask
            self.data[f'{key}_reference'] = reference

            logger.info(f"RMSD calculation complete: {len(rmsd_values)} values")
            return rmsd_values

        except Exception as e:
            logger.error(f"Error calculating RMSD: {e}")
            raise

    def calculate_rmsf(self, mask: str = '@CA', label: str = None) -> np.ndarray:
        """
        Calculate per-atom/residue RMSF (root mean square fluctuation).

        Args:
            mask: AMBER atom selection mask
            label: Storage label (auto-generated if None)

        Returns:
            np.ndarray: RMSF values per atom in selection (Angstroms)
        """
        try:
            logger.info(f"Calculating RMSF with mask '{mask}'")

            # RMSF requires an aligned trajectory
            # Create a copy of the trajectory to avoid modifying the original
            aligned_traj = self.traj[:]

            # Align trajectory to first frame using backbone atoms for alignment
            # This removes overall translation and rotation
            # Note: atomicfluct/rmsf does NOT automatically align, so we must do it explicitly
            alignment_mask = '@CA,C,N'  # Use backbone atoms for alignment
            pt.superpose(aligned_traj, ref=0, mask=alignment_mask)

            # Calculate per-residue RMSF using pytraj with byres option
            # Returns array of shape (n_residues, 2) where columns are [residue_index, rmsf_value]
            # byres calculates mass-weighted average RMSF for each residue
            rmsf_data = pt.rmsf(aligned_traj, mask=mask, options='byres')

            # Extract just the RMSF values (column 1), not the residue indices (column 0)
            rmsf_array = rmsf_data[:, 1]

            # Store with descriptive label
            key = label or f'rmsf_{mask.replace(",", "_").replace(":", "_").replace("@", "")}'
            self.data[key] = rmsf_array.tolist() if hasattr(rmsf_array, 'tolist') else list(rmsf_array)
            self.data[f'{key}_mask'] = mask
            # Store residue indices for proper labeling
            self.data[f'{key}_residue_indices'] = rmsf_data[:, 0].astype(int).tolist()

            logger.info(f"RMSF calculation complete: {len(rmsf_array)} residues, range: {rmsf_array.min():.2f}-{rmsf_array.max():.2f} Å")
            return rmsf_array

        except Exception as e:
            logger.error(f"Error calculating RMSF: {e}")
            raise

    def calculate_distance(self, mask1: str, mask2: str,
                          label: str = None) -> np.ndarray:
        """
        Calculate distance between two atom selections over trajectory.

        Args:
            mask1: First atom selection (e.g., ':42@CA')
            mask2: Second atom selection (e.g., ':89@CA')
            label: Storage label (auto-generated if None)

        Returns:
            np.ndarray: Distance per frame (Angstroms)
        """
        try:
            logger.info(f"Calculating distance between '{mask1}' and '{mask2}'")

            # Calculate distance using pytraj
            # Format: "mask1 mask2"
            distances = pt.distance(self.traj, f"{mask1} {mask2}")

            # Store with descriptive label
            key = label or f'distance_{mask1}_{mask2}'.replace(':', '_').replace('@', '')
            self.data[key] = distances.tolist()
            self.data[f'{key}_mask1'] = mask1
            self.data[f'{key}_mask2'] = mask2

            logger.info(f"Distance calculation complete: {len(distances)} values")
            return distances

        except Exception as e:
            logger.error(f"Error calculating distance: {e}")
            raise

    def calculate_angle(self, mask1: str, mask2: str, mask3: str,
                       label: str = None) -> np.ndarray:
        """
        Calculate angle between three atom selections over trajectory.

        Args:
            mask1: First atom selection
            mask2: Second atom selection (vertex)
            mask3: Third atom selection
            label: Storage label

        Returns:
            np.ndarray: Angle per frame (degrees)
        """
        try:
            logger.info(f"Calculating angle {mask1}-{mask2}-{mask3}")

            # Calculate angle using pytraj
            # Format: "mask1 mask2 mask3"
            angles = pt.angle(self.traj, f"{mask1} {mask2} {mask3}")

            # Store with descriptive label
            key = label or f'angle_{mask1}_{mask2}_{mask3}'.replace(':', '_').replace('@', '')
            self.data[key] = angles.tolist()
            self.data[f'{key}_masks'] = [mask1, mask2, mask3]

            logger.info(f"Angle calculation complete: {len(angles)} values")
            return angles

        except Exception as e:
            logger.error(f"Error calculating angle: {e}")
            raise

    def calculate_dihedral(self, mask1: str, mask2: str, mask3: str, mask4: str,
                          label: str = None) -> np.ndarray:
        """
        Calculate dihedral angle between four atom selections over trajectory.

        Args:
            mask1: First atom selection
            mask2: Second atom selection
            mask3: Third atom selection
            mask4: Fourth atom selection
            label: Storage label

        Returns:
            np.ndarray: Dihedral per frame (degrees, -180 to 180)
        """
        try:
            logger.info(f"Calculating dihedral {mask1}-{mask2}-{mask3}-{mask4}")

            # Calculate dihedral using pytraj
            # Format: "mask1 mask2 mask3 mask4"
            dihedrals = pt.dihedral(self.traj, f"{mask1} {mask2} {mask3} {mask4}")

            # Store with descriptive label
            key = label or f'dihedral_{mask1}_{mask2}_{mask3}_{mask4}'.replace(':', '_').replace('@', '')
            self.data[key] = dihedrals.tolist()
            self.data[f'{key}_masks'] = [mask1, mask2, mask3, mask4]

            logger.info(f"Dihedral calculation complete: {len(dihedrals)} values")
            return dihedrals

        except Exception as e:
            logger.error(f"Error calculating dihedral: {e}")
            raise

    def calculate_hbonds(self, donor_mask: str = None, acceptor_mask: str = None,
                        distance_cutoff: float = 3.0, angle_cutoff: float = 135.0,
                        label: str = None) -> dict:
        """
        Calculate hydrogen bonds over trajectory.

        Args:
            donor_mask: Donor atom selection (None = all)
            acceptor_mask: Acceptor atom selection (None = all)
            distance_cutoff: Maximum donor-acceptor distance (Angstroms)
            angle_cutoff: Minimum donor-H-acceptor angle (degrees)
            label: Storage label

        Returns:
            dict: H-bond analysis results including counts per frame
        """
        try:
            logger.info(f"Calculating H-bonds (distance<{distance_cutoff}, angle>{angle_cutoff})")

            # Build options string for donor/acceptor masks
            options = ""
            if donor_mask:
                options += f"donormask {donor_mask} "
            if acceptor_mask:
                options += f"acceptormask {acceptor_mask}"

            # Calculate H-bonds using pytraj
            # Returns DatasetHBond object with .values[0] = total H-bonds per frame
            hbond_data = pt.hbond(
                self.traj,
                distance=distance_cutoff,
                angle=angle_cutoff,
                options=options.strip(),
                series=True  # Return time series data
            )

            # Extract total H-bonds per frame from first row of values array
            # hbond_data.values shape: (n_hbonds+1, n_frames)
            # Row 0 = total counts, rows 1+ = individual donor-acceptor pairs
            hbond_counts = hbond_data.values[0]

            # Calculate occupancy for each H-bond pair
            # Occupancy = percentage of frames where H-bond is present
            n_frames = len(self.traj)
            pair_occupancies = []

            if hasattr(hbond_data, 'donor_acceptor') and hasattr(hbond_data, 'values'):
                donor_acceptor_pairs = hbond_data.donor_acceptor
                for i, pair_name in enumerate(donor_acceptor_pairs, start=1):
                    # Row i contains 0/1 for each frame (present/absent)
                    pair_data = hbond_data.values[i]
                    occupancy = (np.sum(pair_data) / n_frames) * 100  # Percentage
                    pair_occupancies.append({
                        'pair': pair_name,
                        'occupancy': occupancy,
                        'frames_present': int(np.sum(pair_data))
                    })

                # Sort by occupancy (most persistent first)
                pair_occupancies.sort(key=lambda x: x['occupancy'], reverse=True)

            # Store results
            key = label or 'hbonds'
            self.data[key] = hbond_counts.tolist() if hasattr(hbond_counts, 'tolist') else list(hbond_counts)
            self.data[f'{key}_distance_cutoff'] = distance_cutoff
            self.data[f'{key}_angle_cutoff'] = angle_cutoff
            # Store donor-acceptor pairs with occupancy data
            self.data[f'{key}_pair_occupancies'] = pair_occupancies

            logger.info(f"H-bond calculation complete: {len(pair_occupancies)} unique pairs found")
            return hbond_data

        except Exception as e:
            logger.error(f"Error calculating H-bonds: {e}")
            # For now, store zeros if H-bond calc fails
            key = label or 'hbonds'
            self.data[key] = [0] * len(self.traj)
            logger.warning("H-bond calculation failed, storing zeros")
            return {}

    def calculate_radius_of_gyration(self, mask: str = '@CA',
                                     label: str = None) -> np.ndarray:
        """
        Calculate radius of gyration over trajectory.

        Args:
            mask: Atom selection (typically '@CA' or protein mask)
            label: Storage label

        Returns:
            np.ndarray: Rg per frame (Angstroms)
        """
        try:
            logger.info(f"Calculating radius of gyration with mask '{mask}'")

            # Calculate Rg using pytraj
            rgyr = pt.radgyr(self.traj, mask=mask)

            # Store with descriptive label
            key = label or f'rgyr_{mask.replace(",", "_").replace(":", "_").replace("@", "")}'
            self.data[key] = rgyr.tolist()
            self.data[f'{key}_mask'] = mask

            logger.info(f"Rg calculation complete: {len(rgyr)} values")
            return rgyr

        except Exception as e:
            logger.error(f"Error calculating Rg: {e}")
            raise

    def calculate_water_radial_distribution(self, solute_mask: str,
                                           solvent_mask: str = ':WAT@O',
                                           max_distance: float = 10.0,
                                           bin_spacing: float = 0.1,
                                           label: str = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate radial distribution function (RDF) for water around solute.

        Args:
            solute_mask: Selection for solute (e.g., ':42-45', ':ALA')
            solvent_mask: Selection for solvent (default: water oxygens)
            max_distance: Maximum distance for RDF (Angstroms)
            bin_spacing: Bin width (Angstroms)
            label: Storage label

        Returns:
            tuple: (distances, g(r) values)
        """
        try:
            logger.info(f"Calculating RDF: '{solute_mask}' to '{solvent_mask}'")

            # Calculate RDF using pytraj
            rdf_data = pt.radial(
                self.traj,
                mask=f"{solute_mask} {solvent_mask}",
                spacing=bin_spacing,
                maximum=max_distance
            )

            # rdf_data is typically a 2D array where:
            # column 0 = distances, column 1 = g(r)
            distances = rdf_data[:, 0]
            gr_values = rdf_data[:, 1]

            # Store results
            key = label or f'rdf_{solute_mask}_to_{solvent_mask}'.replace(':', '_').replace('@', '')
            self.data[f'{key}_distances'] = distances.tolist()
            self.data[f'{key}_gr'] = gr_values.tolist()
            self.data[f'{key}_solute_mask'] = solute_mask
            self.data[f'{key}_solvent_mask'] = solvent_mask

            logger.info(f"RDF calculation complete: {len(distances)} bins")
            return distances, gr_values

        except Exception as e:
            logger.error(f"Error calculating RDF: {e}")
            raise

    def calculate_sasa(self, mask: str = '*', probe_radius: float = 1.4,
                      label: str = None) -> np.ndarray:
        """
        Calculate solvent accessible surface area over trajectory.

        Args:
            mask: Atom selection (default: all atoms)
            probe_radius: Probe radius (Angstroms, 1.4 for water)
            label: Storage label

        Returns:
            np.ndarray: SASA per frame (Angstrom^2)
        """
        try:
            logger.info(f"Calculating SASA with mask '{mask}', probe radius {probe_radius}")

            # Calculate SASA using pytraj
            sasa = pt.molsurf(self.traj, mask=mask, probe_radius=probe_radius)

            # Store with descriptive label
            key = label or f'sasa_{mask.replace(",", "_").replace(":", "_").replace("@", "")}'
            self.data[key] = sasa.tolist()
            self.data[f'{key}_mask'] = mask
            self.data[f'{key}_probe_radius'] = probe_radius

            logger.info(f"SASA calculation complete: {len(sasa)} values")
            return sasa

        except Exception as e:
            logger.error(f"Error calculating SASA: {e}")
            raise

    def calculate_dssp(self, mask: str = None, label: str = None) -> dict:
        """
        Calculate secondary structure using DSSP algorithm.

        Args:
            mask: Residue selection (default: all protein residues)
            label: Storage label

        Returns:
            dict: DSSP results including per-residue assignments
        """
        try:
            logger.info("Calculating secondary structure (DSSP)")

            # DSSP assigns secondary structure codes:
            # H = α-helix, E = β-sheet, T = turn, C = coil
            # Additional: G = 3-10 helix, I = π-helix, B = β-bridge, S = bend

            # Calculate DSSP
            # Returns: array of secondary structure codes per frame
            dssp_data = pt.dssp(self.traj, mask=mask if mask else '', simplified=True)

            # Process DSSP output to get per-frame statistics
            # dssp_data format: (n_frames, n_residues) array of SS codes

            # Count secondary structure types per frame
            n_frames = len(self.traj)
            helix_pct = []
            sheet_pct = []
            turn_pct = []
            coil_pct = []

            for frame_ss in dssp_data:
                # frame_ss is array of SS codes for this frame
                total = len(frame_ss)
                if total == 0:
                    helix_pct.append(0.0)
                    sheet_pct.append(0.0)
                    turn_pct.append(0.0)
                    coil_pct.append(0.0)
                    continue

                # Count each type (simplified codes: H=helix, E=sheet, T=turn, C=coil)
                n_helix = np.sum(frame_ss == 'H')
                n_sheet = np.sum(frame_ss == 'E')
                n_turn = np.sum(frame_ss == 'T')
                n_coil = np.sum(frame_ss == 'C')

                helix_pct.append((n_helix / total) * 100)
                sheet_pct.append((n_sheet / total) * 100)
                turn_pct.append((n_turn / total) * 100)
                coil_pct.append((n_coil / total) * 100)

            # Store results
            key = label or 'dssp'
            self.data[f'{key}_helix_pct'] = helix_pct
            self.data[f'{key}_sheet_pct'] = sheet_pct
            self.data[f'{key}_turn_pct'] = turn_pct
            self.data[f'{key}_coil_pct'] = coil_pct
            self.data[f'{key}_per_residue'] = dssp_data.tolist()  # Full per-residue assignments
            self.data[f'{key}_mask'] = mask if mask else 'all protein'

            logger.info(f"DSSP calculation complete: {len(dssp_data[0])} residues analyzed")
            return {
                'helix_pct': helix_pct,
                'sheet_pct': sheet_pct,
                'turn_pct': turn_pct,
                'coil_pct': coil_pct,
                'per_residue': dssp_data
            }

        except Exception as e:
            logger.error(f"Error calculating DSSP: {e}")
            logger.warning("DSSP calculation requires dssp executable in PATH")
            logger.warning("Install via: conda install -c salilab dssp")
            raise

    def calculate_ramachandran(self, residue_selection: str = None,
                              dihedral_type: str = 'phi-psi',
                              label: str = None) -> dict:
        """
        Calculate Ramachandran angles and other backbone/sidechain dihedrals.

        Args:
            residue_selection: Residue range (e.g., '1-156', None = all)
            dihedral_type: 'phi-psi', 'chi1', 'chi2', 'omega'
            label: Storage label

        Returns:
            dict: Dihedral angles per residue over time
        """
        try:
            logger.info(f"Calculating {dihedral_type} dihedral angles")

            # Determine residue range
            if residue_selection is None:
                # Auto-detect protein residues (first 80% if >100 residues)
                n_residues = self.traj.top.n_residues
                if n_residues > 100:
                    # Conservative estimate: first 80% are protein
                    protein_end = int(n_residues * 0.8)
                    resrange = f'1-{protein_end}'
                else:
                    resrange = f'1-{n_residues}'
            else:
                resrange = residue_selection

            results = {}

            if dihedral_type == 'phi-psi':
                # Calculate both phi and psi angles using pytraj calc functions
                logger.info(f"Calculating phi angles for residues {resrange}")
                phi_data = pt.calc_phi(self.traj, resrange=resrange)

                logger.info(f"Calculating psi angles for residues {resrange}")
                psi_data = pt.calc_psi(self.traj, resrange=resrange)

                # phi_data and psi_data are DatasetList objects
                # Convert to numpy arrays: shape (n_residues, n_frames)
                phi_array = np.array([data.values for data in phi_data])
                psi_array = np.array([data.values for data in psi_data])

                # Store results
                key = label or 'ramachandran'
                self.data[f'{key}_phi'] = phi_array.tolist()
                self.data[f'{key}_psi'] = psi_array.tolist()
                self.data[f'{key}_resrange'] = resrange
                self.data[f'{key}_type'] = 'phi-psi'

                results = {
                    'phi': phi_array,
                    'psi': psi_array,
                    'resrange': resrange
                }

                logger.info(f"Ramachandran calculation complete: {len(phi_data)} residues analyzed")

            elif dihedral_type in ['chi1', 'chi2', 'chi3', 'chi4']:
                # Calculate sidechain chi angles
                logger.info(f"Calculating {dihedral_type} angles for residues {resrange}")
                chi_type = int(dihedral_type[-1])  # Extract number from 'chi1', 'chi2', etc.
                chi_data = pt.calc_chi(self.traj, resrange=resrange, chi_type=chi_type)

                # Convert to numpy array
                chi_array = np.array([data.values for data in chi_data])

                # Store results
                key = label or f'dihedral_{dihedral_type}'
                self.data[f'{key}_values'] = chi_array.tolist()
                self.data[f'{key}_resrange'] = resrange
                self.data[f'{key}_type'] = dihedral_type

                results = {
                    'chi': chi_array,
                    'resrange': resrange,
                    'type': dihedral_type
                }

                logger.info(f"{dihedral_type} calculation complete: {len(chi_data)} residues analyzed")

            elif dihedral_type == 'omega':
                # Calculate omega (peptide bond planarity)
                logger.info(f"Calculating omega angles for residues {resrange}")
                omega_data = pt.calc_omega(self.traj, resrange=resrange)

                omega_array = np.array([data.values for data in omega_data])

                key = label or 'omega'
                self.data[f'{key}_values'] = omega_array.tolist()
                self.data[f'{key}_resrange'] = resrange
                self.data[f'{key}_type'] = 'omega'

                results = {
                    'omega': omega_array,
                    'resrange': resrange
                }

                logger.info(f"Omega calculation complete: {len(omega_data)} residues analyzed")

            return results

        except Exception as e:
            logger.error(f"Error calculating {dihedral_type} dihedrals: {e}")
            raise

    def calculate_contacts(self,
                          mask1: str = None,
                          mask2: str = None,
                          distance_cutoff: float = 4.5,
                          reference_frame: int = 0,
                          label: str = None) -> dict:
        """
        Calculate contact maps and native contact analysis.

        Args:
            mask1: First selection mask (default: protein C-alpha atoms)
            mask2: Second selection mask (default: same as mask1 for intra-protein)
            distance_cutoff: Distance cutoff in Angstroms for contact definition
            reference_frame: Frame to use as reference for native contacts (default: 0)
            label: Optional label for data storage

        Returns:
            dict: Contact analysis results including Q-value and contact matrix
        """
        try:
            # Default to C-alpha atoms if not specified
            if mask1 is None:
                mask1 = '@CA'
            if mask2 is None:
                mask2 = mask1

            logger.info(f"Calculating contacts between '{mask1}' and '{mask2}' with cutoff {distance_cutoff} Å")

            # Get atom indices for selections
            atoms1 = self.traj.top.select(mask1)
            atoms2 = self.traj.top.select(mask2)

            n_atoms1 = len(atoms1)
            n_atoms2 = len(atoms2)
            n_frames = self.traj.n_frames

            logger.info(f"Selection 1: {n_atoms1} atoms, Selection 2: {n_atoms2} atoms")

            # Initialize contact matrix (n_atoms1 x n_atoms2 x n_frames)
            # We'll compute contacts frame by frame to avoid memory issues
            contact_matrix = np.zeros((n_atoms1, n_atoms2, n_frames), dtype=bool)

            # Calculate distances for all frames
            for frame_idx in range(n_frames):
                frame = self.traj[frame_idx]

                for i, atom_i in enumerate(atoms1):
                    for j, atom_j in enumerate(atoms2):
                        # Skip self-contacts if same selection
                        if mask1 == mask2 and atom_i == atom_j:
                            continue

                        # Calculate distance using pytraj
                        dist = pt.distance(frame, f'@{atom_i+1} @{atom_j+1}')[0]

                        if dist <= distance_cutoff:
                            contact_matrix[i, j, frame_idx] = True

            # Calculate native contacts from reference frame
            native_contacts = contact_matrix[:, :, reference_frame]
            n_native_contacts = np.sum(native_contacts)

            logger.info(f"Found {n_native_contacts} native contacts in reference frame {reference_frame}")

            # Calculate Q-value for each frame (fraction of native contacts maintained)
            q_values = []
            for frame_idx in range(n_frames):
                frame_contacts = contact_matrix[:, :, frame_idx]
                # Contacts that are native AND present in this frame
                maintained = np.logical_and(native_contacts, frame_contacts)
                n_maintained = np.sum(maintained)

                if n_native_contacts > 0:
                    q_value = n_maintained / n_native_contacts
                else:
                    q_value = 1.0  # No native contacts defined

                q_values.append(q_value)

            # Calculate contact frequency (occupancy) for each pair
            contact_frequency = np.sum(contact_matrix, axis=2) / n_frames

            # Find most persistent contacts (occupancy > 0.5)
            persistent_contacts = []
            for i in range(n_atoms1):
                for j in range(n_atoms2):
                    occupancy = contact_frequency[i, j]
                    if occupancy > 0.5:
                        res_i = self.traj.top.atom(atoms1[i]).resid + 1  # 1-indexed
                        res_j = self.traj.top.atom(atoms2[j]).resid + 1
                        resname_i = self.traj.top.atom(atoms1[i]).resname
                        resname_j = self.traj.top.atom(atoms2[j]).resname
                        atom_i_name = self.traj.top.atom(atoms1[i]).name
                        atom_j_name = self.traj.top.atom(atoms2[j]).name

                        persistent_contacts.append({
                            'res_i': res_i,
                            'res_j': res_j,
                            'resname_i': resname_i,
                            'resname_j': resname_j,
                            'atom_i': atom_i_name,
                            'atom_j': atom_j_name,
                            'occupancy': occupancy,
                            'is_native': bool(native_contacts[i, j])
                        })

            # Sort by occupancy
            persistent_contacts.sort(key=lambda x: x['occupancy'], reverse=True)

            # Store results
            key = label or f'contacts_{mask1}_{mask2}'
            key = key.replace('@', '').replace(',', '_').replace(' ', '_')  # Clean key

            self.data[f'{key}_q_values'] = q_values
            self.data[f'{key}_contact_frequency'] = contact_frequency.tolist()
            self.data[f'{key}_persistent_contacts'] = persistent_contacts
            self.data[f'{key}_n_native_contacts'] = int(n_native_contacts)
            self.data[f'{key}_mask1'] = mask1
            self.data[f'{key}_mask2'] = mask2
            self.data[f'{key}_cutoff'] = distance_cutoff
            self.data[f'{key}_reference_frame'] = reference_frame

            results = {
                'q_values': q_values,
                'contact_frequency': contact_frequency,
                'persistent_contacts': persistent_contacts,
                'n_native_contacts': n_native_contacts,
                'mask1': mask1,
                'mask2': mask2
            }

            logger.info(f"Contact analysis complete: {len(persistent_contacts)} persistent contacts found")

            return results

        except Exception as e:
            logger.error(f"Error calculating contacts: {e}")
            raise

    def calculate_salt_bridges(self,
                               distance_cutoff: float = 4.0,
                               label: str = None) -> dict:
        """
        Calculate salt bridge interactions between charged residues.

        Identifies electrostatic interactions between acidic (ASP, GLU) and
        basic (LYS, ARG, HIS) residues.

        Args:
            distance_cutoff: Maximum distance in Angstroms for salt bridge (default: 4.0)
            label: Optional label for data storage

        Returns:
            dict: Salt bridge analysis results including occupancy and distances
        """
        try:
            logger.info(f"Calculating salt bridges with cutoff {distance_cutoff} Å")

            # Define charged residues
            acidic_residues = ['ASP', 'GLU']
            basic_residues = ['LYS', 'ARG', 'HIS']

            # Define charged atoms for each residue type
            acidic_atoms = {
                'ASP': ['OD1', 'OD2'],
                'GLU': ['OE1', 'OE2']
            }
            basic_atoms = {
                'LYS': ['NZ'],
                'ARG': ['NH1', 'NH2'],
                'HIS': ['ND1', 'NE2']
            }

            # Find all acidic and basic residues in the system
            acidic_res_list = []
            basic_res_list = []

            for res in self.traj.top.residues:
                if res.name in acidic_residues:
                    # Get indices of charged atoms
                    charged_atoms = [atom.index for atom in res.atoms
                                   if atom.name in acidic_atoms.get(res.name, [])]
                    if charged_atoms:
                        acidic_res_list.append({
                            'resid': res.resid + 1,  # 1-indexed
                            'resname': res.name,
                            'atoms': charged_atoms
                        })

                elif res.name in basic_residues:
                    charged_atoms = [atom.index for atom in res.atoms
                                   if atom.name in basic_atoms.get(res.name, [])]
                    if charged_atoms:
                        basic_res_list.append({
                            'resid': res.resid + 1,
                            'resname': res.name,
                            'atoms': charged_atoms
                        })

            logger.info(f"Found {len(acidic_res_list)} acidic and {len(basic_res_list)} basic residues")

            n_frames = self.traj.n_frames

            # Track all potential salt bridge pairs
            salt_bridge_pairs = []
            for acidic in acidic_res_list:
                for basic in basic_res_list:
                    # Skip if same residue
                    if acidic['resid'] == basic['resid']:
                        continue

                    pair_key = f"{acidic['resname']}{acidic['resid']}-{basic['resname']}{basic['resid']}"
                    salt_bridge_pairs.append({
                        'key': pair_key,
                        'acidic_res': acidic,
                        'basic_res': basic,
                        'distances': [],  # min distance per frame
                        'formed': []  # boolean per frame
                    })

            # Calculate distances for each frame
            for frame_idx in range(n_frames):
                frame = self.traj[frame_idx]

                for pair in salt_bridge_pairs:
                    acidic_atoms = pair['acidic_res']['atoms']
                    basic_atoms = pair['basic_res']['atoms']

                    # Find minimum distance between any acidic and basic atom in this pair
                    min_dist = float('inf')
                    for acid_atom in acidic_atoms:
                        for base_atom in basic_atoms:
                            dist = pt.distance(frame, f'@{acid_atom+1} @{base_atom+1}')[0]
                            if dist < min_dist:
                                min_dist = dist

                    pair['distances'].append(min_dist)
                    pair['formed'].append(min_dist <= distance_cutoff)

            # Calculate occupancy and statistics for each salt bridge
            salt_bridge_results = []
            for pair in salt_bridge_pairs:
                occupancy = sum(pair['formed']) / n_frames
                distances = np.array(pair['distances'])

                # Only include if formed at least once
                if occupancy > 0:
                    salt_bridge_results.append({
                        'pair': pair['key'],
                        'acidic_res': pair['acidic_res']['resid'],
                        'acidic_name': pair['acidic_res']['resname'],
                        'basic_res': pair['basic_res']['resid'],
                        'basic_name': pair['basic_res']['resname'],
                        'occupancy': occupancy,
                        'mean_distance': float(np.mean(distances)),
                        'std_distance': float(np.std(distances)),
                        'min_distance': float(np.min(distances)),
                        'distances': distances.tolist()
                    })

            # Sort by occupancy
            salt_bridge_results.sort(key=lambda x: x['occupancy'], reverse=True)

            # Store results
            key = label or 'salt_bridges'
            self.data[f'{key}_results'] = salt_bridge_results
            self.data[f'{key}_cutoff'] = distance_cutoff
            self.data[f'{key}_n_bridges'] = len(salt_bridge_results)

            results = {
                'salt_bridges': salt_bridge_results,
                'cutoff': distance_cutoff,
                'n_bridges': len(salt_bridge_results)
            }

            logger.info(f"Salt bridge analysis complete: {len(salt_bridge_results)} salt bridges found")

            return results

        except Exception as e:
            logger.error(f"Error calculating salt bridges: {e}")
            raise

    def calculate_pca(self,
                     mask: str = None,
                     n_components: int = 3,
                     label: str = None) -> dict:
        """
        Calculate Principal Component Analysis (PCA) for trajectory.

        Identifies collective motions by reducing dimensionality.

        Args:
            mask: AMBER mask for atoms to include (default: C-alpha atoms)
            n_components: Number of principal components to calculate (default: 3)
            label: Optional label for data storage

        Returns:
            dict: PCA results including projections, eigenvalues, and variance explained
        """
        try:
            # Default to C-alpha atoms if not specified
            if mask is None:
                mask = '@CA'

            logger.info(f"Calculating PCA for mask '{mask}' with {n_components} components")

            # Perform PCA using pytraj
            # First, create a dataset for PCA
            pca_data = pt.pca(self.traj, mask=mask, n_vecs=n_components)

            # Get eigenvalues and eigenvectors
            eigenvalues = pca_data[0]  # Array of eigenvalues
            eigenvectors = pca_data[1]  # Matrix of eigenvectors

            # Project trajectory onto principal components
            # Use pytraj projection
            projections = pt.projection(self.traj, mask=mask, eigenvalues=eigenvalues,
                                       eigenvectors=eigenvectors, scalar_type='covar')

            # Calculate variance explained by each component
            total_variance = np.sum(eigenvalues)
            variance_explained = [(ev / total_variance) * 100 for ev in eigenvalues]
            cumulative_variance = np.cumsum(variance_explained)

            # Store projection data for each component
            pc_projections = {}
            for i in range(min(n_components, len(projections))):
                pc_projections[f'PC{i+1}'] = projections[i].tolist() if hasattr(projections[i], 'tolist') else list(projections[i])

            # Store results
            key = label or f'pca_{mask.replace("@", "").replace(",", "_")}'

            self.data[f'{key}_eigenvalues'] = eigenvalues.tolist() if hasattr(eigenvalues, 'tolist') else list(eigenvalues)
            self.data[f'{key}_variance_explained'] = variance_explained
            self.data[f'{key}_cumulative_variance'] = cumulative_variance.tolist()
            self.data[f'{key}_projections'] = pc_projections
            self.data[f'{key}_mask'] = mask
            self.data[f'{key}_n_components'] = n_components

            results = {
                'eigenvalues': eigenvalues,
                'variance_explained': variance_explained,
                'cumulative_variance': cumulative_variance,
                'projections': pc_projections,
                'mask': mask,
                'n_components': n_components
            }

            logger.info(f"PCA complete: PC1 explains {variance_explained[0]:.1f}% of variance")

            return results

        except Exception as e:
            logger.error(f"Error calculating PCA: {e}")
            raise

    def calculate_clustering(self,
                            mask: str = None,
                            n_clusters: int = 5,
                            algorithm: str = 'kmeans',
                            label: str = None) -> dict:
        """
        Calculate conformational clustering of trajectory.

        Groups similar conformations into clusters.

        Args:
            mask: AMBER mask for atoms to include (default: C-alpha atoms)
            n_clusters: Number of clusters (default: 5)
            algorithm: Clustering algorithm - 'kmeans' or 'hierarchical' (default: kmeans)
            label: Optional label for data storage

        Returns:
            dict: Clustering results including assignments, populations, and representatives
        """
        try:
            # Default to C-alpha atoms if not specified
            if mask is None:
                mask = '@CA'

            logger.info(f"Calculating {algorithm} clustering for mask '{mask}' with {n_clusters} clusters")

            # Perform clustering using pytraj
            if algorithm.lower() == 'kmeans':
                # K-means clustering
                cluster_data = pt.cluster.kmeans(
                    self.traj,
                    n_clusters=n_clusters,
                    mask=mask,
                    metric='rms'
                )
            elif algorithm.lower() == 'hierarchical':
                # Hierarchical clustering
                cluster_data = pt.cluster.hierarchical(
                    self.traj,
                    n_clusters=n_clusters,
                    mask=mask,
                    metric='rms'
                )
            else:
                raise ValueError(f"Unknown clustering algorithm: {algorithm}")

            # Extract cluster assignments for each frame
            cluster_assignments = cluster_data.cluster_index

            # Calculate cluster populations
            cluster_counts = {}
            for i in range(n_clusters):
                count = np.sum(cluster_assignments == i)
                cluster_counts[i] = count

            # Get representative frames (centroids)
            representative_frames = cluster_data.representative_frames

            # Calculate RMSD within each cluster
            cluster_rmsd = {}
            for cluster_id in range(n_clusters):
                # Get frames in this cluster
                cluster_frames_idx = np.where(cluster_assignments == cluster_id)[0]

                if len(cluster_frames_idx) > 1:
                    # Calculate average RMSD to centroid
                    centroid_idx = representative_frames[cluster_id]
                    rmsds = []

                    for frame_idx in cluster_frames_idx:
                        rmsd = pt.rmsd(
                            self.traj[frame_idx],
                            self.traj[centroid_idx],
                            mask=mask
                        )[0]
                        rmsds.append(rmsd)

                    cluster_rmsd[cluster_id] = {
                        'mean': float(np.mean(rmsds)),
                        'std': float(np.std(rmsds)),
                        'max': float(np.max(rmsds))
                    }
                else:
                    cluster_rmsd[cluster_id] = {'mean': 0.0, 'std': 0.0, 'max': 0.0}

            # Store results
            key = label or f'clustering_{algorithm}'

            self.data[f'{key}_assignments'] = cluster_assignments.tolist()
            self.data[f'{key}_populations'] = cluster_counts
            self.data[f'{key}_representatives'] = representative_frames.tolist()
            self.data[f'{key}_cluster_rmsd'] = cluster_rmsd
            self.data[f'{key}_n_clusters'] = n_clusters
            self.data[f'{key}_algorithm'] = algorithm
            self.data[f'{key}_mask'] = mask

            results = {
                'assignments': cluster_assignments,
                'populations': cluster_counts,
                'representatives': representative_frames,
                'cluster_rmsd': cluster_rmsd,
                'n_clusters': n_clusters,
                'algorithm': algorithm,
                'mask': mask
            }

            logger.info(f"Clustering complete: {n_clusters} clusters identified")

            return results

        except Exception as e:
            logger.error(f"Error calculating clustering: {e}")
            raise

    def calculate_pairwise_rmsd(self,
                               mask: str = None,
                               subsample: int = None,
                               label: str = None) -> dict:
        """
        Calculate pairwise RMSD matrix for all frames.

        Computes all-vs-all RMSD to identify conformational similarity.

        Args:
            mask: AMBER mask for atoms to include (default: C-alpha atoms)
            subsample: Subsample every Nth frame (default: None, use all frames)
            label: Optional label for data storage

        Returns:
            dict: Pairwise RMSD results including matrix and statistics
        """
        try:
            # Default to C-alpha atoms if not specified
            if mask is None:
                mask = '@CA'

            # Handle subsampling
            if subsample and subsample > 1:
                traj_to_use = self.traj[::subsample]
                logger.info(f"Subsampling trajectory: using every {subsample}th frame ({traj_to_use.n_frames} frames)")
            else:
                traj_to_use = self.traj

            n_frames = traj_to_use.n_frames

            logger.info(f"Calculating pairwise RMSD matrix for {n_frames} frames with mask '{mask}'")

            # Initialize RMSD matrix
            rmsd_matrix = np.zeros((n_frames, n_frames))

            # Calculate pairwise RMSDs using pytraj
            # Use rmsd_nofit=False to include alignment
            for i in range(n_frames):
                for j in range(i, n_frames):
                    if i == j:
                        rmsd_matrix[i, j] = 0.0
                    else:
                        rmsd_val = pt.rmsd(traj_to_use[i], traj_to_use[j], mask=mask)[0]
                        rmsd_matrix[i, j] = rmsd_val
                        rmsd_matrix[j, i] = rmsd_val  # Symmetric

            # Calculate statistics
            # Average RMSD for each frame (to all other frames)
            avg_rmsd_per_frame = np.mean(rmsd_matrix, axis=1)

            # Find most "central" frame (lowest average RMSD to all others)
            most_central_frame = int(np.argmin(avg_rmsd_per_frame))
            central_frame_rmsd = float(avg_rmsd_per_frame[most_central_frame])

            # Find most "extreme" frame (highest average RMSD)
            most_extreme_frame = int(np.argmax(avg_rmsd_per_frame))
            extreme_frame_rmsd = float(avg_rmsd_per_frame[most_extreme_frame])

            # Overall statistics
            upper_triangle = rmsd_matrix[np.triu_indices(n_frames, k=1)]
            overall_mean = float(np.mean(upper_triangle))
            overall_std = float(np.std(upper_triangle))
            overall_max = float(np.max(upper_triangle))

            # Store results
            key = label or f'pairwise_rmsd'

            self.data[f'{key}_matrix'] = rmsd_matrix.tolist()
            self.data[f'{key}_avg_per_frame'] = avg_rmsd_per_frame.tolist()
            self.data[f'{key}_most_central'] = most_central_frame
            self.data[f'{key}_most_extreme'] = most_extreme_frame
            self.data[f'{key}_overall_mean'] = overall_mean
            self.data[f'{key}_overall_std'] = overall_std
            self.data[f'{key}_overall_max'] = overall_max
            self.data[f'{key}_mask'] = mask
            self.data[f'{key}_subsample'] = subsample

            results = {
                'matrix': rmsd_matrix,
                'avg_per_frame': avg_rmsd_per_frame,
                'most_central': most_central_frame,
                'central_rmsd': central_frame_rmsd,
                'most_extreme': most_extreme_frame,
                'extreme_rmsd': extreme_frame_rmsd,
                'overall_mean': overall_mean,
                'overall_std': overall_std,
                'overall_max': overall_max,
                'mask': mask,
                'n_frames': n_frames
            }

            logger.info(f"Pairwise RMSD complete: overall mean = {overall_mean:.2f} Å")

            return results

        except Exception as e:
            logger.error(f"Error calculating pairwise RMSD: {e}")
            raise

    def calculate_bfactors(self,
                          mask: str = None,
                          by_residue: bool = True,
                          label: str = None) -> dict:
        """
        Calculate pseudo B-factors from MD trajectory.

        B-factors are calculated from atomic fluctuations (RMSF):
        B = (8π²/3) * RMSF²

        Args:
            mask: AMBER mask for atoms to include (default: C-alpha atoms)
            by_residue: Calculate per-residue average (default: True)
            label: Optional label for data storage

        Returns:
            dict: B-factor results including values and residue mapping
        """
        try:
            # Default to C-alpha atoms if not specified
            if mask is None:
                mask = '@CA'

            logger.info(f"Calculating B-factors for mask '{mask}'")

            # Calculate RMSF
            rmsf = pt.rmsf(self.traj, mask=mask)

            # Convert RMSF (Å) to B-factors (Å²)
            # B = (8π²/3) * RMSF²
            conversion_factor = (8 * np.pi**2) / 3
            bfactors = conversion_factor * (rmsf ** 2)

            # Get atom/residue information
            atoms = self.traj.top.select(mask)
            residue_ids = []
            residue_names = []

            for atom_idx in atoms:
                atom = self.traj.top.atom(atom_idx)
                residue_ids.append(atom.resid + 1)  # 1-indexed
                residue_names.append(atom.resname)

            # Calculate per-residue if requested
            if by_residue:
                unique_resids = sorted(set(residue_ids))
                residue_bfactors = []
                residue_rmsf = []

                for resid in unique_resids:
                    # Find atoms belonging to this residue
                    res_indices = [i for i, rid in enumerate(residue_ids) if rid == resid]

                    # Average B-factor for this residue
                    res_bf = np.mean([bfactors[i] for i in res_indices])
                    res_rmsf_val = np.mean([rmsf[i] for i in res_indices])

                    residue_bfactors.append(res_bf)
                    residue_rmsf.append(res_rmsf_val)

                # Store residue-level results
                key = label or f'bfactors_{mask.replace("@", "").replace(",", "_")}'

                self.data[f'{key}_residue_ids'] = unique_resids
                self.data[f'{key}_bfactors'] = residue_bfactors
                self.data[f'{key}_rmsf'] = residue_rmsf
                self.data[f'{key}_mask'] = mask
                self.data[f'{key}_by_residue'] = True

                results = {
                    'bfactors': np.array(residue_bfactors),
                    'rmsf': np.array(residue_rmsf),
                    'residue_ids': unique_resids,
                    'residue_names': [residue_names[residue_ids.index(rid)] for rid in unique_resids],
                    'mask': mask,
                    'by_residue': True
                }

            else:
                # Store per-atom results
                key = label or f'bfactors_{mask.replace("@", "").replace(",", "_")}'

                self.data[f'{key}_atom_bfactors'] = bfactors.tolist()
                self.data[f'{key}_atom_rmsf'] = rmsf.tolist()
                self.data[f'{key}_residue_ids'] = residue_ids
                self.data[f'{key}_mask'] = mask
                self.data[f'{key}_by_residue'] = False

                results = {
                    'bfactors': bfactors,
                    'rmsf': rmsf,
                    'residue_ids': residue_ids,
                    'residue_names': residue_names,
                    'mask': mask,
                    'by_residue': False
                }

            logger.info(f"B-factor calculation complete: mean B-factor = {np.mean(results['bfactors']):.2f} Å²")

            return results

        except Exception as e:
            logger.error(f"Error calculating B-factors: {e}")
            raise

    def calculate_water_shells(self,
                               solute_mask: str = None,
                               shell_width: float = 2.0,
                               max_distance: float = 10.0,
                               label: str = None) -> dict:
        """
        Calculate water shell analysis around solute.

        Analyzes hydration layers by counting water molecules in
        concentric shells at different distances from solute.

        Args:
            solute_mask: AMBER mask for solute (default: protein)
            shell_width: Width of each shell in Angstroms (default: 2.0)
            max_distance: Maximum distance for analysis (default: 10.0)
            label: Optional label for data storage

        Returns:
            dict: Water shell results including populations and occupancies
        """
        try:
            # Default to protein if not specified
            if solute_mask is None:
                solute_mask = ':1-500 & !@H='  # Protein heavy atoms (guess protein range)

            logger.info(f"Calculating water shells around '{solute_mask}'")

            # Define shell boundaries
            n_shells = int(max_distance / shell_width)
            shell_boundaries = [(i * shell_width, (i + 1) * shell_width) for i in range(n_shells)]

            # Get solute and water atoms
            solute_atoms = self.traj.top.select(solute_mask)

            # Find water oxygen atoms (common water residue names)
            water_mask = ':WAT,HOH,TIP3,TIP4,TIP5,SPC,OPC@O'
            water_atoms = self.traj.top.select(water_mask)

            if len(water_atoms) == 0:
                logger.warning("No water molecules found with standard names")
                # Try alternative approach - find residues named water-like
                water_atoms = []
                for res in self.traj.top.residues:
                    if res.name in ['WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPC', 'OPC']:
                        # Get oxygen atom
                        for atom in res.atoms:
                            if atom.name in ['O', 'OW', 'OH2']:
                                water_atoms.append(atom.index)
                                break

                water_atoms = np.array(water_atoms)

            logger.info(f"Found {len(solute_atoms)} solute atoms and {len(water_atoms)} water molecules")

            n_frames = self.traj.n_frames

            # Initialize shell populations for each frame
            shell_populations = np.zeros((n_frames, n_shells))

            # Calculate water counts per shell for each frame
            for frame_idx in range(n_frames):
                frame = self.traj[frame_idx]

                for water_idx in water_atoms:
                    # Find minimum distance from this water to any solute atom
                    min_dist = float('inf')

                    for solute_idx in solute_atoms:
                        dist = pt.distance(frame, f'@{water_idx+1} @{solute_idx+1}')[0]
                        if dist < min_dist:
                            min_dist = dist

                    # Assign water to appropriate shell
                    for shell_idx, (lower, upper) in enumerate(shell_boundaries):
                        if lower <= min_dist < upper:
                            shell_populations[frame_idx, shell_idx] += 1
                            break

            # Calculate statistics for each shell
            shell_stats = []
            for shell_idx in range(n_shells):
                lower, upper = shell_boundaries[shell_idx]
                populations = shell_populations[:, shell_idx]

                shell_stats.append({
                    'range': f"{lower:.1f}-{upper:.1f} Å",
                    'mean': float(np.mean(populations)),
                    'std': float(np.std(populations)),
                    'min': int(np.min(populations)),
                    'max': int(np.max(populations))
                })

            # Store results
            key = label or 'water_shells'

            self.data[f'{key}_populations'] = shell_populations.tolist()
            self.data[f'{key}_stats'] = shell_stats
            self.data[f'{key}_shell_boundaries'] = shell_boundaries
            self.data[f'{key}_solute_mask'] = solute_mask
            self.data[f'{key}_n_shells'] = n_shells

            results = {
                'populations': shell_populations,
                'stats': shell_stats,
                'shell_boundaries': shell_boundaries,
                'solute_mask': solute_mask,
                'n_shells': n_shells
            }

            logger.info(f"Water shell analysis complete: {n_shells} shells analyzed")

            return results

        except Exception as e:
            logger.error(f"Error calculating water shells: {e}")
            raise

    def calculate_density_map(self,
                             selection_mask: str,
                             grid_spacing: float = 1.0,
                             label: str = None) -> dict:
        """
        Calculate 3D density map for selected atoms.

        Creates spatial density distribution showing where atoms
        spend most time during simulation.

        Args:
            selection_mask: AMBER mask for atoms to analyze
            grid_spacing: Grid spacing in Angstroms (default: 1.0)
            label: Optional label for data storage

        Returns:
            dict: Density map results including grid and peak locations
        """
        try:
            logger.info(f"Calculating density map for '{selection_mask}' with {grid_spacing} Å spacing")

            # Get selected atoms
            atoms = self.traj.top.select(selection_mask)

            if len(atoms) == 0:
                raise ValueError(f"No atoms selected with mask: {selection_mask}")

            logger.info(f"Selected {len(atoms)} atoms")

            # Determine grid boundaries from trajectory box
            # Get box dimensions (assumes periodic box)
            box_lengths = []
            for frame in self.traj:
                if hasattr(frame, 'box_lengths') and frame.box_lengths is not None:
                    box_lengths.append(frame.box_lengths)

            if len(box_lengths) == 0:
                # No box info, use atom positions to define grid
                logger.warning("No box information found, using atom positions")
                all_coords = []
                for frame in self.traj:
                    coords = frame.xyz[:, atoms, :]
                    all_coords.append(coords)

                all_coords = np.concatenate(all_coords, axis=0)
                min_coords = np.min(all_coords, axis=0)
                max_coords = np.max(all_coords, axis=0)

                # Add padding
                padding = 5.0  # Angstroms
                grid_min = min_coords - padding
                grid_max = max_coords + padding
            else:
                # Use average box dimensions
                avg_box = np.mean(box_lengths, axis=0)
                grid_min = np.array([0.0, 0.0, 0.0])
                grid_max = avg_box[:3]

            # Create 3D grid
            nx = int((grid_max[0] - grid_min[0]) / grid_spacing) + 1
            ny = int((grid_max[1] - grid_min[1]) / grid_spacing) + 1
            nz = int((grid_max[2] - grid_min[2]) / grid_spacing) + 1

            # Limit grid size for memory
            max_grid_dim = 100
            if nx > max_grid_dim or ny > max_grid_dim or nz > max_grid_dim:
                logger.warning(f"Grid too large ({nx}x{ny}x{nz}), adjusting spacing")
                new_spacing = max((grid_max[0] - grid_min[0]) / max_grid_dim,
                                (grid_max[1] - grid_min[1]) / max_grid_dim,
                                (grid_max[2] - grid_min[2]) / max_grid_dim)
                grid_spacing = new_spacing
                nx = int((grid_max[0] - grid_min[0]) / grid_spacing) + 1
                ny = int((grid_max[1] - grid_min[1]) / grid_spacing) + 1
                nz = int((grid_max[2] - grid_min[2]) / grid_spacing) + 1

            logger.info(f"Grid dimensions: {nx}x{ny}x{nz}")

            # Initialize density grid
            density = np.zeros((nx, ny, nz))

            # Accumulate density
            for frame in self.traj:
                coords = frame.xyz[0, atoms, :]  # Shape: (n_atoms, 3)

                for coord in coords:
                    # Convert to grid indices
                    ix = int((coord[0] - grid_min[0]) / grid_spacing)
                    iy = int((coord[1] - grid_min[1]) / grid_spacing)
                    iz = int((coord[2] - grid_min[2]) / grid_spacing)

                    # Check bounds
                    if 0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz:
                        density[ix, iy, iz] += 1

            # Normalize by number of frames and number of atoms
            density = density / (self.traj.n_frames * len(atoms))

            # Find peak density locations
            flat_density = density.flatten()
            sorted_indices = np.argsort(flat_density)[::-1]

            # Get top 5 peak locations
            peaks = []
            for idx in sorted_indices[:5]:
                if flat_density[idx] > 0:
                    # Convert back to 3D indices
                    ix = idx // (ny * nz)
                    iy = (idx % (ny * nz)) // nz
                    iz = idx % nz

                    # Convert to real coordinates
                    x = grid_min[0] + ix * grid_spacing
                    y = grid_min[1] + iy * grid_spacing
                    z = grid_min[2] + iz * grid_spacing

                    peaks.append({
                        'coords': (x, y, z),
                        'density': float(flat_density[idx])
                    })

            # Store results
            key = label or f'density_{selection_mask.replace("@", "").replace(":", "_")}'

            # Store grid info and peaks (not full density for memory)
            self.data[f'{key}_grid_min'] = grid_min.tolist()
            self.data[f'{key}_grid_max'] = grid_max.tolist()
            self.data[f'{key}_grid_spacing'] = grid_spacing
            self.data[f'{key}_grid_dims'] = [nx, ny, nz]
            self.data[f'{key}_peaks'] = peaks
            self.data[f'{key}_max_density'] = float(np.max(density))
            self.data[f'{key}_selection'] = selection_mask

            # Create 2D slice for visualization (middle Z plane)
            mid_z = nz // 2
            density_slice = density[:, :, mid_z]

            results = {
                'density': density,
                'density_slice': density_slice,
                'grid_min': grid_min,
                'grid_max': grid_max,
                'grid_spacing': grid_spacing,
                'grid_dims': (nx, ny, nz),
                'peaks': peaks,
                'max_density': float(np.max(density)),
                'selection': selection_mask
            }

            logger.info(f"Density map complete: max density = {np.max(density):.4f}")

            return results

        except Exception as e:
            logger.error(f"Error calculating density map: {e}")
            raise

    def calculate_vector_analysis(self,
                                  atom1_mask: str,
                                  atom2_mask: str,
                                  label: str = None) -> dict:
        """
        Calculate vector orientation analysis.

        Tracks orientation of vector defined by two atom selections.
        Useful for monitoring helix tilt, bond rotations, molecular axis.

        Args:
            atom1_mask: AMBER mask for vector start
            atom2_mask: AMBER mask for vector end
            label: Optional label for data storage

        Returns:
            dict: Vector analysis results including angles and magnitudes
        """
        try:
            logger.info(f"Calculating vector analysis: '{atom1_mask}' → '{atom2_mask}'")

            # Get atoms
            atoms1 = self.traj.top.select(atom1_mask)
            atoms2 = self.traj.top.select(atom2_mask)

            if len(atoms1) == 0 or len(atoms2) == 0:
                raise ValueError(f"No atoms selected")

            # Calculate center of mass for each selection per frame
            n_frames = self.traj.n_frames
            vectors = []
            magnitudes = []
            polar_angles = []  # Angle from Z axis
            azimuthal_angles = []  # Angle in XY plane from X axis

            for frame in self.traj:
                # Get coordinates
                coords1 = frame.xyz[0, atoms1, :]
                coords2 = frame.xyz[0, atoms2, :]

                # Center of mass (simplified - assume equal mass)
                com1 = np.mean(coords1, axis=0)
                com2 = np.mean(coords2, axis=0)

                # Vector from com1 to com2
                vec = com2 - com1
                vectors.append(vec)

                # Magnitude
                mag = np.linalg.norm(vec)
                magnitudes.append(mag)

                # Polar angle (from Z axis)
                if mag > 0:
                    polar = np.arccos(vec[2] / mag) * (180 / np.pi)
                else:
                    polar = 0.0
                polar_angles.append(polar)

                # Azimuthal angle (in XY plane from X axis)
                azimuthal = np.arctan2(vec[1], vec[0]) * (180 / np.pi)
                azimuthal_angles.append(azimuthal)

            vectors = np.array(vectors)
            magnitudes = np.array(magnitudes)
            polar_angles = np.array(polar_angles)
            azimuthal_angles = np.array(azimuthal_angles)

            # Store results
            key = label or 'vector'

            self.data[f'{key}_magnitudes'] = magnitudes.tolist()
            self.data[f'{key}_polar_angles'] = polar_angles.tolist()
            self.data[f'{key}_azimuthal_angles'] = azimuthal_angles.tolist()
            self.data[f'{key}_vectors'] = vectors.tolist()
            self.data[f'{key}_atom1_mask'] = atom1_mask
            self.data[f'{key}_atom2_mask'] = atom2_mask

            results = {
                'magnitudes': magnitudes,
                'polar_angles': polar_angles,
                'azimuthal_angles': azimuthal_angles,
                'vectors': vectors,
                'atom1_mask': atom1_mask,
                'atom2_mask': atom2_mask
            }

            logger.info(f"Vector analysis complete: mean magnitude = {np.mean(magnitudes):.2f} Å")

            return results

        except Exception as e:
            logger.error(f"Error calculating vector analysis: {e}")
            raise

    def calculate_autocorrelation(self,
                                  data_series: np.ndarray,
                                  max_lag: int = None,
                                  label: str = None) -> dict:
        """
        Calculate autocorrelation function for a time series.

        Identifies periodic motions and relaxation timescales.

        Args:
            data_series: Time series data (e.g., RMSD, distance)
            max_lag: Maximum lag time (default: half of series length)
            label: Optional label for data storage

        Returns:
            dict: Autocorrelation results
        """
        try:
            n = len(data_series)

            if max_lag is None:
                max_lag = n // 2

            logger.info(f"Calculating autocorrelation with max_lag={max_lag}")

            # Calculate autocorrelation
            mean_val = np.mean(data_series)
            var_val = np.var(data_series)

            autocorr = np.zeros(max_lag)

            for lag in range(max_lag):
                if lag < n:
                    autocorr[lag] = np.mean((data_series[:n-lag] - mean_val) *
                                           (data_series[lag:] - mean_val)) / var_val

            # Store results
            key = label or 'autocorrelation'
            self.data[f'{key}_values'] = autocorr.tolist()
            self.data[f'{key}_max_lag'] = max_lag

            results = {
                'autocorr': autocorr,
                'max_lag': max_lag
            }

            logger.info("Autocorrelation calculation complete")

            return results

        except Exception as e:
            logger.error(f"Error calculating autocorrelation: {e}")
            raise

    def calculate_contact_frequency(self,
                                   mask: str = None,
                                   distance_cutoff: float = 4.5,
                                   label: str = None) -> dict:
        """
        Calculate per-residue contact frequency.

        Shows which residues make the most contacts.

        Args:
            mask: AMBER mask for residues (default: protein)
            distance_cutoff: Distance cutoff for contacts (default: 4.5 Å)
            label: Optional label for data storage

        Returns:
            dict: Per-residue contact frequency results
        """
        try:
            import numpy as np
            from scipy.spatial.distance import cdist

            if mask is None:
                mask = '@CA'

            logger.info(f"Calculating per-residue contact frequency")

            atoms = self.traj.top.select(mask)
            n_atoms = len(atoms)
            n_frames = self.traj.n_frames

            # Map atoms to residues
            residue_map = {}
            for atom_idx in atoms:
                atom = self.traj.top.atom(atom_idx)
                resid = atom.resid + 1  # 1-indexed
                if resid not in residue_map:
                    residue_map[resid] = []
                residue_map[resid].append(atom_idx)

            residues = sorted(residue_map.keys())
            n_residues = len(residues)

            logger.info(f"Analyzing {n_residues} residues across {n_frames} frames")

            # Initialize contact counts
            contact_counts = {resid: 0 for resid in residues}

            # Calculate contacts per frame using vectorized operations
            for frame_idx in range(n_frames):
                # Get coordinates for this frame
                coords = self.traj.xyz[frame_idx]  # Shape: (n_atoms_total, 3)

                # Extract coordinates for selected atoms
                selected_coords = coords[atoms]  # Shape: (n_atoms, 3)

                # Calculate distance matrix for all selected atoms at once
                # This is MUCH faster than calling pt.distance repeatedly
                dist_matrix = cdist(selected_coords, selected_coords)

                # For each residue pair, check if any atoms are within cutoff
                for i, res_i in enumerate(residues):
                    atoms_i_indices = [atoms.tolist().index(a) for a in residue_map[res_i]]

                    for j, res_j in enumerate(residues):
                        if i >= j:  # Skip self and duplicates
                            continue

                        atoms_j_indices = [atoms.tolist().index(a) for a in residue_map[res_j]]

                        # Extract sub-matrix for this residue pair
                        pair_distances = dist_matrix[np.ix_(atoms_i_indices, atoms_j_indices)]

                        # Check if any distance is below cutoff
                        if np.any(pair_distances <= distance_cutoff):
                            contact_counts[res_i] += 1
                            contact_counts[res_j] += 1

                # Progress indicator for large trajectories
                if (frame_idx + 1) % max(1, n_frames // 10) == 0:
                    logger.info(f"Processed {frame_idx + 1}/{n_frames} frames ({100*(frame_idx+1)/n_frames:.0f}%)")

            # Normalize by number of frames
            for resid in contact_counts:
                contact_counts[resid] /= n_frames

            # Sort by frequency
            sorted_residues = sorted(contact_counts.items(), key=lambda x: x[1], reverse=True)

            # Store results
            key = label or 'contact_frequency'
            self.data[f'{key}_residues'] = [r[0] for r in sorted_residues]
            self.data[f'{key}_frequencies'] = [r[1] for r in sorted_residues]

            results = {
                'residues': [r[0] for r in sorted_residues],
                'frequencies': [r[1] for r in sorted_residues],
                'contact_counts': contact_counts
            }

            logger.info(f"Contact frequency complete: {len(residues)} residues analyzed")

            return results

        except Exception as e:
            logger.error(f"Error calculating contact frequency: {e}")
            raise

    def get_statistics(self, data_key: str) -> Dict[str, float]:
        """
        Get statistics for a data series.

        Mirrors AMBERMonitor.get_statistics() pattern.

        Args:
            data_key: Key in self.data dictionary

        Returns:
            dict: {'mean', 'std', 'min', 'max', 'median'}
        """
        if data_key not in self.data or not self.data[data_key]:
            return {}

        try:
            values = np.array(self.data[data_key])

            return {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
                'median': float(np.median(values))
            }
        except Exception as e:
            logger.error(f"Error calculating statistics for {data_key}: {e}")
            return {}

    def get_data_keys(self) -> List[str]:
        """
        Get list of all data keys currently stored.

        Returns:
            list: Data keys (excluding metadata keys ending with _mask, _masks, etc.)
        """
        metadata_suffixes = ['_mask', '_masks', '_reference', '_mask1', '_mask2',
                            '_distance_cutoff', '_angle_cutoff', '_probe_radius',
                            '_solute_mask', '_solvent_mask']

        data_keys = []
        for key in self.data.keys():
            if not any(key.endswith(suffix) for suffix in metadata_suffixes):
                data_keys.append(key)

        return sorted(data_keys)

    def export_data(self, output_format: str = 'csv',
                   output_file: Path = None) -> str:
        """
        Export analysis data to file.

        Args:
            output_format: 'csv' or 'json'
            output_file: Output file path (auto-generated if None)

        Returns:
            str: Path to exported file
        """
        import json
        import csv

        if output_file is None:
            output_file = Path(f"trajectory_analysis.{output_format}")

        data_keys = self.get_data_keys()

        if output_format == 'csv':
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)

                # Header row
                header = ['frame', 'time_ps'] + data_keys
                writer.writerow(header)

                # Data rows
                n_frames = self.system_info['n_frames']
                for i in range(n_frames):
                    row = [i, self.frame_times[i] if i < len(self.frame_times) else i * 2.0]

                    for key in data_keys:
                        if key in self.data and i < len(self.data[key]):
                            row.append(self.data[key][i])
                        else:
                            row.append('')

                    writer.writerow(row)

            logger.info(f"Exported data to CSV: {output_file}")

        elif output_format == 'json':
            export_dict = {
                'system_info': self.system_info,
                'frame_times': self.frame_times,
                'data': {}
            }

            for key in data_keys:
                export_dict['data'][key] = {
                    'values': self.data[key],
                    'statistics': self.get_statistics(key)
                }

                # Add metadata if available
                for metadata_key in [f'{key}_mask', f'{key}_masks', f'{key}_reference']:
                    if metadata_key in self.data:
                        export_dict['data'][key]['metadata'] = self.data[metadata_key]

            with open(output_file, 'w') as f:
                json.dump(export_dict, f, indent=2)

            logger.info(f"Exported data to JSON: {output_file}")

        else:
            raise ValueError(f"Unsupported format: {output_format}")

        return str(output_file)
