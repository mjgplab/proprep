"""
Trajectory Analysis for AMBER MD Simulations

Mirrors AMBERMonitor architecture for mdout files but operates on trajectory (.nc) files.
Provides structural analysis to complement energetic analysis.
"""

import logging
import os
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
        """Frame times from the trajectory file, or frame indices if it stores none."""
        times = getattr(self.traj, 'time', None)
        if times is not None:
            self.frame_times = [float(t) for t in times]
            self.has_frame_times = True
        else:
            # No times in the file. Plot against the frame index rather than guess
            # a time step; set_frame_interval() supplies one if the user knows it.
            self.frame_times = [float(i) for i in range(len(self.traj))]
            self.has_frame_times = False

    def set_frame_interval(self, picoseconds: float):
        """Time between frames, for a trajectory file that does not store times."""
        self.frame_times = [i * picoseconds for i in range(len(self.traj))]
        self.has_frame_times = True

    @property
    def time_label(self) -> str:
        """Axis label that matches what frame_times holds."""
        return "Time (ps)" if self.has_frame_times else "Frame"

    def _select(self, mask: str) -> np.ndarray:
        """Atom indices for an AMBER mask.

        pytraj's Topology.select() segfaults on a non-string argument instead of
        raising, which would take the whole ProPrep session down with it.
        """
        if not isinstance(mask, str):
            raise TypeError(f"AMBER mask must be a string, got {type(mask).__name__}: {mask!r}")
        return self.traj.top.select(mask)

    def _superpose(self, traj, alignment_mask: str, reference: int):
        """Fit `traj` (a copy, never self.traj) on alignment_mask to one of its frames."""
        if len(self._select(alignment_mask)) < 3:
            raise ValueError(f"Alignment mask '{alignment_mask}' selects fewer than 3 atoms; choose another")
        if not 0 <= reference < traj.n_frames:
            raise ValueError(f"Alignment reference frame {reference} outside 0-{traj.n_frames - 1}")
        pt.superpose(traj, ref=reference, mask=alignment_mask)

    def _selected_residue_ids(self, mask: str) -> List[int]:
        """Sorted 1-based ids of the residues that a mask touches."""
        top = self.traj.top
        return sorted({top.atom(int(i)).resid + 1 for i in self._select(mask)})

    @staticmethod
    def _residue_ids_to_range(residue_ids: List[int]) -> str:
        """[1, 2, 3, 7, 8] -> '1-3,7-8' (AMBER mask / cpptraj resrange syntax)."""
        ids = sorted(set(residue_ids))
        if not ids:
            return ''
        runs = []
        start = prev = ids[0]
        for rid in ids[1:]:
            if rid != prev + 1:
                runs.append((start, prev))
                start = rid
            prev = rid
        runs.append((start, prev))
        return ','.join(str(a) if a == b else f'{a}-{b}' for a, b in runs)

    def get_protein_residue_ids(self) -> List[int]:
        """1-based ids of residues that have a peptide backbone (N, CA and C atoms).

        Read from the topology, so it is exact whatever the solvent fraction is.
        """
        backbone = [set(self._selected_residue_ids(f'@{name}')) for name in ('N', 'CA', 'C')]
        return sorted(set.intersection(*backbone))

    def get_protein_residue_range(self) -> str:
        """Protein residues as a range string, e.g. '1-10'; '' if there are none."""
        return self._residue_ids_to_range(self.get_protein_residue_ids())

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

            if reference == -1:
                # Average structure. pytraj has no sentinel for it (ref=-1 is the
                # LAST frame), so build it: fit a stripped copy to its first
                # frame, average, and measure against that average.
                fitted = self.traj[mask]
                pt.superpose(fitted, ref=0)
                average = pt.mean_structure(fitted)
                rmsd_values = pt.rmsd(fitted, ref=average)
            else:
                # update_coordinate=False: pt.rmsd otherwise superposes self.traj in
                # place. The box is not rotated with it, so every later analysis that
                # images (RDF, water shells) would silently change with run order.
                rmsd_values = pt.rmsd(self.traj, mask=mask, ref=reference, update_coordinate=False)

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

    def calculate_rmsf(self, mask: str = '@CA', label: str = None,
                       alignment_mask: str = '@CA,C,N', reference: int = 0) -> np.ndarray:
        """
        Calculate per-atom/residue RMSF (root mean square fluctuation).

        Args:
            mask: AMBER atom selection mask
            label: Storage label (auto-generated if None)
            alignment_mask: Atoms the trajectory is fitted on before measuring fluctuations
            reference: Frame the trajectory is fitted to (0-based)

        Returns:
            np.ndarray: RMSF values per atom in selection (Angstroms)
        """
        try:
            logger.info(f"Calculating RMSF with mask '{mask}'")

            # RMSF requires an aligned trajectory. Fit a real copy: self.traj[:] is a
            # VIEW that shares coordinates, so superposing it rotates self.traj too.
            aligned_traj = self.traj.copy()

            # Fitting removes overall translation and rotation.
            # Note: atomicfluct/rmsf does NOT automatically align, so we must do it explicitly
            self._superpose(aligned_traj, alignment_mask, reference)

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
            self.data[f'{key}_alignment_mask'] = alignment_mask
            self.data[f'{key}_alignment_reference'] = reference
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

            if len(self._select(solute_mask)) == 0:
                raise ValueError(f"No atoms selected with solute mask: {solute_mask}")
            if len(self._select(solvent_mask)) == 0:
                raise ValueError(f"No atoms selected with solvent mask: {solvent_mask}")

            # pt.rdf returns (bin centers, g(r)) and images by default
            distances, gr_values = pt.rdf(
                self.traj,
                solvent_mask=solvent_mask,
                solute_mask=solute_mask,
                maximum=max_distance,
                bin_spacing=bin_spacing
            )

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
            sasa = pt.molsurf(self.traj, mask=mask, probe=probe_radius)

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

            # DSSP finds backbone N-H...O=C hydrogen bonds, so it needs whole
            # residues. Handed '@CA' it sees no H-bonds and calls everything coil.
            # Widen whatever was selected to the protein residues it touches.
            protein_ids = self.get_protein_residue_ids()
            touched = set(protein_ids) if not mask else set(protein_ids).intersection(self._selected_residue_ids(mask))
            if not touched:
                raise ValueError(f"No protein residues in selection: {mask or 'all'}")
            reported = sorted(touched)
            residue_mask = ':' + self._residue_ids_to_range(reported)

            # The whole protein is assigned, and the selected residues reported: a
            # residue's H-bond partner may lie outside the selection, and assigned
            # on the selection alone a strand loses its ladder.
            protein_codes = self._dssp_codes_from_cpptraj(protein_ids)
            columns = [protein_ids.index(resid) for resid in reported]
            dssp_data = protein_codes[:, columns]
            # cpptraj saw a stripped copy numbered from 1, so label from the real topology
            residue_names = {res.index + 1: res.name for res in self.traj.top.residues}
            residue_labels = [f"{residue_names[resid]}:{resid}" for resid in reported]
            helix_codes, sheet_codes, turn_codes = ('H', 'G', 'I'), ('B', 'b'), ('T',)

            # Count secondary structure types per frame
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

                n_helix = np.sum(np.isin(frame_ss, helix_codes))
                n_sheet = np.sum(np.isin(frame_ss, sheet_codes))
                n_turn = np.sum(np.isin(frame_ss, turn_codes))
                n_coil = total - n_helix - n_sheet - n_turn  # bend and unassigned

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
            self.data[f'{key}_residue_labels'] = list(residue_labels)
            self.data[f'{key}_mask'] = residue_mask

            logger.info(f"DSSP calculation complete: {len(dssp_data[0])} residues analyzed")
            return {
                'helix_pct': helix_pct,
                'sheet_pct': sheet_pct,
                'turn_pct': turn_pct,
                'coil_pct': coil_pct,
                'per_residue': dssp_data,
                'residue_labels': list(residue_labels),
                'mask': residue_mask
            }

        except Exception as e:
            logger.error(f"Error calculating DSSP: {e}")
            raise

    # cpptraj secstruct writes one integer per residue per frame
    _CPPTRAJ_DSSP_CODES = "0bBGHITS"   # 0 none, 1 parallel, 2 antiparallel, 3 3-10, 4 alpha, 5 pi, 6 turn, 7 bend

    def _dssp_codes_from_cpptraj(self, protein_ids: List[int]) -> np.ndarray:
        """DSSP codes[n_frames, n_protein_residues] from the cpptraj PROGRAM.

        Not pt.dssp: on the linux-64 pytraj shipped with AmberTools 26 it corrupts
        the heap, and a LATER pytraj call aborts the whole process ("corrupted size
        vs. prev_size", "munmap_chunk(): invalid pointer"). Measured on Rocky from
        the 1.21.0 installer: DSSP followed by other analyses aborted mid-run in 11
        of 20 sessions; with no DSSP, 0 of 20. Neither a mask-free call on a
        stripped copy (10 of 20) nor other arguments cure it. The cpptraj
        executable runs the same algorithm without pytraj's wrapper: 40 of 40 clean,
        codes identical to pt.dssp's. macOS never showed the fault.
        """
        import shutil
        import subprocess
        import sys
        import tempfile

        cpptraj = shutil.which("cpptraj") or os.path.join(os.path.dirname(sys.executable), "cpptraj")
        if not os.path.exists(cpptraj):
            raise RuntimeError("cpptraj not found (is AmberTools installed?)")

        protein_only = self.traj[':' + self._residue_ids_to_range(protein_ids)]   # a copy; solvent is of no use here
        with tempfile.TemporaryDirectory(prefix="proprep_dssp_") as work:
            top, nc, out = (os.path.join(work, name) for name in ("protein.parm7", "protein.nc", "dssp.dat"))
            pt.write_parm(top, protein_only.top, overwrite=True)
            pt.write_traj(nc, protein_only, overwrite=True)
            script = f"parm {top}\ntrajin {nc}\nsecstruct out {out}\nrun\nquit\n"
            result = subprocess.run([cpptraj], input=script, capture_output=True, text=True)
            if result.returncode != 0 or not os.path.exists(out):
                tail = "\n".join((result.stdout + result.stderr).splitlines()[-8:])
                raise RuntimeError(f"cpptraj secstruct failed (exit {result.returncode})\n{tail}")
            with open(out) as f:
                rows = [line.split()[1:] for line in f if line.strip() and not line.startswith('#')]

        numbers = np.array([[int(float(value)) for value in row] for row in rows], dtype=int)
        if numbers.shape != (self.traj.n_frames, len(protein_ids)):
            raise RuntimeError(f"cpptraj secstruct returned {numbers.shape}, expected "
                               f"({self.traj.n_frames}, {len(protein_ids)})")
        return np.array(list(self._CPPTRAJ_DSSP_CODES))[numbers]

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
                resrange = self.get_protein_residue_range()
                if not resrange:
                    raise ValueError("No protein residues (N, CA, C backbone) found in topology")
            else:
                resrange = residue_selection

            results = {}

            if dihedral_type == 'phi-psi':
                # Calculate both phi and psi angles using pytraj calc functions
                logger.info(f"Calculating phi angles for residues {resrange}")
                phi_data = pt.calc_phi(self.traj, resrange=resrange)

                logger.info(f"Calculating psi angles for residues {resrange}")
                psi_data = pt.calc_psi(self.traj, resrange=resrange)

                # Datasets are keyed 'phi:<resid>' / 'psi:<resid>'. The first residue
                # of a chain has no phi and the last has no psi, so the two lists
                # cover DIFFERENT residues; pair them by residue id, not by position.
                phi_by_res = {int(d.key.split(':')[1]): d.values for d in phi_data}
                psi_by_res = {int(d.key.split(':')[1]): d.values for d in psi_data}
                residues = sorted(set(phi_by_res) & set(psi_by_res))
                if not residues:
                    raise ValueError(f"No residues with both phi and psi in {resrange}")

                # shape (n_residues, n_frames)
                phi_array = np.array([phi_by_res[r] for r in residues])
                psi_array = np.array([psi_by_res[r] for r in residues])

                # Store results
                key = label or 'ramachandran'
                self.data[f'{key}_phi'] = phi_array.tolist()
                self.data[f'{key}_psi'] = psi_array.tolist()
                self.data[f'{key}_residues'] = residues
                self.data[f'{key}_resrange'] = resrange
                self.data[f'{key}_type'] = 'phi-psi'

                results = {
                    'phi': phi_array,
                    'psi': psi_array,
                    'residues': residues,
                    'resrange': resrange
                }

                logger.info(f"Ramachandran calculation complete: {len(residues)} residues analyzed")

            elif dihedral_type in ['chi1', 'chi2', 'chi3', 'chi4']:
                # Calculate sidechain chi angles
                logger.info(f"Calculating {dihedral_type} angles for residues {resrange}")
                # cpptraj calls the protein chi1 'chip' ('chin' is the nucleic acid chi)
                cpptraj_type = 'chip' if dihedral_type == 'chi1' else dihedral_type
                chi_data = pt.multidihedral(self.traj, dihedral_types=cpptraj_type, resrange=resrange)
                if len(chi_data) == 0:
                    raise ValueError(f"No {dihedral_type} dihedrals found in residues {resrange}")

                # Convert to numpy array
                chi_array = np.array([data.values for data in chi_data])

                # Store results
                key = label or f'dihedral_{dihedral_type}'
                self.data[f'{key}_values'] = chi_array.tolist()
                self.data[f'{key}_resrange'] = resrange
                self.data[f'{key}_type'] = dihedral_type

                results = {
                    'chi': chi_array,
                    'residues': [int(d.key.split(':')[1]) for d in chi_data],
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
                    'residues': [int(d.key.split(':')[1]) for d in omega_data],
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
                          label: str = None,
                          persistence_threshold: float = 0.5) -> dict:
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

            from scipy.spatial.distance import cdist

            atoms1 = self._select(mask1)
            atoms2 = self._select(mask2)

            n_atoms1 = len(atoms1)
            n_atoms2 = len(atoms2)
            n_frames = self.traj.n_frames

            if n_atoms1 == 0 or n_atoms2 == 0:
                raise ValueError(f"No atoms selected (mask1 '{mask1}': {n_atoms1}, mask2 '{mask2}': {n_atoms2})")
            if not 0 <= reference_frame < n_frames:
                raise ValueError(f"Reference frame {reference_frame} outside 0-{n_frames - 1}")

            logger.info(f"Selection 1: {n_atoms1} atoms, Selection 2: {n_atoms2} atoms")

            # An atom is never in contact with itself
            same_atom = atoms1[:, None] == atoms2[None, :]
            xyz = self.traj.xyz

            def contacts_in(frame_idx):
                # Direct (non-imaged) distances, as for pt.distance's default
                within = cdist(xyz[frame_idx, atoms1], xyz[frame_idx, atoms2]) <= distance_cutoff
                within[same_atom] = False
                return within

            native_contacts = contacts_in(reference_frame)
            n_native_contacts = int(np.sum(native_contacts))

            logger.info(f"Found {n_native_contacts} native contacts in reference frame {reference_frame}")

            # One frame at a time: Q-value (fraction of native contacts kept) and a
            # running occupancy count, without holding every frame's matrix
            q_values = []
            contact_counts = np.zeros((n_atoms1, n_atoms2))
            for frame_idx in range(n_frames):
                frame_contacts = contacts_in(frame_idx)
                contact_counts += frame_contacts

                if n_native_contacts > 0:
                    q_value = np.sum(native_contacts & frame_contacts) / n_native_contacts
                else:
                    q_value = 1.0  # No native contacts defined
                q_values.append(float(q_value))

            # Calculate contact frequency (occupancy) for each pair
            contact_frequency = contact_counts / n_frames

            # Persistent contacts: present in more than this fraction of frames
            persistent_contacts = []
            for i in range(n_atoms1):
                for j in range(n_atoms2):
                    occupancy = contact_frequency[i, j]
                    if occupancy > persistence_threshold:
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
                            'occupancy': float(occupancy),
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
                'mask2': mask2,
                'persistence_threshold': persistence_threshold
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
            basic_residues = ['LYS', 'ARG', 'HIS', 'HIP']

            # Define charged atoms for each residue type
            acidic_atoms = {
                'ASP': ['OD1', 'OD2'],
                'GLU': ['OE1', 'OE2']
            }
            basic_atoms = {
                'LYS': ['NZ'],
                'ARG': ['NH1', 'NH2'],
                'HIS': ['ND1', 'NE2'],
                'HIP': ['ND1', 'NE2']
            }

            # Find all acidic and basic residues in the system
            acidic_res_list = []
            basic_res_list = []

            for res in self.traj.top.residues:
                res_atom_indices = range(res.first_atom_index, res.last_atom_index)
                if res.name in acidic_residues:
                    # Get indices of charged atoms
                    charged_atoms = [i for i in res_atom_indices
                                   if self.traj.top.atom(i).name in acidic_atoms.get(res.name, [])]
                    if charged_atoms:
                        acidic_res_list.append({
                            'resid': res.index + 1,  # 1-indexed
                            'resname': res.name,
                            'atoms': charged_atoms
                        })

                elif res.name in basic_residues:
                    charged_atoms = [i for i in res_atom_indices
                                   if self.traj.top.atom(i).name in basic_atoms.get(res.name, [])]
                    if charged_atoms:
                        basic_res_list.append({
                            'resid': res.index + 1,
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
                    })

            # Closest charged-atom approach per frame, all frames at once
            # (direct, non-imaged distances, as for pt.distance's default)
            xyz = self.traj.xyz
            for pair in salt_bridge_pairs:
                acid_xyz = xyz[:, pair['acidic_res']['atoms'], :]   # (n_frames, n_acid, 3)
                base_xyz = xyz[:, pair['basic_res']['atoms'], :]    # (n_frames, n_base, 3)
                separations = np.linalg.norm(acid_xyz[:, :, None, :] - base_xyz[:, None, :, :], axis=-1)
                min_dist = separations.min(axis=(1, 2))

                pair['distances'] = min_dist.tolist()
                pair['formed'] = (min_dist <= distance_cutoff).tolist()

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
                'definitions': {**{r: acidic_atoms[r] for r in acidic_residues},
                                **{r: basic_atoms[r] for r in basic_residues}},
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
                     label: str = None,
                     fit: bool = True) -> dict:
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

            # Work on a stripped copy: with fit=True pt.pca fits its trajectory to the
            # average structure IN PLACE, which must not happen to the shared self.traj.
            selection = self.traj[mask]
            if selection.n_atoms == 0:
                raise ValueError(f"No atoms selected with mask: {mask}")
            n_components = min(n_components, 3 * selection.n_atoms)

            # Returns (projections[n_vecs, n_frames], (eigenvalues, eigenvectors))
            projections, (eigenvalues, eigenvectors) = pt.pca(selection, mask='*', n_vecs=n_components, fit=fit)

            # Variance explained is relative to ALL the motion, not just to the
            # components asked for. `selection` now holds the coordinates pt.pca
            # analysed (fitted to their average if fit=True), so their total variance
            # is the sum of every eigenvalue (the trace of the covariance matrix).
            total_variance = float(selection.xyz.var(axis=0).sum())
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
                'n_components': n_components,
                'fit': fit
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
                            label: str = None,
                            metric: str = 'rms',
                            kseed: int = 1,
                            linkage: str = 'averagelinkage') -> dict:
        """
        Calculate conformational clustering of trajectory.

        Groups similar conformations into clusters.

        Args:
            mask: AMBER mask for atoms to include (default: C-alpha atoms)
            n_clusters: Number of clusters (default: 5)
            algorithm: Clustering algorithm - 'kmeans' or 'hierarchical' (default: kmeans)
            label: Optional label for data storage
            metric: Distance between frames - 'rms' (best-fit RMSD), 'nofit' (RMSD
                without fitting) or 'dme' (distance RMSD)
            kseed: Random seed for k-means' initial points (k-means only)
            linkage: 'averagelinkage', 'linkage' (single) or 'complete' (hierarchical only)

        Returns:
            dict: Clustering results including assignments, populations, and representatives
        """
        try:
            # Default to C-alpha atoms if not specified
            if mask is None:
                mask = '@CA'

            logger.info(f"Calculating {algorithm} clustering for mask '{mask}' with {n_clusters} clusters")

            if len(self._select(mask)) == 0:
                raise ValueError(f"No atoms selected with mask: {mask}")
            if metric not in ('rms', 'nofit', 'dme'):
                raise ValueError(f"Unknown clustering metric: {metric}")
            if linkage not in ('averagelinkage', 'linkage', 'complete'):
                raise ValueError(f"Unknown linkage: {linkage}")

            # Both return a ClusteringDataset: .cluster_index (per frame),
            # .centroids (representative frame of each cluster), .population
            if algorithm.lower() == 'kmeans':
                cluster_data = pt.cluster.kmeans(
                    self.traj,
                    mask=mask,
                    n_clusters=n_clusters,
                    kseed=kseed,
                    metric=metric
                )
            elif algorithm.lower() == 'hierarchical':
                # cpptraj's name for it is hieragglo (hierarchical agglomerative)
                cluster_data = pt.cluster.hieragglo(
                    self.traj,
                    mask=mask,
                    options=f'clusters {n_clusters} {linkage} {metric}'
                )
            else:
                raise ValueError(f"Unknown clustering algorithm: {algorithm}")

            cluster_assignments = np.asarray(cluster_data.cluster_index)
            # .centroids carries cpptraj's 1-based frame numbers (.cluster_index is
            # 0-based); checked against the medoid computed from pairwise RMSDs
            representative_frames = np.asarray(cluster_data.centroids) - 1

            # cpptraj can return fewer clusters than asked for
            n_clusters = len(representative_frames)

            cluster_counts = {i: int(np.sum(cluster_assignments == i)) for i in range(n_clusters)}

            # Spread of each cluster: RMSD of its members to its representative frame
            cluster_rmsd = {}
            for cluster_id in range(n_clusters):
                cluster_frames_idx = np.where(cluster_assignments == cluster_id)[0]

                if len(cluster_frames_idx) > 1:
                    # measured with the same metric the clustering used
                    representative = self.traj[int(representative_frames[cluster_id])]
                    members = cluster_frames_idx.tolist()
                    if metric == 'dme':
                        rmsds = pt.distance_rmsd(self.traj, mask=mask, ref=representative, frame_indices=members)
                    else:
                        rmsds = pt.rmsd(self.traj, mask=mask, ref=representative, frame_indices=members,
                                        nofit=(metric == 'nofit'), update_coordinate=False)

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
                'mask': mask,
                'metric': metric,
                'kseed': kseed if algorithm.lower() == 'kmeans' else None,
                'linkage': linkage if algorithm.lower() == 'hierarchical' else None
            }

            logger.info(f"Clustering complete: {n_clusters} clusters identified")

            return results

        except Exception as e:
            logger.error(f"Error calculating clustering: {e}")
            raise

    def calculate_pairwise_rmsd(self,
                               mask: str = None,
                               subsample: int = None,
                               label: str = None,
                               metric: str = 'rms') -> dict:
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

            if len(self._select(mask)) == 0:
                raise ValueError(f"No atoms selected with mask: {mask}")
            if n_frames < 2:
                raise ValueError(f"Pairwise RMSD needs at least 2 frames, got {n_frames}")

            if metric not in ('rms', 'nofit', 'dme'):
                raise ValueError(f"Unknown pairwise metric: {metric}")

            # Full symmetric (n_frames, n_frames) matrix
            rmsd_matrix = np.asarray(pt.pairwise_rmsd(traj_to_use, mask=mask, metric=metric), dtype=float)

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
                'n_frames': n_frames,
                'metric': metric
            }

            logger.info(f"Pairwise RMSD complete: overall mean = {overall_mean:.2f} Å")

            return results

        except Exception as e:
            logger.error(f"Error calculating pairwise RMSD: {e}")
            raise

    def calculate_bfactors(self,
                          mask: str = None,
                          by_residue: bool = True,
                          label: str = None,
                          alignment_mask: str = '@CA,C,N',
                          reference: int = 0) -> dict:
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

            if len(self._select(mask)) == 0:
                raise ValueError(f"No atoms selected with mask: {mask}")

            # Fluctuations only mean something once overall tumbling is removed,
            # and atomicfluct does not fit. Fit a copy (never self.traj).
            aligned_traj = self.traj.copy()  # [:] would be a view sharing coordinates
            self._superpose(aligned_traj, alignment_mask, reference)

            # Column 0 is the atom number, column 1 the fluctuation (Angstroms)
            rmsf = pt.rmsf(aligned_traj, mask=mask)[:, 1]

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
            # Default to the protein's heavy atoms, read from the topology
            if solute_mask is None:
                protein_range = self.get_protein_residue_range()
                if not protein_range:
                    raise ValueError("No protein residues found; give a custom solute mask")
                solute_mask = f':{protein_range}&!@H='

            logger.info(f"Calculating water shells around '{solute_mask}'")

            if len(self._select(solute_mask)) == 0:
                raise ValueError(f"No atoms selected with solute mask: {solute_mask}")

            # Define shell boundaries
            if shell_width <= 0:
                raise ValueError(f"Shell width must be positive, got {shell_width}")
            n_shells = int(max_distance / shell_width)
            if n_shells < 1:
                raise ValueError(f"Maximum distance {max_distance} is smaller than the shell width {shell_width}")
            shell_boundaries = [(i * shell_width, (i + 1) * shell_width) for i in range(n_shells)]

            # Water oxygens, under whichever water residue name the topology uses
            water_resnames = sorted({res.name for res in self.traj.top.residues}
                                    & {'WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPC', 'OPC'})
            if not water_resnames:
                raise ValueError("No water residues found in topology")
            water_mask = f":{','.join(water_resnames)}@O,OW,OH2"
            n_waters = len(self._select(water_mask))
            if n_waters == 0:
                raise ValueError(f"No water oxygens selected with mask: {water_mask}")

            logger.info(f"Found {len(self._select(solute_mask))} solute atoms and {n_waters} water molecules")

            n_frames = self.traj.n_frames

            # cpptraj's watershell counts, per frame and with periodic imaging, the
            # waters whose oxygen lies within a distance of any solute atom. It takes
            # two distances per call; a shell's population is the difference between
            # the counts within its outer and its inner edge.
            edges = [upper for _, upper in shell_boundaries]
            within = {}
            for i in range(0, len(edges), 2):
                if i + 1 < len(edges):
                    lower, upper = edges[i], edges[i + 1]
                elif i > 0:
                    lower, upper = edges[i - 1], edges[i]  # odd number of edges: pair the last with its neighbour
                else:
                    lower, upper = edges[0] / 2, edges[0]  # a single shell
                counts = pt.watershell(self.traj, solute_mask=solute_mask, solvent_mask=water_mask,
                                       lower=lower, upper=upper, dtype='ndarray')
                within[lower], within[upper] = np.asarray(counts[0]), np.asarray(counts[1])

            shell_populations = np.zeros((n_frames, n_shells))
            inner = np.zeros(n_frames)
            for shell_idx, edge in enumerate(edges):
                shell_populations[:, shell_idx] = within[edge] - inner
                inner = within[edge]

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
                'water_mask': water_mask,
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
                             label: str = None,
                             fit_mask: str = None,
                             max_grid_points: int = 100,
                             n_peaks: int = 5) -> dict:
        """
        Calculate 3D density map for selected atoms.

        Creates spatial density distribution showing where atoms
        spend most time during simulation.

        Args:
            selection_mask: AMBER mask for atoms to analyze
            grid_spacing: Grid spacing in Angstroms (default: 1.0)
            label: Optional label for data storage
            fit_mask: If given, density is accumulated in the frame of these atoms:
                a copy of the trajectory is autoimaged and fitted on them to its
                first frame. If None, coordinates are used as stored (lab frame).
            max_grid_points: Most grid points allowed along an axis; a finer grid
                than that has its spacing coarsened (reported in the result)
            n_peaks: Number of highest-density cells to report

        Returns:
            dict: Density map results including grid and peak locations
        """
        try:
            logger.info(f"Calculating density map for '{selection_mask}' with {grid_spacing} Å spacing")

            # Get selected atoms
            atoms = self._select(selection_mask)

            if len(atoms) == 0:
                raise ValueError(f"No atoms selected with mask: {selection_mask}")

            logger.info(f"Selected {len(atoms)} atoms")

            if grid_spacing <= 0:
                raise ValueError(f"Grid spacing must be positive, got {grid_spacing}")
            requested_spacing = grid_spacing

            if fit_mask:
                # Density around a solute only means something in the solute's own
                # frame: make molecules whole around it, then remove its tumbling.
                # On a copy; the shared trajectory is never moved.
                source = self.traj.copy()
                pt.autoimage(source)
                self._superpose(source, fit_mask, 0)
            else:
                source = self.traj

            # Grid over where the selected atoms actually go. (The box is no guide:
            # it need not start at the origin and may not be orthorhombic.)
            coords = source.xyz[:, atoms, :].reshape(-1, 3)  # every frame, every selected atom
            grid_min = coords.min(axis=0)
            grid_max = coords.max(axis=0)

            # Create 3D grid
            nx = int((grid_max[0] - grid_min[0]) / grid_spacing) + 1
            ny = int((grid_max[1] - grid_min[1]) / grid_spacing) + 1
            nz = int((grid_max[2] - grid_min[2]) / grid_spacing) + 1

            # Limit grid size for memory
            max_grid_dim = max_grid_points
            if nx > max_grid_dim or ny > max_grid_dim or nz > max_grid_dim:
                logger.warning(f"Grid too large ({nx}x{ny}x{nz}), adjusting spacing")
                # n points span n - 1 intervals
                new_spacing = max((grid_max[0] - grid_min[0]) / (max_grid_dim - 1),
                                (grid_max[1] - grid_min[1]) / (max_grid_dim - 1),
                                (grid_max[2] - grid_min[2]) / (max_grid_dim - 1))
                grid_spacing = new_spacing
                nx = int((grid_max[0] - grid_min[0]) / grid_spacing) + 1
                ny = int((grid_max[1] - grid_min[1]) / grid_spacing) + 1
                nz = int((grid_max[2] - grid_min[2]) / grid_spacing) + 1

            logger.info(f"Grid dimensions: {nx}x{ny}x{nz}")

            # Occupancy of each grid cell, as a fraction of all (frame, atom) positions
            edges = [grid_min[axis] + grid_spacing * np.arange(n + 1) for axis, n in enumerate((nx, ny, nz))]
            density, _ = np.histogramdd(coords, bins=edges)
            density = density / (self.traj.n_frames * len(atoms))

            # Find peak density locations
            flat_density = density.flatten()
            sorted_indices = np.argsort(flat_density)[::-1]

            # Get top 5 peak locations
            peaks = []
            for idx in sorted_indices[:n_peaks]:
                if flat_density[idx] > 0:
                    # Convert back to 3D indices
                    ix = idx // (ny * nz)
                    iy = (idx % (ny * nz)) // nz
                    iz = idx % nz

                    # Convert to real coordinates
                    # centre of the grid cell
                    x = grid_min[0] + (ix + 0.5) * grid_spacing
                    y = grid_min[1] + (iy + 0.5) * grid_spacing
                    z = grid_min[2] + (iz + 0.5) * grid_spacing

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
                'selection': selection_mask,
                'fit_mask': fit_mask,
                'requested_spacing': requested_spacing,
                'slice_z': float(grid_min[2] + (mid_z + 0.5) * grid_spacing)
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
            atoms1 = self._select(atom1_mask)
            atoms2 = self._select(atom2_mask)

            if len(atoms1) == 0 or len(atoms2) == 0:
                raise ValueError(f"No atoms selected (start '{atom1_mask}': {len(atoms1)}, end '{atom2_mask}': {len(atoms2)})")

            # Geometric centre of each selection in every frame, and the vector between them
            xyz = self.traj.xyz
            vectors = xyz[:, atoms2, :].mean(axis=1) - xyz[:, atoms1, :].mean(axis=1)  # (n_frames, 3)
            magnitudes = np.linalg.norm(vectors, axis=1)

            # Polar angle from the Z axis (0 for a zero-length vector), and
            # azimuthal angle in the XY plane from the X axis
            safe = np.where(magnitudes > 0, magnitudes, 1.0)
            polar_angles = np.where(magnitudes > 0,
                                    np.degrees(np.arccos(np.clip(vectors[:, 2] / safe, -1.0, 1.0))), 0.0)
            azimuthal_angles = np.degrees(np.arctan2(vectors[:, 1], vectors[:, 0]))

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

            atoms = self._select(mask)
            if len(atoms) == 0:
                raise ValueError(f"No atoms selected with mask: {mask}")
            n_frames = self.traj.n_frames

            # Residue (1-indexed) of each selected atom
            atom_resids = np.array([self.traj.top.atom(int(i)).resid + 1 for i in atoms])
            residues = sorted(set(atom_resids.tolist()))
            n_residues = len(residues)

            logger.info(f"Analyzing {n_residues} residues across {n_frames} frames")

            # membership[r, a] = 1 if selected atom a belongs to residue r
            membership = (np.array(residues)[:, None] == atom_resids[None, :]).astype(float)

            # Two residues are in contact in a frame if ANY of their selected atoms
            # are within the cutoff (direct, non-imaged distances)
            partner_counts = np.zeros(n_residues)
            xyz = self.traj.xyz
            for frame_idx in range(n_frames):
                selected_coords = xyz[frame_idx, atoms]
                atom_contacts = (cdist(selected_coords, selected_coords) <= distance_cutoff).astype(float)
                residue_contacts = (membership @ atom_contacts @ membership.T) > 0
                np.fill_diagonal(residue_contacts, False)  # not with itself
                partner_counts += residue_contacts.sum(axis=1)

            contact_counts = {resid: float(partner_counts[i]) for i, resid in enumerate(residues)}

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
                    row = [i, self.frame_times[i] if i < len(self.frame_times) else i]

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
