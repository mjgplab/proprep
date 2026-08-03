"""
Membrane Configuration Data Structures

Defines the MembraneConfig dataclass that captures all user-specified settings
for a membrane build. This is the central data object passed between UI,
runner, and workspace.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class SoluteConfig:
    """Configuration for a single solute molecule."""
    pdb_file: str
    charge: int = 0
    count: str = "1"  # integer string or e.g. "0.05M"
    in_membrane: bool = False
    prot_distance: Optional[float] = None
    frcmod: Optional[str] = None
    lib: Optional[str] = None
    use_gaff2: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pdb_file": self.pdb_file,
            "charge": self.charge,
            "count": self.count,
            "in_membrane": self.in_membrane,
            "prot_distance": self.prot_distance,
            "frcmod": self.frcmod,
            "lib": self.lib,
            "use_gaff2": self.use_gaff2,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SoluteConfig":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class MembraneConfig:
    """
    Complete configuration for a membrane build.

    Captures every setting the user can configure through the Membrane Builder UI.
    Serializable to/from dict for workspace storage and session replay.
    """

    # --- Protein ---
    protein_pdb: Optional[str] = None
    preoriented: bool = False
    orientation_method: str = "memembed"  # "memembed", "ppm3"
    skip_protonation: bool = False  # Let reduce add H for volume estimation
    n_ter_orientation: Optional[List[str]] = None  # ["in", "out"] per chain

    # --- MEMEMBED options ---
    memembed_algorithm: int = 0  # 0=GA, 1=Grid, 2=Direct, 3=GA×5
    barrel_mode: bool = False
    keep_ligands: bool = False

    # --- Lipid composition ---
    lipids: str = ""  # e.g. "DOPC:CHL1" or "DOPC//DOPE"
    ratio: str = ""  # e.g. "2:1"
    symmetric: bool = True

    # --- Solvent ---
    water_model: str = "auto"  # "auto" selects based on ffprot
    solvents: str = "WAT"
    solvent_ratio: str = ""

    # --- Ions ---
    salt: bool = True
    salt_concentration: float = 0.15
    cation: str = "K+"
    anion: str = "Cl-"
    no_counterions: bool = False
    salt_override: bool = False
    charge_delta: int = 0  # manual charge adjustment
    auto_charge_delta: int = 0  # auto-computed from redox sites

    # --- Geometry ---
    box_padding: float = 15.0
    water_layer: float = 17.5
    leaflet_width: float = 23.0
    fixed_xy: Optional[float] = None
    fixed_dims: Optional[List[float]] = None  # [X, Y, Z]
    lipid_count_method: str = "apl"  # "apl" or "volume"
    apl_offset: Optional[float] = None
    lip_offset: float = 1.0

    # --- Specialized geometry ---
    curvature: Optional[float] = None
    curvature_radius: Optional[float] = None
    gaussian_params: Optional[List[float]] = None  # [C, D, H]
    channel_plug: Optional[float] = None
    head_plane: Optional[float] = None
    tail_plane: Optional[float] = None
    self_assembly: bool = False
    solvate_only: bool = False
    pbc: bool = False
    cubic: bool = False

    # --- Double bilayer (CompEL) ---
    double_bilayer: bool = False
    charge_imbalance: int = 0
    imbalance_ion: str = "cat"  # "cat" or "an"

    # --- PACKMOL settings ---
    tolerance: float = 2.0
    protein_radius: float = 1.5
    nloop: int = 20
    nloop_all: int = 100
    gencan_iterations: int = 20
    move_fraction: float = 0.05
    move_bad_random: bool = False
    short_penalty: bool = False
    pack_all: bool = False
    random_seed: bool = False
    save_trajectory: bool = False
    plot_optimization: bool = False
    no_progress: bool = False

    # --- Force fields (for downstream tleap, informational) ---
    ffprot: str = "ff14SB"
    fflip: str = "lipid21"

    # --- Solutes ---
    solutes: List[SoluteConfig] = field(default_factory=list)

    # --- Output ---
    output_prefix: str = "bilayer"
    keep_temp_files: bool = False
    overwrite: bool = True

    @property
    def effective_water_model(self) -> str:
        """Resolve 'auto' water model based on protein force field."""
        if self.water_model != "auto":
            return self.water_model
        auto_map = {
            "ff14SB": "tip3p",
            "ff15ipq": "spceb",
            "ff19SB": "opc",
        }
        return auto_map.get(self.ffprot, "tip3p")

    @property
    def total_charge_delta(self) -> int:
        """Combined manual + auto charge correction."""
        return self.charge_delta + self.auto_charge_delta

    def get_leaprc_requirements(self) -> List[str]:
        """Determine which leaprc files the membrane system needs."""
        reqs = []
        if self.fflip == "lipid21":
            reqs.append("leaprc.lipid21")
        elif self.fflip == "lipid17":
            reqs.append("leaprc.lipid17")
        # leaprc.lipid_ext is only needed for extended lipids (those after the
        # #EXTENDED LIPIDS marker in memgen.parm). Check via LipidLibrary.
        if self.lipids:
            from .lipid_library import LipidLibrary
            lib = LipidLibrary()
            if lib.load() and lib.needs_extended(self.lipids):
                addpath = lib.get_lipid_ext_addpath()
                if addpath:
                    reqs.append(f"addPath {addpath}")
                reqs.append("leaprc.lipid_ext")
        if self.solvents != "WAT":
            reqs.append("leaprc.extra_solvents")
        return reqs

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dict for workspace storage."""
        d = {}
        for k, v in self.__dataclass_fields__.items():
            val = getattr(self, k)
            if k == "solutes":
                val = [s.to_dict() for s in val]
            d[k] = val
        # Include computed properties
        d["effective_water_model"] = self.effective_water_model
        d["total_charge_delta"] = self.total_charge_delta
        d["leaprc_requirements"] = self.get_leaprc_requirements()
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MembraneConfig":
        """Deserialize from dict."""
        # Remove computed properties that aren't constructor args
        filtered = {k: v for k, v in data.items()
                    if k in cls.__dataclass_fields__}
        if "solutes" in filtered and filtered["solutes"]:
            filtered["solutes"] = [
                SoluteConfig.from_dict(s) if isinstance(s, dict) else s
                for s in filtered["solutes"]
            ]
        return cls(**filtered)

    def to_cli_args(self) -> List[str]:
        """Convert configuration to packmol-memgen CLI arguments."""
        args = []

        # Protein. packmol-memgen accepts the literal string "None" to mean
        # "build an empty bilayer for this slot"; in that case --distxy_fix is
        # required (membrane patch size has nothing to scale around).
        if self.protein_pdb:
            args.extend(["-p", self.protein_pdb])
        else:
            args.extend(["-p", "None"])
        if self.preoriented:
            args.append("--preoriented")
        elif self.orientation_method == "ppm3":
            args.append("--ppm")

        # MEMEMBED options
        if self.memembed_algorithm != 0:
            args.extend(["--mem_opt", str(self.memembed_algorithm)])
        if self.barrel_mode:
            args.append("--barrel")
        if self.keep_ligands:
            args.append("--keepligs")

        # Protonation — reduce adds H for correct volume estimation.
        # ProPrep's Protonation State Analyzer already renamed HIS→HID/HIE/HIP;
        # reduce (without --reducebuild) respects existing residue names.
        # These flags only affect input-protein processing, so skip them for
        # empty bilayers where there's nothing to protonate or trim.
        if self.skip_protonation and self.protein_pdb:
            args.append("--notprotonate")
            args.append("--nottrim")

        # Lipids
        if self.lipids:
            args.extend(["-l", self.lipids])
        if self.ratio:
            args.extend(["-r", self.ratio])

        # N-terminus
        if self.n_ter_orientation:
            for nt in self.n_ter_orientation:
                args.extend(["--n_ter", nt])

        # Solvent
        if self.solvents != "WAT":
            args.extend(["--solvents", self.solvents])
        if self.solvent_ratio:
            args.extend(["--solvent_ratio", self.solvent_ratio])

        # Geometry
        if self.box_padding != 15.0:
            args.extend(["--dist", str(self.box_padding)])
        if self.water_layer != 17.5:
            args.extend(["--dist_wat", str(self.water_layer)])
        if self.leaflet_width != 23.0:
            args.extend(["--leaflet", str(self.leaflet_width)])
        if self.fixed_xy is not None:
            args.extend(["--distxy_fix", str(self.fixed_xy)])
        if self.fixed_dims is not None:
            args.extend(["--dims"] + [str(d) for d in self.fixed_dims])
        if self.lipid_count_method == "volume":
            args.append("--vol")
        if self.apl_offset is not None:
            args.extend(["--apl_offset", str(self.apl_offset)])
        if self.lip_offset != 1.0:
            args.extend(["--lip_offset", str(self.lip_offset)])
        if self.solvate_only:
            args.append("--solvate")
        if self.cubic:
            args.append("--cubic")
        if self.pbc:
            args.append("--pbc")

        # Specialized geometry
        if self.curvature is not None:
            args.extend(["--curvature", str(self.curvature)])
        if self.curvature_radius is not None:
            args.extend(["--curv_radius", str(self.curvature_radius)])
        if self.gaussian_params is not None:
            args.extend(["--xygauss"] + [str(p) for p in self.gaussian_params])
        if self.channel_plug is not None:
            args.extend(["--channel_plug", str(self.channel_plug)])
        if self.head_plane is not None:
            args.extend(["--headplane", str(self.head_plane)])
        if self.tail_plane is not None:
            args.extend(["--tailplane", str(self.tail_plane)])
        if self.self_assembly:
            args.append("--self_assembly")

        # Double bilayer (CompEL)
        if self.double_bilayer:
            args.append("--double")
        if self.charge_imbalance != 0:
            args.extend(["--charge_imbalance", str(self.charge_imbalance)])
        if self.imbalance_ion != "cat":
            args.extend(["--imbalance_ion", self.imbalance_ion])

        # Ions
        if self.salt:
            args.append("--salt")
            if self.salt_concentration != 0.15:
                args.extend(["--saltcon", str(self.salt_concentration)])
            # Always override: if neutralization requires more ions than the
            # specified salt concentration, packmol-memgen should proceed
            # rather than abort.  The actual concentrations are reported in
            # the log and shown to the user.
            args.append("--salt_override")
        if self.cation != "K+":
            args.extend(["--salt_c", self.cation])
        if self.anion != "Cl-":
            args.extend(["--salt_a", self.anion])
        if self.no_counterions:
            args.append("--nocounter")

        # Charge delta (manual + auto from redox sites)
        total_delta = self.total_charge_delta
        if total_delta != 0:
            args.extend(["--charge_pdb_delta", str(total_delta)])

        # PACKMOL settings
        if self.tolerance != 2.0:
            args.extend(["--tolerance", str(self.tolerance)])
        if self.protein_radius != 1.5:
            args.extend(["--prot_rad", str(self.protein_radius)])
        if self.nloop != 20:
            args.extend(["--nloop", str(self.nloop)])
        if self.nloop_all != 100:
            args.extend(["--nloop_all", str(self.nloop_all)])
        if self.gencan_iterations != 20:
            args.extend(["--maxit", str(self.gencan_iterations)])
        if self.move_fraction != 0.05:
            args.extend(["--movefrac", str(self.move_fraction)])
        if self.move_bad_random:
            args.append("--movebadrandom")
        if self.short_penalty:
            args.append("--short_penalty")
        if self.pack_all:
            args.append("--packall")
        if self.random_seed:
            args.append("--random")
        if self.save_trajectory:
            args.append("--traj")
        if self.plot_optimization:
            args.append("--plot")
        if self.no_progress:
            args.append("--noprogress")

        # Solutes
        for solute in self.solutes:
            args.extend(["--solute", solute.pdb_file])
            args.extend(["--solute_charge", str(solute.charge)])
            args.extend(["--solute_con", solute.count])
            if solute.in_membrane:
                args.append("--solute_inmem")
            if solute.prot_distance is not None:
                args.extend(["--solute_prot_dist", str(solute.prot_distance)])

        # Force fields — informational for tleap, but ffprot/fflip affect
        # packmol-memgen's water model auto-selection
        args.extend(["--ffprot", self.ffprot])
        args.extend(["--fflip", self.fflip])

        # Output
        args.extend(["-o", f"{self.output_prefix}.pdb"])
        if self.keep_temp_files:
            args.append("--keep")
        if self.overwrite:
            args.append("--overwrite")

        # NEVER parametrize or minimize — ProPrep handles these
        # (don't add --parametrize or --minimize)

        return args
