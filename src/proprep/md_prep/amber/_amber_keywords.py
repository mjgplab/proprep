"""
Comprehensive AMBER Keyword Database
=====================================
Faithfully parsed from the AMBER Manual:
  - Chapters 22 (sander) and 23 (pmemd): General MD keywords
  - Chapter 4: Generalized Born / Surface Area (GB/SA) Model
  - Chapter 6: Poisson-Boltzmann Surface Area (PBSA) Model
  - Chapter 11: QM/MM Calculations (including Chapters 10-11)
  - Chapter 27: Discrete Constant pH MD (CpHMD)
  - Chapter 28: Constant Redox Potential MD
  - Chapter 29: Continuous Constant pH MD (lambda-dynamics CpHMD)
All descriptions, options, defaults, and notes preserve the manual's exact wording.

Structure:
  - Keyword dataclass: For standard namelist keyword=value entries
    (&cntrl, &ewald, &pol_gauss, &debugf, &qmmm, &pb,
     &adf, &gms, &gau, &orc, &qc, &mrcc, &fb, &quick, &tc,
     &xtb, &dftbplus, &dprc, &sebomd, &vsolv, &adqmmm,
     &phmdin, &phmdparm, &phmdstrt)
  - FileRedirect dataclass: Extended with CpHMD file specifications (cpin, cpout, cprestrt, cein, ceout, cerestrt)
  - WtType dataclass: For &wt varying-conditions TYPE values (Section 22.9)
  - WtParameter dataclass: For &wt shared parameters (ISTEP1, ISTEP2, etc.)
  - FileRedirect dataclass: For file redirection commands (Section 22.10)
"""

from dataclasses import dataclass, field
from typing import Any, Optional, Dict, List

@dataclass
class Keyword:
    name: str
    description: str
    default: Any
    options: Optional[Dict[Any, str]] = None
    notes: Optional[str] = None
    value_type: str = "int"
    min_val: Optional[float] = None
    max_val: Optional[float] = None
    section: str = "cntrl"
    category: str = ""
    related: List[str] = field(default_factory=list)
    commonly_changed: bool = False


# =============================================================================
# &cntrl NAMELIST
# =============================================================================

# ---------------------------------------------------------------------------
# 22.6.1 General flags describing the calculation
# ---------------------------------------------------------------------------

KEYWORDS: List[Keyword] = []

KEYWORDS.append(Keyword(
    name="imin",
    description=(
        "Flag to run minimization."
    ),
    default=0,
    options={
        0: "Run molecular dynamics without any minimization.",
        1: "Perform an energy minimization.",
        5: (
            "Read in a trajectory for analysis using the minimization algorithms. "
            "If imin is set to 5, sander will read a trajectory file (the 'inptraj' argument, "
            "specified using -y on the command line), and will perform the functions described in "
            "the mdin file (e.g., an energy minimization) for each of the structures in this file. "
            "The final structure from each minimization will be written out to the normal mdcrd file. "
            "If you wish to read in a binary (i.e., NetCDF format) trajectory, be sure to set ioutfm "
            "to 1. Note that this will result in the output trajectory having NetCDF format as well. "
            "For example, when imin = 5 and maxcyc = 1000, sander will minimize each structure in the "
            "trajectory for 1000 steps and write a minimized coordinate set for each frame to the mdcrd "
            "file. If maxcyc = 1, the output file can be used to extract the energies of each of the "
            "coordinate sets in the inptraj file. Trajectories containing box coordinates can be "
            "post-processed. In order to read trajectories with box coordinates, ntb should be greater "
            "than 0."
        ),
        6: (
            "Read in a trajectory for analysis using the molecular dynamics driver. "
            "Like imin=5, this option reads a trajectory file for analysis (the 'inptraj' argument, "
            "specified using -y on the command line). Instead of minimizing the potential energy of each "
            "coordinate set, it instead initiates dynamics from each frame as if it were read as a restart "
            "file. If nstlim=0, then this effectively performs a single point energy for each frame. "
            "If the input trajectory file contains velocities (from a previous simulation that set "
            "ntwv=-1), then the stored values are used as the initial velocities of each simulation. "
            "If the input trajectory file does not contain velocities, then each simulation assigns "
            "initial velocities based on the value of tempi. If the random number seed (ig) is set to -1, "
            "then a new random seed value is generated for each frame; otherwise all frames use the same "
            "random number seed."
        ),
        7: (
            "Listen to the selected internet socket and return energies and forces when instructed by "
            "an external server. When this option is set, sander does not perform MD; instead, it listens "
            "for messages from a server instructing it to compute the potential energy and forces of a "
            "system. The server IP address and port number are provided as command line arguments -host "
            "and -port. The default values are -host 127.0.0.1 and -port 31415. The communication pattern "
            "follows the protocol implemented in the i-PI software."
        ),
    },
    notes=(
        "IMPORTANT CAVEAT (for imin=5): The initial coordinates input file used (-c <inpcrd>) should be "
        "the same as the initial coordinates input file used to generate the original trajectory. This is "
        "because sander sets up parameters for PME from the box coordinates in the initial coordinates "
        "input file."
    ),
    value_type="int",
    section="cntrl",
    category="General flags describing the calculation",
    related=["maxcyc", "nstlim", "ntx", "ioutfm"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="nmropt",
    description=(
        "NMR restraint and weight change flag."
    ),
    default=0,
    options={
        0: "No nmr-type analysis will be done.",
        1: "NMR restraints and weight changes will be read.",
        2: (
            "NMR restraints, weight changes, NOESY volumes, chemical shifts and residual "
            "dipolar restraints will be read."
        ),
    },
    value_type="int",
    section="cntrl",
    category="General flags describing the calculation",
    related=["iscale", "noeskp", "ipnlty"],
    commonly_changed=True,
))

# ---------------------------------------------------------------------------
# 22.6.2 Nature and format of the input
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ntx",
    description=(
        "Option to read the initial coordinates, velocities, and box size from the inpcrd file. "
        "Option 1 must be used when one is starting from minimized or model-built coordinates. "
        "If an MD restrt file is specified for inpcrd then option 5 is generally used (unless you "
        "explicitly wish to ignore the velocities that are present)."
    ),
    default=1,
    options={
        1: (
            "Coordinates, but no velocities, will be read; either formatted (ASCII) files or "
            "NetCDF files can be used, as the input file type will be auto-detected."
        ),
        5: (
            "Coordinates and velocities will be read from either a NetCDF or a formatted (ASCII) "
            "coordinate file. Box information will be read if ntb > 0. The velocity information "
            "will only be used if irest = 1."
        ),
    },
    value_type="int",
    section="cntrl",
    category="Nature and format of the input",
    related=["irest", "ntb"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="irest",
    description=(
        "Flag to restart a simulation."
    ),
    default=0,
    options={
        0: (
            "Do not restart the simulation; instead, run as a new simulation. Velocities in the "
            "input coordinate file, if any, will be ignored, and the time step count will be set to 0 "
            "(unless overridden by t)."
        ),
        1: (
            "Restart the simulation, reading coordinates and velocities from a previously saved restart "
            "file. The velocity information is necessary when restarting, so ntx must be 5 if irest = 1."
        ),
    },
    value_type="int",
    section="cntrl",
    category="Nature and format of the input",
    related=["ntx", "t"],
    commonly_changed=True,
))

# ---------------------------------------------------------------------------
# 22.6.3 Nature and format of the output
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ntxo",
    description=(
        "Format of the final coordinates, velocities, and box size (if a constant volume or pressure "
        "run) written to file 'restrt'."
    ),
    default=2,
    options={
        1: "Formatted (ASCII)",
        2: "NetCDF file (recommended, unless you have a workflow that requires the formatted form.)",
    },
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ioutfm"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description=(
        "Every ntpr steps, energy information will be printed in human-readable form to files 'mdout' "
        "and 'mdinfo'. 'mdinfo' is closed and reopened each time, so it always contains the most recent "
        "energy and temperature."
    ),
    default=50,
    value_type="int",
    min_val=1,
    section="cntrl",
    category="Nature and format of the output",
    related=["ntwe", "ntave"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntave",
    description=(
        "Every ntave steps of dynamics, running averages of average energies and fluctuations over the "
        "last ntave steps will be printed out. A value of 0 disables this printout. Setting ntave to a "
        "value 1/2 or 1/4 of nstlim provides a simple way to look at convergence during the simulation."
    ),
    default=0,
    notes=(
        "Avoid setting ntave != 0 on GPU runs. Turning on the printing of running averages results in "
        "the code needing to calculate both energy and forces on every step. This can lead to performance "
        "losses of 20% or more when running on the GPU."
    ),
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Nature and format of the output",
    related=["ntpr", "nstlim"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwr",
    description=(
        "Every ntwr steps during dynamics, the 'restrt' file will be written, ensuring that recovery "
        "from a crash will not be so painful. No matter what the value of ntwr, a restrt file will be "
        "written at the end of the run, i.e., after nstlim steps (for dynamics) or maxcyc steps (for "
        "minimization). If ntwr < 0, a unique copy of the file, 'restrt_<nstep>', is written every "
        "abs(ntwr) steps. This option is useful if for example one wants to run free energy perturbations "
        "from multiple starting points or save a series of restrt files for minimization."
    ),
    default="nstlim",
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["nstlim", "maxcyc"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="iwrap",
    description=(
        "If iwrap = 1, the coordinates written to the restart and trajectory files will be 'wrapped' "
        "into a primary box. This means that for each molecule, its periodic image closest to the middle "
        "of the 'primary box' (with x coordinates between 0 and a, y coordinates between 0 and b, and z "
        "coordinates between 0 and c) will be the one written to the output file. This often makes the "
        "resulting structures look better visually, but has no effect on the energy or forces. Performing "
        "such wrapping, however, can mess up diffusion and other calculations. If iwrap = 0, no wrapping "
        "will be performed, in which case it is typical to use cpptraj as a post-processing program to "
        "translate molecules back to the primary box. For very long runs, setting iwrap = 1 may be required "
        "to keep the coordinate output from overflowing the trajectory and restart file formats, especially "
        "if trajectories are written in ASCII format instead of NetCDF."
    ),
    default=0,
    options={
        0: "No wrapping will be performed.",
        1: "Coordinates will be wrapped into a primary box.",
    },
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ioutfm", "ntwx"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntwx",
    description=(
        "Every ntwx steps, the coordinates will be written to the mdcrd file. If ntwx = 0, no coordinate "
        "trajectory file will be written."
    ),
    default=0,
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Nature and format of the output",
    related=["ioutfm", "ntwprt", "ntwv", "ntwf"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntwv",
    description=(
        "Every ntwv steps, the velocities will be written to the mdvel file. If ntwv = 0, no velocity "
        "trajectory file will be written. If ntwv = -1, velocities will be written to mdcrd, which then "
        "becomes a combined coordinate/velocity trajectory file, at the interval defined by ntwx. This "
        "option is available only for binary NetCDF output (ioutfm = 1). Most users will have no need "
        "for a velocity trajectory file and so can safely leave ntwv at the default."
    ),
    default=0,
    notes=(
        "Dumping velocities frequently, like forces or coordinates, will introduce potentially "
        "significant I/O and communication overhead, hurting both performance and parallel scaling."
    ),
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ntwx", "ioutfm"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ionstepvelocities",
    description=(
        "Controls whether to print the half-step-ahead velocities (0, default) or on-step velocities (1). "
        "The half-step-ahead velocities can potentially be used to restart calculations, but the on-step "
        "velocities correspond to calculated kinetic energy/temperature."
    ),
    default=0,
    options={
        0: "Print half-step-ahead velocities.",
        1: "Print on-step velocities.",
    },
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ntwv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwf",
    description=(
        "Every ntwf steps, the forces will be written to the mdfrc file. If ntwf = 0, no force trajectory "
        "file will be written. If ntwf = -1, forces will be written to the mdcrd, which then becomes a "
        "combined coordinate/force trajectory file, at the interval defined by ntwx. This option is "
        "available only for binary NetCDF output (ioutfm = 1). Most users will have no need for a force "
        "trajectory file and so can safely leave ntwf at the default."
    ),
    default=0,
    notes=(
        "Dumping forces frequently, like velocities or coordinates, will introduce potentially significant "
        "I/O and communication overhead, hurting both performance and parallel scaling."
    ),
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ntwx", "ioutfm"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwe",
    description=(
        "Every ntwe steps, the energies and temperatures will be written to file 'mden' in a compact form. "
        "If ntwe = 0 then no mden file will be written."
    ),
    default=0,
    notes=(
        "Energies in the mden file are not synchronized with coordinates or velocities in the mdcrd or "
        "mdvel file(s). Assuming identical ntwe and ntwx values the energies are one time step before the "
        "coordinates (as well as the velocities which are synchronized with the coordinates). Consequently, "
        "an mden file is rarely written."
    ),
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Nature and format of the output",
    related=["ntpr", "ntwx"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ioutfm",
    description=(
        "The format of coordinate and velocity trajectory files (mdcrd, mdvel and inptraj). As of Amber 9, "
        "the binary format used in previous versions is no longer supported; binary output is now in NetCDF "
        "trajectory format. Binary trajectory files have many advantages: they are smaller, higher precision, "
        "much faster to read and write, and able to accept a wider range of coordinate (or velocity) values "
        "than formatted trajectory files."
    ),
    default=1,
    options={
        0: "Formatted ASCII trajectory",
        1: "Binary NetCDF trajectory",
    },
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ntwx", "ntwv", "iwrap"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwprt",
    description=(
        "The number of atoms to include in trajectory files (mdcrd and mdvel). This flag can be used to "
        "decrease the size of these files, by including only the first part of the system, which is usually "
        "of greater interest (for instance, one might include only the solute and not the solvent). "
        "If ntwprt = 0, all atoms will be included."
    ),
    default=0,
    options={
        0: "Include all atoms of the system when writing trajectories.",
    },
    notes="When > 0, include only atoms 1 to ntwprt when writing trajectories.",
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Nature and format of the output",
    related=["ntwx", "ntwv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="idecomp",
    description=(
        "Perform energy decomposition according to a chosen scheme. In former distributions this option was "
        "only really useful in conjunction with mm_pbsa, where it is turned on automatically if required. Now, "
        "a decomposition of dV/dlambda on a per-residue basis in thermodynamic integration (TI) simulations "
        "is also possible."
    ),
    default=0,
    options={
        0: "Do not decompose energies.",
        1: "Decompose energies on a per-residue basis; 1-4 EEL + 1-4 VDW are added to internal (bond, angle, dihedral) energies.",
        2: "Decompose energies on a per-residue basis; 1-4 EEL + 1-4 VDW are added to EEL and VDW.",
        3: "Decompose energies on a pairwise per-residue basis; otherwise equivalent to idecomp = 1. Not available in TI simulations.",
        4: "Decompose energies on a pairwise per-residue basis; otherwise equivalent to idecomp = 2. Not available in TI simulations.",
    },
    notes="Use of idecomp > 0 is incompatible with ntr > 0 or ibelly > 0.",
    value_type="int",
    section="cntrl",
    category="Nature and format of the output",
    related=["ntr", "ibelly"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.4 Frozen or restrained atoms
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ibelly",
    description=(
        "Flag for belly type dynamics. If set to 1, a subset of the atoms in the system will be allowed "
        "to move, and the coordinates of the rest will be frozen. The moving atoms are specified with "
        "bellymask. This option is not available when igb>0. When belly type dynamics is in use, bonded "
        "energy terms, vdW interactions, and direct space electrostatic interactions are not calculated for "
        "pairs of frozen atoms. Note that this does not provide any significant speed advantage. Freezing "
        "atoms can be useful for some applications but is maintained primarily for backwards compatibility "
        "with older versions of Amber. Most applications should use the ntr variable instead to restrain "
        "parts of the system to stay close to some initial configuration."
    ),
    default=0,
    options={
        0: "No belly type dynamics.",
        1: "Enable belly type dynamics; specify moving atoms with bellymask.",
    },
    value_type="int",
    section="cntrl",
    category="Frozen or restrained atoms",
    related=["bellymask", "ntr"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntr",
    description=(
        "Flag for restraining specified atoms in Cartesian space using a harmonic potential, if ntr = 1. "
        "The restrained atoms are determined by the restraintmask string. The force constant is restraint_wt. "
        "The reference coordinates are read in 'restrt' format from the 'refc' file."
    ),
    default=0,
    options={
        0: "No Cartesian restraints.",
        1: "Apply Cartesian restraints to atoms specified by restraintmask.",
    },
    notes="If ntr=1, you should also set netfrc=0; see the netfrc variable for more information.",
    value_type="int",
    section="cntrl",
    category="Frozen or restrained atoms",
    related=["restraint_wt", "restraintmask", "netfrc"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="restraint_wt",
    description=(
        "The weight (kcal/mol/A^2) for Cartesian restraints when ntr = 1. The restraint is of the form "
        "k*(dx)^2, where k is the force constant of value given by this variable, and dx is the difference "
        "between one of the Cartesian coordinates of a restrained atom and its reference position. There "
        "is a term like this for each Cartesian coordinate of each restrained atom."
    ),
    default=0.0,
    notes=(
        "This variable does not have anything to do with NMR restraints, and there is no way to have "
        "restraint_wt depend upon the time step."
    ),
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Frozen or restrained atoms",
    related=["ntr", "restraintmask"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="restraintmask",
    description=(
        "String that specifies the restrained atoms when ntr = 1. The syntax is given in Section 24.1.1. "
        "Note that these mask strings are limited to a maximum of 256 characters."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Frozen or restrained atoms",
    related=["ntr", "restraint_wt"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="bellymask",
    description=(
        "String that specifies the moving atoms when ibelly=1. The syntax is given in Section 24.1.1. "
        "Note that these mask strings are limited to a maximum of 256 characters."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Frozen or restrained atoms",
    related=["ibelly"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.5 Energy minimization
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="maxcyc",
    description="The maximum number of cycles of minimization.",
    default=1,
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Energy minimization",
    related=["imin", "ncyc", "ntmin", "drms"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ncyc",
    description=(
        "If NTMIN is 1 then the method of minimization will be switched from steepest descent to "
        "conjugate gradient after NCYC cycles."
    ),
    default=10,
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Energy minimization",
    related=["ntmin", "maxcyc"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntmin",
    description="Flag for the method of minimization.",
    default=1,
    options={
        0: "Full conjugate gradient minimization. The first 4 cycles are steepest descent at the start of the run and after every nonbonded pairlist update.",
        1: "For NCYC cycles the steepest descent method is used then conjugate gradient is switched on (default).",
        2: "Only the steepest descent method is used.",
        3: "The XMIN method is used.",
        4: "The LMOD method is used.",
        5: "The DL-Find module is used.",
    },
    value_type="int",
    section="cntrl",
    category="Energy minimization",
    related=["ncyc", "maxcyc", "dx0", "drms"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="dx0",
    description=(
        "The initial step length. If the initial step length is too big then will give a huge energy; "
        "however the minimizer is smart enough to adjust itself."
    ),
    default=0.01,
    value_type="float",
    section="cntrl",
    category="Energy minimization",
    related=["ntmin", "maxcyc"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="drms",
    description=(
        "The convergence criterion for the energy gradient: minimization will halt when the Root-Mean-"
        "Square of the Cartesian elements of the gradient of the energy is less than this."
    ),
    default=1e-4,
    notes="Units: kcal/mol/A.",
    value_type="float",
    section="cntrl",
    category="Energy minimization",
    related=["maxcyc", "ntmin"],
    commonly_changed=True,
))

# ---------------------------------------------------------------------------
# 22.6.6 Molecular dynamics
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="nstlim",
    description="Number of MD-steps to be performed.",
    default=1,
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Molecular dynamics",
    related=["dt", "ntwr", "imin"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="nscm",
    description=(
        "Flag for the removal of translational and rotational center-of-mass (COM) motion at regular "
        "intervals. For non-periodic simulations, after every NSCM steps, translational and rotational "
        "motion will be removed. For periodic systems, just the translational center-of-mass motion will "
        "be removed. This flag is ignored for belly simulations."
    ),
    default=1000,
    notes=(
        "For Langevin dynamics, the position of the center-of-mass of the molecule is reset to zero every "
        "NSCM steps, but the velocities are not affected. Hence there is no change to either the translation "
        "or rotational components of the momenta."
    ),
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Molecular dynamics",
    related=["ntt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="t",
    description=(
        "The time at the start (psec) this is for your own reference and is not critical. Start time "
        "is taken from the coordinate input file if IREST=1."
    ),
    default=0.0,
    value_type="float",
    section="cntrl",
    category="Molecular dynamics",
    related=["irest", "dt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dt",
    description=(
        "The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if it isn't. "
        "Note that for temperatures above 300K, the step size should be reduced since greater temperatures "
        "mean increased velocities and longer distance traveled between each force evaluation, which can "
        "lead to anomalously high energies and system blowup."
    ),
    default=0.001,
    notes=(
        "The use of Hydrogen Mass Repartitioning (HMR), together with SHAKE, allows the time step to "
        "be increased in a stable fashion by about a factor of two (up to .004) by slowing down the high "
        "frequency hydrogen motion in the system. To use HMR, the masses in the topology file need to be "
        "altered before starting the simulation. ParmEd can do this automatically with the HMassRepartition "
        "option."
    ),
    value_type="float",
    min_val=0.0,
    max_val=0.004,
    section="cntrl",
    category="Molecular dynamics",
    related=["ntc", "nstlim"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="nrespa",
    description=(
        "This variable allows the user to evaluate slowly-varying terms in the force field less frequently. "
        "For PME, 'slowly-varying' (now) means the reciprocal sum. For generalized Born runs, the "
        "'slowly-varying' forces are those involving derivatives with respect to the effective radii, and "
        "pair interactions whose distances are greater than the 'inner' cutoff, currently hard-wired at 8 A. "
        "If NRESPA>1 these slowly-varying forces are evaluated every nrespa steps. The forces are adjusted "
        "appropriately, leading to an impulse at that step. If nrespa*dt is less than or equal to 4 fs "
        "then the energy conservation is not seriously compromised. However if nrespa*dt > 4 fs then the "
        "simulation becomes less stable."
    ),
    default=1,
    notes=(
        "Energies and related quantities are only accessible every nrespa steps, since the values at "
        "other times are meaningless. The GPU PME code requires nrespa to be 1."
    ),
    value_type="int",
    min_val=1,
    section="cntrl",
    category="Molecular dynamics",
    related=["dt"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.7 Temperature regulation
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ntt",
    description=(
        "Switch for temperature scaling. Note that setting ntt=0 corresponds to the microcanonical (NVE) "
        "ensemble (which should approach the canonical one for large numbers of degrees of freedom). "
        "Some aspects of the 'weak-coupling ensemble' (ntt=1) have been examined, and roughly interpolate "
        "between the microcanonical and canonical ensembles. The ntt=2 and 3 options correspond to the "
        "canonical (constant T) ensemble."
    ),
    default=0,
    options={
        0: "Constant total energy classical dynamics (assuming that ntb<2, as should probably always be the case when ntt=0).",
        1: (
            "Constant temperature, using the weak-coupling algorithm. A single scaling factor is used for "
            "all atoms. Note that this algorithm just ensures that the total kinetic energy is appropriate "
            "for the desired temperature; it does nothing to ensure that the temperature is even over all "
            "parts of the molecule."
        ),
        2: (
            "Andersen-like temperature coupling scheme, in which imaginary 'collisions' randomize the "
            "velocities to a distribution corresponding to temp0 every vrand steps. Note that in between "
            "these 'massive collisions', the dynamics is Newtonian."
        ),
        3: (
            "Use Langevin dynamics with the collision frequency given by gamma_ln. Note that when gamma_ln "
            "has its default value of zero, this is the same as setting ntt = 0. Since Langevin simulations "
            "are highly susceptible to 'synchronization' artifacts, you should explicitly set the ig variable "
            "to a different value at each restart of a given simulation."
        ),
        9: (
            "Optimized Isokinetic Nose-Hoover chain ensemble (OIN). Constant temperature simulation "
            "utilizing Nose-Hoover chains and an isokinetic constraint on the particle and thermostat "
            "velocities, implemented for use in multiple time-stepping methods, namely for 3D-RISM and RESPA."
        ),
        10: (
            "Stochastic Isokinetic Nose-Hoover RESPA integrator (SINR). A novel isokinetic integrator "
            "that invokes an isokinetic constraint on the particle velocities combined with auxiliary "
            "thermostat velocities v1 and v2."
        ),
        11: (
            "Stochastic version of Berendsen thermostat, also known as Bussi thermostat. This thermostat "
            "samples canonical distribution by scaling all velocities to a random temperature probed from "
            "canonical distribution. Collision frequency with thermostat is controlled by tautp."
        ),
    },
    notes=(
        "Flag 'ntt' is used for the temperature regulation in the default thermostat scheme. The 'middle' "
        "thermostat scheme (Section 22.6.10) is much more efficient than the default scheme to accurately "
        "sample the configuration/conformation space in the molecular dynamics simulation for the NVT "
        "ensemble. Using ntt=1 is especially dangerous for generalized Born simulations, where there are "
        "no collisions with solvent to aid in thermalization. Other temperature coupling options "
        "(especially ntt=3) should be used instead."
    ),
    value_type="int",
    section="cntrl",
    category="Temperature regulation",
    related=["temp0", "gamma_ln", "tautp", "ig", "vrand"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="temp0",
    description=(
        "Reference temperature at which the system is to be kept, if ntt > 0. Note that for temperatures "
        "above 300K, the step size should be reduced since increased distance traveled between evaluations "
        "can lead to SHAKE and other problems."
    ),
    default=300.0,
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "tempi", "dt"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="temp0les",
    description=(
        "This is the target temperature for all LES particles. If temp0les<0, a single temperature bath "
        "is used for all atoms, otherwise separate thermostats are used for LES and non-LES particles."
    ),
    default=-1.0,
    notes="Default is -1, corresponding to a single (weak-coupling) temperature bath.",
    value_type="float",
    section="cntrl",
    category="Temperature regulation",
    related=["temp0"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="tempi",
    description=(
        "Initial temperature. For the initial dynamics run, (ntx = 1) the velocities are assigned from "
        "a Maxwellian distribution at tempi K. If tempi = 0.0, the velocities will be calculated from "
        "the forces instead. tempi has no effect if ntx = 5."
    ),
    default=0.0,
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Temperature regulation",
    related=["ntx", "temp0", "ig"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ig",
    description=(
        "The seed for the pseudo-random number generator. The MD starting velocity is dependent on the "
        "random number generator seed if tempi is nonzero and ntx = 1. The value of this seed also "
        "affects the set of pseudo-random values used for Langevin dynamics or Andersen-like coupling, "
        "and hence should be set to a different value on each restart if ntt = 2 or 3. If ig=-1 (the "
        "default) then the random seed will be based on the current date and time, and hence will be "
        "different for every run."
    ),
    default=-1,
    notes=(
        "Unless you specifically desire reproducibility, it is recommended that you set ig=-1 for all "
        "runs involving ntt = 2 or 3."
    ),
    value_type="int",
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "tempi", "gamma_ln", "vrand"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="tautp",
    description=(
        "Time constant, in ps, for heat bath coupling for the system, if ntt = 1. Generally, values for "
        "TAUTP should be in the range of 0.5-5.0 ps, with a smaller value providing tighter coupling to "
        "the heat bath and, thus, faster heating and a less natural trajectory. Smaller values of TAUTP "
        "result in smaller fluctuations in kinetic energy, but larger fluctuations in the total energy. "
        "Values much larger than the length of the simulation result in a return to constant energy conditions."
    ),
    default=1.0,
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "temp0"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="gamma_ln",
    description=(
        "The collision frequency, in ps^-1, when ntt = 3. A simple Leapfrog integrator is used to "
        "propagate the dynamics, with the kinetic energy adjusted to be correct for the harmonic oscillator "
        "case. Note that it is not necessary that gamma_ln approximate the physical collision frequency, "
        "which is about 50 ps^-1 for liquid water. In fact, it is often advantageous, in terms of sampling "
        "or stability of integration, to use much smaller values, around 2 to 5 ps^-1. For implicit solvent "
        "(GB), even much lower values may be useful: for example, setting gamma_ln to 0.01 ps^-1 can lead "
        "to significant, up to 100-fold in some cases, speedup of conformational sampling."
    ),
    default=0.0,
    notes=(
        "Also used to determine thermostat coupling constant for the Optimized Isokinetic Nose-Hoover "
        "chain integrator (OIN, ntt=9), which is equal to 1/gamma_ln, so the specified gamma_ln must be "
        "> 0. For ntt=10, this is the friction constant associated with the stochastic component of the "
        "integrator and must be > 0."
    ),
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "ig"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="vrand",
    description=(
        "If vrand>0 and ntt=2, the velocities will be randomized to temperature TEMP0 every vrand steps."
    ),
    default=1000,
    value_type="int",
    min_val=0,
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "temp0"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="vlimit",
    description=(
        "If not equal to 0.0, then any component of the velocity that is greater than abs(VLIMIT) will "
        "be reduced to VLIMIT (preserving the sign). This can be used to avoid occasional instabilities "
        "in molecular dynamics runs. VLIMIT should generally be set to a value like 20 (the default), "
        "which is well above the most probable velocity in a Maxwell-Boltzmann distribution at room "
        "temperature. A warning message will be printed whenever the velocities are modified. Runs that "
        "have more than a few such warnings should be carefully examined."
    ),
    default=20.0,
    value_type="float",
    section="cntrl",
    category="Temperature regulation",
    related=["ntt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nkija",
    description=(
        "For use with ntt=9 and ntt=10. For ntt=9, this is the number of substeps of dt when integrating "
        "the thermostat equations of motion, for greater accuracy. For ntt=10, this specifies the number "
        "of additional auxiliary velocity variables v1 and v2, which will total nkija*v1 + nkija*v2."
    ),
    default=1,
    value_type="int",
    min_val=1,
    section="cntrl",
    category="Temperature regulation",
    related=["ntt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="idistr",
    description=(
        "For the isokinetic integrator (ntt=9), the frequency at which the thermostat velocity "
        "distribution functions are accumulated."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "nkija"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="sinrtau",
    description=(
        "For the SINR (Stochastic Isokinetic Nose-Hoover RESPA) integrator (ntt=10), this specifies "
        "the time scale for determining the masses associated with the two auxiliary velocity variables "
        "v1 and v2 (e.g. thermostat velocities) and hence the magnitude of the coupling of the physical "
        "velocities with the auxiliary velocities. Generally this should be related to the time scale of "
        "the system."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Temperature regulation",
    related=["ntt", "nkija"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.8 Pressure regulation
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ntp",
    description=(
        "Flag for constant pressure dynamics. This option should be set to 1 or 2 when Constant Pressure "
        "periodic boundary conditions are used."
    ),
    default=0,
    options={
        0: "No pressure scaling (Default)",
        1: "MD with isotropic position scaling",
        2: (
            "MD with anisotropic (x-,y-,z-) pressure scaling: this should only be used with orthogonal "
            "boxes (i.e. with all angles set to 90 degrees). Anisotropic scaling is primarily intended "
            "for non-isotropic systems, such as membrane simulations, where the surface tensions are "
            "different in different directions. Anisotropic pressure scaling can also be applied to just "
            "one specified direction (x, y or z) with the directional pressure scaling option "
            "(baroscalingdir > 0)."
        ),
        3: (
            "MD with semiisotropic pressure scaling: this is only available with constant surface tension "
            "(csurften > 0) and orthogonal boxes. This links the pressure coupling in the two directions "
            "tangential to the interface."
        ),
        4: (
            "MD towards a targeted volume. This is not for production but for modifying the volume of the "
            "system, particularly useful for preparing REMD simulations where the shape of each replica "
            "needs to be the same."
        ),
    },
    notes=(
        "It is generally wise to equilibrate the temperature to something like the final temperature using "
        "constant volume (ntp=0) before switching on constant pressure simulations to adjust the system to "
        "the correct density. If you fail to do this, the program will try to adjust the density too quickly, "
        "and bad things (such as SHAKE failures) are likely to happen."
    ),
    value_type="int",
    section="cntrl",
    category="Pressure regulation",
    related=["barostat", "pres0", "taup", "ntb", "comp"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="barostat",
    description="Flag used to control which barostat to use in order to control the pressure.",
    default=1,
    options={
        1: "Berendsen (Default)",
        2: "Monte Carlo barostat",
    },
    notes=(
        "While the Berendsen barostat yields the correct target density, it does not strictly sample from "
        "the isothermal-isobaric ensemble and typically yields volume fluctuations that are too low. The "
        "Monte Carlo barostat, on the other hand, samples rigorously from the isothermal-isobaric ensemble. "
        "For GPU performance, use barostat=2: Performance will generally be NVE>NVT>NPT (NVT~NPT for "
        "barostat=2)."
    ),
    value_type="int",
    section="cntrl",
    category="Pressure regulation",
    related=["ntp", "mcbarint", "pres0"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="mcbarint",
    description="Number of steps between volume change attempts performed as part of the Monte Carlo barostat.",
    default=100,
    value_type="int",
    min_val=1,
    section="cntrl",
    category="Pressure regulation",
    related=["barostat", "ntp"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="pres0",
    description=(
        "Reference pressure (in units of bars, where 1 bar ~ 0.987 atm) at which the system is maintained "
        "(when NTP > 0)."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Pressure regulation",
    related=["ntp", "barostat"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="comp",
    description=(
        "Compressibility of the system when NTP > 0. The units are in 1.0 x 10^-6 bar^-1; a value of "
        "44.6 (default) is appropriate for water."
    ),
    default=44.6,
    value_type="float",
    section="cntrl",
    category="Pressure regulation",
    related=["ntp"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="taup",
    description=(
        "Pressure relaxation time (in ps), when NTP > 0. The recommended value is between 1.0 and 5.0 "
        "psec. Larger values may sometimes be necessary (if your trajectories seem unstable)."
    ),
    default=1.0,
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Pressure regulation",
    related=["ntp", "barostat"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="baroscalingdir",
    description=(
        "Flag for pressure scaling direction control. Applicable when using Monte Carlo barostat "
        "(barostat = 2) with anisotropic pressure scaling (ntp = 2)."
    ),
    default=0,
    options={
        0: "Box size scales randomly (x, y or z) each scaling step (default)",
        1: "Box scales only along x-direction, dimensions along y-, z-axes are fixed",
        2: "Box scales only along y-direction, dimensions along x-, z-axes are fixed",
        3: "Box scales only along z-direction, dimensions along x-, y-axes are fixed",
    },
    value_type="int",
    section="cntrl",
    category="Pressure regulation",
    related=["barostat", "ntp"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="csurften",
    description=(
        "Flag for constant surface tension dynamics. Constant surface tension is used in statistical "
        "ensembles for simulating liquid interfaces. This is primarily intended for lipid membrane "
        "simulations with two or more interfaces."
    ),
    default=0,
    options={
        0: "No constant surface tension (default)",
        1: "Constant surface tension with interfaces in the yz plane",
        2: "Constant surface tension with interfaces in the xz plane",
        3: "Constant surface tension with interfaces in the xy plane",
    },
    notes=(
        "In order to use constant surface tension, periodic boundary conditions (ntb = 2), anisotropic "
        "or semiisotropic pressure scaling (ntp = 2 or ntp = 3), and an orthogonal box must be used."
    ),
    value_type="int",
    section="cntrl",
    category="Pressure regulation",
    related=["ntp", "gamma_ten", "ninterface"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="gamma_ten",
    description="Surface tension value in units of dyne/cm.",
    default=0.0,
    value_type="float",
    section="cntrl",
    category="Pressure regulation",
    related=["csurften"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ninterface",
    description=(
        "Number of interfaces in the periodic box. There must be at least two interfaces in the periodic "
        "box. Two interfaces is appropriate for a lipid bilayer system and is the default value."
    ),
    default=2,
    value_type="int",
    min_val=2,
    section="cntrl",
    category="Pressure regulation",
    related=["csurften"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.9 SHAKE bond length constraints
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ntc",
    description=(
        "Flag for SHAKE to perform bond length constraints. (See also NTF in the Potential function "
        "section. In particular, typically NTF = NTC.) The SHAKE option should be used for most MD "
        "calculations. The size of the MD timestep is determined by the fastest motions in the system. "
        "SHAKE removes the bond stretching freedom, which is the fastest motion, and consequently allows "
        "a larger timestep to be used. For water models, a special 'three-point' algorithm is used. "
        "Consequently, to employ TIP3P set NTF = NTC = 2."
    ),
    default=1,
    options={
        1: "SHAKE is not performed (default)",
        2: "Bonds involving hydrogen are constrained",
        3: "All bonds are constrained (not available for parallel or qmmm runs in sander)",
    },
    notes=(
        "Since SHAKE is an algorithm based on dynamics, the minimizer is not aware of what SHAKE is doing; "
        "for this reason, minimizations generally should be carried out without SHAKE. One exception is "
        "short minimizations whose purpose is to remove bad contacts before dynamics can begin."
    ),
    value_type="int",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["ntf", "dt", "tol"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="tol",
    description=(
        "Relative geometrical tolerance for coordinate resetting in SHAKE. Recommended maximum: <0.00005 "
        "Angstrom."
    ),
    default=0.00001,
    value_type="float",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["ntc"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="jfastw",
    description=(
        "Fast water definition flag. By default, the system is searched for water residues, and special "
        "routines are used to SHAKE these systems."
    ),
    default=0,
    options={
        0: "Normal operation. Waters are identified by the default names (given below), unless they are redefined.",
        4: "Do not use the fast SHAKE routines for waters.",
    },
    value_type="int",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["ntc", "WATNAM", "OWTNM", "HWTNM1", "HWTNM2"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="WATNAM",
    description="The residue name the program expects for water.",
    default="WAT ",
    value_type="str",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["jfastw"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="OWTNM",
    description="The atom name the program expects for the oxygen of water.",
    default="O   ",
    value_type="str",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["jfastw"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="HWTNM1",
    description="The atom name the program expects for the 1st H of water.",
    default="H1  ",
    value_type="str",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["jfastw"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="HWTNM2",
    description="The atom name the program expects for the 2nd H of water.",
    default="H2  ",
    value_type="str",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["jfastw"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="noshakemask",
    description=(
        "String that specifies atoms that are not to be shaken (assuming that ntc>1). Any bond that "
        "would otherwise be shaken by virtue of the ntc flag, but which involves an atom flagged here, "
        "will *not* be shaken."
    ),
    default="",
    notes=(
        "If this option is invoked, then all parts of the potential must be evaluated, that is, ntf must "
        "be one. The code enforces this by setting ntf to 1 when a noshakemask string is present in the "
        "input. If you want the noshakemask to apply to all or part of the water molecules, you must also "
        "set jfastw=4, to turn off the special code for water SHAKE."
    ),
    value_type="str",
    section="cntrl",
    category="SHAKE bond length constraints",
    related=["ntc", "ntf", "jfastw"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.10 The "middle" scheme
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ischeme",
    description="Flag for choosing an integration scheme for molecular dynamics.",
    default=0,
    options={
        0: "Conventional scheme in AMBER.",
        1: "'Middle' scheme based on the leapfrog algorithm.",
    },
    notes=(
        "The 'middle' scheme offers a unified framework to develop efficient thermostatting algorithms "
        "for configurational sampling for the canonical ensemble. It can be implemented for performing "
        "molecular dynamics (MD) or path integral molecular dynamics (PIMD), either with or without "
        "holonomic constraints. The 'middle' scheme allows the use of much larger time intervals to "
        "maintain the same accuracy, which significantly improves the configurational sampling efficiency."
    ),
    value_type="int",
    section="cntrl",
    category="Middle scheme",
    related=["ithermostat", "therm_par"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ithermostat",
    description=(
        "Flag for different thermostats when the 'middle' scheme is employed. Two types of thermostats "
        "are currently available."
    ),
    default=1,
    options={
        1: "Langevin dynamics.",
        2: "Andersen thermostat.",
    },
    value_type="int",
    section="cntrl",
    category="Middle scheme",
    related=["ischeme", "therm_par"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="therm_par",
    description=(
        "The parameter used in a thermostatting method of the 'middle' scheme, in the unit of ps^-1, "
        "which should always be a positive number. It refers to the friction coefficient for Langevin "
        "dynamics (ithermostat = 1) or the collision frequency for the Andersen thermostat "
        "(ithermostat = 2)."
    ),
    default=None,
    notes=(
        "The recommended value for therm_par is related to the characteristic frequency of the specific "
        "system. For a liquid water system (216 water molecules in a cell with periodic boundary "
        "conditions) with no holonomic constraints, the thermostat parameter is usually chosen to be "
        "2-50 ps^-1."
    ),
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Middle scheme",
    related=["ischeme", "ithermostat"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.11 Water cap
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ivcap",
    description=(
        "Flag to control cap option. The 'cap' refers to a spherical portion of water centered on a "
        "point in the solute and restrained by a soft half-harmonic potential. For the best physical "
        "realism, this option should be combined with igb=10, in order to include the reaction field of "
        "waters that are beyond the cap radius."
    ),
    default=0,
    options={
        0: "Cap will be in effect if it is in the prmtop file (default).",
        1: (
            "With this option, a cap can be excised from a larger box of water. For this, cutcap (i.e., "
            "the radius of the cap), xcap, ycap, and zcap (i.e., the location of the center of the cap) "
            "need to be specified in the &cntrl namelist."
        ),
        2: "Cap will be inactivated, even if parameters are present in the prmtop file.",
        5: (
            "With this option, a shell of water around a solute can be excised from a larger box of water. "
            "For this, cutcap (i.e., the thickness of the shell) needs to be specified in the &cntrl "
            "namelist. This option only works for a single-step minimization."
        ),
    },
    value_type="int",
    section="cntrl",
    category="Water cap",
    related=["fcap", "cutcap", "xcap", "ycap", "zcap", "igb"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="fcap",
    description="The force constant for the cap restraint potential.",
    default=None,
    value_type="float",
    section="cntrl",
    category="Water cap",
    related=["ivcap"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cutcap",
    description="Radius of the cap, if ivcap=1 is used. Thickness of the shell, if ivcap=5.",
    default=None,
    value_type="float",
    section="cntrl",
    category="Water cap",
    related=["ivcap", "xcap", "ycap", "zcap"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="xcap",
    description="X location of the cap center, if ivcap=1 is used.",
    default=None,
    value_type="float",
    section="cntrl",
    category="Water cap",
    related=["ivcap", "cutcap"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ycap",
    description="Y location of the cap center, if ivcap=1 is used.",
    default=None,
    value_type="float",
    section="cntrl",
    category="Water cap",
    related=["ivcap", "cutcap"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="zcap",
    description="Z location of the cap center, if ivcap=1 is used.",
    default=None,
    value_type="float",
    section="cntrl",
    category="Water cap",
    related=["ivcap", "cutcap"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.12 NMR refinement options
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="iscale",
    description=(
        "Number of additional variables to optimize beyond the 3N structural parameters. (Default = 0). "
        "At present, this is only used with residual dipolar coupling and CSA or pseudo-CSA restraints."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="noeskp",
    description=(
        "The NOESY volumes will only be evaluated if mod(nstep, noeskp) = 0; otherwise the last "
        "computed values for intensities and derivatives will be used."
    ),
    default=1,
    notes="Default = 1, i.e. evaluate volumes at every step.",
    value_type="int",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ipnlty",
    description=(
        "This parameter determines the functional form of the penalty function for NOESY volume and "
        "chemical shift restraints."
    ),
    default=1,
    options={
        1: "The program will minimize the sum of the absolute values of the errors; this is akin to minimizing the crystallographic R-factor (default).",
        2: "The program will optimize the sum of the squares of the errors.",
        3: "For NOESY intensities, the penalty will be of the form awt[I^(1/6)_c - I^(1/6)_o]^2. Chemical shift penalties will be as for ipnlty=1.",
    },
    value_type="int",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mxsub",
    description=(
        "Maximum number of submolecules that will be used. This is used to determine how much space "
        "to allocate for the NOESY calculations."
    ),
    default=1,
    value_type="int",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scalm",
    description=(
        "'Mass' for the additional scaling parameters. Right now they are restricted to all have the "
        "same value. The larger this value, the slower these extra variables will respond to their "
        "environment."
    ),
    default=100.0,
    notes="Units: amu.",
    value_type="float",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt", "iscale"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="pencut",
    description=(
        "In the summaries of the constraint deviations, entries will only be made if the penalty for "
        "that term is greater than PENCUT."
    ),
    default=0.1,
    value_type="float",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="tausw",
    description=(
        "For noesy volume calculations (NMROPT = 2), intensities with mixing times less than TAUSW "
        "(in seconds) will be computed using perturbation theory, whereas those greater than TAUSW will "
        "use a more exact theory. To always use the 'exact' intensities and derivatives, set TAUSW = 0.0; "
        "to always use perturbation theory, set TAUSW to a value larger than the largest mixing time in "
        "the input."
    ),
    default=0.1,
    notes="Default is TAUSW of 0.1 second, which should work pretty well for most systems.",
    value_type="float",
    section="cntrl",
    category="NMR refinement options",
    related=["nmropt"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.6.13 EMAP restraints
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="iemap",
    description=(
        "Turn on EMAP restrained simulation when iemap>0. EMAP restraint information must be input "
        "from &emap namelists in the input file."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="EMAP restraints",
    related=["gammamap"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="gammamap",
    description="Friction constant for the EMAP restraint maps when allowed to move.",
    default=1.0,
    notes="Units: 1/ps.",
    value_type="float",
    section="cntrl",
    category="EMAP restraints",
    related=["iemap"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.7.1 Generic potential function parameters
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ntf",
    description=(
        "Force evaluation. Note: If SHAKE is used (see NTC), it is not necessary to calculate forces "
        "for the constrained bonds."
    ),
    default=1,
    options={
        1: "Complete interaction is calculated (default)",
        2: "Bond interactions involving H-atoms omitted (use with NTC=2)",
        3: "All the bond interactions are omitted (use with NTC=3)",
        4: "Angle involving H-atoms and all bonds are omitted",
        5: "All bond and angle interactions are omitted",
        6: "Dihedrals involving H-atoms and all bonds and all angle interactions are omitted",
        7: "All bond, angle and dihedral interactions are omitted",
        8: "All bond, angle, dihedral and non-bonded interactions are omitted",
    },
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["ntc"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntb",
    description=(
        "This variable controls whether or not periodic boundaries are imposed on the system during the "
        "calculation of non-bonded interactions. Bonds spanning periodic boundaries are not yet supported. "
        "There is no longer any need to set this variable, since it can be determined from igb and ntp "
        "parameters."
    ),
    default=1,
    options={
        0: "No periodicity is applied and PME is off (default when igb > 0)",
        1: "Constant volume (default when igb and ntp are both 0, which are their defaults)",
        2: "Constant pressure (default when ntp > 0)",
    },
    notes=(
        "If NTB is nonzero then there must be a periodic boundary in the topology file. Constant pressure "
        "is not used in minimization (IMIN=1). For a periodic system, constant pressure is the only way "
        "to equilibrate density if the starting state is not correct. Almost every system needs to be "
        "equilibrated at constant pressure (ntb=2, ntp>0) to get to a proper density. But be sure to "
        "equilibrate first (at constant volume) to something close to the final temperature, before "
        "turning on constant pressure."
    ),
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["igb", "ntp", "cut"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="dielc",
    description=(
        "Dielectric multiplicative constant for the electrostatic interactions. Please note this is NOT "
        "related to dielectric constants for generalized Born or Poisson-Boltzmann calculations. It should "
        "only be used for quasi-vacuum simulations."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cut",
    description=(
        "This is used to specify the nonbonded cutoff, in Angstroms. For PME, the cutoff is used to limit "
        "direct space sum, and 8.0 is usually a good value. When igb>0, the cutoff is used to truncate "
        "nonbonded pairs (on an atom-by-atom basis); here a larger value than the default is generally "
        "required."
    ),
    default=8.0,
    notes=(
        "When igb > 0, the default is 9999.0 (effectively infinite). "
        "When igb == 0, the default is 8.0."
    ),
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Potential function - generic parameters",
    related=["ntb", "igb"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="fswitch",
    description=(
        "When off, fswitch<=0, uses a truncation cutoff. When on fswitch>0, sets a force switching region "
        "where the force cutoff smoothly approaches 0 between the region of the fswitch value to the cut "
        "value. Force values below the fswitch value follow the standard Lennard-Jones force."
    ),
    default=-1.0,
    notes=(
        "This option is not supported for use with GB (i.e., only igb=0 and ntb>0), nor is it compatible "
        "with the 12-6-4 and pairwise 12-6-4 Lennard-Jones models (lj1264=1 and plj1264=1). Due to "
        "performance regressions (about 20%) with running with the force switching on, it is recommended "
        "that simulations run with fswitch off unless using a force field that requires or recommends "
        "using the force switch."
    ),
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["cut", "lj1264", "plj1264"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nsnb",
    description=(
        "Determines the frequency of nonbonded list updates when igb=0 and nbflag=0; see the description "
        "of nbflag for more information."
    ),
    default=25,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["nbflag"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ipol",
    description="When set to 1, use a polarizable force field.",
    default=0,
    options={0: "No polarizable force field.", 1: "Use a polarizable force field."},
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["indmeth", "diptol"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ipgm",
    description="When set to 1, use the polarizable Gaussian Multipole force field.",
    default=0,
    options={0: "No pGM force field.", 1: "Use the pGM force field."},
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ifqnt",
    description="Flag for QM/MM run; if set to 1, you must also include a &qmmm namelist.",
    default=0,
    options={0: "No QM/MM.", 1: "QM/MM run."},
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="igb",
    description=(
        "Flag for using the generalized Born implicit solvent models. "
        "Generalized Born simulations can only be run for non-periodic systems, i.e. where ntb=0. "
        "Unlike its use in explicit solvent PME simulations, short nonbonded cutoff values have much "
        "stronger impact on accuracy of the GB calculations. Essentially, any cutoff values other than "
        "cut > structure size can lead to artifacts. Current GPU implementation of the GB can not use cutoffs."
    ),
    default=0,
    options={
        0: "No generalized Born term is used. (Default)",
        1: (
            "The Hawkins, Cramer, Truhlar pairwise generalized Born model is used, with parameters "
            "described by Tsui and Case. This model uses the default radii set up by LEaP. Use LEaP "
            "command 'set default PBRadii mbondi'. Most nucleic acid simulations have used this model."
        ),
        2: (
            "Use a modified GB model developed by A. Onufriev, D. Bashford and D.A. Case (GBOBC). "
            "The effective Born radii are re-scaled to account for interstitial spaces between atom "
            "spheres missed by the GBHCT approximation. Parameters alpha=0.8, beta=0.0, gamma=2.909125 "
            "(model I). Use LEaP command 'set default PBRadii mbondi2'."
        ),
        5: (
            "Same as igb=2, except that alpha=1.0, beta=0.8, gamma=4.85. This corresponds to model II "
            "of the OBC model. Use 'set default PBRadii mbondi2' or 'set default PBRadii bondi'. "
            "The igb=5 variant agrees better with Poisson-Boltzmann treatment in calculating the "
            "electrostatic part of solvation free energy on an extensive set of protein structures."
        ),
        6: (
            "No continuum solvent model used at all; this corresponds to a non-periodic, 'vacuum', model "
            "where the non-bonded interactions are just Lennard-Jones and Coulomb interactions."
        ),
        7: (
            "The GBn model described by Mongan, Simmerling, McCammon, Case and Onufriev is employed. "
            "Uses a pairwise correction term to GBHCT to approximate a molecular surface dielectric "
            "boundary. Use the bondi radii set. Not recommended for systems involving nucleic acids."
        ),
        8: (
            "Same GB functional form as the GBn model (igb=7), but with different parameters. "
            "The offset, overlap screening parameters, and gbneckscale are changed. Individual alpha, "
            "beta, and gamma parameters can be specified for each element H, C, N, O, S, P. "
            "Separate parameters for proteins and nucleic acids. The combination of ff14SBonlysc with "
            "igb=8 gives the best results for proteins, nucleic acids and protein-nucleic acid complexes. "
            "Use LEaP command 'set default PBRadii mbondi3'."
        ),
        10: (
            "Calculate the reaction field and nonbonded interactions using a numerical Poisson-Boltzmann "
            "solver. This option is described in Chapter 6. Note that this is not a generalized Born "
            "simulation, in spite of its use of igb; it is rather an alternative continuum solvent model."
        ),
    },
    notes=(
        "Recommended radii sets: igb=1 -> mbondi, igb=2 -> mbondi2, igb=5 -> mbondi2, "
        "igb=7 -> bondi, igb=8 -> mbondi3. Values 3 and 4 are unused (were used in Amber 7 "
        "for parameter sets no longer supported). GB models are not compatible with polarizable "
        "force fields. If the nonbonded cutoff is used in GB calculations, it should be greater "
        "than that for PME calculations, perhaps cut=16."
    ),
    value_type="int",
    section="cntrl",
    category="Generalized Born",
    related=["ntb", "cut", "intdiel", "extdiel", "saltcon", "rgbmax", "offset", "gbsa", "alpb"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ipb",
    description=(
        "Option to set up a dielectric model for all numerical Poisson-Boltzmann procedures. "
        "IPB=1 corresponds to a classical geometric method, while a level-set based algebraic "
        "method is used when IPB >= 2. Note: in the standalone pbsa program, this is set in &cntrl. "
        "In sander, you must explicitly set IPB to nonzero to invoke pbsa functionalities."
    ),
    default=0,
    options={
        0: "No electrostatic solvation free energy is computed.",
        1: "The dielectric interface between solvent and solute is built with a geometric approach.",
        2: (
            "The dielectric interface is implemented with the level set function. Use of a level set "
            "function simplifies the calculation of the intersection points of the molecular surface "
            "and grid edges and leads to more stable numerical calculations. Default for pbsa."
        ),
        4: (
            "The dielectric interface is implemented with the level set function. The linear equations "
            "on grid points nearby the dielectric boundary are constructed using the IIM. SMOOTHOPT is "
            "useless with this option. Only linear PB equation is supported (NPBOPT=0)."
        ),
        6: (
            "The dielectric interface is implemented analytically with the revised density function "
            "approach (SASOPT=2). The linear equations on irregular points use IIM with the analytical "
            "surface. Otherwise same as IPB=4."
        ),
        7: (
            "The dielectric interface is implemented analytically with the revised density function "
            "approach (SASOPT=2). The linear equations on irregular points use the alpha-factor "
            "harmonic average method."
        ),
        8: (
            "The dielectric interface is implemented analytically with the revised density function "
            "approach (SASOPT=2). The linear equations on irregular points use the second-order "
            "harmonic average method."
        ),
    },
    notes="See Chapter 6 for detailed information. The default IPB is 2 for the standalone pbsa program.",
    value_type="int",
    section="cntrl",
    category="Poisson-Boltzmann",
    related=["igb", "inp", "epsin", "epsout"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="irism",
    description="Flag for 3D-reference interaction site model (RISM) molecular solvation method.",
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ievb",
    description="If set to 1, use the empirical valence bond method to compute energies and forces.",
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="iamoeba",
    description=(
        "Flag for using the amoeba polarizable potentials of Ren and Ponder. When this option is set "
        "to 1, you need to prepare an amoeba namelist with additional parameters."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="lj1264",
    description=(
        "When the Lennard-Jones C-coefficient is present in the prmtop file, the 12-6-4 potential is "
        "active. Setting lj1264 to 0 when the C-coefficient is present will forcibly turn off the 12-6-4 "
        "potential. Setting lj1264 to 1 when no C-coefficient is present will result in a fatal error."
    ),
    default=0,
    notes=(
        "It is currently only compatible with the Particle Mesh Ewald method for long-range "
        "electrostatics."
    ),
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["plj1264", "fswitch"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="plj1264",
    description=(
        "Similar to lj1264 above, if not present in the input file, this keyword will still automatically "
        "turn to 1 as long as D-coefficient is found in the prmtop file."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["lj1264"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="efx",
    description=(
        "This sets the x component of the electric field in kcal/(mol*A*e). Electric fields are "
        "naturally off if efx, efy, efz are 0. It currently only supports pmemd (both the serial "
        "and MPI versions)."
    ),
    default=0.0,
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["efy", "efz", "efn", "efphase", "effreq"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="efy",
    description=(
        "This sets the y component of the electric field in kcal/(mol*A*e). Electric fields are "
        "naturally off if efx, efy, efz are 0."
    ),
    default=0.0,
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["efx", "efz"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="efz",
    description=(
        "This sets the z component of the electric field in kcal/(mol*A*e). Electric fields are "
        "naturally off if efx, efy, efz are 0."
    ),
    default=0.0,
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["efx", "efy"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="efn",
    description=(
        "If efn is on (efn=1), the x, y, z (efx, efy, efz) components are scaled to box size. "
        "This normalizes the electric field charge to your box size. It is off when it is 0."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["efx", "efy", "efz"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="efphase",
    description=(
        "efphase sets the timestep phase for the electric field using the equation "
        "cos((2*pi*effreq/1000)*(dt*step) - (pi*efphase/180))."
    ),
    default=None,
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["effreq", "efx", "efy", "efz"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="effreq",
    description=(
        "effreq sets the timestep frequency for the electric field using the equation "
        "cos((2*pi*effreq/1000)*(dt*step) - (pi*efphase/180))."
    ),
    default=None,
    value_type="float",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["efphase", "efx", "efy", "efz"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mcwat",
    description=(
        "mcwat determines whether the Monte Carlo (MC) water equilibrium feature is used. Set 1 to run, "
        "0 otherwise. Currently only supported on pmemd, pmemd.cuda, and REMD runs of pmemd.cuda.MPI."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - MC water",
    related=["nmd", "nmc", "mcwatmask", "mcligshift"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nmd",
    description="nmd is the number of MD steps per MC/MD cycle.",
    default=1000,
    value_type="int",
    section="cntrl",
    category="Potential function - MC water",
    related=["mcwat", "nmc"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nmc",
    description="nmc is the number of MC steps that are performed per MC/MD cycle.",
    default=100000,
    value_type="int",
    section="cntrl",
    category="Potential function - MC water",
    related=["mcwat", "nmd"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mcwatmask",
    description=(
        "mcwatmask is an amber mask selection of atoms which we calculate the center of mass to define "
        "the center of the MC region. The MC region is the region in which water molecules can do "
        "translational MC moves."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - MC water",
    related=["mcwat", "mcligshift"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mcligshift",
    description=(
        "mcligshift is used to define the length of the grid along each axis starting from the center "
        "of mass of mcwatmask, in Angstroms. For example, a value of 10 means that the grid will extend "
        "10 A from the center along each axis."
    ),
    default=None,
    value_type="float",
    section="cntrl",
    category="Potential function - MC water",
    related=["mcwat", "mcwatmask"],
    commonly_changed=False,
))

# TI decomposition (pmemd specific)
KEYWORDS.append(Keyword(
    name="ntwd",
    description=(
        "(default 0), set to 1 to turn on Thermodynamic Integration (TI) Free Energy Decomposition in "
        "the pmemd.decomp executable. When on, the command line flag -decomp can be set to determine the "
        "log file that outputs decomposition data."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - TI decomposition",
    related=["decompmask", "ligmask", "proteinmask", "cofactormask", "reweight"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="decompmask",
    description=(
        "An amber mask string that sets the atoms to be displayed in the decomp.log file (this is to "
        "help with potential file size issues). You can just set to all atoms to get the decomposition "
        "of all atoms."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - TI decomposition",
    related=["ntwd"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ligmask",
    description=(
        "Amber mask of atoms that are considered the ligands in your calculation. Generally these should "
        "be set to your timask1 and timask2 atoms."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - TI decomposition",
    related=["ntwd", "proteinmask", "cofactormask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="proteinmask",
    description=(
        "Amber mask of atoms that are considered part of your protein. This is used in the calculation "
        "for region decomposition specifically the ligand-protein contribution."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - TI decomposition",
    related=["ntwd", "ligmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cofactormask",
    description=(
        "Amber mask of atoms that are considered part of your cofactor. This is used in the calculation "
        "for region decomposition specifically the ligand-cofactor contribution."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - TI decomposition",
    related=["ntwd", "ligmask"],
    commonly_changed=False,
))

# RAMD
KEYWORDS.append(Keyword(
    name="ramdboost",
    description=(
        "Sets default random boost acceleration for ramd (default 1). This boost is multiplied by the "
        "mass of each atom in the ligand to determine the force each atom receives. This value is in "
        "internal acceleration units. Ramd is a pmemd and pmemd.cuda only feature and does not support MPI."
    ),
    default=1,
    value_type="float",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboostfreq", "ramdboostrate", "ramdint", "ramdmaxdist", "ramdligmask", "ramdprotmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ramdboostfreq",
    description="Sets number of steps between each time ramd boost strength is increased.",
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboost"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ramdboostrate",
    description="Sets the amount to increase the ramdboost acceleration each time ramdboostfreq.",
    default=0,
    value_type="float",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboost", "ramdboostfreq"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ramdint",
    description="Sets the time step interval to apply ramd boost on to the ligand.",
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboost"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ramdmaxdist",
    description=(
        "Determines the end condition for the simulation (ramd terminates when nstlim is reached or when "
        "ramdmaxdist is satisfied). ramdmaxdist is the amount of angstrom displacement from initial center "
        "of mass distance of protein and ligand to when this displacement increases by ramdmaxdist."
    ),
    default=None,
    value_type="float",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboost", "nstlim"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ramdligmask",
    description="Amber selection mask for what is considered the ligand that needs to be boosted in ramd.",
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboost", "ramdprotmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ramdprotmask",
    description=(
        "Amber selection mask for what is considered the protein that is used to calculate the distance "
        "the ligand has moved."
    ),
    default="",
    value_type="str",
    section="cntrl",
    category="Potential function - RAMD",
    related=["ramdboost", "ramdligmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="reweight",
    description=(
        "Allows the re-evaluation of trajectories (usually with a new parameter file). Set 1 to turn on. "
        "When running this command, in the topology command of the run file (-c) place the trajectory "
        "instead of the topology file. This supports netcdf only. Reweight is supported in pmemd and "
        "pmemd.cuda, and does not support MPI."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=["ntwd"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="midpoint",
    description=(
        "Turns on midpoint optimizations (usage of 3-D spatial decomposition). 1 is on, 0 is off "
        "(default). This switch is currently experimental. Currently only supported on pmemd.MPI."
    ),
    default=0,
    value_type="int",
    section="cntrl",
    category="Potential function - generic parameters",
    related=[],
    commonly_changed=False,
))


# ---------------------------------------------------------------------------
# Chapter 27: Discrete Constant pH MD (&cntrl keywords)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="icnstph",
    description=(
        "Flag to turn on constant pH molecular dynamics using discrete protonation states "
        "with Monte Carlo sampling. Requires a cpin file specified via the -cpin command-line "
        "flag. Reference energies and protonation state definitions are pre-computed and stored "
        "in the cpin file generated by cpinutil.py."
    ),
    default=0,
    options={
        0: "No constant pH (default).",
        1: (
            "Constant pH MD in implicit solvent (Generalized Born). Protonation state changes "
            "are attempted every ntcnstph steps using a Monte Carlo scheme. One residue is examined "
            "per MC step. The GB model and parameters used must match those used to derive the "
            "reference energies in the cpin file."
        ),
        2: (
            "Constant pH MD in explicit solvent. Since direct protonation state changes in explicit "
            "solvent are opposed by solvent orientation, a hybrid MD/MC approach is used: MD is run "
            "for ntcnstph steps, solvent and ions are stripped, protonation state changes are attempted "
            "for all titratable residues in random order using a GB implicit solvent model, then "
            "solvent is restored and relaxation dynamics are run for ntrelax steps if any states changed."
        ),
    },
    notes=(
        "When icnstph=1 or 2, you must also specify solvph and ntcnstph in &cntrl. "
        "The cpin file is specified on the command line with -cpin, the protonation state "
        "history is written to the file specified by -cpout, and a restart file is written "
        "to -cprestrt. The cprestrt file should be used as the cpin file when restarting. "
        "For explicit solvent (icnstph=2), also set ntrelax for solvent relaxation steps, and "
        "use cpinutil.py with -op to generate a modified topology file with adjusted GB radii "
        "for AS4/GL4 carboxylate oxygens. Reference energies were derived with: "
        "cut=30.0, igb=#, saltcon=0.1, nrespa=1, temp0=300.0, ntc=2, ntf=2 "
        "using the ff99SB force field. Force fields with the same charge scheme (e.g., ff14SB) "
        "should be valid, but others (e.g., ff03, ff13) require recalculating reference energies."
    ),
    value_type="int",
    section="cntrl",
    category="Constant pH MD",
    related=["solvph", "ntcnstph", "ntrelax", "saltcon", "igb"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="solvph",
    description=(
        "Sets the solvent pH value for constant pH MD simulations. This is an external parameter "
        "that affects the protonation state distribution through the Metropolis Monte Carlo criteria. "
        "The free energy of a protonation state change is evaluated relative to this pH using the "
        "thermodynamic cycle shown in Figure 27.1 of the AMBER manual."
    ),
    default=7.0,
    notes=(
        "The choice of solvph determines the protonation state equilibrium. For calibrating new "
        "titratable residues using finddgref.py, solvph should be set equal to the pKa of the "
        "model compound. For production simulations, set to the desired simulation pH. "
        "For pH-REMD, different replicas use different solvph values."
    ),
    value_type="float",
    section="cntrl",
    category="Constant pH MD",
    related=["icnstph", "ntcnstph"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntcnstph",
    description=(
        "Number of MD steps between protonation state change attempts in constant pH simulations. "
        "For implicit solvent (icnstph=1), one randomly selected residue is examined per MC step. "
        "For explicit solvent (icnstph=2), all titratable residues are examined in random order "
        "at each MC step."
    ),
    default=10,
    notes=(
        "For implicit solvent, decrease ntcnstph as the number of titrating residues increases "
        "to maintain a constant effective step period for each residue. Good results have been "
        "observed with approximately 100 fs effective period per residue (e.g., ntcnstph=5, "
        "dt=0.002 with about 10 residues titrating). For explicit solvent, protonation state "
        "change attempts are done less frequently (e.g., ntcnstph=100) since each attempt "
        "requires costly solvent relaxation dynamics."
    ),
    value_type="int",
    section="cntrl",
    category="Constant pH MD",
    related=["icnstph", "solvph", "ntrelax"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntrelax",
    description=(
        "Number of solvent relaxation MD steps performed following successful protonation state "
        "changes in explicit solvent constant pH simulations (icnstph=2). During relaxation, "
        "the non-solvent atoms are held fixed while the solvent is allowed to relax around the "
        "new protonation state."
    ),
    default=100,
    notes=(
        "While approximately 4 ps of relaxation is needed for a fully relaxed solvent distribution, "
        "200 fs (approximately 100 steps with dt=0.002) has been found sufficient to account for "
        "the bulk of solvent relaxation. Only used when icnstph=2 (explicit solvent CpHMD). "
        "Has no effect for implicit solvent simulations (icnstph=1)."
    ),
    value_type="int",
    section="cntrl",
    category="Constant pH MD",
    related=["icnstph", "ntcnstph"],
    commonly_changed=True,
))


KEYWORDS.append(Keyword(
    name="iphmd",
    description=(
        "Flag to turn on continuous constant pH molecular dynamics based on the lambda-dynamics "
        "method. Unlike the discrete CpHMD (icnstph), this method uses fictitious lambda particles "
        "with coordinates bound between 0 (protonated) and 1 (deprotonated) to represent protonation "
        "states. Lambda particles are propagated with a Langevin integrator using auxiliary theta "
        "variables (lambda = sin^2(theta)). Titratable groups with two competitive protonation sites "
        "(e.g., His, Asp, Glu) use an additional variable x to control tautomeric states."
    ),
    default=0,
    options={
        0: "No continuous constant pH MD (default).",
        1: (
            "Continuous CpHMD in implicit solvent using the GBNeck2 (igb=8) model. Both "
            "conformational and protonation state sampling are performed using GB. Implemented "
            "for both CPUs and GPUs (GPU version available via downloadable patch). Uses ff14SB "
            "with mbondi3 radii. Residues AS2/GL2 replace ASP/GLU, HIS must be HIP."
        ),
        2: (
            "Continuous CpHMD in hybrid solvent mode. Conformational sampling is performed in "
            "explicit solvent, while protonation state sampling (forces on lambda particles) uses "
            "a GB model. Implemented only for CPUs and has not been extensively tested."
        ),
        3: (
            "Continuous CpHMD in all-atom PME mode. Both conformational and protonation state "
            "sampling are performed in explicit solvent with no GB calculations. Implemented only "
            "for GPUs. Supports AMBER ff14SB, ff19SB, and CHARMM22 force fields."
        ),
    },
    notes=(
        "When iphmd=1 or 3, pmemd takes several additional command-line flags: -phmdin (phmd input "
        "file with &phmdin namelist), -phmdparm (parameter file with &phmdparm namelist containing "
        "model compound PMFs), -phmdrstrt (restart file output with &phmdrst namelist), -phmdout "
        "(lambda file output), and optionally -phmdstrt (restart file input with &phmdstrt namelist). "
        "The solution pH is set by solvph in &cntrl. For GB-CpHMD, load ff14SB, set PBRadii to "
        "mbondi3, and load phmd.lib and frcmod.phmd in LEaP. ASP/GLU residues should be AS2/GL2 "
        "and HIS should be HIP. Validated for Asp, Glu, His, Cys, and Lys titration. "
        "pH replica-exchange significantly accelerates sampling and is recommended with 4+ replicas."
    ),
    value_type="int",
    section="cntrl",
    category="Continuous Constant pH MD",
    related=["solvph", "icnstph", "igb"],
    commonly_changed=True,
))


# ---------------------------------------------------------------------------
# Chapter 28: Constant Redox Potential MD (&cntrl keywords)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="icnste",
    description=(
        "Flag to turn on constant Redox Potential molecular dynamics using discrete redox "
        "states with Monte Carlo sampling. Based on the same methodology as constant pH MD "
        "(icnstph) but applied to redox-active titratable residues. Requires a cein file "
        "specified via the -cein command-line flag. Can be run simultaneously with constant "
        "pH simulations (icnstph) for combined constant pH and Redox Potential MD. Residues "
        "can also be simultaneously pH- and redox-active using cpeinutil.py and a cpein file."
    ),
    default=0,
    options={
        0: "No constant Redox Potential MD (default).",
        1: (
            "Constant Redox Potential MD in implicit solvent (Generalized Born). Redox state "
            "changes are attempted every ntcnste steps using a Monte Carlo scheme. One residue "
            "is examined per MC step. The GB model and parameters must match those used to "
            "derive the reference energies in the cein file."
        ),
        2: (
            "Constant Redox Potential MD in explicit solvent. MD runs for ntcnstph steps, "
            "then solvent and ions are stripped, one redox state change attempt is performed "
            "for each redox-active residue in random order using a GB model, solvent is restored, "
            "and solvent relaxation dynamics run for ntrelaxe steps if any states changed."
        ),
    },
    notes=(
        "When icnste=1 or 2, you must also specify solve and ntcnste in &cntrl. "
        "The cein file is specified on the command line with -cein, the redox state history "
        "is written to -ceout, and a restart file is written to -cerestrt. "
        "Currently, Amber provides definitions for titrating a bis-histidine heme group (HEH residue) "
        "as in N-acetylmicroperoxidase-8 or horse heart cytochrome c. The heme propionates (PRN) "
        "are separate pH-active residues. Reference energies were derived with: "
        "cut=1000.0, igb=#, saltcon=0.1, nrespa=1, temp0=300.0, ntc=2, ntf=2 "
        "using the ff10 force field. Use leaprc.conste to load HEH, PRN, HIO, CYO residue "
        "libraries and force field modifications."
    ),
    value_type="int",
    section="cntrl",
    category="Constant Redox Potential MD",
    related=["solve", "ntcnste", "ntrelaxe", "icnstph", "saltcon", "igb"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="solve",
    description=(
        "Sets the solvent Redox Potential value for constant Redox Potential MD simulations, "
        "in units of Volts. This is an external parameter that affects the redox state "
        "distribution through the Metropolis Monte Carlo criteria, analogous to solvph for "
        "constant pH simulations."
    ),
    default=0.0,
    notes=(
        "The choice of solve determines the redox state equilibrium. For example, "
        "solve=-0.203 sets the solvent Redox Potential to -0.203 V. For E-REMD "
        "(Redox Potential Replica Exchange), different replicas use different solve values."
    ),
    value_type="float",
    section="cntrl",
    category="Constant Redox Potential MD",
    related=["icnste", "ntcnste"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntcnste",
    description=(
        "Number of MD steps between redox state change attempts in constant Redox Potential "
        "simulations. For implicit solvent (icnste=1), one randomly selected residue is examined "
        "per MC step. For explicit solvent (icnste=2), all redox-active residues are examined "
        "in random order at each MC step."
    ),
    default=10,
    notes=(
        "For implicit solvent, decrease ntcnste as the number of redox-active residues "
        "increases to maintain a constant effective step period for each residue. For explicit "
        "solvent, redox state change attempts are done less frequently (e.g., ntcnste=100) "
        "since each attempt requires costly solvent relaxation dynamics."
    ),
    value_type="int",
    section="cntrl",
    category="Constant Redox Potential MD",
    related=["icnste", "solve", "ntrelaxe"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ntrelaxe",
    description=(
        "Number of solvent relaxation MD steps performed following successful redox state "
        "changes in explicit solvent constant Redox Potential simulations (icnste=2). During "
        "relaxation, the non-solvent atoms are held fixed while the solvent is allowed to "
        "relax around the new redox state."
    ),
    default=200,
    notes=(
        "A value of 200 steps (200 fs with dt=0.002) is generally sufficient for solvent "
        "relaxation after redox state changes. Only used when icnste=2 (explicit solvent). "
        "Has no effect for implicit solvent simulations (icnste=1)."
    ),
    value_type="int",
    section="cntrl",
    category="Constant Redox Potential MD",
    related=["icnste", "ntcnste"],
    commonly_changed=True,
))

# =============================================================================
# &ewald NAMELIST (22.7.2 Particle Mesh Ewald)
# =============================================================================

KEYWORDS.append(Keyword(
    name="nfft1",
    description=(
        "Size of the charge grid in the X dimension (upon which the reciprocal sums are interpolated). "
        "Higher values lead to higher accuracy (when the DSUM_TOL is also lowered) but considerably slow "
        "the calculation. Generally reasonable results are obtained when NFFT1 is approximately equal to A "
        "(the box dimension), leading to a grid spacing of 1.0 A. Significant performance enhancement is "
        "obtained by having each integer value be a product of powers of 2, 3, and/or 5."
    ),
    default=None,
    notes="If the values are not given, the program will choose values to meet the criteria.",
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["nfft2", "nfft3", "order", "dsum_tol"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nfft2",
    description="Size of the charge grid in the Y dimension. See nfft1 for details.",
    default=None,
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["nfft1", "nfft3"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nfft3",
    description="Size of the charge grid in the Z dimension. See nfft1 for details.",
    default=None,
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["nfft1", "nfft2"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="order",
    description=(
        "The order of the B-spline interpolation. The higher the order, the better the accuracy (unless "
        "the charge grid is too coarse). The minimum order is 3. An order of 4 (the default) implies a "
        "cubic spline approximation which is a good standard value. Note that the cost of the PME goes "
        "as roughly the order to the third power."
    ),
    default=4,
    value_type="int",
    min_val=3,
    section="ewald",
    category="Particle Mesh Ewald",
    related=["nfft1", "nfft2", "nfft3"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="verbose",
    description=(
        "Standard use is to have VERBOSE = 0. Setting VERBOSE to higher values (up to a maximum of 3) "
        "leads to voluminous output of information about the PME run."
    ),
    default=0,
    value_type="int",
    min_val=0,
    max_val=3,
    section="ewald",
    category="Particle Mesh Ewald",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ew_type",
    description=(
        "Standard use is to have EW_TYPE = 0 which turns on the particle mesh ewald (PME) method. When "
        "EW_TYPE = 1, instead of the approximate, interpolated PME, a regular Ewald calculation is run."
    ),
    default=0,
    options={
        0: "Particle mesh Ewald (PME) method.",
        1: "Regular Ewald calculation.",
    },
    notes=(
        "The exact Ewald summation is present mainly to serve as an accuracy check. Although the cost "
        "of the exact Ewald method formally increases with system size at a much higher rate than the PME, "
        "it may be faster for small numbers of atoms (< 500)."
    ),
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["dsum_tol", "rsum_tol"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dsum_tol",
    description=(
        "This relates to the width of the direct sum part of the Ewald sum, requiring that the value of "
        "the direct sum at the Lennard-Jones cutoff value (specified in CUT) be less than DSUM_TOL. "
        "Standard values for DSUM_TOL are in the range of 10^-6 to 10^-5, leading to estimated RMS "
        "deviation force errors of 0.00001 to 0.0005."
    ),
    default=1e-5,
    value_type="float",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["cut", "ew_coeff"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="rsum_tol",
    description=(
        "This serves as a way to generate the number of reciprocal vectors used in an Ewald sum. "
        "Typically the relative RMS reciprocal sum error is about 5-10 times RSUM_TOL."
    ),
    default=5e-5,
    value_type="float",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["ew_type"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ew_coeff",
    description=(
        "Ewald coefficient, in A^-1. Default is determined by dsum_tol and cutoff. If it is explicitly "
        "input then that value is used, and dsum_tol is computed from ew_coeff and cutoff."
    ),
    default=None,
    notes="Determined automatically from dsum_tol and cut if not specified.",
    value_type="float",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["dsum_tol", "cut"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nbflag",
    description=(
        "If nbflag = 0, construct the direct sum nonbonded list in the 'old' way, i.e. update the list "
        "every nsnb steps. If nbflag = 1 (the default when imin = 0 or ntb > 0), nsnb is ignored, and "
        "the list is updated whenever any atom has moved more than 1/2 skinnb since the last list update."
    ),
    default=1,
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["nsnb", "skinnb"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="skinnb",
    description=(
        "Width of the nonbonded 'skin'. The direct sum nonbonded list is extended to cut + skinnb, and "
        "the van der Waals and direct electrostatic interactions are truncated at cut. Use of this "
        "parameter is required for energy conservation, and recommended for all PME runs."
    ),
    default=2.0,
    notes="Units: Angstroms.",
    value_type="float",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["cut", "nbflag"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="skin",
    description=(
        "(pmemd.cuda only) The threshold, as a fraction of skinnb, at which particle migration will "
        "trigger a non-bonded pair list rebuild. Enter values between 0.5 (minimum, default) and 1.0 "
        "(maximum). A setting of 0.75 is recommended for the best tradeoff of performance to safety."
    ),
    default=0.5,
    value_type="float",
    min_val=0.5,
    max_val=1.0,
    section="ewald",
    category="Particle Mesh Ewald",
    related=["skinnb"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nbtell",
    description=(
        "If nbtell = 1, a message is printed when any atom has moved far enough to trigger a list update. "
        "Use only for debugging or analysis."
    ),
    default=0,
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["nbflag"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="netfrc",
    description=(
        "Controls whether the net force arising from the PME mesh calculation is removed. Setting "
        "netfrc=0 gets the 'legacy' operation on pmemd.cuda, which did not calculate or remove the net "
        "force during PME simulations. The nscm setting can be used to remove any net momentum on a much "
        "less frequent time scale."
    ),
    default=1,
    notes=(
        "One should also set netfrc=0 when ntr>0. The performance cost of net force removal is a "
        "fraction of 1% of the total time."
    ),
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["ntr", "nscm"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="vdwmeth",
    description=(
        "Method for treating van der Waals interactions beyond the cutoff. When set to 1, a continuum "
        "model correction is added. When set to 0, no correction is used."
    ),
    default=1,
    options={
        0: "No long-range van der Waals correction.",
        1: "Use continuum model correction for long-range van der Waals (default for periodic).",
    },
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["cut"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="eedmeth",
    description=(
        "Method for the direct sum electrostatic interaction evaluation. Default value of 1 uses a "
        "cubic spline. A value of 2 implies a linear table lookup. A value of 3 implies use of an "
        "'exact' subroutine call."
    ),
    default=1,
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["eedtbdns"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="eedtbdns",
    description="Density of spline or linear lookup table, if eedmeth is 1 or 2.",
    default=500,
    notes="Units: points per unit.",
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=["eedmeth"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="column_fft",
    description=(
        "1 or 0 flag to turn on or off, respectively, column-mode fft for parallel runs. The default "
        "mode is slab mode which is efficient for low processor counts. The column method can be faster "
        "for larger processor counts."
    ),
    default=0,
    value_type="int",
    section="ewald",
    category="Particle Mesh Ewald",
    related=[],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.7.3 IPS (Isotropic Periodic Sum)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ips",
    description=(
        "Flag to control nonbonded interaction calculation method. When IPS is used for electrostatic "
        "interaction, PME will be turned off."
    ),
    default=0,
    options={
        0: "IPS will not be used (default).",
        1: "3D IPS will be used for both electrostatic and L-J interactions.",
        2: "3D IPS will be used only for electrostatic, including all multipole, interactions.",
        3: "3D IPS will be used only for L-J interactions.",
        4: "3D IPS/DFFT will be used for both electrostatic and L-J interactions.",
        5: "3D IPS/DFFT will be used only for electrostatic interactions.",
        6: "3D IPS/DFFT will be used only for L-J interactions.",
    },
    value_type="int",
    section="cntrl",
    category="Isotropic Periodic Sum",
    related=["raips", "cut"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="raips",
    description=(
        "Local region radius. raips is automatically set to cut for 3D IPS calculations (ips<=3) and "
        "should be set larger than cut for 3D IPS/DFFT calculations (ips>=4). A negative value indicates "
        "that it is set to the longest box side of a simulation system."
    ),
    default=-1.0,
    value_type="float",
    section="cntrl",
    category="Isotropic Periodic Sum",
    related=["ips", "cut"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.7.4 Extra point options (&ewald namelist)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="frameon",
    description=(
        "If frameon is set to 1, (default) the bonds, angles and dihedral interactions involving the "
        "lone pairs/extra points are removed except for constraints added during parm. The lone pairs "
        "are kept in ideal geometry relative to local atoms, and resulting torques are transferred to "
        "these atoms. To treat extra points as regular atoms, set frameon=0."
    ),
    default=1,
    value_type="int",
    section="ewald",
    category="Extra point options",
    related=["chngmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="chngmask",
    description=(
        "If chngmask=1 (default), new 1-1, 1-2, 1-3 and 1-4 interactions are calculated for extra "
        "points. An extra point belonging to an atom has a 1-1 interaction with it, and participates "
        "in any 1-2, 1-3 or 1-4 interaction that atom has."
    ),
    default=1,
    value_type="int",
    section="ewald",
    category="Extra point options",
    related=["frameon"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.7.5 Polarizable potentials (&ewald namelist)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="indmeth",
    description=(
        "If indmeth is 0, 1, or 2 then the nonbond force is called iteratively until successive "
        "estimates of the induced dipoles agree to within DIPTOL in the root mean square sense. "
        "If indmeth = 3, use a Car-Parinello scheme wherein dipoles are assigned a fictitious mass "
        "and integrated each time step. This is much more efficient and is the current default."
    ),
    default=3,
    notes="Method 3 is unstable for dt > 1 fs.",
    value_type="int",
    section="ewald",
    category="Polarizable potentials",
    related=["ipol", "diptol", "dipmass"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="diptol",
    description="Convergence criterion for dipoles in the iterative methods.",
    default=0.0001,
    notes="Units: Debye.",
    value_type="float",
    section="ewald",
    category="Polarizable potentials",
    related=["indmeth", "maxiter"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxiter",
    description=(
        "For iterative methods (indmeth<3), this is the maximum number of iterations allowed per "
        "time step."
    ),
    default=20,
    value_type="int",
    section="ewald",
    category="Polarizable potentials",
    related=["indmeth", "diptol"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipmass",
    description=(
        "The fictitious mass assigned to dipoles. Default value is 0.33, which works well for 1 fs time "
        "steps. If dipmass is set much below this, the dynamics are rapidly unstable. If set much above "
        "this the dynamics of the system are affected."
    ),
    default=0.33,
    value_type="float",
    section="ewald",
    category="Polarizable potentials",
    related=["indmeth"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="diptau",
    description=(
        "This is used for temperature control of the dipoles (for indmeth=3). If diptau is greater than "
        "10 (ps units) temperature control of dipoles is turned off."
    ),
    default=11.0,
    notes="Default is 11 ps (i.e. default is turned off).",
    value_type="float",
    section="ewald",
    category="Polarizable potentials",
    related=["indmeth"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="irstdip",
    description=(
        "If indmeth=3, a restart file for dipole positions and velocities is written along with the "
        "restart for atomic coordinates and velocities. If irstdip=1, the dipolar positions and "
        "velocities from the inpdip file are read in. If irstdip=0, an iterative method is used for "
        "step 1, after which Car-Parrinello is used."
    ),
    default=0,
    value_type="int",
    section="ewald",
    category="Polarizable potentials",
    related=["indmeth"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scaldip",
    description=(
        "To scale 1-4 charge-dipole and dipole-dipole interactions the same as 1-4 charge-charge "
        "(i.e. divided by scee) set scaldip=1 (default). If scaldip=0 the 1-4 charge-dipole and "
        "dipole-dipole interactions are treated the same as other dipolar interactions (i.e. divided by 1)."
    ),
    default=1,
    value_type="int",
    section="ewald",
    category="Polarizable potentials",
    related=["ipol"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# 22.7.7 Detailed MPI Timings
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="profile_mpi",
    description=(
        "Adjusts whether detailed per thread timings should be written to a file called profile_mpi "
        "when running sander in parallel. By default only average timings are printed to the output file."
    ),
    default=0,
    options={
        0: "No detailed MPI timings will be written (default).",
        1: "A detailed breakdown of the timings for each MPI thread will be written to the file: profile_mpi.",
    },
    value_type="int",
    section="ewald",
    category="MPI timings",
    related=[],
    commonly_changed=False,
))

# =============================================================================
# PMEMD-specific namelist variables (Chapter 23.3)
# =============================================================================

KEYWORDS.append(Keyword(
    name="mdout_flush_interval",
    description=(
        "In &cntrl, this variable can be used to control the minimum time in integer seconds between "
        "'flushes' of the mdout file. PMEMD does an open/close cycle on mdout at this interval. "
        "The default of 300 seconds provides a good compromise between efficiency and being able to "
        "observe the progress of the simulation."
    ),
    default=300,
    value_type="int",
    min_val=0,
    max_val=3600,
    section="cntrl",
    category="PMEMD-specific",
    related=["mdinfo_flush_interval"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mdinfo_flush_interval",
    description=(
        "In &cntrl, this variable can be used to control the minimum time in integer seconds between "
        "'flushes' of the mdinfo file. PMEMD does an open/close cycle on mdinfo at this interval."
    ),
    default=60,
    value_type="int",
    min_val=0,
    max_val=3600,
    section="cntrl",
    category="PMEMD-specific",
    related=["mdout_flush_interval"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="es_cutoff",
    description=(
        "In &cntrl, controls the cutoff for electrostatic direct force interactions in PME calculations "
        "separately from vdw. If you specify these variables, you should not specify the cut variable, "
        "and there is a requirement that vdw_cutoff >= es_cutoff. This feature is not currently available "
        "on the GPU."
    ),
    default=None,
    value_type="float",
    section="cntrl",
    category="PMEMD-specific",
    related=["vdw_cutoff", "cut"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="vdw_cutoff",
    description=(
        "In &cntrl, controls the cutoff for vdw interactions in PME calculations separately from "
        "electrostatics. If specified, vdw_cutoff must be >= es_cutoff. This feature is not currently "
        "available on the GPU."
    ),
    default=None,
    value_type="float",
    section="cntrl",
    category="PMEMD-specific",
    related=["es_cutoff", "cut"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="no_intermolecular_bonds",
    description=(
        "New variable controlling molecule definition. If 1, any molecules (as defined by the prmtop) "
        "joined by a covalent bond are fused to form a single molecule for purposes of pressure and "
        "virial-related operations; if 0 then the old behaviour (use prmtop molecule definitions) "
        "pertains."
    ),
    default=1,
    notes=(
        "A value of 0 is not supported with force-fields using extra points. If consistency with sander "
        "is more important and you are not using extra points, set to 0."
    ),
    value_type="int",
    section="cntrl",
    category="PMEMD-specific",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ene_avg_sampling",
    description=(
        "New variable controlling the number of steps between energy samples used in energy averages. "
        "If not specified, then ntpr is used (default). To match the behaviour of sander or PMEMD v9 "
        "or earlier, this variable should be set to 1."
    ),
    default=None,
    notes=(
        "This variable is only used for MD, not minimization and will also effectively be turned off "
        "if ntave is in use (non-0) or RESPA is in use (nrespa > 1). Performance can be improved "
        "without really losing anything of value by using the new default for energy average sampling."
    ),
    value_type="int",
    section="cntrl",
    category="PMEMD-specific",
    related=["ntpr", "ntave", "nrespa"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_axis_opt",
    description=(
        "In &ewald. For parallel runs, the most favorable orientation of an orthogonal unit cell is "
        "with the longest side in the Z direction. Axis optimization is only done for mpi runs in which "
        "an orthogonal unit cell has an aspect ratio of at least 3 to 2. It is turned off for all "
        "minimization runs and for runs in which velocities are randomized (ntt = 2 or 3)."
    ),
    default=None,
    notes="Default adapts to run conditions. Set to 1 to force on, 0 to force off.",
    value_type="int",
    section="ewald",
    category="PMEMD-specific",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="fft_grids_per_ang",
    description=(
        "In &ewald. This variable may be used to set the desired reciprocal space fft grid density in "
        "terms of fft grids/angstrom. The nearest grid dimensions that meet or exceed this density will "
        "be used (i.e., nfft1,2,3 are set based on this specification)."
    ),
    default=1.0,
    value_type="float",
    section="ewald",
    category="PMEMD-specific",
    related=["nfft1", "nfft2", "nfft3"],
    commonly_changed=False,
))


# =============================================================================
# &pol_gauss NAMELIST (Section 22.8 - Polarizable Gaussian Multipole Model)
# =============================================================================

KEYWORDS.append(Keyword(
    name="pol_gauss_verbose",
    description=(
        "In addition to the usual sander output, by setting pol_gauss_verbose=1, extra printing of energy "
        "and forces can be found in the output file."
    ),
    default=0,
    options={
        0: "Normal output only.",
        1: "Extra printing of energy and forces in output file.",
    },
    value_type="int",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["ipgm"],
))

KEYWORDS.append(Keyword(
    name="pol_gauss_ips",
    description=(
        "To use isotropic periodic sum method for pGM model, set pol_gauss_ips=1. "
        "ips=1 under &cntrl Namelist is also required."
    ),
    default=0,
    options={
        0: "Do not use IPS for pGM.",
        1: "Use isotropic periodic sum method for pGM model.",
    },
    notes="Requires ips=1 in &cntrl namelist when enabled.",
    value_type="int",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["ipgm", "ips"],
))

KEYWORDS.append(Keyword(
    name="ee_dsum_cut",
    description=(
        "The Ewald direct sum cutoff for pGM. It is recommended to be set to at least 9 Angstroms."
    ),
    default=None,
    notes=(
        "pGM model requires higher accuracy than classic point charge model. "
        "Recommended to be at least 9 Angstroms."
    ),
    value_type="float",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["ew_coeff", "cut"],
))

KEYWORDS.append(Keyword(
    name="dipole_scf_tol",
    description=(
        "The induced dipoles in the pGM model are solutions to a set of linear equations. "
        "These equations are solved iteratively by a linear system solver. dipole_scf_tol is "
        "the convergence criterion for the iterative solution to the linear equations. To achieve "
        "good energy conservation in NVE simulations (i.e. similar to that observed for additive "
        "force fields at otherwise identical conditions), a convergence criterion of 1e-2 is needed. "
        "Starting from 2023, the convergence is measured with the maximum relative error on individual "
        "dipoles instead of overall residue relative error, so the numerical tolerance appears to be "
        "very different, but the convergence quality requirement is similar to previous releases."
    ),
    default=None,
    notes=(
        "For practical simulations, a criterion of 1e-2 achieves satisfactory energy conservation "
        "while significantly reducing computational time. For strict NVE, use 1e-4 or tighter."
    ),
    value_type="float",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["ipgm"],
))

KEYWORDS.append(Keyword(
    name="dipole_solv_opt",
    description=(
        "Set the induction iteration solver."
    ),
    default=3,
    options={
        3: "Preconditioned conjugate gradient solver (recommended).",
        4: "SOR solver (previous default).",
    },
    value_type="int",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["scf_cg_niter", "scf_sor_niter", "scf_sor_coefficient"],
))

KEYWORDS.append(Keyword(
    name="scf_cg_niter",
    description=(
        "The maximum iterations when solving the induction equations with a conjugate gradient solver."
    ),
    default=50,
    value_type="int",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["dipole_solv_opt"],
))

KEYWORDS.append(Keyword(
    name="scf_sor_coefficient",
    description=(
        "This is the successive relaxation parameter in the SOR solver, which can be adjusted "
        "to balance the efficiency and stability of the solver."
    ),
    default=0.65,
    value_type="float",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["dipole_solv_opt", "scf_sor_niter"],
))

KEYWORDS.append(Keyword(
    name="scf_sor_niter",
    description=(
        "The maximum iterations when solving the induction equations with the SOR solver."
    ),
    default=100,
    value_type="int",
    section="pol_gauss",
    category="Polarizable Gaussian Multipole",
    related=["dipole_solv_opt", "scf_sor_coefficient"],
))


# =============================================================================
# &debugf NAMELIST (Section 22.11 - Debugging information)
# =============================================================================

KEYWORDS.append(Keyword(
    name="do_debugf",
    description=(
        "Flag to perform this module. Set to one to turn on debug options. "
        "If the debug options are set, sander will exit after performing the debug tasks set by the user."
    ),
    default=0,
    options={
        0: "Debug module off.",
        1: "Turn on debug options.",
    },
    value_type="int",
    section="debugf",
    category="Debug",
))

KEYWORDS.append(Keyword(
    name="atomn",
    description=(
        "Array of atom numbers to test atomic forces on. Up to 25 atom numbers can be specified, "
        "separated by commas."
    ),
    default=None,
    notes="Up to 25 atom numbers can be specified.",
    value_type="int_array",
    section="debugf",
    category="Debug",
    related=["nranatm", "do_debugf"],
))

KEYWORDS.append(Keyword(
    name="nranatm",
    description=(
        "Number of random atoms to test atomic forces on. Atom numbers are generated via a "
        "random number generator."
    ),
    default=0,
    value_type="int",
    section="debugf",
    category="Debug",
    related=["ranseed", "atomn", "do_debugf"],
))

KEYWORDS.append(Keyword(
    name="ranseed",
    description=(
        "Seed of random number generator used in generating atom numbers."
    ),
    default=71277,
    value_type="int",
    section="debugf",
    category="Debug",
    related=["nranatm"],
))

KEYWORDS.append(Keyword(
    name="neglgdel",
    description=(
        "Negative log of delta used in numerical differentiating; e.g. 4 means delta is 1e-4 Angstroms."
    ),
    default=5,
    notes=(
        "In general it does no good to set neglgdel larger than about 6. This is because the relative "
        "force error is at best the square root of the numerical error in the energy, which ranges from "
        "1e-15 up to 1e-12 for energies involving a large number of terms."
    ),
    value_type="int",
    section="debugf",
    category="Debug",
    related=["do_debugf"],
))

KEYWORDS.append(Keyword(
    name="chkvir",
    description=(
        "Flag to test the atomic and molecular virials numerically."
    ),
    default=0,
    options={
        0: "Do not test virials.",
        1: "Test virials numerically.",
    },
    value_type="int",
    section="debugf",
    category="Debug",
    related=["do_debugf"],
))

KEYWORDS.append(Keyword(
    name="dumpfrc",
    description=(
        "Flag to dump energies, forces and virials, as well as components of forces (bond, angle forces etc.) "
        "to the file 'forcedump.dat'. This produces an ASCII file."
    ),
    default=0,
    options={
        0: "Do not dump forces.",
        1: "Dump forces to forcedump.dat.",
    },
    value_type="int",
    section="debugf",
    category="Debug",
    related=["rmsfrc", "do_debugf"],
))

KEYWORDS.append(Keyword(
    name="rmsfrc",
    description=(
        "Flag to compare energies, forces and virials as well as components of forces (bond, angle forces "
        "etc.) to those in the file 'forcedump.dat'."
    ),
    default=0,
    options={
        0: "Do not compare forces.",
        1: "Compare forces to those in forcedump.dat.",
    },
    notes="Typically used to get an RMS force error for the Ewald method in use.",
    value_type="int",
    section="debugf",
    category="Debug",
    related=["dumpfrc", "do_debugf"],
))

KEYWORDS.append(Keyword(
    name="zerochg",
    description=(
        "Flag to zero all charges before calculating forces."
    ),
    default=0,
    options={
        0: "Keep charges.",
        1: "Remove all charges.",
    },
    value_type="int",
    section="debugf",
    category="Debug",
    related=["zerovdw", "zerodip", "do_debugf"],
))

KEYWORDS.append(Keyword(
    name="zerovdw",
    description=(
        "Flag to remove all van der Waals interactions before calculating forces."
    ),
    default=0,
    options={
        0: "Keep van der Waals interactions.",
        1: "Remove all van der Waals interactions.",
    },
    value_type="int",
    section="debugf",
    category="Debug",
    related=["zerochg", "zerodip", "do_debugf"],
))

KEYWORDS.append(Keyword(
    name="zerodip",
    description=(
        "Flag to remove all atomic dipoles before calculating forces. Only relevant when "
        "polarizability is invoked."
    ),
    default=0,
    options={
        0: "Keep atomic dipoles.",
        1: "Remove all atomic dipoles.",
    },
    value_type="int",
    section="debugf",
    category="Debug",
    related=["zerochg", "zerovdw", "do_debugf"],
))

# do_dir, do_rec, do_adj, do_self, do_bond, do_cbond, do_angle, do_ephi, do_xconst, do_cap
for _flag_name, _flag_desc in [
    ("do_dir", "direct sum interactions (van der Waals as well as electrostatic)"),
    ("do_rec", "reciprocal sum interactions"),
    ("do_adj", "adjusted (1-4 and excluded) interactions"),
    ("do_self", "self-energy interactions"),
    ("do_bond", "bond energy subroutine"),
    ("do_cbond", "constraint bond subroutine"),
    ("do_angle", "angle energy subroutine"),
    ("do_ephi", "dihedral energy subroutine"),
    ("do_xconst", "extra constraints subroutine"),
    ("do_cap", "water cap subroutine"),
]:
    KEYWORDS.append(Keyword(
        name=_flag_name,
        description=(
            f"Flag to turn on or off the {_flag_desc}. "
            f"Set to zero to prevent the subroutine from running. "
            f"These options, as well as the zerochg, zerovdw, zerodip flags, can be used "
            f"to fine tune a test of forces, accuracy, etc."
        ),
        default=1,
        options={
            0: f"Turn off {_flag_desc}.",
            1: f"Turn on {_flag_desc} (default).",
        },
        value_type="int",
        section="debugf",
        category="Debug",
        related=["do_debugf", "zerochg", "zerovdw"],
    ))

# =============================================================================
# &qmmm NAMELIST (Section 22.14 / Chapter 9 / 11.1.6 - QM/MM)
# =============================================================================

# --- Float parameters ---

KEYWORDS.append(Keyword(
    name="qmcut",
    description=(
        "Nonbonded cutoff in Angstroms used for QM/MM nonbonded interactions. Note there is no "
        "such thing as a cutoff within the QM region, since it is the wavefunction of the entire "
        "system we are optimizing. The default value is the MM cutoff being used (i.e., cut from "
        "sander input)."
    ),
    default=None,
    notes="Defaults to the value of 'cut' from &cntrl.",
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["cut", "ifqnt"],
))

KEYWORDS.append(Keyword(
    name="lnk_dis",
    description="Distance in Angstroms of the QM atom to its link atom.",
    default=1.09,
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["lnk_atomic_no", "lnk_method"],
))

KEYWORDS.append(Keyword(
    name="scfconv",
    description=(
        "Controls the convergence of the SCF calculation. The SCF terminates when the energy "
        "difference between the last two steps is smaller than the value given here."
    ),
    default=1e-8,
    notes="The smallest value that can practically be used within limits of double precision is 1e-14.",
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["errconv", "itrmax", "tight_p_conv"],
))

KEYWORDS.append(Keyword(
    name="errconv",
    description=(
        "SCF tolerance on the maximum absolute value of the error matrix (i.e., the commutator "
        "of the Fock matrix with the density matrix). The value is in units of Hartrees. The default "
        "value is large enough that scfconv will always be more strict."
    ),
    default=None,
    notes="Units are Hartrees. Default is set large so scfconv controls convergence.",
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["scfconv"],
))

KEYWORDS.append(Keyword(
    name="dftb_telec",
    description=(
        "Electronic temperature, in K, used to accelerate SCC convergence in DFTB calculations. "
        "The electronic temperature affects the Fermi distribution promoting some HOMO/LUMO mixing, "
        "which can accelerate the convergence in difficult cases. In most cases, a low telec (around "
        "100K) is enough. Should be used only when necessary, and the results checked carefully."
    ),
    default=0.0,
    value_type="float",
    section="qmmm",
    category="QM/MM DFTB",
    related=["dftb_telec_step", "dftb_maxiter", "dftb_disper"],
))

KEYWORDS.append(Keyword(
    name="dftb_telec_step",
    description=(
        "The size of the step to take when reducing the electronic temperature in a DFTB calculation. "
        "The smaller the step, the longer it will take to get the electronic temperature to zero."
    ),
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM DFTB",
    related=["dftb_telec"],
))

KEYWORDS.append(Keyword(
    name="fockp_d1",
    description="First prefactor for the Fock matrix prediction. Changing this is not recommended.",
    default=2.4,
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
    related=["fockp_d2", "fockp_d3", "fockp_d4", "fock_predict"],
))

KEYWORDS.append(Keyword(
    name="fockp_d2",
    description="Second prefactor for the Fock matrix prediction. Changing this is not recommended.",
    default=-1.2,
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
    related=["fockp_d1", "fockp_d3", "fockp_d4", "fock_predict"],
))

KEYWORDS.append(Keyword(
    name="fockp_d3",
    description="Third prefactor for the Fock matrix prediction. Changing this is not recommended.",
    default=-0.8,
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
    related=["fockp_d1", "fockp_d2", "fockp_d4", "fock_predict"],
))

KEYWORDS.append(Keyword(
    name="fockp_d4",
    description="Fourth prefactor for the Fock matrix prediction. Changing this is not recommended.",
    default=0.6,
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
    related=["fockp_d1", "fockp_d2", "fockp_d3", "fock_predict"],
))

KEYWORDS.append(Keyword(
    name="damp",
    description="SCF damping factor. Changing this is not recommended.",
    default=1.0,
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
))

KEYWORDS.append(Keyword(
    name="vshift",
    description=(
        "Controls level shifting for NDDO methods (not DFTB). Virtual orbitals can be shifted up "
        "by vshift (in eV) to improve SCF convergence in cases with a small HOMO/LUMO gap."
    ),
    default=0.0,
    notes="Units are eV. Not applicable to DFTB methods.",
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
))

KEYWORDS.append(Keyword(
    name="kappa",
    description=(
        "Related to the Debye salt concentration for GB models. This is set automatically from "
        "saltcon in the sander input data structure."
    ),
    default=None,
    notes="Set automatically from saltcon; normally should not be set manually.",
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["qmgb"],
))

KEYWORDS.append(Keyword(
    name="pseudo_diag_criteria",
    description=(
        "Controls whether a pseudo-diagonalization of the Fock matrix can be performed "
        "(not applicable for DFTB)."
    ),
    default=0.05,
    value_type="float",
    section="qmmm",
    category="QM/MM SCF",
    related=["pseudo_diag", "diag_routine"],
))

KEYWORDS.append(Keyword(
    name="min_heavy_mass",
    description=(
        "The smallest value, in atomic mass units, that an atomic mass can have and still be "
        "considered a 'heavy-atom' (i.e., anything besides Hydrogen)."
    ),
    default=4.0,
    value_type="float",
    section="qmmm",
    category="QM/MM",
))

KEYWORDS.append(Keyword(
    name="r_switch_hi",
    description=(
        "If qmmm_switch is turned on, this is the distance, in Angstroms, at which the switch "
        "goes to zero. By default, it is the same as qmcut."
    ),
    default=None,
    notes="Defaults to value of qmcut.",
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["r_switch_lo", "qmmm_switch", "qmcut"],
))

KEYWORDS.append(Keyword(
    name="r_switch_lo",
    description=(
        "If qmmm_switch is turned on, this is the distance, in Angstroms, at which the switch "
        "turns on. By default, it is 2 Angstroms smaller than r_switch_hi."
    ),
    default=None,
    notes="Defaults to r_switch_hi - 2.0 Angstroms.",
    value_type="float",
    section="qmmm",
    category="QM/MM",
    related=["r_switch_hi", "qmmm_switch", "qmcut"],
))

# --- Integer parameters ---

KEYWORDS.append(Keyword(
    name="qmgb",
    description="Specifies how the QM region should be treated with Generalized Born.",
    default=2,
    options={
        2: (
            "The electrostatic and 'polarization' fields from the MM charges and the exterior "
            "dielectric, respectively, are included in the Fock matrix for the QM Hamiltonian."
        ),
        3: (
            "Intended for debugging and only useful for single-point calculations. Computes "
            "the GB energy by treating every atom in the QM region as a point charge equal to "
            "its Mulliken charge."
        ),
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["igb", "kappa"],
))

KEYWORDS.append(Keyword(
    name="lnk_atomic_no",
    description="The atomic number of the element you wish to use as the link atom.",
    default=1,
    notes="Default is 1 (Hydrogen).",
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["lnk_dis", "lnk_method"],
))

KEYWORDS.append(Keyword(
    name="ndiis_matrices",
    description="The number of error vectors to use for the DIIS convergence algorithm.",
    default=6,
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["ndiis_attempts"],
))

KEYWORDS.append(Keyword(
    name="ndiis_attempts",
    description=(
        "The number of iterations that DIIS extrapolation will be attempted. Not available for DFTB."
    ),
    default=0,
    notes="Maximum is 1000. Not available for DFTB.",
    max_val=1000,
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["ndiis_matrices"],
))

KEYWORDS.append(Keyword(
    name="lnk_method",
    description=(
        "The method used to define how classical valence terms across the QM/MM boundary "
        "will be treated. See Subsection 11.1.7 for more information."
    ),
    default=1,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["lnk_dis", "lnk_atomic_no"],
))

KEYWORDS.append(Keyword(
    name="qmcharge",
    description="The net charge of the QM region.",
    default=0,
    commonly_changed=True,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["corecharge", "buffercharge", "spin"],
))

KEYWORDS.append(Keyword(
    name="corecharge",
    description="The net charge of the core QM region.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["qmcharge", "buffercharge", "coremask", "core_iqmatoms"],
))

KEYWORDS.append(Keyword(
    name="buffercharge",
    description="The net charge of the buffer QM region.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["qmcharge", "corecharge", "buffermask", "buffer_iqmatoms"],
))

KEYWORDS.append(Keyword(
    name="spin",
    description="Spin multiplicity of the QM region.",
    default=1,
    notes="Default is 1 (singlet).",
    commonly_changed=True,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["qmcharge"],
))

KEYWORDS.append(Keyword(
    name="qmqmdx",
    description=(
        "Controls whether QM-QM derivatives are computed analytically or pseudo-numerically. "
        "The default (and recommended) is to use analytical QM-QM derivatives."
    ),
    default=1,
    options={
        1: "Analytical derivatives (recommended).",
        2: "Pseudo-numerical derivatives.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
))

KEYWORDS.append(Keyword(
    name="verbosity",
    description="Controls verbosity of QM/MM output.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Output",
))

KEYWORDS.append(Keyword(
    name="printcharges",
    description="Controls printing of QM Mulliken charges.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Output",
    related=["printdipole", "print_eigenvalues", "printbondorders"],
))

KEYWORDS.append(Keyword(
    name="printdipole",
    description="Controls printing of QM dipole moment.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Output",
    related=["printcharges"],
))

KEYWORDS.append(Keyword(
    name="print_eigenvalues",
    description="Controls printing of QM eigenvalues.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Output",
    related=["printcharges"],
))

KEYWORDS.append(Keyword(
    name="peptide_corr",
    description=(
        "If set to 0 (default), do not apply a correction to peptide linkages. "
        "If set to 1, apply a MM correction to peptide linkages."
    ),
    default=0,
    options={
        0: "No correction to peptide linkages.",
        1: "Apply MM correction to peptide linkages.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
))

KEYWORDS.append(Keyword(
    name="itrmax",
    description="Maximum number of SCF iterations to perform before deciding that the convergence has failed.",
    default=1000,
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["scfconv", "errconv"],
))

KEYWORDS.append(Keyword(
    name="printbondorders",
    description="Controls printing of QM bond orders.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Output",
    related=["printcharges"],
))

KEYWORDS.append(Keyword(
    name="qmshake",
    description=(
        "Controls whether SHAKE is applied to QM atoms. If 0, no SHAKE. If 1 (default), "
        "SHAKE QM atoms if MM SHAKE is turned on."
    ),
    default=1,
    options={
        0: "No SHAKE on QM atoms.",
        1: "SHAKE QM atoms if MM SHAKE (ntc) is turned on.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["ntc"],
))

KEYWORDS.append(Keyword(
    name="qmmmrij_incore",
    description=(
        "If set to 1 (default), store QM-MM pairs and related equations in memory. "
        "If set to 0, do not."
    ),
    default=1,
    options={
        0: "Do not store QM-MM pairs in memory (calculate on-the-fly).",
        1: "Store QM-MM pairs and related equations in memory.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["qmqm_erep_incore"],
))

KEYWORDS.append(Keyword(
    name="qmqm_erep_incore",
    description=(
        "If set to 1 (default), store QM-QM 1-electron repulsion integrals to memory. "
        "If set to 0, calculate them on-the-fly."
    ),
    default=1,
    options={
        0: "Calculate QM-QM 1-electron repulsion integrals on-the-fly.",
        1: "Store QM-QM 1-electron repulsion integrals in memory.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["qmmmrij_incore"],
))

KEYWORDS.append(Keyword(
    name="pseudo_diag",
    description=(
        "If set to 1 (default), allow the use of pseudo-diagonalization of the Fock matrix "
        "as long as the pseudo_diag_criteria is met."
    ),
    default=1,
    options={
        0: "Do not use pseudo-diagonalization.",
        1: "Allow pseudo-diagonalization if criteria met.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["pseudo_diag_criteria", "diag_routine"],
))

KEYWORDS.append(Keyword(
    name="qm_ewald",
    description=(
        "Specifies how long-range electrostatics for the QM region should be treated."
    ),
    default=None,
    options={
        0: "Use a real-space cutoff for QM-QM and QM-MM long-range interactions. QM atoms do not see their images and QM-MM interactions are truncated at the cutoff. Default for non-periodic simulations.",
        1: "Use PME or Ewald sum for long-range QM-QM and QM-MM electrostatics (default for periodic with PME).",
        2: "Similar to 1 but QM image charges are fixed at Mulliken charges from the previous MD step. Faster SCF convergence with minor loss of accuracy. Not extensively tested.",
    },
    notes="For non-periodic/GB simulations the default is 0. For periodic PME simulations the default is 1.",
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["qm_pme", "kmaxqx", "kmaxqy", "kmaxqz", "ksqmaxsq"],
))

KEYWORDS.append(Keyword(
    name="qm_pme",
    description=(
        "If 0, use a regular Ewald sum for computing QM-QM and QM-MM long-range electrostatic "
        "interactions. If 1 (default), use PME instead."
    ),
    default=1,
    options={
        0: "Use regular Ewald sum for QM long-range electrostatics.",
        1: "Use PME for QM long-range electrostatics.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["qm_ewald"],
))

KEYWORDS.append(Keyword(
    name="kmaxqx",
    description="Number of K-space vectors to use in the Ewald/PME calculations in the X-dimension.",
    default=8,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["kmaxqy", "kmaxqz", "ksqmaxsq", "qm_ewald"],
))

KEYWORDS.append(Keyword(
    name="kmaxqy",
    description="Same as kmaxqx, but in the Y-dimension.",
    default=8,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["kmaxqx", "kmaxqz", "ksqmaxsq"],
))

KEYWORDS.append(Keyword(
    name="kmaxqz",
    description="Same as kmaxqx, but in the Z-dimension.",
    default=8,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["kmaxqx", "kmaxqy", "ksqmaxsq"],
))

KEYWORDS.append(Keyword(
    name="ksqmaxsq",
    description=(
        "Specifies the maximum number of K^2 values for the spherical cutoff in reciprocal space "
        "when doing a QM-MM Ewald sum. The default value of 100 should be optimal for most systems."
    ),
    default=100,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["kmaxqx", "kmaxqy", "kmaxqz", "qm_ewald"],
))

KEYWORDS.append(Keyword(
    name="qmmm_int",
    description=(
        "Controls the way in which QM/MM electrostatic interactions are handled in the direct-space "
        "sum. VDW interactions are always calculated classically using the standard 6-12 potential."
    ),
    default=1,
    options={
        0: "Turn off all electrostatic interaction between QM and MM atoms in the direct-space sum. QM-MM VDW interactions still calculated classically.",
        1: "QM-MM interactions calculated analogously to QM core-core interaction. MM RESP charges included in one-electron Hamiltonian (default).",
        2: "Same as 1 but includes extra Gaussian terms from AM1/PM3 for QM-MM core-core repulsion (as in CHARMM/DYNAMO). Slightly reduces repulsion at small distances.",
        3: "Reformulated QM core-MM charge potential (PM3/MM*). Requires qm_theory=PM3. QM region limited to H, C, N, O atoms.",
        4: "Currently not in use.",
        5: "Mechanical embedding. QM-MM electrostatic interaction treated classically using force field point charges for QM atoms. Electron density not polarized by MM environment. Does not work with GB.",
    },
    notes="With the exception of qmmm_int=0, DFTB calculations always use Mulliken charge - RESP charge interaction regardless of qmmm_int setting.",
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["qm_theory", "qmcut"],
))

KEYWORDS.append(Keyword(
    name="adjust_q",
    description=(
        "Controls how charge is conserved during a QM/MM calculation with respect to link atoms. "
        "See Subsection 11.1.6 for more information."
    ),
    default=2,
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["lnk_method", "lnk_dis"],
))

KEYWORDS.append(Keyword(
    name="tight_p_conv",
    description=(
        "Controls the tightness of the convergence criteria on the density matrix in the SCF. "
        "If 0 (default), the convergence is loose. If set to 1, convergence is tight."
    ),
    default=0,
    options={
        0: "Loose convergence on density matrix.",
        1: "Tight convergence on density matrix.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["scfconv", "errconv"],
))

KEYWORDS.append(Keyword(
    name="diag_routine",
    description=(
        "The diagonalization routine to use to diagonalize the Fock matrix. "
        "By default (diag_routine=0), the fastest routine is chosen."
    ),
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["pseudo_diag"],
))

KEYWORDS.append(Keyword(
    name="density_predict",
    description=(
        "If 1, use the density matrix from the previous MD step. "
        "If 0 (default), do not predict the density matrix."
    ),
    default=0,
    options={0: "Do not predict density matrix.", 1: "Use density matrix from previous MD step."},
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["fock_predict"],
))

KEYWORDS.append(Keyword(
    name="fock_predict",
    description="If set to 0, do not attempt to predict the Fock matrix (default). If set to 1, try to.",
    default=0,
    options={0: "Do not predict Fock matrix.", 1: "Attempt to predict Fock matrix."},
    value_type="int",
    section="qmmm",
    category="QM/MM SCF",
    related=["density_predict", "fockp_d1", "fockp_d2", "fockp_d3", "fockp_d4"],
))

KEYWORDS.append(Keyword(
    name="vsolv",
    description="If set to 1, use variable solvent QM/MM. If set to 0 (default), do not.",
    default=0,
    options={0: "Standard QM/MM (no variable solvent).", 1: "Use variable solvent QM/MM."},
    value_type="int",
    section="qmmm",
    category="QM/MM",
))

KEYWORDS.append(Keyword(
    name="dftb_maxiter",
    description="The maximum number of SCF iterations to be used in SCC-DFTB calculations.",
    default=70,
    value_type="int",
    section="qmmm",
    category="QM/MM DFTB",
    related=["dftb_telec", "dftb_disper"],
))

KEYWORDS.append(Keyword(
    name="dftb_disper",
    description="If set to 1, use a dispersion correction for DFTB/SCC-DFTB. If set to 0 (default), do not.",
    default=0,
    options={0: "No dispersion correction for DFTB.", 1: "Use dispersion correction for DFTB."},
    value_type="int",
    section="qmmm",
    category="QM/MM DFTB",
    related=["dftb_maxiter", "dftb_telec"],
))

KEYWORDS.append(Keyword(
    name="dftb_chg",
    description="Controls printing of DFTB charges.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM DFTB",
))

KEYWORDS.append(Keyword(
    name="abfqmmm",
    description="Toggles the adaptive biased force QM/MM.",
    default=0,
    options={0: "Disable adaptive biased force QM/MM.", 1: "Enable adaptive biased force QM/MM."},
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["hot_spot", "coremask", "buffermask"],
))

KEYWORDS.append(Keyword(
    name="hot_spot",
    description=(
        "If set to 1, activates hot spot-like adaptive calculation in which the forces of atoms "
        "in the buffer region are linear combinations of the forces obtained from the extended and "
        "reduced calculations using a smoothing function."
    ),
    default=0,
    options={0: "Disable hot spot adaptive calculation.", 1: "Enable hot spot adaptive calculation."},
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "coremask", "buffermask"],
))

KEYWORDS.append(Keyword(
    name="qmmm_switch",
    description="If set to 1, use a switching function defined by r_switch_lo and r_switch_hi.",
    default=0,
    options={0: "No switching function.", 1: "Use switching function."},
    value_type="int",
    section="qmmm",
    category="QM/MM",
    related=["r_switch_lo", "r_switch_hi", "qmcut"],
))

# --- String / mask parameters ---

KEYWORDS.append(Keyword(
    name="qmmask",
    description="An Amber selection mask that provides a way of defining the QM region instead of iqmatoms.",
    default="",
    commonly_changed=True,
    notes="Character array, max 8192 characters.",
    value_type="string",
    section="qmmm",
    category="QM/MM",
    related=["iqmatoms", "qmcharge", "spin", "qm_theory"],
))

KEYWORDS.append(Keyword(
    name="iqmatoms",
    description=(
        "List of atom indexes, starting from 1, that will be treated using QM. "
        "This is one way, along with qmmask, of specifying the QM region."
    ),
    default=None,
    notes="Integer array, MAX_QUANTUM_ATOMS (10000).",
    value_type="int_array",
    section="qmmm",
    category="QM/MM",
    related=["qmmask", "qmcharge"],
))

KEYWORDS.append(Keyword(
    name="coremask",
    description=(
        "An Amber selection mask that provides a way of defining the core QM region "
        "in adaptive simulations instead of core_iqmatoms."
    ),
    default="",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["core_iqmatoms", "corecharge", "buffermask", "centermask"],
))

KEYWORDS.append(Keyword(
    name="buffermask",
    description=(
        "An Amber selection mask that provides a way of defining the buffer QM region "
        "in adaptive simulations instead of buffer_iqmatoms."
    ),
    default="",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["buffer_iqmatoms", "buffercharge", "coremask"],
))

KEYWORDS.append(Keyword(
    name="centermask",
    description=(
        "An Amber selection mask that defines the center region. If not set, it defaults to coremask."
    ),
    default="",
    notes="Defaults to coremask if not set.",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["coremask", "buffermask"],
))

KEYWORDS.append(Keyword(
    name="core_iqmatoms",
    description=(
        "A list of atom indices (starting at 1) selected for inclusion in the core "
        "QM/MM region in adaptive simulations."
    ),
    default=None,
    value_type="int_array",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["coremask", "buffer_iqmatoms"],
))

KEYWORDS.append(Keyword(
    name="buffer_iqmatoms",
    description=(
        "A list of atom indices (starting at 1) selected for inclusion in the buffer "
        "QM/MM region in adaptive simulations."
    ),
    default=None,
    value_type="int_array",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["buffermask", "core_iqmatoms"],
))

KEYWORDS.append(Keyword(
    name="dftb_3rd_order",
    description=(
        "Specifies the 3rd-order DFTB correction. Default ('NONE') means no 3rd order "
        "correction is used. See Chapter 9 for more information."
    ),
    default="NONE",
    value_type="string",
    section="qmmm",
    category="QM/MM DFTB",
    related=["dftb_maxiter", "dftb_disper"],
))

KEYWORDS.append(Keyword(
    name="qm_theory",
    description=(
        "String that defines which level of QM theory to use for the QM region. "
        "Controls which Hamiltonian is applied. Must be specified for any QM/MM calculation."
    ),
    default=None,
    commonly_changed=True,
    options={
        "PM3": "PM3 semiempirical Hamiltonian (default if set).",
        "AM1": "AM1 semiempirical Hamiltonian.",
        "MNDO": "MNDO semiempirical Hamiltonian.",
        "PM3-PDDG": "PM3/PDDG semiempirical Hamiltonian.",
        "AM1-PDDG": "AM1/PDDG semiempirical Hamiltonian.",
        "AM1/d": "AM1/d Hamiltonian (d-orbital extension for phosphorus).",
        "AM1-D*": "AM1-D* Hamiltonian with dispersion correction.",
        "PM6": "PM6 semiempirical Hamiltonian.",
        "DFTB2": "SCC-DFTB (DFTB2) density functional tight binding.",
        "DFTB3": "DFTB3 with third-order corrections.",
        "EXTERN": "Use external QM program via file-based interface (requires &adf, &gms, &gau, &orc, &qc, &mrcc, or &tc namelist).",
        "QUICK": "Use QUICK QM code via the API (linked library, recommended for HF/DFT).",
        "TERACHEM": "Use TeraChem via TCPB client/server interface (requires &tc namelist).",
        "XTB": "Use xTB tight-binding code (requires compiled xTB library linked to sander).",
        "DFTBPLUS": "Use DFTB+ code (requires compiled DFTB+ library linked to sander).",
        "SEBOMD": "Full semiempirical Born-Oppenheimer MD (all atoms QM, requires &sebomd namelist).",
    },
    notes=(
        "Must be specified for any QM/MM calculation. No default value. "
        "Character array, max 12 characters. Built-in semiempirical methods (PM3, AM1, etc.) "
        "and DFTB are handled natively. EXTERN uses file-based interface to external QM programs. "
        "QUICK and TERACHEM use API-based interfaces (recommended for performance). "
        "XTB and DFTBPLUS use linked library interfaces."
    ),
    value_type="string",
    section="qmmm",
    category="QM/MM",
    related=["qmmask", "qmcharge", "spin", "ifqnt"],
))



# =============================================================================
# Chapter 4: Generalized Born / Surface Area Model
# =============================================================================
# Keywords from the &cntrl namelist specific to GB implicit solvent.
# Section 4.1: GB/SA input parameters
# Section 4.2: GBION model
# Section 4.3: GTFE (Group Transfer Free Energy) model
# Section 4.4: ALPB (Analytical Linearized Poisson-Boltzmann)

KEYWORDS.append(Keyword(
    name="intdiel",
    description="Sets the interior dielectric constant of the molecule of interest.",
    default=1.0,
    notes="Other values have not been extensively tested.",
    value_type="float",
    section="cntrl",
    category="Generalized Born",
    related=["extdiel", "igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="extdiel",
    description="Sets the exterior or solvent dielectric constant.",
    default=78.5,
    notes="Default of 78.5 corresponds approximately to water.",
    value_type="float",
    section="cntrl",
    category="Generalized Born",
    related=["intdiel", "igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="saltcon",
    description=(
        "Sets the concentration (M) of 1-1 mobile counterions in solution, using a modified "
        "generalized Born theory based on the Debye-Huckel limiting law for ion screening of "
        "interactions. Setting saltcon to a nonzero value does result in some increase in "
        "computation time."
    ),
    default=0.0,
    notes="Unit is Molar (M). Default is 0.0 M (i.e. no Debye-Huckel screening).",
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Generalized Born",
    related=["igb", "intdiel", "extdiel"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="rgbmax",
    description=(
        "Controls the maximum distance between atom pairs that will be considered in carrying out "
        "the pairwise summation involved in calculating the effective Born radii. Atoms whose "
        "associated spheres are farther away than rgbmax from a given atom will not contribute to "
        "that atom's effective Born radius. This is implemented in a 'smooth' fashion, so that when "
        "part of an atom's atomic sphere lies inside rgbmax cutoff, that part contributes to the "
        "low-dielectric region that determines the effective Born radius."
    ),
    default=25.0,
    notes=(
        "The default of 25 Angstroms is usually plenty for single-domain proteins of a few hundred "
        "residues. Even smaller values (10-15 A) are reasonable, changing the functional form of the "
        "generalized Born theory a little bit, in exchange for a considerable speed-up in efficiency. "
        "The rgbmax parameter affects only effective Born radii and their derivatives. The cut parameter "
        "determines the maximum distance for electrostatic, van der Waals and off-diagonal GB terms. "
        "One typically sets rgbmax <= cut."
    ),
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Generalized Born",
    related=["igb", "cut"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="rbornstat",
    description=(
        "If rbornstat = 1, the statistics of the effective Born radii for each atom of the molecule "
        "throughout the molecular dynamics simulation are reported in the output file."
    ),
    default=0,
    options={
        0: "No Born radii statistics reported.",
        1: "Report Born radii statistics in output file.",
    },
    value_type="int",
    section="cntrl",
    category="Generalized Born",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="offset",
    description=(
        "The dielectric radii for generalized Born calculations are decreased by a uniform value "
        "'offset' to give the 'intrinsic radii' used to obtain effective Born radii."
    ),
    default=0.09,
    notes=(
        "Default is 0.09 Angstroms. For igb=8, the default is 0.195141. "
        "Tsui and Case used an offset of 0.13 A with igb=1, which differs from the default."
    ),
    value_type="float",
    section="cntrl",
    category="Generalized Born",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbsa",
    description=(
        "Option to carry out GB/SA (generalized Born/surface area) simulations. Controls whether "
        "and how surface area is computed and included in the solvation term."
    ),
    default=0,
    options={
        0: "Surface area will not be computed and will not be included in the solvation term.",
        1: "Surface area will be computed using the LCPO model.",
        2: (
            "Surface area will be computed by recursively approximating a sphere around an atom, "
            "starting from an icosahedra. No forces are generated; only works for single point "
            "energy calculation, mainly for energy decomposition in MM-GBSA."
        ),
        3: (
            "Surface area will be computed using a fast pairwise approximation suitable for GPU "
            "computing in pmemd.cuda program. About 30 times faster than gbsa=2 on GPU. Not "
            "currently supported in MM-GBSA, QM/MM or libsff. Recommended for use with pmemd.cuda."
        ),
    },
    value_type="int",
    section="cntrl",
    category="Generalized Born",
    related=["igb", "surften"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="surften",
    description=(
        "Surface tension used to calculate the nonpolar contribution to the free energy of solvation "
        "(when gbsa = 1), as Enp = surften*SA."
    ),
    default=0.005,
    notes="Unit is kcal/mol/A^2. For gbsa=3, surften works comparably with gbsa=1 given the same value.",
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Generalized Born",
    related=["gbsa", "igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="rdt",
    description=(
        "This parameter is only used for GB simulations with LES (Locally Enhanced Sampling). "
        "In GB+LES simulations, non-LES atoms require multiple effective Born radii due to "
        "alternate descreening effects of different LES copies. When the multiple radii for a "
        "non-LES atom differ by less than RDT, only a single radius will be used for that atom."
    ),
    default=0.0,
    notes="Unit is Angstroms. See Chapter 32 for more details on LES.",
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Generalized Born",
    related=["igb"],
    commonly_changed=False,
))

# igb=8 element-specific screening parameters (protein)
# These override the defaults when igb=8 is used.
KEYWORDS.append(Keyword(
    name="Sh",
    description="Overlap screening parameter for element H in igb=8 protein GB model.",
    default=1.425952,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="Sc",
    description="Overlap screening parameter for element C in igb=8 protein GB model.",
    default=1.058554,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="Sn",
    description="Overlap screening parameter for element N in igb=8 protein GB model.",
    default=0.733599,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="So",
    description="Overlap screening parameter for element O in igb=8 protein GB model.",
    default=1.061039,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="Ss",
    description="Overlap screening parameter for element S in igb=8 protein GB model.",
    default=-0.703469,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="Sp",
    description="Overlap screening parameter for element P in igb=8 protein GB model.",
    default=0.5,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbneckscale",
    description=(
        "Scaling parameter for the neck correction in GBneck models (igb=7 and igb=8). "
        "The neck correction eliminates interstitial regions of high dielectric smaller "
        "than a solvent molecule."
    ),
    default=0.826836,
    notes="Default of 0.826836 is for igb=8. Recommended for both proteins and nucleic acids.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbalphaH",
    description="The alpha parameter for element H in igb=8 protein GB model.",
    default=0.78844,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbbetaH",
    description="The beta parameter for element H in igb=8 protein GB model.",
    default=0.798699,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbgammaH",
    description="The gamma parameter for element H in igb=8 protein GB model.",
    default=0.437334,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbalphaC",
    description="The alpha parameter for element C in igb=8 protein GB model.",
    default=0.733756,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbbetaC",
    description="The beta parameter for element C in igb=8 protein GB model.",
    default=0.506378,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbgammaC",
    description="The gamma parameter for element C in igb=8 protein GB model.",
    default=0.205844,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbalphaN",
    description="The alpha parameter for element N in igb=8 protein GB model.",
    default=0.503364,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbbetaN",
    description="The beta parameter for element N in igb=8 protein GB model.",
    default=0.316828,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbgammaN",
    description="The gamma parameter for element N in igb=8 protein GB model.",
    default=0.192915,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbalphaOS",
    description="The alpha parameter for elements O and S in igb=8 protein GB model.",
    default=0.867814,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbbetaOS",
    description="The beta parameter for elements O and S in igb=8 protein GB model.",
    default=0.876635,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbgammaOS",
    description="The gamma parameter for elements O and S in igb=8 protein GB model.",
    default=0.387882,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbalphaP",
    description="The alpha parameter for element P in igb=8 protein GB model.",
    default=0.41836,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbbetaP",
    description="The beta parameter for element P in igb=8 protein GB model.",
    default=0.29005,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gbgammaP",
    description="The gamma parameter for element P in igb=8 protein GB model.",
    default=0.10642,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))

# igb=8 nucleic acid parameters (end with 'nu')
KEYWORDS.append(Keyword(
    name="screen_hnu",
    description="Overlap screening parameter for element H in igb=8 nucleic acid GB model.",
    default=1.69654,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="screen_cnu",
    description="Overlap screening parameter for element C in igb=8 nucleic acid GB model.",
    default=1.2689,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="screen_nnu",
    description="Overlap screening parameter for element N in igb=8 nucleic acid GB model.",
    default=1.425974,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="screen_onu",
    description="Overlap screening parameter for element O in igb=8 nucleic acid GB model.",
    default=0.18401,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="screen_pnu",
    description="Overlap screening parameter for element P in igb=8 nucleic acid GB model.",
    default=1.54506,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_alpha_hnu",
    description="The alpha parameter for element H in igb=8 nucleic acid GB model.",
    default=0.53705,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_beta_hnu",
    description="The beta parameter for element H in igb=8 nucleic acid GB model.",
    default=0.36286,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_gamma_hnu",
    description="The gamma parameter for element H in igb=8 nucleic acid GB model.",
    default=0.1167,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_alpha_cnu",
    description="The alpha parameter for element C in igb=8 nucleic acid GB model.",
    default=0.33167,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_beta_cnu",
    description="The beta parameter for element C in igb=8 nucleic acid GB model.",
    default=0.19684,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_gamma_cnu",
    description="The gamma parameter for element C in igb=8 nucleic acid GB model.",
    default=0.09342,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_alpha_nnu",
    description="The alpha parameter for element N in igb=8 nucleic acid GB model.",
    default=0.68631,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_beta_nnu",
    description="The beta parameter for element N in igb=8 nucleic acid GB model.",
    default=0.46319,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_gamma_nnu",
    description="The gamma parameter for element N in igb=8 nucleic acid GB model.",
    default=0.13872,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_alpha_onu",
    description="The alpha parameter for element O in igb=8 nucleic acid GB model.",
    default=0.60634,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_beta_onu",
    description="The beta parameter for element O in igb=8 nucleic acid GB model.",
    default=0.46301,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_gamma_onu",
    description="The gamma parameter for element O in igb=8 nucleic acid GB model.",
    default=0.14226,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_alpha_pnu",
    description="The alpha parameter for element P in igb=8 nucleic acid GB model.",
    default=0.41836,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_beta_pnu",
    description="The beta parameter for element P in igb=8 nucleic acid GB model.",
    default=0.29005,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_gamma_pnu",
    description="The gamma parameter for element P in igb=8 nucleic acid GB model.",
    default=0.10642,
    value_type="float",
    section="cntrl",
    category="Generalized Born - igb8 parameters",
    related=["igb"],
    commonly_changed=False,
))

# Section 4.2.2: GBION input parameters
KEYWORDS.append(Keyword(
    name="gbion",
    description=(
        "Specifies the version of GBION (Implicit Solvent with Explicit Ions) model to be used. "
        "GBION extends the canonical GB model to include explicit ions in simulations, providing "
        "a more nuanced view of ion-solute interactions."
    ),
    default=0,
    options={
        0: "GBION model is turned off (default).",
        1: "Only KGB(a,b) coefficient is available, does not depend on ion type (under development, not recommended).",
        2: "Both KGB(a,b) and K_epsilon coefficients available, but K_epsilon does not depend on ion type (under development, not recommended).",
        3: "Both KGB(a,b) and K_epsilon depend on the respective atom types (currently recommended).",
    },
    notes=(
        "Optimized and tested against explicit solvent distributions of Na+, K+, Cl- and CoHex3+ "
        "ions around DNA. Recommended for simulating DNA duplexes in the presence of these ions. "
        "Parametrized for use with GBneck2 model (igb=8). All KGB and K_epsilon defaults are 1.0, "
        "so GBION is not invoked unless non-default values are specified explicitly."
    ),
    value_type="int",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["igb", "gi_coef_1_p", "gi_coef_1_n", "intdiel_ion_1_p"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gi_coef_1_p",
    description="Sets the KGB(a,b) coefficient for the interaction between solute atoms and cations.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gi_coef_1_n",
    description="Sets the KGB(a,b) coefficient for interactions between solute atoms and anions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gi_coef_2_pp",
    description="Sets the KGB(a,b) coefficient for cation-cation interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gi_coef_2_pn",
    description="Sets the KGB(a,b) coefficient for cation-anion interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gi_coef_2_nn",
    description="Sets the KGB(a,b) coefficient for anion-anion interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="intdiel_ion_1_p",
    description="Sets the internal dielectric constant K_epsilon(a,b) for solute-cation interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="intdiel_ion_1_n",
    description="Sets the internal dielectric constant K_epsilon(a,b) for solute-anion interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="intdiel_ion_2_pp",
    description="Sets the internal dielectric constant K_epsilon(a,b) for cation-cation interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="intdiel_ion_2_pn",
    description="Sets the internal dielectric constant K_epsilon(a,b) for cation-anion interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="intdiel_ion_2_nn",
    description="Sets the internal dielectric constant K_epsilon(a,b) for anion-anion interactions.",
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_neckscale_ion_1_p",
    description=(
        "Sets the KNS(a,b) coefficient for the interaction between solute atoms and cations. These parameters scale the neck correction in GBneck models. "
        "Currently unused (KNS=1 for all atoms and ions) but may change when new ions are parametrized."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion", "gbneckscale"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_neckscale_ion_1_n",
    description=(
        "Sets the KNS(a,b) coefficient for the interaction between solute atoms and anions. These parameters scale the neck correction in GBneck models. "
        "Currently unused (KNS=1 for all atoms and ions) but may change when new ions are parametrized."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion", "gbneckscale"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_neckscale_ion_2_pp",
    description=(
        "Sets the KNS(a,b) coefficient for cation-cation interactions. These parameters scale the neck correction in GBneck models. "
        "Currently unused (KNS=1 for all atoms and ions) but may change when new ions are parametrized."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion", "gbneckscale"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_neckscale_ion_2_pn",
    description=(
        "Sets the KNS(a,b) coefficient for cation-anion interactions. These parameters scale the neck correction in GBneck models. "
        "Currently unused (KNS=1 for all atoms and ions) but may change when new ions are parametrized."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion", "gbneckscale"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gb_neckscale_ion_2_nn",
    description=(
        "Sets the KNS(a,b) coefficient for anion-anion interactions. These parameters scale the neck correction in GBneck models. "
        "Currently unused (KNS=1 for all atoms and ions) but may change when new ions are parametrized."
    ),
    default=1.0,
    value_type="float",
    section="cntrl",
    category="Generalized Born - GBION",
    related=["gbion", "gbneckscale"],
    commonly_changed=False,
))

# Section 4.3: GTFE (Group Transfer Free Energy) model parameters
# These are residue-specific transfer free energies (in kJ/mol) set in &cntrl.
KEYWORDS.append(Keyword(
    name="ala",
    description=(
        "Group transfer free energy for the alanine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="arg",
    description=(
        "Group transfer free energy for the arginine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="asn",
    description=(
        "Group transfer free energy for the asparagine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="asp",
    description=(
        "Group transfer free energy for the aspartate sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="cys",
    description=(
        "Group transfer free energy for the cysteine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gln",
    description=(
        "Group transfer free energy for the glutamine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="glu",
    description=(
        "Group transfer free energy for the glutamate sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="gly",
    description=(
        "Group transfer free energy for the glycine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="his",
    description=(
        "Group transfer free energy for the histidine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="hip",
    description=(
        "Group transfer free energy for the protonated histidine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="ile",
    description=(
        "Group transfer free energy for the isoleucine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="leu",
    description=(
        "Group transfer free energy for the leucine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="lys",
    description=(
        "Group transfer free energy for the lysine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="met",
    description=(
        "Group transfer free energy for the methionine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="phe",
    description=(
        "Group transfer free energy for the phenylalanine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="pro",
    description=(
        "Group transfer free energy for the proline sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="ser",
    description=(
        "Group transfer free energy for the serine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="thr",
    description=(
        "Group transfer free energy for the threonine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="triptophan",
    description=(
        "Group transfer free energy for the tryptophan sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="tyr",
    description=(
        "Group transfer free energy for the tyrosine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="valine",
    description=(
        "Group transfer free energy for the valine sidechain "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="bb",
    description=(
        "Group transfer free energy for the backbone group  "
        "in the GTFE implicit solvent model, in kJ/mol. Used to account for "
        "effects of pressure, temperature, or osmolytes on protein solvation."
    ),
    default=None,
    notes="Set in &cntrl for GTFE simulations. Values are computed using the GTFE python script.",
    value_type="float",
    section="cntrl",
    category="Generalized Born - GTFE",
    related=["igb", "gbsa"],
    commonly_changed=False,
))

# Section 4.4: ALPB (Analytical Linearized Poisson-Boltzmann)
KEYWORDS.append(Keyword(
    name="alpb",
    description=(
        "Flag for using ALPB to handle electrostatic interactions within the implicit solvent model. "
        "The ALPB approximation tends to be more accurate than GB for finite values of solvent dielectric, "
        "with virtually no additional computational overhead."
    ),
    default=0,
    options={
        0: "No ALPB (default).",
        1: (
            "ALPB is turned on. Requires that one of the analytical GB models is also used to compute "
            "the effective Born radii, i.e. one must set igb=1, 2, 5, or 7. The ALPB uses the same "
            "sets of radii as required by the particular GB model."
        ),
    },
    notes=(
        "Statistically significant tests on macromolecular structures have shown that ALPB is more likely "
        "to be a better approximation to PB than the GB. The electrostatic screening effects of monovalent "
        "salt are introduced into the ALPB in the same manner as in the GB, via the saltcon parameter."
    ),
    value_type="int",
    section="cntrl",
    category="Generalized Born - ALPB",
    related=["igb", "arad", "saltcon"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="arad",
    description=(
        "Effective electrostatic size (radius) of the molecule for ALPB. Characterizes the molecule's "
        "overall dimensions and global shape, and is not to be confused with the effective Born radius "
        "of an atom."
    ),
    default=None,
    notes=(
        "An appropriate value should be provided when using ALPB. This is a molecular-level parameter, "
        "not an atomic parameter."
    ),
    value_type="float",
    min_val=0.0,
    section="cntrl",
    category="Generalized Born - ALPB",
    related=["alpb", "igb"],
    commonly_changed=False,
))

# =============================================================================
# Chapter 6: PBSA (Poisson-Boltzmann Surface Area)
# =============================================================================
# The &pb namelist provides detailed control of the numerical PB procedures.
# These keywords are read immediately after the &cntrl namelist.
# Note: Some keywords (imin, ntx, ipb, inp) are in &cntrl and are handled above.
# The 'inp' keyword below controls non-polar solvation in both standalone pbsa
# and in sander PB calculations.

KEYWORDS.append(Keyword(
    name="inp",
    description=(
        "Option to select different methods to compute non-polar solvation free energy "
        "in Poisson-Boltzmann calculations."
    ),
    default=2,
    options={
        0: "No non-polar solvation free energy is computed.",
        1: (
            "The total non-polar solvation free energy is modeled as a single term linearly "
            "proportional to the solvent accessible surface area, as in the PARSE parameter set. "
            "If INP=1, USE_SAV must be equal to 0."
        ),
        2: (
            "The total non-polar solvation free energy is modeled as two terms: the cavity term "
            "and the dispersion term. The dispersion term is computed with a surface-based "
            "integration method. Default."
        ),
    },
    value_type="int",
    section="cntrl",
    category="Poisson-Boltzmann",
    related=["ipb", "use_sav", "cavity_surften", "cavity_offset"],
    commonly_changed=True,
))

# &pb namelist: Physical constants (Section 6.2.3)
KEYWORDS.append(Keyword(
    name="epsin",
    description="Sets the dielectric constant of the solute region. The solute region is defined to be the solvent excluded volume.",
    default=1.0,
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["epsout", "epsmem", "ipb"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="epsout",
    description=(
        "Sets the implicit solvent dielectric constant. The solvent region is defined to be the space "
        "not occupied by the solute region. Only two dielectric regions are allowed in the current release."
    ),
    default=80.0,
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["epsin", "epsmem", "ipb"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="epsmem",
    description=(
        "Sets the membrane dielectric constant. Only used if membraneopt > 0, does nothing otherwise. "
        "Value used should be between epsin and epsout or there may be errors. "
        "Previously spelled as epsmemb, which is being phased out."
    ),
    default=1.0,
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["membraneopt", "epsin", "epsout"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="smoothopt",
    description=(
        "Instructs PB how to set up dielectric values for finite-difference grid edges that are "
        "located across the solute/solvent dielectric boundary."
    ),
    default=1,
    options={
        0: "The dielectric constants of boundary grid edges are always set to the equal-weight harmonic average of EPSIN and EPSOUT.",
        1: "A weighted harmonic average of EPSIN and EPSOUT is used for boundary grid edges. Default.",
        2: "The dielectric constants of boundary grid edges are set to either EPSIN or EPSOUT depending on whether midpoints are inside or outside the solute surface.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["epsin", "epsout", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="istrng",
    description=(
        "Sets the ionic strength (in mM) for the PB equation."
    ),
    default=0.0,
    notes=(
        "Unit is mM, which is different from that (in M) used in the generalized Born methods (saltcon). "
        "Only symmetrical solutions are supported, so the ionic strength equals the square of the valence "
        "of the symmetrical ions times the ion concentration (in mM)."
    ),
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["ipb", "pbtemp", "ivalence"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="pbtemp",
    description="Temperature (in K) used for the PB equation, needed to compute the Boltzmann factor for salt effects.",
    default=300.0,
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["istrng", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="ivalence",
    description="Valence of the symmetrical ions for Debye-Huckel screening in PB calculations.",
    default=1,
    value_type="int",
    min_val=1,
    section="pb",
    category="Poisson-Boltzmann - Physical constants",
    related=["istrng", "ipb"],
    commonly_changed=False,
))

# &pb namelist: Surface and radii options
KEYWORDS.append(Keyword(
    name="radiopt",
    description="Option to set up atomic radii for PB calculations.",
    default=1,
    options={
        0: "Use radii from the prmtop file for both PB and NP calculations.",
        1: (
            "Use atom-type/charge-based radii by Tan and Luo for the PB calculation. "
            "Note that the radii are optimized for Amber atom types as in standard residues. "
            "If GAFF atom types are used, radii from the prmtop file will be used. Default."
        ),
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["ipb", "dprob"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="dprob",
    description="Solvent probe radius for molecular surface used to define the dielectric boundary between solute and solvent.",
    default=1.4,
    notes="Unit is Angstroms.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["iprob", "sasopt", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="iprob",
    description="Mobile ion probe radius for ion accessible surface used to define the Stern layer.",
    default=2.0,
    notes="Unit is Angstroms.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["dprob", "istrng", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="sasopt",
    description="Option to determine which kind of molecular surfaces to be used in the PB implicit solvent model.",
    default=0,
    options={
        0: "Use the solvent excluded surface (SES). Default.",
        1: "Use the solvent accessible surface (SAS). Reduces to VDW surface when dprob is set to zero.",
        2: "Use the smooth surface defined by a revised density function. Must be combined with IPB >= 2.",
        3: "Use the solvent excluded surface inferred by the machine-learned MLSES model. Must be combined with IPB=2.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["ipb", "dprob", "mlses_opt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="saopt",
    description=(
        "Option to compute the surface area of a molecule. Once enabled, the surface area will be "
        "reported in the output file with subtitle 'Total molecular surface'. Only SES and SAS "
        "surface areas are supported."
    ),
    default=0,
    options={
        0: "Do not compute any surface area.",
        1: "Use the field-view method to compute the surface area.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["sasopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="triopt",
    description=(
        "Option to add trimer arc dots for a more accurate and lower memory mapping method "
        "of the analytical solvent excluded surface."
    ),
    default=1,
    options={
        0: "Trimer arc dots are not used.",
        1: "Trimer arc dots are used. Default.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["arcres", "sasopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="arcres",
    description=(
        "Resolution (in Angstroms) of dots used to represent solvent accessible arcs in the numerical "
        "surface computation. Should be reduced to 0.125 A when TRIOPT is turned off to achieve "
        "similar accuracy in reaction field energies."
    ),
    default=0.25,
    notes=(
        "More generally, ARCRES should be set to max(0.125 A, 0.5*h) when TRIOPT is on, "
        "or max(0.0625 A, 0.25*h) when TRIOPT is off, where h is the grid spacing."
    ),
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Surface",
    related=["triopt", "space"],
    commonly_changed=False,
))

# &pb namelist: Implicit membrane options (Section 6.2.4)
KEYWORDS.append(Keyword(
    name="membraneopt",
    description=(
        "Option to turn the implicit membrane on and off. The membrane is implemented as a slab-like "
        "region with a uniform or heterogeneous dielectric constant depth profile."
    ),
    default=0,
    options={
        0: "No implicit membrane used (default).",
        1: "Use a uniform membrane dielectric constant in a slab-like implicit membrane.",
        2: (
            "Use a heterogeneous membrane dielectric constant varying with depth from 1 in the membrane "
            "center to 80 at the periphery, using PCHIP fitting."
        ),
        3: (
            "Use a heterogeneous membrane dielectric constant varying with depth from 1 in the membrane "
            "center to 80 at the periphery, using Spline fitting."
        ),
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Membrane",
    related=["epsmem", "mthick", "mctrdz", "mprob", "poretype"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="mprob",
    description="Membrane probe radius in Angstroms. Used to specify the lipid molecule accessibility versus water.",
    default=2.70,
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Membrane",
    related=["membraneopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="mthick",
    description="Membrane thickness in Angstroms.",
    default=40.0,
    notes="This is different from the previous default of 20 A.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Membrane",
    related=["membraneopt", "mctrdz"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="mctrdz",
    description="Membrane center in Angstroms in the z direction.",
    default=0.0,
    notes="Default of 0 means membrane centered at the center of the protein.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Membrane",
    related=["membraneopt", "mthick"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="poretype",
    description="Turn on and off the automatic depth-first search method to identify pores in membrane proteins.",
    default=0,
    options={
        0: "Do not turn on the pore searching algorithm.",
        1: "Turn on the pore searching algorithm.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Membrane",
    related=["membraneopt", "poreradius"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="poreradius",
    description=(
        "Controls the radius, in Angstroms, of the cylindrical exclusion region for membrane pore. "
        "This is no longer needed given the automatic pore searching algorithm."
    ),
    default=None,
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Membrane",
    related=["poretype", "membraneopt"],
    commonly_changed=False,
))

# &pb namelist: Numerical procedures (Section 6.2.5)
KEYWORDS.append(Keyword(
    name="npbopt",
    description="Option to select the linear or the full nonlinear PB equation.",
    default=0,
    options={
        0: "Linear PB equation is solved. Default.",
        1: "Nonlinear PB equation is solved.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["solvopt", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="solvopt",
    description="Option to select iterative solvers for the PB equation.",
    default=1,
    options={
        1: "Modified ICCG or Periodic (PICCG) if bcopt=10. Default.",
        2: "Geometric multigrid. A four-level v-cycle implementation is applied by default.",
        3: "Conjugate gradient (Periodic version available under bcopt=10). Requires large MAXITN.",
        4: "SOR. Requires large MAXITN to converge.",
        5: "Adaptive SOR. Only compatible with NPBOPT=1. Requires large MAXITN.",
        6: "Damped SOR. Only compatible with NPBOPT=1. Requires large MAXITN.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["npbopt", "maxitn", "accept", "bcopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="accept",
    description="Sets the iteration convergence criterion (relative to the initial residue) for FD solvers.",
    default=0.001,
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["maxitn", "solvopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="maxitn",
    description=(
        "Sets the maximum number of iterations for the finite difference solvers."
    ),
    default=100,
    notes="MAXITN has to be set to a much larger value, e.g. 10000, for less efficient solvers such as CG and SOR.",
    value_type="int",
    min_val=1,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["accept", "solvopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="fillratio",
    description=(
        "The ratio between the longest dimension of the rectangular finite-difference grid and that "
        "of the solute."
    ),
    default=2.0,
    notes=(
        "A larger FILLRATIO (e.g. 4.0) should be used for small solutes such as ligand molecules. "
        "Otherwise, part of the small solute may lie outside the finite-difference grid."
    ),
    value_type="float",
    min_val=1.0,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["space", "nfocus", "nbuffer"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="space",
    description="Sets the grid spacing for the finite difference solver, in Angstroms.",
    default=0.5,
    notes="For PB dynamics/minimization in sander, a finer spacing of 0.25 A is recommended.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["fillratio", "arcres"],
    commonly_changed=True,
))
KEYWORDS.append(Keyword(
    name="nbuffer",
    description=(
        "Sets how far away (in grid units) the boundary of the finite difference grid is from the "
        "solute surface."
    ),
    default=0,
    notes="Default of 0 means automatically set to be at least a solvent probe or ion probe (diameter) away from the solute surface.",
    value_type="int",
    min_val=0,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["fillratio", "space"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="nfocus",
    description=(
        "Set how many successive FD calculations will be used to perform an electrostatic "
        "focussing calculation on a molecule."
    ),
    default=2,
    notes="Maximum is 2. When NFOCUS=1, no focusing is used. Recommended NFOCUS=1 when multigrid solver is used.",
    value_type="int",
    min_val=1,
    max_val=2,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["fscale", "solvopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="fscale",
    description="Set the ratio between the coarse and fine grid spacings in an electrostatic focussing calculation.",
    default=8,
    notes="For PB dynamics in sander, a value of 4 is recommended.",
    value_type="int",
    min_val=1,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["nfocus", "space"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="npbgrid",
    description=(
        "Sets how often the finite-difference grid is regenerated during dynamics."
    ),
    default=1,
    notes=(
        "For molecular dynamics simulations, recommended to be set to at least 100. "
        "The PB solver takes advantage of slowly varying electrostatic potential distributions "
        "by keeping the grid fixed, but molecules move freely so the grid must be regenerated occasionally."
    ),
    value_type="int",
    min_val=1,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["nsnba"],
    commonly_changed=False,
))

# &pb namelist: Energy and force options (Section 6.2.6)
KEYWORDS.append(Keyword(
    name="bcopt",
    description="Boundary condition options for the PB finite-difference grid.",
    default=5,
    options={
        1: "Boundary grid potentials set to zero (conductor). Total electrostatic potentials and energy computed.",
        5: "Boundary grid potentials computed using all grid charges. Total electrostatic potentials and energy computed. Default.",
        6: "Boundary grid potentials computed using all grid charges. Reaction field potentials and energy computed with charge singularity free formalism.",
        10: (
            "Periodic boundary condition. Total electrostatic potentials and energy computed. "
            "Compatible with SOLVOPT=1,2,3,4 and IPB=1 or 2. Should only be used on charge-neutral systems."
        ),
    },
    notes="For PB dynamics in sander, BCOPT=6 is recommended for stability.",
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["eneopt", "solvopt", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="eneopt",
    description="Option to compute total electrostatic energy and forces.",
    default=2,
    options={
        1: (
            "Compute total electrostatic energy and forces with the P3M procedure. "
            "EPB is set to zero, EEL includes both reaction field and Coulombic energy. "
            "Requires nonzero CUTNB and BCOPT=5."
        ),
        2: (
            "Use dielectric boundary surface charges to compute reaction field energy. Default. "
            "Both Coulombic and van der Waals energies computed via pairwise atomic interactions. "
            "EPB gives reaction field energy, EEL gives Coulombic energy."
        ),
        3: "P3M procedure for both solvation and Coulombic energy/forces for larger systems.",
        4: "P3M procedure for the full nonlinear PB equation for both solvation and Coulombic energy/forces.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["frcopt", "bcopt", "cutnb", "cutfd"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="frcopt",
    description="Option to compute and output electrostatic forces to a file named force.dat.",
    default=0,
    options={
        0: "Do not compute or output atomic and total electrostatic forces. Default.",
        1: "Reaction field forces by trilinear interpolation. Dielectric boundary forces using electric field on boundary. Output in kcal/mol*A.",
        2: "Use dielectric boundary surface polarized charges to compute reaction field and dielectric boundary forces. Output in kcal/mol*A.",
        3: "Reaction field forces using dielectric boundary polarized charge. Dielectric boundary forces using electric field on boundary. Output in kcal/mol*A.",
    },
    notes="For PB dynamics in sander, FRCOPT=2 is recommended.",
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["eneopt", "bcopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="scalec",
    description="Option to scale dielectric boundary surface charges before computing reaction field energy and forces.",
    default=0,
    options={
        0: "Do not scale. Default.",
        1: "Scale using Gauss's law.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["eneopt", "frcopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="cutfd",
    description=(
        "Atom-based cutoff distance to remove short-range finite-difference interactions, "
        "and to add pairwise charge-based interactions. Used for both energy and force calculations."
    ),
    default=5.0,
    notes="Unit is Angstroms.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["cutnb", "eneopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="cutnb",
    description=(
        "Atom-based cutoff distance for van der Waals interactions, and pairwise Coulombic "
        "interactions when ENEOPT=2. When set to default 0, no cutoff is used (all pairwise "
        "interactions included). When ENEOPT=1, this is the cutoff for van der Waals only."
    ),
    default=0.0,
    notes="Unit is Angstroms. Default of 0 means infinite cutoff.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["cutfd", "eneopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="cutsa",
    description=(
        "Cutoff distance for SASA-related interactions in PB calculations."
    ),
    default=None,
    notes="Appears in example inputs. Used with eneopt and cutnb.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Energy/Forces",
    related=["cutnb", "cutfd"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="nsnba",
    description="Sets how often atom-based pairlist is generated for PB calculations.",
    default=1,
    notes="For molecular dynamics simulations, a value of 5 is recommended.",
    value_type="int",
    min_val=1,
    section="pb",
    category="Poisson-Boltzmann - Numerical",
    related=["npbgrid"],
    commonly_changed=False,
))

# &pb namelist: Visualization and output options (Section 6.2.7)
KEYWORDS.append(Keyword(
    name="phiout",
    description="Option to output spatial distribution of electrostatic potential for visualization.",
    default=0,
    options={
        0: "No potential file is printed out. Default.",
        1: "Electrostatic potential is printed out in a file named pbsa.phi in the working directory.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Output",
    related=["phiform"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="phiform",
    description="Controls the format of the electrostatic potential file.",
    default=0,
    options={
        0: "Delphi binary format (kT/mol*e). Default.",
        1: "Amber ASCII format (kcal/mol*e).",
        2: "DX volumetric data format (kcal/mol*e) for use with VMD.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Output",
    related=["phiout"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="outlvlset",
    description=(
        "Option to write the total level set (both solute-solvent and membrane level sets combined) "
        "to a DX format volumetric data file named pbsa_lvlset.dx."
    ),
    default="false",
    options={
        "false": "No level set file printed out. Default.",
        "true": "Level set printed out in pbsa_lvlset.dx.",
    },
    value_type="string",
    section="pb",
    category="Poisson-Boltzmann - Output",
    related=["outmlvlset", "membraneopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="outmlvlset",
    description=(
        "Option to write the membrane level set to a separate DX format volumetric data file. "
        "Does nothing if membraneopt is not turned on."
    ),
    default="false",
    options={
        "false": "No membrane level set file printed out. Default.",
        "true": "Membrane level set printed out.",
    },
    value_type="string",
    section="pb",
    category="Poisson-Boltzmann - Output",
    related=["outlvlset", "membraneopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="npbverb",
    description="When set to 1, turns on verbose mode in pbsa; generates detailed information for PB and NP calculations.",
    default=0,
    options={
        0: "Normal output. Default.",
        1: "Verbose mode with detailed PB/NP information.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Output",
    related=["ipb"],
    commonly_changed=False,
))

# &pb namelist: Non-polar solvation options (Section 6.2.8)
KEYWORDS.append(Keyword(
    name="decompopt",
    description=(
        "Option to select different decomposition schemes when INP=2 for non-polar solvation."
    ),
    default=2,
    options={
        1: "The 6/12 decomposition scheme.",
        2: "The sigma decomposition scheme. Default. Best of the three schemes studied.",
        3: "The WCA decomposition scheme.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["inp", "use_rmin", "sprob"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="use_rmin",
    description="Option to set up van der Waals radii for non-polar calculations.",
    default=1,
    options={
        0: "Use atomic van der Waals sigma values.",
        1: "Use atomic van der Waals rmin values. Default. Improves agreement with TIP3P.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["decompopt", "inp"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="sprob",
    description=(
        "Solvent probe radius for solvent accessible surface area (SASA) used to compute "
        "the dispersion term in non-polar solvation calculations."
    ),
    default=0.557,
    notes=(
        "Default of 0.557 A is for the sigma decomposition scheme, optimized with respect to TIP3P "
        "solvent and PME treatment. Recommended values for other schemes are in Table 4 of Ref. [243]. "
        "If USE_SAV=0, SPROB can also be used for the cavity term but with a different recommended value."
    ),
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["vprob", "decompopt", "inp"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="vprob",
    description=(
        "Solvent probe radius for molecular volume (the volume enclosed by SASA) used to compute "
        "non-polar cavity solvation free energy."
    ),
    default=1.300,
    notes="Default of 1.300 A optimized with respect to TIP3P solvent.",
    value_type="float",
    min_val=0.0,
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["sprob", "use_sav", "decompopt"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="rhow_effect",
    description="Effective water density used in the non-polar dispersion term calculation.",
    default=1.129,
    notes="Default of 1.129 is for the sigma decomposition scheme, optimized with respect to TIP3P solvent in PME.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["decompopt", "inp"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="use_sav",
    description=(
        "Option to use molecular volume (enclosed by SASA) or molecular surface (SASA) for cavity "
        "term calculation. Molecular volume transfers better from small molecules to biomacromolecules."
    ),
    default=1,
    options={
        0: "Use SASA to estimate cavity free energy.",
        1: "Use the molecular volume enclosed by SASA. Default.",
    },
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["cavity_surften", "cavity_offset", "vprob", "inp"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="cavity_surften",
    description=(
        "The regression coefficient for the linear relation between the total non-polar solvation "
        "free energy (INP=1) or the cavity free energy (INP=2) and SASA/volume enclosed by SASA."
    ),
    default=None,
    notes="Default value depends on INP and is set to the best tested scheme (DECOMPOPT=2, USE_RMIN=1, USE_SAV=1).",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["cavity_offset", "use_sav", "inp"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="cavity_offset",
    description=(
        "The regression offset for the linear relation between the total non-polar solvation "
        "free energy (INP=1) or the cavity free energy (INP=2) and SASA/volume enclosed by SASA."
    ),
    default=None,
    notes="Default value depends on INP and is set to the best tested scheme.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["cavity_surften", "use_sav", "inp"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="maxsph",
    description=(
        "Approximate number of dots to represent the maximum atomic solvent accessible surface "
        "for SASA computation. Dots are checked against bonded and non-bonded atoms to identify "
        "buried vs. exposed portions."
    ),
    default=400,
    value_type="int",
    min_val=1,
    section="pb",
    category="Poisson-Boltzmann - Non-polar",
    related=["sasopt"],
    commonly_changed=False,
))

# &pb namelist: MLSES options (Section 6.2.9)
KEYWORDS.append(Keyword(
    name="mlses_opt",
    description=(
        "Option to select the runtime for the Machine-Learned Solvent Excluded Surface (MLSES) model. "
        "Requires SASOPT=3 to be set."
    ),
    default=0,
    options={
        0: "Use customized Fortran function/CUDA kernel to run GENIUSES model on CPU/GPU. Default.",
        1: "Use Torch PBSA runtime to run GENIUSES model on CPU or GPU.",
        2: "Use Torch PBSA runtime to run Con2SES-2D model on CPU or GPU.",
        3: "Use Torch PBSA runtime to run Con2SES-3D model on CPU or GPU. Recommended over Con2SES-2D.",
    },
    notes="Options 1-3 require installation of LibTorch.",
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - MLSES",
    related=["sasopt", "ipb"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="mlses_bench",
    description="Benchmarking option for MLSES surface generation.",
    default=0,
    value_type="int",
    section="pb",
    category="Poisson-Boltzmann - MLSES",
    related=["mlses_opt", "sasopt"],
    commonly_changed=False,
))

# &pb namelist: Active site focusing options (Section 6.2.10)
KEYWORDS.append(Keyword(
    name="xmax",
    description="The upper boundary of the active site focusing box in the x direction.",
    default=0.0,
    notes="When all six box boundaries are zero (default), the original electrostatic focusing is used.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Active site focusing",
    related=["nfocus"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="xmin",
    description="The lower boundary of the active site focusing box in the x direction.",
    default=0.0,
    notes="When all six box boundaries are zero (default), the original electrostatic focusing is used.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Active site focusing",
    related=["nfocus"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="ymax",
    description="The upper boundary of the active site focusing box in the y direction.",
    default=0.0,
    notes="When all six box boundaries are zero (default), the original electrostatic focusing is used.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Active site focusing",
    related=["nfocus"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="ymin",
    description="The lower boundary of the active site focusing box in the y direction.",
    default=0.0,
    notes="When all six box boundaries are zero (default), the original electrostatic focusing is used.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Active site focusing",
    related=["nfocus"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="zmax",
    description="The upper boundary of the active site focusing box in the z direction.",
    default=0.0,
    notes="When all six box boundaries are zero (default), the original electrostatic focusing is used.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Active site focusing",
    related=["nfocus"],
    commonly_changed=False,
))
KEYWORDS.append(Keyword(
    name="zmin",
    description="The lower boundary of the active site focusing box in the z direction.",
    default=0.0,
    notes="When all six box boundaries are zero (default), the original electrostatic focusing is used.",
    value_type="float",
    section="pb",
    category="Poisson-Boltzmann - Active site focusing",
    related=["nfocus"],
    commonly_changed=False,
))



# =============================================================================
# Chapter 11: QM/MM Calculations
# =============================================================================
# Updates to existing &qmmm keywords with full option details from Chapter 11.
# New &qmmm keywords for link atoms, switching, QXD, abfQM/MM.
# External QM program namelists: &adf, &gms, &gau, &orc, &qc, &mrcc, &fb
# GPU QM engines: &quick, &tc
# Tight-binding interfaces: &xtb, &dftbplus
# ML corrections: &dprc
# Full QM dynamics: &sebomd
# Adaptive solvent: &vsolv, &adqmmm

# --- New &qmmm namelist keywords from Chapter 11 ---

KEYWORDS.append(Keyword(
    name="writepdb",
    description="Write a PDB file of the QM region on the first step to qmmm_region.pdb for verification.",
    default=0,
    options={
        0: "Do not write a PDB file of the QM region (default).",
        1: "Write a crude PDB file of the atoms in the QM region on the first step.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM Output",
    related=["qmmask", "iqmatoms"],
    commonly_changed=False,
))

# --- QXD charge-dependent exchange-dispersion parameters (Section 11.1.8) ---

KEYWORDS.append(Keyword(
    name="qxd",
    description="Invoke the charge-dependent exchange-dispersion correction for QM/MM van der Waals interactions. Replaces standard Lennard-Jones QM-MM interactions with a charge-dependent model.",
    default=".false.",
    options={".false.": "Do not use charge-dependent vdW model (default).", ".true.": "Use charge-dependent exchange-dispersion correction for QM/MM vdW."},
    notes="Only affects QM/MM interactions in sander, not pure QM calculations through sqm. Default parameters reproduce regular LJ interactions for typical atom types (HC, C*, N, OW). Parameters for F and Cl are also available.",
    value_type="string",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd_s", "qxd_z0", "qxd_zq", "qxd_d0", "qxd_dq", "qxd_q0", "qxd_qq", "qxd_neff"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_s",
    description="QXD overlap scale parameter (s). Controls the scaling of overlap in the charge-dependent vdW model.",
    default=None,
    notes="Default values reproduce regular LJ interactions for typical atom types. Can be modified via external parameter files.",
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_z0", "qxd_neff"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_z0",
    description="QXD zeta-zero parameter. Baseline exponent for the exchange-dispersion model.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_zq"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_zq",
    description="QXD zeta-q parameter. Charge dependence of the exponent for exchange-dispersion.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_z0"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_d0",
    description="QXD dispersion d0 parameter. Baseline C6 dispersion coefficient.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_dq"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_dq",
    description="QXD dispersion dq parameter. Charge dependence of C6 dispersion.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_d0"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_q0",
    description="QXD q0 parameter. Second dispersion baseline coefficient.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_qq"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_qq",
    description="QXD qq parameter. Charge dependence of the second dispersion coefficient.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_q0"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qxd_neff",
    description="QXD effective electron number Neff(0). Controls the effective number of electrons in the model.",
    default=None,
    value_type="float",
    section="qmmm",
    category="QM/MM QXD",
    related=["qxd", "qxd_s"],
    commonly_changed=False,
))

# --- abfQM/MM radii and extension parameters (Section 11.9.5) ---

KEYWORDS.append(Keyword(
    name="r_core_in",
    description="Inner radius (Angstroms) for determining core extension region around user-specified core atoms in abfQM/MM.",
    default=0.0,
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["r_core_out", "abfqmmm", "coremask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="r_core_out",
    description="Outer radius (Angstroms) for determining core extension region. Defaults to r_core_in if not set.",
    default=None,
    notes="If r_core_out < r_core_in then r_core_out = r_core_in.",
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["r_core_in", "abfqmmm", "coremask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="r_qm_in",
    description="Inner radius (Angstroms) for determining qm extension region around core and user-specified qm atoms in abfQM/MM.",
    default=0.0,
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["r_qm_out", "abfqmmm", "qmmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="r_qm_out",
    description="Outer radius (Angstroms) for determining qm extension region. Defaults to r_qm_in if not set.",
    default=None,
    notes="If r_qm_out < r_qm_in then r_qm_out = r_qm_in.",
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["r_qm_in", "abfqmmm", "qmmask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="r_buffer_in",
    description="Inner radius (Angstroms) for determining buffer extension region around qm and core atoms in abfQM/MM.",
    default=0.0,
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["r_buffer_out", "abfqmmm", "buffermask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="r_buffer_out",
    description="Outer radius (Angstroms) for determining buffer extension region. Defaults to r_buffer_in if not set.",
    default=None,
    notes="If r_buffer_out < r_buffer_in then r_buffer_out = r_buffer_in.",
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["r_buffer_in", "abfqmmm", "buffermask"],
    commonly_changed=False,
))

# --- abfQM/MM miscellaneous parameters (Section 11.9.5.2-3) ---

KEYWORDS.append(Keyword(
    name="nchain",
    description="Number of thermostats in each Nose-Hoover chain for adaptive thermostats (ntt=5,7,8) used in abfQM/MM.",
    default=1,
    value_type="int",
    min_val=1,
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "ntt"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="selection_type",
    description="Type of region selection for abfQM/MM adaptive regions.",
    default=1,
    options={
        1: "Atom-atom distance based selection (default). Atom belongs to outer region if distance to any atom in inner region <= criterion.",
        2: "Flexible sphere selection. Region radius calculated from largest distance between centre of mass and any atom in that region.",
        3: "Fixed sphere based selection. Same as 2 but only innermost region edge is calculated; others are concentric spheres.",
    },
    notes="Options 2 and 3 select significantly more atoms than option 1.",
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "centermask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="initial_selection_type",
    description="Controls the initial region selection when not restarting from an abfQM/MM restart file.",
    default=0,
    options={
        -1: "Use the inner radius for the first selection.",
        0: "Middle sphere selection using the mean of inner and outer radii (default).",
        1: "Use the outer radius for the first selection.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "read_idrst_file"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="center_type",
    description="Type of center calculation for selection_type=2 and 3 in abfQM/MM.",
    default=1,
    options={
        1: "Center of mass (default).",
        2: "Geometric center.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "selection_type", "centermask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="gamma_ln_qm",
    description="Collision frequency (ps^-1) of core and qm atoms when adaptive massive Langevin thermostat is used in abfQM/MM.",
    default=None,
    notes="Defaults to the gamma_ln value defined in &cntrl.",
    value_type="float",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "gamma_ln"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mom_cons_type",
    description="Type of force correction for momentum conservation in abfQM/MM.",
    default=1,
    options={
        0: "No momentum conservation applied.",
        1: "Extra force distributed as equal accelerations among selected atoms (default).",
        2: "Equal forces applied on each atom.",
        -1: "Acceleration proportional to absolute value of current acceleration.",
        -2: "Force proportional to absolute value of current force.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "mom_cons_region"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mom_cons_region",
    description="Specifies the region where force correction for momentum conservation is distributed in abfQM/MM.",
    default=1,
    options={
        0: "Distribution applied only among core atoms.",
        1: "Distribution among current core+qm atoms (default).",
        2: "Distribution among current core+qm+buffer atoms.",
        3: "Distribution applied to all atoms.",
    },
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "mom_cons_type"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="fix_atom_list",
    description="Activates fixed atom list method for abfQM/MM. When >0, regions are extended only by solvent molecules satisfying geometrical criteria.",
    default=0,
    notes="Useful when only solvent exchange is expected and no solute atoms should be selected beyond user-specified ones.",
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "solvent_atom_number"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="solvent_atom_number",
    description="Number of atoms in solvent molecule for fixed atom list mode (fix_atom_list > 0) in abfQM/MM.",
    default=3,
    notes="Default of 3 is appropriate for water. Adjust for other solvents.",
    value_type="int",
    min_val=1,
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "fix_atom_list"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="oxidation_number_list_file",
    description="Filename containing oxidation numbers for abfQM/MM. Each line: RES ATOM OXID with hierarchical assignment.",
    default=None,
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ext_coremask_subset",
    description="Possible core extension atom set for abfQM/MM. Only atoms in this set are considered for core extension based on geometrical criteria.",
    default=None,
    notes="If not defined, defaults to the all atom list.",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "coremask", "r_core_in"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ext_qmmask_subset",
    description="Possible qm extension atom set for abfQM/MM. Only atoms in this set are considered for qm extension based on geometrical criteria.",
    default=None,
    notes="If not defined, defaults to the all atom list.",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "qmmask", "r_qm_in"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ext_buffermask_subset",
    description="Possible buffer extension atom set for abfQM/MM. Only atoms in this set are considered for buffer extension based on geometrical criteria.",
    default=None,
    notes="If not defined, defaults to the all atom list.",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "buffermask", "r_buffer_in"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cut_bond_list_file",
    description="Filename of breakable bonds for intelligent termination of regions in abfQM/MM. Format: ATOM1 ARROW ATOM2 per line.",
    default=None,
    notes="ARROW can be '=>', '<=' (one direction) or '<=>' (either direction). ATOM1/ATOM2 are atom types or indices.",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "max_bonds_per_atom"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="max_bonds_per_atom",
    description="Maximum number of ligands around any atom for the intelligent termination scheme in abfQM/MM.",
    default=4,
    notes="Default of 4 is sufficient for most biological systems. Increase if atoms have more than 4 ligands.",
    value_type="int",
    min_val=1,
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "cut_bond_list_file"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="n_max_recursive",
    description="Maximum number of iterations for the recursive intelligent termination subroutine in abfQM/MM.",
    default=10000,
    value_type="int",
    min_val=1,
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "cut_bond_list_file"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="pdb_file",
    description="Filename for the special abfQM/MM PDB file containing region assignments and atom id information.",
    default="abfqmmm.pdb",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "ntwpdb"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwpdb",
    description="Frequency of writing abfQM/MM information to the PDB file. 0 = no printing. Negative values perform a selection test only.",
    default=0,
    notes="With ntwpdb < 0, no dynamics or point calculations are performed; the program terminates after printing the PDB file.",
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "pdb_file"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="read_idrst_file",
    description="Filename of abfQM/MM atom id restart file for restarting simulations with the same region specifications.",
    default=None,
    notes="Avoids natural transient period in region equilibration during restart. Requires same region specifications as previous run.",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "write_idrst_file"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="write_idrst_file",
    description="Filename of abfQM/MM atom id restart file generated during the run.",
    default="abfqmmm.idrst",
    value_type="string",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "read_idrst_file", "ntwidrst"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwidrst",
    description="Frequency of writing the abfQM/MM atom id restart file.",
    default=0,
    value_type="int",
    section="qmmm",
    category="QM/MM Adaptive",
    related=["abfqmmm", "write_idrst_file"],
    commonly_changed=False,
))


# =============================================================================
# &adf NAMELIST (Section 11.2.6.1 - AMBER/ADF interface)
# =============================================================================
# Amsterdam Density Functional interface. Only mechanical embedding supported.

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set type for ADF DFT calculation. Valid types: SZ, DZ, DZP, TZP, TZ2P, TZ2P+, ZORA/QZ4P.",
    default="DZP",
    value_type="string",
    section="adf",
    category="QM/MM External - ADF",
    related=["xc"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="core",
    description="Type of frozen core for ADF. Allowed: None, Small, Medium, Large.",
    default="None",
    value_type="string",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="zlmfit",
    description="Quality of density fit with ZLM fit method for ADF.",
    default="good",
    value_type="string",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="fit_type",
    description="Fit basis set for density fitting with old pair fit method. Empty string uses ZLM fit by default.",
    default="",
    value_type="string",
    section="adf",
    category="QM/MM External - ADF",
    related=["zlmfit"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="xc",
    description="Exchange-correlation functional for ADF. Examples: LDA VWN, GGA BLYP, GGA PBE, HYBRID B3LYP, HYBRID PBE0.",
    default="GGA BLYP",
    value_type="string",
    section="adf",
    category="QM/MM External - ADF",
    related=["basis"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scf_iter",
    description="Maximum number of SCF cycles allowed for ADF.",
    default=50,
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scf_conv",
    description="SCF convergence threshold for ADF. Convergence achieved when max element of Fock/density commutator < scf_conv.",
    default=1e-06,
    value_type="float",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="beckegrid",
    description="Quality of Becke integration grid for ADF. Allowed: Normal, Good, VeryGood.",
    default="Good",
    value_type="string",
    section="adf",
    category="QM/MM External - ADF",
    related=["integration"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="integration",
    description="Numerical integration accuracy for old te Velde-Baerends grid. Set > 0 (recommend >= 5.0) to use old grid. -1 uses Becke grid.",
    default=-1.0,
    value_type="float",
    section="adf",
    category="QM/MM External - ADF",
    related=["beckegrid"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="num_threads",
    description="Number of threads (CPU cores) for ADF. 0 = use all available cores.",
    default=0,
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_dftb",
    description="Use DFTB with ADF's dftb program instead of DFT. Only works with older ADF DFTB versions (prior to 2011).",
    default=0,
    options={
        0: "Regular DFT calculation (default).",
        1: "Use DFTB; only charge and scf_conv variables considered.",
    },
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="exactdensity",
    description="Use exact (not fitted) electron density for XC potential evaluation in ADF.",
    default=0,
    options={
        0: "Use fitted density (default).",
        1: "Use exact density.",
    },
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a user-provided template file (adf_job.tpl) for ADF input.",
    default=0,
    options={
        0: "Do not use template (default).",
        1: "Use template file.",
    },
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole moment to adf_job.dip. Defaults to &cntrl ntpr.",
    default=None,
    notes="Defaults to &cntrl namelist variable ntpr.",
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipole",
    description="Toggle writing of dipole moment to adf_job.dip.",
    default=0,
    options={
        0: "Do not write dipole (default).",
        1: "Write dipole moment.",
    },
    value_type="int",
    section="adf",
    category="QM/MM External - ADF",
    related=["printdipole"],
    commonly_changed=False,
))


# =============================================================================
# &gms NAMELIST (Section 11.2.6.2 - AMBER/GAMESS-US interface)
# =============================================================================
# GAMESS-US interface. Only mechanical embedding supported.
# Available QM models limited to HF, DFT and MP2 (analytical gradients).

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for GAMESS. Supports Pople (STO-3G, 6-31G, etc.), Karlsruhe (KTZV, KTZVP, KTZVPP), and Dunning (CCn, ACCn) types.",
    default="6-31G*",
    value_type="string",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="method",
    description="QM method for GAMESS. HF, MP2, or a DFT functional (BP86, BLYP, PBE, B3LYP, PBE0).",
    default="BP86",
    value_type="string",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dfttyp",
    description="DFT functional type for GAMESS when method='DFT'. Examples: PBE, BLYP, B3LYP.",
    default=None,
    value_type="string",
    section="gms",
    category="QM/MM External - GAMESS",
    related=["method"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nrad",
    description="Number of radial points in Euler-MacLaurin quadrature for XC integration in GAMESS.",
    default=96,
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=["nleb"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nleb",
    description="Number of angular points in Lebedev grids for XC quadrature in GAMESS.",
    default=590,
    notes="GAMESS default of 302 is not accurate enough to conserve energy.",
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scf_conv",
    description="SCF convergence threshold for GAMESS. Convergence when absolute density change < scf_conv.",
    default=1e-06,
    value_type="float",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxit",
    description="Maximum number of SCF iterations for GAMESS.",
    default=50,
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="gms_version",
    description="Version number specified when building GAMESS. Used to find the correct executable.",
    default="00",
    value_type="string",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="num_threads",
    description="Number of threads (CPU cores) for GAMESS. May require special rungms setup.",
    default=1,
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mwords",
    description="Maximum replicated memory per node in units of 1,000,000 64-bit words for GAMESS.",
    default=50,
    notes="Increase if GAMESS crashes due to insufficient memory.",
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a user-provided template file (gms_job.tpl) for GAMESS input.",
    default=0,
    options={
        0: "Do not use template (default).",
        1: "Use template file.",
    },
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole/charges to gms_prop.ext files.",
    default=None,
    notes="Defaults to &cntrl namelist variable ntpr.",
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="chelpg",
    description="Calculate CHELPG charges and write to gms_prop.chg.",
    default=0,
    options={
        0: "Do not calculate CHELPG charges (default).",
        1: "Calculate CHELPG charges.",
    },
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipole",
    description="Toggle writing of dipole moment to gms_prop.dip.",
    default=0,
    options={
        0: "Do not write dipole (default).",
        1: "Write dipole moment.",
    },
    value_type="int",
    section="gms",
    category="QM/MM External - GAMESS",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &gau NAMELIST (Section 11.2.6.3 - AMBER/Gaussian interface)
# =============================================================================
# Gaussian 03/09/16 interface. Supports electrostatic embedding.

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for Gaussian. Any natively supported basis (STO-3G, 3-21G, 6-31G, 6-311G, etc.) with optional +/++ and */** extensions.",
    default="6-31G*",
    value_type="string",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="method",
    description="QM method for Gaussian. WFT methods (RHF, MP2) or DFT functionals (BLYP, PBE, B3LYP).",
    default="BLYP",
    value_type="string",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scf_conv",
    description="SCF convergence threshold for Gaussian in the form 10^-N.",
    default=8,
    value_type="int",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="num_threads",
    description="Number of threads (CPU cores) for Gaussian.",
    default=1,
    value_type="int",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="executable",
    description="Name of the Gaussian executable. If not specified, g16, g09, g03 are tried in order.",
    default=None,
    value_type="string",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a user-provided template file (gau_job.tpl) for Gaussian input.",
    default=0,
    options={
        0: "Do not use template (default).",
        1: "Use template file.",
    },
    value_type="int",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole moment to gau_job.dip.",
    default=None,
    notes="Defaults to &cntrl namelist variable ntpr.",
    value_type="int",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipole",
    description="Toggle writing of dipole moment to gau_job.dip.",
    default=0,
    options={
        0: "Do not write dipole (default).",
        1: "Write dipole moment.",
    },
    value_type="int",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mem",
    description="Memory allocation string for Gaussian.",
    default="256MB",
    value_type="string",
    section="gau",
    category="QM/MM External - Gaussian",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &orc NAMELIST (Section 11.2.6.4 - AMBER/Orca interface)
# =============================================================================
# Orca interface. Supports electrostatic embedding. Requires OpenMPI for parallel.

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for Orca (e.g. svp, 6-31g). See Orca manual for complete list.",
    default="SV(P)",
    value_type="string",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cbasis",
    description="Auxiliary basis set for correlation fitting in Orca.",
    default="NONE",
    value_type="string",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="jbasis",
    description="Auxiliary basis set for Coulomb fitting in Orca.",
    default="NONE",
    value_type="string",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="method",
    description="QM method for Orca (hf, pm3, blyp, mp2, etc.).",
    default="blyp",
    value_type="string",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="convkey",
    description="General SCF convergence keyword for simplified Orca input (TIGHTSCF, VERYTIGHTSCF, etc.).",
    default="VERYTIGHTSCF",
    value_type="string",
    section="orc",
    category="QM/MM External - Orca",
    related=["scfconv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scfconv",
    description="SCF energy convergence threshold for Orca. Energy converges to 10^-N au. -1 disables (uses convkey instead).",
    default=-1,
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=["convkey"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="grid",
    description="Grid type for XC quadrature during SCF in Orca DFT. Grid=4 uses IntAcc=4.34 with 302 Lebedev points.",
    default=4,
    notes="Conservatively chosen together with finalgrid to conserve energy.",
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=["finalgrid"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="finalgrid",
    description="Grid type for final energy/gradient XC quadrature in Orca DFT. Grid=6 uses IntAcc=5.34 with 590 Lebedev points.",
    default=6,
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=["grid"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxiter",
    description="Maximum number of SCF iterations for Orca.",
    default=100,
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxcore",
    description="Global scratch memory (MB) for Orca. Increase for larger jobs.",
    default=1024,
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="num_threads",
    description="Number of threads (CPU cores) for Orca. Note: Orca only supports OpenMPI.",
    default=1,
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a user-provided template file (orc_job.tpl) for Orca input.",
    default=0,
    options={
        0: "Do not use template (default).",
        1: "Use template file.",
    },
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole moment to orc_job.dip.",
    default=None,
    notes="Defaults to &cntrl namelist variable ntpr.",
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipole",
    description="Toggle writing of dipole moment to orc_job.dip.",
    default=0,
    options={
        0: "Do not write dipole (default).",
        1: "Write dipole moment.",
    },
    value_type="int",
    section="orc",
    category="QM/MM External - Orca",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &qc NAMELIST (Section 11.2.6.5 - AMBER/Q-Chem interface)
# =============================================================================
# Q-Chem interface. Supports electrostatic embedding.

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for Q-Chem. Default depends on method: 6-31G* for DFT, cc-pVDZ for MP2.",
    default="6-31G*",
    value_type="string",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="auxbasis",
    description="Auxiliary basis set for RI methods in Q-Chem. Default: rimp2-cc-pVDZ for RI-MP2, otherwise none.",
    default=None,
    value_type="string",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="method",
    description="QM method for Q-Chem (BLYP, MP2, RIMP2, or other functionals). Can use exchange/correlation keywords instead.",
    default="BLYP",
    value_type="string",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=["exchange", "correlation"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="exchange",
    description="Exchange method for Q-Chem. Can be used with correlation keyword in place of method.",
    default="",
    value_type="string",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=["method", "correlation"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="correlation",
    description="Correlation method for Q-Chem. Can be used with exchange keyword in place of method.",
    default="",
    value_type="string",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=["method", "exchange"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scf_conv",
    description="SCF convergence threshold for Q-Chem.",
    default=6,
    value_type="int",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="num_mpi_prcs",
    description="Number of MPI processes for Q-Chem. Total CPUs = num_mpi_prcs * num_threads.",
    default=1,
    value_type="int",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=["num_threads"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="num_threads",
    description="Number of threads per MPI process for Q-Chem. Total CPUs = num_mpi_prcs * num_threads.",
    default=1,
    value_type="int",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=["num_mpi_prcs"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a user-provided template file (qc_job.tpl) for Q-Chem input.",
    default=0,
    options={
        0: "Do not use template (default).",
        1: "Use template file.",
    },
    value_type="int",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole moment to qc_job.dip.",
    default=None,
    notes="Defaults to &cntrl namelist variable ntpr.",
    value_type="int",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipole",
    description="Toggle writing of dipole moment to qc_job.dip. Currently not supported.",
    default=0,
    value_type="int",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="guess",
    description="Toggle use of MOs from previous step as initial guess for SCF convergence in Q-Chem.",
    default="read",
    notes="Any string different from 'read' disables this.",
    value_type="string",
    section="qc",
    category="QM/MM External - Q-Chem",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &mrcc NAMELIST (Section 11.2.6.6 - AMBER/MRCC interface)
# =============================================================================
# MRCC interface. Supports electrostatic embedding and multilayer QM/QM/MM.

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for MRCC (6-31G*, cc-pVDZ, cc-pVTZ, etc.).",
    default="6-31G*",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="calc",
    description="Type of calculation for MRCC (SCF, B3LYP, MP2, CCSD(T), LCCSD(T), etc.).",
    default="SCF",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dft",
    description="DFT method specification for MRCC. 'off' for non-DFT calculations.",
    default="off",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["calc"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mem",
    description="Memory allocation for MRCC calculation.",
    default="256MB",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="verbosity",
    description="Verbosity of MRCC output file.",
    default=2,
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole moment to mrcc_job.dip.",
    default=None,
    notes="Defaults to &cntrl namelist variable ntpr.",
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="do_dipole",
    description="Toggle writing of dipole moment to mrcc_job.dip.",
    default=0,
    options={
        0: "Do not write dipole (default).",
        1: "Write dipole moment.",
    },
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nprintlog",
    description="Frequency of storing MRCC output files during minimization/MD.",
    default=0,
    notes="0 = keep only last output file.",
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="debug",
    description="Toggle debug mode for AMBER/MRCC interface.",
    default=0,
    options={
        0: "No debugging (default).",
        1: "Print subroutine calls and additional info.",
    },
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a template MINP file for MRCC input.",
    default=0,
    options={
        0: "Do not use template (default).",
        1: "Use template file (mrcc_job.tpl).",
    },
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="embed",
    description="Method of embedding QM region (2nd layer) in multilayer QM/QM/MM or 3rd layer in QM/QM/QM/MM.",
    default="off",
    notes="See MRCC manual for available options.",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["embedatoms", "nmo_embed"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="embedatoms",
    description="Active atoms of embedded QM region in multilayer calculations. Comma-separated list of integers.",
    default="0",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["embed", "nmo_embed"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nmo_embed",
    description="Number of active MOs for embedded QM region. 0 = automatic via Boughton-Pulay algorithm.",
    default=0,
    options={
        0: "Automatic MO determination via BP algorithm (default).",
    },
    notes="Values > 0: MOs selected based on Mulliken charges.",
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["embed", "embedatoms"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="corembed",
    description="Low-level correlation method of embedding region in multilayer calculations.",
    default="off",
    notes="See MRCC manual for available options.",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["corembedatoms", "nmo_corembed"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="corembedatoms",
    description="Active atoms for core embedded region. Must be subset of embedatoms for 4-layer calculations.",
    default="0",
    value_type="string",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["corembed", "embedatoms"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nmo_corembed",
    description="Number of active MOs for core embedded region. 0 = automatic via BP algorithm.",
    default=0,
    options={
        0: "Automatic MO determination via BP algorithm (default).",
    },
    value_type="int",
    section="mrcc",
    category="QM/MM External - MRCC",
    related=["corembed", "corembedatoms"],
    commonly_changed=False,
))


# =============================================================================
# &fb NAMELIST (Section 11.2.6.7 - AMBER/Fireball interface)
# =============================================================================
# Fireball DFT tight-binding interface. Requires special sander compilation
# with Intel compilers and MKL linked against libfireball.a.

KEYWORDS.append(Keyword(
    name="basis",
    description="Path to the Fdata directory containing Fireball basis set interactions.",
    default="./Fdata",
    notes="Can be downloaded from fireball-qmd github.",
    value_type="string",
    section="fb",
    category="QM/MM External - Fireball",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="max_scf_iterations",
    description="Maximum number of SCF iterations for self-consistent charges in Fireball.",
    default=70,
    value_type="int",
    section="fb",
    category="QM/MM External - Fireball",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="sigmatol",
    description="Threshold for self-consistency in Fireball electronic structure calculations.",
    default=1e-08,
    value_type="float",
    section="fb",
    category="QM/MM External - Fireball",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="idftd3",
    description="DFTD3 dispersion correction for Fireball.",
    default=0,
    options={
        0: "No dispersion correction (default).",
        1: "Dispersion correction for BLYP.",
    },
    value_type="int",
    section="fb",
    category="QM/MM External - Fireball",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="iwrtcharges",
    description="Write atomic charges in Fireball output.",
    default=0,
    options={
        0: "Do not write charges (default).",
        1: "Write charges.",
    },
    value_type="int",
    section="fb",
    category="QM/MM External - Fireball",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="iwrteigen",
    description="Write energy levels in Fireball output.",
    default=0,
    options={
        0: "Do not write eigenvalues (default).",
        1: "Write eigenvalues.",
    },
    value_type="int",
    section="fb",
    category="QM/MM External - Fireball",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &quick NAMELIST (Section 11.3.1.3 - AMBER/QUICK interface)
# =============================================================================
# QUICK (QUantum Interaction Computational Kernel) GPU-enabled ab initio QM.
# Supports both API (linked library, recommended) and FBI (file-based) interfaces.
# Distributed with AmberTools. Supports HF and DFT with electrostatic/mechanical embedding.

KEYWORDS.append(Keyword(
    name="method",
    description="QM method for QUICK. Can be HF or a supported DFT functional.",
    default="BLYP",
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=["basis"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for QUICK (e.g. 6-31G, def2-svp, cc-pVDZ).",
    default="6-31G",
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=["method"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="executable",
    description="QUICK executable for FBI. Options: quick, quick.MPI, quick.cuda, quick.hip, quick.cuda.MPI, quick.hip.MPI.",
    default="quick",
    notes="FBI only.",
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="do_parallel",
    description="Parallel launch command for FBI (e.g. 'mpirun -np 2'). Only needed for MPI versions.",
    default=None,
    notes="FBI only. The exact command depends on the MPI library installed.",
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=["executable"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scf_cyc",
    description="Number of SCF cycles for QUICK.",
    default=200,
    value_type="int",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="reuse_dmx",
    description="Reuse density matrix from previous MD step (API only).",
    default=1,
    options={
        0: "Do not reuse density matrix.",
        1: "Reuse density matrix (default, API only).",
    },
    value_type="int",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="denserms",
    description="User-defined density matrix maximum RMS for convergence (API only).",
    default=1e-06,
    value_type="float",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="intcutoff",
    description="User-defined integral cutoff (API only).",
    default=1e-08,
    value_type="float",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="xccutoff",
    description="User-defined threshold for grid pruning in exchange-correlation algorithm (API only).",
    default=1e-08,
    value_type="float",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="basiscutoff",
    description="Cutoff for neglecting insignificant basis functions (API only).",
    default=1e-06,
    value_type="float",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="gradcutoff",
    description="User-defined gradient cutoff (API only).",
    default=1e-07,
    value_type="float",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="export",
    description="Export molecular orbitals and QM data. Currently supports 'molden' only (API only).",
    default=None,
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="keywords",
    description="Full QUICK keywords line instead of separate flags (API only). E.g. 'B3LYP BASIS=cc-pVDZ CHARGE=0 MULT=1 GRADIENT EXTCHARGES'.",
    default=None,
    notes="API only. Overrides individual method/basis settings.",
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="outfprefix",
    description="Prefix for QUICK output file (API only). Name will be followed by '.out'.",
    default="quick",
    value_type="string",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="debug",
    description="Debugging information level for QUICK.",
    default=0,
    options={
        0: "No debugging (default).",
        1: "Print debugging information.",
        2: "Extra debugging for FBI.",
    },
    value_type="int",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a template input file for QUICK FBI.",
    default=0,
    options={
        0: "No template (default, FBI only).",
        1: "Use template file (FBI only).",
    },
    value_type="int",
    section="quick",
    category="QM/MM QUICK",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &tc NAMELIST (Section 11.4.1.3 - AMBER/TeraChem interface)
# =============================================================================
# TeraChem GPU-accelerated QM. Supports both client/server (TCPB, recommended)
# and file-based (FBI) interfaces. Supports HF, DFT, CI, CASSCF.

KEYWORDS.append(Keyword(
    name="host",
    description="Address of machine where TeraChem server is hosted (TCPB only).",
    default=None,
    notes="Required for client/server interface.",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=["port"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="port",
    description="Port number used by TeraChem server (TCPB only).",
    default=None,
    notes="Required for client/server interface.",
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=["host"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="tcfile",
    description="TeraChem input file to pass to TCPB-cpp (TCPB only). If specified with method, tcfile becomes output.",
    default=None,
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="method",
    description="QM method for TeraChem (RHF, BLYP, PBE, B3LYP, etc.).",
    default=None,
    notes="Default: none for TCPB, BLYP for FBI.",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=["basis"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="basis",
    description="Basis set for TeraChem (STO-3G, 3-21G, 6-31G, 6-311G, etc.).",
    default="6-31G",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=["method"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dftd",
    description="Dispersion corrections for DFT in TeraChem.",
    default="no",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="precision",
    description="Precision model (single vs double) for TeraChem.",
    default="mixed",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="guess",
    description="Path to initial wavefunction guess file (FBI only).",
    default="scr/c0",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scrdir",
    description="Path to TeraChem scratch directory (FBI only).",
    default="scr",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="keep_scr",
    description="Keep only a single scratch directory (FBI only).",
    default="yes",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=["scrdir"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="threall",
    description="General threshold controlling various cutoffs in TeraChem.",
    default=1e-11,
    value_type="float",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="convthre",
    description="SCF wavefunction convergence threshold for TeraChem.",
    default=3e-05,
    notes="Leads to SCF energy convergence of approximately 10^-7 au.",
    value_type="float",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxit",
    description="Maximum number of SCF iterations for TeraChem.",
    default=100,
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cis",
    description="Perform CIS (Configuration Interaction Singles) calculation.",
    default="no",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=["cisnumstates", "cistarget"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cisnumstates",
    description="Number of CIS excited states to compute.",
    default=1,
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=["cis", "cistarget"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="cistarget",
    description="Target CIS state for forces/gradient.",
    default=1,
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=["cis", "cisnumstates"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dftgrid",
    description="DFT grid for numerical XC quadrature.",
    default=1,
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ngpus",
    description="Number of GPUs for TeraChem (FBI only). 0 = use all available.",
    default=0,
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=["gpuids"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="gpuids",
    description="GPU IDs for TeraChem when ngpus != 0 (FBI only).",
    default=None,
    notes="Default: 0, 1, 2, etc.",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=["ngpus"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="executable",
    description="TeraChem executable name (FBI only).",
    default="terachem",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="use_template",
    description="Use a template input file (tc_job.tpl) for TeraChem FBI.",
    default=0,
    options={
        0: "No template (default, FBI only).",
        1: "Use template file (FBI only).",
    },
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntpr",
    description="Frequency of printing dipole/charges to tc_job.ext files (FBI only).",
    default=None,
    notes="Defaults to &cntrl ntpr.",
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="charge_analysis",
    description="Toggle writing atomic charges to tc_job.chg (FBI only).",
    default="none",
    notes="Options: 'none' or 'Mulliken'.",
    value_type="string",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dipole",
    description="Toggle writing dipole moment to tc_job.dip (FBI only).",
    default=0,
    options={
        0: "Do not write dipole (default).",
        1: "Write dipole.",
    },
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="recycleinitguess",
    description="Reuse wavefunction from previous MD step as initial guess.",
    default=1,
    options={
        0: "Do not reuse wavefunction.",
        1: "Reuse wavefunction (default).",
    },
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="debug",
    description="Debugging information level for TeraChem interface.",
    default=0,
    options={
        0: "No debugging (default).",
        1: "Print debugging info.",
        2: "Extra debugging for FBI.",
    },
    value_type="int",
    section="tc",
    category="QM/MM TeraChem",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &xtb NAMELIST (Section 11.5.1 - AMBER/xTB interface)
# =============================================================================
# xTB tight-binding DFT interface. Supports GFN1-xTB and GFN2-xTB.
# Fully supports QM/MM Ewald electrostatics within SCF.
# Not distributed with AmberTools; requires separate xTB library compilation.

KEYWORDS.append(Keyword(
    name="qm_level",
    description="xTB method to use (GFN1-xTB or GFN2-xTB).",
    default="GFN2-xTB",
    value_type="string",
    section="xtb",
    category="QM/MM xTB",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="tfermi",
    description="Electronic Fermi temperature in Kelvin for xTB.",
    default=300.0,
    value_type="float",
    section="xtb",
    category="QM/MM xTB",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="accuracy",
    description="Controls SCF convergence criteria for xTB. Energy tolerance = 1e-6 * accuracy.",
    default=0.001,
    value_type="float",
    section="xtb",
    category="QM/MM xTB",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxiter",
    description="Maximum number of SCF iterations for xTB.",
    default=250,
    value_type="int",
    section="xtb",
    category="QM/MM xTB",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mmhardness",
    description="Chemical hardness of MM atoms (atomic units) for tuning QM/MM electrostatics in xTB.",
    default=0.0,
    notes="0.0 = hydrogen hardness (0.405771 au). Positive = uniform hardness for all MM atoms. Negative = element-based hardness scaled by abs(mmhardness).",
    value_type="float",
    section="xtb",
    category="QM/MM xTB",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="debug",
    description="Turn on debug printing for xTB.",
    default="F",
    value_type="string",
    section="xtb",
    category="QM/MM xTB",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &dftbplus NAMELIST (Section 11.6.1 - AMBER/DFTB+ interface)
# =============================================================================
# DFTB+ tight-binding DFT interface. Supports DFTB2, DFTB3, DFTB3-D3,
# DFTB3-D3H5, DFTB3-D42B, DFTB3-D43B.
# Supports QM/MM Ewald electrostatics within SCF using Mulliken charges.
# Not distributed with AmberTools; requires separate DFTB+ library compilation.

KEYWORDS.append(Keyword(
    name="qm_level",
    description="DFTB+ method. DFTB3 includes full third-order corrections; DFTB3-D3 adds Grimme D3 dispersion with BJ damping.",
    default="DFTB2",
    notes="DFTB3-D3H5, DFTB3-D42B, DFTB3-D43B require DFTB+ compiled with -DWITH_SDFTD3=TRUE.",
    value_type="string",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="tfermi",
    description="Electronic Fermi temperature in Kelvin for DFTB+.",
    default=0.0,
    value_type="float",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="scftol",
    description="SCF convergence criteria for DFTB+.",
    default=1e-07,
    value_type="float",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="maxiter",
    description="Maximum number of SCF iterations for DFTB+.",
    default=250,
    value_type="int",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="hcorrection",
    description="Hydrogen-bond correction for DFTB+.",
    default=-1,
    options={
        -1: "Auto: on for DFTB3/DFTB3-D3, off for DFTB2 (default).",
        0: "Turn off H-bond correction.",
        1: "Turn on H-bond correction.",
    },
    value_type="int",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="debug",
    description="Turn on debug printing for DFTB+.",
    default="F",
    value_type="string",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="silent",
    description="Suppress detailed DFTB+ output to disk.",
    default="T",
    notes="When False, writes dftb_pin.hsd, charges.bin, and detailed.out at every timestep.",
    value_type="string",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mixer",
    description="SCF iteration mixing algorithm for DFTB+.",
    default="BROYDEN",
    notes="Options: BROYDEN, DIIS, or SIMPLE.",
    value_type="string",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="d3a1",
    description="The a1 parameter (unitless) in Becke-Johnson dispersion damping for DFTB+.",
    default=-1.0,
    notes="If < 0, value chosen automatically based on qm_level.",
    value_type="float",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="d3a2",
    description="The a2 parameter (Bohr) in Becke-Johnson dispersion damping for DFTB+.",
    default=-1.0,
    notes="If < 0, value chosen automatically based on qm_level.",
    value_type="float",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=["d3a1"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="d3s8",
    description="The C8 dispersion scale factor (unitless) in D3 dispersion for DFTB+.",
    default=-1.0,
    notes="If < 0, value chosen automatically based on qm_level.",
    value_type="float",
    section="dftbplus",
    category="QM/MM DFTB+",
    related=["d3a1", "d3a2"],
    commonly_changed=False,
))


# =============================================================================
# &dprc NAMELIST (Section 11.7.1 - DPRc ML corrections)
# =============================================================================
# Deep Potential Range Corrected (DPRc) neural network correction for QM/MM.
# Requires DeePMD-kit library (not distributed with Amber).
# Trains nonelectrostatic corrections to improve QM/MM interactions.

KEYWORDS.append(Keyword(
    name="idprc",
    description="Activate DPRc correction.",
    default=0,
    options={
        0: "DPRc correction disabled (default).",
        1: "DPRc correction enabled.",
    },
    value_type="int",
    section="dprc",
    category="QM/MM DPRc",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="mask",
    description="Amber atom selection for DPRc correction. Normally set to same as qmmask.",
    default="",
    value_type="string",
    section="dprc",
    category="QM/MM DPRc",
    related=["qmmask", "rcut"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="rcut",
    description="Cutoff (Angstroms) for dynamically selecting nearby atoms for DPRc correction.",
    default=0.0,
    notes="QM/MM interactions within rcut of the QM region are corrected. 0.0 = correct only QM/QM interactions. Must be consistent with the neural network model.",
    value_type="float",
    section="dprc",
    category="QM/MM DPRc",
    related=["mask"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="interfile",
    description="List of up to 4 DeePMD-kit neural network parameter files for QM/QM and QM/MM corrections.",
    default=None,
    notes="For inference use interfile(1) only. For active learning provide 4 files; deviation reported every ntwx steps.",
    value_type="string",
    section="dprc",
    category="QM/MM DPRc",
    related=["intrafile", "avg"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="intrafile",
    description="List of up to 4 DeePMD-kit neural network parameter files for QM/QM only corrections.",
    default=None,
    notes="Analogous to interfile with rcut=0.0. Useful for gas-phase trained corrections combined with condensed-phase interfile.",
    value_type="string",
    section="dprc",
    category="QM/MM DPRc",
    related=["interfile", "avg"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="avg",
    description="Average corrections from 4 models instead of using only the first.",
    default="F",
    notes="Only relevant when 4 interfile/intrafile parameter files are provided.",
    value_type="string",
    section="dprc",
    category="QM/MM DPRc",
    related=["interfile", "intrafile"],
    commonly_changed=False,
))


# =============================================================================
# &sebomd NAMELIST (Section 11.10.3 - SEBOMD)
# =============================================================================
# SemiEmpirical Born-Oppenheimer Molecular Dynamics. All atoms treated as QM
# within the NDDO semiempirical approach. No link atoms or QM/MM boundary.
# Supports optional divide-and-conquer linear scaling.
# Max atoms: 1000. Max residues: 1000 (static allocation).

KEYWORDS.append(Keyword(
    name="charge",
    description="Net charge of the full system for SEBOMD. Only closed-shell systems supported.",
    default=0,
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=["hamiltonian"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="method",
    description="Algorithm for SCF computation in SEBOMD.",
    default=0,
    options={
        0: "Standard closed-shell algorithm with full Fock matrix diagonalization at each SCF iteration (default).",
        1: "Divide and conquer linear scaling SCF. Atom-based subsystems. Requires dbuff1, dbuff2.",
        2: "Divide and conquer linear scaling SCF. Residue-based subsystems (recommended over method=1).",
        3: "Divide and conquer linear scaling SCF. Heavy-atom-based subsystems.",
    },
    notes="Parallel execution (sander.MPI) only supported for method > 0.",
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=["ncore", "dbuff1", "dbuff2"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ncore",
    description="Number of subsystems used to build the core in divide-and-conquer SEBOMD (method > 0).",
    default=1,
    value_type="int",
    min_val=1,
    section="sebomd",
    category="SEBOMD",
    related=["method"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dbuff1",
    description="Extent of the first buffer region from the core (Angstroms) in divide-and-conquer SEBOMD.",
    default=6.0,
    value_type="float",
    section="sebomd",
    category="SEBOMD",
    related=["method", "dbuff2"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dbuff2",
    description="Extent of the second buffer region from the core (Angstroms) in divide-and-conquer SEBOMD.",
    default=0.0,
    value_type="float",
    section="sebomd",
    category="SEBOMD",
    related=["method", "dbuff1"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="hamiltonian",
    description="Semiempirical Hamiltonian for SEBOMD energy and force calculations.",
    default="PM3",
    options={
        "MNDO": "MNDO semiempirical Hamiltonian.",
        "AM1": "AM1 semiempirical Hamiltonian.",
        "PM3": "PM3 semiempirical Hamiltonian (default).",
        "PM3PDDG": "PM3/PDDG semiempirical Hamiltonian.",
        "RM1": "RM1 semiempirical Hamiltonian.",
        "AM1/d-PhoT": "AM1/d-PhoT Hamiltonian (P element not yet implemented).",
    },
    notes="Only H, C, N, O, P, S, F, Cl, Br, I elements implemented (no d-orbitals).",
    value_type="string",
    section="sebomd",
    category="SEBOMD",
    related=["modif"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="modif",
    description="Modification/correction to the semiempirical energy in SEBOMD.",
    default="none",
    options={
        "none": "No correction (default).",
        "PIF2": "PM3-PIF intermolecular core-core correction for organic-water systems.",
        "PIF3": "Extended PIF2 distinguishing hydrophilic vs hydrophobic hydrogens.",
        "MAIS1": "MAIS extension replacing intramolecular and intermolecular core-core functions. H and O only.",
        "MAIS2": "Second MAIS version. H, O, and Cl only.",
    },
    notes="PIF and MAIS corrections are only available for the PM3 Hamiltonian.",
    value_type="string",
    section="sebomd",
    category="SEBOMD",
    related=["hamiltonian"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="longrange",
    description="Long-range interaction treatment for periodic boundary conditions in SEBOMD.",
    default=0,
    options={
        0: "No long-range interaction. Minimum image convention only (default).",
        1: "PME using constant atomic charges from topology file.",
        2: "Ewald summation using Mulliken charges from the semiempirical wavefunction. Long-range effects included in Fock matrix.",
    },
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="dpmax",
    description="SCF convergence criteria on the density matrix for SEBOMD.",
    default=1e-07,
    notes="Default ensures energy conservation in NVE. Larger values speed up but may break conservation.",
    value_type="float",
    section="sebomd",
    category="SEBOMD",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="fullscf",
    description="Control pseudo-diagonalization in SEBOMD SCF.",
    default=0,
    options={
        0: "Enable pseudo-diagonalization when possible (faster, default).",
        1: "Full diagonalization of Fock matrix at each SCF iteration.",
    },
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ipolyn",
    description="Polynomial interpolation of guess density matrix in SEBOMD.",
    default=1,
    options={
        0: "Use previous step converged density as guess. Recommended for minimization.",
        1: "Polynomial interpolation from last 3 steps. Recommended for MD (default).",
    },
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="screen",
    description="Verbosity option for SEBOMD calculations.",
    default=0,
    options={
        0: "Minimum output (default).",
        1: "Output semiempirical energy details at each step.",
        2: "Energy details plus subsystem composition when method > 0.",
    },
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="lambda",
    description="Mixing parameter between SEBOMD and MM energy. E = lambda*E(SEBOMD) + (1-lambda)*E(MM).",
    default=1.0,
    notes="Useful for gradual equilibration from MM (lambda=0) to full QM (lambda=1).",
    value_type="float",
    min_val=0.0,
    max_val=1.0,
    section="sebomd",
    category="SEBOMD",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="charge_out",
    description="Filename for saving SEBOMD atomic charges.",
    default="sebomd.chg",
    value_type="string",
    section="sebomd",
    category="SEBOMD",
    related=["ntwc"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="ntwc",
    description="Frequency of writing Mulliken charges to charge_out file. 0 = no output.",
    default=0,
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=["charge_out"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="peptcorr",
    description="Apply force field correction on peptidic bonds to enforce planarity in SEBOMD.",
    default=0,
    options={
        0: "No peptidic correction (default).",
        1: "Apply peptidic correction.",
    },
    value_type="int",
    section="sebomd",
    category="SEBOMD",
    related=["peptk"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="peptk",
    description="Force constant (kcal/mol) for peptidic correction in SEBOMD.",
    default=None,
    notes="Defaults depend on Hamiltonian: AM1=5.9864, PM3=9.8526, MNDO=6.1737.",
    value_type="float",
    section="sebomd",
    category="SEBOMD",
    related=["peptcorr"],
    commonly_changed=False,
))


# =============================================================================
# &vsolv NAMELIST (Section 11.8.2.2 - Adaptive solvent parameters)
# =============================================================================
# Controls which solvent molecules are included in the QM region for
# adaptive solvent QM/MM. Used with vsolv=1,2,3 in &qmmm namelist.

KEYWORDS.append(Keyword(
    name="nearest_qm_solvent_resname",
    description="Residue name of solvent that can exchange between QM and MM regions.",
    default="WAT",
    value_type="string",
    section="vsolv",
    category="QM/MM Adaptive Solvent",
    related=["vsolv", "nearest_qm_solvent"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nearest_qm_solvent",
    description="Number of solvent molecules in the active (A) region.",
    default=0,
    value_type="int",
    min_val=0,
    section="vsolv",
    category="QM/MM Adaptive Solvent",
    related=["vsolv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nearest_qm_solvent_fq",
    description="Frequency of updating the active region. Should be 1 (every step) for adQM/MM.",
    default=1,
    value_type="int",
    min_val=1,
    section="vsolv",
    category="QM/MM Adaptive Solvent",
    related=["vsolv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nearest_qm_solvent_center_id",
    description="Determines which atom of solvent molecules is used for distance calculation to QM region.",
    default=0,
    options={
        0: "Use the closest atom to the QM region (default).",
        -1: "Use the center of mass of the solvent molecule.",
    },
    notes="Values > 0 specify a particular atom number within the solvent residue.",
    value_type="int",
    section="vsolv",
    category="QM/MM Adaptive Solvent",
    related=["vsolv", "qm_center_atom_id"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qm_center_atom_id",
    description="Determines the atom of the permanent QM region used for distance calculation to solvent.",
    default=0,
    options={
        0: "Use the closest atom (not supported for adQM/MM with radii, default).",
        -1: "Use center of mass of permanent QM region.",
    },
    notes="Values > 0 specify an absolute atom number in the permanent QM region. For adQM/MM, 0 is not supported since a common reference point is needed.",
    value_type="int",
    section="vsolv",
    category="QM/MM Adaptive Solvent",
    related=["vsolv", "nearest_qm_solvent_center_id"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="verbosity",
    description="Verbosity of vsolv output in the mdout file.",
    default=0,
    notes="Values > 1 increase verbosity.",
    value_type="int",
    section="vsolv",
    category="QM/MM Adaptive Solvent",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# &adqmmm NAMELIST (Section 11.8.2.2 - Adaptive QM/MM parameters)
# =============================================================================
# Controls adaptive solvent QM/MM simulations with vsolv=2 or 3 in &qmmm.
# Requires multisander with groupfile for parallel execution.

KEYWORDS.append(Keyword(
    name="n_partition",
    description="Number of QM/MM partitions for adQM/MM. For vsolv=2 this also sets NT = n_partition - 1 solvent molecules in the transition region.",
    default=1,
    notes="For vsolv=3, set to the largest number of partitionings expected for chosen RA and RT.",
    value_type="int",
    min_val=1,
    section="adqmmm",
    category="QM/MM Adaptive Solvent",
    related=["vsolv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="RA",
    description="Radius (Angstroms) of the active (A) region. Only relevant for vsolv=3.",
    default=-1.0,
    notes="Must be changed from default. Requires setting RT.",
    value_type="float",
    section="adqmmm",
    category="QM/MM Adaptive Solvent",
    related=["RT", "vsolv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="RT",
    description="Radius (Angstroms) of the transition (T) region. Only relevant for vsolv=3.",
    default=-1.0,
    notes="Must be changed from default. Requires setting RA.",
    value_type="float",
    section="adqmmm",
    category="QM/MM Adaptive Solvent",
    related=["RA", "vsolv"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="calc_wbk",
    description="Controls whether the book-keeping term W is calculated for adQM/MM.",
    default=0,
    options={
        0: "Do not calculate W (default).",
        1: "Calculate W via one-sided difference approximation (not recommended).",
        2: "Calculate W via central-difference approximation (recommended if W desired).",
    },
    value_type="int",
    section="adqmmm",
    category="QM/MM Adaptive Solvent",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="verbosity",
    description="Verbosity of adQM/MM output.",
    default=0,
    options={
        0: "Standard verbosity (default).",
        1: "Write distances to adqmmm_res_distances.dat and weights to adqmmm_weights.dat.",
        2: "Write distances, weights, and lambda values also to mdout file.",
    },
    value_type="int",
    section="adqmmm",
    category="QM/MM Adaptive Solvent",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="print_qm_coords",
    description="Controls writing of QM atom coordinates for each partitioning.",
    default=0,
    options={
        0: "Do not write coordinates (default).",
        1: "Write QM coordinates in xyz format to QM_coords.001 etc.",
    },
    value_type="int",
    section="adqmmm",
    category="QM/MM Adaptive Solvent",
    related=[],
    commonly_changed=False,
))


# =============================================================================
# Chapter 29: Continuous Constant pH MD Namelists
# =============================================================================

# ---------------------------------------------------------------------------
# &phmdin NAMELIST (Section 29.2.1 - Continuous CpHMD input parameters)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="nsolute",
    description=(
        "Number of solute residues in the system. This variable is no longer used "
        "in current versions but may be present in older input files for backward compatibility."
    ),
    default=None,
    notes="Deprecated; not used anymore.",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="phmdcut",
    description=(
        "The cutoff distance in angstroms to use in Generalized Born calculations during "
        "the continuous constant pH simulation. This should in general be very large, if not "
        "large enough to encompass the whole system, to ensure accurate GB energetics for "
        "protonation state evaluation."
    ),
    default=1000.0,
    notes="Only relevant for iphmd=1 (implicit solvent) and iphmd=2 (hybrid solvent) modes.",
    value_type="float",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["iphmd"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qmass_phmd",
    description=(
        "The mass of the virtual (fictitious) particles associated with the titration "
        "coordinates (lambda), in atomic mass units (amu). Should be roughly as large as "
        "the largest masses in the system."
    ),
    default=10.0,
    notes="Controls the inertia of the lambda dynamics. Affects the timescale of protonation state transitions.",
    value_type="float",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["temp_phmd", "phbeta"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="temp_phmd",
    description=(
        "The temperature of the virtual particles associated with the titration coordinates, "
        "in Kelvin. Controls the thermal fluctuations of the lambda particles."
    ),
    default=300.0,
    notes="Should generally match the simulation temperature (temp0) for consistency.",
    value_type="float",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["qmass_phmd", "phbeta", "temp0"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="phbeta",
    description=(
        "The friction constant for the Langevin integrator used to propagate the titration "
        "coordinates (theta variables), in ps^-1."
    ),
    default=5.0,
    notes="Controls the damping of the lambda dynamics Langevin thermostat.",
    value_type="float",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["qmass_phmd", "temp_phmd"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="iphfrq",
    description=(
        "The number of MD steps between updates of the titration coordinates (lambda). "
        "Controls how frequently the lambda particles are propagated."
    ),
    default=1,
    notes="A value of 1 means titration coordinates are updated every MD step, which is the recommended setting.",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=[],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="qphmdstart",
    description=(
        "Controls whether the velocities of the virtual (lambda) particles should be "
        "regenerated from the Boltzmann distribution at the start of the simulation. "
        "If false, velocities are read from the restart file (-phmdstrt)."
    ),
    default=True,
    options={
        True: "Regenerate velocities from Boltzmann distribution (default).",
        False: "Read velocities from the phmdstrt restart file.",
    },
    notes=(
        "Should be set to false (.false.) for thermodynamic integration (TI) simulations "
        "used during parameterization of new titratable residues for CpHMD."
    ),
    value_type="bool",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["phtest"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="nprint_phmd",
    description=(
        "The number of MD steps between prints to the lambda file (-phmdout). Controls "
        "how frequently the lambda values and related data are written to the output file."
    ),
    default=None,
    notes="Default is set to the same as the trajectory printing frequency (ntwx).",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Output",
    related=["prlam", "ntwx"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="prlam",
    description=(
        "Controls whether the lambda values are printed to the lambda file (-phmdout). "
        "Must be set to true to generate lambda output for analysis of protonation states."
    ),
    default=False,
    options={
        True: "Print lambda values to the lambda file.",
        False: "Do not print lambda values (default).",
    },
    notes="Set to .true. for production simulations where protonation state analysis is needed.",
    value_type="bool",
    section="phmdin",
    category="Continuous CpHMD - Output",
    related=["nprint_phmd", "prderiv"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="prderiv",
    description=(
        "Controls whether the theta values and dU/dtheta derivative information are "
        "output to the mdout file. Used for parameterization of new titratable residues."
    ),
    default=False,
    options={
        True: "Print theta and dU/dtheta to mdout.",
        False: "Do not print derivative information (default).",
    },
    notes="Enable during thermodynamic integration (TI) simulations for deriving model compound parameters.",
    value_type="bool",
    section="phmdin",
    category="Continuous CpHMD - Output",
    related=["phtest", "prlam"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="prnlev",
    description=(
        "Print level controlling what gets printed during continuous CpHMD simulations. "
        "Higher values produce more diagnostic output."
    ),
    default=6,
    options={
        0: "Print header information in the output file.",
        2: "Values > 2 generate the full output file.",
        5: "Values >= 5 print additional diagnostic data to the mdout file.",
    },
    notes="The options describe thresholds: >= 0 prints headers, > 2 prints full output, >= 5 prints diagnostics.",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Output",
    related=["outu"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="outu",
    description=(
        "The Fortran unit number for printing continuous CpHMD diagnostic information."
    ),
    default=6,
    notes="Unit 6 is standard output. Typically does not need to be changed.",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Output",
    related=["prnlev"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="phtest",
    description=(
        "If equal to 1, the theta and theta_x variables are held fixed during the simulation. "
        "Used for parameterization of new titratable residues via thermodynamic integration."
    ),
    default=0,
    options={
        0: "Normal dynamics; theta variables are propagated (default).",
        1: "Fix theta and theta_x variables (for TI parameterization).",
    },
    notes="Set to 1 during TI simulations to compute mean forces at fixed lambda values.",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["prderiv", "qphmdstart"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="masktitrres",
    description=(
        "An array specifying the names of the titratable residues in the system. These names "
        "must match residue names defined in the phmdparm file. For AMBER force fields, typical "
        "names are AS2, GL2, HIP, LYS, and CYS. For CHARMM22, use ASP, GLU, HIP, CYS, and LYS."
    ),
    default=None,
    notes=(
        "Specified as a Fortran array, e.g.: MaskTitrRes(:) = 'AS2','GL2','HIP','CYS'. "
        "Only residues listed here will be treated as titratable during the simulation."
    ),
    value_type="string_array",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["masktitrrestypes", "ngt", "res_name"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="masktitrrestypes",
    description=(
        "The number of distinct titratable residue type entries in the masktitrres array. "
        "Must match the number of entries provided in masktitrres."
    ),
    default=None,
    notes="For example, if masktitrres has 4 entries (AS2, GL2, HIP, CYS), set masktitrrestypes = 4.",
    value_type="int",
    section="phmdin",
    category="Continuous CpHMD - Input",
    related=["masktitrres"],
    commonly_changed=True,
))

# ---------------------------------------------------------------------------
# &phmdparm NAMELIST (Section 29.2.2 - Continuous CpHMD parameter file)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ngt",
    description=(
        "The number of titratable residue types defined in the phmdparm file. Each residue "
        "type has its own set of atom names, charges, disappearing atom flags, model pKa values, "
        "model potential parameters, and barrier heights."
    ),
    default=None,
    notes="Must match the number of entries in res_name, numch, and res_type arrays.",
    value_type="int",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["numch", "res_name", "res_type"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="numch",
    description=(
        "An array of the number of atoms in each titratable residue type defined in the "
        "phmdparm file. The length of this array must equal ngt."
    ),
    default=None,
    notes="Example: NUMCH(:) = 14,17,18,22,11 for AS2(14), GL2(17), HIP(18), LYS(22), CYS(11).",
    value_type="int_array",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["ngt", "atom_name"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="res_name",
    description=(
        "An array of the names of the titratable residue types in the phmdparm file. "
        "These names are used to match residues in the topology to their titration parameters."
    ),
    default=None,
    notes="Example: RES_NAME(:) = 'AS2','GL2','HIP','LYS','CYS'.",
    value_type="string_array",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["ngt", "res_type", "masktitrres"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="res_type",
    description=(
        "An array defining the titration type of each titratable residue in the phmdparm file. "
        "The type determines how the lambda and x variables are used to interpolate between "
        "protonation states."
    ),
    default=None,
    options={
        -2: "Coions titrating with linked titratable residues to maintain constant charge (not currently used).",
        0: "Residues with a single titratable hydrogen (e.g., lysine, cysteine).",
        2: (
            "Residues with two deprotonated states (tautomers) and a single protonated state, "
            "with the two deprotonated states having different pKa values (e.g., histidine)."
        ),
        4: (
            "Residues with two protonated states and a single deprotonated state, where the two "
            "protonated states have the same pKa (e.g., aspartic acid, glutamic acid)."
        ),
    },
    notes="Example: RES_TYPE(:) = 4,4,2,0,0 for AS2(4), GL2(4), HIP(2), LYS(0), CYS(0).",
    value_type="int_array",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["ngt", "res_name"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="atom_name",
    description=(
        "A two-dimensional array containing the force field atom names for each titratable "
        "residue type defined in the phmdparm file. ATOM_NAME(i,:) lists the atom names for "
        "residue type i, and must have exactly numch(i) entries."
    ),
    default=None,
    notes="Example: ATOM_NAME(1,:) = 'N','H','CA','HA','CB','HB2','HB3','CG','OD1','OD2','HD2','C','O','HD1'",
    value_type="string_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["numch", "ch", "ch_md", "rad"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ch",
    description=(
        "A two-dimensional array of partial charges used in the dynamics of the titration "
        "coordinates (lambda forces). CH(i,:) contains the charge sets for residue type i. "
        "For type 0/-2 residues: protonated charges then deprotonated charges. "
        "For type 2: protonated charges then the two deprotonated state charges. "
        "For type 4: two protonated state charges then deprotonated charges."
    ),
    default=None,
    notes=(
        "These charges are used to compute forces on lambda particles. They may differ from "
        "ch_md charges which are used for spatial force calculations."
    ),
    value_type="float_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["ch_md", "atom_name", "res_type"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="ch_md",
    description=(
        "A two-dimensional array of partial charges used for the calculation of the spatial "
        "forces (i.e., the forces on the physical atoms). Same structure as the ch array: "
        "CH_MD(i,:) contains the charge sets for residue type i, organized by protonation state "
        "according to the res_type of the residue."
    ),
    default=None,
    notes=(
        "These charges drive the conformational dynamics. They may differ from the ch array "
        "charges which drive the lambda dynamics."
    ),
    value_type="float_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["ch", "atom_name", "res_type"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="rad",
    description=(
        "A two-dimensional array of flags identifying which atoms disappear during protonation "
        "state changes. RAD(i,:) contains the flags for residue type i. Atoms that disappear "
        "in the deprotonated state are flagged with 1.0 in the deprotonated flags and 0.0 in "
        "the protonated flags. Atoms always present are 0.0 in both flag sets. "
        "Organization follows the same state ordering as ch and ch_md based on res_type."
    ),
    default=None,
    notes="Controls scaling of van der Waals interactions for disappearing (dummy) hydrogens.",
    value_type="float_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["ch", "ch_md", "atom_name"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="model_pka",
    description=(
        "A two-dimensional array containing the experimentally determined model pKa values "
        "for each titratable residue type. MODEL_PKA(i,:) contains the pKa(s) for residue type i. "
        "For type 0/-2: a single pKa entry. "
        "For type 2: two pKa entries for the two tautomers. "
        "For type 4: a single pKa entry."
    ),
    default=None,
    notes=(
        "Example: MODEL_PKA(1,:) = 3.5 (AS2), MODEL_PKA(3,:) = 6.1,6.6 (HIP with two tautomer pKas). "
        "These are used to construct the pH-dependent free energy term in the extended Hamiltonian."
    ),
    value_type="float_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["parameters", "bar", "res_type"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="parameters",
    description=(
        "A two-dimensional array containing the model potential parameters for each titratable "
        "residue type, derived from thermodynamic integration (TI) simulations of model compounds. "
        "PARAMETERS(i,:) contains the fitted parameters for residue type i. "
        "For type 0/-2: two entries (A, B). "
        "For types 2 and 4: six entries (A, B, A0, B0, A10, B10). "
        "For type 4: additional entries 7-12 are R1-R6."
    ),
    default=None,
    notes=(
        "These parameters define the model compound PMF (potential of mean force) in lambda "
        "and x space, which is subtracted from the potential energy to yield an approximately "
        "flat PMF at the model pKa."
    ),
    value_type="float_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["model_pka", "bar", "res_type"],
    commonly_changed=True,
))

KEYWORDS.append(Keyword(
    name="bar",
    description=(
        "A two-dimensional array containing the heights of the quadratic barrier (penalty) "
        "potentials in the model potentials for each titratable residue type. These barriers "
        "ensure the system remains near the endpoints of lambda and x (i.e., near 0 or 1). "
        "BAR(i,:) contains the barrier height(s) for residue type i. "
        "For type 0/-2: one entry (barrier in lambda). "
        "For types 2 and 4: two entries (barrier in x, barrier in lambda)."
    ),
    default=None,
    notes="Example: BAR(1,:) = 2.5,2.5 for AS2 (type 4, barriers in x and lambda).",
    value_type="float_array_2d",
    section="phmdparm",
    category="Continuous CpHMD - Parameters",
    related=["parameters", "model_pka"],
    commonly_changed=False,
))

# ---------------------------------------------------------------------------
# &phmdstrt NAMELIST (Section 29.2.3 - Continuous CpHMD restart variables)
# ---------------------------------------------------------------------------

KEYWORDS.append(Keyword(
    name="ph_theta",
    description=(
        "The theta and theta_x values of the titration coordinates. These are the auxiliary "
        "variables from which lambda = sin^2(theta) is computed. Used to restart a continuous "
        "CpHMD simulation from a specific protonation state or to set initial conditions for "
        "TI parameterization runs."
    ),
    default=None,
    notes="Read from the -phmdstrt file. Written to the -phmdrstrt file at the end of a simulation.",
    value_type="float_array",
    section="phmdstrt",
    category="Continuous CpHMD - Restart",
    related=["vph_theta", "iphmd"],
    commonly_changed=False,
))

KEYWORDS.append(Keyword(
    name="vph_theta",
    description=(
        "The velocities of the titration coordinates (theta and theta_x). Used to restart "
        "a continuous CpHMD simulation with the same lambda dynamics state. If qphmdstart is "
        "true, these velocities are regenerated from the Boltzmann distribution instead."
    ),
    default=None,
    notes="Read from the -phmdstrt file when qphmdstart is false.",
    value_type="float_array",
    section="phmdstrt",
    category="Continuous CpHMD - Restart",
    related=["ph_theta", "qphmdstart"],
    commonly_changed=False,
))

# =============================================================================
# &wt NAMELIST (Section 22.9 - Varying Conditions)
# =============================================================================
# The &wt namelist has a fundamentally different structure from &cntrl, &ewald,
# etc. Rather than keyword=value pairs, it is a repeated namelist block:
#   &wt TYPE='...', ISTEP1=..., ISTEP2=..., VALUE1=..., VALUE2=..., /
# read repeatedly until TYPE='END'. The TYPE parameter selects what quantity
# is being varied, while the other parameters control the schedule.
# We model this with separate dataclasses.

@dataclass
class WtParameter:
    """A standard parameter in the &wt namelist block (shared across all TYPEs)."""
    name: str
    description: str
    default: Any
    value_type: str = "int"

@dataclass
class WtType:
    """A valid TYPE value for the &wt namelist varying-conditions block."""
    name: str
    description: str
    category: str = "Energy weighting"
    notes: Optional[str] = None
    special_params: Optional[str] = None  # If ISTEP/VALUE meanings differ from standard


# Standard &wt parameters (shared by all TYPE cards)
WT_PARAMETERS: List[WtParameter] = [
    WtParameter(
        name="TYPE",
        description=(
            "Defines quantity being varied. Must be uppercase. Valid options include energy "
            "term weights, temperature, cutoff, NMR controls, and more. "
            "Terminated by TYPE='END'."
        ),
        default=None,
        value_type="string",
    ),
    WtParameter(
        name="ISTEP1",
        description=(
            "This change is applied starting at step/iteration ISTEP1."
        ),
        default=0,
        value_type="int",
    ),
    WtParameter(
        name="ISTEP2",
        description=(
            "This change is applied through step/iteration ISTEP2. If ISTEP2=0, "
            "this change will remain in effect from step ISTEP1 to the end of the run "
            "at a value of VALUE1 (VALUE2 is ignored in this case)."
        ),
        default=0,
        value_type="int",
    ),
    WtParameter(
        name="VALUE1",
        description=(
            "Value of the change corresponding to ISTEP1."
        ),
        default=None,
        value_type="float",
    ),
    WtParameter(
        name="VALUE2",
        description=(
            "Value of the change corresponding to ISTEP2. If ISTEP2=0, VALUE2 is ignored."
        ),
        default=None,
        value_type="float",
    ),
    WtParameter(
        name="IINC",
        description=(
            "If IINC > 0, then the change is applied as a step function, with IINC steps/iterations "
            "between each change in the target VALUE (ignored if ISTEP2=0). If IINC=0, the change "
            "is done continuously."
        ),
        default=0,
        value_type="int",
    ),
    WtParameter(
        name="IMULT",
        description=(
            "If IMULT=0, then the change will be linearly interpolated from VALUE1 to VALUE2 "
            "as the step number increases from ISTEP1 to ISTEP2 (default). If IMULT=1, then the "
            "change will be effected by a series of multiplicative scalings, using a single factor R "
            "for all scalings: VALUE2 = (R**INCREMENTS) * VALUE1."
        ),
        default=0,
        value_type="int",
    ),
]

# Valid TYPE values
WT_TYPES: List[WtType] = [
    # --- Energy term weighting ---
    WtType(
        name="BOND",
        description="Varies the relative weighting of bond energy terms.",
        category="Energy weighting",
    ),
    WtType(
        name="ANGLE",
        description="Varies the relative weighting of valence angle energy terms.",
        category="Energy weighting",
    ),
    WtType(
        name="TORSION",
        description=(
            "Varies the relative weighting of torsion (and J-coupling) energy terms. "
            "Note that any restraints defined in the input to the PARM program are included. "
            "Improper torsions are handled separately (IMPROP)."
        ),
        category="Energy weighting",
    ),
    WtType(
        name="IMPROP",
        description=(
            "Varies the relative weighting of the 'improper' torsional terms. "
            "These are not included in TORSION."
        ),
        category="Energy weighting",
    ),
    WtType(
        name="VDW",
        description=(
            "Varies the relative weighting of van der Waals energy terms. "
            "This is equivalent to changing the well depth (epsilon) by the given factor."
        ),
        category="Energy weighting",
    ),
    WtType(
        name="HB",
        description="Varies the relative weighting of hydrogen-bonding energy terms.",
        category="Energy weighting",
    ),
    WtType(
        name="ELEC",
        description="Varies the relative weighting of electrostatic energy terms.",
        category="Energy weighting",
    ),
    WtType(
        name="NB",
        description="Varies the relative weights of the non-bonded (VDW, HB, and ELEC) terms.",
        category="Energy weighting",
    ),
    WtType(
        name="ATTRACT",
        description="Varies the relative weights of the attractive parts of the van der Waals and H-bond terms.",
        category="Energy weighting",
    ),
    WtType(
        name="REPULSE",
        description="Varies the relative weights of the repulsive parts of the van der Waals and H-bond terms.",
        category="Energy weighting",
    ),
    WtType(
        name="RSTAR",
        description=(
            "Varies the effective van der Waals radii for the VDW interactions by the given factor. "
            "Note that this is done by changing the relative attractive and repulsive coefficients, so "
            "ATTRACT/REPULSE should not be used over the same step range as RSTAR."
        ),
        category="Energy weighting",
        notes=(
            "Changes result in exponential weighting changes to the attractive and repulsive terms "
            "(proportional to the scale factor**6 and **12, respectively). Scaling RSTAR to a very "
            "small value (e.g. ~0.1) may result in a zeroing-out of the vdw term."
        ),
    ),
    WtType(
        name="INTERN",
        description=(
            "Varies the relative weights of the BOND, ANGLE and TORSION terms. "
            "'Improper' torsions (IMPROP) must be varied separately."
        ),
        category="Energy weighting",
    ),
    WtType(
        name="ALL",
        description=(
            "Varies the relative weights of all the energy terms above (BOND, ANGLE, TORSION, VDW, "
            "HB, and ELEC; does not affect RSTAR or IMPROP)."
        ),
        category="Energy weighting",
    ),
    # --- NMR restraint weighting ---
    WtType(
        name="REST",
        description="Varies the relative weights of *all* the NMR restraint energy terms.",
        category="NMR restraints",
    ),
    WtType(
        name="RESTS",
        description=(
            "Varies the weights of the 'short-range' NMR restraints. Short-range restraints are "
            "defined by the SHORT instruction (see below)."
        ),
        category="NMR restraints",
    ),
    WtType(
        name="RESTL",
        description=(
            "Varies the weights of any NMR restraints which are not defined as 'short range' "
            "by the SHORT instruction. When no SHORT instruction is given, RESTL is equivalent to REST."
        ),
        category="NMR restraints",
    ),
    WtType(
        name="NOESY",
        description=(
            "Varies the overall weight for NOESY volume restraints. Note that this value multiplies "
            "the individual weights read into the 'awt' array. Only if NMROPT=2."
        ),
        category="NMR restraints",
    ),
    WtType(
        name="SHIFTS",
        description=(
            "Varies the overall weight for chemical shift restraints. Note that this value multiplies "
            "the individual weights read into the 'wt' array. Only if NMROPT=2."
        ),
        category="NMR restraints",
    ),
    WtType(
        name="SHORT",
        description=(
            "Defines the short-range restraints. For this instruction, ISTEP1, ISTEP2, VALUE1, "
            "and VALUE2 have different meanings."
        ),
        category="NMR restraints",
        special_params=(
            "Method 1 (sequence proximity): ISTEP1 <= ABS(delta_residue) <= ISTEP2. "
            "Method 2 (distance): VALUE1 <= distance <= VALUE2. "
            "Only one SHORT command can be issued, and values remain fixed throughout the run. "
            "If IINC>0, the short-range interaction list will be re-evaluated every IINC steps."
        ),
    ),
    # --- Simulation conditions ---
    WtType(
        name="TGTRMSD",
        description="Varies the RMSD target value for targeted MD.",
        category="Simulation conditions",
    ),
    WtType(
        name="TEMP0",
        description="Varies the target temperature TEMP0.",
        category="Simulation conditions",
    ),
    WtType(
        name="TEMP0LES",
        description="Varies the LES target temperature TEMP0LES.",
        category="Simulation conditions",
    ),
    WtType(
        name="TAUTP",
        description=(
            "Varies the coupling parameter, TAUTP, used in temperature scaling when "
            "temperature coupling option NTT=1 is used."
        ),
        category="Simulation conditions",
    ),
    WtType(
        name="CUT",
        description="Varies the non-bonded cutoff distance.",
        category="Simulation conditions",
    ),
    # --- Step counter controls ---
    WtType(
        name="NSTEP0",
        description=(
            "If present, this instruction will reset the initial value of the step counter "
            "(against which ISTEP1/ISTEP2 and NSTEP1/NSTEP2 are compared) to the value ISTEP1. "
            "This only affects the way in which NMR weight restraints are calculated. It does not "
            "affect the value of NSTEP printed as part of the dynamics output."
        ),
        category="Step control",
        special_params=(
            "Only ISTEP1 is used; ISTEP2, VALUE1, VALUE2 and IINC are ignored. "
            "An NSTEP0 instruction only has an effect at the beginning of a run. "
            "Useful for simulation restarts, where NSTEP0 is set to the final step of the previous run."
        ),
    ),
    WtType(
        name="STPMLT",
        description=(
            "If present, the NMR step counter will be changed in increments of STPMLT for "
            "each actual dynamics step."
        ),
        category="Step control",
        special_params=(
            "Only VALUE1 is read. ISTEP1, ISTEP2, VALUE2, IINC, and IMULT are ignored. "
            "Default = 1.0."
        ),
    ),
    # --- Time averaging ---
    WtType(
        name="DISAVE",
        description=(
            "If present, then by default time-averaged values (rather than instantaneous values) "
            "for distance restraints will be used."
        ),
        category="Time averaging",
        special_params=(
            "VALUE1 = tau (characteristic time for exponential decay). "
            "VALUE2 = POWER (power used in averaging; nearest integer is used). "
            "The range (ISTEP1-ISTEP2) applies only to TAU; POWER is not changed by "
            "subsequent cards. Default tau (if 0.0) is 1e6, resulting in no exponential decay. "
            "Any tau >= 1e6 results in no exponential decay."
        ),
    ),
    WtType(
        name="ANGAVE",
        description=(
            "If present, then by default time-averaged values (rather than instantaneous values) "
            "for angle restraints will be used."
        ),
        category="Time averaging",
        special_params="Same parameter meanings as DISAVE but for angle data.",
    ),
    WtType(
        name="TORAVE",
        description=(
            "If present, then by default time-averaged values (rather than instantaneous values) "
            "for torsion restraints will be used."
        ),
        category="Time averaging",
        special_params="Same parameter meanings as DISAVE but for torsion data.",
    ),
    WtType(
        name="DISAVI",
        description=(
            "Controls initial values and dump frequency for time-averaged distance restraints."
        ),
        category="Time averaging",
        special_params=(
            "ISTEP1: Ignored. "
            "ISTEP2: Sets IDMPAV; if > 0 and DUMPAVE file specified, time-averaged values "
            "written every IDMPAV steps. "
            "VALUE1: If != 0, resets initial value of internal r. "
            "  -1000 < VALUE1 < 1000: Initial = r_initial + VALUE1. "
            "  VALUE1 <= -1000: Initial = r_target + 1000. "
            "  VALUE1 >= 1000: Initial = r_target - 1000. "
            "VALUE2: If > 0, sets tau for calculating final reported averages. "
            "IINC: If 0, exact force formula. If 1, approximate (non-conservative) forces."
        ),
        notes="Has no effect unless the corresponding DISAVE card is also present.",
    ),
    WtType(
        name="ANGAVI",
        description="Controls initial values and dump frequency for time-averaged angle restraints.",
        category="Time averaging",
        special_params="Same parameter meanings as DISAVI but for angle data.",
        notes="Has no effect unless the corresponding ANGAVE card is also present.",
    ),
    WtType(
        name="TORAVI",
        description="Controls initial values and dump frequency for time-averaged torsion restraints.",
        category="Time averaging",
        special_params="Same parameter meanings as DISAVI but for torsion data.",
        notes="Has no effect unless the corresponding TORAVE card is also present.",
    ),
    WtType(
        name="DUMPFREQ",
        description=(
            "ISTEP1 is the only parameter read, and it sets the frequency at which the coordinates "
            "in the distance or angle restraints are dumped to the file specified by the DUMPAVE "
            "command in the I/O redirection section."
        ),
        category="Time averaging",
        special_params="Only ISTEP1 is read. ISTEP2 and IMULT are ignored.",
    ),
    WtType(
        name="END",
        description="END of the &wt section. This terminates the reading of &wt namelists.",
        category="Control",
    ),
]

# General notes for the &wt section
WT_NOTES = [
    "All weights are relative to a default of 1.0 in the standard force field.",
    "Weights are not cumulative.",
    (
        "For any range where the weight of a term is not modified, the weight reverts to 1.0. "
        "For any range where TEMP0 or CUTOFF is not specified, the value is set to that in the input file."
    ),
    (
        "If a weight is set to 0.0, it is set internally to 1e-7. This can be overridden by setting "
        "the weight to a negative number, in which case exactly 0.0 will be used. However, if any "
        "weight is set to exactly 0.0, it cannot be changed again during the run."
    ),
    (
        "If two or more cards change a particular weight over the same range, the weight given "
        "on the last applicable card will be the one used."
    ),
    (
        "Once any weight change for which NSTEP2=0 becomes active (i.e. one effective for the "
        "remainder of the run), the weight of this term cannot be further modified."
    ),
]


# =============================================================================
# FILE REDIRECTION COMMANDS (Section 22.10)
# =============================================================================
# These follow the &wt section when NMROPT > 0, with format: TYPE = filename

@dataclass
class FileRedirect:
    """A file redirection command for NMR-related I/O (Section 22.10)."""
    name: str
    description: str

FILE_REDIRECTIONS: List[FileRedirect] = [
    FileRedirect(
        name="LISTIN",
        description=(
            "An output listing of the restraints which have been read, and their deviations "
            "from the target distances before the simulation has been run. By default, this "
            "listing is not printed. If POUT is used for the filename, these deviations will "
            "be printed in the normal output file."
        ),
    ),
    FileRedirect(
        name="LISTOUT",
        description=(
            "An output listing of the restraints which have been read, and their deviations "
            "from the target distances after the simulation has finished. By default, this "
            "listing is not printed. If POUT is used for the filename, these deviations will "
            "be printed in the normal output file."
        ),
    ),
    FileRedirect(
        name="DISANG",
        description=(
            "The file from which the distance and angle restraint information (Section 30.1) "
            "will be read."
        ),
    ),
    FileRedirect(
        name="NOESY",
        description="File from which NOESY volume information (Section 30.2) will be read.",
    ),
    FileRedirect(
        name="SHIFTS",
        description="File from which chemical shift information (Section 30.3) will be read.",
    ),
    FileRedirect(
        name="PCSHIFT",
        description="File from which paramagnetic shift information (Section 30.3) will be read.",
    ),
    FileRedirect(
        name="DIPOLE",
        description="File from which residual dipolar couplings (Section 30.5) will be read.",
    ),
    FileRedirect(
        name="CSA",
        description="File from which CSA or pseudo-CSA restraints (Section 30.6) will be read.",
    ),
    FileRedirect(
        name="DUMPAVE",
        description=(
            "File to which the time-averaged values of all restraints will be written. "
            "If DISAVI/ANGAVI/TORAVI has been used to set IDMPAV != 0, then averaged values "
            "will be output. If the DUMPFREQ command has been used, the instantaneous values "
            "will be output."
        ),
    ),

    FileRedirect(
        name="cpin",
        description=(
            "Input protonation state definitions for constant pH MD (Chapter 27). "
            "Generated by cpinutil.py, this file describes which residues titrate, their "
            "possible protonation states, partial charge vectors for each state, and relative "
            "reference energies. When restarting a simulation, use the cprestrt file from the "
            "previous run as the cpin file to preserve final protonation states."
        ),
    ),
    FileRedirect(
        name="cpout",
        description=(
            "Output protonation state history for constant pH MD (Chapter 27). Records are "
            "written at each Monte Carlo step. Full records (listing all residues with solvent pH, "
            "MC step size, and time) are written on the first step and every ntwx steps. Delta "
            "records list only residues that were examined. Analyze with cphstats to compute "
            "pKa values and protonation state populations."
        ),
    ),
    FileRedirect(
        name="cprestrt",
        description=(
            "Protonation state restart file for constant pH MD (Chapter 27). Written in cpin "
            "format, contains the final protonation states from the simulation. Should be used "
            "as the cpin file when restarting a constant pH simulation to maintain protonation "
            "state continuity."
        ),
    ),
    FileRedirect(
        name="cein",
        description=(
            "Input redox state definitions for constant redox potential (Eh) MD. Analogous to "
            "cpin for constant pH, this file describes titratable redox groups, their possible "
            "oxidation states, and reference energies. Generated by ceinutil.py."
        ),
    ),
    FileRedirect(
        name="ceout",
        description=(
            "Output redox state data saved over trajectory for constant redox potential MD. "
            "Analogous to cpout for constant pH simulations."
        ),
    ),
    FileRedirect(
        name="cerestrt",
        description=(
            "Redox state restart file for constant redox potential MD. Contains final redox "
            "states in cein format. Should be used as the cein file when restarting a constant "
            "Eh simulation."
        ),
    ),

    FileRedirect(
        name="phmdin",
        description=(
            "Input file for continuous constant pH MD (Chapter 29) containing the &phmdin "
            "namelist. Specifies titration parameters including lambda particle mass, temperature, "
            "friction, titratable residue masks, and output settings. Required when iphmd=1 or 3."
        ),
    ),
    FileRedirect(
        name="phmdparm",
        description=(
            "Parameter file for continuous constant pH MD (Chapter 29) containing the &phmdparm "
            "namelist. Defines titratable residue types, atom names, charge vectors for each "
            "protonation state, disappearing atom flags, model pKa values, model potential "
            "parameters (from TI), and barrier heights. Required when iphmd=1 or 3."
        ),
    ),
    FileRedirect(
        name="phmdstrt",
        description=(
            "Optional restart input file for continuous constant pH MD (Chapter 29) containing "
            "the &phmdstrt namelist. Provides theta values and velocities of the titration "
            "coordinates to restart from a previous run or to set initial conditions for "
            "thermodynamic integration simulations."
        ),
    ),
    FileRedirect(
        name="phmdrstrt",
        description=(
            "Restart output file for continuous constant pH MD (Chapter 29) containing the "
            "&phmdrst namelist. Written at the end of a simulation with the final theta values "
            "and velocities of the titration coordinates. Should be used as -phmdstrt input "
            "when restarting a continuous CpHMD simulation."
        ),
    ),
    FileRedirect(
        name="phmdout",
        description=(
            "Lambda output file for continuous constant pH MD (Chapter 29). Contains the "
            "lambda values of each titratable residue over the course of the simulation, "
            "written every nprint_phmd steps when prlam is true. Used for post-simulation "
            "analysis of protonation states and pKa computation."
        ),
    ),
]

def get_keyword(name: str) -> Keyword:
    """Look up a keyword by name (case-insensitive)."""
    name_lower = name.lower()
    for kw in KEYWORDS:
        if kw.name.lower() == name_lower:
            return kw
    raise KeyError(f"Keyword '{name}' not found in database.")


def get_by_category(category: str) -> List[Keyword]:
    """Return all keywords matching a category substring (case-insensitive)."""
    cat_lower = category.lower()
    return [kw for kw in KEYWORDS if cat_lower in kw.category.lower()]


def get_by_section(section: str) -> List[Keyword]:
    """Return all keywords in a given namelist section."""
    return [kw for kw in KEYWORDS if kw.section.lower() == section.lower()]


def get_commonly_changed() -> List[Keyword]:
    """Return all keywords marked as commonly changed (boldface in manual)."""
    return [kw for kw in KEYWORDS if kw.commonly_changed]


def get_wt_type(name: str) -> WtType:
    """Look up a &wt TYPE value by name (case-insensitive)."""
    name_upper = name.upper()
    for wt in WT_TYPES:
        if wt.name == name_upper:
            return wt
    raise KeyError(f"&wt TYPE '{name}' not found in database.")


def get_wt_types_by_category(category: str) -> List[WtType]:
    """Return all &wt TYPE values matching a category substring (case-insensitive)."""
    cat_lower = category.lower()
    return [wt for wt in WT_TYPES if cat_lower in wt.category.lower()]


def get_file_redirect(name: str) -> FileRedirect:
    """Look up a file redirection command by name (case-insensitive)."""
    name_upper = name.upper()
    for fr in FILE_REDIRECTIONS:
        if fr.name == name_upper:
            return fr
    raise KeyError(f"File redirect '{name}' not found in database.")


def summary_stats():
    """Print summary statistics about the database."""
    sections = {}
    categories = {}
    for kw in KEYWORDS:
        sections[kw.section] = sections.get(kw.section, 0) + 1
        categories[kw.category] = categories.get(kw.category, 0) + 1

    print(f"Total keywords: {len(KEYWORDS)}")
    print(f"\nBy namelist section:")
    for sec, count in sorted(sections.items()):
        print(f"  &{sec}: {count}")
    print(f"\nBy category:")
    for cat, count in sorted(categories.items()):
        print(f"  {cat}: {count}")
    print(f"\nCommonly changed (boldface in manual): {len(get_commonly_changed())}")

    # &wt summary
    wt_cats = {}
    for wt in WT_TYPES:
        wt_cats[wt.category] = wt_cats.get(wt.category, 0) + 1
    print(f"\n&wt TYPE values: {len(WT_TYPES)}")
    for cat, count in sorted(wt_cats.items()):
        print(f"  {cat}: {count}")

    print(f"\n&wt parameters: {len(WT_PARAMETERS)}")
    print(f"File redirection commands: {len(FILE_REDIRECTIONS)}")


if __name__ == "__main__":
    summary_stats()
    print("\n" + "="*60)
    print("Sample: Commonly changed keywords")
    print("="*60)
    for kw in get_commonly_changed():
        default_str = kw.default if kw.default is not None else "(auto)"
        print(f"  {kw.name:<22s}  default={str(default_str):<12s}  [{kw.category}]")
    print("\n" + "="*60)
    print("&wt TYPE values")
    print("="*60)
    for wt in WT_TYPES:
        special = " *" if wt.special_params else ""
        print(f"  {wt.name:<12s}  [{wt.category}]{special}")
    print("\n" + "="*60)
    print("File redirection commands (Section 22.10)")
    print("="*60)
    for fr in FILE_REDIRECTIONS:
        print(f"  {fr.name:<12s}  {fr.description[:70]}...")
