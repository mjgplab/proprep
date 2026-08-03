"""
AMBER Parameter Educational Database
=====================================

Pure data module containing structured educational metadata for all AMBER mdin
parameters exposed by the wizard. No UI, no Rich imports, no prompts.

Each parameter includes:
- brief: One-line summary
- what_it_does: Physical/computational meaning
- options: Valid choices with descriptions (where applicable)
- manual_notes: Guidance the Amber manual actually states (often None)
- common: Whether this is commonly or rarely changed
- related: Cross-references to related parameters

Provenance
----------
Descriptions here are written *about* the parameters; they are not quotations
from the Amber manual. Defaults, option sets and namelist assignments are
checked against it (chapters 23 sander / 24 pmemd) -- see
docs/mdin_help_audit.md for the parameter-by-parameter record.

This module previously carried ``why_default`` and ``when_to_change`` fields.
Both were removed: the manual rarely explains why a default was chosen or when
to deviate, so the fields could only be filled with unsourceable advice, and an
audit found invented specifics sitting indistinguishably beside sourced ones
(gamma_ln's "~50 ps^-1 for water" is the manual's; therm_par's "2-50 ps^-1" was
not). Anything the manual does say belongs in ``manual_notes``.

Do not reintroduce advisory fields. If a claim cannot be pointed at a line in
the manual, it does not belong here -- no help is better than wrong help.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# =============================================================================
# Data Structures
# =============================================================================
@dataclass
class Parameter:
    """Educational metadata for a single AMBER mdin parameter."""
    name: str
    brief: str
    what_it_does: str
    default: Any
    param_type: str  # "int", "float", "str"
    options: Optional[Dict[Any, str]] = None
    # Guidance the Amber manual actually states, quoted or closely paraphrased.
    # Leave as None when the manual says nothing -- an empty note is the honest
    # outcome for most parameters, and inventing advice to fill it is the exact
    # failure this field replaced. See docs/mdin_help_audit.md.
    manual_notes: Optional[str] = None
    min_val: Optional[float] = None
    max_val: Optional[float] = None
    related: List[str] = field(default_factory=list)
    common: bool = True
    namelist: str = "cntrl"
    category: str = "general"


# =============================================================================
# Complete Parameter Database
# =============================================================================
def build_parameter_database() -> Dict[str, Parameter]:
    """Build the complete educational parameter database.

    Returns:
        Dictionary mapping parameter name to Parameter instance.
    """
    db: Dict[str, Parameter] = {}

    # =========================================================================
    # SIMULATION CONTROL
    # =========================================================================
    db["imin"] = Parameter(
        name="imin",
        brief="Minimization vs. molecular dynamics",
        what_it_does=(
            "Controls whether to run energy minimization or molecular dynamics. "
            "Minimization finds a local energy minimum by moving atoms 'downhill' on the "
            "potential energy surface -- no kinetic energy, no time evolution. "
            "MD simulates actual time evolution with Newton's equations."
        ),
        default=0,
        param_type="int",
        options={
            0: "Molecular dynamics (time evolution with velocities)",
            1: "Energy minimization (find local energy minimum)",
            5: "Trajectory post-processing with minimization",
            6: "Trajectory post-processing with MD driver (single points)",
            7: "Socket server mode (for external drivers like i-PI)",
        },
        related=["maxcyc", "ncyc", "ntmin", "nstlim", "dt"],
        common=True,
        category="Simulation Control",
    )

    db["nmropt"] = Parameter(
        name="nmropt",
        brief="NMR restraints and varying conditions",
        what_it_does=(
            "Enables reading of NMR-style restraints and/or varying simulation conditions. "
            "Despite the name, this is useful even without NMR data: distance/angle/torsion "
            "restraints, simulated annealing (changing temperature during simulation), "
            "and umbrella sampling."
        ),
        default=0,
        param_type="int",
        options={
            0: "No restraints or varying conditions",
            1: "Read restraints and/or weight changes (&wt namelist)",
            2: "Full NMR: restraints + NOESY + chemical shifts + dipolar couplings",
        },
        related=["temp0"],
        common=False,
        category="Simulation Control",
    )

    # =========================================================================
    # INPUT / RESTART
    # =========================================================================
    db["ntx"] = Parameter(
        name="ntx",
        brief="What to read from input coordinate file",
        what_it_does=(
            "Controls what information is read from the inpcrd/restart file. "
            "ntx=1 reads coordinates only (for new simulations). "
            "ntx=5 reads coordinates and velocities (for continuing a simulation). "
            "Velocities encode the momentum of each atom; reading them allows seamless continuation."
        ),
        default=1,
        param_type="int",
        options={
            1: "Read coordinates only (no velocities)",
            5: "Read coordinates AND velocities (for restarts)",
        },
        related=["irest"],
        common=True,
        category="Input/Restart",
    )

    db["irest"] = Parameter(
        name="irest",
        brief="Restart vs. new simulation",
        what_it_does=(
            "Determines whether this is a continuation or a fresh start. "
            "irest=0: new simulation -- time starts at 0, velocities ignored or regenerated. "
            "irest=1: restart -- time continues from restart file, velocities used. "
            "Common patterns: new simulation (ntx=1, irest=0), restart (ntx=5, irest=1)."
        ),
        default=0,
        param_type="int",
        options={
            0: "New simulation (fresh start, time = 0)",
            1: "Restart (continue from restart file)",
        },
        related=["ntx", "ig"],
        common=True,
        category="Input/Restart",
    )

    # =========================================================================
    # MINIMIZATION
    # =========================================================================
    db["maxcyc"] = Parameter(
        name="maxcyc",
        brief="Maximum minimization cycles",
        what_it_does=(
            "The maximum number of steps the minimizer will take. "
            "Each cycle moves atoms to reduce the potential energy. "
            "The minimizer may stop early if convergence (drms) is reached."
        ),
        default=1,
        param_type="int",
        min_val=1,
        related=["ncyc", "drms", "ntmin"],
        common=True,
        category="Minimization",
    )

    db["ncyc"] = Parameter(
        name="ncyc",
        brief="Steepest descent cycles before conjugate gradient",
        what_it_does=(
            "Number of steepest descent steps before switching to conjugate gradient. "
            "Steepest descent is robust but slow -- good for bad starting structures with clashes. "
            "Conjugate gradient converges much faster but can fail on very bad structures."
        ),
        default=10,
        param_type="int",
        min_val=0,
        manual_notes=(
            "Applies only when ntmin=1: \"If NTMIN is 1 then the method of minimization will be switched from steepest descent to conjugate gradient after NCYC cycles.\""
        ),
        related=["maxcyc", "ntmin"],
        common=True,
        category="Minimization",
    )

    db["ntmin"] = Parameter(
        name="ntmin",
        brief="Minimization algorithm",
        what_it_does=(
            "Selects the minimization method. "
            "Standard choice is ntmin=1 (steepest descent then conjugate gradient)."
        ),
        default=1,
        param_type="int",
        options={
            0: "Conjugate gradient (with brief SD after pairlist updates)",
            1: "Steepest descent for ncyc, then conjugate gradient",
            2: "Steepest descent only (most robust, slowest)",
            3: "XMIN method",
            4: "LMOD method",
            5: "Truncated Newton conjugate gradient",
        },
        related=["ncyc", "maxcyc"],
        common=False,
        category="Minimization",
    )

    db["drms"] = Parameter(
        name="drms",
        brief="RMS gradient convergence criterion",
        what_it_does=(
            "Minimization stops when the root-mean-square of the energy gradient "
            "(force magnitude) falls below this value. Units: kcal/mol/A."
        ),
        default=1.0e-4,
        param_type="float",
        min_val=0.0,
        related=["maxcyc"],
        common=True,
        category="Minimization",
    )

    db["dx0"] = Parameter(
        name="dx0",
        brief="Initial minimization step size",
        what_it_does="Initial step length for the minimizer. Usually auto-adjusted.",
        default=0.01,
        param_type="float",
        min_val=0.0,
        common=False,
        category="Minimization",
    )

    # =========================================================================
    # MOLECULAR DYNAMICS
    # =========================================================================
    db["nstlim"] = Parameter(
        name="nstlim",
        brief="Number of MD steps",
        what_it_does=(
            "Total number of MD steps to perform. "
            "Total simulation time = nstlim x dt. "
            "Example: 500,000 steps x 0.002 ps = 1000 ps = 1 ns."
        ),
        default=1,
        param_type="int",
        min_val=0,
        related=["dt"],
        common=True,
        category="Molecular Dynamics",
    )

    db["dt"] = Parameter(
        name="dt",
        brief="Time step (picoseconds)",
        what_it_does=(
            "Time between force evaluations. Must be small enough to accurately "
            "integrate the fastest motions. Without SHAKE: ~1 fs (0.001 ps). "
            "With SHAKE on H-bonds: 2 fs (0.002 ps). "
            "With SHAKE + hydrogen mass repartitioning: 4 fs (0.004 ps)."
        ),
        default=0.001,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "\"Recommended MAXIMUM is .002 if SHAKE is used, or .001 if it isn't.\" Above 300 K the step size should be reduced. Hydrogen mass repartitioning with SHAKE allows up to .004."
        ),
        related=["ntc", "nstlim"],
        common=True,
        category="Molecular Dynamics",
    )

    db["t"] = Parameter(
        name="t",
        brief="Initial time (ps)",
        what_it_does="Starting time for the simulation clock. Mainly for bookkeeping.",
        default=0.0,
        param_type="float",
        manual_notes=(
            "The start time \"is for your own reference and is not critical\". It is taken from the coordinate input file when irest=1."
        ),
        common=False,
        category="Molecular Dynamics",
    )

    db["nrespa"] = Parameter(
        name="nrespa",
        brief="Multiple time stepping factor",
        what_it_does=(
            "For PME, evaluates the reciprocal-space (slow) forces every nrespa steps "
            "instead of every step. Can speed up CPU simulations modestly."
        ),
        default=1,
        param_type="int",
        min_val=1,
        manual_notes=(
            "Energies \"are only accessible every nrespa steps, since the values at other times are meaningless.\""
        ),
        related=["dt"],
        common=False,
        category="Molecular Dynamics",
    )

    db["nscm"] = Parameter(
        name="nscm",
        brief="Center-of-mass motion removal frequency",
        what_it_does=(
            "How often to remove center-of-mass translation (and rotation for non-periodic). "
            "Prevents the 'flying ice cube' effect where kinetic energy accumulates "
            "in overall translation/rotation instead of internal motion."
        ),
        default=1000,
        param_type="int",
        min_val=0,
        common=False,
        category="Molecular Dynamics",
    )

    # =========================================================================
    # TEMPERATURE CONTROL
    # =========================================================================
    db["ntt"] = Parameter(
        name="ntt",
        brief="Thermostat type",
        what_it_does=(
            "Selects the method for temperature control. "
            "NVE (ntt=0): true Newtonian dynamics, total energy conserved. "
            "Langevin (ntt=3): adds friction + random forces, mimics solvent collisions, "
            "very stable and widely used but affects dynamics. "
            "Bussi (ntt=11): stochastic rescaling, correct ensemble, less friction than Langevin. "
            "Berendsen (ntt=1): scales velocities toward target, produces wrong fluctuations -- "
            "NOT recommended for production."
        ),
        default=0,
        param_type="int",
        options={
            0: "NVE -- constant energy (no thermostat)",
            1: "Berendsen weak coupling (NOT for production)",
            2: "Andersen (velocity randomization)",
            3: "Langevin dynamics (recommended for most simulations)",
            9: "Optimized isokinetic Nose-Hoover (RESPA compatible)",
            10: "Stochastic isokinetic Nose-Hoover",
            11: "Bussi thermostat (good alternative to Langevin)",
        },
        manual_notes=(
            "The manual warns ntt=1 is \"especially dangerous for generalized Born simulations, where there are no collisions with solvent to aid in thermalization. Other temperature coupling options (especially ntt=3) should be used instead.\""
        ),
        related=["temp0", "gamma_ln", "tautp", "vrand"],
        common=True,
        category="Temperature Control",
    )

    db["temp0"] = Parameter(
        name="temp0",
        brief="Target temperature (Kelvin)",
        what_it_does=(
            "The temperature the thermostat tries to maintain. "
            "For NVT/NPT simulations, this is your simulation temperature."
        ),
        default=300.0,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "The manual notes that for temperatures above 300 K the step size should be reduced, since increased distance travelled between evaluations \"can lead to SHAKE and other problems\"."
        ),
        related=["ntt", "tempi"],
        common=True,
        category="Temperature Control",
    )

    db["tempi"] = Parameter(
        name="tempi",
        brief="Initial temperature for velocity assignment (K)",
        what_it_does=(
            "When starting a new simulation (ntx=1), velocities are assigned from a "
            "Maxwell-Boltzmann distribution at this temperature. "
            "If tempi=0, velocities are calculated from forces instead (cold start). "
            "IGNORED when ntx=5 (velocities read from restart file)."
        ),
        default=0.0,
        param_type="float",
        min_val=0.0,
        related=["ntx", "temp0"],
        common=True,
        category="Temperature Control",
    )

    db["ig"] = Parameter(
        name="ig",
        brief="Random number seed",
        what_it_does=(
            "Seed for the random number generator used for initial velocity assignment, "
            "Langevin dynamics random forces (ntt=3), and Andersen thermostat collisions (ntt=2). "
            "ig=-1 uses the current time to generate a unique seed each run."
        ),
        default=-1,
        param_type="int",
        related=["ntt"],
        common=True,
        category="Temperature Control",
    )

    db["gamma_ln"] = Parameter(
        name="gamma_ln",
        brief="Langevin collision frequency (ps^-1)",
        what_it_does=(
            "Friction coefficient for Langevin dynamics (ntt=3). "
            "Higher = stronger coupling to heat bath, more damped dynamics. "
            "Physical water collision rate is ~50 ps^-1, but lower values work "
            "better for sampling in simulations."
        ),
        default=0.0,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "The physical collision frequency for liquid water is about 50 ps^-1, but the manual notes it \"is often advantageous, in terms of sampling or stability of integration, to use much smaller values, around 2 to 5 ps^-1\"."
        ),
        related=["ntt"],
        common=True,
        category="Temperature Control",
    )

    db["tautp"] = Parameter(
        name="tautp",
        brief="Berendsen/Bussi coupling time constant (ps)",
        what_it_does=(
            "For Berendsen (ntt=1) or Bussi (ntt=11): how quickly temperature adjusts "
            "toward temp0. Smaller = faster adjustment, less natural dynamics."
        ),
        default=1.0,
        param_type="float",
        min_val=0.0,
        related=["ntt"],
        common=False,
        category="Temperature Control",
    )

    db["vrand"] = Parameter(
        name="vrand",
        brief="Andersen thermostat collision interval (steps)",
        what_it_does="For ntt=2: steps between velocity randomization.",
        default=1000,
        param_type="int",
        min_val=1,
        related=["ntt"],
        common=False,
        category="Temperature Control",
    )

    db["vlimit"] = Parameter(
        name="vlimit",
        brief="Maximum velocity magnitude (A/ps)",
        what_it_does=(
            "Safety valve: any velocity component exceeding vlimit is reduced to vlimit. "
            "Warning messages indicate problems with the simulation."
        ),
        default=20.0,
        param_type="float",
        manual_notes=(
            "VLIMIT 20 is \"well above the most probable velocity in a Maxwell-Boltzmann distribution at room temperature\". GPU executables such as pmemd.cuda default to -1 and do NOT support vlimit. \"Runs that have more than a few such warnings should be carefully examined.\""
        ),
        common=False,
        category="Temperature Control",
    )

    db["nkija"] = Parameter(
        name="nkija",
        brief="Isokinetic integrator substeps / auxiliary variables (ntt=9,10)",
        what_it_does=(
            "Means different things per integrator. For ntt=9 it is the number of substeps "
            "of dt used when integrating the thermostat equations of motion. For ntt=10 it "
            "is the number of additional auxiliary velocity variables."
        ),
        default=1,
        param_type="int",
        manual_notes=(
            "\"For ntt=9, this the number of substeps of dt when integrating the thermostat "
            "equations of motion, for greater accuracy. For ntt=10, this specifies the number "
            "of additional auxiliary velocity variables v1 and v2.\" Default is 1 for both."
        ),
        min_val=1,
        related=["ntt"],
        common=False,
        category="Temperature Control",
    )

    db["idistr"] = Parameter(
        name="idistr",
        brief="Thermostat distribution accumulation frequency (ntt=9)",
        what_it_does=(
            "For the optimized isokinetic Nose-Hoover integrator (ntt=9), how often the "
            "thermostat velocity distribution functions are accumulated. A step count, "
            "not a mode selector."
        ),
        default=0,
        param_type="int",
        manual_notes=(
            "\"For the isokinetic integrator (ntt=9), the frequency at which the "
            "thermostat velocity distribution functions are accumulated.\""
        ),
        related=["ntt", "nkija"],
        common=False,
        category="Temperature Control",
    )

    db["sinrtau"] = Parameter(
        name="sinrtau",
        brief="Stochastic time constant (ps) for ntt=10",
        what_it_does="Stochastic coupling time constant for the stochastic isokinetic Nose-Hoover thermostat.",
        default=1.0,
        param_type="float",
        min_val=0.0,
        related=["ntt", "nkija"],
        common=False,
        category="Temperature Control",
    )

    db["temp0les"] = Parameter(
        name="temp0les",
        brief="LES particle temperature (K)",
        what_it_does=(
            "Target temperature for locally enhanced sampling (LES) particles. "
            "-1 means use the same temperature as regular atoms (temp0)."
        ),
        default=-1.0,
        param_type="float",
        related=["temp0"],
        common=False,
        category="Temperature Control",
    )

    # =========================================================================
    # INTEGRATOR (MIDDLE SCHEME)
    # =========================================================================
    db["ischeme"] = Parameter(
        name="ischeme",
        brief="Integration scheme",
        what_it_does=(
            "Selects the integration scheme. The 'middle' scheme (ischeme=1) is based on "
            "the leapfrog algorithm and places the thermostat step between the two half-step "
            "velocity updates."
        ),
        default=0,
        param_type="int",
        options={
            0: "Conventional scheme in Amber",
            1: "'Middle' scheme, based on the leapfrog algorithm",
        },
        manual_notes=(
            "\"=0 (default) Conventional scheme in AMBER. =1 'middle' scheme based on the "
            "leapfrog algorithm.\" The manual notes the middle scheme is \"much more "
            "efficient than the default scheme to accurately sample the "
            "configuration/conformation space\" for NVT."
        ),
        related=["ithermostat", "therm_par"],
        common=False,
        category="Integrator",
    )

    db["ithermostat"] = Parameter(
        name="ithermostat",
        brief="Thermostat type for middle scheme",
        what_it_does=(
            "When using the middle scheme (ischeme=1), selects the thermostat type. "
            "Different from ntt -- only applies when ischeme=1."
        ),
        default=1,
        param_type="int",
        options={
            1: "Langevin dynamics (therm_par is the friction coefficient)",
            2: "Andersen thermostat (therm_par is the collision frequency)",
        },
        manual_notes=(
            "\"Flag for different thermostats when the 'middle' scheme is employed. Two "
            "types of thermostats are currently available.\" therm_par supplies the "
            "coupling in both cases -- not gamma_ln."
        ),
        related=["ischeme", "therm_par"],
        common=False,
        category="Integrator",
    )

    db["therm_par"] = Parameter(
        name="therm_par",
        brief="Thermostat parameter for middle scheme (ps^-1)",
        what_it_does=(
            "With the middle scheme (ischeme=1), the thermostat coupling constant: the "
            "friction coefficient for Langevin (ithermostat=1) or the collision frequency "
            "for Andersen (ithermostat=2)."
        ),
        default=5.0,
        param_type="float",
        manual_notes=(
            "The manual states therm_par \"must be set in the input file and should always "
            "be a positive number\", and ties the recommended value to the characteristic "
            "frequency of the specific system rather than giving a fixed range."
        ),
        min_val=0.0,
        related=["ischeme", "ithermostat"],
        common=False,
        category="Integrator",
    )

    # =========================================================================
    # PRESSURE CONTROL
    # =========================================================================
    db["ntp"] = Parameter(
        name="ntp",
        brief="Pressure coupling method",
        what_it_does=(
            "Controls constant pressure (NPT) simulations. "
            "The box volume fluctuates to maintain average pressure. "
            "Essential for equilibrating system density, production at experimental "
            "conditions (1 atm), and membrane simulations (anisotropic)."
        ),
        default=0,
        param_type="int",
        options={
            0: "No pressure control (NVT or NVE)",
            1: "Isotropic scaling (same in all directions) -- MOST COMMON",
            2: "Anisotropic (x,y,z independent) -- for membranes",
            3: "Semi-isotropic (xy coupled, z separate) -- membrane with surface tension",
            4: "Semi-isotropic membrane (standard lipid bilayer)",
        },
        related=["barostat", "pres0", "taup", "ntb"],
        common=True,
        category="Pressure Control",
    )

    db["barostat"] = Parameter(
        name="barostat",
        brief="Barostat algorithm",
        what_it_does=(
            "Algorithm for pressure control. "
            "Berendsen scales coordinates toward target pressure -- fast equilibration "
            "but produces wrong volume fluctuations. "
            "Monte Carlo uses random volume moves with Metropolis criterion -- correct NPT ensemble."
        ),
        default=1,
        param_type="int",
        options={
            1: "Berendsen (weak coupling) -- OK for equilibration only",
            2: "Monte Carlo (correct NPT ensemble) -- RECOMMENDED",
        },
        related=["ntp", "mcbarint"],
        common=True,
        category="Pressure Control",
    )

    db["baro_stochastic"] = Parameter(
        name="baro_stochastic",
        brief="Stochastic cell rescaling barostat",
        what_it_does=(
            "Adds a stochastic term to the Berendsen barostat, giving the stochastic cell "
            "rescaling (SCR) method. Requires ntp>0 and barostat=1."
        ),
        default=0,
        param_type="int",
        options={
            0: "Standard Berendsen barostat",
            1: "Add the stochastic term (stochastic cell rescaling)",
        },
        manual_notes=(
            "\"It must be used with ntp > 0, barostat = 1. This method improves on the "
            "Berendsen barostat and produces the correct isothermal-isobaric ensemble.\" "
            "Note the manual spells the definition heading \"baro_stochasitc\"; its index "
            "and all worked examples use baro_stochastic."
        ),
        related=["ntp", "barostat", "taup"],
        common=False,
        category="Pressure Control",
    )

    db["pres0"] = Parameter(
        name="pres0",
        brief="Target pressure (bar)",
        what_it_does="Reference pressure for NPT simulations. 1 bar ~ 1 atm.",
        default=1.0,
        param_type="float",
        related=["ntp"],
        common=True,
        category="Pressure Control",
    )

    db["taup"] = Parameter(
        name="taup",
        brief="Pressure relaxation time (ps)",
        what_it_does=(
            "How quickly box size adjusts toward target pressure. "
            "Smaller = faster adjustment but potentially less stable."
        ),
        default=1.0,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "\"The recommended value is between 1.0 and 5.0 ps... larger values may sometimes be necessary (if your trajectories seem unstable).\""
        ),
        related=["ntp", "barostat"],
        common=True,
        category="Pressure Control",
    )

    db["comp"] = Parameter(
        name="comp",
        brief="System compressibility (10^-6 bar^-1)",
        what_it_does="Compressibility for pressure coupling. Affects equilibration speed.",
        default=44.6,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "The units are 10^-6 bar^-1; \"a value of 44.6 (default) is appropriate for water\"."
        ),
        related=["ntp"],
        common=False,
        category="Pressure Control",
    )

    db["mcbarint"] = Parameter(
        name="mcbarint",
        brief="MC barostat volume move interval (steps)",
        what_it_does="Steps between Monte Carlo volume change attempts.",
        default=100,
        param_type="int",
        min_val=1,
        related=["barostat"],
        common=False,
        category="Pressure Control",
    )

    db["baroscalingdir"] = Parameter(
        name="baroscalingdir",
        brief="Monte Carlo barostat scaling direction",
        what_it_does=(
            "Restricts which box axis the Monte Carlo barostat rescales. Applies only "
            "with barostat=2 and anisotropic pressure scaling (ntp=2); it has no effect "
            "otherwise."
        ),
        default=0,
        param_type="int",
        options={
            0: "Scale a randomly chosen axis (x, y or z) each step",
            1: "Scale along x only; y and z fixed",
            2: "Scale along y only; x and z fixed",
            3: "Scale along z only; x and y fixed",
        },
        manual_notes=(
            "\"Flag for pressure scaling direction control. Applicable when using Monte "
            "Carlo barostat (barostat = 2) with anisotropic pressure scaling (ntp = 2).\""
        ),
        related=["ntp", "barostat"],
        common=False,
        category="Pressure Control",
    )

    db["csurften"] = Parameter(
        name="csurften",
        brief="Constant surface tension (membrane simulations)",
        what_it_does=(
            "Enables constant surface tension control for membrane simulations. "
            "Specifies which plane contains the membrane interface."
        ),
        default=0,
        param_type="int",
        options={
            0: "No surface tension control",
            1: "Interface in yz plane",
            2: "Interface in xz plane",
            3: "Interface in xy plane (most common for membranes)",
        },
        related=["ntp", "gamma_ten"],
        common=False,
        category="Pressure Control",
    )

    db["ninterface"] = Parameter(
        name="ninterface",
        brief="Number of interfaces in the periodic box",
        what_it_does=(
            "How many interfaces the constant surface tension calculation should assume. "
            "Used with csurften>0."
        ),
        default=2,
        param_type="int",
        min_val=2,
        manual_notes=(
            "\"There must be at least two interfaces in the periodic box. Two interfaces "
            "is appropriate for a lipid bilayer system and is the default value.\""
        ),
        related=["csurften", "gamma_ten"],
        common=False,
        category="Pressure Control",
    )

    db["netfrc"] = Parameter(
        name="netfrc",
        brief="Remove net force each step (PME)",
        what_it_does=(
            "Smooth PME does not strictly conserve momentum, so by default the total force "
            "on the system is removed every step. Set 0 when positional restraints are in "
            "use, since the restraint forces are real external forces."
        ),
        default=1,
        param_type="int",
        options={
            0: "Do not remove the net force",
            1: "Artificially remove the total force each step",
        },
        manual_notes=(
            "\"If netfrc = 1, (the default) the total force on the system is artificially "
            "removed at every step.\" It is set to 0 automatically for minimization. \"It "
            "should also be set to 0 if coordinate restraints are used (ntr>0); note: "
            "setting netfrc=0 when ntr>0 is not now the default, but will be in future "
            "versions.\" NOTE: in ProPrep, positional restraints are normally configured by "
            "the MD Manager's restraint manager rather than here."
        ),
        namelist="ewald",
        related=["ntr", "restraintmask", "restraint_wt"],
        common=False,
        category="Restraints",
    )

    db["nbflag"] = Parameter(
        name="nbflag",
        brief="Non-bonded list update trigger",
        what_it_does=(
            "Chooses how the direct-sum non-bonded list is rebuilt: on a fixed step interval "
            "(nsnb), or whenever an atom has moved more than half the skin width. With the "
            "default, nsnb is ignored."
        ),
        default=1,
        param_type="int",
        options={
            0: "Rebuild every nsnb steps (the 'old' way)",
            1: "Rebuild when any atom has moved more than skinnb/2",
        },
        manual_notes=(
            "\"If nbflag = 1 (the default when imin = 0 or ntb > 0), nsnb is ignored, and "
            "the list is updated whenever any atom has moved more than 1/2 skinnb since the "
            "last list update.\""
        ),
        namelist="ewald",
        related=["nsnb", "skinnb", "cut"],
        common=False,
        category="Non-bonded",
    )

    db["skinnb"] = Parameter(
        name="skinnb",
        brief="Non-bonded list skin width (A)",
        what_it_does=(
            "Extra margin beyond cut used when building the pair list. The list extends to "
            "cut + skinnb while interactions are still truncated at cut, so pairs that drift "
            "inward between rebuilds are already present."
        ),
        default=2.0,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "\"Default is 2.0 A. Use of this parameter is required for energy conservation, "
            "and recommended for all PME runs.\""
        ),
        namelist="ewald",
        related=["cut", "nbflag", "nsnb"],
        common=False,
        category="Non-bonded",
    )

    db["lj1264"] = Parameter(
        name="lj1264",
        brief="12-6-4 Lennard-Jones for divalent metal ions",
        what_it_does=(
            "Enables the 12-6-4 Lennard-Jones form, which adds an r^-4 ion-induced-dipole "
            "term for metal ions. Activated automatically when the topology carries the "
            "C-coefficient; set 0 to force it off."
        ),
        default="auto (from prmtop)",
        param_type="int",
        options={0: "12-6-4 potential inactive", 1: "12-6-4 potential active"},
        manual_notes=(
            "\"In general, you should rarely have to set this variable.\" Default is 1 when "
            "the Lennard-Jones C-coefficient is present in the prmtop and 0 when it is not. "
            "\"Setting lj1264 to 1 when no C-coefficient is present will result in a fatal "
            "error.\" Incompatible with fswitch."
        ),
        related=["plj1264", "fswitch"],
        common=False,
        category="Non-bonded",
    )

    db["plj1264"] = Parameter(
        name="plj1264",
        brief="Pairwise 12-6-4 Lennard-Jones",
        what_it_does=(
            "Pairwise variant of the 12-6-4 form. Activated automatically when the topology "
            "carries the D-coefficient."
        ),
        default="auto (from prmtop)",
        param_type="int",
        options={0: "Pairwise 12-6-4 inactive", 1: "Pairwise 12-6-4 active"},
        manual_notes=(
            "\"Similar to lj1264 above, if not present in the input file, this keyword will "
            "still automatically turn to 1 as long as D-coefficient is found in the prmtop "
            "file.\" Incompatible with fswitch."
        ),
        related=["lj1264", "fswitch"],
        common=False,
        category="Non-bonded",
    )

    db["xcap"] = Parameter(
        name="xcap",
        brief="Water cap centre, x (A)",
        what_it_does="x coordinate of the cap centre. Required when ivcap=1.",
        default=0.0,
        param_type="float",
        manual_notes='"xcap,ycap,zcap Location of the cap center, if ivcap=1 is used." With ivcap=1 the manual says the cap parameters "must be chosen such that the whole solute is covered by solvent".',
        related=["ivcap", "cutcap", "fcap", "ycap", "zcap"],
        common=False,
        category="Water Cap",
    )

    db["ycap"] = Parameter(
        name="ycap",
        brief="Water cap centre, y (A)",
        what_it_does="y coordinate of the cap centre. Required when ivcap=1.",
        default=0.0,
        param_type="float",
        manual_notes='"xcap,ycap,zcap Location of the cap center, if ivcap=1 is used." With ivcap=1 the manual says the cap parameters "must be chosen such that the whole solute is covered by solvent".',
        related=["ivcap", "cutcap", "fcap", "xcap", "zcap"],
        common=False,
        category="Water Cap",
    )

    db["zcap"] = Parameter(
        name="zcap",
        brief="Water cap centre, z (A)",
        what_it_does="z coordinate of the cap centre. Required when ivcap=1.",
        default=0.0,
        param_type="float",
        manual_notes='"xcap,ycap,zcap Location of the cap center, if ivcap=1 is used." With ivcap=1 the manual says the cap parameters "must be chosen such that the whole solute is covered by solvent".',
        related=["ivcap", "cutcap", "fcap", "xcap", "ycap"],
        common=False,
        category="Water Cap",
    )

    db["efn"] = Parameter(
        name="efn",
        brief="Normalise the electric field to box size",
        what_it_does=(
            "When 1, the efx/efy/efz components are divided by the corresponding box "
            "dimension, so the field is expressed relative to the box rather than absolutely."
        ),
        default=0,
        param_type="int",
        options={0: "Use efx/efy/efz as given", 1: "Scale each component to the box dimension"},
        manual_notes=(
            '"If efn is on (efn=1), the x, y, z (efx, efy, efz) components are scaled to '
            'box size... This normalizes the electric field charge to your box size." '
            "pmemd only."
        ),
        related=["efx", "efy", "efz"],
        common=False,
        category="Electric Field",
    )

    db["effreq"] = Parameter(
        name="effreq",
        brief="Electric field oscillation frequency",
        what_it_does=(
            "Frequency term for a time-varying electric field. Left at 0 the field is static."
        ),
        default=0.0,
        param_type="float",
        manual_notes='The oscillating field follows cos((2*pi*effreq/1000)(dt*step) - (pi*efphase/180)). "It currently only supports pmemd (both the serial and MPI versions)."',
        related=["efx", "efy", "efz", "efphase"],
        common=False,
        category="Electric Field",
    )

    db["efphase"] = Parameter(
        name="efphase",
        brief="Electric field phase (degrees)",
        what_it_does="Phase offset, in degrees, for a time-varying electric field.",
        default=0.0,
        param_type="float",
        manual_notes='The oscillating field follows cos((2*pi*effreq/1000)(dt*step) - (pi*efphase/180)). "It currently only supports pmemd (both the serial and MPI versions)."',
        related=["efx", "efy", "efz", "effreq"],
        common=False,
        category="Electric Field",
    )

    db["mdinfo_flush_interval"] = Parameter(
        name="mdinfo_flush_interval",
        brief="mdinfo rewrite interval (seconds)",
        what_it_does=(
            "Minimum interval, in SECONDS, between rewrites of the mdinfo file by pmemd. "
            "Companion to mdout_flush_interval, with a shorter default."
        ),
        default=60,
        param_type="int",
        manual_notes=(
            'pmemd "does an open/close cycle on mdinfo at a default minimum interval of '
            '60 seconds".'
        ),
        min_val=0,
        related=["mdout_flush_interval", "ntpr"],
        common=False,
        category="PMEMD",
    )

    db["rsum_tol"] = Parameter(
        name="rsum_tol",
        brief="Ewald reciprocal sum tolerance",
        what_it_does=(
            "Sets how many reciprocal vectors an Ewald sum uses. Companion to dsum_tol, "
            "which governs the direct sum."
        ),
        default=5e-5,
        param_type="float",
        manual_notes=(
            '"Typically the relative RMS reciprocal sum error is about 5-10 times '
            'RSUM_TOL. Default is 5 x 10^-5."'
        ),
        namelist="ewald",
        related=["dsum_tol", "ew_type", "order"],
        common=False,
        category="Electrostatics",
    )

    db["gamma_ten"] = Parameter(
        name="gamma_ten",
        brief="Surface tension (dyne/cm)",
        what_it_does="Target surface tension for membrane simulations.",
        default=0.0,
        param_type="float",
        min_val=0.0,
        related=["csurften"],
        common=False,
        category="Pressure Control",
    )

    # =========================================================================
    # SYSTEM SETUP
    # =========================================================================
    db["ntb"] = Parameter(
        name="ntb",
        brief="Periodic boundary conditions",
        what_it_does=(
            "Controls whether the system has periodic boundaries (repeating box). "
            "ntb=0: no periodicity (for implicit solvent or gas phase). "
            "ntb=1: periodic, constant volume (NVE or NVT). "
            "ntb=2: periodic, constant pressure (NPT). "
            "Periodic boundaries simulate bulk solution by repeating the box infinitely."
        ),
        default=1,
        param_type="int",
        options={
            0: "No periodicity (implicit solvent, gas phase)",
            1: "Constant volume (NVE, NVT)",
            2: "Constant pressure (NPT)",
        },
        manual_notes=(
            "The manual states there is \"no longer any need to set this variable, since it can be determined from igb and ntp parameters\" (ntb=0 when igb>0, ntb=2 when ntp>0, ntb=1 otherwise), and that overriding it \"is discouraged to prevent errors\"."
        ),
        related=["igb", "ntp", "cut"],
        common=True,
        category="System Setup",
    )

    db["iwrap"] = Parameter(
        name="iwrap",
        brief="Wrap coordinates into primary box",
        what_it_does=(
            "In periodic simulations, molecules can drift outside the primary box. "
            "iwrap=1 shifts them back. This is purely a coordinate convention -- "
            "it doesn't affect the physics since the box repeats infinitely."
        ),
        default=0,
        param_type="int",
        options={
            0: "No wrapping (true displacements preserved)",
            1: "Wrap into primary box (for visualization)",
        },
        common=False,
        category="System Setup",
    )

    # =========================================================================
    # IMPLICIT SOLVENT (GENERALIZED BORN)
    # =========================================================================
    db["igb"] = Parameter(
        name="igb",
        brief="Generalized Born implicit solvent model",
        what_it_does=(
            "Selects the Generalized Born model for implicit solvent. "
            "GB approximates aqueous solvation by treating water as a dielectric continuum. "
            "Much faster than explicit water but less accurate for some properties."
        ),
        default=0,
        param_type="int",
        options={
            0: "No GB (explicit solvent or vacuum)",
            1: "HCT model (older)",
            2: "OBC model (widely used, recommended for general use)",
            5: "OBC with modified parameters",
            6: "Vacuum (no solvent at all)",
            7: "GBneck model",
            8: "GBneck2 (recommended for proteins)",
            10: "Poisson-Boltzmann (see the PB chapter)",
        },
        related=["extdiel", "intdiel", "saltcon", "cut"],
        common=True,
        category="Implicit Solvent",
    )

    db["saltcon"] = Parameter(
        name="saltcon",
        brief="Salt concentration (Molar)",
        what_it_does=(
            "Salt concentration for Debye-Huckel ionic screening in GB. "
            "Affects long-range interactions between charged groups."
        ),
        default=0.0,
        param_type="float",
        min_val=0.0,
        related=["igb"],
        common=True,
        category="Implicit Solvent",
    )

    db["rgbmax"] = Parameter(
        name="rgbmax",
        brief="Maximum distance for GB radii calculation (A)",
        what_it_does="Cutoff for pairwise summation when calculating effective Born radii.",
        default=25.0,
        param_type="float",
        min_val=0.0,
        related=["igb"],
        common=False,
        category="Implicit Solvent",
    )

    db["gbsa"] = Parameter(
        name="gbsa",
        brief="Surface area term for nonpolar solvation",
        what_it_does=(
            "Adds a surface area-dependent term to the GB energy to model nonpolar "
            "(hydrophobic) solvation effects."
        ),
        default=0,
        param_type="int",
        options={
            0: "No surface area term",
            1: "LCPO algorithm (CPU)",
            2: "Recursive algorithm (no forces -- single point only)",
            3: "GPU-optimized GBSA",
        },
        related=["igb", "surften"],
        common=False,
        category="Implicit Solvent",
    )

    db["surften"] = Parameter(
        name="surften",
        brief="Surface tension (kcal/mol/A^2)",
        what_it_does=(
            "Surface tension coefficient for the nonpolar solvation term, applied as "
            "Enp = surften * SA when gbsa=1."
        ),
        manual_notes=(
            "\"Surface tension used to calculate the nonpolar contribution to the free "
            "energy of solvation (when gbsa = 1), as Enp = surften*SA. The default is "
            "0.005 kcal/mol/A2.\" For gbsa = 3 it \"works comparably with gbsa = 1 given "
            "the same value.\""
        ),
        default=0.005,
        param_type="float",
        min_val=0.0,
        related=["gbsa"],
        common=False,
        category="Implicit Solvent",
    )

    db["extdiel"] = Parameter(
        name="extdiel",
        brief="External (solvent) dielectric constant",
        what_it_does="Dielectric constant of the solvent in GB calculations.",
        default=78.5,
        param_type="float",
        min_val=0.0,
        related=["igb", "intdiel"],
        common=False,
        category="Implicit Solvent",
    )

    db["intdiel"] = Parameter(
        name="intdiel",
        brief="Internal (solute) dielectric constant",
        what_it_does="Dielectric constant inside the solute in GB calculations.",
        default=1.0,
        param_type="float",
        min_val=0.0,
        related=["igb", "extdiel"],
        common=False,
        category="Implicit Solvent",
    )

    # =========================================================================
    # CONSTRAINTS (SHAKE)
    # =========================================================================
    db["ntc"] = Parameter(
        name="ntc",
        brief="SHAKE bond constraints",
        what_it_does=(
            "Constrains bond lengths using the SHAKE algorithm. Bond stretching is the "
            "fastest motion in the system, so removing it is what allows a larger time step."
        ),
        default=1,
        param_type="int",
        options={
            1: "No SHAKE (bonds free)",
            2: "Constrain bonds involving hydrogen (STANDARD)",
            3: "Constrain ALL bonds",
        },
        manual_notes=(
            "The manual notes \"typically NTF = NTC\", and that to employ TIP3P you should set NTF = NTC = 2."
        ),
        related=["ntf", "dt", "tol"],
        common=True,
        category="Constraints",
    )

    db["ntf"] = Parameter(
        name="ntf",
        brief="Force evaluation (which interactions to calculate)",
        what_it_does=(
            "Controls which bonded interactions are calculated. "
            "If you constrain bonds with SHAKE (ntc>1), there's no point calculating "
            "forces on those bonds -- SHAKE will enforce the constraint anyway. "
            "So ntf should match ntc to avoid wasted computation."
        ),
        default=1,
        param_type="int",
        options={
            1: "Complete interaction is calculated",
            2: "Omit bond interactions involving H (use with ntc=2)",
            3: "Omit all bond interactions (use with ntc=3)",
            4: "Omit angles involving H, and all bonds",
            5: "Omit all bond and angle interactions",
            6: "Omit dihedrals involving H, and all bonds and angles",
            7: "Omit all bond, angle and dihedral interactions",
            8: "Omit all bond, angle, dihedral and non-bonded interactions",
        },
        manual_notes=(
            "\"If SHAKE is used (see NTC), it is not necessary to calculate forces for "
            "the constrained bonds.\" The manual notes that typically NTF = NTC, and that "
            "TIP3P requires NTF = NTC = 2."
        ),
        related=["ntc"],
        common=True,
        category="Constraints",
    )

    db["tol"] = Parameter(
        name="tol",
        brief="SHAKE tolerance",
        what_it_does=(
            "Relative tolerance for SHAKE convergence. "
            "SHAKE iteratively adjusts positions until bond lengths are within tolerance."
        ),
        default=0.00001,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "\"Recommended maximum: <0.00005 Angstrom.\""
        ),
        related=["ntc"],
        common=False,
        category="Constraints",
    )

    db["jfastw"] = Parameter(
        name="jfastw",
        brief="Fast water SHAKE method",
        what_it_does=(
            "By default the system is searched for water residues and special fast SHAKE "
            "routines are used for them. Set to 4 to disable those routines."
        ),
        default=0,
        param_type="int",
        options={
            0: "Use SETTLE for water (faster, default)",
            4: "Use iterative SHAKE for water",
        },
        related=["ntc"],
        common=False,
        category="Constraints",
    )

    db["noshakemask"] = Parameter(
        name="noshakemask",
        brief="Atoms to exclude from SHAKE",
        what_it_does=(
            "Amber mask selecting atoms that should NOT be constrained even when ntc>1. "
            "Useful for QM/MM or free energy calculations. "
            "When used, ntf is automatically set to 1."
        ),
        default="",
        param_type="str",
        related=["ntc"],
        common=False,
        category="Constraints",
    )

    # =========================================================================
    # NON-BONDED INTERACTIONS
    # =========================================================================
    db["cut"] = Parameter(
        name="cut",
        brief="Non-bonded cutoff distance (A)",
        what_it_does=(
            "The distance beyond which non-bonded interactions are truncated or handled differently. "
            "For PME (explicit solvent): this is the direct-space cutoff; electrostatics "
            "beyond this are handled by reciprocal space (very accurate). "
            "For GB (implicit solvent): this actually truncates interactions, so use a very large value."
        ),
        default=8.0,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "For PME \"the cutoff is used to limit direct space sum, and 8.0 is usually a good value\". When igb>0 the cutoff truncates non-bonded pairs and \"a larger value than the default is generally required\"."
        ),
        related=["igb", "ntb"],
        common=True,
        category="Non-bonded",
    )

    db["fswitch"] = Parameter(
        name="fswitch",
        brief="Force switching start distance (A)",
        what_it_does=(
            "Enables smooth switching of van der Waals forces near the cutoff. "
            "Force switching smoothly reduces forces to zero between fswitch and cut. "
            "CHARMM force fields are parameterized with force switching. "
            "fswitch < 0: no switching (hard cutoff) -- AMBER default."
        ),
        default=-1.0,
        param_type="float",
        manual_notes=(
            "Not supported with GB (only igb=0 and ntb>0), and incompatible with the 12-6-4 and pairwise 12-6-4 Lennard-Jones models. The manual notes about a 20% performance cost when force switching is on."
        ),
        related=["cut"],
        common=False,
        category="Non-bonded",
    )

    db["dielc"] = Parameter(
        name="dielc",
        brief="Dielectric constant multiplier",
        what_it_does=(
            "Multiplies all electrostatic interactions by 1/dielc. "
            "This is a CRUDE approximation -- not related to GB dielectric constants."
        ),
        default=1.0,
        param_type="float",
        min_val=0.0,
        manual_notes=(
            "Not related to the GB or PB dielectric constants. The manual says it \"should only be used for quasi-vacuum simulations\"."
        ),
        common=False,
        category="Non-bonded",
    )

    db["nsnb"] = Parameter(
        name="nsnb",
        brief="Neighbor list update frequency",
        what_it_does=(
            "How often the non-bonded pair list is rebuilt. Only consulted when igb=0 "
            "and nbflag=0; otherwise the list is updated automatically."
        ),
        default=25,
        param_type="int",
        manual_notes=(
            "\"Determines the frequency of nonbonded list updates when igb=0 and "
            "nbflag=0.\" Default is 25."
        ),
        min_val=0,
        common=False,
        category="Non-bonded",
    )

    # =========================================================================
    # ELECTROSTATICS / PME
    # =========================================================================
    db["ew_type"] = Parameter(
        name="ew_type",
        brief="Ewald method type",
        what_it_does=(
            "Selects particle mesh Ewald (the interpolated, approximate method used "
            "for essentially all production work) or an exact Ewald summation."
        ),
        default=0,
        param_type="int",
        options={
            0: "Particle mesh Ewald (PME) -- standard",
            1: "Exact Ewald summation -- accuracy check only",
        },
        manual_notes=(
            "\"Standard use is to have EW_TYPE = 0 which turns on the particle mesh "
            "ewald (PME) method.\" The exact summation \"is present mainly to serve as "
            "an accuracy check\"; it may be faster below ~500 atoms, but for larger "
            "systems PME is significantly faster."
        ),
        namelist="ewald",
        related=["nfft1", "nfft2", "nfft3", "order"],
        common=False,
        category="Electrostatics",
    )

    db["order"] = Parameter(
        name="order",
        brief="PME interpolation order",
        what_it_does="Order of B-spline interpolation for PME. 4 = cubic spline.",
        default=4,
        param_type="int",
        min_val=3,
        namelist="ewald",
        manual_notes=(
            "Minimum order is 3; order 4 \"implies a cubic spline approximation which is a good standard value\". The cost of PME goes roughly as the order cubed."
        ),
        common=False,
        category="Electrostatics",
    )

    db["nfft1"] = Parameter(
        name="nfft1",
        brief="PME grid points in X",
        what_it_does="Number of grid points for PME in the X direction.",
        default=0,
        param_type="int",
        min_val=0,
        namelist="ewald",
        manual_notes=(
            "Reasonable results are obtained when nfft1/2/3 are approximately equal to the box dimensions A, B and C respectively."
        ),
        related=["nfft2", "nfft3", "order"],
        common=False,
        category="Electrostatics",
    )

    db["nfft2"] = Parameter(
        name="nfft2",
        brief="PME grid points in Y",
        what_it_does="Number of grid points for PME in the Y direction.",
        default=0,
        param_type="int",
        min_val=0,
        namelist="ewald",
        related=["nfft1", "nfft3"],
        common=False,
        category="Electrostatics",
    )

    db["nfft3"] = Parameter(
        name="nfft3",
        brief="PME grid points in Z",
        what_it_does="Number of grid points for PME in the Z direction.",
        default=0,
        param_type="int",
        min_val=0,
        namelist="ewald",
        related=["nfft1", "nfft3"],
        common=False,
        category="Electrostatics",
    )

    db["dsum_tol"] = Parameter(
        name="dsum_tol",
        brief="PME direct sum tolerance",
        what_it_does="Controls accuracy of direct space sum in PME.",
        default=1.0e-5,
        param_type="float",
        min_val=0.0,
        namelist="ewald",
        manual_notes=(
            "\"Standard values for DSUM_TOL are in the range of 10^-6 to 10^-5.\" The relative RMS force error is roughly 10-50 times dsum_tol."
        ),
        common=False,
        category="Electrostatics",
    )

    # =========================================================================
    # ELECTRIC FIELD
    # =========================================================================
    db["efx"] = Parameter(
        name="efx",
        brief="External electric field X component (kcal/mol/A/e)",
        what_it_does="Applies a constant electric field force on charged atoms in the X direction.",
        default=0.0,
        param_type="float",
        manual_notes=(
            "Electric fields are off when efx, efy and efz are all 0. The manual notes this \"currently only supports pmemd (both the serial and MPI versions)\"."
        ),
        related=["efy", "efz"],
        common=False,
        category="Electric Field",
    )

    db["efy"] = Parameter(
        name="efy",
        brief="External electric field Y component (kcal/mol/A/e)",
        what_it_does="Applies a constant electric field force on charged atoms in the Y direction.",
        default=0.0,
        param_type="float",
        manual_notes=(
            "Electric fields are off when efx, efy and efz are all 0. The manual notes this \"currently only supports pmemd (both the serial and MPI versions)\"."
        ),
        related=["efx", "efz"],
        common=False,
        category="Electric Field",
    )

    db["efz"] = Parameter(
        name="efz",
        brief="External electric field Z component (kcal/mol/A/e)",
        what_it_does="Applies a constant electric field force on charged atoms in the Z direction.",
        default=0.0,
        param_type="float",
        manual_notes=(
            "Electric fields are off when efx, efy and efz are all 0. The manual notes this \"currently only supports pmemd (both the serial and MPI versions)\"."
        ),
        related=["efx", "efy"],
        common=False,
        category="Electric Field",
    )

    # =========================================================================
    # OUTPUT CONTROL
    # =========================================================================
    db["ntpr"] = Parameter(
        name="ntpr",
        brief="Energy output frequency (steps)",
        what_it_does=(
            "How often to print energy information to mdout and mdinfo. "
            "This is your main way to monitor the simulation's progress and health. "
            "Note that on GPUs small ntpr values cost performance, because energies must be "
            "copied from GPU memory."
        ),
        default=50,
        param_type="int",
        min_val=1,
        related=["ntwe", "ntwx"],
        common=True,
        category="Output Control",
    )

    db["ntwx"] = Parameter(
        name="ntwx",
        brief="Coordinate trajectory output frequency",
        what_it_does=(
            "How often to save coordinates to the trajectory file (mdcrd). "
            "The trajectory is your record of the simulation -- used for all analysis. "
            "Consider storage: 1 ns at dt=0.002, ntwx=500 = 1000 frames."
        ),
        default=0,
        param_type="int",
        min_val=0,
        related=["ntpr", "ioutfm"],
        common=True,
        category="Output Control",
    )

    db["ntwr"] = Parameter(
        name="ntwr",
        brief="Restart file write frequency",
        what_it_does=(
            "How often to write the restrt file during the run. A restrt file is "
            "always written at the end of the run regardless of this value. Negative "
            "values write a uniquely named restrt_<nstep> file every abs(ntwr) steps."
        ),
        default="nstlim",
        param_type="int",
        manual_notes=(
            "\"Default = nstlim.\" \"No matter what the value of ntwr, a restrt file "
            "will be written at the end of the run.\""
        ),
        related=["nstlim"],
        common=True,
        category="Output Control",
    )

    db["ntwe"] = Parameter(
        name="ntwe",
        brief="Energy file (mden) output",
        what_it_does=(
            "Write energies and temperatures to the compact mden file every ntwe steps. "
            "0 disables it. Note that mden energies are one time step ahead of the "
            "coordinates in mdcrd/mdvel, so they are not directly comparable frame by frame."
        ),
        default=0,
        param_type="int",
        manual_notes=(
            "mden energies \"are not synchronized with coordinates or velocities\"; "
            "assuming identical ntwe and ntwx \"the energies are one time step before the "
            "coordinates\". \"Consequently, an mden file is rarely written.\""
        ),
        min_val=0,
        common=False,
        category="Output Control",
    )

    db["ntwv"] = Parameter(
        name="ntwv",
        brief="Velocity trajectory output",
        what_it_does=(
            "How often to save velocities. Most users never need velocity trajectories. "
            "Useful for kinetic energy analysis and velocity autocorrelation functions."
        ),
        default=0,
        param_type="int",
        min_val=-1,
        manual_notes=(
            "ntwv=-1 writes velocities into mdcrd, making it a combined coordinate/velocity file, at the ntwx interval; this needs ioutfm=1. \"Most users will have no need for a velocity trajectory file.\""
        ),
        common=False,
        category="Output Control",
    )

    db["ntwf"] = Parameter(
        name="ntwf",
        brief="Force trajectory output",
        what_it_does="How often to write forces to the mdfrc file. 0 disables it.",
        default=0,
        param_type="int",
        min_val=-1,
        manual_notes=(
            "ntwf=-1 writes forces into mdcrd at the ntwx interval; this needs ioutfm=1. \"Most users will have no need for a force trajectory file.\""
        ),
        common=False,
        category="Output Control",
    )

    db["ntave"] = Parameter(
        name="ntave",
        brief="Running average frequency",
        what_it_does=(
            "Controls printing of running energy averages and fluctuations. "
            "If ntave > 0, prints averages every ntave steps. "
            "Note that on GPUs a non-zero ntave forces an energy calculation every step."
        ),
        default=0,
        param_type="int",
        min_val=0,
        manual_notes=(
            "\"Setting ntave to a value 1/2 or 1/4 of nstlim provides a simple way to look at convergence during the simulation.\""
        ),
        related=["ntpr"],
        common=False,
        category="Output Control",
    )

    db["ioutfm"] = Parameter(
        name="ioutfm",
        brief="Trajectory file format",
        what_it_does=(
            "Format for the coordinate and velocity trajectory files (mdcrd, mdvel and "
            "inptraj). NetCDF is a self-describing binary format."
        ),
        default=1,
        param_type="int",
        options={
            0: "ASCII formatted (legacy, larger, slower)",
            1: "NetCDF binary (recommended)",
        },
        common=False,
        category="Output Control",
    )

    db["ntxo"] = Parameter(
        name="ntxo",
        brief="Restart file format",
        what_it_does=(
            "Format for the restart file. "
            "NetCDF is smaller, faster, and maintains full precision."
        ),
        default=2,
        param_type="int",
        options={
            1: "ASCII formatted (human-readable, larger)",
            2: "NetCDF binary (recommended: smaller, faster, full precision)",
        },
        common=False,
        category="Output Control",
    )

    db["ntwprt"] = Parameter(
        name="ntwprt",
        brief="Number of atoms to include in trajectory",
        what_it_does=(
            "Controls which atoms are written to trajectory files. "
            "ntwprt=0 writes ALL atoms. Setting to a positive integer writes only atoms 1 to ntwprt. "
            "Very useful for excluding solvent from trajectories to reduce file size."
        ),
        default=0,
        param_type="int",
        min_val=0,
        manual_notes=(
            "ntwprt=0 includes all atoms; a positive value includes only atoms 1 to ntwprt."
        ),
        related=["ntwx"],
        common=False,
        category="Output Control",
    )

    db["ionstepvelocities"] = Parameter(
        name="ionstepvelocities",
        brief="On-step velocities vs half-step",
        what_it_does="Whether to output on-step velocities (1) or half-step velocities (0) in trajectory.",
        default=0,
        param_type="int",
        options={0: "Half-step velocities", 1: "On-step velocities"},
        common=False,
        category="Output Control",
    )

    # =========================================================================
    # RESTRAINTS
    # =========================================================================
    db["ntr"] = Parameter(
        name="ntr",
        brief="Positional restraints",
        what_it_does=(
            "Enables harmonic positional restraints on selected atoms. "
            "Atoms are restrained to reference coordinates from a separate file. "
            "Energy penalty: E = k x (displacement)^2."
        ),
        default=0,
        param_type="int",
        options={
            0: "No positional restraints",
            1: "Apply harmonic restraints to selected atoms",
        },
        manual_notes=(
            "The manual advises that \"If ntr=1, you should also set netfrc=0\"."
        ),
        related=["restraint_wt", "restraintmask"],
        common=True,
        category="Restraints",
    )

    db["restraint_wt"] = Parameter(
        name="restraint_wt",
        brief="Restraint force constant (kcal/mol/A^2)",
        what_it_does="Strength of positional restraints. Higher = atoms stay closer to reference.",
        default=0.0,
        param_type="float",
        min_val=0.0,
        related=["ntr", "restraintmask"],
        common=True,
        category="Restraints",
    )

    db["restraintmask"] = Parameter(
        name="restraintmask",
        brief="Atom selection for restraints",
        what_it_does="Amber mask string specifying which atoms to restrain.",
        default="",
        param_type="str",
        manual_notes=(
            "Mask strings are limited to a maximum of 256 characters."
        ),
        related=["ntr", "restraint_wt"],
        common=True,
        category="Restraints",
    )

    # =========================================================================
    # ADVANCED
    # =========================================================================
    db["ibelly"] = Parameter(
        name="ibelly",
        brief="Belly dynamics (freeze atoms)",
        what_it_does=(
            "Freezes atoms NOT in bellymask -- only bellymask atoms move. "
            "Legacy feature -- positional restraints (ntr=1) are usually preferred."
        ),
        default=0,
        param_type="int",
        options={
            0: "Normal dynamics (all atoms move)",
            1: "Only bellymask atoms move",
        },
        manual_notes=(
            "Not available when igb>0. The manual notes belly dynamics \"does not provide any significant speed advantage\"."
        ),
        related=["bellymask"],
        common=False,
        category="Advanced",
    )

    db["bellymask"] = Parameter(
        name="bellymask",
        brief="Moving atoms for belly dynamics",
        what_it_does="Mask specifying atoms that ARE allowed to move (opposite of restraintmask).",
        default="",
        param_type="str",
        manual_notes=(
            "Mask strings are limited to a maximum of 256 characters."
        ),
        related=["ibelly"],
        common=False,
        category="Advanced",
    )

    db["idecomp"] = Parameter(
        name="idecomp",
        brief="Energy decomposition",
        what_it_does=(
            "Enables per-residue or pairwise per-residue energy decomposition. "
            "Useful for identifying key residue-residue interactions."
        ),
        default=0,
        param_type="int",
        options={
            0: "No decomposition",
            1: "Per-residue (1-4 terms with internal)",
            2: "Per-residue (1-4 terms with non-bonded)",
            3: "Pairwise per-residue (1-4 with internal)",
            4: "Pairwise per-residue (1-4 with non-bonded)",
        },
        common=False,
        category="Advanced",
    )

    # =========================================================================
    # WATER CAP
    # =========================================================================
    db["ivcap"] = Parameter(
        name="ivcap",
        brief="Water cap type",
        what_it_does="Cap type for spherical boundary water simulations.",
        default=0,
        param_type="int",
        options={
            0: "Cap in effect if present in the prmtop file",
            1: "Excise a cap from a larger water box (needs cutcap, xcap, ycap, zcap)",
            2: "Cap inactivated, even if parameters are present in the prmtop",
            5: "Excise a shell of water around the solute (cutcap = shell thickness)",
        },
        manual_notes=(
            "\"= 0 Cap will be in effect if it is in the prmtop file (default).\" For best "
            "physical realism the manual recommends combining the cap with igb=10, to "
            "include the reaction field of waters beyond the cap radius."
        ),
        common=False,
        category="Water Cap",
    )

    db["fcap"] = Parameter(
        name="fcap",
        brief="Cap force constant (kcal/mol/A^2)",
        what_it_does="Force constant for the water cap restraining potential.",
        default=1.5,
        param_type="float",
        min_val=0.0,
        related=["ivcap", "cutcap"],
        common=False,
        category="Water Cap",
    )

    db["cutcap"] = Parameter(
        name="cutcap",
        brief="Cap radius (A)",
        what_it_does="Radius of the spherical water cap.",
        default=15.0,
        param_type="float",
        min_val=0.0,
        related=["ivcap", "fcap"],
        common=False,
        category="Water Cap",
    )

    # =========================================================================
    # PMEMD-SPECIFIC
    # =========================================================================
    db["mdout_flush_interval"] = Parameter(
        name="mdout_flush_interval",
        brief="PMEMD output flush interval",
        what_it_does=(
            "Minimum interval, in SECONDS, between rewrites of the mdout file by pmemd. "
            "pmemd does an open/close cycle rather than a flush. Setting 0 reopens mdout "
            "for every printed step."
        ),
        default=300,
        param_type="int",
        manual_notes=(
            "pmemd \"does an open/close cycle on mdout at a default minimum interval of "
            "300 seconds\", changeable in the range 0-3600."
        ),
        min_val=0,
        common=False,
        category="PMEMD",
    )

    db["es_cutoff"] = Parameter(
        name="es_cutoff",
        brief="PMEMD electrostatic cutoff (A)",
        what_it_does="Separate cutoff for electrostatic interactions (PMEMD only).",
        default=0.0,
        param_type="float",
        manual_notes=(
            "If es_cutoff/vdw_cutoff are specified you \"should not specify the cut variable\", and \"there is a requirement that vdw_cutoff >= es_cutoff\"."
        ),
        related=["cut", "vdw_cutoff"],
        common=False,
        category="PMEMD",
    )

    db["vdw_cutoff"] = Parameter(
        name="vdw_cutoff",
        brief="PMEMD van der Waals cutoff (A)",
        what_it_does="Separate cutoff for van der Waals interactions (PMEMD only).",
        default=0.0,
        param_type="float",
        manual_notes=(
            "If es_cutoff/vdw_cutoff are specified you \"should not specify the cut variable\", and \"there is a requirement that vdw_cutoff >= es_cutoff\"."
        ),
        related=["cut", "es_cutoff"],
        common=False,
        category="PMEMD",
    )

    db["fft_grids_per_ang"] = Parameter(
        name="fft_grids_per_ang",
        brief="PMEMD FFT grids per angstrom",
        what_it_does=(
            "Requested reciprocal-space FFT grid density in grids per angstrom. nfft1/2/3 "
            "are chosen automatically to meet or exceed it, subject to the FFT's supported "
            "prime factors."
        ),
        default=1.0,
        param_type="float",
        manual_notes=(
            "\"The default value is 1.0 grids/angstrom and gives very reasonable accuracy.\""
        ),
        namelist="ewald",
        min_val=0.1,
        common=False,
        category="PMEMD",
    )

    return db


# =============================================================================
# Utility Functions
# =============================================================================
def get_parameters_by_category(db: Dict[str, Parameter]) -> Dict[str, list]:
    """Group parameters by category, preserving a logical order.

    Returns:
        Ordered dict of category_name -> list of Parameter instances.
    """
    category_order = [
        "Simulation Control",
        "Input/Restart",
        "Minimization",
        "Molecular Dynamics",
        "Temperature Control",
        "Integrator",
        "Pressure Control",
        "System Setup",
        "Implicit Solvent",
        "Constraints",
        "Non-bonded",
        "Electrostatics",
        "Electric Field",
        "Output Control",
        "Restraints",
        "Advanced",
        "Water Cap",
        "PMEMD",
    ]

    grouped: Dict[str, list] = {}
    for cat in category_order:
        grouped[cat] = []

    for param in db.values():
        cat = param.category
        if cat not in grouped:
            grouped[cat] = []
        grouped[cat].append(param)

    # Remove empty categories
    return {cat: params for cat, params in grouped.items() if params}
