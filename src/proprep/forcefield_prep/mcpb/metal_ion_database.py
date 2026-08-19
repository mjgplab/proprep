"""
Metal Ion Database

Auto-configuration database for metal ions based on MCPB/PymSMT data.
Provides element, mass, and VDW parameters for standard metal ions without
requiring user input.

Data sourced from:
- METAL_PDB: PymSMT element.py
- IOD Parameters: Li et al. JCTC, 2013, 9, 2733; JCTC, 2015, 11, 1645; JCTC, 2020, 16, 4429

© 2024 ProPrep Developer. All rights reserved.
"""

from typing import Optional, Tuple, Dict
from dataclasses import dataclass
import logging


# Water models with their own metal LJ set. Anything else has to borrow one.
SUPPORTED_WATER_MODELS = ('tip3p', 'spce', 'tip4pew', 'opc3', 'opc', 'fb3', 'fb4')

# Models ProPrep offers that have no LJ set of their own, and the closest one
# that does. tip5p has no counterpart at all, so it takes the tip3p default.
_WATER_MODEL_SUBSTITUTES = {
    'opc3pol': 'opc3',   # polarizable OPC3
    'spceb': 'spce',     # SPC/E-b
}


def water_model_from_leaprc(leaprc: str, logger=None) -> str:
    """The LJ-table water model named by a ``leaprc.water.<model>`` string.

    Metal LJ parameters are fitted per water model -- Fe3+ is (1.386, 0.01357)
    against TIP3P but (1.400, 0.01571) against OPC -- so the model the system
    will actually be solvated in decides which set is correct.

    Matched on the exact leaprc suffix. A substring test would send
    ``leaprc.water.opc3`` to the OPC set, since "opc" is inside "opc3".
    """
    log = logger or logging.getLogger(__name__)
    name = (leaprc or '').strip().lower().rsplit('.', 1)[-1]

    if not name:
        return 'tip3p'
    if name in SUPPORTED_WATER_MODELS:
        return name
    if name in _WATER_MODEL_SUBSTITUTES:
        substitute = _WATER_MODEL_SUBSTITUTES[name]
        log.warning(
            f"No metal Lennard-Jones set is published for the {name} water "
            f"model; using the {substitute} set.")
        return substitute

    log.warning(
        f"No metal Lennard-Jones set is published for the {name} water model; "
        f"using tip3p parameters.")
    return 'tip3p'


@dataclass
class MetalConfig:
    """Complete configuration for a metal ion"""
    resname: str           # PDB residue name (e.g., 'FE', 'FE2')
    atom_name: str         # PDB atom name (e.g., 'FE', 'FE2')
    element: str           # Element symbol (e.g., 'Fe')
    atomic_number: int     # Atomic number
    mass: float            # Atomic mass (amu)
    atom_type: str         # AMBER atom type (same as resname for metals)
    charge: int            # Formal charge (e.g., +2, +3)
    spin: int              # Number of unpaired electrons (+alpha, -beta)
    vdw_radius: float      # VDW radius (Å) - will be set based on water model
    vdw_epsilon: float     # VDW well depth (kcal/mol) - will be set based on water model
    vdw_source: str        # Citation for VDW parameters
    water_model: str       # Water model used for VDW params


class MetalIonDatabase:
    """
    Auto-configuration database for metal ions.

    Recognizes metals from (resname, atomname) pairs and provides:
    - Element symbol and atomic properties
    - AMBER atom type assignment
    - VDW (non-bonded) parameters for different water models

    No user input required for standard metals.
    """

    # Metal PDB mappings: (resname, atomname) -> (element, atomic_num, mass)
    # Sourced from PymSMT element.py METAL_PDB
    METAL_PDB = {
        ('LI', 'LI'):    ('Li', 3, 6.94),
        ('BE', 'BE'):    ('Be', 4, 9.01),
        ('F', 'F'):      ('F', 9, 19.00),
        ('NA', 'NA'):    ('Na', 11, 22.99),
        ('MG', 'MG'):    ('Mg', 12, 24.305),
        ('AL', 'AL'):    ('Al', 13, 26.98),
        ('SI', 'SI'):    ('Si', 14, 28.09),
        ('CL', 'CL'):    ('Cl', 17, 35.45),
        ('K', 'K'):      ('K', 19, 39.10),
        ('CA', 'CA'):    ('Ca', 20, 40.08),
        ('SC', 'SC'):    ('Sc', 21, 44.96),
        ('TI', 'TI'):    ('Ti', 22, 47.87),
        ('V', 'V'):      ('V', 23, 50.94),
        ('CR', 'CR'):    ('Cr', 24, 52.00),
        ('MN', 'MN'):    ('Mn', 25, 54.94),
        ('FE2', 'FE2'):  ('Fe', 26, 55.85),  # Fe2+
        ('FE', 'FE'):    ('Fe', 26, 55.85),  # Fe3+
        ('CO', 'CO'):    ('Co', 27, 58.93),
        ('NI', 'NI'):    ('Ni', 28, 58.69),
        ('CU1', 'CU'):   ('Cu', 29, 63.55),  # Cu+
        ('CU', 'CU'):    ('Cu', 29, 63.55),  # Cu2+
        ('ZN', 'ZN'):    ('Zn', 30, 65.4),
        ('GA', 'GA'):    ('Ga', 31, 69.72),
        ('GE', 'GE'):    ('Ge', 32, 72.63),
        ('ARS', 'AS'):   ('As', 33, 74.92),
        ('SE', 'SE'):    ('Se', 34, 78.97),
        ('BR', 'BR'):    ('Br', 35, 79.90),
        ('RB', 'RB'):    ('Rb', 37, 85.47),
        ('SR', 'SR'):    ('Sr', 38, 87.62),
        ('Y', 'Y'):      ('Y', 39, 88.91),
        ('ZR', 'ZR'):    ('Zr', 40, 91.22),
        ('NB', 'NB'):    ('Nb', 41, 92.91),
        ('MO', 'MO'):    ('Mo', 42, 95.96),
        ('TC', 'TC'):    ('Tc', 43, 98.0),
        ('RU', 'RU'):    ('Ru', 44, 101.07),
        ('RH1', 'RH'):   ('Rh', 45, 102.91),  # Rh+
        ('RH', 'RH'):    ('Rh', 45, 102.91),
        ('RH3', 'RH'):   ('Rh', 45, 102.91),  # Rh3+
        ('PD', 'PD'):    ('Pd', 46, 106.42),
        ('AG', 'AG'):    ('Ag', 47, 107.87),
        ('CD', 'CD'):    ('Cd', 48, 112.41),
        ('IN', 'IN'):    ('In', 49, 114.82),
        ('SN', 'SN'):    ('Sn', 50, 118.7),
        ('SB', 'SB'):    ('Sb', 51, 121.8),
        ('TE', 'TE'):    ('Te', 52, 127.6),
        ('IOD', 'I'):    ('I', 53, 126.9),
        ('CS', 'CS'):    ('Cs', 55, 132.91),
        ('BA', 'BA'):    ('Ba', 56, 137.33),
        ('LA', 'LA'):    ('La', 57, 138.91),
        ('CE', 'CE'):    ('Ce', 58, 140.12),
        ('PR', 'PR'):    ('Pr', 59, 140.91),
        ('ND', 'ND'):    ('Nd', 60, 144.2),
        ('PM', 'PM'):    ('Pm', 61, 145.0),
        ('SM', 'SM'):    ('Sm', 62, 150.36),
        ('EU', 'EU'):    ('Eu', 63, 151.96),  # Eu2+
        ('EU3', 'EU'):   ('Eu', 63, 151.96),  # Eu3+
        ('GD3', 'GD'):   ('Gd', 64, 157.25),  # Gd3+
        ('GD', 'GD'):    ('Gd', 64, 157.25),
        ('TB', 'TB'):    ('Tb', 65, 158.93),
        ('DY', 'DY'):    ('Dy', 66, 162.5),
        ('HO3', 'HO'):   ('Ho', 67, 164.93),  # Ho3+
        ('HO', 'HO'):    ('Ho', 67, 164.93),
        ('ER', 'ER'):    ('Er', 68, 167.3),
        ('TM', 'TM'):    ('Tm', 69, 168.9),
        ('YB2', 'YB2'):  ('Yb', 70, 173.05),  # Yb2+
        ('YB', 'YB'):    ('Yb', 70, 173.05),  # Yb3+
        ('LU', 'LU'):    ('Lu', 71, 174.97),
        ('HF', 'HF'):    ('Hf', 72, 178.5),
        ('TA', 'TA'):    ('Ta', 73, 180.9),
        ('W', 'W'):      ('W', 74, 183.84),
        ('RE', 'RE'):    ('Re', 75, 186.21),
        ('OS4', 'OS'):   ('Os', 76, 190.23),  # Os4+
        ('OS', 'OS'):    ('Os', 76, 190.23),
        ('IR3', 'IR'):   ('Ir', 77, 192.22),  # Ir3+
        ('IR', 'IR'):    ('Ir', 77, 192.22),
        ('PT', 'PT'):    ('Pt', 78, 195.08),
        ('AU', 'AU'):    ('Au', 79, 197.0),
        ('HG', 'HG'):    ('Hg', 80, 200.59),
        ('TL', 'TL'):    ('Tl', 81, 204.38),
        ('PB', 'PB'):    ('Pb', 82, 207.2),
        ('BI', 'BI'):    ('Bi', 83, 209.0),
        ('PO', 'PO'):    ('Po', 84, 209.0),
        ('AT', 'AT'):    ('At', 85, 210.0),
        ('FR', 'FR'):    ('Fr', 87, 223.0),
        ('RA', 'RA'):    ('Ra', 88, 226.0),
        ('AC', 'AC'):    ('Ac', 89, 227.0),
        ('TH', 'TH'):    ('Th', 90, 232.04),
        ('PA', 'PA'):    ('Pa', 91, 231.0),
        ('U', 'U'):      ('U', 92, 238.0),
        ('NP', 'NP'):    ('Np', 93, 237.0),
        ('PU', 'PU'):    ('Pu', 94, 244.0),
        ('AM', 'AM'):    ('Am', 95, 243.0),
        ('CM', 'CM'):    ('Cm', 96, 247.0),
        ('BK', 'BK'):    ('Bk', 97, 247.0),
        ('CF', 'CF'):    ('Cf', 98, 251.0),
        ('ES', 'ES'):    ('Es', 99, 252.0),
        ('FM', 'FM'):    ('Fm', 100, 257.0),
        ('MD', 'MD'):    ('Md', 101, 258.0),
        ('NO', 'NO'):    ('No', 102, 259.0),
        ('LR', 'LR'):    ('Lr', 103, 262.0),
    }

    # IOD Parameters for monovalent ions
    # Format: {ion_type: (radius, epsilon, citation)}
    IOD_MONOVALENT = {
        'tip3p': {
            'Li1': (1.315, 0.00594975, 'IOD set for Li+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'Na1': (1.465, 0.02909167, 'IOD set for Na+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'K1': (1.745, 0.17018074, 'IOD set for K+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'Rb1': (1.820, 0.22962229, 'IOD set for Rb+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'Cs1': (2.000, 0.38943250, 'IOD set for Cs+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'Tl1': (1.870, 0.27244486, 'IOD set for Tl+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'Cu1': (1.214, 0.00139196, 'IOD set for Cu+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'Ag1': (1.500, 0.03899838, 'IOD set for Ag+ ion from Li et al. JCTC, 2015, 11, 1645'),
            'F-1': (1.739, 0.16573832, 'IOD set for F- ion from Li et al. JCTC, 2015, 11, 1645'),
            'Cl-1': (2.162, 0.53154665, 'IOD set for Cl- ion from Li et al. JCTC, 2015, 11, 1645'),
            'Br-1': (2.331, 0.65952968, 'IOD set for Br- ion from Li et al. JCTC, 2015, 11, 1645'),
            'I-1': (2.590, 0.80293907, 'IOD set for I- ion from Li et al. JCTC, 2015, 11, 1645'),
        }
    }

    # IOD Parameters for divalent ions
    IOD_DIVALENT = {
        'tip3p': {
            'Be2': (1.168, 0.00063064, 'IOD set for Be2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Cu2': (1.409, 0.01721000, 'IOD set for Cu2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Ni2': (1.373, 0.01179373, 'IOD set for Ni2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Zn2': (1.395, 0.01491700, 'IOD set for Zn2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Co2': (1.404, 0.01636246, 'IOD set for Co2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Cr2': (1.388, 0.01386171, 'IOD set for Cr2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Fe2': (1.409, 0.01721000, 'IOD set for Fe2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Mg2': (1.395, 0.01491700, 'IOD set for Mg2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'V2': (1.476, 0.03198620, 'IOD set for V2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Mn2': (1.467, 0.02960343, 'IOD set for Mn2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Hg2': (1.575, 0.06751391, 'IOD set for Hg2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Cd2': (1.506, 0.04090549, 'IOD set for Cd2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Ca2': (1.608, 0.08337961, 'IOD set for Ca2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Sn2': (1.738, 0.16500296, 'IOD set for Sn2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Sr2': (1.753, 0.17618319, 'IOD set for Sr2+ ion from Li et al. JCTC, 2013, 9, 2733'),
            'Ba2': (1.913, 0.31060194, 'IOD set for Ba2+ ion from Li et al. JCTC, 2013, 9, 2733'),
        }
    }

    # IOD Parameters for trivalent ions
    # From Li et al. JPCB, 2015, 119, 883
    IOD_TRIVALENT = {
        'tip3p': {
            'Al3': (1.297, 0.00471279, 'IOD set for Al3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Fe3': (1.386, 0.01357097, 'IOD set for Fe3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Cr3': (1.344, 0.00848000, 'IOD set for Cr3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'In3': (1.461, 0.02808726, 'IOD set for In3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Tl3': (1.513, 0.04321029, 'IOD set for Tl3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Y3': (1.602, 0.08034231, 'IOD set for Y3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'La3': (1.718, 0.15060822, 'IOD set for La3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Ce3': (1.741, 0.16721338, 'IOD set for Ce3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Pr3': (1.733, 0.16134811, 'IOD set for Pr3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Nd3': (1.681, 0.12564307, 'IOD set for Nd3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Sm3': (1.659, 0.11189491, 'IOD set for Sm3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Eu3': (1.666, 0.11617738, 'IOD set for Eu3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Gd3': (1.623, 0.09126804, 'IOD set for Gd3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Tb3': (1.630, 0.09509276, 'IOD set for Tb3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Dy3': (1.609, 0.08389240, 'IOD set for Dy3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Er3': (1.602, 0.08034231, 'IOD set for Er3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Tm3': (1.602, 0.08034231, 'IOD set for Tm3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Lu3': (1.588, 0.07351892, 'IOD set for Lu3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
        }
    }

    # IOD Parameters for tetravalent ions
    # From Li et al. JPCB, 2015, 119, 883
    IOD_TETRAVALENT = {
        'tip3p': {
            'Hf4': (1.499, 0.03868661, 'IOD set for Hf4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Zr4': (1.519, 0.04525501, 'IOD set for Zr4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Ce4': (1.684, 0.12758274, 'IOD set for Ce4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'U4': (1.684, 0.12758274, 'IOD set for U4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Pu4': (1.662, 0.11371963, 'IOD set for Pu4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
            'Th4': (1.708, 0.14364160, 'IOD set for Th4+ ion TIP3P water from Li et al. JPCB, 2015, 119, 883'),
        }
    }

    # Map element+charge to IOD parameter key
    ELEMENT_TO_IOD = {
        # Monovalent
        ('Li', 1): 'Li1',
        ('Na', 1): 'Na1',
        ('K', 1): 'K1',
        ('Rb', 1): 'Rb1',
        ('Cs', 1): 'Cs1',
        ('Tl', 1): 'Tl1',
        ('Cu', 1): 'Cu1',
        ('Ag', 1): 'Ag1',
        ('F', -1): 'F-1',
        ('Cl', -1): 'Cl-1',
        ('Br', -1): 'Br-1',
        ('I', -1): 'I-1',
        # Divalent
        ('Be', 2): 'Be2',
        ('Cu', 2): 'Cu2',
        ('Ni', 2): 'Ni2',
        ('Zn', 2): 'Zn2',
        ('Co', 2): 'Co2',
        ('Cr', 2): 'Cr2',
        ('Fe', 2): 'Fe2',
        ('Mg', 2): 'Mg2',
        ('V', 2): 'V2',
        ('Mn', 2): 'Mn2',
        ('Hg', 2): 'Hg2',
        ('Cd', 2): 'Cd2',
        ('Ca', 2): 'Ca2',
        ('Sn', 2): 'Sn2',
        ('Sr', 2): 'Sr2',
        ('Ba', 2): 'Ba2',
        # Trivalent
        ('Al', 3): 'Al3',
        ('Fe', 3): 'Fe3',
        ('Cr', 3): 'Cr3',
        ('In', 3): 'In3',
        ('Tl', 3): 'Tl3',
        ('Y', 3): 'Y3',
        ('La', 3): 'La3',
        ('Ce', 3): 'Ce3',
        ('Pr', 3): 'Pr3',
        ('Nd', 3): 'Nd3',
        ('Sm', 3): 'Sm3',
        ('Eu', 3): 'Eu3',
        ('Gd', 3): 'Gd3',
        ('Tb', 3): 'Tb3',
        ('Dy', 3): 'Dy3',
        ('Er', 3): 'Er3',
        ('Tm', 3): 'Tm3',
        ('Lu', 3): 'Lu3',
        # Tetravalent
        ('Hf', 4): 'Hf4',
        ('Zr', 4): 'Zr4',
        ('Ce', 4): 'Ce4',
        ('U', 4): 'U4',
        ('Pu', 4): 'Pu4',
        ('Th', 4): 'Th4',
    }

    def __init__(self, water_model: str = 'tip3p', logger=None):
        """
        Initialize metal ion database.

        Args:
            water_model: Water model for VDW parameters (default: tip3p)
            logger: Optional logger instance
        """
        self.water_model = water_model.lower()
        self.logger = logger or logging.getLogger(__name__)

        # Validate water model
        if self.water_model not in SUPPORTED_WATER_MODELS:
            self.logger.warning(
                f"Water model '{water_model}' not fully supported. "
                f"Falling back to tip3p parameters."
            )
            self.water_model = 'tip3p'

    def is_metal(self, resname: str, atom_name: str) -> bool:
        """
        Check if (resname, atomname) corresponds to a known metal.

        Args:
            resname: PDB residue name
            atom_name: PDB atom name

        Returns:
            True if this is a recognized metal
        """
        key = (resname.upper(), atom_name.upper())
        return key in self.METAL_PDB

    def is_metal_element(self, element: str) -> bool:
        """
        Check if an element symbol corresponds to a metal.

        This is useful for identifying metal atoms in complex residues
        like iron-sulfur clusters (SF4) where the resname doesn't match
        the standard metal residue names.

        Args:
            element: Element symbol (e.g., 'Fe', 'Zn', 'Cu')

        Returns:
            True if this element is a recognized metal
        """
        element_upper = element.upper() if element else ""
        # Check all METAL_PDB entries for this element
        for (res, atom), (elem, _, _) in self.METAL_PDB.items():
            if elem.upper() == element_upper:
                return True
        return False

    def get_metal_config(self, resname: str, atom_name: str,
                        charge: int, spin: int) -> Optional[MetalConfig]:
        """
        Get complete metal configuration.

        Args:
            resname: PDB residue name (e.g., 'FE', 'FE2')
            atom_name: PDB atom name (e.g., 'FE', 'FE2')
            charge: Formal charge (REQUIRED - must be provided by user)
            spin: Number of unpaired electrons (REQUIRED - must be provided by user)
                  Positive for alpha spins, negative for beta spins

        Returns:
            MetalConfig object or None if not found
        """
        key = (resname.upper(), atom_name.upper())

        if key not in self.METAL_PDB:
            return None

        element, atomic_num, mass = self.METAL_PDB[key]

        # Get VDW parameters using user-provided charge
        vdw_radius, vdw_epsilon, vdw_source = self._get_vdw_params(element, charge)

        return MetalConfig(
            resname=resname.upper(),
            atom_name=atom_name.upper(),
            element=element,
            atomic_number=atomic_num,
            mass=mass,
            atom_type=resname.upper(),  # For metals, atom type = resname
            charge=charge,  # Store user-provided charge
            spin=spin,  # Store user-provided spin
            vdw_radius=vdw_radius,
            vdw_epsilon=vdw_epsilon,
            vdw_source=vdw_source,
            water_model=self.water_model
        )

    def get_metal_config_by_element(self, element: str, atom_name: str,
                                    resname: str, charge: int, spin: int) -> Optional[MetalConfig]:
        """
        Get metal configuration by element symbol.

        This is useful for metal atoms in complex residues (like SF4 iron-sulfur clusters)
        where the resname doesn't match standard metal residue names.

        Args:
            element: Element symbol (e.g., 'Fe')
            atom_name: Original PDB atom name (e.g., 'FE1')
            resname: Original PDB residue name (e.g., 'SF4')
            charge: Formal charge (REQUIRED - must be provided by user)
            spin: Number of unpaired electrons (REQUIRED - must be provided by user)
                  Positive for alpha spins, negative for beta spins

        Returns:
            MetalConfig object with generic resname/atomtype, or None if element not found
        """
        element_upper = element.upper() if element else ""

        # Find any entry with this element to get atomic properties
        for (res, atom), (elem, atomic_num, mass) in self.METAL_PDB.items():
            if elem.upper() == element_upper:
                # Get VDW parameters using user-provided charge
                vdw_radius, vdw_epsilon, vdw_source = self._get_vdw_params(elem, charge)

                # Use element symbol as atom type (e.g., "FE" for all iron atoms)
                # This ensures consistency for iron-sulfur clusters
                atom_type = elem.upper()

                return MetalConfig(
                    resname=resname.upper(),  # Keep original resname (e.g., "SF4")
                    atom_name=atom_name.upper(),  # Keep original atom name (e.g., "FE1")
                    element=elem,
                    atomic_number=atomic_num,
                    mass=mass,
                    atom_type=atom_type,  # Use element as atom type (e.g., "FE")
                    charge=charge,  # Store user-provided charge
                    spin=spin,  # Store user-provided spin
                    vdw_radius=vdw_radius,
                    vdw_epsilon=vdw_epsilon,
                    vdw_source=vdw_source,
                    water_model=self.water_model
                )

        return None

    def _get_vdw_params(self, element: str, charge: int) -> Tuple[float, float, str]:
        """
        Get VDW parameters (radius, epsilon, source) for metal.

        Args:
            element: Element symbol
            charge: Formal charge (user-provided)

        Returns:
            (radius, epsilon, source) tuple
        """
        # Look up IOD parameters using user-provided charge
        iod_key = self.ELEMENT_TO_IOD.get((element, charge))

        if iod_key:
            # Try monovalent first
            if abs(charge) == 1:
                params = self.IOD_MONOVALENT.get(self.water_model, {})
                if iod_key in params:
                    return params[iod_key]

            # Try divalent
            elif abs(charge) == 2:
                params = self.IOD_DIVALENT.get(self.water_model, {})
                if iod_key in params:
                    return params[iod_key]

            # Try trivalent
            elif abs(charge) == 3:
                params = self.IOD_TRIVALENT.get(self.water_model, {})
                if iod_key in params:
                    return params[iod_key]

            # Try tetravalent
            elif abs(charge) == 4:
                params = self.IOD_TETRAVALENT.get(self.water_model, {})
                if iod_key in params:
                    return params[iod_key]

        # The IOD tables above cover ions with a measured ion-oxygen distance,
        # which excludes metals that do not exist as free aqueous ions. MCPB.py
        # carries a wider table that fills those in from UFF -- Mo(VI), for one,
        # is there as 'Mo6' and would otherwise reach the placeholder below.
        mcpb = self._mcpb_lj_params(element, charge)
        if mcpb:
            return mcpb

        # No parameters available for this element/charge combination
        self.logger.warning(
            f"No IOD parameters found for {element} with charge {charge:+d}. "
            f"Using generic fallback values (VERIFY MANUALLY)."
        )
        return (1.5, 0.01, f'Generic fallback for {element}{charge:+d} (VERIFY MANUALLY)')

    # MCPB.py's LJ table, per water model. Keyed by water model because the
    # IOD entries in it differ between models; the UFF-derived ones do not.
    _LJ_TABLES: Dict[str, Dict[str, tuple]] = {}

    @classmethod
    def _mcpb_lj_table(cls, water_model: str) -> Dict[str, tuple]:
        """MCPB.py's IonLJParaDict for one water model, or {} without pymsmt."""
        if water_model not in cls._LJ_TABLES:
            try:
                from pymsmt.mol.element import get_ionljparadict
                # Table keys are padded to a fixed width ('K1 ', 'Mg2').
                cls._LJ_TABLES[water_model] = {
                    k.strip(): v for k, v in get_ionljparadict(water_model).items()
                }
            except Exception:  # noqa: BLE001 — AmberTools may not be installed
                cls._LJ_TABLES[water_model] = {}
        return cls._LJ_TABLES[water_model]

    def _mcpb_lj_params(self, element: str, charge: int):
        """(radius, epsilon, source) from MCPB.py's table, or None.

        Keys are the element symbol followed by the charge magnitude, with a
        leading minus for anions: 'Fe3', 'Mo6', 'Cl-1'.
        """
        if not element:
            return None
        symbol = element.strip().capitalize()
        key = f"{symbol}-{abs(charge)}" if charge < 0 else f"{symbol}{charge}"
        return self._mcpb_lj_table(self.water_model).get(key)

    # MCPB.py's own MK radius table, loaded once. Class-level so the import is
    # attempted a single time per process.
    _MK_RADII: Optional[Dict[str, tuple]] = None

    @classmethod
    def _mcpb_mk_radii(cls) -> Dict[str, tuple]:
        """MCPB.py's vdwRadiiDict2023, or {} if pymsmt is not importable."""
        if cls._MK_RADII is None:
            try:
                from pymsmt.mol.element import vdwRadiiDict2023
                cls._MK_RADII = dict(vdwRadiiDict2023)
            except Exception:  # noqa: BLE001 — AmberTools may not be installed
                cls._MK_RADII = {}
        return cls._MK_RADII

    def get_vdw_radius(self, element: str, charge: int) -> Optional[float]:
        """
        Radius for Gaussian's Merz-Kollman ESP fit, Pop(MK,ReadRadii).

        This is NOT the force-field Lennard-Jones radius. It is the sampling
        radius that decides how close to the nucleus ESP grid points may fall,
        so it is looked up separately from the IOD parameters that feed the
        frcmod NONBON section (``_get_vdw_params``, unchanged).

        Preferred source is MCPB.py's own ``vdwRadiiDict2023`` — the table
        MCPB.py writes into its ReadRadii block, from Smith et al. JCTC 2023,
        19, 2064, with 267 entries carrying their provenance. The IOD tables
        used for force-field parameters only run to tetravalent, so a Mo(VI)
        cofactor could never resolve there and fell to a generic 1.5 A with no
        cited source.

        Args:
            element: Element symbol (e.g., 'Zn', 'Cu', 'Fe')
            charge: Formal charge (e.g., 2 for Zn2+)

        Returns:
            Radius in Angstroms, or None if neither source has one.
        """
        entry = self._mcpb_mk_radii().get(f"{element}{charge}")
        if entry:
            try:
                return float(entry[0])
            except (TypeError, ValueError, IndexError):
                pass

        radius, _, _ = self._get_vdw_params(element, charge)
        return radius

    def get_all_metals(self) -> Dict[Tuple[str, str], Tuple[str, int, float]]:
        """Get all known metal definitions."""
        return self.METAL_PDB.copy()

    def __repr__(self):
        return f"MetalIonDatabase(water_model='{self.water_model}', n_metals={len(self.METAL_PDB)})"
