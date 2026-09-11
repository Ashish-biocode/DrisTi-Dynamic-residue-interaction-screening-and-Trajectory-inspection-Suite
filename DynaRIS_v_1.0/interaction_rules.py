import numpy as np
from rdkit import Chem
import tempfile
import os
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

# ==========================================
# IMPORT LIGAND PROCESSOR
# ==========================================
from chemistry.ligand_processor import (
    create_rdkit_mol_from_ligand,
    generate_ligand_smiles
)


# ==========================================
# BASIC GEOMETRY
# ==========================================
def atom_distance(a, b):

    return np.linalg.norm(
        np.array(a.position) -
        np.array(b.position)
    )


def centroid(atoms):

    coords = np.array(
        [a.position for a in atoms]
    )

    return coords.mean(axis=0)


def ring_normal(atoms):

    coords = np.array(
        [a.position for a in atoms]
    )

    if len(coords) < 3:
        return np.array([0.0, 0.0, 1.0])

    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]

    normal = np.cross(v1, v2)

    norm = np.linalg.norm(normal)

    if norm == 0:
        return np.array([0.0, 0.0, 1.0])

    return normal / norm


def angle_between(v1, v2):

    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)

    if norm1 == 0 or norm2 == 0:
        return 0.0

    v1 = v1 / norm1
    v2 = v2 / norm2

    dot = np.clip(
        np.dot(v1, v2),
        -1.0,
        1.0
    )

    return np.degrees(
        np.arccos(dot)
    )


def centroid_distance(group1, group2):

    c1 = centroid(group1)
    c2 = centroid(group2)

    return np.linalg.norm(c1 - c2)


def closest_atom_distance(
        atoms1,
        atoms2):

    min_dist = 999.0

    for a1 in atoms1:
        for a2 in atoms2:

            d = atom_distance(a1, a2)

            if d < min_dist:
                min_dist = d

    return min_dist


# ==========================================
# ELEMENT UTILITIES
# ==========================================
def get_element(atom):

    if hasattr(atom, "element") and atom.element:

        return atom.element.strip().upper()

    name = atom.name.strip()

    two_letter = name[:2].upper()

    if two_letter in [
        "CL", "BR", "NA",
        "MG", "ZN", "FE",
        "CA", "CU", "MN"
    ]:
        return two_letter

    return name[0].upper()


def is_heavy(atom):

    return get_element(atom) != "H"


def filter_heavy(atoms):

    return [
        a for a in atoms
        if is_heavy(a)
    ]


# ==========================================
# RESIDUE DEFINITIONS
# ==========================================
HYDROPHOBIC_RESIDUES = [

    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PRO"
]

ALKYL_RESIDUES = [

    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PRO",
    "LYS"
]

AROMATIC_RESIDUES = [

    "PHE",
    "TYR",
    "TRP",
    "HIS"
]

NEGATIVE_RESIDUES = [

    "ASP",
    "GLU"
]

POSITIVE_RESIDUES = [

    "LYS",
    "ARG",
    "HIS"
]

BACKBONE_ATOMS = [

    "N",
    "CA",
    "C",
    "O"
]

NEGATIVE_SIDECHAIN_ATOMS = [

    "OD1",
    "OD2",
    "OE1",
    "OE2"
]

POSITIVE_SIDECHAIN_ATOMS = [

    "NZ",
    "NE",
    "NH1",
    "NH2"
]

AROMATIC_RING_ATOMS = {

    "PHE": [

        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ"
    ],

    "TYR": [

        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ"
    ],

    "HIS": [

        "CG",
        "ND1",
        "CD2",
        "CE1",
        "NE2"
    ],

    "TRP": [

        "CD2",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2"
    ]
}


# ==========================================
# H-BOND DEFINITIONS
# ==========================================
def get_ligand_donor_acceptor_atoms(
        ligand_atoms,
        ligand_sdf_file,
        sdf_molecule_index=0):

    donor_indices = set()
    acceptor_indices = set()

    mol = create_rdkit_mol_from_ligand(
        ligand_sdf_file,
        molecule_index=sdf_molecule_index
    )

    if mol is None:
        return donor_indices, acceptor_indices

    fdef = os.path.join(
        RDConfig.RDDataDir,
        "BaseFeatures.fdef"
    )

    factory = ChemicalFeatures.BuildFeatureFactory(
        fdef
    )

    feats = factory.GetFeaturesForMol(
        mol
    )

    for feat in feats:

        fam = feat.GetFamily()

        atom_ids = feat.GetAtomIds()

        atom_names = []

        for idx in atom_ids:

            if idx < len(ligand_atoms):

                atom_names.append(
                    ligand_atoms[idx].name
                )

        pass

        if fam == "Donor":

            donor_indices.update(
                atom_ids
            )

        elif fam == "Acceptor":

            acceptor_indices.update(
                atom_ids
            )

    return donor_indices, acceptor_indices


HBOND_DONOR_ATOMS = {

    "ARG": ["NE", "NH1", "NH2"],
    "ASN": ["ND2"],
    "GLN": ["NE2"],
    "HIS": [],
    "LYS": ["NZ"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TRP": ["NE1"],
    "TYR": ["OH"],

    # backbone donor
    "BACKBONE": ["N"]
}


HBOND_ACCEPTOR_ATOMS = {

    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "ASN": ["OD1"],
    "GLN": ["OE1"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
    "HIS": ["NE2"],

    # backbone acceptor
    "BACKBONE": ["O"]
}


def is_hbond_donor(atom):

    if atom.name.startswith("N"):
        return True

    atoms = HBOND_DONOR_ATOMS.get(
        atom.resname,
        []
    )

    return atom.name in atoms


def is_hbond_acceptor(atom):

    if atom.name.startswith("O"):
        return True

    atoms = HBOND_ACCEPTOR_ATOMS.get(
        atom.resname,
        []
    )

    return atom.name in atoms

# ==========================================
# H-BOND GEOMETRY
# ==========================================
def get_attached_hydrogens(atom):

    hydrogens = []

    try:
        for bonded in atom.bonded_atoms:

            if get_element(bonded) == "H":
                hydrogens.append(bonded)

    except Exception:
        pass

    return hydrogens


# ==========================================
# CARBON H-BOND DONOR CHECK
# ==========================================
def is_carbon_hbond_donor(atom):

    # Must be carbon
    if get_element(atom) != "C":
        return False

    # Must have at least one attached hydrogen
    if len(get_attached_hydrogens(atom)) == 0:
        return False

    try:

        rd_atom = atom.rdkit_atom

    except AttributeError:
        # RDKit atom unavailable
        return True

    # Aromatic carbons are good donors
    if rd_atom.GetIsAromatic():
        return True

    # sp2 carbons are preferred donors
    if str(rd_atom.GetHybridization()) == "SP2":
        return True

    # Carbon attached to O/N also becomes a better donor
    for nbr in rd_atom.GetNeighbors():

        if nbr.GetSymbol() in ["O", "N"]:
            return True

    # Ordinary alkyl carbon
    return False

def estimate_hbond_angle(
        protein_atom,
        ligand_atom,
        protein_is_donor=False,
        ligand_is_donor=False):

    donor_atom = None
    acceptor_atom = None

    if protein_is_donor:
        donor_atom = protein_atom
        acceptor_atom = ligand_atom

    elif ligand_is_donor:
        donor_atom = ligand_atom
        acceptor_atom = protein_atom

    else:
        return 180.0

    hydrogens = get_attached_hydrogens(donor_atom)

    if not hydrogens:
        # Crystal structures often lack explicit hydrogens.
        return 120.0

    best_angle = 0.0

    for h in hydrogens:

        v1 = donor_atom.position - h.position
        v2 = acceptor_atom.position - h.position

        angle = angle_between(v1, v2)

        if angle > best_angle:
            best_angle = angle

    return best_angle

# ==========================================
# LIGAND CHARGED ATOMS
# ==========================================
def get_ligand_charged_atoms(
        ligand_file,
        sdf_molecule_index=0):

    positive_atoms = set()
    negative_atoms = set()

    mol = create_rdkit_mol_from_ligand(
        ligand_file,
        molecule_index=sdf_molecule_index
    )

    if mol is None:
        return positive_atoms, negative_atoms

    for atom in mol.GetAtoms():

        charge = atom.GetFormalCharge()

        if charge > 0:
            positive_atoms.add(atom.GetIdx())

        elif charge < 0:
            negative_atoms.add(atom.GetIdx())

    return positive_atoms, negative_atoms


# ==========================================
# H-BOND DETECTION
# ==========================================
def detect_hbond(
        patom,
        latom,
        ligand_index,
        ligand_donors,
        ligand_acceptors):

    interactions = []

    pe = get_element(patom)
    le = get_element(latom)

    if pe == "H" or le == "H":
        return interactions

    protein_donor = is_hbond_donor(
    patom
    )

    protein_acceptor = is_hbond_acceptor(
    patom
    )

    ligand_donor = (
    ligand_index
    in ligand_donors
    )

    ligand_acceptor = (
    ligand_index
    in ligand_acceptors
    )



    valid_pair = (
    (
        protein_donor
        and
        ligand_acceptor
    )
    or
    (
        ligand_donor
        and
        protein_acceptor
    )
)

    if not valid_pair:
        return interactions

    dist = atom_distance(
        patom,
        latom
    )


    if dist > 3.8:
        return interactions

    pass

    angle = estimate_hbond_angle(
        patom,
        latom,
        protein_is_donor=protein_donor,
        ligand_is_donor=ligand_donor
    )

    # Reject H-bonds whose angle could not
    # be determined from real hydrogens

    hydrogens = []

    if ligand_donor:
        hydrogens = get_attached_hydrogens(latom)

    elif protein_donor:
        hydrogens = get_attached_hydrogens(patom)

    # If hydrogens are missing, skip this interaction.
    # A future version can estimate hydrogen positions.
    if len(hydrogens) == 0:

        return interactions
    
    pass



    if angle < 100:
        return interactions

    if dist <= 3.0 and angle >= 150:

        interaction_type = "H-Bond (Strong)"

    elif dist <= 3.3 and angle >= 130:

        interaction_type = "H-Bond (Moderate)"

    else:

        interaction_type = "H-Bond (Weak)"
    


    pass

    interactions.append({

        "protein_atom":
            patom,

        "ligand_atom":
            latom,

        "type":
            interaction_type,

        "distance":
            float(dist)
    })

    return interactions

# ==========================================
# SALT BRIDGE DETECTION
# ==========================================
def detect_salt_bridge(
        patom,
        latom,
        ligand_index,
        ligand_positive,
        ligand_negative):

    interactions = []

    dist = atom_distance(
        patom,
        latom
    )

    if dist > 4.0:
        return interactions

    # -------------------------------
    # Protein positive
    # -------------------------------
    protein_positive = (

        patom.resname in POSITIVE_RESIDUES

        and

        patom.name in POSITIVE_SIDECHAIN_ATOMS
    )

    # -------------------------------
    # Protein negative
    # -------------------------------
    protein_negative = (

        patom.resname in NEGATIVE_RESIDUES

        and

        patom.name in NEGATIVE_SIDECHAIN_ATOMS
    )

    # -------------------------------
    # Ligand charge
    # -------------------------------
    ligand_is_positive = (
        ligand_index in ligand_positive
    )

    ligand_is_negative = (
        ligand_index in ligand_negative
    )

    # -------------------------------
    # Opposite charges required
    # -------------------------------
    valid_pair = (

        (
            protein_positive
            and
            ligand_is_negative
        )

        or

        (
            protein_negative
            and
            ligand_is_positive
        )
    )

    if not valid_pair:
        return interactions

    pass

    interactions.append({

        "protein_atom":
            patom,

        "ligand_atom":
            latom,

        "type":
            "Salt Bridge",

        "distance":
            float(dist)
    })

    return interactions


# ==========================================
# CARBON H-BOND DETECTION
# ==========================================
def detect_carbon_hbond(
        patom,
        latom):

    interactions = []

    protein_acceptor = is_hbond_acceptor(
        patom
    )

    if not protein_acceptor:
        return interactions

    if get_element(latom) != "C":
        return interactions


    if get_element(latom) == "C":
        rd = latom.rdkit_atom
        pass

    if not is_carbon_hbond_donor(latom):
        return interactions

    # Find hydrogens attached to ligand carbon
    hydrogens = get_attached_hydrogens(latom)

    if len(hydrogens) == 0:
        return interactions

    best_h_distance = None
    best_angle = None
    best_hydrogen = None

    for hydrogen in hydrogens:

        h_dist = atom_distance(
            hydrogen,
            patom
        )

        v1 = (
            latom.position
            - hydrogen.position
        )

        v2 = (
            patom.position
            - hydrogen.position
        )

        angle = angle_between(
            v1,
            v2
        )

        if (
            best_h_distance is None
            or
            h_dist < best_h_distance
        ):

            best_h_distance = h_dist
            best_angle = angle
            best_hydrogen = hydrogen

    # Carbon H-bond cutoff (H···acceptor)
    if best_h_distance > 2.8:
        return interactions

    if best_angle < 110:
        return interactions

    # Heavy-atom distance (for reporting only)
    dist = atom_distance(
        patom,
        latom
    )

    interactions.append({

        "protein_atom":
            patom,

        "ligand_atom":
            latom,

        "type":
            "Carbon H-Bond",

        "distance":
            float(dist)
    })

    return interactions


# ==========================================
# HYDROPHOBIC DETECTION
# ==========================================
def is_hydrophobic(patom, latom):

    if get_element(latom) != "C":
        return False

    if not is_hydrophobic_sidechain_carbon(patom):
        return False

    if patom.resname not in HYDROPHOBIC_RESIDUES:
        return False

    return True


# ==========================================
# HYDROPHOBIC HELPERS
# ==========================================
def is_hydrophobic_sidechain_carbon(atom):

    if get_element(atom) != "C":
        return False

    if atom.name in BACKBONE_ATOMS:
        return False

    return True


# ==========================================
# PROTEIN AROMATIC RINGS
# ==========================================
def get_protein_rings(protein_atoms):

    rings = []

    residues = {}

    for atom in protein_atoms:

        key = (
            atom.resname,
            atom.resid
        )

        residues.setdefault(
            key,
            []
        ).append(atom)

    for (resname, resid), atoms in residues.items():

        if resname not in AROMATIC_RESIDUES:
            continue

        needed = AROMATIC_RING_ATOMS.get(
            resname,
            []
        )

        ring_atoms = [

            a for a in atoms
            if a.name in needed
        ]

        if len(ring_atoms) < 5:
            continue

        rings.append({

            "resname":
                resname,

            "resid":
                resid,

            "atoms":
                ring_atoms,

            "centroid":
                centroid(ring_atoms),

            "normal":
                ring_normal(ring_atoms)
        })

    return rings


# ==========================================
# LIGAND AROMATIC RINGS
# ==========================================
def get_ligand_aromatic_rings(
        ligand_atoms,
        ligand_pdb_file=None,
        sdf_molecule_index=0):

    rings = []

    # ======================================
    # CREATE RDKit MOLECULE
    # ======================================
    mol = create_rdkit_mol_from_ligand(
        ligand_pdb_file,
        molecule_index=sdf_molecule_index
    )

    if mol is None:

        pass

        return rings

    # ======================================
    # AUTO GENERATE SMILES
    # ======================================
    try:

        smiles = generate_ligand_smiles(mol)

        pass

    except Exception:

        smiles = None

    try:
        Chem.SanitizeMol(mol)

    except Exception:
        pass

    ring_info = mol.GetRingInfo()

    atom_rings = ring_info.AtomRings()

    for ring in atom_rings:

        if len(ring) not in [5, 6]:
            continue

        rd_ring_atoms = [
            mol.GetAtomWithIdx(i)
            for i in ring
        ]

        if not all(
            a.GetIsAromatic()
            for a in rd_ring_atoms
        ):
            continue

        ring_atoms = []

        for idx in ring:

            if idx >= len(ligand_atoms):
                continue

            ring_atoms.append(
                ligand_atoms[idx]
            )

        coords = np.array(
            [a.position for a in ring_atoms]
        )

        centroid_coords = coords.mean(axis=0)

        shifted = coords - centroid_coords

        try:
            _, _, vh = np.linalg.svd(shifted)

        except Exception:
            continue

        normal = vh[2]

        distances = np.dot(
            shifted,
            normal
        )

        max_dev = np.max(
            np.abs(distances)
        )

        if max_dev > 0.30:
            continue

        rings.append({

            "atoms":
                ring_atoms,

            "centroid":
                centroid(ring_atoms),

            "normal":
                ring_normal(ring_atoms)
        })

    return rings


# ==========================================
# PI-PI DETECTION
# ==========================================
def detect_pi_pi(
        protein_rings,
        ligand_rings):

    interactions = []

    seen = set()

    for pring in protein_rings:

        for lring in ligand_rings:

            center_dist = np.linalg.norm(
                pring["centroid"]
                -
                lring["centroid"]
            )

            centroid_vector = (
                lring["centroid"]
                -
                pring["centroid"]
            )

            vertical = abs(
                np.dot(
                    centroid_vector,
                    pring["normal"]
                )
            )
            horizontal = np.sqrt(
                max(
                    center_dist**2 - vertical**2,
                    0.0
                )
            )

            if horizontal > 5.0:
                continue

            if center_dist > 6.0:
                continue

            closest_dist = closest_atom_distance(
                pring["atoms"],
                lring["atoms"]
            )

            if closest_dist > 4.5:
                continue

            dot = np.dot(
                pring["normal"],
                lring["normal"]
            )

            dot = np.clip(
                np.abs(dot),
                -1.0,
                1.0
            )

            angle = np.degrees(
                np.arccos(dot)
            )

            interaction_type = None

            if angle <= 35:
                interaction_type = "Pi-Pi Stacked"

            elif angle >= 55:
                interaction_type = "Pi-Pi T-Shaped"

            else:
                continue

            key = (
                pring["resname"],
                pring["resid"],
                interaction_type
            )

            if key in seen:
                continue

            seen.add(key)

            pass

            interactions.append({

                "protein_atom":
                    pring["atoms"][0],

                "ligand_atom":
                    lring["atoms"][0],

                "type":
                    interaction_type,

                "distance":
                    float(center_dist)
            })

    return interactions


# ==========================================
# PI-ANION DETECTION
# ==========================================
def detect_pi_anion(
        protein_atoms,
        ligand_rings):

    interactions = []

    for lring in ligand_rings:

        ligand_centroid = lring["centroid"]
        ligand_normal = lring["normal"]

        best_interaction = None

        for patom in protein_atoms:

            if patom.resname not in NEGATIVE_RESIDUES:
                continue

            if patom.name not in NEGATIVE_SIDECHAIN_ATOMS:
                continue

            vec = (
                np.array(patom.position)
                -
                ligand_centroid
            )

            dist = np.linalg.norm(vec)

            if dist > 4.2:
                continue

            angle = angle_between(
                vec,
                ligand_normal
            )

            deviation = min(
                angle,
                180 - angle
            )

            if deviation > 40:
                continue

            if (
                best_interaction is None
                or
                dist < best_interaction["distance"]
            ):

                best_interaction = {

                    "protein_atom":
                        patom,

                    "ligand_atom":
                        lring["atoms"][0],

                    "type":
                        "Pi-Anion",

                    "distance":
                        float(dist)
                }

        if best_interaction is not None:
            interactions.append(
                best_interaction
            )

    return interactions


# ==========================================
# PI-CATION DETECTION
# ==========================================
def detect_pi_cation(
        protein_atoms,
        ligand_rings):

    interactions = []

    seen = set()

    for lring in ligand_rings:

        ligand_centroid = lring["centroid"]
        ligand_normal = lring["normal"]

        for patom in protein_atoms:

            if patom.resname not in POSITIVE_RESIDUES:
                continue

            if patom.name not in POSITIVE_SIDECHAIN_ATOMS:
                continue

            vec = (
                np.array(patom.position)
                -
                ligand_centroid
            )

            dist = np.linalg.norm(vec)

            if dist > 6.0:
                continue

            angle = angle_between(
                vec,
                ligand_normal
            )

            deviation = min(
                angle,
                180 - angle
            )

            if deviation > 40:
                continue

            key = (
                patom.resname,
                patom.resid
            )

            if key in seen:
                continue

            seen.add(key)

            interactions.append({

                "protein_atom":
                    patom,

                "ligand_atom":
                    lring["atoms"][0],

                "type":
                    "Pi-Cation",

                "distance":
                    float(dist)
            })

    return interactions

# ==========================================
# PI-ALKYL / HYdrophobic DETECTION
# ==========================================
def detect_pi_alkyl(
        protein_atoms,
        ligand_rings):

    interactions = []

    best_pi_alkyl = {}

    for lring in ligand_rings:

        ligand_centroid = np.array(
            lring["centroid"]
        )

        ring_normal = np.array(
            lring["normal"]
        )

        for patom in protein_atoms:

            if patom.resname not in ALKYL_RESIDUES:
                continue

            if get_element(patom) != "C":
                continue

            if patom.name in BACKBONE_ATOMS:
                continue

            v = (
                np.array(patom.position)
                - ligand_centroid
            )

            dist = np.linalg.norm(v)

            # NEW
            closest = closest_atom_distance(
                [patom],
                lring["atoms"]
            )

            # Reject only if BOTH criteria fail
            if dist > 5.5 and closest > 4.5:
                continue

            vertical = abs(
                np.dot(v, ring_normal)
            )

            horizontal = np.sqrt(
                max(
                    0.0,
                    dist * dist
                    - vertical * vertical
                )
            )


            pass
            # Lateral offset filter
            if horizontal > 3.3 and closest > 4.5:
                continue

            interaction_type = "Pi-Alkyl"


            key = (
                patom.resname,
                patom.resid
            )

            candidate = {

                "protein_atom": patom,

                "ligand_atom": lring["atoms"][0],

                "type": interaction_type,

                "distance": float(dist)
            }

            if (
                key not in best_pi_alkyl
                or
                dist < best_pi_alkyl[key]["distance"]
            ):
                best_pi_alkyl[key] = candidate

    interactions.extend(best_pi_alkyl.values())

    return interactions


# ==========================================
# AMIDE-PI DETECTION
# ==========================================
def detect_amide_pi(
        protein_atoms,
        ligand_rings):

    interactions = []
    seen = set()

    # build lookup
    residue_atoms = {}

    for atom in protein_atoms:

        residue_atoms.setdefault(
            atom.resid,
            []
        ).append(atom)

    for resid in sorted(residue_atoms):

        current_residue = residue_atoms[resid]
        previous_residue = residue_atoms.get(
            resid - 1,
            None
        )

        if previous_residue is None:
            continue

        n_atom = next(
            (
                a for a in current_residue
                if a.name == "N"
            ),
            None
        )

        c_atom = next(
            (
                a for a in previous_residue
                if a.name == "C"
            ),
            None
        )

        o_atom = next(
            (
                a for a in previous_residue
                if a.name == "O"
            ),
            None
        )

        if (
            n_atom is None
            or c_atom is None
            or o_atom is None
        ):
            continue

        amide_centroid = (
            np.array(c_atom.position)
            + np.array(o_atom.position)
            + np.array(n_atom.position)
        ) / 3.0

        # Amide plane normal
        v1 = np.array(o_atom.position) - np.array(c_atom.position)
        v2 = np.array(n_atom.position) - np.array(c_atom.position)

        amide_normal = np.cross(v1, v2)

        norm = np.linalg.norm(amide_normal)

        if norm < 1e-6:
            continue

        amide_normal = amide_normal / norm
        
        for lring in ligand_rings:

            ring_centroid = lring["centroid"]
            ring_normal = lring["normal"]


            offset_vec = amide_centroid - ring_centroid

            dist = np.linalg.norm(offset_vec)
            
            if dist > 4.5:
                continue

            # Plane angle
            angle = np.degrees(
                np.arccos(
                    np.clip(
                        abs(np.dot(
                            amide_normal,
                            ring_normal
                        )),
                        -1.0,
                        1.0
                    )
                )
            )

            # Vertical distance above ring plane
            vertical = abs(
                np.dot(
                    offset_vec,
                    ring_normal
                )
            )

            # Horizontal displacement
            horizontal = np.sqrt(
                max(
                    dist**2 - vertical**2,
                    0.0
                )
            )

            if horizontal > 3.5:
                continue

            key = (
                c_atom.resname,
                c_atom.resid,
                n_atom.resname,
                n_atom.resid
            )
            interactions.append({

                "protein_atom": c_atom,

                "ligand_atom": lring["atoms"][0],

                "type": "Amide-Pi Stacked",

                "distance": float(dist),

                "amide_prev_res": f"{c_atom.resname}-{c_atom.resid}",

                "amide_curr_res": f"{n_atom.resname}-{n_atom.resid}",

                "amide_atoms": "C,O,N"
            })

    return interactions

# ==========================================
# MASTER INTERACTION ENGINE
# ==========================================
def detect_interactions(
        protein_atoms,
        ligand_atoms,
        ligand_pdb_file=None,
        sdf_molecule_index=0):

    interactions = []

    protein_all_atoms = protein_atoms
    ligand_all_atoms = ligand_atoms

    ligand_hydrogens = [
    	a for a in ligand_all_atoms
    	if get_element(a) == "H"
    ]

    pass

    protein_atoms = filter_heavy(
        protein_atoms
    )

    ligand_atoms = filter_heavy(
        ligand_atoms
    )

    # ======================================
    # AROMATIC RINGS
    # ======================================
    protein_rings = get_protein_rings(
        protein_atoms
    )

    ligand_rings = (
        get_ligand_aromatic_rings(
            ligand_atoms,
            ligand_pdb_file,
            sdf_molecule_index=sdf_molecule_index
        )
    )


    # ======================================
    # PI INTERACTIONS
    # ======================================
    interactions.extend(

        detect_pi_pi(
            protein_rings,
            ligand_rings
        )
    )


    interactions.extend(

        detect_pi_anion(
            protein_atoms,
            ligand_rings
        )
    )

    interactions.extend(

        detect_pi_cation(
            protein_atoms,
            ligand_rings
        )
    )

    interactions.extend(

        detect_pi_alkyl(
            protein_atoms,
            ligand_rings
        )
    )

    interactions.extend(

        detect_amide_pi(
            protein_atoms,
            ligand_rings
        )
    )
    # ======================================
    # H-BONDS
    # ======================================
    ligand_donors, ligand_acceptors = \
    get_ligand_donor_acceptor_atoms(
        ligand_all_atoms,
        ligand_pdb_file,
        sdf_molecule_index=sdf_molecule_index
    )
    
    ligand_positive, ligand_negative = \
    get_ligand_charged_atoms(
        ligand_pdb_file,
        sdf_molecule_index=sdf_molecule_index
    )
    hbond_seen = set()

    # ======================================
    # CARBON H-BONDS
    # ======================================
    for patom in protein_all_atoms:

        for latom in ligand_all_atoms:

            carbon_hbond_interactions = (
                detect_carbon_hbond(
                    patom,
                    latom
                )
            )

            for interaction in carbon_hbond_interactions:

                interactions.append(
                    interaction
                )

    # ======================================
    # SALT BRIDGES + H-BONDS
    # ======================================
    for patom in protein_all_atoms:

        for ligand_index, latom in enumerate(
                ligand_all_atoms):

            # -----------------------------
            # Salt Bridges
            # -----------------------------
            salt_interactions = detect_salt_bridge(
                patom,
                latom,
                ligand_index,
                ligand_positive,
                ligand_negative
            )

            interactions.extend(
                salt_interactions
            )

            # -----------------------------
            # Hydrogen Bonds
            # -----------------------------
            hbond_interactions = detect_hbond(
                patom,
                latom,
                ligand_index,
                ligand_donors,
                ligand_acceptors
            )

            for interaction in hbond_interactions:

                key = (

                    interaction[
                        "protein_atom"
                    ].resname,

                    interaction[
                        "protein_atom"
                    ].resid,

                    interaction[
                        "ligand_atom"
                    ].name,

                    interaction[
                        "type"
                    ]
                )

                if key in hbond_seen:
                    continue

                hbond_seen.add(key)

                interactions.append(
                    interaction
                )

    # ======================================
    # HYDROPHOBIC INTERACTIONS
    # ======================================
    best_hydrophobic = {}

    for patom in protein_atoms:

        for latom in ligand_atoms:

            if not is_hydrophobic(
                    patom,
                    latom):
                continue

            dist = atom_distance(
                patom,
                latom
            )

            if dist > 4.2:
                continue

            key = (
                patom.resname,
                patom.resid
            )

            if (
                key not in best_hydrophobic
                or
                dist < best_hydrophobic[key]["distance"]
            ):

                best_hydrophobic[key] = {

                    "protein_atom":
                        patom,

                    "ligand_atom":
                        latom,

                    "type":
                        "Hydrophobic",

                    "distance":
                        float(dist)
                }

    # ======================================
    # PLIP-LIKE PRIORITIZATION
    # Hydrophobic contacts are removed if the
    # same residue already forms a stronger
    # interaction with the ligand
    # ======================================

    priority_residues = set()

    for interaction in interactions:

        if interaction["type"] in [

            "Pi-Alkyl",
            "Pi-Cation",
            "Pi-Anion",
            "Pi-Pi Stacked",
            "Pi-Pi T-Shaped",
            "Amide-Pi Stacked",
            "H-Bond (Strong)",
            "H-Bond (Moderate)",
            "H-Bond (Weak)"
        ]:

            patom = interaction["protein_atom"]

            priority_residues.add(
                (
                    patom.resname,
                    patom.resid
                )
            )

    filtered_hydrophobic = []

    for interaction in best_hydrophobic.values():

        patom = interaction["protein_atom"]

        key = (
            patom.resname,
            patom.resid
        )

        if key in priority_residues:
            continue

        filtered_hydrophobic.append(
            interaction
        )

    interactions.extend(
        filtered_hydrophobic
    )

    return interactions