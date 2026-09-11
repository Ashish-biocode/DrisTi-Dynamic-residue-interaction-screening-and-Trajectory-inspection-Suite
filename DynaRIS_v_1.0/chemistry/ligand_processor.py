from rdkit import Chem
import tempfile
import os
import numpy as np


_SDF_LIGAND_CACHE = {}


# ==========================================
# ELEMENT DETECTION
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


# ==========================================
# LOAD LIGAND(S) FROM SDF
# ==========================================
def load_ligands_from_sdf(sdf_file):
    """
    Load all valid molecules from an SDF.

    Returns a list of dictionaries containing:
        index       : zero-based SDF molecule index
        name        : SDF molecule name (_Name) or generated name
        mol         : RDKit molecule
    """
    cache_key = (
        os.path.abspath(sdf_file),
        os.path.getmtime(sdf_file)
    )

    if cache_key in _SDF_LIGAND_CACHE:
        return _SDF_LIGAND_CACHE[cache_key]

    # Remove stale cache entries for the same file.
    file_path = os.path.abspath(sdf_file)
    stale_keys = [
        key for key in _SDF_LIGAND_CACHE
        if key[0] == file_path and key != cache_key
    ]

    for key in stale_keys:
        del _SDF_LIGAND_CACHE[key]

    supplier = Chem.SDMolSupplier(
        sdf_file,
        removeHs=False
    )

    ligands = []

    for index, mol in enumerate(supplier):
        if mol is None:
            pass
            continue

        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""

        if not name:
            name = f"Ligand_{index + 1}"

        ligands.append({
            "index": index,
            "name": name,
            "mol": mol
        })

    if not ligands:
        raise ValueError(
            f"No valid molecules found in ligand SDF: {sdf_file}"
        )

    _SDF_LIGAND_CACHE[cache_key] = ligands

    return ligands


def create_rdkit_mol_from_ligand(
        sdf_file,
        molecule_index=0):

    ligands = load_ligands_from_sdf(sdf_file)

    matches = [
        item for item in ligands
        if item["index"] == molecule_index
    ]

    if not matches:
        raise IndexError(
            f"SDF molecule index {molecule_index} "
            f"was not found in {sdf_file}."
        )

    mol = matches[0]["mol"]

    pass
    pass
    pass
    pass
    pass

    aromatic_atoms = sum(
        atom.GetIsAromatic()
        for atom in mol.GetAtoms()
    )

    pass
    pass
    pass

    return mol


# ==========================================
# GENERATE SMILES
# ==========================================
def generate_ligand_smiles(mol):

    return Chem.MolToSmiles(
        mol,
        canonical=True
    )

# ==========================================
# AROMATIC RING DETECTION
# ==========================================
def get_aromatic_rings(mol):

    ring_data = []

    conf = mol.GetConformer()

    ring_info = mol.GetRingInfo()

    for ring in ring_info.AtomRings():

        # keep only fully aromatic rings

        if not all(
            mol.GetAtomWithIdx(idx).GetIsAromatic()
            for idx in ring
        ):
            continue

        coords = []

        for idx in ring:

            pos = conf.GetAtomPosition(idx)

            coords.append(
                [
                    pos.x,
                    pos.y,
                    pos.z
                ]
            )

        coords = np.array(coords)

        centroid = coords.mean(axis=0)

        # calculate ring normal

        if len(coords) >= 3:

            v1 = coords[1] - coords[0]

            v2 = coords[2] - coords[0]

            normal = np.cross(v1, v2)

            norm = np.linalg.norm(normal)

            if norm > 0:

                normal = normal / norm

            else:

                normal = np.array(
                    [0.0, 0.0, 1.0]
                )

        else:

            normal = np.array(
                [0.0, 0.0, 1.0]
            )

        ring_data.append(
            {
                "atoms": list(ring),
                "centroid": centroid,
                "normal": normal
            }
        )

    return ring_data


# ==========================================
# RDKIT ATOM WRAPPER
# ==========================================
def get_ligand_atoms_from_sdf(
        sdf_file,
        molecule_index=0):

    mol = create_rdkit_mol_from_ligand(
        sdf_file,
        molecule_index=molecule_index
    )

    conf = mol.GetConformer()

    atoms = []

    for atom in mol.GetAtoms():

        pos = conf.GetAtomPosition(
            atom.GetIdx()
        )

        wrapped = type(
            "LigandAtom",
            (),
            {
                "position": np.array(
                    [
                        pos.x,
                        pos.y,
                        pos.z
                    ],
                    dtype=float
                ),

                "element":
                    atom.GetSymbol(),

                "name":
                    f"{atom.GetSymbol()}"
                    f"{atom.GetIdx()+1}",

                "resname":
                    "LIG",

                "resid":
                    1,

                "rdkit_atom":
                    atom,

                "bonded_atoms":
                    []
            }
        )

        atoms.append(wrapped)

    for bond in mol.GetBonds():

        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        atoms[a1].bonded_atoms.append(
            atoms[a2]
        )

        atoms[a2].bonded_atoms.append(
            atoms[a1]
        )

    pass

    return atoms
