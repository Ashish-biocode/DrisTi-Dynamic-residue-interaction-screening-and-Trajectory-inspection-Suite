import MDAnalysis as mda
import numpy as np


class StructureParser:

    def __init__(
        self,
        protein_file,
        ligand_file=None):

        self.protein_file = protein_file
        self.ligand_file = ligand_file

        self.universe = None

    # ==========================================
    # LOAD STRUCTURE
    # ==========================================
    def load_structure(self):

        self.universe = mda.Universe(
            self.protein_file
        )

        return self.universe


    # ==========================================
    # CLEAN LIGAND SELECTION
    # ==========================================
    def _get_clean_ligand_selection(self):

        exclude_resnames = [

            "HOH",
            "WAT",

            "NA",
            "CL",
            "K",
            "MG",
            "CA",

            "ZN",
            "FE",
            "CU",
            "MN"
        ]

        selection = (

            "not protein "
            "and not resname "
            +
            " ".join(exclude_resnames)
        )

        return self.universe.select_atoms(
            selection
        )

    # ==========================================
    # LIGAND DETECTION
    # ==========================================
    def get_ligands(self):

        ligands = self._get_clean_ligand_selection()

        return ligands.residues

    # ==========================================
    # GET LIGAND SMILES
    # ==========================================
    def get_ligand_smiles(self):

        """
        Return ligand SMILES.

        Add ligand mappings here.
        Later this can be upgraded to:
        - CCD lookup
        - PubChem lookup
        - SDF parser
        """

        ligand_smiles_map = {

            # ==================================
            # EXAMPLE MAPPINGS
            # ==================================

            # "ATP": "O=P(O)(O)OP(=O)(O)OP(=O)(O)OC",

            # ==================================
            # CURRENT LIGAND
            # ==================================

            "UNK": "YOUR_SMILES_HERE"
        }

        ligands = self.get_ligands()

        if not ligands:
            return None

        ligand_name = (
            ligands[0]
            .resname
            .strip()
        )

        smiles = ligand_smiles_map.get(
            ligand_name,
            None
        )


        return smiles

    # ==========================================
    # ATOM WRAPPER
    # ==========================================
    def _wrap_atom(self, atom):

        element = "C"

        try:

            if (
                hasattr(atom, "element")
                and atom.element
            ):

                element = (
                    atom.element
                    .strip()
                    .upper()
                )

            else:

                name = atom.name.strip()

                if name[:2].upper() in [
                    "CL", "BR", "NA",
                    "MG", "ZN", "FE",
                    "CA", "CU", "MN"
                ]:

                    element = name[:2].upper()

                else:

                    element = name[0].upper()

        except Exception:

            element = atom.name[0].upper()

        return type(

            "AtomObj",
            (),
            {

                "position":
                    np.array(
                        atom.position,
                        dtype=float
                    ),

                "element":
                    element,

                "resid":
                    int(atom.resid),

                "resname":
                    atom.resname,

                "name":
                    atom.name,


                # NEW
                "mda_atom": atom
            }
        )

    # ==========================================
    # GET PROTEIN ATOMS
    # ==========================================
    def get_protein_atoms(self):

        protein = self.universe.select_atoms(
            "protein"
        )

        wrapped_atoms = []

        for atom in protein:

            try:

                if atom.element == "H":
                    continue

            except Exception:
                pass

            wrapped_atoms.append(
                self._wrap_atom(atom)
            )

        return wrapped_atoms

    # ==========================================
    # GET LIGAND ATOMS
    # ==========================================
    def get_ligand_atoms(
            self,
            molecule_index=0):

        from chemistry.ligand_processor import (
            get_ligand_atoms_from_sdf
        )

        return get_ligand_atoms_from_sdf(
            self.ligand_file,
            molecule_index=molecule_index
        )