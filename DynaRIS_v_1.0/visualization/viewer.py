# visualization/viewer.py

import py3Dmol


class InteractionViewer:
    """Build the embedded 3Dmol.js molecular scene."""

    # Distinct colours for the interaction classes.
    
    INTERACTION_COLORS = {
        "h-bond": "#18864b",
        "hydrophobic": "#e08a00",
        "pi-alkyl": "#c13dbb",
        "pi-pi": "#6b42c1",
        "pi-cation": "#2d6cdf",
        "pi-anion": "#d64545",
        "salt bridge": "#b08b00",
        "salt": "#b08b00",
        "carbon h-bond": "#159a9c",
        "amide-pi": "#7a4fb3",
        "polar": "#159a9c",
        "halogen": "#8b5a2b",
    }

    def __init__(
        self,
        pdb_file,
        ligand_file=None,
        interactions=None,
        ligand_resname="UNK",
        ligand_index=0,
        show_labels=True,
    ):
        self.pdb_file = pdb_file
        self.ligand_file = ligand_file
        self.ligand_resname = ligand_resname
        self.ligand_index = ligand_index
        self.interactions = interactions or []
        self.show_labels = bool(show_labels)

    def get_interaction_color(self, interaction_type):
        text = str(interaction_type).lower()
        for key, color in self.INTERACTION_COLORS.items():
            if key in text:
                return color
        return "#6b7280"

    def show_structure(self):
        with open(
            self.pdb_file, "r", encoding="utf-8", errors="replace"
        ) as handle:
            pdb_data = handle.read()

        viewer = py3Dmol.view(width=1400, height=900)
        viewer.addModel(pdb_data, "pdb")
        viewer.setBackgroundColor("#ffffff")

        # Protein
        viewer.setStyle(
            {"model": 0},
            {"cartoon": {"color": "#c9d1db", "opacity": 0.8}},
        )

        if self.ligand_file:
            from rdkit import Chem

            supplier = Chem.SDMolSupplier(
                self.ligand_file, removeHs=False
            )
            selected_mol = None
            for index, mol in enumerate(supplier):
                if index == self.ligand_index:
                    selected_mol = mol
                    break

            if selected_mol is None:
                raise ValueError(
                    f"Could not load SDF molecule index {self.ligand_index}."
                )

            sdf_data = Chem.MolToMolBlock(selected_mol)
            viewer.addModel(sdf_data, "sdf")

            viewer.setStyle(
                {"model": 1},
                {
                    "stick": {
                        "colorscheme": "greenCarbon",
                        "radius": 0.25,
                    }
                },
            )
            viewer.addStyle(
                {"model": 1},
                {
                    "sphere": {
                        "scale": 0.28,
                        "colorscheme": "greenCarbon",
                    }
                },
            )

        # Highlight residues involved in an interaction.
        interacting_residues = set()
        for interaction in self.interactions:
            try:
                interacting_residues.add(
                    int(interaction["protein_atom"].resid)
                )
            except Exception:
                pass

        for resid in interacting_residues:
            viewer.addStyle(
                {"model": 0, "resi": resid},
                {
                    "stick": {
                        "colorscheme": "whiteCarbon",
                        "radius": 0.20,
                    }
                },
            )

        # Draw one short, coloured dotted/dashed connector per interaction.
        labeled_residues = set()

        for interaction in self.interactions:
            try:
                protein_atom = interaction["protein_atom"]
                ligand_atom = interaction["ligand_atom"]
                interaction_type = interaction["type"]
                color = self.get_interaction_color(interaction_type)

                x1, y1, z1 = map(float, protein_atom.position)
                x2, y2, z2 = map(float, ligand_atom.position)

                viewer.addLine(
                    {
                        "start": {"x": x1, "y": y1, "z": z1},
                        "end": {"x": x2, "y": y2, "z": z2},
                        "color": color,
                        "linewidth": 5.0,
                        "dashed": True,
                        "dashLength": 0.16,
                        "gapLength": 0.18,
                    }
                )

                residue_id = (
                    f"{protein_atom.resname}{int(protein_atom.resid)}"
                )

                if self.show_labels and residue_id not in labeled_residues:
                    viewer.addLabel(
                        residue_id,
                        {
                            "position": {
                                "x": x1,
                                "y": y1,
                                "z": z1,
                            },
                            "fontColor": "#1f2937",
                            "fontSize": 14,
                            "showBackground": True,
                            "backgroundColor": "#ffffff",
                            "backgroundOpacity": 0.72,
                            "borderThickness": 0,
                        },
                    )
                    labeled_residues.add(residue_id)

            except Exception as exc:
                pass

        # Fit the complete protein + selected ligand scene.
        viewer.zoomTo()
        viewer.render()
        return viewer
