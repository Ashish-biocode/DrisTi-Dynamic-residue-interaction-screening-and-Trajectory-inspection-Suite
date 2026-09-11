from analysis.parser import StructureParser
# Scientific interaction engine
from interaction_rules import detect_interactions
from chemistry.ligand_processor import load_ligands_from_sdf



def analyze_ligands(protein_path, ligand_path):
    """
    Analyze every valid molecule in a multi-molecule SDF.

    Returns:
        protein_atoms
        ligand_results: list of dictionaries containing ligand
                        metadata and its detected interactions.
    """

    parser = StructureParser(
        protein_path,
        ligand_path
    )

    parser.load_structure()
    protein_atoms = parser.get_protein_atoms()

    ligands = load_ligands_from_sdf(
        ligand_path
    )

    ligand_results = []
    for position, ligand in enumerate(ligands, start=1):

        sdf_index = ligand["index"]
        ligand_name = ligand["name"]
        ligand_atoms = parser.get_ligand_atoms(
            molecule_index=sdf_index
        )

        interactions = detect_interactions(
            protein_atoms,
            ligand_atoms,
            ligand_path,
            sdf_molecule_index=sdf_index
        )

        # Add ligand metadata to every interaction.
        for interaction in interactions:
            interaction["ligand_index"] = sdf_index
            interaction["ligand_name"] = ligand_name

        ligand_results.append({
            "index": sdf_index,
            "name": ligand_name,
            "interactions": interactions
        })
    return protein_atoms, ligand_results



if __name__ == "__main__":
    from GUI.main_window import launch_gui
    launch_gui()
