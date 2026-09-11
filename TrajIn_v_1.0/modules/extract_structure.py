def extract_structure_at_time(
    topology,
    trajectory,
    time_ns,
    selection_type="protein_ligand",
    ligand_resname=None,
    output_dir="outputs"
):

    import MDAnalysis as mda
    import os

    try:
        print(f"[INFO] Extracting structure at {time_ns} ns | Mode: {selection_type}")

        # -----------------------------
        # Output directory
        # -----------------------------
        output_dir = os.path.join(output_dir, "structure_extraction")
        os.makedirs(output_dir, exist_ok=True)

        # -----------------------------
        # Load Universe
        # -----------------------------
        u = mda.Universe(topology, trajectory)

        # -----------------------------
        # Convert ns → ps
        # -----------------------------
        target_time_ps = time_ns * 1000

        closest_frame = None
        min_diff = float("inf")

        # -----------------------------
        # Finding closest frame
        # -----------------------------
        for ts in u.trajectory:
            diff = abs(ts.time - target_time_ps)
            if diff < min_diff:
                min_diff = diff
                closest_frame = ts.frame

        if closest_frame is None:
            print("[ERROR] No frame found!")
            return None

        print(f"[INFO] Closest frame: {closest_frame}, time diff: {min_diff:.2f} ps")

        # Move to frame
        u.trajectory[closest_frame]

        # -----------------------------
        # Selection
        # -----------------------------
        if selection_type == "full":
            atoms = u.atoms
            filename_tag = "full_system"

        elif selection_type == "protein":
            atoms = u.select_atoms("protein")
            filename_tag = "protein_only"

        elif selection_type == "protein_ligand":

            # If ligand name provided → use it
            if ligand_resname:
                ligand = u.select_atoms(f"resname {ligand_resname}")
            else:
                # Try auto-detect (exclude common solvent)
                ligand = u.select_atoms("not protein and not resname SOL")

            if len(ligand) == 0:
                print("[WARNING] No ligand found, saving protein only")
                atoms = u.select_atoms("protein")
                filename_tag = "protein_only"
            else:
                atoms = u.select_atoms(f"protein or resname {ligand.residues[0].resname}")
                filename_tag = "protein_ligand"

        else:
            raise ValueError("Invalid selection_type")

        # -----------------------------
        # Save structure
        # -----------------------------
        output_file = os.path.join(
            output_dir,
            f"{filename_tag}_{time_ns}ns.pdb"
        )

        atoms.write(output_file)

        print(f"[INFO] Structure saved at: {output_file}")

        # -----------------------------
        # RETURN
        # -----------------------------
        return {
            "type": "structure",
            "file": output_file,
            "selection": selection_type
        }

    except Exception as e:
        print(f"[ERROR] {str(e)}")
        return None