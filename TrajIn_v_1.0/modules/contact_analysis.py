def run_contact_analysis(topology, trajectory, ligand_resname,
                        distance_cutoff=3.5, top_residues=20, output_dir="outputs"):

    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import pandas as pd
    from collections import Counter
    from MDAnalysis.lib.distances import distance_array
    import os

    try:
        print("[INFO] Running Contact Analysis...")

        os.makedirs(output_dir, exist_ok=True)

        output_plot = os.path.join(output_dir, "ligand_contact_frequency.png")

        # -----------------------------
        # Load system
        # -----------------------------
        u = mda.Universe(topology, trajectory)

        ligand = u.select_atoms(f"resname {ligand_resname}")
        protein = u.select_atoms("protein")

        if len(ligand) == 0:
            print("[ERROR] Ligand not found")
            return None

        interacting_residues = []

        # -----------------------------
        # TRAJECTORY LOOP
        # -----------------------------
        for ts in u.trajectory:
            dist_matrix = distance_array(protein.positions, ligand.positions)
            mask = dist_matrix < distance_cutoff

            indices = mask.any(axis=1).nonzero()[0]

            for idx in indices:
                atom = protein[idx]
                interacting_residues.append(f"{atom.resname}{atom.resid}")

        # -----------------------------
        # COUNT FREQUENCY
        # -----------------------------
        residue_counts = Counter(interacting_residues)

        df = pd.DataFrame.from_dict(
            residue_counts, orient="index", columns=["count"]
        )

        df = df.sort_values(by="count", ascending=False).head(top_residues)

        residues = df.index.tolist()
        counts = df["count"].values

        # -----------------------------
        # SAVE STATIC PLOT
        # -----------------------------
        plt.figure(figsize=(10, 6))
        plt.bar(residues, counts)

        plt.xticks(rotation=45)
        plt.xlabel("Residue")
        plt.ylabel("Contact Count")
        plt.title("Top Interacting Residues")

        plt.tight_layout()
        plt.savefig(output_plot, dpi=300)
        plt.close()

        print("[INFO] Contact analysis completed")

        # -----------------------------
        # RETURN FOR GUI
        # -----------------------------
        return {
            "x": residues,
            "y": counts,
            "residues": residues,
            "counts": counts,
            "file": output_plot,
            "xlabel": "Residue",
            "ylabel": "Contact Count",
            "title": "Top Interacting Residues"
        }

    except Exception as e:
        print(f"[ERROR] {str(e)}")
        return None