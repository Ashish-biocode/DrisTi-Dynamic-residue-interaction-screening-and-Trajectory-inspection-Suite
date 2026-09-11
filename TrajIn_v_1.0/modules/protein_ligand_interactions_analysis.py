def run_protein_ligand_interactions_analysis(topology_file, trajectory_file, ligand_resname):

    import os
    import MDAnalysis as mda
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
    from MDAnalysis.lib.distances import distance_array

    try:
        print("[INFO] Loading trajectory...")

        output_dir = os.path.join("outputs", "protein_ligand_interactions_analysis")
        os.makedirs(output_dir, exist_ok=True)

        u = mda.Universe(topology_file, trajectory_file)

        protein = u.select_atoms("protein")
        ligand = u.select_atoms(f"resname {ligand_resname}")

        if len(ligand) == 0:
            raise ValueError("Ligand not found!")

        n_frames = len(u.trajectory)
        residues = protein.residues
        residue_labels = [f"{r.resname}{r.resid}" for r in residues]

        # 🔥 Adaptive sampling (max ~1000 frames)
        step = max(1, n_frames // 1000)
        sampled_frames = len(u.trajectory[::step])

        print(f"[INFO] Total frames: {n_frames}")
        print(f"[INFO] Using every {step} frame(s) → {sampled_frames} frames")

        contact_matrix = np.zeros((len(residues), sampled_frames))

        # Time series arrays
        hbond_frame_counts = np.zeros(sampled_frames)
        hydrophobic_frame_counts = np.zeros(sampled_frames)
        pipi_frame_counts = np.zeros(sampled_frames)

        # -----------------------------
        # Hydrogen bonds
        # -----------------------------
        hbond_PL = HydrogenBondAnalysis(
            universe=u,
            donors_sel="protein",
            acceptors_sel=f"resname {ligand_resname}",
        )

        hbond_LP = HydrogenBondAnalysis(
            universe=u,
            donors_sel=f"resname {ligand_resname}",
            acceptors_sel="protein",
        )

        hbond_PL.run()
        hbond_LP.run()

        hbonds = np.vstack((hbond_PL.results.hbonds, hbond_LP.results.hbonds))

        hbond_residues = []

        for h in hbonds:
            frame = int(h[0])

            # Map full trajectory frame → sampled frame
            mapped_frame = frame // step
            if mapped_frame >= sampled_frames:
                continue

            resname = h[4]
            resid = int(h[5])

            hbond_frame_counts[mapped_frame] += 1
            hbond_residues.append(f"{resname}{resid}")

        # -----------------------------
        # Interaction selections
        # -----------------------------
        hydrophobic_atoms = protein.select_atoms("resname ALA VAL LEU ILE MET PHE TRP PRO TYR")
        aromatic_residues = protein.select_atoms("resname PHE TYR TRP HIS").residues

        hydrophobic_contacts = []
        pi_contacts = []

        # -----------------------------
        # OPTIMIZED TRAJECTORY LOOP
        # -----------------------------
        for frame_index, ts in enumerate(u.trajectory[::step]):

            # 🔥 Vectorized distance (ALL protein atoms vs ligand)
            dist_all = distance_array(protein.positions, ligand.positions)
            contact_atoms = (dist_all < 4.5).any(axis=1)

            # Map atom → residue index
            for atom_idx, atom in enumerate(protein.atoms):
                if contact_atoms[atom_idx]:
                    res_index = atom.residue.ix
                    contact_matrix[res_index, frame_index] = 1

            # -----------------------------
            # Hydrophobic interactions
            # -----------------------------
            if len(hydrophobic_atoms) > 0:
                dist = distance_array(hydrophobic_atoms.positions, ligand.positions)
                mask = dist < 4.5
                indices = mask.any(axis=1).nonzero()[0]

                hydrophobic_frame_counts[frame_index] = len(indices)

                for i in indices:
                    atom = hydrophobic_atoms[i]
                    hydrophobic_contacts.append(f"{atom.resname}{atom.resid}")

            # -----------------------------
            # Pi-Pi interactions
            # -----------------------------
            ligand_centroid = ligand.positions.mean(axis=0)

            pipi_count = 0
            for residue in aromatic_residues:
                ring_atoms = residue.atoms.select_atoms("name CG CD1 CD2 CE1 CE2 CZ")
                if len(ring_atoms) == 0:
                    continue

                ring_centroid = ring_atoms.positions.mean(axis=0)
                dist = np.linalg.norm(ring_centroid - ligand_centroid)

                if dist < 5.0:
                    pipi_count += 1
                    pi_contacts.append(f"{residue.resname}{residue.resid}")

            pipi_frame_counts[frame_index] = pipi_count

        # -----------------------------
        # SAVE DATA
        # -----------------------------
        summary = pd.DataFrame()
        summary["Hydrogen_Bonds"] = pd.Series(Counter(hbond_residues))
        summary["Hydrophobic"] = pd.Series(Counter(hydrophobic_contacts))
        summary["Pi_Stacking"] = pd.Series(Counter(pi_contacts))

        summary.fillna(0, inplace=True)
        summary.to_csv(os.path.join(output_dir, "interaction_summary.csv"))

        contact_frequency = contact_matrix.sum(axis=1) / sampled_frames * 100

        persistence_df = pd.DataFrame({
            "Residue": residue_labels,
            "Contact_%": contact_frequency
        }).sort_values("Contact_%", ascending=False)

        persistence_df.to_csv(os.path.join(output_dir, "contact_persistence.csv"), index=False)

        # -----------------------------
        # PLOTS
        # -----------------------------
        plots = {}

        # Time series
        plt.figure()
        plt.plot(hbond_frame_counts, label="H-Bonds")
        plt.plot(hydrophobic_frame_counts, label="Hydrophobic")
        plt.plot(pipi_frame_counts, label="Pi-Pi")
        plt.legend()
        plt.xlabel("Frame")
        plt.ylabel("Count")
        plt.title("Interaction Time Series")

        path1 = os.path.join(output_dir, "interaction_timeseries.png")
        plt.savefig(path1, dpi=300)
        plt.close()
        plots["timeseries"] = path1

        # Top residues
        top_res = persistence_df.head(15)

        plt.figure()
        plt.bar(top_res["Residue"], top_res["Contact_%"])
        plt.xticks(rotation=45)
        plt.xlabel("Residue")
        plt.ylabel("Contact %")
        plt.title("Top Interacting Residues")

        path2 = os.path.join(output_dir, "top_residues.png")
        plt.savefig(path2, dpi=300)
        plt.close()
        plots["top_residues"] = path2

        # Heatmap
        plt.figure()
        plt.imshow(contact_matrix, aspect='auto')
        plt.xlabel("Frame")
        plt.ylabel("Residue Index")
        plt.title("Contact Heatmap")

        path3 = os.path.join(output_dir, "contact_heatmap.png")
        plt.savefig(path3, dpi=300)
        plt.close()
        plots["heatmap"] = path3

        print("[INFO] Analysis complete")

        return {
            "plots": plots,
            "data": {
                "time_series": {
                    "hbonds": hbond_frame_counts,
                    "hydrophobic": hydrophobic_frame_counts,
                    "pipi": pipi_frame_counts
                },
                "top_residues": persistence_df.head(20)
            }
        }

    except Exception as e:
        print(f"[ERROR] {e}")
        return None