def run_hbond_analysis(topology, trajectory, ligand_resname, output_dir="outputs"):

    import MDAnalysis as mda
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import os

    try:
        os.makedirs(output_dir, exist_ok=True)

        u = mda.Universe(topology, trajectory)

        ligand = u.select_atoms(f"resname {ligand_resname}")
        if len(ligand) == 0:
            print("Ligand not found!")
            return None

        # -----------------------------
        # TIME
        # -----------------------------
        time_ps = np.array([ts.time for ts in u.trajectory])
        time_ps -= time_ps[0]
        time_ns = time_ps / 1000.0

        u.trajectory.rewind()

        # -----------------------------
        # H-BOND ANALYSIS
        # -----------------------------
        h = HBA(
            universe=u,
            donors_sel="protein",
            acceptors_sel=f"resname {ligand_resname}"
        )

        h.d_a_cutoff = 3.5
        h.d_h_a_angle_cutoff = 130.0

        h.run()

        hbonds = h.results.hbonds

        if hbonds is None or len(hbonds) == 0:
            print("No H-bonds found")
            return None

        df = pd.DataFrame(hbonds, columns=[
            'frame', 'donor_idx', 'hydrogen_idx',
            'acceptor_idx', 'distance', 'angle'
        ])

        hbond_counts = df.groupby("frame").size()

        counts_full = np.zeros(len(time_ns))
        counts_full[hbond_counts.index.astype(int)] = hbond_counts.values

        # -----------------------------
        # SAVE DATA FILE
        # -----------------------------
        data_path = os.path.join(output_dir, "hbond_counts.dat")
        np.savetxt(
            data_path,
            np.column_stack((time_ns, counts_full)),
            header="Time(ns) Hbond_Count",
            fmt="%.4f"
        )

        print("H-bond data saved to:", data_path)

        # -----------------------------
        # SAVE PLOT
        # -----------------------------
        plt.figure()
        plt.plot(time_ns, counts_full)
        plt.xlabel("Time (ns)")
        plt.ylabel("Number of H-bonds")
        plt.title("Protein-Ligand H-bonds Over Time")

        plot_path = os.path.join(output_dir, "hbond_plot.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()

        print("H-bond plot saved to:", plot_path)

        # -----------------------------
        # RETURN
        # -----------------------------
        return {
            "x": time_ns,
            "y": counts_full,
            "xlabel": "Time (ns)",
            "ylabel": "Number of H-Bonds"
        }

    except Exception as e:
        print(f"HBA Error: {e}")
        return None