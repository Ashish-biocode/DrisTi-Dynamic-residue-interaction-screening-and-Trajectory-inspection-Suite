def run_rmsf(topology, trajectory, output_dir="outputs"):

    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    from MDAnalysis.analysis import align
    import matplotlib.pyplot as plt
    import os
    import numpy as np

    # -----------------------------
    # Create output directory
    # -----------------------------
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------
    # Load trajectory
    # -----------------------------
    u = mda.Universe(topology, trajectory)

    # -----------------------------
    # Align trajectory
    # -----------------------------
    aligner = align.AlignTraj(u, u, select="name CA", in_memory=True)
    aligner.run()

    # -----------------------------
    # Select C-alpha atoms
    # -----------------------------
    calpha = u.select_atoms("name CA")

    if len(calpha) == 0:
        raise ValueError("No C-alpha atoms found. Check topology!")

    # -----------------------------
    # RMSF calculation
    # -----------------------------
    rmsf_analysis = rms.RMSF(calpha)
    rmsf_analysis.run()

    rmsf_values = rmsf_analysis.rmsf
    residues = calpha.resids

    # -----------------------------
    # Save data
    # -----------------------------
    data_file = os.path.join(output_dir, "C_alpha_RMSF.dat")

    np.savetxt(
        data_file,
        np.column_stack((residues, rmsf_values)),
        header="Residue   RMSF(Å)"
    )

    # -----------------------------
    # Plot
    # -----------------------------
    plt.figure()
    plt.plot(residues, rmsf_values)

    plt.xlabel("Residue Number")
    plt.ylabel("RMSF (Å)")
    plt.title("C-alpha RMSF")

    plot_file = os.path.join(output_dir, "C_alpha_RMSF.jpeg")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # -----------------------------
    # RETURN
    # -----------------------------
    return {
        "x": residues,
        "y": rmsf_values,
        "xlabel": "Residue Number",
        "ylabel": "RMSF (Å)",
        "title": "C-alpha RMSF"
    }