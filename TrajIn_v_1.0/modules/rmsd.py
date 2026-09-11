def run_rmsd(topology, trajectory, output_dir="outputs"):

    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    os.makedirs(output_dir, exist_ok=True)

    u = mda.Universe(topology, trajectory)

    R = rms.RMSD(u, u, select="backbone")
    R.run()

    time = R.rmsd[:, 1] / 1000
    rmsd_values = R.rmsd[:, 2]

    # -----------------------------
    # SAVE DATA FILE
    # -----------------------------
    data_path = os.path.join(output_dir, "rmsd.dat")
    np.savetxt(
        data_path,
        np.column_stack((time, rmsd_values)),
        header="Time(ns) RMSD(A)",
        fmt="%.4f"
    )

    print("RMSD data saved to:", data_path)

    # -----------------------------
    # SAVE PLOT
    # -----------------------------
    plt.figure()
    plt.plot(time, rmsd_values)
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("Backbone RMSD vs Time")

    plot_path = os.path.join(output_dir, "Backbone_RMSD.jpeg")
    plt.savefig(plot_path, dpi=300)
    plt.close()

    print("RMSD plot saved to:", plot_path)

    # -----------------------------
    # RETURN
    # -----------------------------
    return {
        "x": time,
        "y": rmsd_values,
        "xlabel": "Time (ns)",
        "ylabel": "RMSD (Å)",
        "title": "Backbone RMSD vs Time"
    }