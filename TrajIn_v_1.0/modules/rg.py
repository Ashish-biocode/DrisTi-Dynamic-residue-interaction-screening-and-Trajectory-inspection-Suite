def run_rg(topology, trajectory, output_dir="outputs"):

    import MDAnalysis as mda
    import numpy as np
    import matplotlib.pyplot as plt
    import os

    # -----------------------------
    # Create outputs folder
    # -----------------------------
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------
    # Load trajectory
    # -----------------------------
    u = mda.Universe(topology, trajectory)

    # Select protein
    protein = u.select_atoms("protein")

    time = []
    rgyr = []

    # -----------------------------
    # Calculate Radius of Gyration
    # -----------------------------
    for ts in u.trajectory:

        # ps → ns
        time.append(ts.time / 1000)

        # Å → nm
        rg = protein.radius_of_gyration() / 10
        rgyr.append(rg)

    time = np.array(time)
    rgyr = np.array(rgyr)

    # -----------------------------
    # Save data
    # -----------------------------
    data_file = os.path.join(output_dir, "radius_of_gyration.dat")
    plot_file = os.path.join(output_dir, "radius_of_gyration.jpeg")

    np.savetxt(
        data_file,
        np.column_stack((time, rgyr)),
        header="Time(ns)   Rg(nm)"
    )

    # -----------------------------
    # Save plot
    # -----------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(time, rgyr)

    plt.xlabel("Time (ns)")
    plt.ylabel("Radius of Gyration (nm)")
    plt.title("Radius of Gyration vs Time")

    plt.grid()

    plt.savefig(plot_file, dpi=300)
    plt.close()

    print("Radius of Gyration calculation completed.")
    print(f"Plot saved at: {plot_file}")
    print(f"Data saved at: {data_file}")

    # -----------------------------
    # RETURN
    # -----------------------------
    return {
        "x": time,
        "y": rgyr,
        "xlabel": "Time (ns)",
        "ylabel": "Radius of Gyration (nm)",
        "title": "Radius of Gyration vs Time"
    }