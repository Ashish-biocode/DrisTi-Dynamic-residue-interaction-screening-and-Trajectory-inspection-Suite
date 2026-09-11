def run_dssp_selected(topology, trajectory, residue_input, output_dir="outputs"):

    import os
    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    print("Running DSSP for selected residues...")

    os.makedirs(output_dir, exist_ok=True)

    output_csv = os.path.join(output_dir, "dssp_selected_residues.csv")
    heatmap_jpeg = os.path.join(output_dir, "dssp_selected_residues_heatmap.jpg")

    # -----------------------------
    # Residue input
    # -----------------------------
    if isinstance(residue_input, str):
        binding_site_resids = [int(x.strip()) for x in residue_input.split(",")]
    else:
        binding_site_resids = list(residue_input)

    # -----------------------------
    # Load Universe
    # -----------------------------
    u = mda.Universe(topology, trajectory)
    protein = u.select_atoms("protein")

    # Time
    time_ps = np.array([ts.time for ts in u.trajectory])
    time_ps -= time_ps[0]
    time_ns = time_ps / 1000.0
    u.trajectory.rewind()

    # -----------------------------
    # DSSP
    # -----------------------------
    dssp = DSSP(u, protein).run()

    ss_array = dssp.results.dssp
    residue_ids = list(dssp.results.resids)

    valid_res = [r for r in binding_site_resids if r in residue_ids]
    if not valid_res:
        raise ValueError("Invalid residues")

    df = pd.DataFrame(ss_array, columns=residue_ids)[valid_res]
    df.insert(0, "Time (ns)", time_ns)
    df.to_csv(output_csv, index=False)

    # -----------------------------
    # Mapping
    # -----------------------------
    ss_map = {'-':0,'H':1,'G':2,'I':3,'E':4,'B':5,'T':6,'S':7}
    ss_numeric = df.drop(columns=["Time (ns)"]).replace(ss_map).to_numpy()

    colors = [
        "#d3d3d3","#ff4d4d","#ffa64d","#ffcc66",
        "#4da6ff","#80ccff","#66ff66","#ccff99"
    ]
    cmap = ListedColormap(colors)

    # -----------------------------
    # CREATE FIGURE
    # -----------------------------
    fig, ax = plt.subplots(figsize=(10, 5))

    sns.heatmap(
        ss_numeric.T,
        cmap=cmap,
        cbar=False,
        ax=ax,
        xticklabels=False
    )

    # Axis
    total_time = time_ns[-1]
    num_frames = len(time_ns)

    if total_time > 0:
        tick_times = np.arange(0, int(total_time)+1, 10)
        tick_pos = [int(t/total_time*(num_frames-1)) for t in tick_times]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_times)

    ax.set_yticks(np.arange(len(valid_res)) + 0.5)
    ax.set_yticklabels(valid_res)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Residue ID")
    ax.set_title("DSSP - Selected Residues")

    # Legend
    legend_elements = [
        Patch(facecolor=colors[i], label=label)
        for i, label in enumerate([
            "Coil","Alpha Helix","3-10 Helix","Pi Helix",
            "Beta Strand","Beta Bridge","Turn","Bend"
        ])
    ]

    ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()

    fig.savefig(heatmap_jpeg, dpi=300, bbox_inches='tight')

    print("DSSP selected completed")

    # RETURN FIGURE
    return fig