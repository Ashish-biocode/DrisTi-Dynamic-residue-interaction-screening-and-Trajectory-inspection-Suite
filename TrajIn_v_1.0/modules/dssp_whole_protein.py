def run_dssp(topology, trajectory):

    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    try:
        print("Loading trajectory...")
        u = mda.Universe(topology, trajectory)

        # -----------------------------
        # TIME
        # -----------------------------
        time_ps = np.array([ts.time for ts in u.trajectory])
        time_ps = time_ps - time_ps[0]
        time_ns = time_ps / 1000.0

        u.trajectory.rewind()

        # -----------------------------
        # OUTPUT DIR
        # -----------------------------
        base_dir = os.path.dirname(os.path.abspath(trajectory))
        output_dir = os.path.join(base_dir, "outputs")
        os.makedirs(output_dir, exist_ok=True)

        # -----------------------------
        # DSSP
        # -----------------------------
        print("Running DSSP...")
        dssp = DSSP(u).run()
        ss_array = dssp.results.dssp

        # -----------------------------
        # SAVE CSV
        # -----------------------------
        df = pd.DataFrame(ss_array)
        df.insert(0, "Time (ns)", time_ns)
        df.to_csv(os.path.join(output_dir, "dssp_per_frame.csv"), index=False)

        # -----------------------------
        # Mapping
        # -----------------------------
        ss_labels_map = {
            '-': 'Coil', 'H': 'Alpha Helix', 'G': '3-10 Helix',
            'I': 'Pi Helix', 'E': 'Beta Strand', 'B': 'Beta Bridge',
            'T': 'Turn', 'S': 'Bend'
        }

        ss_label_to_int = {
            label: i for i, label in enumerate(sorted(set(ss_labels_map.values())))
        }

        ss_numeric = np.vectorize(
            lambda x: ss_label_to_int.get(ss_labels_map.get(x, "Unknown"), -1)
        )(ss_array)

        # -----------------------------
        # COLORS
        # -----------------------------
        label_colors = {
            'Coil': "#d3d3d3",
            'Alpha Helix': "#ff4d4d",
            '3-10 Helix': "#ffa64d",
            'Pi Helix': "#ffcc66",
            'Beta Strand': "#4da6ff",
            'Beta Bridge': "#80ccff",
            'Turn': "#66ff66",
            'Bend': "#ccff99"
        }

        cmap = ListedColormap(
            [label_colors[label] for label in sorted(label_colors.keys())]
        )

        # -----------------------------
        # PLOT
        # -----------------------------
        fig = plt.figure(figsize=(14, 6))
        ax = fig.add_subplot(111)

        sns.heatmap(
            ss_numeric.T,
            cmap=cmap,
            cbar=False,
            ax=ax,
            xticklabels=False
        )

        # Prevent squeezing
        ax.set_aspect('auto')

        # -----------------------------
        # AXIS
        # -----------------------------
        total_time = time_ns[-1]
        num_frames = len(time_ns)

        if total_time > 0:
            tick_times = np.arange(0, int(total_time) + 1, 10)
            tick_positions = [
                int(t / total_time * (num_frames - 1)) for t in tick_times
            ]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_times)

        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Residue Index")
        ax.set_title("Secondary Structure Evolution (DSSP)")

        # -----------------------------
        # LEGENDS
        # -----------------------------
        legend_elements = [
            Patch(facecolor=color, label=label)
            for label, color in label_colors.items()
        ]

        ax.legend(
            handles=legend_elements,
            bbox_to_anchor=(1.02, 1),
            loc='upper left'
        )

        
        fig.subplots_adjust(right=0.78)

        
        fig.tight_layout(rect=[0, 0, 0.78, 1])

        # -----------------------------
        # SAVE FIGURE
        # -----------------------------
        fig.savefig(
            os.path.join(output_dir, "dssp_heatmap.jpeg"),
            dpi=300,
            bbox_inches='tight'
        )

        print("DSSP completed successfully.")

        # RETURN FIGURE
        return fig

    except Exception as e:
        print(f"Error in DSSP: {str(e)}")
        return None