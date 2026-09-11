import warnings
warnings.filterwarnings("ignore")

import MDAnalysis as mda
from MDAnalysis.analysis import align, pca
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt


def run_pca_fel_analysis(topology, trajectory, output_dir="outputs"):

    try:
        print("[INFO] Starting PCA + FEL Analysis...")

        os.makedirs(output_dir, exist_ok=True)

        # File paths
        pca_csv = os.path.join(output_dir, "PCA_projection_data.csv")
        variance_plot = os.path.join(output_dir, "PCA_variance.png")
        projection_plot = os.path.join(output_dir, "PCA_projection.png")
        fel_plot = os.path.join(output_dir, "Free_Energy_Landscape_PCA.png")

        # -----------------------------
        # LOAD TRAJECTORY
        # -----------------------------
        u = mda.Universe(topology, trajectory)
        selection = u.select_atoms("protein and backbone")

        # -----------------------------
        # ALIGN
        # -----------------------------
        align.AlignTraj(
            u,
            u,
            select="protein and backbone",
            in_memory=True
        ).run()

        # -----------------------------
        # PCA
        # -----------------------------
        PCA = pca.PCA(
            u,
            select="protein and backbone",
            align=True,
            n_components=10
        )
        PCA.run()

        # -----------------------------
        # PROJECTION
        # -----------------------------
        proj = PCA.transform(selection, n_components=2)

        pc1_values = proj[:, 0]
        pc2_values = proj[:, 1]

        # -----------------------------
        # SAVE DATA
        # -----------------------------
        df = pd.DataFrame({
            "PC1": pc1_values,
            "PC2": pc2_values
        })
        df.to_csv(pca_csv, index=False)

        # -----------------------------
        # REPRESENTATIVE STRUCTURES
        # -----------------------------
        frames_to_save = {
            "min_pc1": np.argmin(pc1_values),
            "max_pc1": np.argmax(pc1_values),
            "mean_pc1": np.argmin(np.abs(pc1_values - np.mean(pc1_values)))
        }

        for label, frame_index in frames_to_save.items():
            u.trajectory[frame_index]
            pdb_path = os.path.join(output_dir, f"frame_{label}.pdb")

            with mda.Writer(pdb_path, selection.n_atoms) as W:
                W.write(selection)

        # -----------------------------
        # VARIANCE
        # -----------------------------
        variance = PCA.variance[:10]

        # -----------------------------
        # PLOT: VARIANCE
        # -----------------------------
        plt.figure()
        plt.plot(range(1, len(variance) + 1), variance, marker='o')
        plt.xlabel("Principal Component")
        plt.ylabel("Variance")
        plt.title("PCA Variance")
        plt.savefig(variance_plot)
        plt.close()

        # -----------------------------
        # PLOT: PROJECTION
        # -----------------------------
        plt.figure()
        plt.scatter(pc1_values, pc2_values, s=10)
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("PCA Projection")
        plt.savefig(projection_plot)
        plt.close()

        # -----------------------------
        # FEL CALCULATION
        # -----------------------------
        H, xedges, yedges = np.histogram2d(
            pc1_values,
            pc2_values,
            bins=50,
            density=True
        )

        H = np.where(H == 0, np.nan, H)

        F = -np.log(H)
        F = F - np.nanmin(F)

        # -----------------------------
        # PLOT: FEL
        # -----------------------------
        plt.figure()
        plt.imshow(
            F.T,
            origin='lower',
            aspect='auto',
            extent=[
                xedges[0], xedges[-1],
                yedges[0], yedges[-1]
            ]
        )
        plt.colorbar(label="Free Energy")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("Free Energy Landscape (PCA)")
        plt.savefig(fel_plot)
        plt.close()

        print("[INFO] PCA + FEL completed")

        # -----------------------------
        # RETURN
        # -----------------------------
        return {
            "type": "pca_fel",
            "plots": {
                "variance": variance_plot,
                "projection": projection_plot,
                "fel": fel_plot
            }
        }

    except Exception as e:
        print(f"[ERROR] {str(e)}")
        return None