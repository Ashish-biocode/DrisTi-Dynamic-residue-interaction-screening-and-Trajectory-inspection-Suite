import os
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt


def find_optimal_clusters(data, max_k=8, output_dir="."):
    scores = []
    K = range(2, max_k + 1)

    print("Finding optimal clusters (Silhouette)...")

    for k in K:
        model = KMeans(n_clusters=k, random_state=42)
        labels = model.fit_predict(data)
        score = silhouette_score(data, labels)
        scores.append(score)

    best_k = K[np.argmax(scores)]

    # Save silhouette plot
    plt.figure()
    plt.plot(K, scores, marker='o')
    plt.xlabel("Number of clusters")
    plt.ylabel("Silhouette Score")
    plt.title("Optimal Cluster Selection")
    plot_path = os.path.join(output_dir, "silhouette_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

    print(f"[INFO] Optimal clusters: {best_k}")
    print(f"[INFO] Silhouette plot saved: {plot_path}")

    return best_k


def run_clustering(topology, trajectory,
                   n_clusters=4,
                   method="kmeans",
                   auto_k=False,
                   output_dir="outputs"):

    try:
        print("Loading trajectory...")
        u = mda.Universe(topology, trajectory)

        # ==============================
        # CREATE OUTPUT DIRECTORY
        # ==============================
        output_dir = os.path.join(output_dir, "clustering")
        os.makedirs(output_dir, exist_ok=True)

        # ==============================
        # ALIGN TRAJECTORY
        # ==============================
        print("Aligning trajectory...")
        align.AlignTraj(
            u, u,
            select="protein and name CA",
            in_memory=True
        ).run()

        # ==============================
        # SELECT CA ATOMS
        # ==============================
        selection = u.select_atoms("protein and name CA")

        # ==============================
        # EXTRACT COORDINATES
        # ==============================
        print("Extracting coordinates...")
        coords = []

        for ts in u.trajectory:
            coords.append(selection.positions.flatten())

        coords = np.array(coords)

        print(f"Frames: {coords.shape[0]}")
        print(f"Features: {coords.shape[1]}")

        # ==============================
        # PCA (for clustering + plotting)
        # ==============================
        print("Running PCA...")
        pca = PCA(n_components=10)
        coords_reduced = pca.fit_transform(coords)

        # ==============================
        # AUTO CLUSTER SELECTION
        # ==============================
        if auto_k:
            n_clusters = find_optimal_clusters(coords_reduced, output_dir=output_dir)

        print(f"Using {method} with k={n_clusters}")

        # ==============================
        # CLUSTERING
        # ==============================
        if method == "kmeans":
            model = KMeans(n_clusters=n_clusters, random_state=42)
            labels = model.fit_predict(coords_reduced)

        elif method == "agglomerative":
            model = AgglomerativeClustering(n_clusters=n_clusters)
            labels = model.fit_predict(coords_reduced)

        else:
            raise ValueError("Invalid method. Use 'kmeans' or 'agglomerative'.")

        # ==============================
        # PCA PLOT
        # ==============================
        print("Saving cluster plot...")
        plot_path = os.path.join(output_dir, "cluster_pca_plot.png")

        plt.figure()
        plt.scatter(coords_reduced[:, 0], coords_reduced[:, 1], c=labels)
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("MD Clustering (PCA)")
        plt.savefig(plot_path, dpi=300)
        plt.close()

        print(f"[INFO] Plot saved: {plot_path}")

        # ==============================
        # REPRESENTATIVE STRUCTURES
        # ==============================
        print("Saving representative structures...")

        for i in range(n_clusters):
            cluster_indices = np.where(labels == i)[0]

            if len(cluster_indices) == 0:
                continue

            if method == "kmeans":
                cluster_center = model.cluster_centers_[i]
                distances = np.linalg.norm(
                    coords_reduced[cluster_indices] - cluster_center,
                    axis=1
                )
                rep_frame_index = cluster_indices[np.argmin(distances)]

            else:
                # For agglomerative → choose middle frame
                rep_frame_index = cluster_indices[len(cluster_indices) // 2]

            u.trajectory[rep_frame_index]

            pdb_path = os.path.join(
                output_dir,
                f"cluster_{i}_representative.pdb"
            )

            with mda.Writer(pdb_path) as W:
                W.write(u)

            print(f"Cluster {i}: Frame {rep_frame_index} saved → {pdb_path}")

        # ==============================
        # SAVE LABELS
        # ==============================
        labels_path = os.path.join(output_dir, "cluster_labels.txt")
        np.savetxt(labels_path, labels, fmt="%d")

        print(f"[INFO] Labels saved: {labels_path}")

        # ==============================
        # RETURN
        # ==============================
        return coords_reduced, labels

    except Exception as e:
        print(f"Error: {str(e)}")
        return None