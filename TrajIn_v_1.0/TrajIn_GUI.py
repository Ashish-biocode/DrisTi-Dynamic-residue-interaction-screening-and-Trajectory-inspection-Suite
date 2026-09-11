import sys
import os

from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton,
    QVBoxLayout, QFileDialog, QLabel,
    QInputDialog, QGridLayout, QProgressBar,
    QScrollArea
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QPixmap, QFont

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

# Modules
from modules.rmsd import run_rmsd
from modules.rg import run_rg
from modules.rmsf import run_rmsf
from modules.dssp_whole_protein import run_dssp
from modules.dssp_selected_residues import run_dssp_selected
from modules.HBA_complete import run_hbond_analysis
from modules.contact_analysis import run_contact_analysis
from modules.pca_fel_analysis import run_pca_fel_analysis
from modules.clustering_and_extract_representative_structure import run_clustering
from modules.protein_ligand_interactions_analysis import run_protein_ligand_interactions_analysis
from modules.extract_structure import extract_structure_at_time


# -----------------------------

class Worker(QThread):
    finished = pyqtSignal(object)

    def __init__(self, func, *args):
        super().__init__()
        self.func = func
        self.args = args

    def run(self):
        try:
            result = self.func(*self.args)
        except Exception as e:
            print(f"Error: {e}")
            result = None
        self.finished.emit(result)


class MDToolGUI(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("TrajIn - Trajectory Inspection Platform")
        self.setGeometry(200, 200, 600, 900)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        container = QWidget()
        layout = QVBoxLayout(container)

        base_dir = os.path.dirname(os.path.abspath(__file__))

        # TITLE
        self.title = QLabel("TrajIn")
        self.title.setFont(QFont("Segoe UI", 30, QFont.Bold))
        self.title.setAlignment(Qt.AlignCenter)

        self.subtitle = QLabel("Trajectory Inspection Platform")
        self.subtitle.setAlignment(Qt.AlignCenter)

        layout.addWidget(self.title)
        layout.addWidget(self.subtitle)

        # IMAGE
        for ext in ["jpeg", "jpg", "png"]:
            path = os.path.join(base_dir, f"main.{ext}")
            if os.path.exists(path):
                img = QLabel()
                pixmap = QPixmap(path)
                img.setPixmap(pixmap.scaled(350, 220, Qt.KeepAspectRatio))
                img.setAlignment(Qt.AlignCenter)
                layout.addWidget(img)
                break

        # STATUS
        self.label = QLabel("Load TPR/GRO and XTC files to begin")
        self.label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.label)

        # PROGRESS
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        # PLOT AREA
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        # NAV BUTTONS
        self.btn_prev = QPushButton("⬅ Previous")
        self.btn_next = QPushButton("Next ➡")

        self.btn_prev.setVisible(False)
        self.btn_next.setVisible(False)

        layout.addWidget(self.btn_prev)
        layout.addWidget(self.btn_next)

        self.btn_prev.clicked.connect(self.show_prev_image)
        self.btn_next.clicked.connect(self.show_next_image)

        # FILES
        self.btn_tpr = QPushButton("Load TPR File")
        self.btn_gro = QPushButton("Load GRO File")
        self.btn_xtc = QPushButton("Load Trajectory")

        layout.addWidget(self.btn_tpr)
        layout.addWidget(self.btn_gro)
        layout.addWidget(self.btn_xtc)

        # ANALYSIS
        grid = QGridLayout()

        self.btn_rmsd = QPushButton("RMSD")
        self.btn_rg = QPushButton("Radius of Gyration")
        self.btn_rmsf = QPushButton("RMSF")
        self.btn_dssp = QPushButton("DSSP (Whole Protein)")
        self.btn_dssp_selected = QPushButton("DSSP (Selected Residues)")
        self.btn_hbond = QPushButton("H-Bond Analysis")
        self.btn_contact = QPushButton("Contact Analysis")
        self.btn_pca = QPushButton("PCA and PCA based FEL Analysis")
        self.btn_clustering = QPushButton("Clustering")
        self.btn_interactions = QPushButton("Protein-Ligand Interactions")

        buttons = [
            self.btn_rmsd, self.btn_rg, self.btn_rmsf,
            self.btn_dssp, self.btn_dssp_selected,
            self.btn_hbond, self.btn_contact,
            self.btn_pca, self.btn_clustering,
            self.btn_interactions
        ]

        for i, btn in enumerate(buttons):
            grid.addWidget(btn, i % 5, i // 5)

        layout.addLayout(grid)

        # EXTRACTION
        self.btn_extract = QPushButton("Extract Structure at Time")
        layout.addWidget(self.btn_extract)

        scroll.setWidget(container)
        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll)

        # VARIABLES
        self.topology = None
        self.trajectory = None

        self.image_list = []
        self.image_index = 0

        # CONNECTIONS
        self.btn_tpr.clicked.connect(self.load_tpr)
        self.btn_gro.clicked.connect(self.load_gro)
        self.btn_xtc.clicked.connect(self.load_xtc)

        self.btn_rmsd.clicked.connect(lambda: self.run_analysis(run_rmsd, self.plot_basic))
        self.btn_rg.clicked.connect(lambda: self.run_analysis(run_rg, self.plot_basic))
        self.btn_rmsf.clicked.connect(lambda: self.run_analysis(run_rmsf, self.plot_basic))
        self.btn_dssp.clicked.connect(lambda: self.run_analysis(run_dssp, self.plot_figure))
        self.btn_dssp_selected.clicked.connect(self.run_dssp_selected_with_input)
        self.btn_hbond.clicked.connect(self.run_hbond_with_input)
        self.btn_contact.clicked.connect(self.run_contact_with_input)
        self.btn_pca.clicked.connect(lambda: self.run_analysis(run_pca_fel_analysis, self.plot_pca_fel))

        self.btn_clustering.clicked.connect(self.run_clustering_with_options)

        self.btn_interactions.clicked.connect(self.run_protein_ligand_with_input)

        self.btn_extract.clicked.connect(self.run_extraction)

    def run_extraction(self):
        if self.topology and self.trajectory:

            t, ok_time = QInputDialog.getDouble(self, "Time", "Enter time (ns):")
            if not ok_time:
                return

            mode_choice, ok_mode = QInputDialog.getItem(
                self,
                "Extraction Type",
                "Select what to extract:",
                ["Protein only", "Protein + Ligand", "Full system"],
                0,
                False
            )
            if not ok_mode:
                return

            if mode_choice == "Protein only":
                mode = "protein"
                ligand = None

            elif mode_choice == "Protein + Ligand":
                mode = "protein_ligand"

                ligand_name, ok_lig = QInputDialog.getText(
                    self,
                    "Ligand Name",
                    "Enter ligand residue name (e.g., LIG, ATP):"
                )

                if not ok_lig or not ligand_name.strip():
                    return

                ligand = ligand_name.strip()
            else:
                mode = "full"
                ligand = None

            self.label.setText("Extracting structure...")
            self.progress.setVisible(True)

            self.worker = Worker(
                extract_structure_at_time,
                self.topology,
                self.trajectory,
                t,
                mode,
                ligand
            )

            self.worker.finished.connect(self.analysis_done)
            self.worker.start()

    def run_protein_ligand_with_input(self):
        if self.topology and self.trajectory:
            ligand, ok = QInputDialog.getText(
                self, "Ligand Name",
                "Enter ligand residue name (e.g., LIG, ATP):"
            )
            if not ok or not ligand.strip():
                return
            self.label.setText("Running Protein-Ligand Analysis...")
            self.progress.setVisible(True)
            self.worker = Worker(
                run_protein_ligand_interactions_analysis,
                self.topology,
                self.trajectory,
                ligand.strip()
            )
            self.worker.finished.connect(self.plot_protein_ligand)
            self.worker.start()

    def plot_protein_ligand(self, result):
        self.progress.setVisible(False)
        if result is None:
            self.label.setText("Error ❌")
            return
        plots = result.get("plots", {})
        images = []
        for key in ["timeseries", "top_residues", "heatmap"]:
            if key in plots and os.path.exists(plots[key]):
                images.append(plots[key])
        if not images:
            self.label.setText("No plots found ❌")
            return
        self.image_list = images
        self.image_index = 0
        self.btn_prev.setVisible(True)
        self.btn_next.setVisible(True)
        self.show_image()
        self.label.setText("Completed ✅")

    def run_clustering_with_options(self):
        if self.topology and self.trajectory:
            method, ok1 = QInputDialog.getItem(
                self, "Clustering Method",
                "Select algorithm:",
                ["kmeans", "agglomerative"], 0, False
            )
            if not ok1:
                return
            auto_k_choice, ok2 = QInputDialog.getItem(
                self, "Auto Cluster Selection",
                "Use automatic cluster detection?",
                ["Yes", "No"], 0, False
            )
            if not ok2:
                return
            auto_k = True if auto_k_choice == "Yes" else False
            self.label.setText("Running Clustering...")
            self.progress.setVisible(True)
            self.worker = Worker(
                run_clustering,
                self.topology,
                self.trajectory,
                4,
                method,
                auto_k
            )
            self.worker.finished.connect(self.plot_clustering)
            self.worker.start()

    def plot_clustering(self, result):
        self.progress.setVisible(False)
        output_dir = os.path.join("outputs", "clustering")
        images = []
        for f in ["cluster_pca_plot.png", "silhouette_plot.png"]:
            p = os.path.join(output_dir, f)
            if os.path.exists(p):
                images.append(p)
        if not images:
            self.label.setText("No plots found ❌")
            return
        self.image_list = images
        self.image_index = 0
        self.btn_prev.setVisible(True)
        self.btn_next.setVisible(True)
        self.show_image()
        self.label.setText("Completed ✅")

    def show_image(self):
        import matplotlib.image as mpimg
        img = mpimg.imread(self.image_list[self.image_index])
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.imshow(img)
        ax.axis("off")
        self.canvas.draw_idle()

    def show_next_image(self):
        if self.image_index < len(self.image_list) - 1:
            self.image_index += 1
            self.show_image()

    def show_prev_image(self):
        if self.image_index > 0:
            self.image_index -= 1
            self.show_image()

    def run_analysis(self, func, plot_func=None):
        if self.topology and self.trajectory:
            self.label.setText("Running analysis...")
            self.progress.setVisible(True)
            self.btn_prev.setVisible(False)
            self.btn_next.setVisible(False)
            self.worker = Worker(func, self.topology, self.trajectory)
            if plot_func:
                self.worker.finished.connect(plot_func)
            else:
                self.worker.finished.connect(self.analysis_done)
            self.worker.start()

    def analysis_done(self, result):
        self.progress.setVisible(False)
        self.label.setText("Completed ✅" if result is not None else "Error ❌")

    # Plot_basic
    def plot_basic(self, result):
        self.progress.setVisible(False)

        if result is None:
            self.label.setText("Error ❌")
            return

        # Handle tuple/list with extra values
        if isinstance(result, (tuple, list)):
            if len(result) < 2:
                self.label.setText("Invalid result format ❌")
                return
            x, y = result[0], result[1]
            xlabel = "Frame"
            ylabel = "Value"

        elif isinstance(result, dict):
            x = result.get("x")
            y = result.get("y")
            xlabel = result.get("xlabel", "Frame")
            ylabel = result.get("ylabel", "Value")

            if x is None or y is None:
                self.label.setText("Invalid dict format ❌")
                return

        else:
            self.label.setText("Unsupported result format ❌")
            return

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        ax.plot(x, y)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        self.figure.tight_layout()  # ensures labels are visible properly


        self.canvas.draw_idle()

        self.label.setText("Completed ✅")

    def plot_figure(self, fig):
        self.progress.setVisible(False)
        if fig is None:
            self.label.setText("Error ❌")
            return
        self.figure = fig
        self.canvas.figure = self.figure
        self.toolbar.update()
        self.canvas.draw_idle()
        self.label.setText("Completed ✅")

    def plot_pca_fel(self, result):
        self.progress.setVisible(False)
        if result is None:
            self.label.setText("Error ❌")
            return
        plots = result.get("plots", {})
        self.image_list = [
            plots["variance"],
            plots["projection"],
            plots["fel"]
        ]
        self.image_index = 0
        self.btn_prev.setVisible(True)
        self.btn_next.setVisible(True)
        self.show_image()
        self.label.setText("Completed ✅")

    def run_dssp_selected_with_input(self):
        if self.topology and self.trajectory:
            text, ok = QInputDialog.getText(self, "Select Residues",
                                            "Enter residue numbers (comma-separated):")
            if ok and text.strip():
                self.label.setText("Running DSSP...")
                self.progress.setVisible(True)
                self.worker = Worker(run_dssp_selected, self.topology, self.trajectory, text)
                self.worker.finished.connect(self.plot_figure)
                self.worker.start()

    def run_hbond_with_input(self):
        if self.topology and self.trajectory:
            text, ok = QInputDialog.getText(self, "Ligand Name",
                                            "Enter ligand residue name:")
            if ok and text.strip():
                self.label.setText("Running H-bond...")
                self.progress.setVisible(True)
                self.worker = Worker(run_hbond_analysis, self.topology, self.trajectory, text.strip())
                self.worker.finished.connect(self.plot_basic)
                self.worker.start()

    def run_contact_with_input(self):
        if self.topology and self.trajectory:
            ligand, ok = QInputDialog.getText(self, "Ligand Name",
                                              "Enter ligand resname:")
            if ok and ligand.strip():
                self.label.setText("Running Contact Analysis...")
                self.progress.setVisible(True)
                self.worker = Worker(run_contact_analysis, self.topology, self.trajectory, ligand.strip())
                self.worker.finished.connect(self.plot_basic)
                self.worker.start()

    def load_tpr(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select TPR File", "", "*.tpr")
        if f:
            self.topology = f

    def load_gro(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select GRO File", "", "*.gro")
        if f:
            self.topology = f

    def load_xtc(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Trajectory", "", "*.xtc *.trr")
        if f:
            self.trajectory = f


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MDToolGUI()
    w.show()
    sys.exit(app.exec_())