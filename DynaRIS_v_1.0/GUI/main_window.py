"""DynaRIS — Dynamic Residue Interaction Screening Tool.

The scientific interaction engine is intentionally kept separate from the GUI
and is imported through main.py.
"""

import os
import subprocess
import sys

from PyQt5.QtCore import QObject, QThread, QUrl, pyqtSignal, pyqtSlot, Qt
from PyQt5.QtGui import QFont, QPixmap
from PyQt5.QtWidgets import (
    QApplication,
    QAction,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QComboBox,
    QProgressBar,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
    QAbstractItemView,
    QSizePolicy,
)
from PyQt5.QtWebEngineWidgets import QWebEngineView

from main import analyze_ligands
from GUI.viewer_widget import (
    build_viewer_html,
    collect_interaction_rows,
    export_complex_pdb,
    export_ligand_pdb,
    export_ligand_sdf,
    export_protein_pdb,
    export_ligand_package,
    export_viewer_html,
    export_interactions_csv,
    export_interactions_txt,
)


class AnalysisWorker(QObject):
    finished = pyqtSignal(object, object)
    failed = pyqtSignal(str)

    def __init__(self, protein, ligand):
        super().__init__()
        self.protein = protein
        self.ligand = ligand

    @pyqtSlot()
    def run(self):
        try:
            protein_atoms, ligand_results = analyze_ligands(
                self.protein, self.ligand
            )
            self.finished.emit(protein_atoms, ligand_results)
        except Exception as exc:
            self.failed.emit(str(exc))


class InteractionToolGUI(QMainWindow):
    BG = "#f4f7fb"
    PANEL = "#ffffff"
    PANEL_ALT = "#f8fafc"
    TEXT = "#17202a"
    MUTED = "#667085"
    ACCENT = "#2457a6"
    ACCENT_DARK = "#173e78"
    BORDER = "#d9e1eb"
    SUCCESS = "#16834a"

    def __init__(self):
        super().__init__()
        self.setWindowTitle("DynaRIS — Protein–Ligand Interaction Analysis")
        self.resize(1580, 1000)
        self.setMinimumSize(1250, 800)

        self.protein_path = ""
        self.ligand_path = ""
        self.protein_atoms = None
        self.ligand_results = []
        self.current_ligand = None
        self.worker_thread = None
        self.worker = None
        self.current_viewer_html = ""
        self.show_interaction_labels = True

        self._build_menu_bar()
        self._build_ui()
        self._update_action_states()
        self._set_status("Ready — load a protein PDB and ligand SDF to begin.")

    # ------------------------------------------------------------------
    # Menu bar / toolbar
    # ------------------------------------------------------------------
    def _build_menu_bar(self):
        menubar = self.menuBar()

        self.file_menu = menubar.addMenu("&File")
        self.tools_menu = menubar.addMenu("&Tools")
        help_menu = menubar.addMenu("&Help")

        # File
        self.open_protein_action = QAction("Open &Protein PDB…", self)
        self.open_protein_action.triggered.connect(self._browse_protein)
        self.file_menu.addAction(self.open_protein_action)

        self.open_ligand_action = QAction("Open &Ligand SDF…", self)
        self.open_ligand_action.triggered.connect(self._browse_ligand)
        self.file_menu.addAction(self.open_ligand_action)

        self.file_menu.addSeparator()

        self.save_complex_action = QAction("Save &Complex PDB…", self)
        self.save_complex_action.triggered.connect(self._export_complex)
        self.file_menu.addAction(self.save_complex_action)

        self.save_protein_action = QAction("Save &Protein PDB…", self)
        self.save_protein_action.triggered.connect(self._export_protein)
        self.file_menu.addAction(self.save_protein_action)

        self.save_ligand_action = QAction("Save &Ligand PDB…", self)
        self.save_ligand_action.triggered.connect(self._export_ligand)
        self.file_menu.addAction(self.save_ligand_action)

        self.save_ligand_sdf_action = QAction("Save Selected Ligand SDF…", self)
        self.save_ligand_sdf_action.triggered.connect(self._export_ligand_sdf)
        self.file_menu.addAction(self.save_ligand_sdf_action)


        self.export_package_action = QAction("Export Selected Ligand Package (.zip)…", self)
        self.export_package_action.triggered.connect(self._export_ligand_package)
        self.file_menu.addAction(self.export_package_action)

        self.save_html_action = QAction("Save 3D Viewer HTML…", self)
        self.save_html_action.triggered.connect(self._export_html)
        self.file_menu.addAction(self.save_html_action)

        self.file_menu.addSeparator()
        exit_action = QAction("E&xit", self)
        exit_action.triggered.connect(self.close)
        self.file_menu.addAction(exit_action)

        # Tools
        self.analyze_action = QAction("&Run Interaction Analysis", self)
        self.analyze_action.triggered.connect(self._start_analysis)
        self.tools_menu.addAction(self.analyze_action)

        self.tools_menu.addSeparator()

        self.export_current_csv_action = QAction("Current Ligand → CSV…", self)
        self.export_current_csv_action.triggered.connect(lambda: self._export_interactions("current", "csv"))
        self.tools_menu.addAction(self.export_current_csv_action)

        self.export_all_csv_action = QAction("All Ligands → CSV…", self)
        self.export_all_csv_action.triggered.connect(lambda: self._export_interactions("all", "csv"))
        self.tools_menu.addAction(self.export_all_csv_action)


        self.export_current_txt_action = QAction("Current Ligand → Text Table…", self)
        self.export_current_txt_action.triggered.connect(lambda: self._export_interactions("current", "txt"))
        self.tools_menu.addAction(self.export_current_txt_action)

        self.export_all_txt_action = QAction("All Ligands → Text Table…", self)
        self.export_all_txt_action.triggered.connect(lambda: self._export_interactions("all", "txt"))
        self.tools_menu.addAction(self.export_all_txt_action)

        self.export_png_action = QAction("Export High-Resolution 3D PNG…", self)
        self.export_png_action.triggered.connect(self._export_png)
        self.tools_menu.addAction(self.export_png_action)

        self.capture_interface_action = QAction(
            "Capture DynaRIS Interface…", self
        )
        self.capture_interface_action.triggered.connect(
            self._capture_interface
        )
        self.tools_menu.addAction(self.capture_interface_action)

        self.tools_menu.addSeparator()

        self.launch_trajin_action = QAction("Launch &TrajIn", self)
        self.launch_trajin_action.triggered.connect(self._launch_trajin)
        self.tools_menu.addAction(self.launch_trajin_action)

        self.tools_menu.addSeparator()

        self.reset_view_action = QAction("&Reset 3D View", self)
        self.reset_view_action.triggered.connect(self._reset_view)
        self.tools_menu.addAction(self.reset_view_action)

        self.tools_menu.addSeparator()
        self.labels_action = QAction("Show interacting amino-acid labels", self)
        self.labels_action.setCheckable(True)
        self.labels_action.setChecked(True)
        self.labels_action.triggered.connect(self._toggle_interaction_labels)
        self.tools_menu.addAction(self.labels_action)

        # Help
        about_action = QAction("&About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        central.setObjectName("centralWidget")
        self.setCentralWidget(central)

        root = QVBoxLayout(central)
        root.setContentsMargins(18, 14, 18, 12)
        root.setSpacing(10)

        header = QFrame()
        header.setObjectName("headerCard")
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(18, 12, 18, 12)

        title_box = QVBoxLayout()
        title = QLabel("DynaRIS")
        title.setObjectName("appTitle")
        title_box.addWidget(title)

        subtitle = QLabel(
            "Protein–ligand interaction analysis  •  multi-ligand SDF support  •  embedded 3D visualization"
        )
        subtitle.setObjectName("appSubtitle")
        title_box.addWidget(subtitle)

        header_layout.addLayout(title_box, 1)

        brand_meta = QVBoxLayout()
        brand_meta.setSpacing(2)
        product_label = QLabel("Dynamic Residue Interaction Screening Tool")
        product_label.setObjectName("productLabel")
        product_label.setAlignment(Qt.AlignRight)
        version_label = QLabel("Version 1.0")
        version_label.setObjectName("versionLabel")
        version_label.setAlignment(Qt.AlignRight)
        brand_meta.addWidget(product_label)
        brand_meta.addWidget(version_label)
        header_layout.addLayout(brand_meta)
        header_layout.addSpacing(16)

        self.analysis_badge = QLabel("READY")
        self.analysis_badge.setObjectName("statusBadge")
        self.analysis_badge.setAlignment(Qt.AlignCenter)
        self.analysis_badge.setMinimumWidth(92)
        header_layout.addWidget(self.analysis_badge)

        root.addWidget(header)
        root.addWidget(self._build_input_panel())

        # Main split: ligand list on the left; viewer ABOVE tables on right.
        main_split = QSplitter(Qt.Horizontal)
        main_split.setChildrenCollapsible(False)

        main_split.addWidget(self._build_ligand_panel())
        main_split.addWidget(self._build_right_panel())
        main_split.setSizes([245, 1280])

        root.addWidget(main_split, 1)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        self.progress.setMaximumWidth(230)

        status_row = QHBoxLayout()
        status_row.setContentsMargins(4, 0, 4, 0)
        status_row.addWidget(self.progress)

        self.status_label = QLabel()
        self.status_label.setObjectName("statusLabel")
        status_row.addWidget(self.status_label)
        status_row.addStretch()

        root.addLayout(status_row)

        self.setStyleSheet(
            f"""
            QMainWindow, QWidget#centralWidget {{
                background: {self.BG};
                color: {self.TEXT};
            }}

            QFrame#headerCard {{
                background: {self.PANEL};
                border: 1px solid #cfd9e6;
                border-radius: 12px;
            }}

            QLabel#productLabel {{
                color: {self.ACCENT};
                font: 600 10pt "Segoe UI";
            }}

            QLabel#versionLabel {{
                color: {self.MUTED};
                font: 9pt "Segoe UI";
            }}

            QLabel#appTitle {{
                color: {self.TEXT};
                font: 700 25px "Segoe UI";
            }}

            QLabel#appSubtitle {{
                color: {self.MUTED};
                font: 10pt "Segoe UI";
                margin-top: 2px;
            }}

            QLabel#statusBadge {{
                background: #e9f7ef;
                color: {self.SUCCESS};
                border: 1px solid #ccebd8;
                border-radius: 14px;
                padding: 6px 12px;
                font-weight: 700;
            }}

            QFrame#card {{
                background: {self.PANEL};
                border: 1px solid {self.BORDER};
                border-radius: 8px;
            }}

            QLabel#sectionTitle {{
                color: {self.TEXT};
                font: 700 11pt "Segoe UI";
            }}

            QLabel#muted {{
                color: {self.MUTED};
            }}

            QPushButton {{
                background: #f4f7fb;
                color: #253449;
                border: 1px solid #cfd9e6;
                border-radius: 7px;
                padding: 8px 13px;
                font-weight: 600;
            }}

            QPushButton:hover {{
                background: #e8eff8;
                border-color: #b8c8dc;
            }}

            QPushButton:pressed {{
                background: #dce7f5;
            }}

            QPushButton:disabled {{
                color: #9aa5b4;
                background: #f1f3f5;
            }}

            QPushButton#primaryButton {{
                background: {self.ACCENT};
                color: white;
                border: none;
                border-radius: 7px;
                font-weight: 700;
                padding: 9px 20px;
            }}

            QPushButton#primaryButton:hover {{
                background: {self.ACCENT_DARK};
            }}

            QListWidget {{
                background: #fbfcfe;
                border: 1px solid {self.BORDER};
                border-radius: 6px;
                padding: 3px;
                outline: none;
            }}

            QListWidget::item {{
                padding: 9px 8px;
                border-radius: 5px;
            }}

            QListWidget::item:selected {{
                background: #e7effb;
                color: #173e78;
                font-weight: 600;
            }}

            QComboBox {{
                background: white;
                border: 1px solid {self.BORDER};
                border-radius: 5px;
                padding: 6px 9px;
                min-width: 180px;
            }}

            QTableWidget {{
                background: white;
                alternate-background-color: #f8fafc;
                border: 1px solid {self.BORDER};
                border-radius: 6px;
                gridline-color: #e8edf3;
                selection-background-color: #dfeafb;
                selection-color: #17202a;
            }}

            QHeaderView::section {{
                background: #edf2f7;
                color: #344054;
                padding: 8px;
                border: none;
                border-bottom: 1px solid {self.BORDER};
                font-weight: 700;
            }}

            QToolButton {{
                background: transparent;
                border-radius: 5px;
                padding: 6px 9px;
                color: #344054;
            }}

            QToolButton:hover {{
                background: #eef3f8;
            }}

            QMenuBar {{
                background: {self.PANEL};
                border-bottom: 1px solid {self.BORDER};
                padding: 2px 5px;
            }}

            QMenuBar::item {{
                padding: 6px 10px;
            }}

            QMenuBar::item:selected {{
                background: #edf2f7;
                border-radius: 4px;
            }}

            QMenu {{
                background: white;
                border: 1px solid {self.BORDER};
                padding: 5px;
            }}

            QMenu::item {{
                padding: 7px 24px 7px 12px;
                border-radius: 4px;
            }}

            QMenu::item:selected {{
                background: #edf2f7;
            }}
            """
        )

    def _panel(self, object_name="card"):
        frame = QFrame()
        frame.setObjectName(object_name)
        return frame

    def _button(self, text, slot, primary=False):
        b = QPushButton(text)
        b.clicked.connect(slot)
        if primary:
            b.setObjectName("primaryButton")
        return b

    def _build_input_panel(self):
        panel = self._panel()
        layout = QHBoxLayout(panel)
        layout.setContentsMargins(12, 9, 12, 9)
        layout.setSpacing(8)

        protein_caption = QLabel("PROTEIN")
        protein_caption.setObjectName("muted")
        self.protein_label = QLabel("Not selected")
        self.protein_label.setToolTip("Protein PDB file")
        self.protein_label.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Preferred
        )

        ligand_caption = QLabel("LIGANDS")
        ligand_caption.setObjectName("muted")
        self.ligand_label = QLabel("Not selected")
        self.ligand_label.setToolTip("Ligand SDF file")
        self.ligand_label.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Preferred
        )

        layout.addWidget(protein_caption)
        layout.addWidget(self.protein_label, 2)
        layout.addWidget(self._button("Browse PDB", self._browse_protein))
        layout.addSpacing(12)
        layout.addWidget(ligand_caption)
        layout.addWidget(self.ligand_label, 2)
        layout.addWidget(self._button("Browse SDF", self._browse_ligand))
        layout.addSpacing(12)

        analyze_button = self._button(
            "Run Analysis", self._start_analysis, primary=True
        )
        layout.addWidget(analyze_button)

        return panel

    def _build_ligand_panel(self):
        panel = self._panel()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        heading = QLabel("Ligands")
        heading.setObjectName("sectionTitle")
        layout.addWidget(heading)

        hint = QLabel(
            "Select a ligand to inspect its interactions and 3D structure."
        )
        hint.setWordWrap(True)
        hint.setObjectName("muted")
        layout.addWidget(hint)

        self.ligand_list = QListWidget()
        self.ligand_list.currentRowChanged.connect(self._ligand_selected)
        layout.addWidget(self.ligand_list, 1)

        self.open_button = self._button(
            "Show Selected Ligand in 3D", self._focus_current_ligand
        )
        layout.addWidget(self.open_button)

        return panel

    def _build_right_panel(self):
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(9)

        summary = self._panel()
        sl = QHBoxLayout(summary)
        sl.setContentsMargins(12, 9, 12, 9)

        self.summary_label = QLabel("No analysis loaded")
        self.summary_label.setObjectName("sectionTitle")
        sl.addWidget(self.summary_label)
        sl.addStretch()

        filter_label = QLabel("Interaction filter:")
        filter_label.setObjectName("muted")
        sl.addWidget(filter_label)

        self.filter_combo = QComboBox()
        self.filter_combo.addItem("All interactions")
        self.filter_combo.currentTextChanged.connect(self._refresh_table)
        sl.addWidget(self.filter_combo)

        layout.addWidget(summary)

        # Viewer.
        content_split = QSplitter(Qt.Vertical)
        content_split.setChildrenCollapsible(False)

        viewer_frame = self._panel()
        viewer_layout = QVBoxLayout(viewer_frame)
        viewer_layout.setContentsMargins(0, 0, 0, 0)
        viewer_layout.setSpacing(0)

        viewer_header = QHBoxLayout()
        viewer_title = QLabel("3D Structure & Interaction Viewer")
        viewer_title.setObjectName("sectionTitle")
        viewer_header.addWidget(viewer_title)
        viewer_header.addStretch()

        self.viewer_state = QLabel("No ligand loaded")
        self.viewer_state.setObjectName("muted")
        viewer_header.addWidget(self.viewer_state)

        self.reset_button = self._button("Reset View", self._reset_view)
        viewer_header.addWidget(self.reset_button)

        viewer_layout.addLayout(viewer_header)

        legend = QLabel(
            "<b>Interaction legend:</b> "
            "<span style='color:#18864b'>● H-Bond</span> &nbsp; "
            "<span style='color:#e08a00'>● Hydrophobic</span> &nbsp; "
            "<span style='color:#c13dbb'>● Pi-Alkyl</span> &nbsp; "
            "<span style='color:#2d6cdf'>● Pi-Cation</span> &nbsp; "
            "<span style='color:#d64545'>● Pi-Anion</span> &nbsp; "
            "<span style='color:#b08b00'>● Salt Bridge</span> &nbsp; "
            "<span style='color:#6b42c1'>● Pi-Pi</span> &nbsp; "
            "<span style='color:#159a9c'>● Carbon H-Bond</span>"
        )
        legend.setObjectName("muted")
        legend.setTextFormat(Qt.RichText)
        viewer_layout.addWidget(legend)

        self.viewer = QWebEngineView()
        self.viewer.setMinimumHeight(360)
        self.viewer.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Expanding
        )
        viewer_layout.addWidget(self.viewer, 1)

        content_split.addWidget(viewer_frame)

        table_frame = self._panel()
        table_layout = QVBoxLayout(table_frame)
        table_layout.setContentsMargins(9, 9, 9, 9)
        table_layout.setSpacing(6)

        table_header = QHBoxLayout()
        table_title = QLabel("Detected Interactions")
        table_title.setObjectName("sectionTitle")
        table_header.addWidget(table_title)
        table_header.addStretch()

        self.table_count_label = QLabel("0 interactions")
        self.table_count_label.setObjectName("muted")
        table_header.addWidget(self.table_count_label)
        table_layout.addLayout(table_header)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["Ligand", "Interaction", "Protein atom", "Ligand atom", "Distance (Å)"]
        )
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setSortingEnabled(True)

        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)

        table_layout.addWidget(self.table, 1)
        content_split.addWidget(table_frame)

        # Viewer in the upper pane.
        content_split.setSizes([760, 170])
        content_split.setStretchFactor(0, 5)
        content_split.setStretchFactor(1, 1)
        layout.addWidget(content_split, 1)

        self._set_viewer_empty()
        return container

    # ------------------------------------------------------------------
    # File handling
    # ------------------------------------------------------------------
    def _browse_protein(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Protein PDB",
            "",
            "PDB files (*.pdb *.ent);;All files (*)",
        )
        if path:
            self.protein_path = path
            self.protein_label.setText(os.path.basename(path))
            self.protein_label.setToolTip(path)
            self._set_status(f"Protein loaded: {os.path.basename(path)}")
            self._update_action_states()

    def _browse_ligand(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Ligand SDF",
            "",
            "SDF files (*.sdf *.sd);;All files (*)",
        )
        if path:
            self.ligand_path = path
            self.ligand_label.setText(os.path.basename(path))
            self.ligand_label.setToolTip(path)
            self._set_status(f"Ligand library loaded: {os.path.basename(path)}")
            self._update_action_states()

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    def _start_analysis(self):
        if not self.protein_path or not os.path.isfile(self.protein_path):
            QMessageBox.warning(
                self, "Missing protein",
                "Please select a valid protein PDB file."
            )
            return

        if not self.ligand_path or not os.path.isfile(self.ligand_path):
            QMessageBox.warning(
                self, "Missing ligand",
                "Please select a valid ligand SDF file."
            )
            return

        self.analyze_action.setEnabled(False)
        self.progress.setVisible(True)
        self._set_status("Running interaction analysis…")
        self.analysis_badge.setText("ANALYZING")
        self.analysis_badge.setStyleSheet(
            "background:#fff4d6;color:#9a6700;border:1px solid #f0d99a;"
            "border-radius:14px;padding:6px 12px;font-weight:700;"
        )

        self.worker_thread = QThread()
        self.worker = AnalysisWorker(self.protein_path, self.ligand_path)
        self.worker.moveToThread(self.worker_thread)
        self.worker_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self._analysis_finished)
        self.worker.failed.connect(self._analysis_failed)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.failed.connect(self.worker_thread.quit)
        self.worker_thread.finished.connect(self._analysis_thread_cleanup)
        self.worker_thread.start()

    @pyqtSlot(object, object)
    def _analysis_finished(self, protein_atoms, ligand_results):
        self.protein_atoms = protein_atoms
        self.ligand_results = ligand_results
        self.current_ligand = None

        total = sum(len(x["interactions"]) for x in ligand_results)
        self.summary_label.setText(
            f"{len(ligand_results)} ligand(s) analyzed  •  "
            f"{total} interaction(s) detected"
        )

        self.ligand_list.clear()
        for item in ligand_results:
            self.ligand_list.addItem(
                f"{item['name']}   ({len(item['interactions'])})"
            )

        types = sorted(
            {
                interaction["type"]
                for ligand in ligand_results
                for interaction in ligand["interactions"]
            }
        )
        self.filter_combo.blockSignals(True)
        self.filter_combo.clear()
        self.filter_combo.addItem("All interactions")
        self.filter_combo.addItems(types)
        self.filter_combo.blockSignals(False)

        self._refresh_table()
        if ligand_results:
            self.ligand_list.setCurrentRow(0)

        self.analysis_badge.setText("COMPLETE")
        self.analysis_badge.setStyleSheet(
            "background:#e9f7ef;color:#16834a;border:1px solid #ccebd8;"
            "border-radius:14px;padding:6px 12px;font-weight:700;"
        )
        self._set_status("Analysis completed successfully.")
        self._update_action_states()

    @pyqtSlot(str)
    def _analysis_failed(self, message):
        self.analysis_badge.setText("ERROR")
        self.analysis_badge.setStyleSheet(
            "background:#fff0f0;color:#b42318;border:1px solid #f4c7c3;"
            "border-radius:14px;padding:6px 12px;font-weight:700;"
        )
        self._set_status("Analysis failed.")
        QMessageBox.critical(
            self,
            "Analysis error",
            f"The analysis could not be completed.\n\n{message}",
        )

    def _analysis_thread_cleanup(self):
        self.progress.setVisible(False)
        self.analyze_action.setEnabled(True)
        if self.worker is not None:
            self.worker.deleteLater()
        if self.worker_thread is not None:
            self.worker_thread.deleteLater()
        self.worker = None
        self.worker_thread = None
        self._update_action_states()

    # ------------------------------------------------------------------
    # Ligand(s) / table
    # ------------------------------------------------------------------
    def _ligand_selected(self, row):
        if row < 0 or row >= len(self.ligand_results):
            self.current_ligand = None
            self._set_viewer_empty()
            self._refresh_table()
            self._update_action_states()
            return

        self.current_ligand = self.ligand_results[row]
        self._refresh_table()
        self._load_current_viewer()
        self._update_action_states()

    def _focus_current_ligand(self):
        if self.current_ligand is None:
            QMessageBox.information(
                self, "No ligand selected",
                "Select a ligand first."
            )
            return
        self._load_current_viewer()

    def _refresh_table(self, *_):
        selected_type = self.filter_combo.currentText()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)

        rows = collect_interaction_rows(
            self.ligand_results,
            current_ligand=self.current_ligand,
            interaction_type=selected_type,
        )

        for values in rows:
            row = self.table.rowCount()
            self.table.insertRow(row)
            for col, value in enumerate(values):
                if col == 4:
                    item = QTableWidgetItem(f"{float(value):.2f}")
                else:
                    item = QTableWidgetItem(str(value))
                self.table.setItem(row, col, item)

        self.table.setSortingEnabled(True)
        self.table_count_label.setText(
            f"{len(rows)} interaction{'s' if len(rows) != 1 else ''}"
        )

    # ------------------------------------------------------------------
    # 3D viewer
    # ------------------------------------------------------------------
    def _set_viewer_empty(self):
        self.current_viewer_html = ""
        self.viewer_state.setText("No ligand loaded")
        self.viewer.setHtml(
            """
            <html>
            <body style="margin:0;background:#fbfcfe;
            font-family:Segoe UI,Arial;color:#667085;
            display:flex;align-items:center;justify-content:center;
            height:100vh;">
            <div style="text-align:center;">
              <div style="font-size:28px;margin-bottom:10px;">3D</div>
              <div>Select a ligand to load its embedded structure.</div>
            </div>
            </body>
            </html>
            """,
            QUrl("https://3dmol.org/"),
        )

    def _load_current_viewer(self):
        if self.current_ligand is None:
            return

        try:
            html = build_viewer_html(
                self.protein_path,
                self.ligand_path,
                self.current_ligand,
                show_labels=self.show_interaction_labels,
            )
            self.current_viewer_html = html
            self.viewer.setHtml(html, QUrl("https://3dmol.org/"))
            self.viewer_state.setText(
                f"Showing: {self.current_ligand['name']}"
            )
            self._set_status(
                f"3D viewer loaded for {self.current_ligand['name']}."
            )
        except Exception as exc:
            self.viewer_state.setText("Viewer error")
            QMessageBox.critical(
                self,
                "3D viewer error",
                f"Could not load the embedded 3D viewer.\n\n{exc}",
            )

    def _reset_view(self):
        if not self.current_ligand:
            return
        self.viewer.page().runJavaScript(
            """
            (function() {
                var candidates = Object.keys(window).filter(function(k) {
                    return (k.indexOf('viewer_') === 0 ||
                            k === 'viewer') &&
                           window[k] &&
                           typeof window[k].zoomTo === 'function';
                });
                candidates.forEach(function(k) {
                    try {
                        window[k].zoomTo();
                        window[k].render();
                    } catch (e) {}
                });
            })();
            """
        )
        self._set_status("3D view reset.")

    # ------------------------------------------------------------------
    # Exports
    # ------------------------------------------------------------------
    def _require_current(self):
        if self.current_ligand is None:
            QMessageBox.information(
                self, "No ligand selected",
                "Select a ligand first."
            )
            return False
        return True

    def _export_complex(self):
        if not self._require_current():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Complex PDB", "", "PDB files (*.pdb)"
        )
        if not path:
            return
        try:
            export_complex_pdb(
                self.protein_path, self.ligand_path,
                self.current_ligand["index"], path
            )
            self._set_status(f"Saved complex PDB: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))

    def _export_protein(self):
        if not self.protein_path:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Protein PDB", "", "PDB files (*.pdb)"
        )
        if not path:
            return
        try:
            export_protein_pdb(self.protein_path, path)
            self._set_status(f"Saved protein PDB: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))

    def _export_ligand(self):
        if not self._require_current():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Selected Ligand PDB", "", "PDB files (*.pdb)"
        )
        if not path:
            return
        try:
            export_ligand_pdb(
                self.ligand_path, self.current_ligand["index"], path
            )
            self._set_status(f"Saved ligand PDB: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))

    def _export_ligand_sdf(self):
        if not self._require_current():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Selected Ligand SDF", "", "SDF files (*.sdf)"
        )
        if not path:
            return
        try:
            export_ligand_sdf(
                self.ligand_path, self.current_ligand["index"], path
            )
            self._set_status(f"Saved ligand SDF: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))


    def _export_ligand_package(self):
        if not self._require_current():
            return
        default_name = str(self.current_ligand.get("name", "ligand")).strip() or "ligand"
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Selected Ligand Analysis Package",
            f"{default_name}_interaction_package.zip",
            "ZIP packages (*.zip)",
        )
        if not path:
            return

        # Capture a high-resolution image using the exact current camera.
        js = r"""
        (function() {
            function findViewer() {
                var keys = Object.keys(window);
                for (var i = 0; i < keys.length; i++) {
                    try {
                        var obj = window[keys[i]];
                        if (obj && typeof obj.pngURI === 'function' &&
                            typeof obj.getView === 'function' &&
                            typeof obj.setView === 'function') {
                            return obj;
                        }
                    } catch (e) {}
                }
                return null;
            }

            var v = findViewer();
            if (!v) return '';

            try {
                var savedView = v.getView();

                // Keep exactly the same aspect ratio as the live GUI viewer.
                // Changing the aspect ratio changes the projection and makes
                // the same camera appear shifted in the exported PNG.
                var canvas = v.getCanvas ? v.getCanvas() : null;
                var oldW = v.getWidth ? v.getWidth() : (canvas ? canvas.clientWidth : 1400);
                var oldH = v.getHeight ? v.getHeight() : (canvas ? canvas.clientHeight : 900);
                var cssW = canvas && canvas.clientWidth ? canvas.clientWidth : oldW;
                var cssH = canvas && canvas.clientHeight ? canvas.clientHeight : oldH;
                var liveAspect = cssW / Math.max(1, cssH);

                var exportW = 2400;
                var exportH = Math.max(1, Math.round(exportW / liveAspect));

                v.setWidth(exportW);
                v.setHeight(exportH);
                v.setView(savedView);
                v.render();

                var uri = v.pngURI();

                v.setWidth(oldW);
                v.setHeight(oldH);
                v.setView(savedView);
                v.render();

                return uri;
            } catch (e) {
                return '';
            }
        })();
        """


        def finish_package(data):
            try:
                import base64
                png_bytes = None
                if data and str(data).startswith("data:image/png;base64,"):
                    png_bytes = base64.b64decode(str(data).split(",", 1)[1])
                export_ligand_package(
                    self.protein_path,
                    self.ligand_path,
                    self.current_ligand,
                    path,
                    viewer_html=self.current_viewer_html,
                    png_bytes=png_bytes,
                    show_labels=self.show_interaction_labels,
                )
                self._set_status(f"Exported ligand package: {path}")
            except Exception as exc:
                QMessageBox.critical(self, "Package export error", str(exc))

        self.viewer.page().runJavaScript(js, finish_package)

    def _toggle_interaction_labels(self, checked):
        """Toggle labels in-place without rebuilding or changing the camera."""
        self.show_interaction_labels = bool(checked)

        if self.current_ligand is not None:
            js = (
                "if (window.setInteractionLabels) "
                "{setInteractionLabels(%s);}" % (
                    "true" if self.show_interaction_labels else "false"
                )
            )

            def after_toggle(_result):
                # Keep the saved HTML representation in sync for package export.
                try:
                    self.current_viewer_html = build_viewer_html(
                        self.protein_path,
                        self.ligand_path,
                        self.current_ligand,
                        show_labels=self.show_interaction_labels,
                    )
                except Exception:
                    pass

            self.viewer.page().runJavaScript(js, after_toggle)

        self._set_status(
            "Interacting amino-acid labels enabled." if checked
            else "Interacting amino-acid labels hidden."
        )

    def _export_html(self):
        if not self._require_current():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save 3D Viewer HTML", "", "HTML files (*.html)"
        )
        if not path:
            return
        try:
            export_viewer_html(self.current_viewer_html, path)
            self._set_status(f"Saved 3D viewer HTML: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))

    def _export_interactions(self, scope, fmt):
        """Export interactions for the selected ligand or all SDF ligands."""
        if not self.ligand_results:
            QMessageBox.information(
                self, "No results",
                "Run the interaction analysis before exporting."
            )
            return

        if scope == "current":
            if not self._require_current():
                return
            rows = collect_interaction_rows(
                self.ligand_results,
                current_ligand=self.current_ligand,
                interaction_type=self.filter_combo.currentText(),
            )
            scope_label = self.current_ligand["name"]
        else:
            rows = collect_interaction_rows(
                self.ligand_results,
                current_ligand=None,
                interaction_type=self.filter_combo.currentText(),
            )
            scope_label = "all ligands"

        if fmt == "csv":
            path, _ = QFileDialog.getSaveFileName(
                self, f"Export {scope_label} interactions",
                "", "CSV files (*.csv)"
            )
            if not path:
                return
            export_interactions_csv(rows, path)
        else:
            path, _ = QFileDialog.getSaveFileName(
                self, f"Export {scope_label} interactions",
                "", "Text files (*.txt)"
            )
            if not path:
                return
            export_interactions_txt(
                rows, path,
                title=f"Interaction Analysis Results — {scope_label}"
            )

        self._set_status(
            f"Exported {len(rows)} interaction(s) for {scope_label}: {path}"
        )

    def _export_png(self):
        if not self._require_current():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export High-Resolution 3D Viewer Image",
            "", "PNG images (*.png)"
        )
        if not path:
            return

        # Render a larger canvas while preserving the user's exact current
        # camera (rotation, translation and zoom).  This avoids the old
        # behaviour where pngURI() used a different canvas/projection and
        # shifted the complex toward one side of the exported image.
        js = r"""
        (function() {
            function findViewer() {
                var keys = Object.keys(window);
                for (var i = 0; i < keys.length; i++) {
                    try {
                        var obj = window[keys[i]];
                        if (obj && typeof obj.pngURI === 'function' &&
                            typeof obj.getView === 'function' &&
                            typeof obj.setView === 'function') {
                            return obj;
                        }
                    } catch (e) {}
                }
                return null;
            }

            var v = findViewer();
            if (!v) return '';

            try {
                var savedView = v.getView();
                var canvas = v.getCanvas ? v.getCanvas() : null;
                var oldW = v.getWidth ? v.getWidth() : (canvas ? canvas.clientWidth : 1400);
                var oldH = v.getHeight ? v.getHeight() : (canvas ? canvas.clientHeight : 900);
                var cssW = canvas && canvas.clientWidth ? canvas.clientWidth : oldW;
                var cssH = canvas && canvas.clientHeight ? canvas.clientHeight : oldH;
                var liveAspect = cssW / Math.max(1, cssH);

                // Export at higher resolution without changing the live
                // viewer's aspect ratio. This preserves the exact framing
                // and centering chosen by the user.
                var exportW = 2400;
                var exportH = Math.max(1, Math.round(exportW / liveAspect));

                v.setWidth(exportW);
                v.setHeight(exportH);
                v.setView(savedView);
                v.render();

                var uri = v.pngURI();

                // Restore the live viewer exactly as the user left it.
                v.setWidth(oldW);
                v.setHeight(oldH);
                v.setView(savedView);
                v.render();

                return uri;
            } catch (e) {
                return '';
            }
        })();
        """

        def receive_png(data):
            try:
                if not data or not str(data).startswith("data:image/png;base64,"):
                    raise RuntimeError(
                        "The 3D viewer could not provide a high-resolution PNG."
                    )
                import base64
                payload = str(data).split(",", 1)[1]
                with open(path, "wb") as handle:
                    handle.write(base64.b64decode(payload))
                self._set_status(f"Exported high-resolution 3D image: {path}")
            except Exception as exc:
                QMessageBox.critical(self, "Image export error", str(exc))

        self.viewer.page().runJavaScript(js, receive_png)

    def _capture_interface(self):
        """Capture the complete DynaRIS application window as a high-resolution PNG."""

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Capture DynaRIS Interface",
            "DynaRIS_Interface_HighRes.png",
            "PNG images (*.png)"
        )

        if not path:
            return

        try:
            # Capture the complete DynaRIS window
            pixmap = self.grab()

            if pixmap.isNull():
                raise RuntimeError(
                    "Could not capture the DynaRIS interface."
                )

            # Increase output resolution
            scale_factor = 3

            high_res_pixmap = pixmap.scaled(
                pixmap.width() * scale_factor,
                pixmap.height() * scale_factor,
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )

            if not high_res_pixmap.save(path, "PNG"):
                raise RuntimeError(
                    "Could not save the high-resolution interface image."
                )

            self._set_status(
                f"High-resolution interface captured: {path}"
            )

        except Exception as exc:
            QMessageBox.critical(
                self,
                "Interface capture error",
                str(exc)
            )

    def _launch_trajin(self):
        """Launch the bundled TrajIn application as an independent process."""

        # main_window.py
        #     ↓
        # GUI/
        #     ↓
        # DynaRIS_1.0/
        #     ↓
        # ProLIG-Suite/
        suite_dir = os.path.dirname(
            os.path.dirname(
                os.path.dirname(os.path.abspath(__file__))
            )
        )

        trajin_script = os.path.join(
            suite_dir,
            "TrajIn_v_1.0",
            "TrajIn_GUI.py",
        )

        if not os.path.isfile(trajin_script):
            QMessageBox.warning(
                self,
                "TrajIn not found",
                "Could not locate the bundled TrajIn application.\n\n"
                f"Expected location:\n{trajin_script}",
            )
            return

        try:
            subprocess.Popen(
                [sys.executable, trajin_script],
                cwd=os.path.dirname(trajin_script),
            )

            self._set_status("TrajIn launched successfully.")

        except Exception as exc:
            QMessageBox.critical(
                self,
                "Could not launch TrajIn",
                f"TrajIn could not be started.\n\n{exc}",
            )


    # ------------------------------------------------------------------
    # State / utility
    # ------------------------------------------------------------------
    def _update_action_states(self):
        has_analysis = bool(self.ligand_results)
        has_current = self.current_ligand is not None

        self.save_complex_action.setEnabled(has_current)
        self.save_ligand_action.setEnabled(has_current)
        self.save_ligand_sdf_action.setEnabled(has_current)
        self.export_package_action.setEnabled(has_current)
        self.save_html_action.setEnabled(has_current)
        self.export_png_action.setEnabled(has_current)
        self.reset_view_action.setEnabled(has_current)

        self.save_protein_action.setEnabled(bool(self.protein_path))
        for action in (
            self.export_current_csv_action,
            self.export_current_txt_action,
        ):
            action.setEnabled(has_current)

        for action in (
            self.export_all_csv_action,
            self.export_all_txt_action,
        ):
            action.setEnabled(has_analysis)
        self.open_button.setEnabled(has_current)

    def _set_status(self, text):
        self.status_label.setText(text)

    def _show_about(self):
        QMessageBox.about(
            self,
            "About DynaRIS",
            "<div style='font-size:15px'><b>DynaRIS</b></div>"
            "<div style='color:#667085'>Dynamic Residue Interaction Screening Tool</div><br>"
            "Research-oriented analysis of protein–ligand interactions with "
            "multi-ligand SDF support and an embedded interactive 3D viewer.<br><br>"
            "<b>Version:</b> 1.0<br>"
            "<b>Analysis engine:</b> DynaRIS interaction rules<br><br>"
            "Designed for structural bioinformatics and computer-aided drug discovery.",
        )


def launch_gui():
    app = QApplication.instance()
    owns_app = app is None
    if owns_app:
        app = QApplication([])
        app.setStyle("Fusion")

    app.setApplicationName("DynaRIS")
    app.setApplicationDisplayName("DynaRIS")
    app.setOrganizationName("DynaRIS")

    window = InteractionToolGUI()
    window.show()

    if owns_app:
        app.exec()


if __name__ == "__main__":
    launch_gui()
