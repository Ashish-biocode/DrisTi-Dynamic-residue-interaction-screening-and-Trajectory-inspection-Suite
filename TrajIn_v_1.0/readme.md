# 🧬 TrajAn-Trajectory Analyzer

![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)
![Framework](https://img.shields.io/badge/GUI-PyQt5-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/Status-Active-success.svg)

---

## 📌 Overview

This repository presents a **TrajAn-Trajectory Analyzer**, developed to simplify and automate complex computational workflows into an intuitive graphical interface.

The tool integrates multiple structural and dynamical analysis methods within a single platform, eliminating the need for extensive command-line operations and improving accessibility for researchers in computational biology and drug discovery.

---

## 🚀 Key Features

* Integrated platform for MD trajectory analysis
* User-friendly GUI built using PyQt5
* Supports RMSD, RMSF, Rg, hydrogen bonds
* Advanced analysis: PCA, FEL, clustering
* DSSP-based secondary structure evaluation
* Automated workflow execution
* Publication-quality plots and outputs
* Modular and extensible design

---

## 🧪 Modules Included

* **Structural Analysis**: RMSD, RMSF, Radius of Gyration
* **Interaction Analysis**: Hydrogen Bonds, Contact Analysis, Detailed Protein-ligand interactions with key residues
* **Conformational Analysis**: PCA & Free Energy Landscape (FEL)
* **Clustering**: K-means / Agglomerative clustering
* **Secondary Structure**: DSSP
* **Structure Extraction**

---

## ⚙️ Installation

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/MD-Trajectory-GUI.git
cd MD-Trajectory-GUI
```

### 2. Create Environment (Recommended)

```bash
conda create -n md_gui python=3.9
conda activate md_gui
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

---

## ▶️ Usage

Run the main GUI application:

```bash
python TrajAn_GUI.py
```

Then:

* Load your trajectory, .tpr and .gro files
* Select analysis type
* Generate plots and results

---

## 📊 Output

* High-quality plots (Images)
* Numerical data files
* Analysis summaries suitable for publication

---


## 📖 Citation

If you use this tool in your research, please cite:

```text
Gupta, A, Purohit R. (2026). TrajAn: (Trajectory Analyser): An interactive python based graphical user interface (GUI) to analyse molecular dynamics simulation trajectories. [Journal Name]. DOI: XXXXX
```

You can also use the GitHub citation feature.

---

## 📜 License

This project is licensed under the **MIT License**.

---

## ⚠️ Disclaimer

This software is developed for **research and educational purposes only**.
The author is not responsible for misuse or misinterpretation of results.

---

## 👨‍💻 Author

**Ashish Gupta**
Computational Drug Discovery Researcher

---

## 🔗 Acknowledgements

* MDAnalysis
* PyEMMA
* PyQt5
* Scientific Python ecosystem

---

## ⭐ Support

If you find this work useful:

* ⭐ Star this repository
* 📄 Cite the associated paper
* 🤝 Contribute to further development

---
