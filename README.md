# DrisTi — Dynamic Residue Interaction Screening and Trajectory Inspection Suite

**DrisTi** is a Python-based software suite for integrated **protein–ligand interaction analysis and molecular dynamics (MD) trajectory analysis**. It combines two complementary tools, **DynaRIS** and **TrajIn**, into a unified workflow for structural and dynamic characterization of protein–ligand systems.

## Overview

Molecular docking and molecular dynamics simulations are widely used in structure-based drug discovery, but interaction analysis and trajectory analysis are often performed using separate tools and workflows. DrisTi brings these analyses together to facilitate the systematic investigation of protein–ligand binding interactions and their evolution during molecular dynamics simulations.

### Components

| Tool        | Description                                                                                                                                   |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| **DynaRIS** | Protein–ligand interaction analysis with support for multiple ligands from a single SDF file and embedded 3D visualization.                   |
| **TrajIn**  | Analysis and visualization of molecular dynamics trajectories, including structural, conformational, and protein–ligand interaction analyses. |

## Key Features

### DynaRIS

* Automated protein–ligand interaction detection
* Supports **multiple ligands within a single SDF file**
* Hydrogen-bond analysis
* Hydrophobic interactions
* π-related interactions
* Salt-bridge and related electrostatic interactions
* Embedded interactive **3D molecular visualization**
* Interaction labels and visual interaction lines
* Individual and multi-ligand analysis
* Export of interaction results
* Suitable for rapid analysis of docked protein–ligand complexes

### DynaRIS Interface

<img width="2813" height="1503" alt="sc_1" src="https://github.com/user-attachments/assets/5d4247f3-0873-424e-942c-bf6d28b09606" />

### Key Capabilities

<img width="3047" height="1627" alt="sc_2" src="https://github.com/user-attachments/assets/ec4bf98b-1144-49a1-87b3-a50a4759acde" />

<img width="2608" height="1989" alt="sc_3" src="https://github.com/user-attachments/assets/bf169ba5-7ce9-4c58-8451-8804d1fd60c9" />

<img width="2939" height="1933" alt="sc_4" src="https://github.com/user-attachments/assets/f00709af-64f8-4760-8147-cd4c033769f7" />


### TrajIn

* Molecular dynamics trajectory analysis
* RMSD analysis
* RMSF analysis
* Radius of gyration
* Protein–ligand contact analysis
* Hydrogen-bond analysis
* Secondary-structure analysis
* Principal component analysis (PCA)
* Free-energy landscape analysis
* Clustering and representative-structure extraction
* Structural extraction from trajectories

## DrisTi Workflow

```text
Protein–Ligand Complex
          │
          ▼
       DynaRIS
          │
          ├── Protein–ligand interactions
          ├── Interaction visualization
          └── Multi-ligand analysis
          │
          ▼
   Molecular Dynamics
       Trajectory
          │
          ▼
        TrajIn
          │
          ├── Structural analysis
          ├── Conformational analysis
          ├── Contact analysis
          └── Dynamic interaction analysis
```

## Repository Structure

```text
DrisTi/
│
├── DynaRIS_v_1.0/
│   ├── GUI/
│   ├── analysis/
│   ├── chemistry/
│   ├── visualization/
│   ├── test_files/
│   ├── interaction_rules.py
│   ├── main.py
│   ├── README.md
│   └── INSTALLATION.txt
│
├── TrajIn_v_1.0/
│   ├── modules/
│   ├── outputs/
│   ├── TrajIn_GUI.py
│   ├── requirements.txt
│   └── readme.md
│
├── LICENSE
└── README.md
```

## Installation

Clone the repository:

```bash
git clone https://github.com/Ashish-biocode/DrisTi-Dynamic-residue-interaction-screening-and-Trajectory-inspection-Suite.git
```

Enter the repository:

```bash
cd DrisTi-Dynamic-residue-interaction-screening-and-Trajectory-inspection-Suite
```

Installation instructions and dependencies for each component are provided in their respective folders.

### DynaRIS

See:

```text
DynaRIS_v_1.0/INSTALLATION.txt
```

### TrajIn

See:

```text
TrajIn_v_1.0/requirements.txt
```

## Usage

### DynaRIS

DynaRIS can be launched from its directory using:

```bash
python main.py
```

The tool accepts a protein structure and ligand SDF file and performs automated protein–ligand interaction analysis.

A single SDF file can contain multiple ligands, allowing sequential analysis of multiple protein–ligand complexes within the same workflow.

### TrajIn

TrajIn can be launched using:

```bash
python TrajIn_GUI.py
```

Users can select the required trajectory and topology files and perform the available structural, conformational, and interaction analyses through the graphical interface.

## Test Data

Example protein and ligand files are provided with DynaRIS in:

```text
DynaRIS_v_1.0/test_files/
```

These files can be used to test the installation and basic functionality of the software.

## Applications

DrisTi is intended for computational studies involving:

* Protein–ligand binding analysis
* Structure-based drug discovery
* Molecular docking analysis
* Molecular dynamics simulations
* Binding-site characterization
* Interaction persistence and evolution
* Comparative analysis of ligand binding modes

## Version

**DrisTi v1.0**

Components:

* **DynaRIS v1.0**
* **TrajIn v1.0**

## Citation

If you use DrisTi in your research, please cite the associated publication:

**Citation information will be added following publication.**

## License

DrisTi is distributed under the license provided in this repository.

See the [`LICENSE`](LICENSE) file for details.

## Contact

For questions, suggestions, bug reports, or collaboration related to DrisTi, please use the **GitHub Issues** section of this repository.

---

**DrisTi — Dynamic Residue Interaction Screening and Trajectory Inspection Suite**

*A unified workflow for protein–ligand interaction and molecular dynamics trajectory analysis.*
