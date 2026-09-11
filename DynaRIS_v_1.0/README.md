# DynaRIS 1.0

**Dynamic Residue Interaction Screening Tool**

DynaRIS is a PyQt5 desktop application for protein–ligand interaction analysis. It supports multi-molecule SDF input, interaction filtering, an embedded 3D viewer, and export of structures, interaction tables, viewer HTML, and high-resolution images.

## Core workflow
1. Load a protein PDB file.
2. Load an SDF containing one or more ligand molecules.
3. Run the interaction analysis.
4. Select individual ligands from the results list.
5. Inspect interactions and the embedded 3D structure.
6. Export structures, tables, HTML, images, or a selected-ligand package.

## Run
```text
python main.py
```

## Main dependencies
- Python 3.x
- PyQt5
- PyQtWebEngine
- MDAnalysis
- RDKit
- NumPy

The 3D viewer uses 3Dmol.js.
