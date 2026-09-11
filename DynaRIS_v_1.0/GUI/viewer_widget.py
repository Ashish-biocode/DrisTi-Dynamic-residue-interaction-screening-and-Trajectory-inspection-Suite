"""Embedded 3D viewer helpers and export utilities."""

from pathlib import Path
import json
import csv
import shutil
import subprocess
import tempfile

from rdkit import Chem


def build_viewer_html(protein_path, ligand_path, ligand_result, show_labels=True):
    """Build a standalone 3Dmol.js scene for QWebEngineView."""
    from visualization.viewer import InteractionViewer

    viewer = InteractionViewer(
        protein_path,
        ligand_path,
        ligand_result.get("interactions", []),
        ligand_resname="LIG",
        ligand_index=ligand_result["index"],
        show_labels=show_labels,
    )
    view = viewer.show_structure()

    if not hasattr(view, "_make_html"):
        raise RuntimeError(
            "Installed py3Dmol does not expose the HTML builder required "
            "for the embedded viewer."
        )

    core_html = view._make_html()

    # Keep label coordinates able to be toggled in-place.
    label_data = []
    seen = set()
    for interaction in ligand_result.get("interactions", []):
        try:
            atom = interaction["protein_atom"]
            key = f"{atom.resname}{int(atom.resid)}"
            if key in seen:
                continue
            seen.add(key)
            x, y, z = map(float, atom.position)
            label_data.append({
                "text": key,
                "position": {"x": x, "y": y, "z": z},
            })
        except Exception:
            continue

    labels_json = json.dumps(label_data)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Interaction 3D Viewer</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
html, body {{
    width: 100%;
    height: 100%;
    margin: 0;
    padding: 0;
    overflow: hidden;
    background: #ffffff;
}}
</style>
</head>
<body>
{core_html}
<script>
window.__interactionLabelData = {labels_json};

function __findInteractionViewer() {{
    var keys = Object.keys(window);
    for (var i = 0; i < keys.length; i++) {{
        try {{
            var obj = window[keys[i]];
            if (obj &&
                typeof obj.removeAllLabels === 'function' &&
                typeof obj.addLabel === 'function' &&
                typeof obj.render === 'function') {{
                return obj;
            }}
        }} catch (e) {{}}
    }}
    return null;
}}

window.setInteractionLabels = function(show) {{
    var v = __findInteractionViewer();
    if (!v) return false;

    // Change labels only. Never change the camera.
    v.removeAllLabels();

    if (show) {{
        (__interactionLabelData || []).forEach(function(item) {{
            v.addLabel(item.text, {{
                position: item.position,
                fontColor: "#1f2937",
                fontSize: 14,
                showBackground: true,
                backgroundColor: "#ffffff",
                backgroundOpacity: 0.72,
                borderThickness: 0
            }});
        }});
    }}

    v.render();
    return true;
}};
</script>
</body>
</html>
"""


def _load_ligand_mol(ligand_path, ligand_index):
    supplier = Chem.SDMolSupplier(ligand_path, removeHs=False)
    for index, mol in enumerate(supplier):
        if index == ligand_index:
            if mol is None:
                raise ValueError(
                    f"Could not read SDF molecule index {ligand_index}."
                )
            return mol
    raise ValueError(f"SDF molecule index {ligand_index} was not found.")


def export_protein_pdb(protein_path, output_path):
    Path(output_path).write_bytes(Path(protein_path).read_bytes())


def export_ligand_pdb(ligand_path, ligand_index, output_path):
    mol = _load_ligand_mol(ligand_path, ligand_index)
    Path(output_path).write_text(Chem.MolToPDBBlock(mol), encoding="utf-8")


def export_complex_pdb(protein_path, ligand_path, ligand_index, output_path):
    protein = Path(protein_path).read_text(
        encoding="utf-8", errors="replace"
    )
    mol = _load_ligand_mol(ligand_path, ligand_index)
    ligand = Chem.MolToPDBBlock(mol)

    protein_lines = [
        line for line in protein.splitlines()
        if line.strip() and line.strip() != "END"
    ]
    ligand_lines = [
        line for line in ligand.splitlines()
        if line.strip() and line.strip() != "END"
    ]

    combined = "\n".join(protein_lines + ligand_lines) + "\nEND\n"
    Path(output_path).write_text(combined, encoding="utf-8")


def export_ligand_sdf(ligand_path, ligand_index, output_path):
    mol = _load_ligand_mol(ligand_path, ligand_index)
    writer = Chem.SDWriter(str(output_path))
    try:
        writer.write(mol)
    finally:
        writer.close()


def export_viewer_html(html, output_path):
    Path(output_path).write_text(html, encoding="utf-8")


def _interaction_row(ligand, interaction):
    p = interaction.get("protein_atom")
    l = interaction.get("ligand_atom")

    if p is not None:
        if interaction.get("type") == "Amide-Pi Stacked":
            protein_text = (
                f"{interaction.get('amide_prev_res', 'Unknown')} + "
                f"{interaction.get('amide_curr_res', 'Unknown')}"
            )
        else:
            protein_text = f"{p.resname}-{p.resid}:{p.name}"
    else:
        protein_text = "Unknown"

    ligand_text = str(getattr(l, "name", "Unknown"))
    return [
        ligand.get("name", ""),
        interaction.get("type", ""),
        protein_text,
        ligand_text,
        float(interaction.get("distance", 0.0)),
    ]


def collect_interaction_rows(ligand_results, current_ligand=None,
                             interaction_type="All interactions"):
    rows = []
    for ligand in ligand_results:
        if current_ligand is not None and ligand is not current_ligand:
            continue
        for interaction in ligand.get("interactions", []):
            if (
                interaction_type != "All interactions"
                and interaction.get("type") != interaction_type
            ):
                continue
            rows.append(_interaction_row(ligand, interaction))
    return rows


def export_interactions_csv(rows, output_path):
    with open(output_path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["Ligand", "Interaction", "Protein atom",
             "Ligand atom", "Distance (Å)"]
        )
        writer.writerows(rows)


def export_interactions_txt(rows, output_path, title="Interaction Analysis Results"):
    """Write interactions as a human-readable fixed-width tabular text file."""
    headers = ["Ligand", "Interaction", "Protein atom", "Ligand atom", "Distance (Å)"]
    widths = [28, 28, 30, 22, 14]

    def fmt_row(values):
        return "  ".join(
            str(v)[:widths[i]].ljust(widths[i])
            for i, v in enumerate(values)
        )

    lines = [
        title,
        "=" * sum(widths + [8] * (len(widths) - 1)),
        fmt_row(headers),
        "-" * sum(widths + [8] * (len(widths) - 1)),
    ]

    for row in rows:
        values = [
            row[0],
            row[1],
            row[2],
            row[3],
            f"{float(row[4]):.2f}",
        ]
        lines.append(fmt_row(values))

    lines.extend([
        "-" * sum(widths + [8] * (len(widths) - 1)),
        f"Total interactions: {len(rows)}",
    ])
    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")



def export_ligand_package(
    protein_path,
    ligand_path,
    ligand_result,
    output_path,
    viewer_html=None,
    png_bytes=None,
    show_labels=True,
):
    """Create a self-contained ZIP package for one selected ligand."""
    from pathlib import Path
    import zipfile

    ligand_name = str(ligand_result.get("name", "Ligand")).strip() or "Ligand"
    safe_name = "".join(
        ch if ch.isalnum() or ch in "-_ ." else "_"
        for ch in ligand_name
    ).strip().replace(" ", "_") or "Ligand"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / f"{safe_name}_Interaction_Analysis"
        root.mkdir(parents=True, exist_ok=True)

        complex_pdb = root / "complex.pdb"
        ligand_sdf = root / "ligand.sdf"
        interactions_txt = root / "interactions.txt"
        interactions_csv = root / "interactions.csv"
        viewer_file = root / "viewer.html"
        image_file = root / "interaction_view.png"
        readme = root / "README.txt"

        export_complex_pdb(
            protein_path, ligand_path, ligand_result["index"], complex_pdb
        )
        export_ligand_sdf(
            ligand_path, ligand_result["index"], ligand_sdf
        )

        rows = collect_interaction_rows([ligand_result])
        export_interactions_txt(
            rows,
            interactions_txt,
            title=f"Interaction Analysis Results — {ligand_name}",
        )
        export_interactions_csv(rows, interactions_csv)

        if viewer_html:
            viewer_file.write_text(viewer_html, encoding="utf-8")

        if png_bytes:
            image_file.write_bytes(png_bytes)

        readme.write_text(
            "DynaRIS - Selected Ligand Package\n"
            "===================================================\n\n"
            f"Ligand: {ligand_name}\n"
            f"SDF molecule index: {ligand_result.get('index', '')}\n"
            f"Interactions: {len(rows)}\n\n"
            "Files:\n"
            "  complex.pdb              Protein + selected ligand\n"
            "  ligand.sdf               Selected ligand structure\n"
            "  interactions.txt         Human-readable interaction table\n"
            "  interactions.csv         Tabular interaction data\n"

            "  viewer.html              Embedded 3D viewer page\n"
            "  interaction_view.png     High-resolution 3D image (if captured)\n\n"
            "Open complex.pdb in Discovery Studio or PyMOL for structural analysis.\n",
            encoding="utf-8",
        )

        output = Path(output_path)
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for file in root.iterdir():
                if file.is_file():
                    archive.write(file, arcname=f"{root.name}/{file.name}")


