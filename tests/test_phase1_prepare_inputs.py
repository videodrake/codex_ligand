from __future__ import annotations

import csv
from pathlib import Path

from egfr_pipeline.phase1 import prepare_inputs as module


def _atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resnum: int,
    *,
    element: str | None = None,
    occupancy: float = 1.0,
    altloc: str = " ",
) -> str:
    element = element or atom_name.strip().lstrip("0123456789")[0]
    x = float(serial)
    y = float(serial) + 0.5
    z = float(serial) + 1.0
    return (
        f"ATOM  {serial:5d} {atom_name:>4}{altloc}{resname:>3} {chain}{resnum:4d}"
        f"    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{0.00:6.2f}          {element:>2}\n"
    )


def _write_pdb(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "TER\nEND\n", encoding="utf-8")


def test_normalize_rosetta_input_lines_handles_charmm_aliases() -> None:
    lines = [
        _atom_line(1, "N", "ILE", "X", 699, element="N"),
        _atom_line(2, "CA", "ILE", "X", 699, element="C"),
        _atom_line(3, "C", "ILE", "X", 699, element="C"),
        _atom_line(4, "O", "ILE", "X", 699, element="O"),
        _atom_line(5, "CD", "ILE", "X", 699, element="C"),
        _atom_line(6, "HN", "ILE", "X", 699, element="H"),
        _atom_line(7, "N", "HSD", "X", 700, element="N"),
        _atom_line(8, "CA", "HSD", "X", 700, element="C"),
        _atom_line(9, "C", "HSD", "X", 700, element="C"),
        _atom_line(10, "O", "HSD", "X", 700, element="O"),
        _atom_line(11, "CB", "HSD", "X", 700, element="C"),
        _atom_line(12, "HN", "HSD", "X", 700, element="H"),
        _atom_line(13, "CB", "SER", "X", 701, element="C", altloc="B"),
        _atom_line(14, "OG", "SER", "X", 701, element="O", occupancy=0.0),
    ]

    normalized, stats = module.normalize_rosetta_input_lines(lines)
    normalized_text = "".join(normalized)

    assert " HSD " not in normalized_text
    assert " HIS " in normalized_text
    assert " CD1" in normalized_text
    assert " CD " not in normalized_text
    assert not any(module._is_hydrogen_atom(line) for line in normalized)
    assert stats["hydrogens_removed"] == 2
    assert stats["alternate_location_removed"] == 1
    assert stats["zero_occupancy_removed"] == 1
    assert stats["residue_renames"]["HSD->HIS"] == 5
    assert stats["atom_renames"]["ILE:CD->CD1"] == 1


def test_prepare_phase1_inputs_writes_rosetta_normalization_metadata(
    tmp_path: Path,
    monkeypatch,
) -> None:
    receptor_path = tmp_path / "input" / "receptors" / "demo_receptor.pdb"
    partner_path = tmp_path / "input" / "PPI" / "demo_partner.pdb"

    receptor_lines = [
        _atom_line(1, "N", "ILE", "X", 699, element="N"),
        _atom_line(2, "CA", "ILE", "X", 699, element="C"),
        _atom_line(3, "C", "ILE", "X", 699, element="C"),
        _atom_line(4, "O", "ILE", "X", 699, element="O"),
        _atom_line(5, "CD", "ILE", "X", 699, element="C"),
        _atom_line(6, "HN", "ILE", "X", 699, element="H"),
        _atom_line(7, "N", "HSD", "X", 700, element="N"),
        _atom_line(8, "CA", "HSD", "X", 700, element="C"),
        _atom_line(9, "C", "HSD", "X", 700, element="C"),
        _atom_line(10, "O", "HSD", "X", 700, element="O"),
        _atom_line(11, "CB", "HSD", "X", 700, element="C"),
        _atom_line(12, "HN", "HSD", "X", 700, element="H"),
    ]
    partner_lines = [
        _atom_line(1, "N", "SER", "A", 955, element="N"),
        _atom_line(2, "CA", "SER", "A", 955, element="C"),
        _atom_line(3, "C", "SER", "A", 955, element="C"),
        _atom_line(4, "O", "SER", "A", 955, element="O"),
        _atom_line(5, "CB", "SER", "A", 955, element="C"),
        _atom_line(6, "HN", "SER", "A", 955, element="H"),
        _atom_line(7, "N", "VAL", "A", 956, element="N"),
        _atom_line(8, "CA", "VAL", "A", 956, element="C"),
        _atom_line(9, "C", "VAL", "A", 956, element="C"),
        _atom_line(10, "O", "VAL", "A", 956, element="O"),
        _atom_line(11, "CB", "VAL", "A", 956, element="C"),
        _atom_line(12, "HN", "VAL", "A", 956, element="H"),
    ]
    _write_pdb(receptor_path, receptor_lines)
    _write_pdb(partner_path, partner_lines)

    monkeypatch.setattr(module, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        module,
        "RECEPTOR_STATES",
        {
            "demo_state": {
                "pdb": "input/receptors/demo_receptor.pdb",
                "source_chain": "X",
                "description": "demo receptor",
            }
        },
    )
    monkeypatch.setattr(module, "KINASE_DOMAIN_START", 699)
    monkeypatch.setattr(module, "KINASE_DOMAIN_END", 700)
    monkeypatch.setattr(module, "NLOBE_CLOBE_BOUNDARY", 700)
    monkeypatch.setattr(module, "BETA_MEANDER_SOURCE", "input/PPI/demo_partner.pdb")
    monkeypatch.setattr(module, "BETA_MEANDER_START", 955)
    monkeypatch.setattr(module, "BETA_MEANDER_END", 956)
    monkeypatch.setattr(
        module,
        "SHEET_DEFINITIONS",
        {
            "sheet_demo": {
                "residues": [955, 956],
                "role": "demo",
            }
        },
    )
    monkeypatch.setattr(module, "MEMBRANE_PROXIMAL", "699")

    output_dir = tmp_path / "runtime_inputs"
    result = module.prepare_phase1_inputs(output_dir=output_dir)

    receptor_output = output_dir / "receptor_demo_state.pdb"
    receptor_text = receptor_output.read_text(encoding="utf-8")
    assert " HSD " not in receptor_text
    assert " HIS " in receptor_text
    assert " CD1" in receptor_text
    assert " HN " not in receptor_text

    receptor_rows = list(csv.DictReader((output_dir / "receptor_metadata.csv").open(encoding="utf-8")))
    assert len(receptor_rows) == 1
    receptor_row = receptor_rows[0]
    assert receptor_row["rosetta_removed_hydrogens"] == "2"
    assert receptor_row["rosetta_residue_renames"] == "HSD->HIS:5"
    assert receptor_row["rosetta_atom_renames"] == "ILE:CD->CD1:1"

    report_text = (output_dir / "phase1_input_validation_report.md").read_text(encoding="utf-8")
    assert "## Rosetta Compatibility Normalization" in report_text
    assert "HSD->HIS:5" in report_text
    assert "ILE:CD->CD1:1" in report_text
    assert result["status_code"] == 0


def test_register_pilot_data_uses_historical_namespaced_fields(tmp_path: Path) -> None:
    output_path = module.register_pilot_data(tmp_path)

    rows = list(csv.DictReader(output_path.open(encoding="utf-8")))
    assert len(rows) == 2
    for row in rows:
        assert row["status"] == "historical_reference_only"
        assert row["historical_only"] == "true"
        assert row["downstream_filter"] == "historical_only == true AND status == historical_reference_only"
        assert "results_dir" not in row
        assert "runtime_inputs" not in row

    assert rows[0]["historical_reference.results_dir"]
    assert rows[1]["legacy_reference.results_dir"]
