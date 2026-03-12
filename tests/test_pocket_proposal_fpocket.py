from egfr_pipeline.phase2.pocket_proposal import parse_fpocket_output


def test_fpocket_zero_padded_pocket_file_is_parsed(tmp_path):
    state = "3GT8_raw"
    state_dir = tmp_path / state
    fpocket_dir = state_dir / f"receptor_{state}_out"
    pockets_dir = fpocket_dir / "pockets"
    pockets_dir.mkdir(parents=True)

    (fpocket_dir / f"receptor_{state}_info.txt").write_text(
        "\n".join(
            [
                "Pocket 1 :",
                "    Score :                0.4521",
                "    Druggability Score :   0.2340",
                "    Number of Alpha Spheres :    35",
                "    Volume :               456.78",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (pockets_dir / "pocket0001_atm.pdb").write_text(
        "\n".join(
            [
                "ATOM      1  C   ALA A 123      10.000  20.000  30.000  1.00 20.00           C",
                "ATOM      2  C   GLY A 124      12.000  22.000  32.000  1.00 20.00           C",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    pockets = parse_fpocket_output(state, state_dir)

    assert len(pockets) == 1
    assert pockets[0]["residue_ids"] == "ALA123;GLY124"
    assert pockets[0]["centroid_x"] == "11.000"
    assert pockets[0]["centroid_y"] == "21.000"
    assert pockets[0]["centroid_z"] == "31.000"
