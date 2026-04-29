from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"C:\Users\dudwn\OneDrive\문서\EGFR_MYO1D\analysis_v0_5")
PDB = ROOT / "03_receptor_pack" / "receptor" / "EGFR_plus10_final_model_EGFR_rot10_best.pdb"
ASSETS = ROOT / "reports" / "assets"


COLORS = {
    "A": (0, 168, 200),
    "B": (230, 144, 58),
    "tm": (19, 32, 38),
    "jm": (196, 145, 39),
    "kd": (88, 126, 153),
    "bg": (248, 251, 252),
    "ink": (19, 32, 38),
    "muted": (91, 112, 123),
    "membrane": (186, 225, 228),
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\malgunbd.ttf" if bold else r"C:\Windows\Fonts\malgun.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


def parse_ca(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            chain = line[21].strip() or "?"
            try:
                resi = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            rows.append((chain, resi, np.array([x, y, z], dtype=float)))
    if not rows:
        raise RuntimeError(f"No CA atoms parsed from {path}")
    return rows


def domain(resi: int) -> str:
    if 634 <= resi <= 670:
        return "tm"
    if 671 <= resi <= 704:
        return "jm"
    return "kd"


def pca_basis(points: np.ndarray) -> np.ndarray:
    centered = points - points.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return vh


def project(rows, axes: tuple[int, int], basis: np.ndarray):
    points = np.array([row[2] for row in rows])
    centered = points - points.mean(axis=0)
    transformed = centered @ basis.T
    coords = transformed[:, list(axes)]
    return coords


def fit_coords(coords: np.ndarray, width: int, height: int, margin: int):
    min_xy = coords.min(axis=0)
    max_xy = coords.max(axis=0)
    span = np.maximum(max_xy - min_xy, 1e-6)
    scale = min((width - margin * 2) / span[0], (height - margin * 2) / span[1])
    fitted = (coords - min_xy) * scale
    fitted[:, 0] += margin + (width - margin * 2 - span[0] * scale) / 2
    fitted[:, 1] = height - margin - fitted[:, 1] - (height - margin * 2 - span[1] * scale) / 2
    return fitted


def draw_trace(draw: ImageDraw.ImageDraw, screen_rows):
    by_chain: dict[str, list] = {}
    for item in screen_rows:
        by_chain.setdefault(item["chain"], []).append(item)

    for chain, items in by_chain.items():
        items.sort(key=lambda row: row["resi"])
        for left, right in zip(items, items[1:]):
            color = COLORS["A"] if chain == "A" else COLORS["B"]
            draw.line([left["xy"], right["xy"]], fill=color, width=5)
        for row in items:
            d = domain(row["resi"])
            color = COLORS[d] if d in {"tm", "jm"} else (COLORS["A"] if chain == "A" else COLORS["B"])
            r = 4 if d != "kd" else 3
            x, y = row["xy"]
            draw.ellipse((x - r, y - r, x + r, y + r), fill=color, outline=(255, 255, 255), width=1)


def draw_legend(draw: ImageDraw.ImageDraw, x: int, y: int):
    entries = [
        ("Chain A", COLORS["A"]),
        ("Chain B", COLORS["B"]),
        ("TM 634-670", COLORS["tm"]),
        ("JM 671-704", COLORS["jm"]),
        ("KD 705+", COLORS["kd"]),
    ]
    f = font(22)
    for i, (label, color) in enumerate(entries):
        yy = y + i * 34
        draw.rounded_rectangle((x, yy, x + 24, yy + 18), radius=6, fill=color)
        draw.text((x + 34, yy - 4), label, fill=COLORS["muted"], font=f)


def render(rows, axes: tuple[int, int], title: str, subtitle: str, output: Path, side: bool = False):
    width, height = 1400, 900
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.rectangle((0, 0, width, height), fill=(255, 255, 255))
    draw.text((70, 54), title, fill=COLORS["ink"], font=font(44, True))
    draw.text((70, 110), subtitle, fill=COLORS["muted"], font=font(24))

    plot_box = (70, 165, 1010, 820)
    draw.rounded_rectangle(plot_box, radius=18, fill=COLORS["bg"])

    if side:
        y_mid = 675
        draw.rectangle((plot_box[0] + 26, y_mid - 48, plot_box[2] - 26, y_mid + 48), fill=COLORS["membrane"])
        draw.line((plot_box[0] + 26, y_mid, plot_box[2] - 26, y_mid), fill=(130, 191, 197), width=3)
        draw.text((plot_box[0] + 36, y_mid + 56), "approx. membrane plane / TM anchor zone", fill=COLORS["muted"], font=font(18))

    points = np.array([r[2] for r in rows])
    basis = pca_basis(points)
    coords = project(rows, axes, basis)
    fitted = fit_coords(coords, plot_box[2] - plot_box[0] - 70, plot_box[3] - plot_box[1] - 70, 30)
    fitted[:, 0] += plot_box[0] + 35
    fitted[:, 1] += plot_box[1] + 35

    screen_rows = []
    for (chain, resi, _), xy in zip(rows, fitted):
        screen_rows.append({"chain": chain, "resi": resi, "xy": (float(xy[0]), float(xy[1]))})
    draw_trace(draw, screen_rows)

    draw.rounded_rectangle((1045, 165, 1330, 465), radius=18, fill=COLORS["bg"])
    draw.text((1074, 198), "Data source", fill=COLORS["ink"], font=font(24, True))
    source_lines = [
        "EGFR_plus10_final_model",
        "_EGFR_rot10_best.pdb",
        "CA trace only",
        "domain ranges:",
        "TM 634-670",
        "JM 671-704",
        "KD 705+",
    ]
    for i, line in enumerate(source_lines):
        draw.text((1074, 244 + i * 27), line, fill=COLORS["muted"], font=font(18))

    draw_legend(draw, 1050, 520)
    draw.text((70, 842), "Generated from local PDB coordinates for presentation overview; not a publication-quality molecular rendering.", fill=COLORS["muted"], font=font(18))
    output.parent.mkdir(parents=True, exist_ok=True)
    img.save(output)


def main():
    rows = parse_ca(PDB)
    render(
        rows,
        axes=(0, 2),
        title="+10 EGFR receptor overview: side projection",
        subtitle="Actual CA-coordinate projection from the selected final model PDB",
        output=ASSETS / "slide6_plus10_side_view.png",
        side=True,
    )
    render(
        rows,
        axes=(0, 1),
        title="+10 EGFR receptor overview: top projection",
        subtitle="Actual CA-coordinate projection highlighting protomer/domain organization",
        output=ASSETS / "slide6_plus10_top_view.png",
        side=False,
    )


if __name__ == "__main__":
    main()
