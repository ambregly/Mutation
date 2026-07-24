#!/usr/bin/env python3
"""Carte des duplications FLT3-ITD : positions (axe Y) x patients (axe X).

Genere une figure SVG (aucune dependance) a partir de la table
`variants_par_echantillon.tsv` produite par match_mutations.py.

- axe Y (a gauche) : positions du chromosome 13, de la plus basse (en HAUT) a la
  plus haute (en BAS) ;
- axe X (en HAUT)  : les patients (JB_01 ... JB_20) ;
- chaque duplication est un segment [pos, pos+svlen] avec un point a son debut ;
- les ITD **detectees dans les references** (connu=oui) sont mises en lumiere
  (rouge, epais) et leur etendue est surlignee dans la colonne du patient, ce
  qui permet de voir si les nouveaux variants (bleu) sont des sous-groupes de la
  duplication principale.

Exemple :
    python3 plot_itd.py --input resultats/variants_par_echantillon.tsv \\
        --out resultats/carte_itd.svg
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import re
import shutil
import subprocess


def _read(path: str, dup_only: bool):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    h = {name: i for i, name in enumerate(rows[0])}
    out = []
    for r in rows[1:]:
        if not r or len(r) <= h["pos"]:
            continue
        if dup_only and r[h["DUP"]] != "oui":
            continue
        try:
            pos = int(r[h["pos"]])
            svlen = abs(int(r[h["svlen"]])) if r[h["svlen"]] else 0
        except ValueError:
            continue
        try:
            vaf = float(r[h["VAF_max"]]) if r[h["VAF_max"]] else 0.0
        except ValueError:
            vaf = 0.0
        try:
            m = int(r[h["M_total"]]) if h.get("M_total") is not None and \
                r[h["M_total"]] else 0
        except (ValueError, KeyError):
            m = 0
        out.append({
            "sample": r[h["sample"]],
            "pos": pos,
            "svlen": svlen,
            "vaf": vaf,
            "m": m,
            "ref": r[h["ref"]] if h.get("ref") is not None else "",
            "alt": r[h["alt"]] if h.get("alt") is not None else "",
            "known": r[h["connu"]] == "oui",
        })
    return out


def read_exons(path: str) -> list[dict]:
    """Lit les regions d'exons a colorier.

    Accepte :
      * un FASTA (ex. la reference FiLT3r sequence-FLT3-Ex13-14-15-20.fa) : les
        coordonnees sont extraites des en-tetes '>...13:28033760-28034429:-1...' ;
      * un BED / TSV : colonnes chrom, start, end, [nom], [couleur].

    Renvoie une liste de {chrom, start, end, name, color?}.
    """
    exons: list[dict] = []
    with open(path, encoding="utf-8") as fh:
        first = fh.read(1)
        fh.seek(0)
        if first == ">":  # FASTA
            coord = re.compile(r"(?:chr)?(\w+)\s*:\s*(\d+)\s*-\s*(\d+)")
            for line in fh:
                if not line.startswith(">"):
                    continue
                m = coord.search(line)
                if not m:
                    continue
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                exons.append({"chrom": m.group(1),
                              "start": int(m.group(2)), "end": int(m.group(3)),
                              "name": name})
        else:  # BED / TSV
            for line in fh:
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                cols = line.split("\t") if "\t" in line else line.split()
                if len(cols) < 3:
                    continue
                try:
                    start, end = int(cols[1]), int(cols[2])
                except ValueError:
                    continue
                e = {"chrom": cols[0], "start": min(start, end),
                     "end": max(start, end),
                     "name": cols[3] if len(cols) > 3 else ""}
                if len(cols) > 4 and cols[4].startswith("#"):
                    e["color"] = cols[4]
                exons.append(e)
    return exons


def _radius(value, size_by="vaf") -> float:
    """Rayon d'un point selon la metrique choisie (meme formule partout).

    - 'vaf' : proportionnel a la VAF (ratio) ;
    - 'm'   : proportionnel au log10 du nombre de reads (M).
    """
    if size_by == "m":
        return 2.0 + min(6.0, 1.7 * math.log10(float(value or 0.0) + 1.0))
    return 2.0 + min(6.0, 260.0 * float(value or 0.0))


def _lanes(variants, gap_px, ytop):
    """Repartit les variants d'un patient en 'couloirs' pour eviter le chevauchement."""
    ordered = sorted(variants, key=lambda v: v["pos"])
    lane_last = []  # derniere position basse (y) de chaque couloir
    for v in ordered:
        y1 = ytop(v["pos"])
        y2 = ytop(v["pos"] + v["svlen"])
        placed = False
        for li in range(len(lane_last)):
            if lane_last[li] <= y1 - gap_px:
                v["lane"] = li
                lane_last[li] = max(y2, y1 + 2)
                placed = True
                break
        if not placed:
            v["lane"] = len(lane_last)
            lane_last.append(max(y2, y1 + 2))
    n_lanes = max((v["lane"] for v in ordered), default=0) + 1
    return ordered, n_lanes


def build_svg(variants, dup_only=True, size_by="vaf", exons=None):
    samples = sorted({v["sample"] for v in variants})
    if not variants:
        raise SystemExit("Aucune duplication a tracer.")

    for v in variants:
        v["cat"] = "known" if v["known"] else "new"

    pos_lo = min(v["pos"] for v in variants)
    pos_hi = max(v["pos"] + v["svlen"] for v in variants)
    span = max(pos_hi - pos_lo, 1)
    pad = span * 0.04
    pos_lo -= pad
    pos_hi += pad
    span = pos_hi - pos_lo

    # Geometrie
    ML, MR, MT, MB = 104, 300, 104, 46
    col_w = 56
    plot_w = col_w * len(samples)
    plot_h = 660
    W = ML + plot_w + MR
    H = MT + plot_h + MB

    def y(pos):
        return MT + (pos - pos_lo) / span * plot_h

    def col_x(i):
        return ML + i * col_w

    C_KNOWN = "#E4572E"      # ITD de reference (mise en lumiere)
    C_NEW = "#3A7BD5"        # nouveau variant
    C_AXIS = "#444"
    C_GRID = "#e6e6e6"
    C_BAND = "#E4572E"
    STYLE = {
        "known": (C_KNOWN, 4.0, 1.0),
        "new": (C_NEW, 2.0, 0.8),
    }

    s = []
    s.append(
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
        'viewBox="0 0 %d %d" font-family="Helvetica,Arial,sans-serif">' % (W, H, W, H))
    s.append('<rect width="%d" height="%d" fill="white"/>' % (W, H))

    # Titre
    s.append('<text x="%d" y="34" font-size="20" font-weight="bold" fill="#222">'
             'Carte des duplications détectées par FiLT3r</text>' % ML)
    metric_label = "nombre de reads (M)" if size_by == "m" else "VAF"
    s.append('<text x="%d" y="56" font-size="13" fill="#666">'
             'position chr13 (basse en haut, haute en bas) x patient ; '
             'rouge = ITD de reference, bleu = nouveau variant ; '
             'taille du point ~ %s%s</text>'
             % (ML, metric_label, " ; DUP uniquement" if dup_only else ""))

    # Grille + graduations Y (positions)
    n_ticks = 8
    for t in range(n_ticks + 1):
        pos = pos_lo + span * t / n_ticks
        yy = y(pos)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" '
                 'stroke-width="1"/>' % (ML, yy, ML + plot_w, yy, C_GRID))
        s.append('<text x="%d" y="%.1f" font-size="11" fill="%s" '
                 'text-anchor="end">%d</text>'
                 % (ML - 8, yy + 4, C_AXIS, round(pos)))
    # Axe Y
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
             'stroke-width="1.5"/>' % (ML, MT, ML, MT + plot_h, C_AXIS))
    s.append('<text transform="translate(24,%d) rotate(-90)" font-size="13" '
             'fill="#333" text-anchor="middle">Position chr13</text>'
             % (MT + plot_h / 2))

    # Colonnes + labels patients (en HAUT)
    for i, sample in enumerate(samples):
        x0 = col_x(i)
        if i % 2 == 0:
            s.append('<rect x="%.1f" y="%d" width="%d" height="%d" '
                     'fill="#fafafa"/>' % (x0, MT, col_w, plot_h))
        cx = x0 + col_w / 2
        s.append('<text x="%.1f" y="%d" font-size="12" fill="#333" '
                 'text-anchor="start" transform="rotate(-55 %.1f %d)">%s</text>'
                 % (cx, MT - 10, cx, MT - 10, html.escape(sample)))

    # Grille verticale : une ligne a chaque separation de patient
    for i in range(len(samples) + 1):
        xx = col_x(i)
        s.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" '
                 'stroke-width="1"/>' % (xx, MT, xx, MT + plot_h, "#d0d0d0"))
    # Cadre de la zone de trace
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
             'stroke-width="1.2"/>' % (ML, MT, ML + plot_w, MT, C_AXIS))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
             'stroke-width="1.2"/>'
             % (ML, MT + plot_h, ML + plot_w, MT + plot_h, C_AXIS))

    # Bandes d'exons (colorees) sur toute la largeur, si fournies.
    exon_palette = ["#F6C177", "#9CCFD8", "#C4A7E7", "#A3BE8C",
                    "#EBBCBA", "#B8C0E0"]
    drawn_exons = []
    for k, ex in enumerate(exons or []):
        y1, y2 = y(ex["start"]), y(ex["end"])
        top, bot = min(y1, y2), max(y1, y2)
        # ignore les exons entierement hors de la fenetre affichee
        if bot < MT or top > MT + plot_h:
            continue
        top = max(top, MT)
        bot = min(bot, MT + plot_h)
        color = ex.get("color") or exon_palette[k % len(exon_palette)]
        s.append('<rect x="%d" y="%.1f" width="%d" height="%.1f" fill="%s" '
                 'opacity="0.22"/>' % (ML, top, plot_w, max(bot - top, 1.5), color))
        label = ex.get("name") or ("exon %d" % (k + 1))
        s.append('<text x="%d" y="%.1f" font-size="11" fill="#333" '
                 'text-anchor="end">%s</text>'
                 % (ML + plot_w + 6, top + 12, html.escape(label)))
        drawn_exons.append((label, color))

    # Bandes de surlignage : etendue de l'ITD connu de chaque patient
    for i, sample in enumerate(samples):
        x0 = col_x(i)
        for v in variants:
            if v["sample"] == sample and v["known"]:
                y1, y2 = y(v["pos"]), y(v["pos"] + v["svlen"])
                s.append('<rect x="%.1f" y="%.1f" width="%d" height="%.1f" '
                         'fill="%s" opacity="0.13"/>'
                         % (x0 + 3, min(y1, y2), col_w - 6, abs(y2 - y1) + 2, C_BAND))

    # Variants (segments + points), repartis en couloirs par patient
    for i, sample in enumerate(samples):
        x0 = col_x(i)
        vs = [v for v in variants if v["sample"] == sample]
        ordered, n_lanes = _lanes(vs, gap_px=2, ytop=y)
        inner = col_w - 16
        for v in ordered:
            if n_lanes <= 1:
                x = x0 + col_w / 2
            else:
                x = x0 + 8 + inner * v["lane"] / (n_lanes - 1)
            y1, y2 = y(v["pos"]), y(v["pos"] + v["svlen"])
            r = _radius(v[size_by], size_by)
            col, wdt, op = STYLE[v["cat"]]
            rr = max(r, 4.5) if v["cat"] == "known" else r
            s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                     'stroke-width="%.1f" opacity="%.2f" stroke-linecap="round"/>'
                     % (x, y1, x, y2, col, wdt, op))
            s.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" '
                     'opacity="%.2f"/>' % (x, y1, rr, col, min(op + 0.15, 1.0)))

    # Legende
    lx = ML + plot_w + 26
    ly = MT + 6
    s.append('<text x="%d" y="%d" font-size="13" font-weight="bold" '
             'fill="#333">Legende</text>' % (lx, ly))
    items = [
        (C_KNOWN, 4.0, 5.0, "ITD de reference (connu=oui)"),
        (C_NEW, 2.0, 3.5, "nouveau variant (connu=non)"),
    ]
    for j, (col, wdt, rr, label) in enumerate(items):
        yy = ly + 24 + j * 24
        s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
                 'stroke-width="%.1f" stroke-linecap="round"/>'
                 % (lx, yy, lx + 26, yy, col, wdt))
        s.append('<circle cx="%d" cy="%d" r="%.1f" fill="%s"/>'
                 % (lx + 13, yy, rr, col))
        s.append('<text x="%d" y="%d" font-size="12" fill="#333">%s</text>'
                 % (lx + 34, yy + 4, label))
    yb = ly + 24 + len(items) * 24
    s.append('<rect x="%d" y="%d" width="26" height="14" fill="%s" '
             'opacity="0.13"/>' % (lx, yb, C_BAND))
    s.append('<text x="%d" y="%d" font-size="12" fill="#333">etendue de '
             "l'ITD connu</text>" % (lx + 34, yb + 11))
    s.append('<text x="%d" y="%d" font-size="11" fill="#666">segment vertical = '
             'etendue [pos, pos+svlen]</text>' % (lx, yb + 40))
    s.append('<text x="%d" y="%d" font-size="11" fill="#666">point = debut ; '
             'taille ~ %s</text>' % (lx, yb + 58, metric_label))
    s.append('<text x="%d" y="%d" font-size="11" fill="#666">Un variant proche '
             "(bleu) de l'ITD rouge</text>" % (lx, yb + 84))
    s.append('<text x="%d" y="%d" font-size="11" fill="#666">= sous-groupe '
             'probable de l\'ITD.</text>' % (lx, yb + 100))

    # Bandes d'exons dans la legende
    ys = yb + 132
    if drawn_exons:
        s.append('<text x="%d" y="%d" font-size="12" font-weight="bold" '
                 'fill="#333">Exons (region)</text>' % (lx, ys))
        ys += 22
        for label, color in drawn_exons:
            s.append('<rect x="%d" y="%d" width="26" height="13" fill="%s" '
                     'opacity="0.5"/>' % (lx, ys - 11, color))
            s.append('<text x="%d" y="%d" font-size="11" fill="#333">%s</text>'
                     % (lx + 34, ys, html.escape(label)))
            ys += 22
        ys += 12

    # Echelle de taille (VAF ou M)
    s.append('<text x="%d" y="%d" font-size="12" font-weight="bold" '
             'fill="#333">Taille du point = %s</text>' % (lx, ys, metric_label))
    if size_by == "m":
        scale = [(10, "10 reads"), (50, "50 reads"),
                 (200, "200 reads"), (1000, "1000 reads")]
    else:
        scale = [(0.001, "VAF 0,1 %"), (0.005, "VAF 0,5 %"),
                 (0.01, "VAF 1 %"), (0.02, "VAF 2 %")]
    cx = lx + 13
    row_y = ys + 24
    for val, label in scale:
        r = _radius(val, size_by)
        s.append('<circle cx="%d" cy="%.1f" r="%.1f" fill="%s" '
                 'opacity="0.85"/>' % (cx, row_y, r, C_NEW))
        s.append('<text x="%d" y="%.1f" font-size="11" fill="#333">'
                 '%s</text>' % (lx + 34, row_y + 4, label))
        row_y += 26

    s.append('</svg>')
    return "\n".join(s)


def svg_to_png(svg_path: str, png_path: str, scale: float = 2.0) -> str | None:
    """Convertit un SVG en PNG avec le premier outil disponible.

    Essaie dans l'ordre : cairosvg (module Python), rsvg-convert, inkscape,
    un navigateur Chromium/Chrome, puis Firefox (headless). Renvoie le nom de
    l'outil utilise, ou None si aucun n'est disponible.
    """
    svg_path = os.path.abspath(svg_path)
    png_path = os.path.abspath(png_path)

    # 1. cairosvg (pip install cairosvg)
    try:
        import cairosvg  # type: ignore
        cairosvg.svg2png(url=svg_path, write_to=png_path, scale=scale)
        return "cairosvg"
    except Exception:
        pass

    # 2. rsvg-convert (paquet librsvg2-bin)
    if shutil.which("rsvg-convert"):
        subprocess.run(["rsvg-convert", "-z", str(scale), "-o", png_path,
                        svg_path], check=True)
        return "rsvg-convert"

    # 3. inkscape
    if shutil.which("inkscape"):
        subprocess.run(["inkscape", svg_path, "--export-type=png",
                        "--export-filename=" + png_path], check=True)
        return "inkscape"

    # Dimensions du SVG (pour dimensionner la fenetre du navigateur).
    w_svg, h_svg = _svg_size(svg_path)

    # 4. navigateur Chromium / Chrome (headless)
    for exe in ("chromium", "chromium-browser", "google-chrome",
                "google-chrome-stable", "chrome", "headless_shell"):
        path = shutil.which(exe)
        if path:
            subprocess.run([path, "--headless", "--no-sandbox", "--disable-gpu",
                            "--force-device-scale-factor=%s" % scale,
                            "--screenshot=" + png_path,
                            "--window-size=%d,%d" % (w_svg + 4, h_svg + 4),
                            "file://" + svg_path], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return os.path.basename(path)

    # 5. Firefox (headless) - resolution 1x
    ff = shutil.which("firefox")
    if ff:
        subprocess.run([ff, "--headless",
                        "--window-size=%d,%d" % (w_svg + 4, h_svg + 4),
                        "--screenshot", png_path, "file://" + svg_path],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        if os.path.exists(png_path):
            return "firefox"

    return None


def _svg_size(svg_path: str) -> tuple[int, int]:
    """Lit les attributs width/height du SVG (repli 1600x900)."""
    try:
        with open(svg_path, encoding="utf-8") as fh:
            head = fh.read(600)
        import re
        w = re.search(r'width="(\d+)"', head)
        h = re.search(r'height="(\d+)"', head)
        return (int(w.group(1)) if w else 1600,
                int(h.group(1)) if h else 900)
    except Exception:
        return 1600, 900


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True,
                   help="Chemin de variants_par_echantillon.tsv")
    p.add_argument("--out", default="carte_itd.svg",
                   help="Fichier de sortie. Extension .svg ou .png (le PNG "
                        "necessite cairosvg, rsvg-convert, inkscape ou un "
                        "navigateur Chromium/Chrome).")
    p.add_argument("--png", action="store_true",
                   help="Ecrire aussi un PNG a cote du SVG.")
    p.add_argument("--scale", type=float, default=2.0,
                   help="Facteur de resolution du PNG (2 = ~2x, plus net).")
    p.add_argument("--all", action="store_true",
                   help="Tracer tous les variants (pas seulement DUP=oui)")
    p.add_argument("--size-by", choices=["vaf", "m"], default="vaf",
                   help="Metrique codee par la taille des points : 'vaf' "
                        "(defaut) ou 'm' (nombre de reads).")
    p.add_argument("--exons", metavar="FICHIER",
                   help="Fichier definissant les exons a colorier : soit la "
                        "reference FiLT3r (.fa, coordonnees lues dans les "
                        "en-tetes), soit un BED (chrom start end [nom] [couleur]).")
    args = p.parse_args(argv)

    variants = _read(args.input, dup_only=not args.all)
    exons = read_exons(args.exons) if args.exons else None
    svg = build_svg(variants, dup_only=not args.all,
                    size_by=args.size_by, exons=exons)

    # Le SVG est toujours ecrit (c'est la source). Si --out finit par .png, on
    # ecrit le SVG a cote et on rend le PNG demande.
    want_png = args.png or args.out.lower().endswith(".png")
    svg_path = args.out[:-4] + ".svg" if args.out.lower().endswith(".png") else args.out
    with open(svg_path, "w", encoding="utf-8") as fh:
        fh.write(svg)
    print("Figure ecrite : %s (%d duplications, %d patients)"
          % (svg_path, len(variants), len({v["sample"] for v in variants})))

    if want_png:
        png_path = (args.out if args.out.lower().endswith(".png")
                    else svg_path[:-4] + ".png")
        tool = svg_to_png(svg_path, png_path, scale=args.scale)
        if tool:
            print("PNG ecrit : %s (via %s)" % (png_path, tool))
        else:
            print("PNG non genere : aucun convertisseur trouve "
                  "(cairosvg, rsvg-convert, inkscape, chromium/chrome, firefox).\n"
                  "  Sans droits admin, le plus simple est un environnement "
                  "virtuel :\n"
                  "    python3 -m venv ~/venv-mut\n"
                  "    ~/venv-mut/bin/pip install cairosvg\n"
                  "    ~/venv-mut/bin/python %s --input ... --out %s\n"
                  "  Avec sudo : apt install librsvg2-bin (fournit rsvg-convert).\n"
                  "  Ou ouvrez simplement le SVG dans un navigateur :\n"
                  "    %s" % (os.path.basename(__file__), png_path, svg_path))


if __name__ == "__main__":
    main()
