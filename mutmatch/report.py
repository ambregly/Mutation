"""Generation des rapports de sortie (TSV, et .ods optionnel)."""

from __future__ import annotations

import csv
import os

from . import ods
from .match import MergedVariant, KnownMatch


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "oui" if v else "non"
    if isinstance(v, float):
        # evite la notation scientifique pour les petites VAF
        return ("%f" % v).rstrip("0").rstrip(".") if v else "0"
    return str(v)


def _pair(mv: MergedVariant, pair: str, field: str):
    d = mv.per_pair.get(pair)
    return d[field] if d else None


# --- Tables ---------------------------------------------------------------

def variants_table(per_sample: dict[str, list[MergedVariant]]) -> list[list[str]]:
    header = [
        "sample", "chrom", "pos", "ref", "alt", "svlen", "DUP",
        "trouve_dans_paires", "nb_paires",
        "M_R1", "VAF_R1", "WT_R1", "M_R2", "VAF_R2", "WT_R2",
        "M_total", "VAF_max", "connu", "niveau_correspondance",
        "query_connue", "patient_id",
    ]
    rows = [header]
    for sample in sorted(per_sample):
        for mv in per_sample[sample]:
            queries = ";".join(k.query for k in mv.known_refs if k.query)
            patients = ";".join(k.patient_id for k in mv.known_refs if k.patient_id)
            rows.append([_fmt(x) for x in [
                mv.sample, mv.chrom, mv.pos, mv.ref, mv.alt, mv.svlen, mv.is_dup,
                mv.found_in_pairs, mv.n_pairs,
                _pair(mv, "1", "m"), _pair(mv, "1", "vaf"), _pair(mv, "1", "wt"),
                _pair(mv, "2", "m"), _pair(mv, "2", "vaf"), _pair(mv, "2", "wt"),
                mv.total_m, mv.max_vaf, mv.known, mv.match_level,
                queries, patients,
            ]])
    return rows


def known_matches_table(matches: list[KnownMatch]) -> list[list[str]]:
    header = [
        "source", "query", "chrom", "pos", "ref", "alt",
        "patient_id", "sample_id", "position_brute",
        "detecte", "niveau_correspondance",
        "echantillons_trouves", "VAF_max_detectee", "M_total_detecte",
    ]
    rows = [header]
    for km in matches:
        k = km.known
        detected = "oui" if km.match_level != "none" else "non"
        vaf = max([mv.max_vaf for mv in km.matched_variants], default="")
        mtot = sum([mv.total_m for mv in km.matched_variants]) if km.matched_variants else ""
        rows.append([_fmt(x) for x in [
            k.source, k.query, k.chrom, k.pos, k.ref, k.alt,
            k.patient_id, k.sample_id, k.position_raw,
            detected, km.match_level,
            ";".join(km.matched_samples), vaf, mtot,
        ]])
    return rows


def novel_table(per_sample: dict[str, list[MergedVariant]]) -> list[list[str]]:
    header = [
        "sample", "chrom", "pos", "ref", "alt", "svlen", "DUP",
        "trouve_dans_paires", "M_total", "VAF_max",
    ]
    rows = [header]
    for sample in sorted(per_sample):
        for mv in per_sample[sample]:
            if mv.known:
                continue
            rows.append([_fmt(x) for x in [
                mv.sample, mv.chrom, mv.pos, mv.ref, mv.alt, mv.svlen, mv.is_dup,
                mv.found_in_pairs, mv.total_m, mv.max_vaf,
            ]])
    return rows


# --- Ecriture -------------------------------------------------------------

def write_tsv(path: str, rows: list[list[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(rows)


def write_outputs(outdir: str,
                  per_sample: dict[str, list[MergedVariant]],
                  matches: list[KnownMatch],
                  also_ods: bool = False) -> list[str]:
    """Ecrit les trois tables dans `outdir`. Renvoie la liste des fichiers crees."""
    os.makedirs(outdir, exist_ok=True)
    tables = {
        "variants_par_echantillon": variants_table(per_sample),
        "correspondance_mutations_connues": known_matches_table(matches),
        "variants_nouveaux": novel_table(per_sample),
    }
    written = []
    for name, rows in tables.items():
        tsv = os.path.join(outdir, name + ".tsv")
        write_tsv(tsv, rows)
        written.append(tsv)
        if also_ods:
            odspath = os.path.join(outdir, name + ".ods")
            ods.write_table(odspath, rows, table_name=name[:31])
            written.append(odspath)
    return written
