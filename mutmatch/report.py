"""Generation des rapports de sortie (TSV, et .ods optionnel)."""

from __future__ import annotations

import csv
import os

from . import ods
from .match import MergedVariant, KnownMatch
from .references import ClinicalMutation, norm_hgvs
from .vcf import norm_chrom


def _norm_sample(s: str) -> str:
    return str(s).strip().lower().replace("-", "_") if s else ""


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


def synthesis_table(clinical: list[ClinicalMutation],
                    matches: list[KnownMatch],
                    vcf_chroms: set) -> list[list[str]]:
    """Synthese par patient : mutations cliniques (HGVS) x detection VCF.

    Les mutations cliniques (Results_patientJB, sans coordonnees) sont croisees
    avec les coordonnees + le statut de detection issus du Fichier_CHU, par
    (echantillon, GENE:c.xxx).
    """
    # Table de correspondance : (echantillon, hgvs) -> infos de detection.
    lut = {}
    for km in matches:
        k = km.known
        sample = k.sample_id or k.patient_id
        key = (_norm_sample(sample), norm_hgvs(k.query))
        vaf = max([mv.max_vaf for mv in km.matched_variants], default=None)
        mtot = (sum(mv.total_m for mv in km.matched_variants)
                if km.matched_variants else None)
        lut[key] = {
            "level": km.match_level, "chrom": k.chrom, "pos": k.pos,
            "ref": k.ref, "alt": k.alt, "vaf": vaf, "m": mtot,
        }

    header = [
        "sample", "n_echantillon", "caryotype", "gene",
        "hgvs_c", "hgvs_p", "chrom", "pos", "ref", "alt",
        "couvert_par_ce_vcf", "detecte", "niveau_correspondance",
        "VAF_max", "M_total",
    ]
    rows = [header]
    for c in sorted(clinical, key=lambda m: (_norm_sample(m.sample_id), m.gene)):
        info = lut.get((_norm_sample(c.sample_id), norm_hgvs(c.hgvs)))
        if info is not None:
            covered = "oui" if norm_chrom(info["chrom"]) in vcf_chroms else "non"
            detecte = "oui" if info["level"] != "none" else "non"
            chrom, pos, ref, alt = info["chrom"], info["pos"], info["ref"], info["alt"]
            level, vaf, mtot = info["level"], info["vaf"], info["m"]
        else:
            # mutation clinique sans coordonnee correspondante dans le Fichier_CHU
            covered = "inconnu (pas de coordonnee)"
            detecte = "?"
            chrom = pos = ref = alt = ""
            level = ""
            vaf = mtot = ""
        rows.append([_fmt(x) for x in [
            c.sample_id, c.n_echantillon, c.caryotype, c.gene,
            c.hgvs_c, c.hgvs_p, chrom, pos, ref, alt,
            covered, detecte, level, vaf, mtot,
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
                  clinical: list[ClinicalMutation] | None = None,
                  also_ods: bool = False) -> list[str]:
    """Ecrit les tables dans `outdir`. Renvoie la liste des fichiers crees."""
    os.makedirs(outdir, exist_ok=True)
    tables = {
        "variants_par_echantillon": variants_table(per_sample),
        "correspondance_mutations_connues": known_matches_table(matches),
        "variants_nouveaux": novel_table(per_sample),
    }
    if clinical:
        vcf_chroms = {norm_chrom(mv.chrom)
                      for variants in per_sample.values() for mv in variants}
        tables["synthese_par_patient"] = synthesis_table(
            clinical, matches, vcf_chroms)
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
