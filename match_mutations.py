#!/usr/bin/env python3
"""Outil en ligne de commande.

Lit les fichiers VCF de resultats (JB_*_*_fast.gz.results.vcf), applique des
filtres, fusionne les deux paires de lecture par echantillon, puis compare les
variants detectes aux mutations connues des patients (Fichier_CHU.ods,
Results_patientJB.ods).

Exemple :

    python match_mutations.py \\
        --vcf-dir /data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/fastq \\
        --reference Fichier_CHU.ods --reference Results_patientJB.ods \\
        --out resultats \\
        --min-vaf 0.001 --min-m 3

Aucune dependance externe : bibliotheque standard Python 3 uniquement.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

from mutmatch.vcf import read_vcf
from mutmatch.references import read_any
from mutmatch.match import FilterConfig, merge_all, match_known
from mutmatch.report import write_outputs


def find_vcf_files(vcf_dir: str, pattern: str) -> list[str]:
    files = sorted(glob.glob(os.path.join(vcf_dir, pattern)))
    return files


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Correspondance des variants VCF FLT3 avec les mutations connues.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--vcf-dir", help="Dossier contenant les fichiers VCF de resultats.")
    p.add_argument("--vcf", action="append", default=[],
                   help="Fichier VCF explicite (repetable). Complete --vcf-dir.")
    p.add_argument("--pattern", default="*results.vcf*",
                   help="Motif de recherche des VCF dans --vcf-dir "
                        "(gere les formes _fast et .fast, compressees ou non).")
    p.add_argument("--reference", action="append", default=[],
                   help="Fichier de reference .ods/.csv/.tsv (repetable).")
    p.add_argument("--out", default="resultats",
                   help="Dossier de sortie pour les rapports.")
    p.add_argument("--ods-output", action="store_true",
                   help="Ecrire aussi les rapports au format .ods.")

    g = p.add_argument_group("filtres")
    g.add_argument("--min-vaf", type=float, default=0.0,
                   help="VAF minimale (ratio) pour conserver un variant.")
    g.add_argument("--min-m", type=int, default=0,
                   help="Nombre minimal de lectures mutees (M).")
    g.add_argument("--min-wt", type=int, default=0,
                   help="Nombre minimal de lectures sauvages (WT).")
    g.add_argument("--only-dup", action="store_true",
                   help="Ne garder que les evenements DUP (duplications/ITD).")
    g.add_argument("--require-both-pairs", action="store_true",
                   help="Ne garder que les variants presents dans R1 ET R2.")
    g.add_argument("--no-require-pass", action="store_true",
                   help="Ne pas exiger FILTER=PASS.")

    t = p.add_argument_group("triage des variants nouveaux")
    t.add_argument("--candidat-vaf-min", type=float, default=0.01,
                   help="VAF minimale pour qu'un variant nouveau soit un "
                        "'candidat ITD' plutot que du bruit.")
    t.add_argument("--artefact-nb-echantillons", type=int, default=5,
                   help="Un variant present dans au moins ce nombre "
                        "d'echantillons est classe 'artefact probable'.")

    p.add_argument("--no-restrict-to-sample", action="store_true",
                   help="Chercher chaque mutation connue dans TOUS les "
                        "echantillons (par defaut : uniquement dans "
                        "l'echantillon dont le sample_id/patient correspond, "
                        "ce qui evite les faux positifs inter-patients).")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    vcf_paths = list(args.vcf)
    if args.vcf_dir:
        vcf_paths += find_vcf_files(args.vcf_dir, args.pattern)
    vcf_paths = sorted(set(vcf_paths))

    if not vcf_paths:
        print("Aucun fichier VCF trouve. Utilisez --vcf-dir ou --vcf.",
              file=sys.stderr)
        return 2

    cfg = FilterConfig(
        require_pass=not args.no_require_pass,
        min_vaf=args.min_vaf,
        min_m=args.min_m,
        min_wt=args.min_wt,
        only_dup=args.only_dup,
        require_both_pairs=args.require_both_pairs,
    )

    print("Lecture de %d fichier(s) VCF..." % len(vcf_paths))
    vcf_files = []
    for path in vcf_paths:
        try:
            vf = read_vcf(path)
        except Exception as exc:  # noqa: BLE001 - on continue sur les autres
            print("  ! echec de lecture %s : %s" % (path, exc), file=sys.stderr)
            continue
        vcf_files.append(vf)
        print("  %-40s echantillon=%s paire=%s variants=%d"
              % (os.path.basename(path), vf.sample, vf.pair or "?", len(vf.variants)))

    per_sample = merge_all(vcf_files, cfg)
    n_variants = sum(len(v) for v in per_sample.values())
    print("\nFiltres appliques : %s" % cfg.describe())
    print("%d echantillon(s), %d variant(s) apres filtrage/fusion."
          % (len(per_sample), n_variants))

    known = []       # mutations avec coordonnees (Fichier_CHU) -> matching VCF
    clinical = []    # mutations cliniques HGVS par patient (Results_patientJB)
    for ref in args.reference:
        try:
            kind, items = read_any(ref)
        except Exception as exc:  # noqa: BLE001
            print("  ! echec de lecture reference %s : %s" % (ref, exc),
                  file=sys.stderr)
            continue
        if kind == "clinical":
            clinical += items
            print("Reference %s (clinique) : %d mutation(s) sur %d patient(s)."
                  % (os.path.basename(ref), len(items),
                     len({m.sample_id for m in items})))
        else:
            known += items
            print("Reference %s (coordonnees) : %d mutation(s) connue(s)."
                  % (os.path.basename(ref), len(items)))

    matches = match_known(per_sample, known,
                          restrict_to_sample=not args.no_restrict_to_sample)

    if known:
        n_detected = sum(1 for m in matches if m.match_level != "none")
        print("\nMutations connues (coordonnees) detectees : %d / %d"
              % (n_detected, len(matches)))
        n_novel = sum(1 for v in _flat(per_sample) if not v.known)
        print("Variants detectes non presents dans les references : %d" % n_novel)
    if clinical:
        vcf_chroms = sorted({v.chrom for v in _flat(per_sample)})
        print("Chromosome(s) couvert(s) par les VCF : %s"
              % (", ".join(vcf_chroms) if vcf_chroms else "aucun"))

    written = write_outputs(args.out, per_sample, matches,
                            clinical=clinical, also_ods=args.ods_output,
                            candidate_min_vaf=args.candidat_vaf_min,
                            artifact_min_samples=args.artefact_nb_echantillons)
    print("\nRapports ecrits dans %s/ :" % args.out)
    for w in written:
        print("  %s" % w)

    return 0


def _flat(per_sample):
    for variants in per_sample.values():
        yield from variants


if __name__ == "__main__":
    raise SystemExit(main())
