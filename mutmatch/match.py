"""Filtrage, fusion des paires R1/R2 et correspondance avec les mutations connues."""

from __future__ import annotations

from dataclasses import dataclass, field

from .vcf import Variant, VcfFile, norm_chrom
from .references import KnownMutation


@dataclass
class FilterConfig:
    """Seuils de filtrage appliques aux variants VCF."""

    require_pass: bool = True          # ne garder que FILTER == PASS
    min_vaf: float = 0.0               # VAF minimale (ratio)
    min_m: int = 0                     # nombre minimal de lectures mutees
    min_wt: int = 0                    # nombre minimal de lectures sauvages
    only_dup: bool = False             # ne garder que les evenements DUP (ITD)
    require_both_pairs: bool = False   # variant present dans R1 ET R2

    def describe(self) -> str:
        parts = []
        if self.require_pass:
            parts.append("FILTER=PASS")
        if self.min_vaf > 0:
            parts.append("VAF>=%g" % self.min_vaf)
        if self.min_m > 0:
            parts.append("M>=%d" % self.min_m)
        if self.min_wt > 0:
            parts.append("WT>=%d" % self.min_wt)
        if self.only_dup:
            parts.append("DUP uniquement")
        if self.require_both_pairs:
            parts.append("present dans R1 et R2")
        return ", ".join(parts) if parts else "aucun"


def passes_filters(v: Variant, cfg: FilterConfig) -> bool:
    """Teste un variant individuel (hors regle inter-paires)."""
    if cfg.require_pass and v.filt not in ("PASS", ".", ""):
        # '.' et '' sont toleres car certains VCF ne renseignent pas FILTER.
        if v.filt != "PASS":
            return False
    if cfg.only_dup and not v.is_dup:
        return False
    vaf = v.vaf if isinstance(v.vaf, (int, float)) else 0.0
    if vaf < cfg.min_vaf:
        return False
    m = v.m if isinstance(v.m, (int, float)) else 0
    if m < cfg.min_m:
        return False
    wt = v.wt if isinstance(v.wt, (int, float)) else 0
    if wt < cfg.min_wt:
        return False
    return True


@dataclass
class MergedVariant:
    """Un variant unique pour un echantillon, fusionnant les deux paires."""

    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str
    svlen: object = None
    is_dup: bool = False
    # metriques par paire : pair -> {"m":..,"wt":..,"vaf":..,"filter":..}
    per_pair: dict = field(default_factory=dict)
    # correspondance avec le referentiel (rempli plus tard)
    known: bool = False
    match_level: str = ""              # "exact", "position" ou ""
    known_refs: list = field(default_factory=list)

    def key(self) -> tuple:
        return (norm_chrom(self.chrom), self.pos, self.ref.upper(), self.alt.upper())

    def pos_key(self) -> tuple:
        return (norm_chrom(self.chrom), self.pos)

    @property
    def found_in_pairs(self) -> str:
        return ",".join(sorted(self.per_pair.keys()))

    @property
    def n_pairs(self) -> int:
        return len(self.per_pair)

    @property
    def max_vaf(self) -> float:
        vals = [d["vaf"] for d in self.per_pair.values() if isinstance(d["vaf"], (int, float))]
        return max(vals) if vals else 0.0

    @property
    def total_m(self) -> int:
        vals = [d["m"] for d in self.per_pair.values() if isinstance(d["m"], (int, float))]
        return sum(vals) if vals else 0


def merge_sample(vcf_files: list[VcfFile], cfg: FilterConfig) -> list[MergedVariant]:
    """Fusionne les variants (filtres) de plusieurs fichiers d'un meme echantillon."""
    merged: dict[tuple, MergedVariant] = {}
    for vf in vcf_files:
        for v in vf.variants:
            if not passes_filters(v, cfg):
                continue
            k = v.key()
            mv = merged.get(k)
            if mv is None:
                mv = MergedVariant(
                    sample=v.sample,
                    chrom=v.chrom,
                    pos=v.pos,
                    ref=v.ref,
                    alt=v.alt,
                    svlen=v.svlen,
                    is_dup=v.is_dup,
                )
                merged[k] = mv
            pair = v.pair or "?"
            mv.per_pair[pair] = {
                "m": v.m,
                "wt": v.wt,
                "vaf": v.vaf,
                "filter": v.filt,
            }
            mv.is_dup = mv.is_dup or v.is_dup

    result = list(merged.values())
    if cfg.require_both_pairs:
        result = [mv for mv in result if len(mv.per_pair) >= 2]
    # tri par chromosome puis position
    result.sort(key=lambda m: (norm_chrom(m.chrom), m.pos))
    return result


def group_by_sample(vcf_files: list[VcfFile]) -> dict[str, list[VcfFile]]:
    groups: dict[str, list[VcfFile]] = {}
    for vf in vcf_files:
        groups.setdefault(vf.sample, []).append(vf)
    return groups


def merge_all(vcf_files: list[VcfFile], cfg: FilterConfig) -> dict[str, list[MergedVariant]]:
    """Renvoie {echantillon: [MergedVariant, ...]} apres filtrage/fusion."""
    out: dict[str, list[MergedVariant]] = {}
    for sample, files in group_by_sample(vcf_files).items():
        out[sample] = merge_sample(files, cfg)
    return out


@dataclass
class KnownMatch:
    """Resultat de correspondance pour une mutation connue."""

    known: KnownMutation
    match_level: str = "none"            # "exact", "position", "none"
    matched_samples: list = field(default_factory=list)
    matched_variants: list = field(default_factory=list)  # MergedVariant


def match_known(
    per_sample: dict[str, list[MergedVariant]],
    known: list[KnownMutation],
    restrict_to_sample: bool = False,
) -> list[KnownMatch]:
    """Confronte les mutations connues aux variants detectes.

    - remplit `known`/`match_level`/`known_refs` sur les MergedVariant
    - renvoie la liste des KnownMatch (une par mutation connue)

    Si `restrict_to_sample` est vrai, une mutation connue n'est cherchee que
    dans l'echantillon dont le sample_id correspond.
    """
    # index des variants detectes
    by_key: dict[tuple, list[MergedVariant]] = {}
    by_pos: dict[tuple, list[MergedVariant]] = {}
    for variants in per_sample.values():
        for mv in variants:
            by_key.setdefault(mv.key(), []).append(mv)
            by_pos.setdefault(mv.pos_key(), []).append(mv)

    def sample_ok(mv: MergedVariant, km: KnownMutation) -> bool:
        if not restrict_to_sample:
            return True
        want = _norm_sample(km.sample_id) or _norm_sample(km.patient_id)
        if not want:
            return True
        return _norm_sample(mv.sample) == want

    results: list[KnownMatch] = []
    for km in known:
        res = KnownMatch(known=km)
        full = km.key()
        candidates: list[MergedVariant] = []
        level = "none"

        if full is not None and full in by_key:
            candidates = [mv for mv in by_key[full] if sample_ok(mv, km)]
            if candidates:
                level = "exact"

        if not candidates:
            pk = km.pos_key()
            if pk is not None and pk in by_pos:
                candidates = [mv for mv in by_pos[pk] if sample_ok(mv, km)]
                if candidates:
                    level = "position"

        res.match_level = level
        res.matched_variants = candidates
        res.matched_samples = sorted({mv.sample for mv in candidates})
        for mv in candidates:
            mv.known = True
            mv.known_refs.append(km)
            # on garde le niveau le plus fort
            if mv.match_level != "exact":
                mv.match_level = level
        results.append(res)

    return results


def _norm_sample(s: str) -> str:
    return str(s).strip().lower().replace("-", "_") if s else ""
