"""Lecture des fichiers VCF de resultats (JB_*_*_fast.gz.results.vcf).

Ces fichiers sont produits par l'outil de detection de duplications/ITD sur
FLT3. Chaque echantillon a deux fichiers : un par paire de lecture (R1 = _1,
R2 = _2). Le module gere les fichiers en clair et compresses (.gz).
"""

from __future__ import annotations

import gzip
import io
import os
import re
from dataclasses import dataclass, field

# Nom de fichier attendu : JB_01_1_fast.gz.results.vcf ou JB_01_1.fast.gz...
#   groupe 1 = echantillon (JB_01), groupe 2 = paire (1 ou 2)
#   le separateur avant 'fast' peut etre '_' ou '.' (les deux formes existent).
_FNAME_RE = re.compile(r"(?P<sample>.+?)_(?P<pair>[12])[._]fast", re.IGNORECASE)

# Ligne ##sample=.../JB_01_1.fastq.gz,.../JB_01_2.fastq.gz
_SAMPLE_FASTQ_RE = re.compile(r"(?P<sample>[A-Za-z0-9]+_\d+)_[12]\.fastq", re.IGNORECASE)


@dataclass
class Variant:
    """Un variant VCF avec ses champs INFO deja decodes."""

    chrom: str
    pos: int
    vid: str
    ref: str
    alt: str
    qual: str
    filt: str
    info: dict = field(default_factory=dict)

    # Renseignes lors de la lecture (contexte du fichier).
    sample: str = ""
    pair: str = ""
    source_file: str = ""

    # Champs INFO frequents, exposes en attributs pour plus de commodite.
    @property
    def svlen(self):
        return self.info.get("SVLEN")

    @property
    def m(self):
        """Nombre de lectures portant le variant."""
        return self.info.get("M")

    @property
    def wt(self):
        """Nombre de lectures sauvages (wild type)."""
        return self.info.get("WT")

    @property
    def vaf(self):
        """Variant Allele Frequency (ratio)."""
        return self.info.get("VAF")

    @property
    def is_dup(self) -> bool:
        return bool(self.info.get("DUP"))

    def key(self) -> tuple:
        """Cle d'identite d'un variant : (chrom, pos, ref, alt)."""
        return (norm_chrom(self.chrom), self.pos, self.ref.upper(), self.alt.upper())


def norm_chrom(chrom: str) -> str:
    """Normalise un nom de chromosome (chr13 -> 13, insensible a la casse)."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def _parse_info(info_field: str) -> dict:
    """Decode le champ INFO en dictionnaire. Convertit les nombres."""
    out: dict = {}
    if info_field in (".", ""):
        return out
    for item in info_field.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = _coerce(val)
        else:
            out[item] = True  # drapeau (ex: DUP)
    return out


def _coerce(val: str):
    """Convertit une chaine en int ou float si possible."""
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        return val


def _open_text(path: str):
    """Ouvre un fichier texte, gere .gz de maniere transparente."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def sample_pair_from_filename(path: str) -> tuple[str, str]:
    """Deduit (echantillon, paire) du nom de fichier."""
    base = os.path.basename(path)
    m = _FNAME_RE.search(base)
    if m:
        return m.group("sample"), m.group("pair")
    return base, ""


@dataclass
class VcfFile:
    """Un fichier VCF lu : ses metadonnees et ses variants."""

    path: str
    sample: str
    pair: str
    variants: list[Variant] = field(default_factory=list)
    meta: dict = field(default_factory=dict)


def read_vcf(path: str) -> VcfFile:
    """Lit un fichier VCF de resultats et renvoie un objet VcfFile."""
    sample_fn, pair = sample_pair_from_filename(path)
    sample_meta = None
    meta: dict = {}
    variants: list[Variant] = []

    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                # On tente d'extraire le nom d'echantillon de la ligne ##sample=
                if line.startswith("##sample="):
                    meta["sample_line"] = line[len("##sample=") :]
                    m = _SAMPLE_FASTQ_RE.search(line)
                    if m:
                        sample_meta = m.group("sample")
                elif line.startswith("##reference="):
                    meta["reference"] = line[len("##reference=") :]
                continue
            if line.startswith("#"):
                continue  # ligne d'en-tete des colonnes (#CHROM ...)

            cols = line.split("\t")
            if len(cols) < 8:
                # Tolerance : certains fichiers utilisent des espaces.
                cols = line.split()
            if len(cols) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            variants.append(
                Variant(
                    chrom=chrom,
                    pos=pos_i,
                    vid=vid,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    filt=filt,
                    info=_parse_info(info),
                )
            )

    # Le nom d'echantillon issu de la ligne ##sample= fait autorite ;
    # sinon on retombe sur le nom de fichier.
    sample = sample_meta or sample_fn
    for v in variants:
        v.sample = sample
        v.pair = pair
        v.source_file = os.path.basename(path)

    return VcfFile(path=path, sample=sample, pair=pair, variants=variants, meta=meta)
