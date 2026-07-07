"""Lecture des fichiers de reference decrivant les mutations connues.

Deux formats sont geres (fichiers .ods, .csv ou .tsv) :

  * Fichier "coordonnees" (ex. Fichier_CHU) : une mutation par ligne, avec des
    colonnes Query, #chr, start, ID, ref, alt, patient id, sample_id. Ce format
    porte les coordonnees genomiques et sert au matching avec les VCF.

  * Fichier "clinique large" (ex. Results_patientJB) : une ligne par patient,
    avec N Echantillon, Sample, Caryotype, puis une ou plusieurs colonnes de
    mutations en nomenclature HGVS ('GENE : c.xxx; p.xxx'), sans coordonnees.

La detection des colonnes est souple (insensible a la casse, tolere les
espaces, '#', accents et quelques synonymes). Le type de fichier est reconnu
automatiquement (voir read_any / is_clinical_wide).
"""

from __future__ import annotations

import csv
import os
import re
import unicodedata
from dataclasses import dataclass

from . import ods
from .vcf import norm_chrom


@dataclass
class KnownMutation:
    """Une mutation connue issue d'un fichier de reference."""

    query: str = ""
    chrom: str = ""
    pos: int | None = None
    ref: str = ""
    alt: str = ""
    patient_id: str = ""
    sample_id: str = ""
    position_raw: str = ""
    source: str = ""
    row: dict = None

    def key(self) -> tuple | None:
        """Cle (chrom, pos, ref, alt) si assez d'info, sinon None."""
        if self.pos is None:
            return None
        return (norm_chrom(self.chrom), self.pos, self.ref.upper(), self.alt.upper())

    def pos_key(self) -> tuple | None:
        """Cle de position seule (chrom, pos)."""
        if self.pos is None:
            return None
        return (norm_chrom(self.chrom), self.pos)


def _slug(name: str) -> str:
    """Normalise un nom de colonne pour la comparaison."""
    s = unicodedata.normalize("NFKD", str(name))
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower().strip()
    s = s.replace("#", "").replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


# Synonymes acceptes pour chaque champ logique (apres _slug).
_COLUMN_SYNONYMS = {
    "query": {"query", "requete"},
    "chrom": {"chr", "chrom", "chromosome"},
    "pos": {"start", "pos", "position", "debut"},
    "vid": {"id", "identifiant"},
    "ref": {"ref", "reference", "allele ref"},
    "alt": {"alt", "alternate", "allele alt", "mutation"},
    "patient_id": {"patient id", "patient", "idpatient", "id patient"},
    "sample_id": {"sample id", "sample", "echantillon", "id sample"},
}


def _map_columns(header: list[str]) -> dict[str, int]:
    """Associe chaque champ logique a l'indice de colonne correspondant."""
    slugs = [_slug(h) for h in header]
    mapping: dict[str, int] = {}
    for field, names in _COLUMN_SYNONYMS.items():
        for i, s in enumerate(slugs):
            if s in names and field not in mapping:
                mapping[field] = i
                break
    return mapping


def _read_rows(path: str) -> list[list[str]]:
    """Lit un fichier .ods / .csv / .tsv en liste de lignes."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".ods":
        return ods.read_table(path)
    if ext in (".tsv", ".txt"):
        with open(path, "r", encoding="utf-8-sig", newline="") as fh:
            return [row for row in csv.reader(fh, delimiter="\t")]
    if ext == ".csv":
        with open(path, "r", encoding="utf-8-sig", newline="") as fh:
            sample = fh.read(4096)
            fh.seek(0)
            delim = ";" if sample.count(";") > sample.count(",") else ","
            return [row for row in csv.reader(fh, delimiter=delim)]
    raise ValueError("Extension non prise en charge : %s" % ext)


_BASE_RE = re.compile(r"^[ACGTNacgtn]+$")
_CHROM_RE = re.compile(r"^[0-9]+$|^[XYMxym]+$")
_INT_RE = re.compile(r"^[0-9][0-9,]*$")


def _parse_coord(raw: str) -> tuple[str, int | None, str, str]:
    """Decode une valeur de coordonnee -> (chrom, pos, ref, alt).

    Formats reconnus :
      * 'chr-pos-ref-alt'  ex. '19-33301387-C-CGGAAGATGCCCCG'
      * 'chr:pos:ref:alt', 'chr pos ref alt', melanges de -, :, _ , espace
      * 'chr:pos' / 'chr-pos'  (ref/alt vides)
      * 'chr13:...' (prefixe 'chr' ignore)
      * une position entiere seule
    """
    raw = str(raw).strip()
    if not raw:
        return "", None, "", ""

    parts = [p for p in re.split(r"[\s:_\-]+", raw) if p != ""]
    if parts and parts[0].lower().startswith("chr"):
        parts[0] = parts[0][3:]

    # chr-pos-ref-alt : chrom, pos entier, puis deux chaines de bases
    if (len(parts) >= 4 and _CHROM_RE.match(parts[0]) and _INT_RE.match(parts[1])
            and _BASE_RE.match(parts[2]) and _BASE_RE.match(parts[3])):
        return (parts[0], int(parts[1].replace(",", "")),
                parts[2].upper(), parts[3].upper())

    # chr:pos
    if len(parts) >= 2 and _CHROM_RE.match(parts[0]) and _INT_RE.match(parts[1]):
        return parts[0], int(parts[1].replace(",", "")), "", ""

    # position seule
    if len(parts) == 1 and _INT_RE.match(parts[0]):
        return "", int(parts[0].replace(",", "")), "", ""

    return "", None, "", ""


def _find_coord_in_row(row: list[str]) -> tuple[str, int | None, str, str]:
    """Repli : cherche dans toute la ligne une cellule ressemblant a une coordonnee."""
    for cell in row:
        c, p, r, a = _parse_coord(cell)
        if p is not None:
            return c, p, r, a
    return "", None, "", ""


def _cell(row: list[str], idx: int | None) -> str:
    if idx is None or idx >= len(row):
        return ""
    return str(row[idx]).strip()


def read_reference(path: str) -> list[KnownMutation]:
    """Lit un fichier de reference et renvoie la liste des mutations connues."""
    rows = _read_rows(path)
    if not rows:
        return []

    header = rows[0]
    cols = _map_columns(header)
    source = os.path.basename(path)
    out: list[KnownMutation] = []

    for raw_row in rows[1:]:
        if not any(str(c).strip() for c in raw_row):
            continue  # ligne vide

        chrom = _cell(raw_row, cols.get("chrom"))
        ref = _cell(raw_row, cols.get("ref"))
        alt = _cell(raw_row, cols.get("alt"))
        pos_raw = _cell(raw_row, cols.get("pos"))
        pos: int | None = None
        position_raw = pos_raw

        # La colonne 'start'/'position' peut etre une simple position ou une
        # coordonnee complete 'chr-pos-ref-alt'.
        c2, p2, r2, a2 = _parse_coord(pos_raw)
        if p2 is not None:
            pos = p2
            chrom = chrom or c2
            ref = ref or r2
            alt = alt or a2

        # Repli : si aucune position trouvee dans la colonne attendue (cas ou
        # l'en-tete ne correspond pas au nombre de colonnes), on scanne la ligne.
        if pos is None:
            c2, p2, r2, a2 = _find_coord_in_row(raw_row)
            if p2 is not None:
                pos = p2
                chrom = chrom or c2
                ref = ref or r2
                alt = alt or a2
                if not position_raw:
                    position_raw = "%s-%s-%s-%s" % (c2, p2, r2, a2)

        km = KnownMutation(
            query=_cell(raw_row, cols.get("query")),
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            patient_id=_cell(raw_row, cols.get("patient_id")),
            sample_id=_cell(raw_row, cols.get("sample_id")),
            position_raw=position_raw,
            source=source,
            row={header[i] if i < len(header) else "col%d" % i: _cell(raw_row, i)
                 for i in range(len(raw_row))},
        )
        out.append(km)

    return out


# --- Format clinique large (Results_patientJB) ---------------------------

@dataclass
class ClinicalMutation:
    """Une mutation connue en nomenclature HGVS, rattachee a un patient."""

    sample_id: str = ""
    n_echantillon: str = ""
    caryotype: str = ""
    gene: str = ""
    hgvs: str = ""      # 'GENE:c.xxx'
    hgvs_c: str = ""    # 'c.xxx'
    hgvs_p: str = ""    # 'p.xxx' (ou vide)
    raw: str = ""
    source: str = ""


# 'GENE : c.xxxx; p.yyyy'  (le ' ; p...' est optionnel)
_MUT_RE = re.compile(r"^\s*([A-Za-z0-9]+)\s*:\s*(c\.[^;]+?)\s*(?:;\s*(.*))?$")


def _parse_mutation_cell(cell: str):
    """Decode 'GENE : c.xxx; p.yyy' -> (gene, hgvs_c, hgvs_p) ou None."""
    cell = str(cell).strip()
    if not cell:
        return None
    m = _MUT_RE.match(cell)
    if not m:
        return None
    gene = m.group(1).strip()
    hgvs_c = m.group(2).strip()
    hgvs_p = (m.group(3) or "").strip()
    return gene, hgvs_c, hgvs_p


def is_clinical_wide(rows: list[list[str]]) -> bool:
    """Reconnait le format clinique large (une ligne par patient)."""
    if not rows:
        return False
    slugs = {_slug(h) for h in rows[0]}
    if any("mutations" in s for s in slugs):
        return True
    return "sample" in slugs and any("caryotype" in s for s in slugs)


def _first_col(slugs: list[str], predicate) -> int | None:
    for i, s in enumerate(slugs):
        if predicate(s):
            return i
    return None


def read_clinical_wide(path: str) -> list[ClinicalMutation]:
    """Lit un fichier clinique large -> liste de ClinicalMutation (a plat)."""
    rows = _read_rows(path)
    if not rows:
        return []
    header = rows[0]
    slugs = [_slug(h) for h in header]
    source = os.path.basename(path)

    i_sample = _first_col(slugs, lambda s: s == "sample")
    i_caryo = _first_col(slugs, lambda s: "caryotype" in s)
    i_num = _first_col(slugs, lambda s: "echantillon" in s or s.startswith("n "))
    i_mut = _first_col(slugs, lambda s: "mutations" in s)
    if i_mut is None:
        i_mut = (i_caryo + 1) if i_caryo is not None else 3

    out: list[ClinicalMutation] = []
    for row in rows[1:]:
        if not any(str(c).strip() for c in row):
            continue
        sample = _cell(row, i_sample)
        caryo = _cell(row, i_caryo)
        nech = _cell(row, i_num)
        # Les mutations occupent la colonne 'Mutations NGS' et toutes les
        # colonnes suivantes (format irregulier : autant de colonnes que de
        # mutations).
        for j in range(i_mut, len(row)):
            parsed = _parse_mutation_cell(_cell(row, j))
            if not parsed:
                continue
            gene, hgvs_c, hgvs_p = parsed
            out.append(ClinicalMutation(
                sample_id=sample,
                n_echantillon=nech,
                caryotype=caryo,
                gene=gene,
                hgvs="%s:%s" % (gene, hgvs_c),
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                raw=_cell(row, j),
                source=source,
            ))
    return out


def mutation_type(hgvs_c: str) -> str:
    """Classe une mutation d'apres sa nomenclature HGVS 'c.'.

    Renvoie : 'duplication', 'insertion', 'delins', 'deletion', 'SNV' ou 'autre'.
    Sert a savoir ce que FiLT3r peut detecter (uniquement les duplications).
    """
    s = str(hgvs_c).lower()
    if "delins" in s:
        return "delins"
    if "dup" in s:
        return "duplication"
    if "ins" in s:
        return "insertion"
    if "del" in s:
        return "deletion"
    if ">" in s:
        return "SNV"
    return "autre"


def norm_hgvs(text: str) -> str:
    """Normalise une chaine HGVS pour comparer Fichier_CHU et Results.

    Retire les espaces, met en majuscules et ne garde que la partie
    'GENE:c.xxx' (avant un eventuel '; p...').
    """
    s = str(text).split(";")[0]
    s = re.sub(r"\s+", "", s).upper()
    return s


def read_any(path: str):
    """Lit un fichier de reference et devine son type.

    Renvoie ('clinical', [ClinicalMutation]) ou ('coord', [KnownMutation]).
    """
    rows = _read_rows(path)
    if is_clinical_wide(rows):
        return "clinical", read_clinical_wide(path)
    return "coord", read_reference(path)
