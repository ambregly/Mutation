"""Lecture des fichiers de reference decrivant les mutations connues.

Deux fichiers sont attendus (formats .ods, .csv ou .tsv acceptes) :

  * Fichier_CHU        : colonnes Query, #chr, start, ID, ref, alt,
                         patient id, sample_id
  * Results_patientJB  : colonnes Query, Position

La detection des colonnes est souple (insensible a la casse, tolere les
espaces, '#', accents et quelques synonymes), de sorte que de petites
variations d'en-tete ne cassent pas la lecture.
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
