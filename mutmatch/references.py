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


def _parse_position(raw: str) -> tuple[str, int | None]:
    """Decode une valeur 'Position' -> (chrom, pos).

    Gere 'chr13:28034317', '13:28034317', '13-28034317' ou un simple entier.
    """
    raw = str(raw).strip()
    if not raw:
        return "", None
    m = re.match(r"^\s*(?:chr)?([0-9XYMxym]+)\s*[:\-_ ]\s*([0-9,]+)", raw)
    if m:
        pos = int(m.group(2).replace(",", ""))
        return m.group(1), pos
    m = re.match(r"^\s*([0-9,]+)\s*$", raw)
    if m:
        return "", int(m.group(1).replace(",", ""))
    return "", None


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
        pos_raw = _cell(raw_row, cols.get("pos"))
        pos: int | None = None
        position_raw = pos_raw

        if pos_raw:
            # 'start' peut deja etre une position ou une valeur 'chrom:pos'.
            c2, p2 = _parse_position(pos_raw)
            if p2 is not None:
                pos = p2
                if not chrom and c2:
                    chrom = c2
            else:
                try:
                    pos = int(str(pos_raw).replace(",", ""))
                except ValueError:
                    pos = None

        # Fichier de type Results_patientJB : colonne 'Position' unique.
        if pos is None and "pos" not in cols:
            # cherche une colonne dont le nom ressemble a 'position'
            for i, h in enumerate(header):
                if _slug(h) in {"position", "pos"}:
                    c2, p2 = _parse_position(_cell(raw_row, i))
                    chrom = chrom or c2
                    pos = p2
                    position_raw = _cell(raw_row, i)
                    break

        km = KnownMutation(
            query=_cell(raw_row, cols.get("query")),
            chrom=chrom,
            pos=pos,
            ref=_cell(raw_row, cols.get("ref")),
            alt=_cell(raw_row, cols.get("alt")),
            patient_id=_cell(raw_row, cols.get("patient_id")),
            sample_id=_cell(raw_row, cols.get("sample_id")),
            position_raw=position_raw,
            source=source,
            row={header[i] if i < len(header) else "col%d" % i: _cell(raw_row, i)
                 for i in range(len(raw_row))},
        )
        out.append(km)

    return out
