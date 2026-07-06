"""Lecture et ecriture de fichiers OpenDocument Spreadsheet (.ods) en Python pur.

Un fichier .ods est une archive ZIP contenant un `content.xml`. On n'a donc
besoin d'aucune dependance externe (ni pandas, ni odfpy) pour lire les
fichiers de reference du laboratoire, ce qui est pratique sur un serveur ou
l'on ne peut pas forcement installer de paquets.

Fonctions publiques :
    read_table(path)  -> liste de lignes, chaque ligne etant une liste de
                         chaines (cellules). La premiere table du classeur est
                         renvoyee.
    read_rows(path)   -> alias de read_table.
    write_table(path, rows) -> ecrit une liste de lignes dans un .ods valide.
"""

from __future__ import annotations

import zipfile
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

# Espaces de noms OpenDocument.
NS = {
    "office": "urn:oasis:names:tc:opendocument:xmlns:office:1.0",
    "table": "urn:oasis:names:tc:opendocument:xmlns:table:1.0",
    "text": "urn:oasis:names:tc:opendocument:xmlns:text:1.0",
}

# Au dela de cette valeur, un attribut "repeated" est considere comme du
# remplissage de fin de ligne/table et est ignore (LibreOffice met souvent
# 1000+ pour les cellules vides finales).
_MAX_REPEAT = 4096


def _qn(prefix_tag: str) -> str:
    prefix, tag = prefix_tag.split(":")
    return "{%s}%s" % (NS[prefix], tag)


def _cell_text(cell: ET.Element) -> str:
    """Concatene le texte de tous les <text:p> d'une cellule."""
    parts = []
    for p in cell.iter(_qn("text:p")):
        parts.append("".join(p.itertext()))
    return "\n".join(parts)


def read_table(path: str) -> list[list[str]]:
    """Renvoie la premiere table du .ods sous forme de liste de lignes."""
    with zipfile.ZipFile(path) as zf:
        with zf.open("content.xml") as fh:
            tree = ET.parse(fh)

    root = tree.getroot()
    table = root.find(".//" + _qn("table:table"))
    if table is None:
        return []

    rows: list[list[str]] = []
    for row in table.findall(_qn("table:table-row")):
        row_repeat = int(row.get(_qn("table:number-rows-repeated"), "1"))
        cells: list[str] = []
        for cell in row.findall(_qn("table:table-cell")):
            repeat = int(cell.get(_qn("table:number-columns-repeated"), "1"))
            if repeat > _MAX_REPEAT:
                repeat = 1  # remplissage de fin de ligne
            value = _cell_text(cell)
            cells.extend([value] * repeat)

        # Supprime les cellules vides en fin de ligne.
        while cells and cells[-1] == "":
            cells.pop()

        if row_repeat > _MAX_REPEAT:
            row_repeat = 1
        for _ in range(row_repeat):
            rows.append(list(cells))

    # Supprime les lignes entierement vides en fin de table.
    while rows and not any(c.strip() for c in rows[-1]):
        rows.pop()

    return rows


# Alias plus explicite.
read_rows = read_table


_CONTENT_HEAD = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<office:document-content '
    'xmlns:office="urn:oasis:names:tc:opendocument:xmlns:office:1.0" '
    'xmlns:table="urn:oasis:names:tc:opendocument:xmlns:table:1.0" '
    'xmlns:text="urn:oasis:names:tc:opendocument:xmlns:text:1.0" '
    'office:version="1.3">'
    "<office:body><office:spreadsheet>"
)
_CONTENT_TAIL = "</office:spreadsheet></office:body></office:document-content>"

_MANIFEST = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<manifest:manifest '
    'xmlns:manifest="urn:oasis:names:tc:opendocument:xmlns:manifest:1.0" '
    'manifest:version="1.3">'
    '<manifest:file-entry manifest:full-path="/" '
    'manifest:media-type="application/vnd.oasis.opendocument.spreadsheet"/>'
    '<manifest:file-entry manifest:full-path="content.xml" '
    'manifest:media-type="text/xml"/>'
    "</manifest:manifest>"
)

_MIMETYPE = "application/vnd.oasis.opendocument.spreadsheet"


def write_table(path: str, rows: list[list[str]], table_name: str = "Sheet1") -> None:
    """Ecrit `rows` (liste de listes de chaines) dans un .ods valide."""
    body = [_CONTENT_HEAD, '<table:table table:name="%s">' % escape(table_name)]
    for row in rows:
        body.append("<table:table-row>")
        for cell in row:
            body.append(
                '<table:table-cell office:value-type="string">'
                "<text:p>%s</text:p></table:table-cell>" % escape(str(cell))
            )
        body.append("</table:table-row>")
    body.append("</table:table>")
    body.append(_CONTENT_TAIL)
    content = "".join(body)

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        # Le mimetype doit etre le premier membre et non compresse.
        zf.writestr(
            zipfile.ZipInfo("mimetype"), _MIMETYPE, compress_type=zipfile.ZIP_STORED
        )
        zf.writestr("META-INF/manifest.xml", _MANIFEST)
        zf.writestr("content.xml", content)
