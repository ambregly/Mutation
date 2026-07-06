#!/usr/bin/env python3
"""Genere les fichiers .ods de reference d'exemple.

Reproductible : relancer ce script recree examples/Fichier_CHU.ods et
examples/Results_patientJB.ods.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mutmatch import ods  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

ALT_ITD = ("CAAATTAGCAGGGTGTTGTGACAGGTGCCACCCAGCCTGGCCACCGTGGTG"
           "AAACTCCGTCTCTACTAAAAATACAA")

# Fichier_CHU : mutations connues avec coordonnees completes.
fichier_chu = [
    ["Query", "#chr", "start", "ID", "ref", "alt", "patient id", "sample_id"],
    # correspondance exacte avec le variant detecte dans JB_01
    ["FLT3-ITD-76", "13", "28034317", ".", "C", ALT_ITD, "P001", "JB_01"],
    # mutation connue NON detectee (pour illustrer un "non detecte")
    ["FLT3-D835", "13", "28035000", ".", "A", "T", "P001", "JB_01"],
]

# Results_patientJB : Query / Position.
# La colonne Position encode la coordonnee complete 'chr-pos-ref-alt'
# (ex. reel : 19-33301387-C-CGGAAGATGCCCCG). On ajoute volontairement une
# colonne intermediaire, comme dans les fichiers reels, pour verifier le repli
# de detection de coordonnee (le lecteur scanne la ligne si besoin).
results_jb = [
    ["Query", "Position"],
    ["FLT3-ITD-76;p.Xaa", "FLT3-ITD-76", "13-28034317-C-" + ALT_ITD],
    # position seule -> correspondance par position
    ["FLT3-ins36", "13:28034314"],
]

ods.write_table(os.path.join(HERE, "Fichier_CHU.ods"), fichier_chu)
ods.write_table(os.path.join(HERE, "Results_patientJB.ods"), results_jb)
print("Fichiers .ods d'exemple generes.")
