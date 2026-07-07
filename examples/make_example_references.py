#!/usr/bin/env python3
"""Genere les fichiers .ods de reference d'exemple.

Reproductible : relancer ce script recree
  * examples/Fichier_CHU.ods       (format "coordonnees", pour le matching VCF)
  * examples/Results_patientJB.ods (format "clinique large", par patient)
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mutmatch import ods  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

ALT_ITD = ("CAAATTAGCAGGGTGTTGTGACAGGTGCCACCCAGCCTGGCCACCGTGGTG"
           "AAACTCCGTCTCTACTAAAAATACAA")

# Fichier_CHU : une mutation par ligne, avec coordonnees genomiques.
# La colonne Query est en nomenclature 'GENE:c.xxx' (comme dans les vraies
# donnees), ce qui permet de la relier au fichier clinique.
fichier_chu = [
    ["Query", "#chr", "start", "ID", "ref", "alt", "patient id", "sample_id"],
    # ITD FLT3 presente dans le VCF -> correspondance exacte
    ["FLT3:c.1747_1794dup", "13", "28034317", ".", "C", ALT_ITD, "JB_01", "JB_01"],
    # variant FLT3 dont seule la position coincide -> correspondance 'position'
    ["FLT3:c.pos", "13", "28034314", ".", "G", "GA", "JB_01", "JB_01"],
    # mutation CEBPA (chr19) : hors de la portee du VCF FLT3 -> non couverte
    ["CEBPA:c.1015_1027dup", "19", "33301387", ".", "C", "CGGAAGATGCCCCG",
     "JB_01", "JB_01"],
    # SNV FLT3 absent du VCF -> non detecte
    ["FLT3:c.2503G>T", "13", "28035000", ".", "A", "T", "JB_01", "JB_01"],
]

# Results_patientJB : format clinique large, une ligne par patient. Les
# colonnes de mutations (a partir de 'Mutations NGS') sont en nombre variable.
results_jb = [
    ["N° Echantillon", "Sample", "Caryotype", "Mutations NGS"],
    ["SMO000001 X Y", "JB_01", "46,XX[40]",
     "FLT3 : c.1747_1794dup; p.G583_E598dup",
     "CEBPA : c.1015_1027dup; p.Arg343ProfsTer64",
     "FLT3 : c.2503G>T; p.Asp835Tyr"],
]

ods.write_table(os.path.join(HERE, "Fichier_CHU.ods"), fichier_chu)
ods.write_table(os.path.join(HERE, "Results_patientJB.ods"), results_jb)
print("Fichiers .ods d'exemple generes.")
