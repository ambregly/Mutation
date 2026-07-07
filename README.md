# Mutation — correspondance des variants VCF FLT3 avec les mutations connues

Outil en **Python 3 pur (aucune dependance externe)** pour :

1. lire les fichiers VCF de resultats `JB_*_*_fast.gz.results.vcf` (detection de
   duplications / ITD sur FLT3) ;
2. appliquer des **filtres** (VAF, nombre de lectures, PASS, DUP...) ;
3. **fusionner les deux paires** de lecture (R1 = `_1`, R2 = `_2`) par
   echantillon ;
4. **comparer** les variants detectes aux mutations deja connues des patients
   (`Fichier_CHU.ods`, `Results_patientJB.ods`) ;
5. produire des rapports **TSV** (et `.ods` en option) ouvrables dans
   LibreOffice / Excel.

Aucune installation de paquet n'est requise : seule la bibliotheque standard de
Python 3.8+ est utilisee (pratique sur un serveur comme `ella` ou l'on ne peut
pas toujours faire `pip install`).

## Utilisation rapide

```bash
python3 match_mutations.py \
    --vcf-dir /data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/fastq \
    --reference Fichier_CHU.ods \
    --reference Results_patientJB.ods \
    --out resultats \
    --min-vaf 0.001 --min-m 3
```

Cela produit dans `resultats/` :

| Fichier | Contenu |
|---|---|
| `synthese_par_patient.tsv` | **Vue clinique par patient** : chaque mutation connue (HGVS + caryotype de `Results_patientJB`), croisee avec sa detection dans le VCF. Colonne `couvert_par_ce_vcf` = le gene est-il dans la portee du VCF. |
| `correspondance_mutations_connues.tsv` | Une ligne par mutation connue **a coordonnees** (`Fichier_CHU`), avec `detecte = oui/non` et le niveau de correspondance. |
| `variants_par_echantillon.tsv` | Tous les variants retenus, un par echantillon, avec les metriques R1/R2 et le statut « connu / nouveau ». |
| `variants_nouveaux.tsv` | Les variants detectes **absents** des fichiers de reference (candidats a verifier). |

> **Portee des VCF** : ces fichiers proviennent d'un detecteur d'ITD/duplications
> **FLT3** (reference `sequence-FLT3-Ex13-14-15-20.fa`) et ne contiennent que des
> variants du chromosome 13 (region FLT3). Les mutations connues d'autres genes
> (CEBPA, JAK2, IDH1...) ne peuvent donc pas y etre retrouvees : la colonne
> `couvert_par_ce_vcf` vaut alors `non`. Pour les valider, il faut les VCF du
> panel correspondant.

Ajouter `--ods-output` pour ecrire aussi ces tables au format `.ods`.

## Essayer sur les donnees d'exemple

Le depot contient un jeu d'exemple (`examples/`) reproduisant vos formats :

```bash
python3 examples/make_example_references.py   # (re)genere les .ods d'exemple
python3 match_mutations.py \
    --vcf-dir examples/vcf \
    --reference examples/Fichier_CHU.ods \
    --reference examples/Results_patientJB.ods \
    --out /tmp/resultats
```

## Filtres disponibles

| Option | Effet | Defaut |
|---|---|---|
| `--min-vaf X` | VAF minimale (ratio) | `0.0` |
| `--min-m N` | nombre minimal de lectures mutees (`M`) | `0` |
| `--min-wt N` | nombre minimal de lectures sauvages (`WT`) | `0` |
| `--only-dup` | ne garder que les evenements `DUP` (duplications / ITD) | off |
| `--require-both-pairs` | variant present dans **R1 ET R2** (plus fiable) | off |
| `--no-require-pass` | ne pas exiger `FILTER=PASS` | (PASS exige) |

Par defaut aucun seuil de VAF/lectures n'est impose (seul `FILTER=PASS` l'est) :
les ITD FLT3 pertinentes peuvent avoir une VAF tres faible. Ajustez les seuils
selon votre pratique.

## Comment se fait la correspondance

Chaque variant detecte est identifie par la cle `(chromosome, position, ref,
alt)`. Pour chaque mutation connue :

- **`exact`** — le referentiel fournit `chr`, `start`, `ref` et `alt` et les
  quatre coincident avec un variant detecte ;
- **`position`** — seule la position `(chr, pos)` coincide (ref/alt differents
  ou absents) ;
- **`none`** — mutation connue non retrouvee dans les VCF.

**Par defaut**, une mutation connue n'est cherchee que dans l'echantillon dont
le `sample_id` / `patient id` correspond : cela evite les faux positifs
inter-patients (une position qui coincide par hasard chez un autre patient).
Utilisez `--no-restrict-to-sample` pour chercher dans tous les echantillons.

### Detection des colonnes de reference

La lecture des `.ods` / `.csv` / `.tsv` reconnait les colonnes de maniere
souple (insensible a la casse, tolere `#`, accents, espaces et quelques
synonymes) :

Le **type** de chaque fichier `--reference` est reconnu automatiquement :

- **Format coordonnees** (ex. `Fichier_CHU`) : une mutation par ligne, colonnes
  `Query`, `#chr`, `start`, `ID`, `ref`, `alt`, `patient id`, `sample_id`. Une
  colonne `start`/`Position` peut aussi encoder la coordonnee complete
  `chr-pos-ref-alt` (ex. `19-33301387-C-CGGAAGATGCCCCG`). C'est ce format qui
  porte les coordonnees et pilote le matching avec les VCF.

- **Format clinique large** (ex. `Results_patientJB`) : **une ligne par
  patient**, colonnes `N° Echantillon`, `Sample`, `Caryotype`, puis une ou
  plusieurs colonnes de mutations en nomenclature HGVS
  (`GENE : c.xxx; p.yyy`). Ce format n'a pas de coordonnees ; il est relie au
  fichier coordonnees par `(echantillon, GENE:c.xxx)` pour produire la
  `synthese_par_patient.tsv`.

> **Multi-genes** : le pipeline produit un VCF par gene/amplicon (FLT3 sur le
> chr13, CEBPA sur le chr19, etc.). Pour retrouver une mutation d'un gene donne,
> le dossier `--vcf-dir` doit contenir les VCF correspondants ; l'outil regroupe
> automatiquement tous les VCF d'un meme echantillon (`JB_01`), quel que soit le
> gene.

## Hypotheses / points a valider

Ces choix sont volontairement explicites — dites-moi si l'un ne correspond pas a
votre besoin, ils sont faciles a ajuster :

1. **Coordonnees** : on suppose que `start` du `Fichier_CHU` est la position VCF
   (1-based). Si vos coordonnees sont 0-based ou proviennent d'un autre
   referentiel (hg19/hg38), la correspondance exacte peut echouer et retomber en
   `position`. A verifier sur un cas connu.
2. **Fusion R1/R2** : par defaut un variant present dans **au moins une** paire
   est conserve (union). `--require-both-pairs` bascule en intersection.
3. **Nom d'echantillon** : deduit de la ligne `##sample=` du VCF (ex. `JB_01`),
   avec repli sur le nom de fichier.

## Structure du depot

```
match_mutations.py        # interface en ligne de commande
mutmatch/
  vcf.py                  # lecture des VCF (gere .gz)
  ods.py                  # lecture/ecriture .ods en Python pur
  references.py           # lecture des fichiers de reference (ods/csv/tsv)
  match.py                # filtres, fusion des paires, correspondance
  report.py               # ecriture des rapports
examples/                 # jeu de donnees d'exemple
tests/                    # tests (python3 -m unittest discover -s tests)
```

## Tests

```bash
python3 -m unittest discover -s tests -v
```
