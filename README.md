# Mutation — correspondance des variants VCF FLT3 avec les mutations connues

Les fichiers VCF proviennent de **FiLT3r** ([Baudry *et al*,
2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04983-6)),
un detecteur de **duplications internes (ITD) de FLT3**. Il ne detecte donc
**que les FLT3-ITD** : ni les autres genes, ni les mutations ponctuelles (SNV)
de FLT3. L'outil ci-dessous en tient compte (colonnes `cible_filt3r` /
`commentaire`).

Outil en **Python 3 pur (aucune dependance externe)** pour :

1. lire les fichiers VCF de resultats `JB_*_*_fast.gz.results.vcf` (FiLT3r,
   duplications / ITD de FLT3) ;
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
| `synthese_par_patient.tsv` | **Vue clinique par patient** : chaque mutation connue (HGVS + caryotype de `Results_patientJB`), croisee avec sa detection dans le VCF. Colonnes cle : `type_mutation`, `cible_filt3r` (est-ce une ITD FLT3, seule chose que FiLT3r detecte), `detecte`, `commentaire` (interpretation). |
| `correspondance_mutations_connues.tsv` | Une ligne par mutation connue **a coordonnees** (`Fichier_CHU`), avec `detecte = oui/non` et le niveau de correspondance. |
| `variants_par_echantillon.tsv` | Tous les variants retenus, un par echantillon, avec les metriques R1/R2 et le statut « connu / nouveau ». |
| `variants_nouveaux.tsv` | Les variants detectes **absents** des references, tries par VAF decroissante et **pre-classes** (colonne `triage`) pour distinguer un vrai ITD du bruit de fond. Voir plus bas. |

> **Portee de FiLT3r** : ces VCF ne contiennent que des variants du chromosome 13
> (region FLT3) et **uniquement des duplications (ITD)**. Dans
> `synthese_par_patient.tsv`, seules les lignes `cible_filt3r = oui` (= une ITD
> FLT3) peuvent legitimement etre `detecte = oui`. Tout le reste est marque hors
> perimetre dans `commentaire` :
>
> - **autre gene** (CEBPA, JAK2, IDH1...) -> hors perimetre FiLT3r ;
> - **SNV de FLT3** (ex. D835, `c.2503G>T`) -> FiLT3r ne detecte pas les
>   mutations ponctuelles.
>
> Autrement dit, `detecte = non` sur une ligne `cible_filt3r = oui` est un vrai
> signal (ITD attendue mais non retrouvee, a investiguer) ; ailleurs, c'est
> simplement hors perimetre.

Ajouter `--ods-output` pour ecrire aussi ces tables au format `.ods`.

## Carte des duplications (figure SVG)

`plot_itd.py` produit une figure a partir de `variants_par_echantillon.tsv` :

```bash
python3 plot_itd.py --input resultats/variants_par_echantillon.tsv \
    --out resultats/carte_itd.svg
```

- **axe Y (a gauche)** : positions du chromosome 13, de la plus basse (en haut)
  a la plus haute (en bas) ;
- **axe X (en haut)** : les patients ;
- chaque duplication (`DUP=oui` par defaut ; `--all` pour tous les variants) est
  un segment `[pos, pos+svlen]` avec un point a son debut, taille ~ VAF ;
- couleurs : **rouge** = ITD de reference (`connu=oui`, etendue surlignee),
  **bleu** = nouveau variant (`connu=non`).

Un variant bleu proche (voire inclus dans l'etendue) d'une ITD rouge est un
**sous-groupe probable** de la duplication principale. Le SVG est vectoriel
(zoom sans perte, ouvrable dans un navigateur ou LibreOffice Draw).

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
| `--drop-recurrent N` | supprimer les variants presents dans **>= N echantillons** (artefacts systematiques) | `0` (off) |
| `--no-require-pass` | ne pas exiger `FILTER=PASS` | (PASS exige) |

Par defaut aucun seuil de VAF/lectures n'est impose (seul `FILTER=PASS` l'est) :
les ITD FLT3 pertinentes peuvent avoir une VAF tres faible. Ajustez les seuils
selon votre pratique.

### Reduire le bruit de fond

Un seul seuil de VAF ne suffit pas (un artefact recurrent peut passer au-dessus,
un vrai ITD minoritaire en dessous). **Combinez** plutot les filtres. Effet
mesure sur un jeu reel de ~7000 variants non connus :

| Filtres | Variants restants |
|---|---|
| aucun | 7062 |
| `--min-m 10` | 2956 |
| `--min-m 10 --require-both-pairs` | 1590 |
| `--min-m 10 --require-both-pairs --only-dup` | 117 |
| `--min-m 10 --require-both-pairs --only-dup --drop-recurrent 5` | 0 |

`--drop-recurrent` est la cle contre les artefacts **systematiques** (memes
positions chez de nombreux patients) que les seuils de VAF/M ne retirent pas.
Un ITD FLT3 etant propre a un patient, un seuil eleve (>= 5) ne supprime que du
bruit. Pour de la recherche d'ITD **minoritaire** (MRD), gardez un `--min-m`
plus bas mais conservez `--require-both-pairs` et `--drop-recurrent`.

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

## Distinguer un vrai ITD d'un bruit de fond

FiLT3r ressort **beaucoup** de variants a tres basse frequence (bruit
d'alignement/sequencage). `variants_nouveaux.tsv` est donc trie par VAF
decroissante et annote pour faciliter le tri :

| Colonne | Aide au tri |
|---|---|
| `VAF_max`, `M_total` | un vrai ITD a une VAF et un nombre de reads **nettement au-dessus** du bruit du meme echantillon. |
| `trouve_dans_paires` | `1,2` (present dans R1 **et** R2) = plus fiable. |
| `DUP` | un vrai ITD est une **duplication** (`oui`). |
| `cadre_lecture` | un ITD FLT3 est generalement **in-frame** (taille multiple de 3). |
| `nb_echantillons_meme_variant` | **critere cle** : un variant present chez **beaucoup** de patients a basse VAF est un **artefact systematique**, pas un ITD propre a un patient. |
| `triage` | pre-classement synthetique (voir ci-dessous). |

La colonne `triage` combine ces criteres :

- **`artefact probable`** — variant recurrent (>= `--artefact-nb-echantillons`,
  defaut 5 echantillons) : artefact systematique ;
- **`bruit probable`** — VAF < `--candidat-vaf-min` (defaut 0.01 = 1 %) ;
- **`candidat ITD`** — duplication, presente dans R1+R2, in-frame et VAF >= seuil ;
- **`a verifier`** — VAF suffisante mais un critere manque (une seule paire,
  non-duplication, hors-cadre).

> C'est une **aide au tri, pas un diagnostic** : ajustez les seuils avec
> `--candidat-vaf-min` et `--artefact-nb-echantillons`. Un vrai ITD deja connu
> apparait dans `synthese_par_patient.tsv` (il n'est pas "nouveau") ; cette table
> sert surtout a reperer une **ITD non encore connue**.

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

Ces choix sont volontairement explicites :

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
