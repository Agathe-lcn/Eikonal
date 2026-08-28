# Eikonal2D_MPI

Le dossier Eikonal2D_MPI est un projet en C permettant de calculer sur une grille 2D la distance entre chaque point de la grille et un ensemble de points sources à l'aide de l'équation eikonale et de la Fast Iterative Method (FIM) implémentée en parallèle. Il est également possible de retrouver des fichiers python permettant de visualiser les résultats obtenus.

## Auteur

Agathe Luciani


### Visualisation (Python)

Installation des dépendances Python :

```bash
pip install numpy matplotlib pandas
```

## Compilation

Compiler tous les exécutables:

```bash
make all
```

### Exécution

Compiler et exécuter manuellement un programme spécifique:

```bash
make fim-mpi
make benchmark-overlap
make benchmark-optim
```

### Nettoyage

```bash
make clean
```

## Utilisation

### Génération d'un maillage et application de la FIM : `fim`

```bash
./bin/fim_mpi_run NP=nb_proc mon_fichier.txt
```

ou

```bash
make fim-mpi NP=nb_proc ARGS=mon_fichier.txt
```

Si aucun fichier n'est fourni en argument, `fim-mpi` cherche automatiquement un fichier nommé `config.txt` dans le dossier courant.

Le fichier de configuration `.txt` doit respecter la structure suivante :

```
n = ...
m = ...
h = ...
max_depth = ...
overlap = ...
OPTIM_COM_MPI = ...


sources:
... ... 
... ...

walls:
... ... ... ...
```

- `n`, `m` : dimensions de la grille
- `h` : pas de la grille
- `max_depth` : profondeur maximale
- `overlap` : taille de la zone de recouvrement
- `OPTIM_COM_MPI` : booléen activant ou non une optimisation sur les communications
- `sources` : liste des coordonnées des points sources
- `walls` : liste des coordonnées des obstacles (colonne de départ, colonne de fin, ligne de départ, ligne de fin)

Tous les fichiers générés (grilles, résultats) sont au format `.txt`.

### Profondeur maximale (`max_depth`)

Le paramètre `max_depth` de la fonction `fim_solve` (appelée dans `fim_run.c`) contrôle la zone parcourue par la FIM autour des sources:

- **`max_depth` égale à -1** : la FIM parcourt toute la grille.
- **`max_depth` strictement positif** : la valeur définit un seuil, et seuls les points situés à une distance inférieure à ce seuil des sources sont traités.

## Visualisation

Les scripts du dossier `visualization/` permettent d'analyser les résultats produits par la FIM. Ils génèrent des images au format `.png`.

- **`benchmark_overlap.py`**: visualisation du temps de calcul et du nombre d'échanges en fonction de overlap.

  ```bash
  python3 visualization/benchmark_overlap.py
  ```

- **`erreur_euclid.py`**: visualisation de l'erreur absolue entre les résultats de la FIM parallèle et la distance euclidienne. Utilisable uniquement dans les cas avec une ou plusieurs sources et **aucun obstacle**.

  ```bash
  python3 visualization/erreur_euclid.py
  ```

- **`fim_compare.py`**: visualisation de l'erreur absolue entre les résultats de la FIM séquentielle et de la FIM parallèle.

  ```bash
  python3 visualization/fim_compare.py resultats_sequentiel.txt resultats_mpi.txt
  ```

- **`visu_fim_global.py`** : visualisation du maillage totale obtenue après l'exécution de la FIM parallèle.

  ```bash
  python3 visualization/visu_fim_global.py
  ```

- **`visualisation_fim_local.py`** : visualisation des sous-domaines de chacun des processus après l'exécution de la FIM parallèle.

  ```bash
  python3 visualization/visu_fim_local.py
  ```