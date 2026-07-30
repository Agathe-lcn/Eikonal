# Eikonal2D

Le dossier Eikonal2D est un projet en C permettant de calculer sur une grille 2D la distance entre chaque point de la grille et un ensemble de points sources à l'aide de l'équation eikonale et de la Fast Iterative Method (FIM). Il est également possible de retrouver des fichiers python permettant de visualiser les résultats obtenus (carte de distance, vitesse, source la plus proche, erreur par rapport à la distance euclidienne).

## Auteur

Agathe Luciani


### Visualisation (Python)

Installation des dépendances Python :

```bash
pip install numpy matplotlib pandas
```

## Compilation

Compiler tours les exécutables:

```bash
make all
```

### Exécution

Compiler et exécuter manuellement un programme spécifique:

```bash
make run-generate
make run-1source
make run-5sources
make benchmark
```

### Nettoyage

```bash
make clean
```

## Utilisation

### Génération d'une grille : `generate_grid`

```bash
./bin/generate_grid mon_fichier.txt
```

Si aucun fichier n'est fourni en argument, `generate_grid` cherche automatiquement un fichier nommé `config.txt` dans le dossier courant.

Le fichier de configuration `.txt` doit respecter la structure suivante :

```
n = ...
m = ...
h = ...

sources:
... ... 
... ...

walls:
... ... ... ...
```

- `n`, `m` : dimensions de la grille
- `h` : pas de la grille
- `sources` : liste des coordonnées des points sources
- `walls` : liste des coordonnées des obstacles (colonne de départ, colonne de fin, ligne de départ, ligne de fin)

Tous les fichiers générés (grilles, résultats) sont au format `.txt`.

### Vitesse variable

Par défaut, la vitesse de propagation `F` est constante égale à 1. Pour utiliser une vitesse variable dans `generate_grid.c` :

1. Implémenter la fonction de vitesse souhaitée dans `set_variable_speed` (fichier `grids/generate_grid.c`).
2. Dans le `main()` de `generate_grid.c` :
   - commenter la ligne :
     ```c
     eikonal_grid_set_speed_constant(g, 1.0);
     ```
   - décommenter la ligne :
     ```c
     set_variable_speed(g, cfg.n, cfg.m);
     ```

Pour revenir à une vitesse constante, il suffit d'inverser ces deux étapes: décommenter `eikonal_grid_set_speed_constant(g, 1.0)` et recommenter `set_variable_speed(g, cfg.n, cfg.m)`.

### Rayon maximal (`max_radius`)

Le paramètre `max_radius` de la fonction `fim_solve` (appelée dans `generate_grid.c`) contrôle la zone parcourue par la FIM autour des sources:

- **`max_radius` négatif ou nul** : la FIM parcourt toute la grille.
- **`max_radius` strictement positif** : la valeur définit un seuil, et seuls les points situés à une distance inférieure à ce seuil des sources sont traités.

## Visualisation

Les scripts du dossier `visualization/` permettent d'analyser les résultats produits par la FIM. Ils génèrent des images au format `.png`.

- **`visu_fim.py`**: visualisation simple de la grille.

  ```bash
  python3 visualization/visu_fim.py
  ```

- **`visu_fim_speed.py`**: visualisation des résultats de la FIM et de la vitesse `F` en tout point de la grille.

  ```bash
  python3 visualization/visu_fim_speed.py
  ```

- **`closest_source.py`**: détermine et visualise, pour chaque point de la grille, la source la plus proche.

  ```bash
  python3 visualization/closest_source.py
  ```

- **`error_fim.py`** : visualisation de l'erreur absolue entre les résultats de la FIM et la distance euclidienne. Utilisable uniquement dans les cas avec une ou plusieurs sources et **aucun obstacle**.

  ```bash
  python3 visualization/error_fim.py
  ```

- **`benchmark.py`** : visualisation du temps de calcul de la FIM en fonction de la taille de la grille, avec et sans seuil.

  ```bash
  python3 visualization/benchmark.py
  ```