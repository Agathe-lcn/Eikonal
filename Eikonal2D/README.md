
##### A ECRIRE PROPREMENT ####




Pour avoir le plot avec 5 sources et un mur vertical au centre:
dans dossier Eikonal2D:
gcc -O2 -o generate_grid grids/5_sources_1_wall.c src/SolveEikonal2D.c src/FIM2D.c -lm
./generate_grid
python3 visualization/visu_fim.py

pour avoir plot avec 1 source et l'erreur par rapport à la distance euclidienne:
dans dossier Eikonal2D:
gcc -O2 -o error grids/1_source.c src/SolveEikonal2D.c src/FIM2D.c -lm 
./error
python3 visualization/error_fim_1_source.py

pr avoir plot avec 1 source au centre:
gcc -O2 -o generate_grid grids/1_source.c src/SolveEikonal2D.c src/FIM2D.c -lm 
./generate_grid
python3 visualization/visu_fim.py