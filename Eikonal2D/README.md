Pour avoir le plot avec 5 sources et un mur vertical au centre:
dans dossier Eikonal2D:
gcc -o generate_grid grids/5_sources_1_wall.c src/SolveEikonal2D.c src/FIM2D.c -lm
./generate_grid
python3 visualization/test.py

pour avoir plot avec 1 source et l'erreur par rapport à la distance euclidienne:
dans dossier Eikonal2D:
gcc -o generate_grid grids/1_source.c src/SolveEikonal2D.c src/FIM2D.c -lm 
./generate_grid
python3 visualization.visu_fim.py