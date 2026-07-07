
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






pour mettre vitesse non cste:
modifier la fonction set_variable_speed du fichier generate_grid.c avec la fonction de vitesse voulue
dans le main du fichier generate_grid.c, mettre en commentaire la ligne "eikonal_grid_set_speed_constant(g,1.0) et décommenter la ligne set_variable_speed(g, cfg.n, cfg.m)




dans les params de fim_solve, si on veut parcourir toute la grille on met un max_radius négatif et si on veut un seuil on met la valeur du seuil