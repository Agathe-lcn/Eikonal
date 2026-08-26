CC = mpicc
CFLAGS = -O2 -Iinclude -I../../Eikonal2D/include -Wall
LDFLAGS = -lm

SRC_DIR = src
GRID_DIR = grids
BIN_DIR = bin
SEQ_DIR = ../Eikonal2D

SRC_COMMON = $(SRC_DIR)/FIM2D_mpi.c $(SRC_DIR)/FIM2D_io.c $(SRC_DIR)/config_mpi.c $(SEQ_DIR)/src/SolveEikonal2D.c $(SEQ_DIR)/src/FIM2D.c $(SEQ_DIR)/src/config.c

EXECS = $(BIN_DIR)/generate_grid_mpi

.PHONY: all dirs clean run-generate benchmark

all: dirs $(EXECS)

dirs:
	mkdir -p $(BIN_DIR)

# Compilation des programmes MPI avec chargement du module
$(BIN_DIR)/generate_grid_mpi: $(GRID_DIR)/generate_grid_mpi.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Cible pour exécuter le programme
run-generate: $(BIN_DIR)/generate_grid_mpi
	@if [ -z "$(NP)" ]; then \
		echo "Veuillez spécifier le nombre de processus avec NP=nombre"; \
		echo "Par exemple, make run-generate NP=4"; \
		echo "Ou avec un fichier de configuration en argument"; \
		echo "Par exemple, make run-generate NP=4 ARGS='config.txt"; \
	elif [ -z "$(ARGS)" ]; then \
		mpirun -np $(NP) ./$(BIN_DIR)/generate_grid_mpi config.txt; \
	else \
		mpirun -np $(NP) ./$(BIN_DIR)/generate_grid_mpi $(ARGS); \
	fi

clean:
	rm -rf $(BIN_DIR)
	rm -f coords_source.txt matrix_fim.txt.chkpt result_rank0.txt result_rank1.txt result_rank2.txt result_rank3.txt
