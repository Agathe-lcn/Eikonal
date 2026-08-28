CC = mpicc
CFLAGS = -O2 -Iinclude -I../../Eikonal2D/include -Wall
LDFLAGS = -lm

SRC_DIR = src
GRID_DIR = runners
BIN_DIR = bin
SEQ_DIR = ../Eikonal2D

SRC_COMMON = $(SRC_DIR)/FIM2D_mpi.c $(SRC_DIR)/FIM2D_io.c $(SRC_DIR)/config_mpi.c $(SEQ_DIR)/src/SolveEikonal2D.c $(SEQ_DIR)/src/FIM2D.c $(SEQ_DIR)/src/config.c

EXECS = $(BIN_DIR)/fim_mpi_run $(BIN_DIR)/benchmark_overlap $(BIN_DIR)/benchmark_optim

.PHONY: all dirs clean fim-mpi benchmark-overlap bencharm-optim

all: dirs $(EXECS)

dirs:
	mkdir -p $(BIN_DIR)

# Compilation des programmes MPI avec chargement du module
$(BIN_DIR)/fim_mpi_run: $(GRID_DIR)/fim_mpi_run.c $(SRC_COMMON) | dirs
	module load mpi/openmpi-x86_64 && $(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/benchmark_overlap: $(GRID_DIR)/benchmark_overlap.c $(SRC_COMMON) | dirs
	module load mpi/openmpi-x86_64 && $(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/benchmark_optim: $(GRID_DIR)/benchmark_optim_com.c $(SRC_COMMON) | dirs
	module load mpi/openmpi-x86_64 && $(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Cible pour exécuter les programmes
fim-mpi: $(BIN_DIR)/fim_mpi_run
	@if [ -z "$(NP)" ]; then \
		echo "Veuillez spécifier le nombre de processus avec NP=nombre"; \
		echo "Par exemple, make run-generate NP=4"; \
		echo "Ou avec un fichier de configuration en argument"; \
		echo "Par exemple, make run-generate NP=4 ARGS='config.txt"; \
	elif [ -z "$(ARGS)" ]; then \
		module load mpi/openmpi-x86_64 && mpirun -np $(NP) ./$(BIN_DIR)/fim_mpi_run config.txt; \
	else \
		module load mpi/openmpi-x86_64 && mpirun -np $(NP) ./$(BIN_DIR)/fim_run_mpi $(ARGS); \
	fi

benchmark-overlap: $(BIN_DIR)/benchmark_overlap
	@if [ -z "$(NP)" ]; then \
		echo "Veuillez spécifier le nombre de processus avec NP=nombre"; \
		echo "Par exemple, make run-generate NP=4"; \
		echo "Ou avec un fichier de configuration en argument"; \
		echo "Par exemple, make run-generate NP=4 ARGS='config.txt"; \
	elif [ -z "$(ARGS)" ]; then \
		module load mpi/openmpi-x86_64 && mpirun -np $(NP) ./$(BIN_DIR)/benchmark_overlap config.txt; \
	else \
		module load mpi/openmpi-x86_64 && mpirun -np $(NP) ./$(BIN_DIR)/benchmark_overlap $(ARGS); \
	fi

benchmark-optim: $(BIN_DIR)/benchmark_optim
	@if [ -z "$(NP)" ]; then \
		echo "Veuillez spécifier le nombre de processus avec NP=nombre"; \
		echo "Par exemple, make run-generate NP=4"; \
		echo "Ou avec un fichier de configuration en argument"; \
		echo "Par exemple, make run-generate NP=4 ARGS='config.txt"; \
	elif [ -z "$(ARGS)" ]; then \
		module load mpi/openmpi-x86_64 && mpirun -np $(NP) ./$(BIN_DIR)/benchmark_optim config.txt; \
	else \
		module load mpi/openmpi-x86_64 && mpirun -np $(NP) ./$(BIN_DIR)/benchmark_optim $(ARGS); \
	fi

clean:
	rm -rf $(BIN_DIR)
	rm -f coords_source.txt matrix_fim.txt.chkpt result_rank0.txt result_rank1.txt result_rank2.txt result_rank3.txt