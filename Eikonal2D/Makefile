CC = gcc
CFLAGS = -O2 -Iinclude -Wall
LDFLAGS = -lm

SRC_DIR = src
GRID_DIR = runners
BIN_DIR = bin

SRC_COMMON = $(SRC_DIR)/FIM2D.c $(SRC_DIR)/SolveEikonal2D.c $(SRC_DIR)/config.c

EXECS = $(BIN_DIR)/fim_run $(BIN_DIR)/benchmark_fim

fim: $(BIN_DIR)/fim_run
benchmark: $(BIN_DIR)/benchmark_fim

.PHONY: all dirs clean run benchmark

all: dirs $(EXECS)

dirs:
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/fim_run: $(GRID_DIR)/fim_run.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/benchmark_fim: $(GRID_DIR)/benchmark_fim.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Cible pour exécuter le programme
fim: $(BIN_DIR)/fim_run
	@if [ -z "$(ARGS)" ]; then \
		./$(BIN_DIR)/fim_run config.txt; \
	else \
		./$(BIN_DIR)/fim_run $(ARGS); \
	fi

# Exécution du benchmark
benchmark: $(BIN_DIR)/benchmark_fim
	./$(BIN_DIR)/benchmark_fim


clean:
	rm -rf $(BIN_DIR)
	rm -f coords_source.txt matrix_fim.txt source_tags.txt speed.txt results_fim.txt