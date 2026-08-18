CC = gcc
CFLAGS = -O2 -Iinclude -Wall
LDFLAGS = -lm

SRC_DIR = src
GRID_DIR = grids
BIN_DIR = bin

SRC_COMMON = $(SRC_DIR)/FIM2D.c $(SRC_DIR)/SolveEikonal2D.c $(SRC_DIR)/config.c

EXECS = $(BIN_DIR)/generate_grid $(BIN_DIR)/benchmark_fim

generate: $(BIN_DIR)/generate_grid
benchmark: $(BIN_DIR)/benchmark_fim

.PHONY: all dirs clean run benchmark

all: dirs $(EXECS)

dirs:
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/generate_grid: $(GRID_DIR)/generate_grid.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/benchmark_fim: $(GRID_DIR)/benchmark_fim.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Cible pour exécuter le programme
run-generate: $(BIN_DIR)/generate_grid
	@if [ -z "$(ARGS)" ]; then \
		./$(BIN_DIR)/generate_grid config.txt; \
	else \
		./$(BIN_DIR)/generate_grid $(ARGS); \
	fi

# Exécution du benchmark
benchmark: $(BIN_DIR)/benchmark_fim
	./$(BIN_DIR)/benchmark_fim


clean:
	rm -rf $(BIN_DIR)
	rm -f coords_source.txt matrix_fim.txt source_tags.txt speed.txt results_fim.txt