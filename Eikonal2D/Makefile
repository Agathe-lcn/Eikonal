CC = gcc
CFLAGS = -O2 -Iinclude -Wall
LDFLAGS = -lm

SRC_DIR = src
GRID_DIR = grids
BIN_DIR = bin

SRC_COMMON = $(SRC_DIR)/FIM2D.c $(SRC_DIR)/SolveEikonal2D.c $(SRC_DIR)/config.c

EXECS = $(BIN_DIR)/generate_grid $(BIN_DIR)/1_source $(BIN_DIR)/5_sources_1_wall $(BIN_DIR)/benchmark_fim

generate: $(BIN_DIR)/generate_grid
1source: $(BIN_DIR)/1_source
5sources: $(BIN_DIR)/5_sources_1_wall
benchmark: $(BIN_DIR)/benchmark_fim

.PHONY: all dirs clean run benchmark

all: dirs $(EXECS)

dirs:
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/generate_grid: $(GRID_DIR)/generate_grid.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/1_source: $(GRID_DIR)/1_source.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/5_sources_1_wall: $(GRID_DIR)/5_sources_1_wall.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/benchmark_fim: $(GRID_DIR)/benchmark_fim.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Cibles pour exécuter chaque programme
run-generate: $(BIN_DIR)/generate_grid
	./$(BIN_DIR)/generate_grid config.txt

run-1source: $(BIN_DIR)/1_source
	./$(BIN_DIR)/1_source config.txt

run-5sources: $(BIN_DIR)/5_sources_1_wall
	./$(BIN_DIR)/5_sources_1_wall config.txt

# Exécution du benchmark
benchmark: $(BIN_DIR)/benchmark_fim
	./$(BIN_DIR)/benchmark_fim

# Exécution par défaut (generate_grid)
run: all
	./$(BIN_DIR)/generate_grid config.txt

clean:
	rm -rf $(BIN_DIR)
	rm -f coords_source.txt matrix_fim.txt source_tags.txt speed.txt results_fim.txt