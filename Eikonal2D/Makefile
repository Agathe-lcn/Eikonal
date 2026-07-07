CC = gcc
CFLAGS = -O2 -Iinclude -Wall
LDFLAGS = -lm

SRC_DIR = src
GRID_DIR = grids
BIN_DIR = bin

SRC_COMMON = $(SRC_DIR)/FIM2D.c $(SRC_DIR)/SolveEikonal2D.c

GRID_BINS = $(BIN_DIR)/generate_grid $(BIN_DIR)/1_source $(BIN_DIR)/5_sources_1_wall

.PHONY: all dirs tests clean run

all: dirs $(BIN_DIR)/generate_grid

dirs:
	mkdir -p $(BIN_DIR) $(OUTPUT_DIR)

$(BIN_DIR)/generate_grid: $(GRID_DIR)/generate_grid.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/1_source: $(GRID_DIR)/1_source.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/5_sources_1_wall: $(GRID_DIR)/5_sources_1_wall.c $(SRC_COMMON) | dirs
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

run: all
	./$(BIN_DIR)/generate_grid config.txt

clean:
	rm -rf $(BIN_DIR)
	rm -f coords_source.txt matrix_fim.txt source_tags.txt speed.txt