ifeq ($(OS),Windows_NT)
RM = del /q
else
RM = rm -rf
endif

SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
TST_DIR := test

EXE := $(BIN_DIR)/correlation.exe
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CC			 := g++
CPPFLAGS := -Iinclude
CFLAGS   := -Wall
LDFLAGS  := -Llib
LDLIBS   := -lm -static-libgcc -static-libstdc++

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir $@

clean:
	$(RM) $(OBJ_DIR)

tests: test1 test2 test3

test1: $(TST_DIR)/test_1/aPd.cell
	$(EXE) $(TST_DIR)/test_1/aPd.cell

test2: $(TST_DIR)/test_2/Si.car
	$(EXE) $(TST_DIR)/test_2/Si.car

test3: $(TST_DIR)/test_3/PdH.car
	$(EXE) $(TST_DIR)/test_3/PdH.car
