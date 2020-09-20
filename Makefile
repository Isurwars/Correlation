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

tests: test1 test2 test3 test4

test1: $(TST_DIR)/test_1/aPd.cell
	$(EXE) $(TST_DIR)/test_1/aPd.cell

test2: $(TST_DIR)/test_2/xSi.car
	$(EXE) $(TST_DIR)/test_2/xSi.car

test3: $(TST_DIR)/test_3/aPdH.car
	$(EXE) -i $(TST_DIR)/test_3/aPdH.dat $(TST_DIR)/test_3/aPdH.car

test4: $(TST_DIR)/test_4/graphite.car
	$(EXE) $(TST_DIR)/test_4/graphite.car

test5: $(TST_DIR)/test_5/aPdH_3144.car
	$(EXE) $(TST_DIR)/test_5/aPdH_3144.car

test6: $(TST_DIR)/test_6/aPdH_25152.car
	$(EXE) $(TST_DIR)/test_6/aPdH_25152.car
