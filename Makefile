ifeq ($(OS),Windows_NT)
RM = del
RMDIR = rmdir /S /Q
else
RM = rm -rf
RMDIR = rm -rf
endif

SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
TST_DIR := test

EXE := $(BIN_DIR)/correlation.exe
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CC			 := g++
CPPFLAGS := -Iinclude -std=c++17
CFLAGS   := -Wall
LDFLAGS  := -Llib
LDLIBS   := -lm -static-libstdc++

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir $@

clean:
	$(RMDIR) $(OBJ_DIR)

tests: test1 test2 test3 test4 test5

test1: $(TST_DIR)/test_1/xSi.car
	$(EXE) $(TST_DIR)/test_1/xSi.car

test2: $(TST_DIR)/test_2/Graphene.car
	$(EXE) $(TST_DIR)/test_2/Graphene.car

test3: $(TST_DIR)/test_3/aPd.cell
	$(EXE) $(TST_DIR)/test_3/aPd.cell

test4: $(TST_DIR)/test_4/aPdH.dat
	$(EXE) -i $(TST_DIR)/test_4/aPdH.txt $(TST_DIR)/test_4/aPdH.dat

test5: $(TST_DIR)/test_5/lBi.car
	$(EXE) $(TST_DIR)/test_5/lBi.car

test6: $(TST_DIR)/test_6/aPdH_3144.car
	$(EXE) $(TST_DIR)/test_6/aPdH_3144.car
