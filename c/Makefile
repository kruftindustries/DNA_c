# Build for the C port. No external libraries: only libc and libm.
#
#   make            build bin/genetic-health
#   make test       build and run the C unit tests
#   make debug      build with sanitizers
#   make clean

CC      ?= cc
CSTD    ?= -std=c11
WARN    := -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes \
           -Wmissing-prototypes -Wconversion -Wsign-conversion \
           -Wno-unused-parameter
OPT     ?= -O2
# SAN must reach both the compile and the link step, or the sanitizer runtime
# is never linked in.
SAN     ?=
VERSION ?= $(shell sed -n 's/^VERSION = //p' engine.pro)
CFLAGS  += $(CSTD) $(WARN) $(OPT) $(SAN) -Iinclude -DGH_VERSION='"$(VERSION)"'
LDFLAGS += $(SAN)
LDLIBS  += -lm

SRC_DIR   := src
OBJ_DIR   := build
BIN_DIR   := bin
TEST_DIR  := tests

# Library sources: everything except the CLI entry point.
LIB_SRCS := $(filter-out $(SRC_DIR)/main.c,$(wildcard $(SRC_DIR)/*.c))
LIB_OBJS := $(LIB_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
MAIN_OBJ := $(OBJ_DIR)/main.o

TEST_SRCS := $(wildcard $(TEST_DIR)/*.c)
TEST_BINS := $(TEST_SRCS:$(TEST_DIR)/%.c=$(BIN_DIR)/%)

TARGET := $(BIN_DIR)/genetic-health

.PHONY: all test debug clean

all: $(TARGET)

$(TARGET): $(LIB_OBJS) $(MAIN_OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BIN_DIR)/%: $(TEST_DIR)/%.c $(LIB_OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LIB_OBJS) $(LDLIBS) -o $@

$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

test: $(TEST_BINS)
	@fail=0; \
	for t in $(TEST_BINS); do \
		echo "--- $$t"; \
		./$$t || fail=1; \
	done; \
	exit $$fail

# Sanitizer build: catches the classes of bug a manual port is most likely to
# introduce (bounds, lifetime, integer overflow).
debug:
	$(MAKE) clean
	$(MAKE) OPT="-O1 -g" SAN="-fsanitize=address,undefined -fno-omit-frame-pointer" all test

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

-include $(OBJ_DIR)/*.d

