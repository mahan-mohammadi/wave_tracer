# Wave renderer - Phase 1 numerical foundation.
# Usage:  make          build libwr.a and the tests
#         make test     build and run all tests
#         make clean
# Never add -ffast-math: compensated summation needs strict IEEE evaluation.
CC      ?= gcc
CFLAGS  ?= -std=c17 -O2 -Wall -Wextra -Wpedantic -fopenmp
CPPFLAGS += -Iinclude
LDLIBS  += -lm

LIB_SRC := $(wildcard src/*.c)
LIB_OBJ := $(LIB_SRC:src/%.c=build/%.o)
TESTS   := $(patsubst tests/%.c,build/%,$(wildcard tests/test_*.c))

all: build/libwr.a $(TESTS)

build:
	mkdir -p build

build/%.o: src/%.c | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

build/libwr.a: $(LIB_OBJ)
	ar rcs $@ $^

build/test_%: tests/test_%.c build/libwr.a tests/wr_test.h | build
	$(CC) $(CPPFLAGS) -Itests $(CFLAGS) $< build/libwr.a $(LDLIBS) -o $@

test: all
	@fail=0; for t in $(TESTS); do ./$$t || fail=1; done; \
	if [ $$fail -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "TESTS FAILED"; exit 1; fi

clean:
	rm -rf build

.PHONY: all test clean
