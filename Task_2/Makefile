# Makefile for N-Body OpenMP Simulation

CXX = g++
CXXFLAGS = -O3 -march=native -fopenmp -Wall -Wextra
LDFLAGS = -fopenmp -lm
TARGET = nbody_simulation
SRCDIR = src
BINDIR = bin
SRC = $(SRCDIR)/nbody_openmp.cpp

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(SRC)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(BINDIR)/*

run: $(BINDIR)/$(TARGET)
	@echo "Running with example parameters..."
	@./$(BINDIR)/$(TARGET) 1.0 example_input.txt

test: $(BINDIR)/$(TARGET)
	@echo "Creating test data..."
	@python3 generate_input.py 10 test_input.txt
	@echo "\nRunning test simulation..."
	@./$(BINDIR)/$(TARGET) 0.1 test_input.txt
	@echo "\nVisualizing results..."
	@python3 visualize_results.py trajectories.csv

.PHONY: all clean run test
