#!/bin/bash

# Скрипт запуска симуляции N тел

if [ $# -lt 2 ]; then
    echo "Usage: $0 <tend> <input_file>"
    echo "Example: $0 10.0 data/input.txt"
    echo ""
    echo "To generate input data:"
    echo "  python3 generate_input.py <num_particles> <output_file>"
    echo "  python3 generate_input.py 100 input.txt"
    exit 1
fi

TEND=$1
INPUT_FILE=$2

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    echo "Generate one with: python3 generate_input.py <N> $INPUT_FILE"
    exit 1
fi

# Проверяем, скомпилирована ли программа
if [ ! -f "bin/nbody_simulation" ]; then
    echo "Program not compiled. Building..."
    ./build.sh
    if [ ! -f "bin/nbody_simulation" ]; then
        echo "Build failed!"
        exit 1
    fi
fi

echo "=== Starting N-Body OpenMP Simulation ==="
echo "End time: $TEND seconds"
echo "Input file: $INPUT_FILE"
echo ""

# Запускаем симуляцию
time ./bin/nbody_simulation "$TEND" "$INPUT_FILE"

echo ""
echo "=== Simulation Complete ==="
echo "Results saved to: trajectories.csv"
echo ""
echo "To visualize results:"
echo "  python3 visualize_results.py trajectories.csv"
echo ""
echo "For different visualization types:"
echo "  python3 visualize_results.py trajectories.csv all"
echo "  python3 visualize_results.py trajectories.csv simple"
echo "  python3 visualize_results.py trajectories.csv animation"
