#!/bin/bash

echo "=== Building N-Body OpenMP Simulation ==="

# Создаем директории
mkdir -p bin

# Проверяем компилятор
echo "Checking compiler..."
if ! command -v g++ &> /dev/null; then
    echo "Error: g++ not found. Installing..."
    sudo apt update
    sudo apt install g++ -y
fi

# Проверяем OpenMP
echo "Checking OpenMP support..."
g++ -fopenmp -dM -E - < /dev/null | grep -i openmp > /dev/null
if [ $? -ne 0 ]; then
    echo "Warning: OpenMP not supported by compiler"
    echo "Trying to install g++ with OpenMP support..."
    sudo apt install g++-11 -y
    if command -v g++-11 &> /dev/null; then
        echo "Using g++-11"
        sed -i 's/CXX = g++/CXX = g++-11/' Makefile
    fi
fi

# Компилируем
echo -e "\nCompiling..."
make clean
make

if [ $? -eq 0 ]; then
    echo -e "\nBuild successful!"
    echo "Executable: bin/nbody_simulation"
    echo -e "\nUsage:"
    echo "  ./bin/nbody_simulation <tend> <input_file>"
    echo "  ./run_simulation.sh <tend> <input_file>"
    echo -e "\nTo test:"
    echo "  make test"
else
    echo -e "\nBuild failed!"
    exit 1
fi
