#!/usr/bin/env python3
import random
import sys
import numpy as np

def generate_input(filename, n=100):
    """Генерация тестовых данных для N тел"""
    with open(filename, 'w') as f:
        f.write(f"{n}\n")
        
        # Солнце в центре
        f.write(f"1.989e30 0.0 0.0 0.0 0.0 0.0 0.0\n")
        
        # Остальные частицы
        for i in range(1, n):
            # Массы в диапазоне 1e22 - 1e26 кг (астероиды, маленькие планеты)
            mass = 10**random.uniform(22, 26)
            
            # Позиции на кольце вокруг Солнца
            angle = random.uniform(0, 2 * np.pi)
            distance = random.uniform(1e11, 5e11)  # 0.1-5 а.е.
            
            x = distance * np.cos(angle)
            y = distance * np.sin(angle)
            z = random.uniform(-1e10, 1e10)  # Немного вне плоскости
            
            # Орбитальные скорости (приблизительно круговые орбиты)
            orbital_speed = np.sqrt(6.67430e-11 * 1.989e30 / distance)
            vx = -orbital_speed * np.sin(angle) * random.uniform(0.9, 1.1)
            vy = orbital_speed * np.cos(angle) * random.uniform(0.9, 1.1)
            vz = random.uniform(-1e3, 1e3)
            
            f.write(f"{mass:.6e} {x:.6e} {y:.6e} {z:.6e} ")
            f.write(f"{vx:.6e} {vy:.6e} {vz:.6e}\n")
    
    print(f"Generated input file: {filename}")
    print(f"   Number of particles: {n}")
    print(f"   First particle (Sun): mass = 1.989e30 kg")
    print(f"   Other particles: masses ~10^22-10^26 kg")
    print(f"   Orbital distances: 0.1-5 AU")

def generate_simple_test(filename="example_input.txt"):
    """Создание простого тестового файла"""
    with open(filename, 'w') as f:
        f.write("5\n")
        # Солнце и 4 планеты
        f.write("1.989e30 0.0 0.0 0.0 0.0 0.0 0.0\n")
        f.write("5.972e24 1.5e11 0.0 0.0 0.0 29780.0 0.0\n")  # Земля
        f.write("4.867e24 0.0 1.08e11 0.0 -35020.0 0.0 0.0\n")  # Венера
        f.write("6.39e23 -2.28e11 0.0 0.0 0.0 -24130.0 0.0\n")  # Марс
        f.write("1.898e27 0.0 -7.78e11 0.0 13070.0 0.0 0.0\n")  # Юпитер
    
    print(f"Generated simple test file: {filename}")
    print("   Contains: Sun + Earth + Venus + Mars + Jupiter")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        # По умолчанию создаем простой тест
        generate_simple_test()
    elif len(sys.argv) == 2:
        n = int(sys.argv[1])
        generate_input(f"input_{n}.txt", n)
    elif len(sys.argv) == 3:
        n = int(sys.argv[1])
        filename = sys.argv[2]
        generate_input(filename, n)
    else:
        print("Usage:")
        print("  python3 generate_input.py                    # Create simple test (5 particles)")
        print("  python3 generate_input.py <N>               # Create N particles")
        print("  python3 generate_input.py <N> <filename>    # Create N particles in file")
        print("\nExamples:")
        print("  python3 generate_input.py 100               # 100 particles")
        print("  python3 generate_input.py 50 test.txt       # 50 particles in test.txt")
