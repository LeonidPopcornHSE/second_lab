#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def plot_simple_trajectories(df, output_file="trajectories_plot.png"):
    """Простая визуализация траекторий"""
    n_particles = (len(df.columns) - 1) // 2
    
    plt.figure(figsize=(12, 10))
    
    # Ограничиваем количество частиц для читаемости
    max_particles = min(n_particles, 20)
    colors = plt.cm.rainbow(np.linspace(0, 1, max_particles))
    
    for i in range(max_particles):
        x_col = f'x{i+1}'
        y_col = f'y{i+1}'
        
        plt.plot(df[x_col], df[y_col], '-', color=colors[i], alpha=0.7, linewidth=1.5)
        plt.plot(df[x_col].iloc[0], df[y_col].iloc[0], 'o', 
                color=colors[i], markersize=10, label=f'Particle {i+1}')
        plt.plot(df[x_col].iloc[-1], df[y_col].iloc[-1], 's', 
                color=colors[i], markersize=8, alpha=0.8)
    
    plt.xlabel('X coordinate (m)', fontsize=12)
    plt.ylabel('Y coordinate (m)', fontsize=12)
    plt.title(f'Trajectories of {n_particles} Particles', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"Plot saved to {output_file}")

def plot_energy(df, output_file="energy_plot.png"):
    """Построение энергии системы"""
    # Для простоты будем считать только кинетическую энергию
    n_particles = (len(df.columns) - 1) // 2
    
    # Нужны массы частиц - их нет в выводе, используем фиктивные
    masses = [1.989e30] + [10**24] * (n_particles - 1)  # Солнце + планеты
    
    kinetic_energy = []
    for idx in range(len(df)):
        total_ke = 0
        for i in range(n_particles):
            vx_col = None  # Скоростей нет в выводе, пропускаем
            # Для реального расчета нужны скорости
        kinetic_energy.append(total_ke)
    
    plt.figure(figsize=(10, 6))
    plt.plot(df['t'], kinetic_energy, 'b-', linewidth=2)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Kinetic Energy (J)', fontsize=12)
    plt.title('System Kinetic Energy vs Time', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=150)
    plt.show()

def create_summary(df, csv_file):
    """Создание текстового отчета"""
    n_particles = (len(df.columns) - 1) // 2
    
    print("=" * 60)
    print("SIMULATION SUMMARY")
    print("=" * 60)
    print(f"Results file: {csv_file}")
    print(f"Number of particles: {n_particles}")
    print(f"Simulation time: {df['t'].iloc[-1]:.2e} s")
    print(f"Time steps recorded: {len(df)}")
    
    print("\nInitial positions (first 5 particles):")
    for i in range(min(5, n_particles)):
        x0 = df[f'x{i+1}'].iloc[0]
        y0 = df[f'y{i+1}'].iloc[0]
        print(f"  Particle {i+1}: x = {x0:.2e} m, y = {y0:.2e} m")
    
    print("\nFinal positions (first 5 particles):")
    for i in range(min(5, n_particles)):
        xf = df[f'x{i+1}'].iloc[-1]
        yf = df[f'y{i+1}'].iloc[-1]
        dx = xf - df[f'x{i+1}'].iloc[0]
        dy = yf - df[f'y{i+1}'].iloc[0]
        displacement = np.sqrt(dx**2 + dy**2)
        print(f"  Particle {i+1}: x = {xf:.2e} m, y = {yf:.2e} m")
        print(f"               Displacement: {displacement:.2e} m")
    
    print("\nPosition ranges:")
    for i in range(min(3, n_particles)):
        x_min = df[f'x{i+1}'].min()
        x_max = df[f'x{i+1}'].max()
        y_min = df[f'y{i+1}'].min()
        y_max = df[f'y{i+1}'].max()
        print(f"  Particle {i+1}: x ∈ [{x_min:.2e}, {x_max:.2e}], "
              f"y ∈ [{y_min:.2e}, {y_max:.2e}]")
    
    print("=" * 60)
def main():
    if len(sys.argv) < 2:
        print("Usage: python3 visualize_results.py <csv_file> [plot_type]")
        print("\nPlot types:")
        print("  simple    - Simple trajectory plot (default)")
        print("  summary   - Text summary only")
        print("  all       - All plots")
        print("\nExample: python3 visualize_results.py trajectories.csv simple")
        return
    
    csv_file = sys.argv[1]
    plot_type = sys.argv[2] if len(sys.argv) > 2 else "simple"
    
    if not os.path.exists(csv_file):
        print(f"Error: File '{csv_file}' not found!")
        return
    
    print(f"Reading data from {csv_file}...")
    df = pd.read_csv(csv_file)
    
    print(f"   Loaded: {len(df)} time steps, {(len(df.columns)-1)//2} particles")
    
    # Создаем директорию для графиков
    os.makedirs("plots", exist_ok=True)
    
    if plot_type == "simple" or plot_type == "all":
        plot_simple_trajectories(df, "plots/trajectories.png")
    
    if plot_type == "summary" or plot_type == "all":
        create_summary(df, csv_file)
    
    print("\nVisualization complete!")
    print(f"   Plots saved to 'plots/' directory")

if __name__ == "__main__":
    main()
