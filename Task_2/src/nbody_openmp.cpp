#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <omp.h>

const double G = 6.67430e-11;  // гравитационная постоянная

struct Particle {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    
    Particle() : mass(0), x(0), y(0), z(0), 
                 vx(0), vy(0), vz(0), ax(0), ay(0), az(0) {}
};

// Чтение входных данных
bool read_input(const std::string& filename, 
                std::vector<Particle>& particles, 
                int& n) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    file >> n;
    particles.resize(n);
    
    for (int i = 0; i < n; i++) {
        file >> particles[i].mass
             >> particles[i].x >> particles[i].y >> particles[i].z
             >> particles[i].vx >> particles[i].vy >> particles[i].vz;
        particles[i].ax = particles[i].ay = particles[i].az = 0.0;
    }
    
    file.close();
    return true;
}

// Запись результатов в CSV
void write_output(const std::string& filename, 
                  const std::vector<double>& time_points,
                  const std::vector<std::vector<Particle>>& trajectories,
                  int n) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create output file" << std::endl;
        return;
    }
    
    // Заголовок
    file << "t";
    for (int i = 0; i < n; i++) {
        file << ",x" << i+1 << ",y" << i+1;
    }
    file << "\n";
    
    // Данные
    for (size_t t_idx = 0; t_idx < time_points.size(); t_idx++) {
        file << std::scientific << std::setprecision(6) << time_points[t_idx];
        for (int i = 0; i < n; i++) {
            file << "," << trajectories[t_idx][i].x 
                 << "," << trajectories[t_idx][i].y;
        }
        file << "\n";
    }
    
    file.close();
}

// Оптимизированное вычисление ускорений с OpenMP
void compute_accelerations(std::vector<Particle>& particles, int n) {
    // Обнуляем ускорения
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        particles[i].ax = 0.0;
        particles[i].ay = 0.0;
        particles[i].az = 0.0;
    }
    
    // Вычисляем силы попарно
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        double ax_i = 0.0, ay_i = 0.0, az_i = 0.0;
        double xi = particles[i].x;
        double yi = particles[i].y;
        double zi = particles[i].z;
        
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dx = particles[j].x - xi;
                double dy = particles[j].y - yi;
                double dz = particles[j].z - zi;
                
                double dist_sq = dx*dx + dy*dy + dz*dz + 1e-10; // Регуляризация
                double dist = sqrt(dist_sq);
                double force_magnitude = G * particles[j].mass / (dist_sq * dist);
                
                ax_i += force_magnitude * dx;
                ay_i += force_magnitude * dy;
                az_i += force_magnitude * dz;
            }
        }
        
        particles[i].ax = ax_i;
        particles[i].ay = ay_i;
        particles[i].az = az_i;
    }
}
// Метод Верле (Velocity Verlet) для интегрирования
void verlet_integration(std::vector<Particle>& particles, double dt, int n) {
    // Шаг 1: обновление позиций
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        particles[i].x += particles[i].vx * dt + 0.5 * particles[i].ax * dt * dt;
        particles[i].y += particles[i].vy * dt + 0.5 * particles[i].ay * dt * dt;
        particles[i].z += particles[i].vz * dt + 0.5 * particles[i].az * dt * dt;
    }
    
    // Сохраняем старые ускорения
    std::vector<double> old_ax(n), old_ay(n), old_az(n);
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        old_ax[i] = particles[i].ax;
        old_ay[i] = particles[i].ay;
        old_az[i] = particles[i].az;
    }
    
    // Шаг 2: вычисляем новые ускорения
    compute_accelerations(particles, n);
    
    // Шаг 3: обновление скоростей
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        particles[i].vx += 0.5 * (particles[i].ax + old_ax[i]) * dt;
        particles[i].vy += 0.5 * (particles[i].ay + old_ay[i]) * dt;
        particles[i].vz += 0.5 * (particles[i].az + old_az[i]) * dt;
    }
}

int main(int argc, char* argv[]) {
    // Установка количества потоков по умолчанию
    omp_set_num_threads(omp_get_max_threads());
    
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <tend> <input_file>" << std::endl;
        std::cerr << "Example: " << argv[0] << " 10.0 input_data.txt" << std::endl;
        return 1;
    }
    
    double tend = std::stod(argv[1]);
    std::string input_file = argv[2];
    
    int n;
    std::vector<Particle> particles;
    
    if (!read_input(input_file, particles, n)) {
        return 1;
    }
    
    std::cout << "=== N-Body Simulation (OpenMP) ===" << std::endl;
    std::cout << "Number of particles: " << n << std::endl;
    std::cout << "Simulation time: " << tend << " seconds" << std::endl;
    std::cout << "Number of OpenMP threads: " << omp_get_max_threads() << std::endl;
    
    // Параметры симуляции
    double dt = 0.01;  // шаг по времени
    int steps = static_cast<int>(tend / dt);
    int output_freq = std::max(1, steps / 1000);  // выводим примерно 1000 точек
    
    std::vector<double> time_points;
    std::vector<std::vector<Particle>> trajectories;
    
    // Начальные ускорения
    compute_accelerations(particles, n);
    
    std::cout << "\nStarting simulation..." << std::endl;
    std::cout << "Total steps: " << steps << std::endl;
    std::cout << "Output frequency: every " << output_freq << " steps" << std::endl;
    
    // Основной цикл симуляции
    double start_time = omp_get_wtime();
    
    for (int step = 0; step <= steps; step++) {
        double t = step * dt;
        
        // Сохраняем состояние
        if (step % output_freq == 0) {
            time_points.push_back(t);
            trajectories.push_back(particles);
            
            // Прогресс
            if (step % (steps / 10) == 0 && step > 0) {
                double progress = 100.0 * step / steps;
                double elapsed = omp_get_wtime() - start_time;
                std::cout << "Progress: " << std::fixed << std::setprecision(1) 
                          << progress << "%, Time: " << t << "s, "
                          << "Elapsed: " << elapsed << "s" << std::endl;
            }
        }
        
        // Интегрируем
        verlet_integration(particles, dt, n);
    }
    
    double end_time = omp_get_wtime();
    double total_time = end_time - start_time;
    
    // Запись результатов
    std::string output_file = "trajectories.csv";
    write_output(output_file, time_points, trajectories, n);
    
    // Статистика
    std::cout << "\n=== Simulation Complete ===" << std::endl;
    std::cout << "Total simulation time: " << total_time << " seconds" << std::endl;
    std::cout << "Performance: " << steps * n * n / total_time / 1e9 << " G interactions/sec" << std::endl;
    std::cout << "Output points: " << time_points.size() << std::endl;
    std::cout << "Results saved to: " << output_file << std::endl;
return 0;
}
