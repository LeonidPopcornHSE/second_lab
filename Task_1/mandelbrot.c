#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define MAX_ITER 1000

int in_mandelbrot(double cx, double cy) {
    double zx = 0.0, zy = 0.0;
    for (int i = 0; i < MAX_ITER; i++) {
        double zx_new = zx * zx - zy * zy + cx;
        double zy_new = 2.0 * zx * zy + cy;
        zx = zx_new;
        zy = zy_new;

        if (zx * zx + zy * zy >= 4.0)
            return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s nthreads npoints\n", argv[0]);
        return 1;
    }

    int nthreads = atoi(argv[1]);
    int npoints  = atoi(argv[2]);

    omp_set_num_threads(nthreads);

    double xmin = -2.0, xmax = 1.0;
    double ymin = -1.5, ymax = 1.5;

    FILE *out = fopen("mandelbrot.csv", "w");
    if (!out) {
        perror("fopen");
        return 1;
    }

    fprintf(out, "x,y\n");

    #pragma omp parallel
    {
        unsigned int seed = omp_get_thread_num();

        #pragma omp for schedule(static)
        for (int i = 0; i < npoints; i++) {
            double x = xmin + (xmax - xmin) * rand_r(&seed) / RAND_MAX;
            double y = ymin + (ymax - ymin) * rand_r(&seed) / RAND_MAX;

            if (in_mandelbrot(x, y)) {
                #pragma omp critical
                {
                    fprintf(out, "%.8f,%.8f\n", x, y);
                }
            }
        }
    }

    fclose(out);
    return 0;
}
