#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#define MAX_CITIES 100

typedef struct {
    double x;
    double y;
} City;

double calculateDistance(City city1, City city2) {
    double dx = city1.x - city2.x;
    double dy = city1.y - city2.y;
    return sqrt(dx * dx + dy * dy);
}

double calculateTotalDistance(int* path, City* cities, int numCities) {
    double totalDistance = 0.0;
    for (int i = 0; i < numCities - 1; i++) {
        int city1 = path[i];
        int city2 = path[i + 1];
        totalDistance += calculateDistance(cities[city1], cities[city2]);
    }
    int lastCity = path[numCities - 1];
    int firstCity = path[0];
    totalDistance += calculateDistance(cities[lastCity], cities[firstCity]);
    return totalDistance;
}

void twoOptSwap(int* path, int i, int j) {
    while (i < j) {
        int temp = path[i];
        path[i] = path[j];
        path[j] = temp;
        i++;
        j--;
    }
}

int main() {
    int numCities;
    printf("Enter the number of cities: ");
    scanf("%d", &numCities);

    if (numCities <= 0 || numCities > MAX_CITIES) {
        printf("Invalid number of cities. Exiting...\n");
        return 1;
    }

    City cities[MAX_CITIES];
    srand(time(NULL));

    // Generate random cities
    for (int i = 0; i < numCities; i++) {
        cities[i].x = (double)rand() / RAND_MAX;
        cities[i].y = (double)rand() / RAND_MAX;
    }

    // Select a random starting city
    int startCity = rand() % numCities;

    // Print the input cities
    printf("Input Cities:\n");
    for (int i = 0; i < numCities; i++) {
        printf("City %d: (%.2f, %.2f)\n", i + 1, cities[i].x, cities[i].y);
    }
    printf("\n");

    // Solve TSP using serial algorithm
    printf("----- Serial Algorithm -----\n");
    int* path = (int*)malloc(numCities * sizeof(int));
    for (int i = 0; i < numCities; i++) {
        path[i] = i;
    }
    double serialStart = omp_get_wtime();
    double bestSerialDistance = calculateTotalDistance(path, cities, numCities);
    for (int i = 1; i < numCities - 1; i++) {
        for (int j = i + 1; j < numCities; j++) {
            twoOptSwap(path, i, j);
            double newDistance = calculateTotalDistance(path, cities, numCities);
            if (newDistance < bestSerialDistance) {
                bestSerialDistance = newDistance;
            }
        }
    }
    double serialEnd = omp_get_wtime();
    double serialTime = serialEnd - serialStart;
    printf("Serial Minimum Distance: %.2f\n", bestSerialDistance);
        printf("Serial Elapsed Time: %.6f seconds\n\n", serialTime);

    // Solve TSP using parallel algorithm
    printf("----- Parallel Algorithm -----\n");
    srand(time(NULL)); // Set the seed for random number generation
    double bestParallelDistance = bestSerialDistance;
    int* parallelPath = (int*)malloc(numCities * sizeof(int));
    memcpy(parallelPath, path, numCities * sizeof(int));
    double parallelStart = omp_get_wtime();
    int* newPath; // Declare newPath outside the parallel region

    #pragma omp parallel num_threads(2) shared(bestParallelDistance, parallelPath) private(newPath)
    {
        #pragma omp for
        for (int i = 1; i < numCities - 1; i++) {
            newPath = (int*)malloc(numCities * sizeof(int)); // Assign newPath inside the parallel region
            memcpy(newPath, parallelPath, numCities * sizeof(int));
            for (int j = i + 1; j < numCities; j++) {
                twoOptSwap(newPath, i, j);
                double newDistance = calculateTotalDistance(newPath, cities, numCities);
                #pragma omp critical
                {
                    if (newDistance < bestParallelDistance) {
                        bestParallelDistance = newDistance;
                        memcpy(parallelPath, newPath, numCities * sizeof(int));
                    }
                }
            }
            free(newPath);
        }
    }

    double parallelEnd = omp_get_wtime();
    double parallelTime = parallelEnd - parallelStart;
    printf("Parallel Minimum Distance: %.2f\n", bestParallelDistance);
    printf("Parallel Elapsed Time: %.6f seconds\n\n", parallelTime);

    // Print speedup
    double speedup = serialTime / parallelTime;
    printf("Speedup: %.2f\n", speedup);

    free(path);
    free(parallelPath);

    return 0;
}

