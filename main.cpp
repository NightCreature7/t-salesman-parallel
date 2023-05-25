#include <iostream>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
#include <omp.h>
#include <algorithm>


#define MAX_CITIES 100

struct City {
    double x;
    double y;
};

double calculateDistance(City city1, City city2) {
    double dx = city1.x - city2.x;
    double dy = city1.y - city2.y;
    return std::sqrt(dx * dx + dy * dy);
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

double calculateSerialBruteForce(int* path, City* cities, int numCities) {
    double bestDistance = calculateTotalDistance(path, cities, numCities);

    do {
        double newDistance = calculateTotalDistance(path, cities, numCities);
        if (newDistance < bestDistance) {
            bestDistance = newDistance;
        }
    } while (std::next_permutation(path, path + numCities));

    return bestDistance;
}

double calculateParallelBruteForce(int* path, City* cities, int numCities) {
    double bestDistance = calculateTotalDistance(path, cities, numCities);
    double newDistance;

#pragma omp parallel for shared(bestDistance) private(newDistance)
    for (int i = 0; i < numCities; i++) {
        int localPath[MAX_CITIES];
        std::memcpy(localPath, path, numCities * sizeof(int));
        std::next_permutation(localPath + i, localPath + numCities);

        do {
            newDistance = calculateTotalDistance(localPath, cities, numCities);
#pragma omp critical
            {
                if (newDistance < bestDistance) {
                    bestDistance = newDistance;
                }
            }
        } while (std::next_permutation(localPath, localPath + numCities));
    }

    return bestDistance;
}

int main() {
    int maxCities;
    std::cout << "Enter the maximum number of cities for the table: ";
    std::cin >> maxCities;

    if (maxCities <= 0 || maxCities > MAX_CITIES) {
        std::cout << "Invalid maximum number of cities. Exiting..." << std::endl;
        return 1;
    }

    int choice;
    std::cout << "Choose an option:\n";
    std::cout << "1. Calculate both Brute Force and Two Opt distances.\n";
    std::cout << "2. Calculate only Brute Force distance.\n";
    std::cout << "3. Calculate only Two Opt distance.\n";
    std::cout << "Enter your choice: ";
    std::cin >> choice;

    if (choice < 1 || choice > 3) {
        std::cout << "Invalid choice. Exiting..." << std::endl;
        return 1;
    }

    std::cout << "Cities\t\t\tBrute Force\t\t\tTwo-Opt" << std::endl;

    for (int numCities = 1; numCities <= maxCities; numCities++) {
        City cities[MAX_CITIES];
        std::srand(std::time(NULL));

        // Generate random cities
        for (int i = 0; i < numCities; i++) {
            cities[i].x = (double)std::rand() / RAND_MAX;
            cities[i].y = (double)std::rand() / RAND_MAX;
        }

        // Select a random starting city
        int startCity = std::rand() % numCities;

        int* path = new int[numCities];
        for (int i = 0; i < numCities; i++) {
            path[i] = i;
        }

        double bruteForceDistance = 0.0;
        double twoOptDistance = 0.0;

        if (choice == 1 || choice == 2) {
            bruteForceDistance = calculateSerialBruteForce(path, cities, numCities);
        }

        if (choice == 1 || choice == 3) {
            twoOptDistance = calculateTotalDistance(path, cities, numCities);
            for (int i = 1; i < numCities - 1; i++) {
                for (int j = i + 1; j < numCities; j++) {
                    twoOptSwap(path, i, j);
                    double newDistance = calculateTotalDistance(path, cities, numCities);
                    if (newDistance < twoOptDistance) {
                        twoOptDistance = newDistance;
                    }
                }
            }
        }

        std::cout << numCities << "\t\t\t" << bruteForceDistance << "\t\t\t" << twoOptDistance << std::endl;

        delete[] path;
    }

    return 0;
}
