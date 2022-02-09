#include "simulation.h"

constexpr int N = 64;
constexpr double boxSize = 16;
constexpr double dt = 1;
constexpr double eps = 1;
constexpr double sigma = 1;

int main() {
    Simulation({ N, boxSize, dt, eps, sigma }).run();
}
