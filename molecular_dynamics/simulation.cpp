#include "simulation.h"
#include "linalg.h"
#include <cmath>
#include <algorithm>

// Always make sure N is a cube of a natural number.
// It is important for initial placement of the molecules.
constexpr size_t N = 64;

constexpr double boxSize = 16;
constexpr double dt = 1;
constexpr double eps = 1;
constexpr double sigma = 1;

Simulation::Simulation() :
    x_(N), y_(N), z_(N),
    vx_(N), vy_(N), vz_(N),
    ax_(N), ay_(N), az_(N)
{
    size_t N1D = std::lround(std::cbrt(N));
    double cellSize = boxSize / N1D;

    // Initializing coordinates.
    auto x1D = linalg::linspace(cellSize / 2, boxSize - cellSize / 2, N1D);
    auto y1D = linalg::linspace(cellSize / 2, boxSize - cellSize / 2, N1D);
    auto z1D = linalg::linspace(cellSize / 2, boxSize - cellSize / 2, N1D);

    size_t count = 0;
    for (size_t i = 0; i < N1D; ++i) {
        double x = x1D[i];
        for (size_t j = 0; j < N1D; ++j) {
            double y = y1D[j];
            for (size_t k = 0; k < N1D; ++k) {
                double z = z1D[k];
                x_[count] = x;
                y_[count] = y;
                z_[count++] = z;
            }
        }
    }
    
    // Generating velocities.
    vx_ = linalg::randomUniform(-boxSize, boxSize, N);
    vy_ = linalg::randomUniform(-boxSize, boxSize, N);
    vz_ = linalg::randomUniform(-boxSize, boxSize, N);

    // Center of masses velocity.

    // Calculate accelerations.


    //// Calculate initial energy.
    //E_.push_back(getEnergy());
}

void Simulation::run() {
    step();
}

void Simulation::step() {
    // Calculate current accelerations.

    // Calculate next step coordinates.

    // Calculate next step acceleration.
    
    // Calculate next step velocities.
}

void Simulation::save() const {

}

double Simulation::U(const linalg::Vector &r1, const linalg::Vector &r2) {
    linalg::Vector dr = r1 - r2;
    double dx = std::abs(dr[0]);
    double dy = std::abs(dr[1]);
    double dz = std::abs(dr[2]);

    linalg::Vector r(
        std::min(dx, boxSize - dx),
        std::min(dy, boxSize - dy),
        std::min(dz, boxSize - dz)
    );

    return U(linalg::norm(r));
}

double Simulation::U(double r) {
    double sigma_over_r6 = pow(sigma / r, 6);
    return 4 * eps * sigma_over_r6 * (sigma_over_r6 - 1);
}

//double Simulation::getEnergy() const {
//    double E = 0;
//
//    // Kinetic energy.
//    for (size_t i = 0; i < N; ++i) {
//        E += (vx_[i] * vx_[i] + vy_[i] * vy_[i] + vz_[i]) / 2;
//    }
//
//    // Potential energy.
//    for (size_t i = 0; i < N; ++i) {
//        linalg::Vector ri(x_[i], y_[i], z_[i]);
//        for (size_t j = i + 1; j < N; ++j) {
//            linalg::Vector rj(x_[j], y_[j], z_[j]);
//            double r = linalg::norm(rj - ri);
//        }
//    }
//
//    return E;
//}

void Simulation::getAccelerations(
    linalg::Vector &ax,
    linalg::Vector &ay,
    linalg::Vector &az
) const {
    for (size_t i = 0; i < N; ++i) {
        linalg::Vector ri(x_[i], y_[i], z_[i]);
        ax[i] = 0; ay[i] = 0; az[i] = 0;

        for (size_t j = i + 1; j < N; ++j) {
            linalg::Vector rj(x_[j], y_[j], z_[j]);

        }
    }
}
