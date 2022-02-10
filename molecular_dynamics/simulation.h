#pragma once
#include "linalg.h"
#include <vector>

class Simulation {
public:
    Simulation(const Simulation &) = delete;
    Simulation &operator=(const Simulation &) = delete;

    Simulation();
    void run();

    // N-sized vectors of coordinates.
private:
    linalg::Vector x_;
    linalg::Vector y_;
    linalg::Vector z_;

    // N-sized vectors of velocities.
private:
    linalg::Vector vx_;
    linalg::Vector vy_;
    linalg::Vector vz_;

    // N-sized vectors of accelerations.
private:
    linalg::Vector ax_;
    linalg::Vector ay_;
    linalg::Vector az_;

//    // History of energies of the system.
//private:
//    linalg::Vector E_;

private:
    void step();
    void save() const;

    static double U(const linalg::Vector &r1, const linalg::Vector &r2);
    static double U(double r);

    //double getEnergy() const;

    // Calculate accelerations from forces based on current coordinates.
    void getAccelerations(
        linalg::Vector &ax,
        linalg::Vector &ay,
        linalg::Vector &az
    ) const;
};
