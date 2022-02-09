#pragma once
#include <vector>

class Simulation {
public:
    struct SimulationParams {
        int N_;
        double boxSize_;
        double dt_;
        double eps_;
        double sigma_;
    } params_;

    Simulation() = delete;
    Simulation(const Simulation &) = delete;
    Simulation &operator=(const Simulation &) = delete;

    Simulation(SimulationParams params);
    void run();

    // N-sized vectors of coordinates.
private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;

    // N-sized vectors of velocities.
private:
    std::vector<double> vx_;
    std::vector<double> vy_;
    std::vector<double> vz_;

private:
    void step();
    void save() const;
    double U(double r);
};
