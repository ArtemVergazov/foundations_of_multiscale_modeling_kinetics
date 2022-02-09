#include "simulation.h"
#include <cmath>

Simulation::Simulation(SimulationParams params) :
    params_(params),
    x_(params_.N_),
    y_(params_.N_),
    z_(params_.N_),
    vx_(params_.N_),
    vy_(params_.N_),
    vz_(params_.N_)
{
    // Initialize coordinates and velocities.
    
}

void Simulation::run() {
    step();
}

void Simulation::step() {

}

void Simulation::save() const {

}

double Simulation::U(double r) {
    return 4 * params_.eps_ * (pow(params_.sigma_ / r, 12) - pow(params_.sigma_ / r, 6));
}
