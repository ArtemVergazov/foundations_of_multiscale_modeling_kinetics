#include "linalg.h"
#include <stdexcept>
#include <string>
#include <random>

linalg::Vector::Vector(double x, double y, double z) : vector_() {
    vector_.push_back(x);
    vector_.push_back(y);
    vector_.push_back(z);
}

linalg::Vector &linalg::Vector::operator+=(const Vector &other) {
    if (vector_.size() != other.vector_.size()) {
        throw std::invalid_argument("Wrong size for vector addition!");
    }

    for (size_t i = 0; i < vector_.size(); ++i) {
        vector_[i] += other.vector_[i];
    }

    return *this;
}

linalg::Vector linalg::operator*(double val, Vector vec) {
    Vector res(vec);
    for (size_t i = 0; i < res.size(); ++i) {
        res[i] *= val;
    }
    return res;
}

linalg::Vector &linalg::Vector::operator-=(const Vector &other) {
    return (*this) += (-1) * other;
}

linalg::Vector linalg::operator+(Vector left, const Vector &right) {
    return left += right;
}

linalg::Vector linalg::operator-(Vector left, const Vector &right) {
    return left -= right;
}

linalg::Vector linalg::linspace(double start, double stop, size_t num) {
    if (num < 2) throw std::logic_error(
        "Cannot create linspace with " + std::to_string(num) + "elements!"
    );
    
    if (stop - start < 0) throw std::logic_error("Cannot create linspace with stop < start!");

    double delta = (stop - start) / (num - 1);
    linalg::Vector v(num);
    v[0] = start; // make sure start is exactly the same as the input

    for (size_t i = 1; i < num - 1; ++i) {
        v[i] = start + delta * (i - 1);
    }

    v[num - 1] = stop; // make sure start is exactly the same as the input
    return v;
}

linalg::Vector linalg::randomUniform(double low, double high, size_t size) {
    if (!size) return Vector();

    std::random_device rd; // TODO: reproducibility
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(low, high);

    Vector v(size);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = distribution(gen);
    }

    return v;
}

double linalg::norm(const Vector &v) {
    double res = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        res += v[i] * v[i];
    }
    return sqrt(res);
}
