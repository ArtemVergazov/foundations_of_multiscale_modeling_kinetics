#pragma once
#include <vector>

// Imitating useful numpy-like functions.
namespace linalg {

    // np.array
    class Vector {
    public:
        explicit Vector(size_t count = 0, double val = 0) : vector_(count, val) {}

        // 3D vector - specific case.
        Vector(double x, double y, double z);

        size_t size() const { return vector_.size(); }
        void push_back(double val) { vector_.push_back(val); }

        double operator[](size_t i) const { return vector_[i]; }
        double &operator[](size_t i) { return vector_[i]; }

        Vector &operator+=(const Vector &other);
        Vector &operator-=(const Vector &other);

    private:
        std::vector<double> vector_;
    };

    Vector operator*(double val, Vector vec);
    Vector operator+(Vector left, const Vector &right);
    Vector operator-(Vector left, const Vector &right);

    // np.linspace
    Vector linspace(double start, double stop, size_t num = 50);

    // np.random.uniform
    Vector randomUniform(double low = 0, double high = 1, size_t size = 1);

    // np.linalg.norm
    double norm(const Vector &v);
}
