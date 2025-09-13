#ifndef LIE_HPP_
#define LIE_HPP_

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <ostream>
#include <utility>
#include <vector>

namespace lie {

template<typename T, int R>
struct Vector;

template<typename T, int D>
struct SO;

template<typename T>
struct Array
{
    explicit Array(int size) { data_.resize(size); }
    Array<T> (const Array<T>& other);
    Array<T> (Array<T>&& other);
    Array<T>& operator=(const Array<T>& other);
    Array<T>& operator=(Array<T>&& other);

    inline T& at(int i) { return data_.at(i); }
    inline T at(int i) const { return data_.at(i); }
    inline int size() const { return data_.size(); }

    Array<T> operator*(const Array<T>& other) const;
    Array<T> operator-(const Array<T>& other) const;
    Array<T> operator+(const Array<T>& other) const;
private:
    std::vector<T> data_;
}; // struct Array

template<typename T, int R, int C>
struct Matrix
{
private:
    struct MatrixInitializer
    {
        explicit MatrixInitializer(Matrix<T, R, C> *m) : m(m), idx(1) {}
        MatrixInitializer& operator,(T v);
        Matrix<T, R, C> *m;
        int idx;
    }; // struct MatrixInitializer
public:
    MatrixInitializer operator<<(T v);
    Matrix() : array_(R * C), row_(R), col_(C) {}
    Matrix(const Array<T>& a) : array_(a), row_(R), col_(C) { assert(a.size() == R * C); }

    inline Array<T>& array() { return array_; }
    inline const Array<T>& array() const { return array_; }
    inline int row() const { return row_; }
    inline int col() const { return col_; }
    inline int size() const { return array_.size(); }
    inline T at(int r, int c) const { return array_.at(r * C + c); }
    inline T& at(int r, int c) { return array_.at(r * C + c); }

    virtual inline Matrix operator-(const Matrix& other) const { return Matrix<T, R, C>(array_ - other.array_); }
    virtual inline Matrix operator+(const Matrix& other) const { return Matrix<T, R, C>(array_ + other.array_); }
    template<int C2> Matrix<T, R, C2> operator*(const Matrix<T, C, C2>& other) const;

    Matrix<T, C, R> t() const;
    bool is_skew_sym() const;
    Vector<T, R> vee() const;
    SO<T, R> exp() const;
private:
    Array<T> array_;
    int row_;
    int col_;
}; // struct Matrix

template<int R, int C>
using Matrixf = Matrix<float, R, C>;
template<int R, int C>
using Matrixd = Matrix<double, R, C>;
using Matrix2d = Matrix<double, 2, 2>;
using Matrix3d = Matrix<double, 3, 3>;
using Matrix4d = Matrix<double, 4, 4>;
using Matrix2f = Matrix<float, 2, 2>;
using Matrix3f = Matrix<float, 3, 3>;
using Matrix4f = Matrix<float, 4, 4>;

template<typename T, int D>
struct Vector: public Matrix<T, D, 1>
{
    Vector() : Matrix<T, D, 1>() {}
    Matrix<T, D, D> hat() const;
    SO<T, D> Exp() const;
}; // struct Vector

template<int R>
using Vectorf = Vector<float, R>;
template<int R>
using Vectord = Vector<double, R>;
using Vector2f = Vector<float, 2>;
using Vector2d = Vector<double, 2>;
using Vector3f = Vector<float, 3>;
using Vector3d = Vector<double, 3>;
using Vector4f = Vector<float, 4>;
using Vector4d = Vector<double, 4>;

template<typename T, int D>
struct SO: public Matrix<T, D, D>
{
    SO() : Matrix<T, D, D>() {}
    Vector<T, D> operator*(const Vector<T, D>& v) const;
    Vector<T, D> Log() const;
    Matrix<T, D, D> log() const;
}; // struct SO

using SO2f = SO<float, 2>;
using SO2d = SO<double, 2>;
using SO3f = SO<float, 3>;
using SO3d = SO<double, 3>;


template<typename T>
Array<T>::Array(const Array<T>& other)
    : data_(other.data_)
{}

template<typename T>
Array<T>::Array(Array<T>&& other)
    : data_(std::move(other.data_))
{}

template<typename T>
Array<T>& Array<T>::operator=(const Array<T>& other)
{
    data_ = other.data_;
}

template<typename T>
Array<T>& Array<T>::operator=(Array<T>&& other)
{
    data_ = std::move(other.data_);
}

template<typename T>
Array<T> Array<T>::operator*(const Array<T>& other) const
{
    assert(size() == other.size());
    Array<T> out(*this);
    for (int i = 0; i < size(); ++i) {
        out.at(i) *= other.at(i);
    }

    return out;
}

template<typename T>
Array<T> Array<T>::operator-(const Array<T>& other) const
{
    assert(size() == other.size());
    Array<T> out(*this);
    for (int i = 0; i < size(); ++i) {
        out.at(i) -= other.at(i);
    }

    return out;
}

template<typename T>
Array<T> Array<T>::operator+(const Array<T>& other) const
{
    assert(size() == other.size());
    Array<T> out(*this);
    for (int i = 0; i < size(); ++i) {
        out.at(i) += other.at(i);
    }

    return out;
}

template<typename T, int R, int C>
typename Matrix<T, R, C>::MatrixInitializer Matrix<T, R, C>::operator<<(T v)
{
    array_.at(0) = v;
    return MatrixInitializer(this);    
}

template<typename T, int R, int C>
typename Matrix<T, R, C>::MatrixInitializer& Matrix<T, R, C>::MatrixInitializer::operator,(T v)
{
    m->array().at(idx++) = v;
    return *this;
}

template<typename T, int R, int C>
template<int C2>
Matrix<T, R, C2> Matrix<T, R, C>::operator*(const Matrix<T, C, C2>& other) const
{
    Matrix<T, R, C2> out;
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C2; ++j) {
            T val{0};
            for (int k = 0; k < C; ++k) {
                val += at(i, k) * other.at(k, j);
            }
            out.at(i, j) = val;
        }
    }
    return out;
}

template<typename T, int R, int C>
Matrix<T, C, R> Matrix<T, R, C>::t() const
{
    Matrix<T, C, R> out;
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            out.at(j, i) = at(i, j);
        }
    }

    return out;
}

template<typename T, int R, int C>
bool Matrix<T, R, C>::is_skew_sym() const
{
    if (R != C) return false;

    for (int i = 0; i < R; ++i) {
        for (int j = i; j < C; ++j) {
            if (i == j) if (at(i, j) != 0) {
                return false;
            }
            if (std::abs(at(i, j) + at(j, i)) > std::numeric_limits<T>::epsilon()) {
                return false;
            }
        }
    }

    return true;
}

template<typename T, int R, int C>
Vector<T, R> Matrix<T, R, C>::vee() const
{

}

template<typename T, int R, int C>
SO<T, R> Matrix<T, R, C>::exp() const
{

}

template<typename T, int D>
Matrix<T, D, D> Vector<T, D>::hat() const
{

}

template<typename T, int D>
SO<T, D> Vector<T, D>::Exp() const
{

}

template<typename T, int D>
Vector<T, D> SO<T, D>::Log() const
{

}

template<typename T, int D>
Matrix<T, D, D> SO<T, D>::log() const
{

}

template<typename T, int R, int C>
std::ostream& operator<<(std::ostream& os, const Matrix<T, R, C>& m)
{
    os << "Matrix" << typeid(T).name() << "<" << m.row() << ", " << m.col() << ">{";
    os << std::fixed << std::setprecision(3);
    for (int i = 0; i < m.size(); ++i) {
        if (i % m.col() == 0) os << "\n  ";
        os << m.array().at(i) << ", ";
    }
    os << "\n}";
    return os;
}

template<typename T, int D>
std::ostream& operator<<(std::ostream& os, const SO<T, D>& so)
{
    os << "SO" << D << typeid(T).name() << " {";
    os << std::fixed << std::setprecision(3);
    for (int i = 0; i < so.size(); ++i) {
        if (i % so.col() == 0) os << "\n  ";
        os << so.array().at(i) << ", ";
    }
    os << "\n}";
    return os;
}

} // namespace lie

#endif // LIE_HPP_


#define LIE_IMPLEMENTATION // delete me


#ifdef LIE_IMPLEMENTATION
#ifndef LIE_CPP_
#define LIE_CPP_

namespace lie {

}; // namespace lie

#endif // LIE_CPP_
#endif // LIE_IMPLEMENTATION