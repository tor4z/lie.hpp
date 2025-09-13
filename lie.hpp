#ifndef LIE_HPP_
#define LIE_HPP_

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <utility>
#include <vector>

namespace lie {

template<typename T, int D>
struct SO;

template<typename T, int R, int C>
struct Matrix;

template<typename T, int D>
using Vector = Matrix<T, D, 1>;

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
    inline bool operator!=(const Array<T>& other) const { return !(*this == other); }

    bool operator==(const Array<T>& other) const;
    Array<T> operator*(T v) const;
    Array<T> operator-(T v) const;
    Array<T> operator+(T v) const;
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
    inline T at(int i) const { return array_.at(i); }
    inline T& at(int i) { return array_.at(i); }
    inline T at(int r, int c) const { return array_.at(r * C + c); }
    inline T& at(int r, int c) { return array_.at(r * C + c); }

    virtual bool operator!=(const Matrix& other) const { return array_ != other.array_; }
    virtual bool operator==(const Matrix& other) const { return array_ == other.array_; }
    virtual inline Matrix operator-(const Matrix& other) const { return Matrix<T, R, C>(array_ - other.array_); }
    virtual inline Matrix operator+(const Matrix& other) const { return Matrix<T, R, C>(array_ + other.array_); }
    virtual inline Matrix operator*(T v) const { return Matrix<T, R, C>(array_ * v); }
    template<int C2> Matrix<T, R, C2> operator*(const Matrix<T, C, C2>& other) const;

    bool is_skew_sym() const;
    T tr() const;
    T norm2() const;
    T det() const;
    SO<T, R> Exp() const;
    SO<T, R> exp() const;
    Matrix<T, R, R> hat() const;
    Matrix<T, C, R> t() const;
    Vector<T, R> vee() const;
    Vector<T, R> null_space() const;

    static Matrix<T, R, C> eye();
private:
    Array<T> array_;
    int row_;
    int col_;
}; // struct Matrix

// invalid matrix
template<typename T>
struct Matrix<T, 0, 0>
{
    T at(int i) const;
    T& at(int i);
    T at(int r, int c) const;
    T& at(int r, int c);
    T det() const;
}; // struct Matrix<T, 0, 0>

template<typename T, int D>
struct SO: public Matrix<T, D, D>
{
    using Matrix<T, D, D>::array;
    using Matrix<T, D, D>::tr;

    SO() : Matrix<T, D, D>() {}
    explicit SO(const Array<T>& a) : Matrix<T, D, D>(a) {}
    inline SO<T, D> t() const { return SO<T, D>(mat().t().array()); }
    inline Matrix<T, D, D> mat() const { return Matrix<T, D, D>(array()); }
    /**
     * @brief Composition
     * 
     * @param other 
     * @return SO<T, D> 
     */
    inline SO<T, D> operator*(const SO<T, D>& other) const { return SO<T, D>((mat() * other.mat()).array()); }

    SO r_plus(const Vector<T, D>& v) const;
    Vector<T, D> r_minus(const SO so) const;
    Vector<T, D> Log() const;
    Matrix<T, D, D> log() const;
}; // struct SO

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

template<int D>
using Vectorf = Vector<float, D>;
template<int D>
using Vectord = Vector<double, D>;
using Vector2f = Vectorf<2>;
using Vector2d = Vectord<2>;
using Vector3f = Vectorf<3>;
using Vector3d = Vectord<3>;
using Vector4f = Vectorf<4>;
using Vector4d = Vectord<4>;

template<int D>
using SOf = SO<float, D>;
template<int D>
using SOd = SO<double, D>;
using SO2f = SOf<2>;
using SO2d = SOd<2>;
using SO3f = SOf<3>;
using SO3d = SOd<3>;

template<typename T>
bool is_0(T v)
{
    return std::abs(v) <= std::numeric_limits<T>::epsilon();
}

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
bool Array<T>::operator==(const Array<T>& other) const
{
    if (size() != other.size()) return false;
    for (int i = 0; i < size(); ++i) {
        if (!is_0(at(i) - other.at(i))) {
            return false;
        }
    }

    return true;
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

template<typename T>
Array<T> Array<T>::operator*(T v) const
{
    Array<T> out(*this);
    for (int i = 0; i < size(); ++i) {
        out.at(i) *= v;
    }

    return out;
}

template<typename T>
Array<T> Array<T>::operator-(T v) const
{
    Array<T> out(*this);
    for (int i = 0; i < size(); ++i) {
        out.at(i) -= v;
    }

    return out;
}

template<typename T>
Array<T> Array<T>::operator+(T v) const
{
    Array<T> out(*this);
    for (int i = 0; i < size(); ++i) {
        out.at(i) += v;
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
        if (!is_0(at(i, i))) {
            return false;
        }
        for (int j = i; j < C; ++j) {
            if (!is_0(at(i, j) + at(j, i))) {
                return false;
            }
        }
    }

    return true;
}

template<typename T, int R, int C>
Vector<T, R> Matrix<T, R, C>::vee() const
{
    static_assert(R == C);
    static_assert(R == 3);

    Vector<T, R> out;
    out.at(2) = -at(0, 1);
    out.at(1) = at(0, 2);
    out.at(0) = -at(1, 2);
    return out;
}

template<typename T, int R, int C>
Vector<T, R> Matrix<T, R, C>::null_space() const
{
    // gaussian elimination
    Matrix tmp{*this};

    auto swap_row{[] (Matrix<T, R, C>& m, int r1, int r2) -> void {
        for (int c = 0; c < C; ++c) {
            std::swap(m.at(r1, c), m.at(r2, c));
        }
    }};

    auto has_not0_before_c{[&tmp] (int r, int c) -> bool {
        for (int i = 0; i < c; ++i) {
            if (!is_0(tmp.at(r, i))) {
                return true;
            }
        }
        return false;
    }};

    auto has_not0_after_c{[&tmp] (int r, int c) -> bool {
        for (int i = c + 1; i < C; ++i) {
            if (!is_0(tmp.at(r, i))) {
                return true;
            }
        }
        return false;
    }};

    T pivot{};
    // downword elimination
    for (int c = 0; c < C; ++c) {
        pivot = tmp.at(c, c);
        if (is_0(pivot)) {
            // try swap row
            for (int r = 0; r < R; ++r) {
                if (c != r && !is_0(tmp.at(r, c))) {
                    if (r < c && has_not0_before_c(r, c)) { continue; }
                    swap_row(tmp, r, c);
                    break;
                }
            }
            pivot = tmp.at(c, c);   // update pivot
        }

        if (is_0(pivot)) break; // finished
        for (int r = 0; r < R; ++r) {
            if (is_0(tmp.at(r, c)) || r == c) {
                continue;
            }

            const auto ceof{tmp.at(r, c) / pivot};
            for (int c1 = 0; c1 < C; ++c1) {
                tmp.at(r, c1) -= tmp.at(c, c1) * ceof;
            }
        }
    }

    // upword elimination and solve
    Vector<T, R> out;
    for (int r = 0; r < R; ++r) {
        if (!is_0(tmp.at(r, r))) {
            if (!has_not0_after_c(r, C)) {
                for (int r1 = 0; r1 < r; ++r1) {
                    tmp.at(r1, r) = 0;
                }
                out.at(r) = 0;
            } else {
                out.at(r) = 1;
                out.at(r) = 1;
            }
        }
    }

    std::cout << "tmp: " << tmp << "\n";
    return out;
}

template<typename T, int R, int C>
SO<T, R> Matrix<T, R, C>::exp() const
{
    static_assert(R == C);

    auto vee_v{vee()};
    T vee_norm2{vee_v.norm2()};
    auto hat_m{(vee_v * (1 / vee_norm2)).hat()};
    auto m{eye() + (hat_m * hat_m) * (1 - std::cos(vee_norm2)) + hat_m * std::sin(vee_norm2)};
    return SO<T, R>{m.array()};
}

template<typename T, int R, int C>
Matrix<T, R, C> Matrix<T, R, C>::eye()
{
    static_assert(R == C);
    Matrix<T, R, C> out;
    for (int i = 0; i < R; ++i) {
        out.at(i, i) = 1;
    }
    return out;
}

template<typename T, int R, int C>
Matrix<T, R, R> Matrix<T, R, C>::hat() const
{
    static_assert(R == 3);
    static_assert(C == 1);
    Matrix<T, R, R> out;
    out.at(0, 1) = -at(2);
    out.at(0, 2) = at(1);
    out.at(1, 2) = -at(0);
    out.at(1, 0) = at(2);
    out.at(2, 0) = -at(1);
    out.at(2, 1) = at(0);
    return out;
}

template<typename T, int R, int C>
SO<T, R> Matrix<T, R, C>::Exp() const
{
    static_assert(C == 1);
    return hat().exp();
}

template<typename T, int R, int C>
T Matrix<T, R, C>::det() const
{
    static_assert(R == C);
    if (R == 1) {
        return array_.at(0);
    } else if (R == 2) {
        return at(0, 0) * at(1, 1) - at(1, 0) * at(0, 1);
    }

    // calculate determinant with cofactor expansions
    T result{0};
    Matrix<T, R - 1, C - 1> minor;
    for (int c = 0; c < C; ++c) {
        T sign{(1 + 1 + c) % 2 == 0 ? static_cast<T>(1) : static_cast<T>(-1)};
        T pivot{at(0, c)};

        // expansion for first row
        for (int mr = 1; mr < R; ++mr) {
            int cc{0};
            for (int mc = 0; mc < C; ++mc) {
                if (c == mc) continue;
                minor.at(mr - 1, cc++) = at(mr, mc);
            }
        }

        result += sign * pivot * minor.det();
    }
    return result;
}

template<typename T, int R, int C>
T Matrix<T, R, C>::norm2() const
{
    static_assert(C == 1, "Undefined");

    if (C == 1) {
        T sq_sum{0};
        for (int i = 0; i < R; ++i) {
            sq_sum += at(i) * at(i);
        }
        return std::sqrt(sq_sum);
    }
}

template<typename T, int R, int C>
T Matrix<T, R, C>::tr() const
{
    static_assert(R == C);
    T result{0};
    for (int i = 0; i < R; ++i) {
        result += at(i, i);
    }
    return result;
}

template<typename T, int D>
Vector<T, D> SO<T, D>::Log() const
{
    auto theta{std::acos((tr() - 1) / 2)};
    Vector<T, D> norm_v{(*this - Matrix<T, D, D>::eye()).null_space()};   
    return norm_v * theta;
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

template<typename T, int D>
std::ostream& operator<<(std::ostream& os, const Vector<T, D>& v)
{
    os << "Vector" << D << typeid(T).name() << " {";
    os << std::fixed << std::setprecision(3);
    for (int i = 0; i < v.size(); ++i) {
        os << "\n  " << v.array().at(i) << ", ";
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