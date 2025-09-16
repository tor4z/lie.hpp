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

template<typename T>
struct SOX;

template<typename T>
struct SEX;

template<typename T>
struct MatrixX;

template<typename T>
struct VectorX;

template<typename T>
struct Array
{
    explicit Array(int size) { data_.resize(size, 0); }
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

template<typename T>
struct MatrixX
{
private:
    struct MatrixXInitializer
    {
        explicit MatrixXInitializer(MatrixX<T> *m) : m(m), idx(1) {}
        MatrixXInitializer& operator,(T v);
        MatrixX<T> *m;
        int idx;
    }; // struct MatrixXInitializer
public:
    MatrixXInitializer operator<<(T v);
    MatrixX(int row, int col) : array_(row * col), row_(row), col_(col) {}
    MatrixX(int row, int col, const Array<T>& a) : array_(a), row_(row), col_(col) { assert(a.size() == row * col); }

    inline Array<T>& array() { return array_; }
    inline const Array<T>& array() const { return array_; }
    inline int row() const { return row_; }
    inline int col() const { return col_; }
    inline int size() const { return array_.size(); }
    inline T at(int i) const { return array_.at(i); }
    inline T& at(int i) { return array_.at(i); }
    inline T at(int r, int c) const { return array_.at(r * col_ + c); }
    inline T& at(int r, int c) { return array_.at(r * col_ + c); }

    virtual bool operator!=(const MatrixX& other) const { return array_ != other.array_; }
    virtual bool operator==(const MatrixX& other) const { return array_ == other.array_; }
    virtual inline MatrixX operator-(const MatrixX& other) const { return MatrixX(row_, col_, array_ - other.array_); }
    virtual inline MatrixX operator+(const MatrixX& other) const { return MatrixX(row_, col_, array_ + other.array_); }
    inline MatrixX operator*(T v) const { return MatrixX(row_, col_, array_ * v); }
    MatrixX operator*(const MatrixX<T>& other) const;
    MatrixX sub_mat(int sr, int sc, int r, int c) const;
    void sub_mat_assign(int sr, int sc, int r, int c, const MatrixX& other);

    bool is_skew_sym() const;
    T tr() const;
    T norm2() const;
    T det() const;
    SOX<T> SO_Exp() const;
    SOX<T> SO_exp() const;
    SEX<T> SE_Exp() const;
    SEX<T> SE_exp() const;
    MatrixX<T> hat() const;
    MatrixX<T> t() const;
    VectorX<T> vee() const;
    VectorX<T> null_space() const;
    MatrixX normalize() const;
    MatrixX inv() const;

    static MatrixX eye(int row, int col);
private:
    Array<T> array_;
    const int row_;
    const int col_;
}; // struct MatrixX

template<typename T>
struct VectorX: public MatrixX<T>
{
    using MatrixX<T>::array;
    using MatrixX<T>::at;
    using MatrixX<T>::row;

    explicit VectorX(int dim) : MatrixX<T>(dim, 1) {}
    VectorX(const MatrixX<T>& m) : MatrixX<T>(m) { assert(m.col() == 1); }
    VectorX(int dim, const Array<T>& a) : MatrixX<T>(dim, 1, a) {}
    inline MatrixX<T> mat() const { return MatrixX<T>(dim(), 1, array()); }
    inline int dim() const { return row(); }
    inline VectorX operator*(T v) const { return VectorX(dim(), array() * v); }
}; // struct VectorX

template<typename T>
struct SOX: public MatrixX<T>
{
    using MatrixX<T>::array;
    using MatrixX<T>::tr;
    using MatrixX<T>::at;
    using MatrixX<T>::row;
    using MatrixX<T>::col;

    inline static SOX eye(int dim) { return SOX(MatrixX<T>(dim).eye().array()); }

    explicit SOX(int dim) : MatrixX<T>(dim, dim) {}
    SOX(const MatrixX<T>& m) : MatrixX<T>(m) { assert(m.row() == m.col()); }
    SOX(int dim, const Array<T>& a) : MatrixX<T>(dim, dim, a) {}
    inline MatrixX<T> mat() const { return MatrixX<T>(dim(), dim(), array()); }
    inline int dim() const { return row(); }

    inline MatrixX<T> Adj() { return mat(); }
    inline VectorX<T> adj(const VectorX<T>& v) { return (mat() * v.hat() * mat().inv()).vee(); }
    inline SOX t() const { return SOX(dim(), MatrixX<T>::t().array()); }

    /**
     * @brief Composition
     * 
     * @param other 
     * @return SOX<T>
     */
    inline SOX operator*(const SOX<T>& other) const { return SOX(dim(), (this->mat() * other.mat()).array()); }

    /**
     * @brief action on a vector
     * 
     * @param other 
     * @return Vector<T>
     */
    inline VectorX<T> operator*(const VectorX<T>& other) const { return this->mat() * other; }

    /**
     * @brief inv
     * 
     * @param other 
     * @return VectorX<T> 
     */
    inline SOX inv() const { return SOX<T>(dim(), MatrixX<T>::inv().array()); }

    /**
     * @brief right circle plus
     * 
     * @param v 
     * @return SO 
     */
    SOX plus(const VectorX<T>& v) const;

    /**
     * @brief right circle minus
     * 
     * @param so 
     * @return Vector<T, D> 
     */
    VectorX<T> minus(const SOX& so) const;
    VectorX<T> Log() const;
    MatrixX<T> log() const;
}; // struct SOX

template<typename T>
struct SEX: public MatrixX<T>
{
    using MatrixX<T>::array;
    using MatrixX<T>::tr;
    using MatrixX<T>::at;
    using MatrixX<T>::col;
    using MatrixX<T>::row;
    using MatrixX<T>::sub_mat;
    using MatrixX<T>::sub_mat_assign;

    inline static int to_m_dim(int dim) { return dim + 1; }
    inline static SEX eye(int dim) { return SEX(MatrixX<T>::eye(to_m_dim(dim), to_m_dim(dim))); }
    
    explicit SEX(int dim) : MatrixX<T>(to_m_dim(dim), to_m_dim(dim)) {}
    SEX(int dim, const Array<T>& a) : MatrixX<T>(to_m_dim(dim), to_m_dim(dim), a) {}
    SEX(const MatrixX<T>& m) : MatrixX<T>(m) { assert(m.row() == m.col()); }
    SEX(const SOX<T>& so, const VectorX<T>& off);
    inline int dim() const { return row() - 1; }
    inline MatrixX<T> mat() const { return MatrixX<T>(to_m_dim(dim()), to_m_dim(dim()), array()); }
    inline VectorX<T> adj(const VectorX<T>& v) { return (mat() * v.hat() * mat().inv()).vee(); }
    
    MatrixX<T> Adj();
    VectorX<T> offset() const;
    SOX<T> rot() const;

    /**
     * @brief Composition
     * 
     * @param other 
     * @return SEX<T>
     */
    SEX<T> operator*(const SEX<T>& other) const;

    /**
     * @brief inverse
     * 
     * @return SEX<T> 
     */
    SEX<T> inv() const;

    /**
     * @brief action on a vector
     * 
     * @param other 
     * @return VectorX<T>
     */
    VectorX<T> operator*(const VectorX<T>& other) const;

    /**
     * @brief right circle plus
     * 
     * @param v 
     * @return SEX<T>
     */
    SEX plus(const VectorX<T>& v) const;

    /**
     * @brief right circle minus
     * 
     * @param se 
     * @return VectorX<T> 
     */
    VectorX<T> minus(const SEX& se) const;
    VectorX<T> Log() const;
    MatrixX<T> log() const;
}; // struct SEX

template<typename T, int R, int C>
struct Matrix: public MatrixX<T>
{
    Matrix() : MatrixX<T>(R, C) {}
    explicit Matrix(const Array<T>& a) : MatrixX<T>(R, C, a) {}
    static Matrix eye();
}; // struct Matrix

template<typename T, int D>
struct Vector: public VectorX<T>
{
    Vector() : VectorX<T>(D) {}
    Vector(const VectorX<T>& v) : VectorX<T>(v) {}
    explicit Vector(const Array<T>& a) : MatrixX<T>(D, 1, a) {}
}; // struct Vector

template<typename T, int D>
struct SO: public SOX<T>
{
    SO() : SOX<T>(D) {}
    explicit SO(const Array<T>& a) : SOX<T>(D, a) {}
    static SO eye();
}; // struct SO

template<typename T, int D>
struct SE: public SEX<T>
{
    SE() : SEX<T>(D) {}
    explicit SE(const Array<T>& a) : SEX<T>(D, a) {}
    static SE eye();
}; // struct SE

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
using Vector6f = Vectorf<6>;
using Vector6d = Vectord<6>;

template<int D>
using SOf = SO<float, D>;
template<int D>
using SOd = SO<double, D>;
using SO2f = SOf<2>;
using SO2d = SOd<2>;
using SO3f = SOf<3>;
using SO3d = SOd<3>;

template<int D>
using SEf = SE<float, D>;
template<int D>
using SEd = SE<double, D>;
using SE2f = SEf<2>;
using SE2d = SEd<2>;
using SE3f = SEf<3>;
using SE3d = SEd<3>;

template<typename T>
bool is_0(T v)
{
    return std::abs(v) <= std::numeric_limits<T>::epsilon();
}

template<typename T>
void swap_row(MatrixX<T>& m, int r1, int r2)
{
    for (int c = 0; c < m.col(); ++c) {
        std::swap(m.at(r1, c), m.at(r2, c));
    }
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
    return *this;
}

template<typename T>
Array<T>& Array<T>::operator=(Array<T>&& other)
{
    data_ = std::move(other.data_);
    return *this;
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

template<typename T>
typename MatrixX<T>::MatrixXInitializer MatrixX<T>::operator<<(T v)
{
    array_.at(0) = v;
    return MatrixXInitializer(this);    
}

template<typename T>
typename MatrixX<T>::MatrixXInitializer& MatrixX<T>::MatrixXInitializer::operator,(T v)
{
    m->array().at(idx++) = v;
    return *this;
}

template<typename T>
MatrixX<T> MatrixX<T>::operator*(const MatrixX<T>& other) const
{
    MatrixX<T> out(row_, other.col_);
    for (int i = 0; i < out.row_; ++i) {
        for (int j = 0; j < out.col_; ++j) {
            T val{0};
            for (int k = 0; k < col_; ++k) {
                val += at(i, k) * other.at(k, j);
            }
            out.at(i, j) = val;
        }
    }
    return out;
}

template<typename T>
MatrixX<T> MatrixX<T>::eye(int row, int col)
{
    assert(row == col);
    MatrixX<T> out(row, col);
    for (int i = 0; i < row; ++i) {
        out.at(i, i) = 1;
    }
    return out;
}

template<typename T>
MatrixX<T> MatrixX<T>::t() const
{
    MatrixX<T> out(col_, row_);
    for (int i = 0; i < row_; ++i) {
        for (int j = 0; j < col_; ++j) {
            out.at(j, i) = at(i, j);
        }
    }

    return out;
}

template<typename T>
bool MatrixX<T>::is_skew_sym() const
{
    if (row_ != col_) return false;

    for (int i = 0; i < row_; ++i) {
        if (!is_0(at(i, i))) {
            return false;
        }
        for (int j = i; j < col_; ++j) {
            if (!is_0(at(i, j) + at(j, i))) {
                return false;
            }
        }
    }

    return true;
}

template<typename T>
VectorX<T> MatrixX<T>::vee() const
{
    assert(col_ == col_);
    assert(row_ == 3 || row_ == 4);

    if (row_ == 3) {
        VectorX<T> out(3);
        out.at(2) = -at(0, 1);
        out.at(1) = at(0, 2);
        out.at(0) = -at(1, 2);
        return out;
    } else if (row_ == 4) {
        VectorX<T> out(6);
        out.at(5) = -at(0, 1);
        out.at(4) = at(0, 2);
        out.at(3) = -at(1, 2);

        out.at(0) = at(0, 3);
        out.at(1) = at(1, 3);
        out.at(2) = at(2, 3);
        return out;
    }

    return VectorX<T>(3); // unreachable, return something to supress compiler warning
}

template<typename T>
MatrixX<T> MatrixX<T>::hat() const
{
    assert(row_ == 3 || row_ == 6);
    assert(col_ == 1);

    if (row_ == 3) {
        MatrixX<T> out(3, 3);
        out.at(0, 1) = -at(2);
        out.at(0, 2) = at(1);
        out.at(1, 2) = -at(0);
        out.at(1, 0) = at(2);
        out.at(2, 0) = -at(1);
        out.at(2, 1) = at(0);
        return out;
    } else if (row_ == 6) {
        MatrixX<T> out(4, 4);
        out.at(0, 1) = -at(5);
        out.at(0, 2) = at(4);
        out.at(1, 2) = -at(3);
        out.at(1, 0) = at(5);
        out.at(2, 0) = -at(4);
        out.at(2, 1) = at(3);

        out.at(0, 3) = at(0);
        out.at(1, 3) = at(1);
        out.at(2, 3) = at(2);
        return out;
    }
    return MatrixX<T>(3, 3); // unreachable, return something to supress compiler warning
}

template<typename T>
VectorX<T> MatrixX<T>::null_space() const
{
    // gaussian elimination
    MatrixX tmp{*this};

    auto has_not0_before_c{[&tmp] (int r, int c) -> bool {
        for (int i = 0; i < c; ++i) {
            if (!is_0(tmp.at(r, i))) {
                return true;
            }
        }
        return false;
    }};

    auto has_not0_after_c{[&tmp, this] (int r, int c) -> bool {
        for (int i = c + 1; i < col_; ++i) {
            if (!is_0(tmp.at(r, i))) {
                return true;
            }
        }
        return false;
    }};

    T pivot{};
    // downword elimination
    for (int c = 0; c < col_; ++c) {
        pivot = tmp.at(c, c);
        if (is_0(pivot)) {
            // try swap row
            for (int r = 0; r < row_; ++r) {
                if (c != r && !is_0(tmp.at(r, c))) {
                    if (r < c && has_not0_before_c(r, c)) { continue; }
                    swap_row(tmp, r, c);
                    break;
                }
            }
            pivot = tmp.at(c, c);       // update pivot
        }

        if (is_0(pivot)) continue;      // finished
        for (int r = 0; r < row_; ++r) {
            if (is_0(tmp.at(r, c)) || r == c) {
                continue;
            }

            const auto ceof{tmp.at(r, c) / pivot};
            for (int c1 = 0; c1 < col_; ++c1) {
                tmp.at(r, c1) -= tmp.at(c, c1) * ceof;
            }
        }
    }

    // upword elimination
    for (int r = 0; r < row_; ++r) {
        if (!is_0(tmp.at(r, r))) {
            if (!has_not0_after_c(r, col_)) {
                for (int r1 = 0; r1 < r; ++r1) {
                    tmp.at(r1, r) = 0;
                }
            }
        }
    }

    // solve
    VectorX<T> out(col_);
    int not0_cs[2];  // not 0 col index on row r
    for (int r = 0; r < row_; ++r) {
        int num_not0{0};
        for (int c = 0; c < col_; ++c) {
            if (!is_0(tmp.at(r, c))) {
                not0_cs[num_not0++] = c;
            }
        }
        if (num_not0 == 1) {
            out.at(not0_cs[0]) = 0;
        } else if (num_not0 == 2) {
            if (is_0(out.at(not0_cs[1]))) {
                out.at(not0_cs[0]) = 1;
                out.at(not0_cs[1]) = -tmp.at(r, not0_cs[0]) / tmp.at(r, not0_cs[1]);
            } else {
                out.at(not0_cs[0]) = -out.at(not0_cs[1]) * tmp.at(r, not0_cs[1]) / tmp.at(r, not0_cs[0]);
            }
        } else if (num_not0 == 0) {
            if (is_0(out.at(r)))
                out.at(r) = 1;
        } else {
            assert(false && "Undefined");
        }
    }

    return out;
}

template<typename T>
T MatrixX<T>::norm2() const
{
    assert(col_ == 1 && "Undefined");

    T sq_sum{0};
    if (col_ == 1) {
        for (int i = 0; i < row_; ++i) {
            sq_sum += at(i) * at(i);
        }
    }
    return std::sqrt(sq_sum);
}

template<typename T>
MatrixX<T> MatrixX<T>::normalize() const
{
    assert(col_ == 1 && "Undefined");

    auto norm{norm2()};
    MatrixX<T> out(row_, col_);
    if (is_0(norm)) {
        return out;
    }

    if (col_ == 1) {
        for (int i = 0; i < size(); ++i) {
            out.at(i) = at(i) / norm;
        }
    }
    return out;
}

template<typename T>
SOX<T> MatrixX<T>::SO_Exp() const
{
    assert(col_ == 1);

    T theta{norm2()};
    assert(!is_0(theta) && "Invalid Exp operator");

    auto hat_m{(*this * (1 / theta)).hat()};
    auto eye_m{MatrixX<T>::eye(row_, row_)};
    auto m{eye_m + (hat_m * hat_m) * (1 - std::cos(theta)) + hat_m * std::sin(theta)};
    return SOX<T>(row_, m.array());
}

template<typename T>
SOX<T> MatrixX<T>::SO_exp() const
{
    assert(row_ == col_);
    return vee().SO_Exp();
}

template<typename T>
SEX<T> MatrixX<T>::SE_Exp() const
{
    assert(row_ % 2 == 0 && col_ == 1 && "Bad vector shape");

    VectorX<T> offset_vec{sub_mat(0, 0, row_ / 2, 1)};
    VectorX<T> rot_vec{sub_mat(row_ / 2, 0, row_ / 2, 1)};

    T theta{rot_vec.norm2()};
    assert(!is_0(theta) && "Invalid Exp operator");

    SOX<T> so{rot_vec.SO_Exp()};
    auto eye_m{MatrixX<T>::eye(row_ / 2, row_ / 2)};
    auto hat_m{(rot_vec * (1 / theta)).hat()};

    MatrixX<T> J{eye_m + hat_m * ((1 - std::cos(theta)) / theta)
        + hat_m * hat_m * ((theta - std::sin(theta)) / theta)};

    return SEX<T>(so, J * offset_vec);
}

template<typename T>
SEX<T> MatrixX<T>::SE_exp() const
{
    assert(row_ == col_);
    assert(row_ == 3 /* for SE2 */ || row_ == 4 /* for SE3 */);
    return vee().SE_Exp();
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

template<typename T>
T MatrixX<T>::det() const
{
    assert(row_ == col_);
    if (row_ == 1) {
        return array_.at(0);
    } else if (row_ == 2) {
        return at(0, 0) * at(1, 1) - at(1, 0) * at(0, 1);
    }

    // calculate determinant with cofactor expansions
    T result{0};
    MatrixX<T> minor(row_ - 1, col_ - 1);
    for (int c = 0; c < col_; ++c) {
        T sign{(1 + 1 + c) % 2 == 0 ? static_cast<T>(1) : static_cast<T>(-1)};
        T pivot{at(0, c)};

        // expansion for first row
        for (int mr = 1; mr < row_; ++mr) {
            int cc{0};
            for (int mc = 0; mc < col_; ++mc) {
                if (c == mc) continue;
                minor.at(mr - 1, cc++) = at(mr, mc);
            }
        }

        result += sign * pivot * minor.det();
    }
    return result;
}

template<typename T>
T MatrixX<T>::tr() const
{
    assert(row_ == col_);

    T result{0};
    for (int i = 0; i < row_; ++i) {
        result += at(i, i);
    }
    return result;
}

template<typename T>
MatrixX<T> MatrixX<T>::inv() const
{
    assert(row_ == col_ && "Invalid matrix shape");
    
    if (row_ == 1) {
        assert(!is_0(at(0, 0)) && "This matrix not invertible");
        
        MatrixX<T> out(row_, col_);
        out << 1 / at(0, 0);
        return out;
    }

    // Gauss Jordan algorithm
    auto out{MatrixX<T>::eye(row_, col_)};
    auto tmp_m{*this};

    for (int c = 0; c < out.col(); ++c) {
        if (tmp_m.at(c, c) == 0) {
            // swap row
            for (int r1 = c + 1; r1 < out.row(); ++r1) {
                if (tmp_m.at(r1, c) != 0) {
                    swap_row(tmp_m, c, r1);
                    swap_row(out, c, r1);
                }
            }
        }

        auto tmp{tmp_m.at(c, c)};
        assert(!is_0(tmp) && "This matrix not invertible");

        for (int c1 = 0; c1 < out.col(); ++c1) {
            tmp_m.at(c, c1) /= tmp;
            out.at(c, c1) /= tmp;
        }

        for (int r = 0; r < out.row(); ++r) {
            if (r == c) continue;
            const auto t{tmp_m.at(r, c)};
            for (int c1 = 0; c1 < out.col(); ++c1) {
                out.at(r, c1) -= out.at(c, c1) * t;
                tmp_m.at(r, c1) -= t * tmp_m.at(c, c1);
            }
        }
    }

    return out;
}

template<typename T>
MatrixX<T> MatrixX<T>::sub_mat(int sr, int sc, int r, int c) const
{
    assert(sr >= 0 && sc >= 0);
    assert((r + sr) <= row_ && (c + sc) <= col_ && "Bad sub matrix");

    MatrixX<T> out(r, c);
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            out.at(i, j) = at(sr + i, sc + j);
        }
    }
    return out;
}

template<typename T>
void MatrixX<T>::sub_mat_assign(int sr, int sc, int r, int c, const MatrixX& other)
{
    assert(sr >= 0 && sc >= 0);
    assert(r == other.row_ && c == other.col_);
    assert((r + sr) <= row_ && (c + sc) <= col_ && "Bad sub matrix");

    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            at(sr + i, sc + j) = other.at(i, j);
        }
    }
}

template<typename T>
VectorX<T> SOX<T>::Log() const
{
    auto theta{std::acos((tr() - 1) / 2)};
    auto eye_m{MatrixX<T>::eye(dim(), dim())};
    VectorX<T> norm_v{dim(), (*this - eye_m).null_space().normalize().array()};
    return norm_v * theta;
}

template<typename T>
MatrixX<T> SOX<T>::log() const
{
    return Log().hat();
}

template<typename T>
SOX<T> SOX<T>::plus(const VectorX<T>& v) const
{
    assert(dim() == v.row());
    return (*this) * v.SO_Exp();
}

template<typename T>
VectorX<T> SOX<T>::minus(const SOX& so) const
{
    assert(dim() == so.dim());
    return (so.inv() * (*this)).Log();
}

template<typename T>
SEX<T>::SEX(const SOX<T>& so, const VectorX<T>& off)
    : MatrixX<T>(so.dim() + 1, so.dim() + 1)
{
    // copy so
    sub_mat_assign(0, 0, so.dim(), so.dim(), so.mat());
    // copy offset
    sub_mat_assign(0, dim(), off.row(), off.col(), off.mat());
    // set const
    at(dim(), dim()) = 1;
}

template<typename T>
SEX<T> SEX<T>::operator*(const SEX<T>& other) const
{
    const auto this_rot{rot()};
    return SEX<T>(this_rot * other.rot(), offset() + this_rot * other.offset());
}

template<typename T>
MatrixX<T> SEX<T>::Adj()
{
    MatrixX<T> out(2 * dim(), 2 * dim());
    out.sub_mat_assign(0, 0, dim(), dim(), rot());
    out.sub_mat_assign(0, dim(), dim(), dim(), offset().hat() * rot());
    out.sub_mat_assign(dim(), dim(), dim(), dim(), rot());
    return out;
}

template<typename T>
SEX<T> SEX<T>::inv() const
{
    const auto this_rot_t{rot().t()};
    return SEX<T>(this_rot_t, this_rot_t * offset() * -1);
}

template<typename T>
VectorX<T> SEX<T>::offset() const
{
    return sub_mat(0, dim(), dim(), 1);
}

template<typename T>
SOX<T> SEX<T>::rot() const
{
    return SOX<T>(dim(), sub_mat(0, 0, dim(), dim()).array());
}

template<typename T>
VectorX<T> SEX<T>::operator*(const VectorX<T>& other) const
{
    return offset() + (rot() * other);
}

template<typename T>
VectorX<T> SEX<T>::Log() const
{
    VectorX<T> rot_vec(rot().Log());
    
    T theta{rot_vec.norm2()};
    assert(!is_0(theta) && "Invalid Exp operator");
    
    auto eye_m{MatrixX<T>::eye(dim(), dim())};
    auto hat_m{(rot_vec * (1 / theta)).hat()};

    MatrixX<T> J{eye_m + hat_m * ((1 - std::cos(theta)) / theta)
        + hat_m * hat_m * ((theta - std::sin(theta)) / theta)};

    VectorX<T> offset_vec(J.inv() * offset());

    // concat vector
    VectorX<T> out(offset_vec.row() + offset_vec.row());
    out.sub_mat_assign(0, 0, offset_vec.row(), offset_vec.col(), offset_vec);
    out.sub_mat_assign(offset_vec.row(), 0, rot_vec.row(), rot_vec.col(), rot_vec);
    return out;
}

template<typename T>
MatrixX<T> SEX<T>::log() const
{
    return Log().hat();
}

template<typename T>
SEX<T> SEX<T>::plus(const VectorX<T>& v) const
{
    assert(dim() == v.dim() / 2);
    auto v_exp{v.SE_Exp()};

    SOX<T> rot_m(rot() * v_exp.rot());
    VectorX<T> off_v(offset() + rot().mat() * v_exp.offset());
    return SEX<T>(rot_m, off_v);
}

template<typename T>
VectorX<T> SEX<T>::minus(const SEX<T>& other) const
{
    assert(dim() == other.dim());

    auto other_rot_t{other.rot().t()};
    VectorX<T> rot_vec((other_rot_t * rot()).Log());
    
    T theta{rot_vec.norm2()};
    assert(!is_0(theta) && "Invalid Exp operator");
    
    auto eye_m{MatrixX<T>::eye(dim(), dim())};
    auto hat_m{(rot_vec * (1 / theta)).hat()};
    MatrixX<T> J{eye_m + hat_m * ((1 - std::cos(theta)) / theta)
        + hat_m * hat_m * ((theta - std::sin(theta)) / theta)};

    VectorX<T> offset_vec(J.inv() * other_rot_t * (offset() - other.offset()));

    VectorX<T> out(offset_vec.dim() + rot_vec.dim());
    out.sub_mat_assign(0, 0, offset_vec.row(), offset_vec.col(), offset_vec);
    out.sub_mat_assign(offset_vec.row(), 0, rot_vec.row(), rot_vec.col(), rot_vec);
    return out;
}

template<typename T, int D>
SO<T, D> SO<T, D>::eye()
{
    SO<T, D> out;
    for (int i = 0; i < D; ++i) {
        out.at(i, i) = 1;
    }
    return out;
}

template<typename T, int D>
SE<T, D> SE<T, D>::eye()
{
    return SE<T, D>(SEX<T>::eye(D).array());
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const MatrixX<T>& m)
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

template<typename T>
std::ostream& operator<<(std::ostream& os, const SOX<T>& so)
{
    os << "SO" << so.dim() << typeid(T).name() << " {";
    os << std::fixed << std::setprecision(3);
    for (int i = 0; i < so.size(); ++i) {
        if (i % so.col() == 0) os << "\n  ";
        os << so.array().at(i) << ", ";
    }
    os << "\n}";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const SEX<T>& se)
{
    os << "SE" << se.dim() << typeid(T).name() << " {";
    os << std::fixed << std::setprecision(3);
    for (int i = 0; i < se.size(); ++i) {
        if (i % se.col() == 0) os << "\n  ";
        os << se.array().at(i) << ", ";
    }
    os << "\n}";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const VectorX<T>& v)
{
    os << "Vector" << v.dim() << typeid(T).name() << " {";
    os << std::fixed << std::setprecision(3);
    for (int i = 0; i < v.size(); ++i) {
        if (i % v.col() == 0) os << "\n  ";
        os << v.array().at(i) << ", ";
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