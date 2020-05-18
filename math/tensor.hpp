//
// Created by janos on 16.05.20.
//

#pragma once

#include <Magnum/Magnum.h>
#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Math.h>
#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>

#include <cblas.h>


namespace Mg = Magnum;
namespace Cr = Corrade;

template<class T, class U = void>
struct choose_tag {
    using type = Cr::Containers::ValueInitT;
};

template<class T>
struct choose_tag<T, std::enable_if_t<std::is_trivial_v<T>>> {
    using type = Cr::Containers::NoInitT;
};

template<class T>
using choose_tag_t = typename choose_tag<T>::type;

template<typename Scalar>
class Matrix {

    template<class TagType>
    explicit Matrix(TagType, Mg::Int rows, Mg::Int cols) :
        m_data(TagType{}, rows * cols),
        m_rows(rows),
        m_cols(cols)
    {
    }

    explicit Matrix(Mg::Int rows, Mg::Int cols) :
        m_data(choose_tag_t<Scalar>{}, rows * cols),
        m_rows(rows),
        m_cols(cols)
    {
    }

    Matrix& operator=(Matrix other) noexcept {
        other.swap(*this);
        return *this;
    }

    Matrix(Matrix const& other) :
        m_rows(other.m_rows),
        m_cols(other.m_cols),
        m_data(choose_tag_t<Scalar>{}, other.size())
    {
        Cr::Utility::copy(other.m_data, m_data);
    }

    Matrix(Matrix&& other) noexcept :
        m_data(std::move(other.m_data)),
        m_rows(other.m_rows),
        m_cols(other.m_cols)
    {
    }

    void swap(Matrix& other){
        auto rows = other.m_rows;
        auto cols = other.m_cols;
        auto temp = std::move(other.m_data);
        other.m_data = m_data;
        other.m_rows = m_rows;
        other.m_cols = m_cols;
        m_data = std::move(temp);
        m_rows = rows;
        m_cols = cols;
    }

    Scalar& operator()(Mg::Int x, Mg::Int y){
        return m_data[x * m_rows + y];
    };

    friend Matrix operator*(Matrix const&, Matrix const&);

    [[nodiscard]] constexpr inline Mg::Int rows() const { return m_rows; }
    [[nodiscard]] constexpr inline Mg::Int cols() const { return m_cols; }
    [[nodiscard]] constexpr inline Mg::Int size() const { return m_cols * m_rows; }

    Matrix inverted(){

    }

    Scalar determinant(){

    }


private:
    Mg::Int m_rows = 0, m_cols = 0;
    Cr::Containers::Array<Scalar> m_data;
};

template<class Scalar>
Matrix<Scalar> operator*(Matrix<Scalar> const& A, Matrix<Scalar> const& B){
    auto rows = A.m_rows;
    auto n = A.m_cols;
    auto cols = B.m_cols;
    CORRADE_ASSERT(n == B.m_rows, "Matrices not compatible for product", {});

    Matrix<Scalar> C(Cr::Containers::ValueInit, rows, cols);
    for (int x = 0; x < cols; ++x) {
        for (int y = 0; y < rows; ++y) {
            for (int i = 0; i < n; ++i) {
                C(y,x) += A(y, i) * B(i, x);
            }
        }
    }

    return C;
}

template<unsigned dimension, class Scalar>
struct Map
{
    Map(Cr::Containers::StridedArrayView2D<Scalar> const& view) :
        m_view(view)
    {
    }

    [[nodiscard]] constexpr inline Mg::Int rows() const { return m_view.size()[0]; }
    [[nodiscard]] constexpr inline Mg::Int cols() const { return m_view.size()[1]; }
    [[nodiscard]] constexpr inline Mg::Int size() const { return m_view.size().sum(); }

private:
    Cr::Containers::StridedArrayView2D<Scalar> m_view;
};


template<class Scalar>
Matrix<Scalar> operator+(Matrix<Scalar> const& A, Matrix<Scalar> const& B){
    auto rows = A.rows();
    auto cols = A.cols();
    CORRADE_ASSERT(rows == B.m_rows, "Matrix Sum : Not the same number of rows", {});
    CORRADE_ASSERT(cols == B.m_cols, "Matrix Sum : Not the same number of cols", {});

    Matrix<Scalar> C(choose_tag_t<Scalar>{}, rows, cols);
    for (int x = 0; x < cols; ++x) {
        for (int y = 0; y < rows; ++y) {
            C(x,y) = A(x,y) + B(x,y);
        }
    }
    return C;
}

template<class Scalar>
struct HouseholderQR {

    HouseholderQR(Matrix<Scalar> const& matrix){
        Matrix copy(matrix);

    }

    Scalar determinant(){

    }


};