//
// Created by janos on 30.04.20.
//

#include "fem.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/GrowableArray.h>

using namespace Corrade;
using namespace Magnum;

auto computeArea(Vector3d const& a, Vector3d const& b, Vector3d const& c){
    auto area = Math::cross(b - a, c - a).length() * .5;
    return area;
}

Containers::Array<Double> computeAreas(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices){
    Containers::Array<Double> areas(Containers::NoInit, triangles.size());
    for (int j = 0; j < triangles.size(); ++j) {
        auto& t = triangles[j];
        areas[j] = computeArea(vertices[t[0]],vertices[t[1]],vertices[t[2]]);
    }
    return areas;
}


Containers::Array<Double> computeIntegralOperator(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices){
    Containers::Array<Double> integral(Containers::ValueInit, vertices.size());
    for (int i = 0; i < triangles.size(); ++i) {
        auto const& t = triangles[i];
        auto area = computeArea(vertices[t[0]],vertices[t[1]],vertices[t[2]]);
        for (int j = 0; j < 3; ++j)
            integral[t[j]] += area / 3.;
    }
    return integral;
}

Eigen::SparseMatrix<Mg::Double> computeMassMatrix(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices) {

    Containers::Array<Eigen::Triplet<Double>> triplets;
    Containers::arrayReserve(triplets, 9 * triangles.size());
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        auto const& t = triangles[i];

        auto a = vertices[t[1]] - vertices[t[2]];
        auto b = vertices[t[2]] - vertices[t[0]];

        auto area = Math::cross(a, -b).length() * .5;

        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                auto factor = (j == k) ? 2. : 1.;
                Containers::arrayAppend(triplets, Containers::InPlaceInit, t[j], t[k], factor * area);
            }
        }
    }
    Eigen::SparseMatrix<Double> massMatrix(vertices.size(), vertices.size());
    massMatrix.setFromTriplets(triplets.begin(), triplets.end());
    massMatrix.makeCompressed();
    return massMatrix;
}


Eigen::SparseMatrix<Double> gradient(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices,
        bool uniform)
{
    // Number of faces
    const auto m = triangles.size();
    // Number of vertices
    const auto nv = vertices.size();

    Containers::Array<Eigen::Triplet<Double>> triplets;
    Containers::arrayReserve(triplets, 4*3*m);

    for (std::size_t i=0; i < m; ++i)
    {
        auto i1 = triangles[i][0];
        auto i2 = triangles[i][1];
        auto i3 = triangles[i][2];

        // #F x 3 matrices of triangle edge vectors, named after opposite vertices
        auto v32 = vertices[i3] - vertices[i2];
        auto v13 = vertices[i1] - vertices[i3];
        auto v21 = vertices[i2] - vertices[i1];
        auto n = Math::cross(v32, v13);

        Double dblA = n.length();
        Vector3d u(0,0,1);
        if (!uniform) {
            // now normalize normals to get unit normals
            u = n / dblA;
        } else {
            // Abstract equilateral triangle v1=(0,0), v2=(h,0), v3=(h/2, (sqrt(3)/2)*h)
            // get h (by the area of the triangle)
            Float h = sqrt(dblA/sin(Math::Constants<Float>::pi() / 3.0)); // (h^2*sin(60))/2. = Area => h = sqrt(2*Area/sin_60)

            Vector3d v1(0,0,0);
            Vector3d v2(h,0,0);
            Vector3d v3(h/2.,(sqrt(3)/2.)*h,0);

            // now fix v32,v13,v21 and the normal
            v32 = v3-v2;
            v13 = v1-v3;
            v21 = v2-v1;
            n = Math::cross(v32, v13);
        }

        // rotate each vector 90 degrees around normal
        auto eperp21 = Math::cross(u, v21).normalized() * v21.dot() / dblA;
        auto eperp13 = Math::cross(u, v13).normalized() * v13.dot() / dblA;

        for(int d = 0; d < 3; d++)
        {
            Containers::arrayAppend(triplets, Containers::InPlaceInit, m+d*m, triangles[m][1], eperp13[d]);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, m+d*m, triangles[m][0], -eperp13[d]);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, m+d*m, triangles[m][2], eperp21[d]);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, m+d*m, triangles[m][1], -eperp13[d]);
        }
    }

    Eigen::SparseMatrix<Double> G(3*m, nv);
    G.setFromTriplets(triplets.begin(), triplets.end());
    G.makeCompressed();
    return G;
}

Eigen::SparseMatrix<Double> computeStiffnessMatrix(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices) {
    auto m = triangles.size();
    auto n = vertices.size();
    Containers::Array<Eigen::Triplet<Double>> triplets;
    Containers::arrayReserve(triplets, 12 * m);

    // run over all faces
    for (std::size_t i = 0; i < m; ++i) {
        auto const &t = triangles[i];

        auto a = vertices[t[1]] - vertices[t[2]];
        auto b = vertices[t[2]] - vertices[t[0]];
        auto c = vertices[t[0]] - vertices[t[1]];

        auto area = Math::cross(a, -b).length() * .5;
        Vector3d lengthSqr(a.dot(), b.dot(), c.dot());

        // add local contribution
        for (int j = 0; j < 3; j++) {
            auto entry = 0.125 * (lengthSqr[j] - lengthSqr[(j + 1) % 3] - lengthSqr[(j + 2) % 3]) / area;
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 1) % 3], t[(j + 2) % 3], entry);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 2) % 3], t[(j + 1) % 3], entry);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 1) % 3], t[(j + 1) % 3], -entry);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 2) % 3], t[(j + 2) % 3], -entry);
        }
    }

    Eigen::SparseMatrix<Double> stiffnessMatrix(n, n);
    stiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());
    stiffnessMatrix.makeCompressed();
    return stiffnessMatrix;
}