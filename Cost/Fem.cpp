//
// Created by janos on 30.04.20.
//

#include "Fem.h"

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Math/Vector3.h>
#include <Corrade/Containers/StridedArrayView.h>

namespace Phasefield {

using namespace Corrade;
using namespace Magnum;

Double computeArea(Vector3d const& a, Vector3d const& b, Vector3d const& c) {
    auto area = Math::cross(b - a, c - a).length()*.5;
    return area;
}

Containers::Array<Double> computeAreas(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices) {
    Containers::Array<Double> areas(Containers::NoInit, triangles.size());
    for(std::size_t j = 0; j < triangles.size(); ++j) {
        auto& t = triangles[j];
        areas[j] = computeArea(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
    }
    return areas;
}


Containers::Array<Double> computeIntegralOperator(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices) {
    Containers::Array<Double> integral(Containers::ValueInit, vertices.size());
    for(std::size_t i = 0; i < triangles.size(); ++i) {
        auto const& t = triangles[i];
        auto area = computeArea(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
        for(int j = 0; j < 3; ++j)
            integral[t[j]] += area/3.;
    }
    return integral;
}

Containers::Array<Triplet> computeMassMatrix(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices) {

    Containers::Array<Triplet> triplets;
    Containers::arrayReserve(triplets, 9*triangles.size());
    for(std::size_t i = 0; i < triangles.size(); ++i) {
        auto const& t = triangles[i];

        auto a = vertices[t[1]] - vertices[t[2]];
        auto b = vertices[t[2]] - vertices[t[0]];

        auto area = Math::cross(a, -b).length()*.5;

        for(int j = 0; j < 3; j++) {
            for(int k = 0; k < 3; k++) {
                auto factor = (j == k) ? 2. : 1.;
                Containers::arrayAppend(triplets, Containers::InPlaceInit, t[j], t[k], factor*area);
            }
        }
    }

    return triplets;
}


Containers::Array<Mg::Vector3d> gradient(
        const Containers::ArrayView<const Vector3ui>& triangles,
        const Containers::ArrayView<const Vector3d>& vertices) {
    const std::size_t m = triangles.size();

    Containers::Array<Mg::Vector3d> data(Containers::NoInit, 3*m);
    Containers::StridedArrayView2D<Mg::Vector3d> edgeNormals{data, {m, 3}};

    for(std::size_t i = 0; i < m; ++i) {
        auto i1 = triangles[i][0];
        auto i2 = triangles[i][1];
        auto i3 = triangles[i][2];

        auto e32 = vertices[i3] - vertices[i2];
        auto e13 = vertices[i1] - vertices[i3];
        auto e21 = vertices[i2] - vertices[i1];
        auto n = Math::cross(e32, e13);

        Double dblA = n.length();

        /* rotate each vector 90 degrees around normal */
        edgeNormals[i][0] = Math::cross(n, e32).normalized()*e32.length()/dblA;
        edgeNormals[i][1] = Math::cross(n, e13).normalized()*e13.length()/dblA;
        edgeNormals[i][2] = Math::cross(n, e21).normalized()*e21.length()/dblA;
    }

    return data;
}

Containers::Array<Triplet> computeStiffnessMatrix(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices) {

    auto m = triangles.size();
    auto n = vertices.size();
    Containers::Array<Triplet> triplets;
    Containers::arrayReserve(triplets, 12*m);

    // run over all faces
    for(std::size_t i = 0; i < m; ++i) {
        auto const& t = triangles[i];

        auto a = vertices[t[1]] - vertices[t[2]];
        auto b = vertices[t[2]] - vertices[t[0]];
        auto c = vertices[t[0]] - vertices[t[1]];

        auto area = Math::cross(a, -b).length()*.5;
        Vector3d lengthSqr(a.dot(), b.dot(), c.dot());

        // add local contribution
        for(int j = 0; j < 3; j++) {
            auto entry = 0.125*(lengthSqr[j] - lengthSqr[(j + 1)%3] - lengthSqr[(j + 2)%3])/area;
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 1)%3], t[(j + 2)%3], entry);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 2)%3], t[(j + 1)%3], entry);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 1)%3], t[(j + 1)%3], -entry);
            Containers::arrayAppend(triplets, Containers::InPlaceInit, t[(j + 2)%3], t[(j + 2)%3], -entry);
        }
    }

    return triplets;
}

}