//
// Created by janos on 8/17/20.
//

#include "Mesh.h"

#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Trade/MeshData.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Matrix4.h>

namespace Phasefield {

using namespace Magnum;
using namespace Corrade;

Mesh::Mesh(Mg::Trade::MeshData meshData) {
    m_vertexCount = meshData.vertexCount();
    Containers::arrayResize(m_attributes, Containers::NoInit, stride*m_vertexCount);

}

void Mesh::requireGaussianCurvature() {
    if(m_gaussianCurvature.size() != vertexCount()) {
        requireInternalAngles();
        Containers::arrayResize(m_angles, 0);
        Containers::arrayReserve(m_angles, m_indices.size());
        for(auto face : faces()) {
            for(std::size_t i = 0; i < face.size(); ++i) {
                Vector3 p1 = vertices()[face[(i + face.size() - 1)%face.size()]];
                Vector3 p2 = vertices()[face[i]];
                Vector3 p3 = vertices()[face[(i + 1)%face.size()]];

                Containers::arrayAppend(m_angles, Math::angle((p1 - p2).normalized(), (p3 - p2).normalized()));
            }
        }
    }
}

void Mesh::requireInternalAngles() {
    if(m_angles.size() != vertexCount()) {
        Containers::arrayResize(m_angles, 0);
        Containers::arrayReserve(m_angles, m_indices.size());
        for(Face face : faces()) {
            Vector3 A = (face.c() - face.b()).normalized();
            Vector3 B = (face.a() - face.c()).normalized();
            Vector3 C = (face.b() - face.a()).normalized();
            Containers::arrayAppend(m_angles, {Math::angle(C, -B), Math::angle(A, -C), Math::angle(B, -A)});
        }
    }
}

StridedArrayView1D<Vector3> Mesh::positions() {
    return m_attributes.slice(&Implementation::Attributes::position);
}

StridedArrayView1D<Vector3> Mesh::normals() {
    return m_attributes.slice(&Implementation::Attributes::normal);
}

StridedArrayView1D<Float> Mesh::scalarField() {
    return m_attributes.slice(&Implementation::Attributes::scalar);
}

Corrade::Containers::StridedArrayView2D<Magnum::UnsignedInt> Mesh::faces() {
    return {m_indices, {m_faceCount, m_primitiveSize}};
}

Corrade::Containers::StridedArrayView2D<Rad> Mesh::angles() {
    return {m_angles, {m_faceCount, m_primitiveSize}};
}

Mesh::Mesh(Mg::Trade::MeshData meshData) {

}

bool Mesh::isManifold() {
    for (Edge e : edges()) {
        if (!e.isManifold()) return false;
    }
    for (Vertex v : vertices()) {
        if (!v.isManifold()) return false;
    }
    return true;
}

UnsignedInt Mesh::faceCount() const {
    return m_faceCount;
}

UnsignedInt Mesh::vertexCount() const {
    return m_vertexCount;
}

}
