//
// Created by janos on 10/11/20.
//

#pragma once

#include "Mesh.h"

#include "FastBVH.h"

#include <Magnum/Math/Range.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/Optional.h>

namespace Phasefield {

//why do people write such interfaces ?? I just dont get it!
struct BVHAdapter : MeshFeature {

    Mesh& mesh;

    struct Converter {
        Mesh& mesh;

        FastBVH::BBox<float> operator()(Face const& triangle) {
            auto bb = triangle.bb();
            Vector3 min = bb.min();
            Vector3 max = bb.max();
            return {{min.x(), min.y(), min.z()}, {max.x(), max.y(), max.z()}};
        }
    };

    struct Intersector {
        Mesh& mesh;

        FastBVH::Intersection<float, Face> operator()(const Face& triangle, FastBVH::Ray<float> const& ray) const noexcept {
            constexpr static float EPSILON = 0.0000001;
            auto pos = triangle.positions();
            Vector3 vertex0 = pos[0];
            Vector3 vertex1 = pos[1];
            Vector3 vertex2 = pos[2];
            Vector3 rayVector{ray.d.x, ray.d.y, ray.d.z};
            Vector3 rayOrigin{ray.o.x, ray.o.y, ray.o.z};
            Vector3 edge1, edge2, h, s, q;
            float a,f,u,v;
            edge1 = vertex1 - vertex0;
            edge2 = vertex2 - vertex0;
            h = Math::cross(rayVector,edge2);
            a = Math::dot(edge1,h);
            if (a > -EPSILON && a < EPSILON)
                return {};    // This ray is parallel to this triangle.
            f = 1.0f/a;
            s = rayOrigin - vertex0;
            u = f * Math::dot(s,h);
            if (u < 0.0 || u > 1.0)
                return {};
            q = Math::cross(s,edge1);
            v = f * Math::dot(rayVector,q);
            if (v < 0.0 || u + v > 1.0)
                return {};
            // At this stage we can compute t to find out where the intersection point is on the line.
            float t = f * Math::dot(edge2,q);
            if (t > EPSILON) // ray intersection
            {
                return {t, &triangle};
            }
            else // This means that there is a line intersection but not a ray intersection.
                return {};
        }
    };

    Converter converter;
    FaceData<Face> primitives;

    Optional<FastBVH::BVH<float, Face>> bvh;
    Optional<FastBVH::Traverser<float, Face, Intersector>> traverser;

    explicit BVHAdapter(Mesh& m) : MeshFeature(m, false),  mesh(m), converter{m} {
        updateImpl();
    }

    size_t computeIntersection(Vector3 const& p, Vector3 const& dir) {
        if(primitives.empty()) return Invalid;
        FastBVH::Ray<float> ray{{p.x(), p.y(), p.z()}, {dir.x(), dir.y(), dir.z()}};
        auto isect = traverser->traverse(ray);
        Face const* f = isect.object;
        if(!f) return Invalid;
        Vertex v = *(f->vertices().begin());
        return v.idx;
    }

    void updateImpl() {
        if(mesh.vertexCount() == 0) return; // FastBVH does not handle empty mesh
        arrayResize(primitives, mesh.faceCount());
        for(Face f : mesh.faces()) {
            primitives[f] = f;
        }
        FastBVH::Iterable<Face> view{primitives.data(), primitives.size()};
        FastBVH::BuildStrategy<float, 1> strategy;
        bvh.emplace(strategy(view, converter));

        Intersector intersector{mesh};
        traverser.emplace(*bvh, intersector);
    }

    void update() override { updateImpl(); }

    FEATURE_NAME("BVH")

};

}

