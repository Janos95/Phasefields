//
// Created by janos on 10/11/20.
//

#include "Bvh.h"
#include "Mesh.h"

#include <bvh/bvh.hpp>
#include <bvh/vector.hpp>
#include <bvh/triangle.hpp>
#include <bvh/ray.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/primitive_intersectors.hpp>

using Triangle = bvh::Triangle<float>;
using Bvh      = bvh::Bvh<float>;
using Ray      = bvh::Ray<float>;

namespace Phasefield {

struct BVHAdapter::Impl {
    Mesh& mesh;
    FaceData<Triangle> triangles;
    Bvh bvh;

    bool computeIntersection(Vector3 const& p, Vector3 const& dir, Intersection&);
    void update();

};

void BVHAdapter::Impl::update() {
    arrayResize(triangles, mesh.faceCount());
    for(Face f : mesh.faces()) {
        bvh::Vector3<float> ps[3];
        HalfEdge he = f.halfEdge();
        size_t i = 0;
        for(Vertex v : f.vertices()) {
            auto p = v.position();
            ps[i++] = bvh::Vector3<float>{p.x(), p.y(), p.z()};
        }
        triangles[f] = Triangle{ps[0], ps[1], ps[2]};
    }

    bvh::SweepSahBuilder<Bvh> builder(bvh);
    auto[bboxes, centers] = bvh::compute_bounding_boxes_and_centers(triangles.data(), triangles.size());
    auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), triangles.size());
    builder.build(global_bbox, bboxes.get(), centers.get(), triangles.size());
}


bool BVHAdapter::Impl::computeIntersection(Vector3 const& p, Vector3 const& dir, Intersection& intersection) {
    bvh::ClosestPrimitiveIntersector<Bvh, Triangle> intersector{bvh, triangles.data()};
    bvh::SingleRayTraverser<Bvh> traverser(bvh);

    Ray ray(
            {p.x(), p.y(), p.z()}, // origin
            {dir.x(), dir.y(), dir.z()}, // direction
            0.0,                    // minimum distance
            100.0                   // maximum distance
    );
    auto hit = traverser.traverse(ray, intersector);
    if (hit) {
        intersection.t = hit->intersection.t;
        intersection.u = hit->intersection.t;
        intersection.v = hit->intersection.v;
        return true;

    } else return false;
}

BVHAdapter::BVHAdapter(Mesh& mesh) : MeshFeature(mesh, false), impl(new Impl{mesh}) { impl->update(); }

bool BVHAdapter::computeIntersection(Vector3 const& p, Vector3 const& dir, Intersection& intersection) {
    return impl->computeIntersection(p, dir, intersection);
}

void BVHAdapter::update() {
    impl->update();
}

BVHAdapter::~BVHAdapter() {
    delete impl;
}

}
