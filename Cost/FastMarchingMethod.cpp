//
// Created by janos on 8/17/20.
//

#include "FastMarchingMethod.h"
#include "Heap.h"
#include "GraphCommon.hpp"

#include <Magnum/Math/Functions.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Corrade/Containers/Array.h>

#include <cstring>

namespace Phasefield {

using namespace Corrade;
using namespace Magnum;

namespace {

// The super fun quadratic distance function in the Fast Marching Method on triangle meshes
// TODO parameter c isn't actually defined in paper, so I guessed that it was an error
double eikonalDistanceSubroutine(double a, double b, double theta, double dA, double dB) {
    if(theta <= Math::Constants<Double>::pi()/2.0) {
        double u = dB - dA;
        double cTheta = std::cos(theta);
        double sTheta2 = 1.0 - cTheta*cTheta;

        // Quadratic equation
        double quadA = a*a + b*b - 2*a*b*cTheta;
        double quadB = 2*b*u*(a*cTheta - b);
        double quadC = b*b*(u*u - a*a*sTheta2);
        double sqrtVal = std::sqrt(quadB*quadB - 4*quadA*quadC);
        // double tVals[] = {(-quadB + sqrtVal) / (2*quadA),        // seems to always be the first one
        //                   (-quadB - sqrtVal) / (2*quadA)};

        double t = (-quadB + sqrtVal)/(2*quadA);
        if(u < t && a*cTheta < b*(t - u)/t && b*(t - u)/t < a/cTheta) {
            return dA + t;
        } else {
            return Math::min(b + dA, a + dB);
        }

    }
        // Custom by Nick to get acceptable results in obtuse triangles without fancy unfolding
    else {

        double maxDist = Math::max(dA, dB); // all points on base are less than this far away, by convexity
        double c = Math::sqrt(a*a + b*b - 2*a*b*std::cos(theta));
        double area = 0.5*std::sin(theta)*a*b;
        double altitude = 2*area/c; // distance to base, must be inside triangle since obtuse
        double baseDist = maxDist + altitude;

        return Math::min({b + dA, a + dB, baseDist});
    }
}


} // namespace

struct FastMarchingMethod::Impl {

    struct Entry {
        UnsignedInt idx;
        Double distance;
    };

    explicit Impl(Mesh& m) :
        mesh(m),
        distances(Containers::NoInit, mesh.vertexCount()),
        finalized(Containers::NoInit, mesh.vertexCount()) {

        mesh.requireCornerAngles();
        mesh.requireEdgeLengths();

        reset();
    }

    void reset(){
        for(auto& x: distances) x = 0.;
        std::memset(finalized.data(), 0, finalized.size());
        frontierPQ.clear();
        nFound = 0;
    }

    void setSource(UnsignedInt idx){
        frontierPQ.emplace(idx, 0.);
    }

    bool step(UnsignedInt& idx, Double& distance);

    Mesh& mesh;

    Containers::Heap<Entry> frontierPQ;

    // Initialize
    Containers::Array<double> distances;
    Containers::Array<char> finalized;

    std::size_t nFound;
};

bool FastMarchingMethod::Impl::step(UnsignedInt& idx, Double& distance) {

    // TODO this could handle nonmanifold geometry with a few small tweaks
    CORRADE_ASSERT(mesh.isManifold(), "handling of nonmanifold mesh not yet implemented",{});

    size_t nVert = mesh.vertexCount();

    // Search
    while(nFound < nVert && !frontierPQ.empty()) {

        // Pop the nearest element
        auto item = frontierPQ.extractMin();
        idx = item.idx;
        distance = item.distance;

        // Accept it if not stale
        if(finalized[idx]) {
            continue;
        }
        distances[idx] = distance;
        finalized[idx] = true;
        nFound++;

        // Add any eligible neighbors
        for(Halfedge he : currV.incomingHalfedges()) {
            Vertex neighVert = he.vertex();

            // Add with length
            if(!finalized[neighVert]) {
                double newDist = currDist + geometry.edgeLengths[he.edge()];
                if(newDist < distances[neighVert]) {
                    frontierPQ.push(std::make_pair(currDist + geometry.edgeLengths[he.edge()], neighVert));
                    distances[neighVert] = newDist;
                }
                continue;
            }

            // Check the third point of the "left" triangle straddling this edge
            if(he.isInterior()) {
                Vertex newVert = he.next().next().vertex();
                if(!finalized[newVert]) {

                    // Compute the distance
                    double lenB = geometry.edgeLengths[he.next().next().edge()];
                    double distB = currDist;
                    double lenA = geometry.edgeLengths[he.next().edge()];
                    double distA = distances[neighVert];
                    double theta = geometry.cornerAngles[he.next().next().corner()];
                    double newDist = eikonalDistanceSubroutine(lenA, lenB, theta, distA, distB);

                    if(newDist < distances[newVert]) {
                        frontierPQ.push(std::make_pair(newDist, newVert));
                        distances[newVert] = newDist;
                    }
                }
            }

            // Check the third point of the "right" triangle straddling this edge
            Halfedge heT = he.twin();
            if(heT.isInterior()) {
                Vertex newVert = heT.next().next().vertex();
                if(!finalized[newVert]) {

                    // Compute the distance
                    double lenB = geometry.edgeLengths[heT.next().edge()];
                    double distB = currDist;
                    double lenA = geometry.edgeLengths[heT.next().next().edge()];
                    double distA = distances[neighVert];
                    double theta = geometry.cornerAngles[heT.next().next().corner()];
                    double newDist = eikonalDistanceSubroutine(lenA, lenB, theta, distA, distB);

                    if(newDist < distances[newVert]) {
                        frontierPQ.push(std::make_pair(newDist, newVert));
                        distances[newVert] = newDist;
                    }
                }
            }
        }
    }

    return distances;
}


FastMarchingMethod::FastMarchingMethod(Mesh& mesh) : m_impl(Containers::pointer<Impl>(mesh)) {}

FastMarchingMethod::Result FastMarchingMethod::tick() { return m_impl->tick(); }

}