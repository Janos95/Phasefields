//
// Created by janos on 8/17/20.
// taken from geometry central

#include "FastMarchingMethod.h"
#include "Types.h"
#include "MeshElements.h"

#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Math/Functions.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Magnum/Math/Constants.h>

#include <cstring>
#include <limits>

namespace Phasefield {

using namespace Corrade;
using namespace Magnum;

namespace {

// The super fun quadratic distance function in the Fast Marching Method on triangle meshes
// TODO parameter c isn't actually defined in paper, so I guessed that it was an error
double eikonalDistanceSubroutine(double a, double b, Radd theta, double dA, double dB) {
    if(theta <= Radd{Math::Constants<double>::piHalf() + 0.1}) {
        double u = dB - dA;
        double cTheta = Math::cos(theta);
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
        //Debug{} << "Obtuse triangle";
        double maxDist = Math::max(dA, dB); // all points on base are less than this far away, by convexity
        double c = Math::sqrt(a*a + b*b - 2*a*b*Math::cos(theta));
        double area = 0.5*Math::sin(theta)*a*b;
        double altitude = 2*area/c; // distance to base, must be inside triangle since obtuse
        double baseDist = maxDist + altitude;

        return Math::min({b + dA, a + dB, baseDist});
    }
}


} // namespace


FastMarchingMethod::FastMarchingMethod(Mesh& mesh) :
    m_mesh(mesh),
    m_distances(NoInit, mesh.vertexCount()),
    m_finalized(NoInit, mesh.vertexCount())
{
    //CORRADE_ASSERT(mesh.isManifold(), "Mesh needs to be manifold", );
    mesh.requireAngles();
    mesh.requireEdgeLengths();
    reset();
}

void FastMarchingMethod::reset(){
    for(auto& x: m_distances) x = std::numeric_limits<double>::infinity();
    std::memset(m_finalized.data(), 0, m_finalized.size());
    m_frontier.clear();
}

void FastMarchingMethod::update() {

    m_mesh.requireAngles();
    m_mesh.requireEdgeLengths();

    arrayResize(m_finalized, NoInit, m_mesh.vertexCount());
    arrayResize(m_distances, NoInit, m_mesh.vertexCount());

    reset();
}

void FastMarchingMethod::setSource(Vertex v){
    m_frontier.emplace(v, 0.);
}

bool FastMarchingMethod::step(Vertex& vertex, Double& distance) {

    size_t nVert = m_mesh.vertexCount();

    /* get next non finished item */
    do {
        if(m_frontier.empty()) return false;

        // Pop the nearest element
        auto item = m_frontier.extractMin();
        vertex = item.node;
        distance = item.distance;

    } while(m_finalized[vertex]);

    m_distances[vertex] = distance;
    m_finalized[vertex] = true;

    // Add any eligible neighbors
    for(HalfEdge he : vertex.incomingHalfEdges()) {
        Vertex neighbor = he.tail();
        // Add with length
        if(!m_finalized[neighbor]) {
            double newDist = distance + m_mesh.edgeLength[he.edge()];
            if(newDist < m_distances[neighbor]) {
                m_frontier.emplace(neighbor, distance + m_mesh.edgeLength[he.edge()]);
                m_distances[neighbor] = newDist;
            }
            continue;
        }

        // Check the third point of the "left" triangle straddling this edge
        if(he.isInterior()) {
            Vertex newVert = he.next().next().tail();
            if(!m_finalized[newVert]) {

                // Compute the distance
                double lenB = m_mesh.edgeLength[he.next().next().edge()];
                double distB = distance;
                double lenA = m_mesh.edgeLength[he.next().edge()];
                double distA = m_distances[neighbor];
                Radd theta = m_mesh.angle[he.next().next().corner()];
                double newDist = eikonalDistanceSubroutine(lenA, lenB, Radd{theta}, distA, distB);

                if(newDist < m_distances[newVert]) {
                    m_frontier.emplace(newVert, newDist);
                    m_distances[newVert] = newDist;
                }
            }
        }

        // Check the third point of the "right" triangle straddling this edge
        HalfEdge heT = he.twin();
        if(heT.isInterior()) {
            Vertex newVert = heT.next().next().tail();
            if(!m_finalized[newVert]) {

                // Compute the distance
                double lenB = m_mesh.edgeLength[heT.next().edge()];
                double distB = distance;
                double lenA = m_mesh.edgeLength[heT.next().next().edge()];
                double distA = m_distances[neighbor];
                Radd theta = m_mesh.angle[heT.next().next().corner()];
                double newDist = eikonalDistanceSubroutine(lenA, lenB, Radd{theta}, distA, distB);

                if(newDist < m_distances[newVert]) {
                    m_frontier.emplace(newVert, newDist);
                    m_distances[newVert] = newDist;
                }
            }
        }
    }

    return true;
}

}