//
// Created by janos on 8/25/20.
//

#include "EigenHelper.h"
#include "LinearizedElasticity.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Tree.h"
#include "Functional.hpp"
#include "FEM.hpp"
#include "VisualizationProxy.h"
#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Matrix3.h>

#include <imgui.h>

namespace Phasefield {

LinearizedElasticity::LinearizedElasticity(Mesh& m) : mesh(m),
    displacementsX(m.vertexCount()), displacementsY(m.vertexCount()), positions(m.vertexCount())
{
    mesh.requireIntegralOperator();
    mesh.requireGradientOperator();

    for(Vertex v : mesh.vertices())
        positions[v] = Vector2d{v.position().xy()};

    for(size_t i = 0; i < neumannBoundary.size(); ++i) {
        Vertex v = neumannBoundary[i];
        HalfEdge hes[2];
        size_t j = 0;
        for(HalfEdge he : v.outgoingHalfEdges()) {
            if(he.edge().onBoundaryLoop())
                hes[j++] = he;
        }
        CORRADE_ASSERT(j == 3, "More than two boundary edge at vertex",);
        double boundaryElement = 0.5*(hes[0].asVector().length() + hes[1].asVector().length());
        neumannElements[i] = boundaryElement;
    }
}

double LinearizedElasticity::surfaceLoad() {
    const Vector2d force{0., -forceNorm};

    double load = 0;

    for(size_t i = 0; i < neumannBoundary.size(); ++i) {
        Vertex v = neumannBoundary[i];
        Vector2d displacement{displacementsX[v], displacementsY[v]};
        load += neumannElements[i]*Math::dot(force, displacement);
    }

    return load;
}

template<class Scalar>
void LinearizedElasticity::operator()(
         ArrayView<const Scalar> parameters, ArrayView<const Scalar> weights, Scalar& out,
         ArrayView<Scalar> gradP, ArrayView<Scalar> gradW) {


}

size_t LinearizedElasticity::numParameters() const { return mesh.vertexCount(); }

void LinearizedElasticity::drawImGuiOptions(VisualizationProxy& proxy) {
    bool redraw = false;

    if(ImGui::Checkbox("Apply Displacements", &applyDisplacements)) {
        if(applyDisplacements) {
            proxy.setCallbacks(
                    [this, &proxy](Node) {
                        for(Vertex v : mesh.vertices())
                            mesh.position(v).xy() += Vector2(displacementsX[v], displacementsY[v]);
                    },
                    [this]{ applyDisplacements = false; });

        } else {
            for(Vertex v : mesh.vertices())
                mesh.position(v).xy() += Vector2{positions[v]};
        }
        redraw = true;
    }

    static const double min = 0.0001;
    static const double max = 1;
    if(ImGui::DragScalar("Poisson Ratio", ImGuiDataType_Double, &poissonRatio, 1.f, &min, &max, "%f", 1)) {
        if(applyDisplacements) {
            computeDisplacements();
            redraw = true;
        }
    }

    if(ImGui::DragScalar("Youngs Modulus", ImGuiDataType_Double, &youngModulus, 1.f, &min, &max, "%f", 1)) {
        if(applyDisplacements) {
            computeDisplacements();
            redraw = true;
        }
    }

    if(redraw) proxy.redraw();
}

DEFINE_FUNCTIONAL_CONSTRUCTOR(LinearizedElasticity)
DEFINE_FUNCTIONAL_OPERATOR(LinearizedElasticity, double)

#ifdef PHASEFIELD_WITH_ADOLC
DEFINE_FUNCTIONAL_OPERATOR(LinearizedElasticity, adouble)
#endif

}
