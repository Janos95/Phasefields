//
// Created by janos on 09.12.19.
//

#pragma once

#include "interpolation.hpp"
#include "quadrature_ref_triangle.hpp"
#include "dijkstra.hpp"
#include "utility.hpp"
#include "union_find.hpp"
#include "bfs.hpp"
#include "detach.hpp"
#include "small_array.hpp"

#include <scoped_timer/scoped_timer.hpp>

#include <Magnum/Magnum.h>
#include <Corrade/Containers/EnumSet.h>

#include <ceres/first_order_function.h>

#include <Eigen/Core>
#include <igl/doublearea.h>
#include <igl/edges.h>

#include <ratio>
#include <numeric>


enum class GradientFlag : Magnum::UnsignedByte {
    Analytic = 1 << 0,
    Automatic = 1 << 1,
};


template<class Scalar>
class ConnectednessConstraint : public ceres::FirstOrderFunction
{
public:
    ConnectednessConstraint() noexcept;
    
    ConnectednessConstraint(ConnectednessConstraint&&) noexcept;
    
    ConnectednessConstraint& operator=(ConnectednessConstraint&&);
    
    //deleted copy constructor
    ConnectednessConstraint(ConnectednessConstraint&) = delete;
    
    //deleted copy assigment
    ConnectednessConstraint& operator=(ConnectednessConstraint&) = delete;
    
    ConnectednessConstraint(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F,
            const double epsilon,
            const double a,
            const double b);

    void setPreimageInterval(double a, double b);

    bool Evaluate(double const* parameters, double* cost, double* jacobian) const override ;

    int NumParameters() const override;

private:

    struct Neighbor {
        Neighbor(int v, double w) : vertex(v), weight(w) {}

        int vertex;
        Scalar weight;
    };

    struct Edge {
        Edge(const int v1_, const int v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

        int v1, v2;

#ifdef __cpp_impl_three_way_comparison
        auto operator<=>(const Edge&) const = default;
#else

        bool operator<(const Edge& other) const{
            return std::tie(v1, v2) < std::tie(other.v1, other.v2);
        }

        bool operator==(const Edge& other) const{
            return std::tie(v1, v2) == std::tie(other.v1, other.v2);
        }
#endif
    };

    const Eigen::MatrixXd &m_V;
    const Eigen::MatrixXi& m_F;

    std::vector<Edge> m_dualEdges;
    mutable std::vector<SmallArray<3, Neighbor>> m_adjacencyList;

    std::vector<double> m_lineElements;
    Eigen::VectorXd m_areas;
    std::vector<double> m_diams;

    double m_epsilon;
    double m_a, m_b;
};

template<class Scalar>
ConnectednessConstraint<Scalar>::ConnectednessConstraint(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const double epsilon,
        const double a,
        const double b):
        m_V(V),
        m_F(F),
        m_epsilon(epsilon), m_a(a), m_b(b) {

    CORRADE_INTERNAL_ASSERT(m_F.colwise().maxCoeff().maxCoeff() < m_V.size());
    //compute edges in dual graph
    std::map<Edge, int> edgeMap;
    m_dualEdges.reserve(3 * F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            auto[it, inserted] = edgeMap.try_emplace(Edge{m_F(i, j), m_F(i, (j + 1) % 3)}, i);
            if (!inserted) {
                //both faces share edge F(i,j) - F(i,j+1mod3)
                m_dualEdges.emplace_back(it->second, i);
            }
        }
    }

    m_lineElements.resize(m_dualEdges.size(), 0);
    m_adjacencyList.resize(m_F.rows(), SmallArray<3, Neighbor>{Containers::NoInit});
    for (auto[v1, v2] : m_dualEdges) {
        m_adjacencyList[v1].emplace_back(v2, .0);
        m_adjacencyList[v2].emplace_back(v1, .0);
    }

    m_diams.resize(m_F.rows());
    for (int i = 0; i < m_F.rows(); ++i) {
        auto f = m_F.row(i);
        m_diams[i] = triangleDiameter(m_V.row(f[0]), m_V.row(f[2]), m_V.row(f[2]));
    }

    igl::doublearea(m_V, m_F, m_areas);
    m_areas *= 0.5;

#ifndef NODEBUG
    BreadthFirstSearch bfs(m_adjacencyList, 0);
    bfs.run();
    CORRADE_INTERNAL_ASSERT(bfs.isConnected());
#endif
}

template<class Scalar>
bool ConnectednessConstraint<Scalar>::Evaluate(double const *params,
                            double *cost,
                            double *jacobian) const {
    using VectorT = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatrixT = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    ScopedTimer t("Connectedness", true);

    F weight(m_a, m_b);
    FGrad weightGrad(m_a, m_b);
    W bump(m_a, m_b);
    WGrad bumpGrad(m_a, m_b);

    int numFaces = m_F.rows();
    int numVertices = m_V.rows();

    VectorT U(numVertices);
    std::copy_n(params, numVertices, U.begin());

    if constexpr(detail::IsDiffArray<Scalar>) {
        if (jacobian)
            for (auto &u : U)
                enoki::set_requires_gradient(u);
    }

    UnionFind set(numFaces);
    std::vector inInterface(numFaces, false);
    std::vector ws(numFaces, Scalar(-1.));
    std::vector uT(numFaces, Scalar(0.));

    for (int i = 0; i < numFaces; i++) {
        auto f = m_F.row(i);
        uT[i] = 1. / 3. * (U[f[0]] + U[f[1]] + U[f[2]]);
        bool in = m_a <= detail::detach(uT[i]) && detail::detach(uT[i]) <= m_b;
        if (in) {
            ws[i] = bump(uT[i]);
            inInterface[i] = true;
        }
    }

    for (auto[dualV1, dualV2] : m_dualEdges) {
        auto fuT1 = weight(uT[dualV1]);
        auto fuT2 = weight(uT[dualV2]);
        auto w = .5 * (m_diams[dualV1] + m_diams[dualV2]) * .5 * (fuT1 + fuT2);
        updateWeight(dualV1, w, m_adjacencyList[dualV2]);
        updateWeight(dualV2, w, m_adjacencyList[dualV1]);
        if (std::abs(detail::detach(w)) < std::numeric_limits<double>::epsilon())
            set.unite(dualV1, dualV2);
    }

    std::vector<int> components(set.size());
    std::vector<int> roots;
    for (std::size_t i = 0; i < components.size(); ++i) {
        if (inInterface[i]) {
            components[i] = set.find(i);
            roots.push_back(components[i]);
        } else
            components[i] = -1;
    }

    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());

    auto numComponents = roots.size();
    if (numComponents <= 1){
        *cost = 0.;
        if(jacobian)
            std::fill_n(jacobian, numVertices, 0.);
        return true;
    }

    Debug{} << "Found" << numComponents <<  "connencted components. Enforcing Connectedness!\n";

    std::vector W(numComponents, Scalar(0.));

    for (int i = 0; i < numFaces; ++i) {
        auto &c = components[i];
        if (c != -1) {
            auto it = std::lower_bound(roots.begin(), roots.end(), c);
            CORRADE_INTERNAL_ASSERT(it != roots.end() && *it == c);
            auto k = std::distance(roots.begin(), it);
            CORRADE_INTERNAL_ASSERT(std::abs(-1 - detail::detach(ws[i])) > 1e-6);
            W[k] += ws[i] * m_areas[i];
            c = k;
        }
    }

    //run dijkstra from each connected component except last one
    std::vector dijkstras(numComponents - 1, Dijkstra(m_adjacencyList));
    std::vector stops(numComponents - 1, StoppingCriteria{});
    {
        ScopedTimer t("dijkstra", true);
        for (std::size_t i = 0; i < numComponents - 1; ++i) {
            StoppingCriteria stop(roots[i], numComponents, components);
            dijkstras[i].setSource(roots[i]);
            dijkstras[i].run({stop});
            CORRADE_INTERNAL_ASSERT(stop.foundAll());
            stops[i] = std::move(stop);
        }
    }

    Scalar d = 0;
    MatrixT D;
    if constexpr(!detail::IsDiffArray<Scalar>)
        D = MatrixT(numComponents, numComponents);
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(numVertices);

    {
        ScopedTimer tDiff("connectedness diff", true);

        for (std::size_t i = 0; i < numComponents; ++i) {
            for (std::size_t j = i + 1; j < numComponents; ++j) {

                Scalar dij = 0;
                auto Wij = W[i] * W[j];

                for (auto&&[a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))) {

                    auto &av = m_adjacencyList[a];
                    auto it = std::find_if(av.begin(), av.end(), [b = b](const auto &n) { return b == n.vertex; });
                    dij += it->weight;

                    if (jacobian && !detail::IsDiffArray<Scalar>) {
                        auto lineElement = .5 * (m_diams[a] + m_diams[b]);
                        const double fgrada = weightGrad(detail::detach(uT[a]));
                        const double fgradb = weightGrad(detail::detach(uT[b]));
                        auto weightedFGrad = detail::detach(Wij) * lineElement * .5 * (1. / 3.);
                        for (auto v : m_F.row(a))
                            gradient[v] += weightedFGrad * fgrada;
                        for (auto v : m_F.row(b))
                            gradient[v] += weightedFGrad * fgradb;
                    }
                }

                d += dij * Wij;
                if (jacobian && !detail::IsDiffArray<Scalar>) {
                    D(i, j) = dij;
                    D(j, i) = dij;
                }
            }
        }

        if constexpr (!detail::IsDiffArray<Scalar>) {
            if (jacobian) {
                for (int k = 0; k < numFaces; ++k) {
                    if (components[k] < 0) //not in interface
                        continue;
                    auto wgrad = bumpGrad(uT[k]);
                    if (std::abs(wgrad) < std::numeric_limits<double>::epsilon())
                        continue;
                    std::size_t i = components[k];
                    for (std::size_t j = 0; j < numComponents; ++j) {
                        if (i == j)
                            continue;
                        double weightedGrad = D(i, j) * W[j] * wgrad * m_areas[k] / 3.;
                        //each incident vertex has the same influence
                        for (auto v : m_F.row(k)) {
                            gradient[v] += weightedGrad;
                        }
                    }
                }
            }
        }
    }

    auto scaleFactor = (1. / std::pow(m_epsilon, 2)) * 2.; // 1/(eps * eps)
    d *= scaleFactor;
    *cost = detail::detach(d);

    if (jacobian) {
        if constexpr(detail::IsDiffArray<Scalar>) {
            Eigen::VectorXd autodiffgrad = Eigen::VectorXd::Zero(numVertices);
            Scalar::simplify_graph_();
            enoki::backward(d);
            for (int i = 0; i < numVertices; ++i) {
                jacobian[i] = enoki::gradient(U[i]);
            }
        } else {
            Eigen::Map<Eigen::VectorXd> map(jacobian, numVertices);
            map = scaleFactor * gradient;
        }
    }
    return true;
}

template<class Scalar>
int ConnectednessConstraint<Scalar>::NumParameters() const {
    return m_V.rows();
}
