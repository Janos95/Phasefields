//
// Created by janos on 7/16/20.
//

#include "RecursiveProblem.h"
#include "../Cost/SparseMatrix.h"
#include "c1_functions.hpp"

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Math/Vector3.h>

using namespace Magnum;

namespace Phasefield::Solver {

RecursiveProblem::RecursiveProblem(Tree& t) : tree(t) {
    updateInternalDataStructures();
}

Mg::UnsignedInt RecursiveProblem::numParameters() const {
    return tree.phasefieldData.size();
}

Mg::UnsignedInt RecursiveProblem::numConstraints() const {
    return tree.numLeafs*constraints.size();
}

void RecursiveProblem::updateInternalDataStructures() {
    for(Functional& f : objectives)
        f.updateInternalDataStructures();
    for(Functional& f : constraints)
        f.updateInternalDataStructures();
}

void RecursiveProblem::determineSparsityStructure(SparseMatrix& jacobian) const {
    jacobian.clear();
    jacobian.numCols = numParameters();
    jacobian.numRows = numConstraints()*tree.numLeafs;

    auto m = numConstraints();
    auto n = numParameters();
    UnsignedInt row = 0;
    tree.traverse([&](Node const& leaf) {
        if(leaf.isLeaf()){
            Node const* node = &leaf;
            while(true){
                for(UnsignedInt j = 0; j < m; ++j, ++row){
                    for(UnsignedInt i = 0; i < n; ++i){
                        Containers::arrayAppend(jacobian.rows, row);
                        auto col = node->idx*n + i;
                        Containers::arrayAppend(jacobian.cols, col);
                    }
                }
                if(node->parent != Node::None){
                    break;
                } else {
                    node = &tree.nodes[node->parent];
                }
            }
        }
    });
    jacobian.nnz = jacobian.rows.size();
    Containers::arrayResize(jacobian.values, jacobian.nnz);
}

void RecursiveProblem::operator()(
        Containers::ArrayView<const double> data,
        double& cost,
        Containers::ArrayView<double> gradData,
        Containers::ArrayView<double> constr,
        SparseMatrix* jacobian) const {

    SmootherStep smoothStep;

    auto& nodes = tree.nodes;
    auto n = tree.phasefieldSize;
    auto m = numConstraints();
    auto size = nodes.size();

    //auto gradients = tree.gradients();
    auto& jacValues = jacobian->values;
    auto prefixes = tree.temps();
    Containers::StridedArrayView2D<const double> phasefields{data, {size, n}};
    Containers::StridedArrayView2D<double> gradients{gradData, {size, n}};

    cost = 0.;
    for(auto& x : prefixes[0]) x = 1.;

    int jacIdx = 0;
    tree.traverse([&](Node& node) -> void {
        if(!node.isLeaf()){
            for(UnsignedInt i = 0; i < n; ++i){
                if(node.leftChild != Node::None){
                    double pos = smoothStep.eval(phasefields[node.leftChild][i]);
                    prefixes[node.leftChild][i] = pos*prefixes[node.idx][i];
                }
                if(node.rightChild != Node::None){
                    double neg = smoothStep.eval(-phasefields[node.leftChild][i]);
                    prefixes[node.rightChild][i] = neg*prefixes[node.idx][i];
                }
            }
        } else {
            Node const& leaf = node;

            Containers::Array<double> productsGrad(Containers::ValueInit, n);
            Containers::Array<double> postfix(Containers::NoInit, n);
            Containers::Array<double> productsJacData(Containers::ValueInit, n*m);
            Containers::StridedArrayView2D<double> productsJac(productsJacData, {m, n});

            for(Functional const& functional : objectives){
                functional(phasefields[node.idx].asContiguous(),
                           prefixes[node.idx].asContiguous(), cost,
                           gradients[leaf.idx].asContiguous(), productsGrad);
            }

            for(std::size_t i = 0; i < constraints.size(); ++i){
                auto jacSlice = jacValues.slice(jacIdx, jacIdx + n);
                constraints[i](phasefields[node.idx].asContiguous(),
                               prefixes[node.idx].asContiguous(), cost,
                               jacSlice, productsJac[i].asContiguous());
                jacIdx += n;
            }

            for(auto& x : postfix) x = 1.;

            Node const* node1 = &leaf;
            while(true) {
                auto idx = node1->idx;
                double sign = tree.isLeftChild(*node1);
                for(UnsignedInt i = 0; i < n; ++i){
                    gradients[idx][i] += prefixes[idx][i]*smoothStep.grad(sign*phasefields[idx][i])*postfix[i]*productsGrad[i];
                    postfix[i] *= phasefields[idx][i];
                }

                for(std::size_t j = 0; j < m; ++j) {
                    for(UnsignedInt i = 0; i < n; ++i) {
                        jacValues[jacIdx++] += prefixes[idx][i]*smoothStep.grad(sign*phasefields[idx][i])*postfix[i]*productsJac[j][i];
                    }
                }

                if(node1->parent != Node::None){
                    break;
                } else {
                    node1 = &tree.nodes[node1->parent];
                }
            }
        }
    });
}

}
