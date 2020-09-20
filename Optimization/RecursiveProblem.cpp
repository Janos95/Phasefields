//
// Created by janos on 7/16/20.
//

#include "RecursiveProblem.h"
#include "SparseMatrix.h"
#include "C1Functions.h"
#include "Types.h"

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Math/Vector3.h>

namespace Phasefield::Solver {


RecursiveProblem::RecursiveProblem(Tree& t) : tree(t) {}

size_t RecursiveProblem::numParameters() const { return tree.level(levelToOptimize).size(); }

size_t RecursiveProblem::numConstraints() const { return tree.numLeafs*constraints.size(); }

void RecursiveProblem::determineSparsityStructure(SparseMatrix& jacobian) const {
    //jacobian.clear();
    //jacobian.numCols = numParameters();
    //jacobian.numRows = numConstraints()*tree.numLeafs;

    //size_t m = numConstraints();
    //size_t n = numParameters();
    //size_t row = 0;
    //tree.traverse([&](Node const& leaf) {
    //    if(leaf.isLeaf()){
    //        Node const* node = &leaf;
    //        while(true){
    //            for(size_t j = 0; j < m; ++j, ++row){
    //                for(size_t i = 0; i < n; ++i){
    //                    arrayAppend(jacobian.rows, row);
    //                    size_t col = node->idx*n + i;
    //                    arrayAppend(jacobian.cols, col);
    //                }
    //            }
    //            if(node->m_parent != Invalid){
    //                break;
    //            } else {
    //                node = &tree.nodes[node->m_parent];
    //            }
    //        }
    //    }
    //    return true;
    //});
    //jacobian.nnz = jacobian.rows.size();
    //arrayResize(jacobian.values, jacobian.nnz);
}

void RecursiveProblem::operator()(
        Containers::ArrayView<const double> data,
        double& cost,
        Containers::ArrayView<double> gradData,
        Containers::ArrayView<double> constr,
        SparseMatrix* jacobian) const {

    //SmootherStep smoothStep;

    //auto& nodes = tree.nodes;
    //auto n = tree.vertexCount();
    //auto m = numConstraints();
    //auto size = nodes.size();

    ////auto gradients = tree.gradients();
    //auto& jacValues = jacobian->values;
    //auto prefixes = tree.temporaryData();
    //StridedArrayView2D<const double> phasefields{data, {size, n}};
    //StridedArrayView2D<double> gradients{gradData, {size, n}};

    //cost = 0.;
    //for(auto& x : prefixes[0]) x = 1.;

    //int jacIdx = 0;
    //tree.traverse([&](Node& node) -> bool {
    //    if(!node.isLeaf()){
    //        for(UnsignedInt i = 0; i < n; ++i){
    //            if(node.m_leftChild != Invalid){
    //                double pos = smoothStep.eval(phasefields[node.m_leftChild][i]);
    //                prefixes[node.m_leftChild][i] = pos*prefixes[node.idx][i];
    //            }
    //            if(node.m_rightChild != Invalid){
    //                double neg = smoothStep.eval(-phasefields[node.m_leftChild][i]);
    //                prefixes[node.m_rightChild][i] = neg*prefixes[node.idx][i];
    //            }
    //        }
    //    } else {
    //        Node const& leaf = node;

    //        Array<double> productsGrad(Containers::ValueInit, n);
    //        Array<double> postfix(Containers::NoInit, n);
    //        Array<double> productsJacData(Containers::ValueInit, n*m);
    //        Containers::StridedArrayView2D<double> productsJac(productsJacData, {m, n});

    //        for(Functional const& functional : objectives){
    //            functional(phasefields[node.idx].asContiguous(),
    //                       prefixes[node.idx].asContiguous(), cost,
    //                       gradients[leaf.idx].asContiguous(), productsGrad);
    //        }

    //        for(std::size_t i = 0; i < constraints.size(); ++i){
    //            auto jacSlice = jacValues.slice(jacIdx, jacIdx + n);
    //            constraints[i](phasefields[node.idx].asContiguous(),
    //                           prefixes[node.idx].asContiguous(), cost,
    //                           jacSlice, productsJac[i].asContiguous());
    //            jacIdx += n;
    //        }

    //        for(double& x : postfix) x = 1.;

    //        Node const* node1 = &leaf;
    //        while(true) {
    //            size_t idx = node1->idx;
    //            double sign = tree.isLeftChild(*node1);
    //            for(UnsignedInt i = 0; i < n; ++i){
    //                gradients[idx][i] += prefixes[idx][i]*smoothStep.grad(sign*phasefields[idx][i])*postfix[i]*productsGrad[i];
    //                postfix[i] *= phasefields[idx][i];
    //            }

    //            for(std::size_t j = 0; j < m; ++j) {
    //                for(UnsignedInt i = 0; i < n; ++i) {
    //                    jacValues[jacIdx++] += prefixes[idx][i]*smoothStep.grad(sign*phasefields[idx][i])*postfix[i]*productsJac[j][i];
    //                }
    //            }

    //            if(node1->m_parent != Invalid){
    //                break;
    //            } else {
    //                node1 = &tree.nodes[node1->m_parent];
    //            }
    //        }
    //    }
    //    return true;
    //});
}



void RecursiveProblem::operator()(ArrayView<const double> parameters, double& cost, ArrayView<double> gradient) const {
    size_t n = tree.vertexCount();
    SmootherStep smootherStep;

    size_t nodeCountOnLevel = tree.nodeCountOnLevel(levelToOptimize);
    StridedArrayView2D<double> gradientsOnLevel{gradient, {nodeCountOnLevel, n}};
    StridedArrayView2D<const double> phasesOnLevel{parameters, {nodeCountOnLevel, n}};

    for(double& x : tree.root().temporary()) x = 1.;
    if(gradient)
        for(double& x : gradient) x = 0.;
    cost = 0.;

    for(Node node : tree.ancestorsOfLevel(levelToOptimize)) {
        auto weights = node.temporary();
        for(size_t i = 0; i < n; ++i) {
            if(node.hasLeftChild()) {
                Node leftChild = node.leftChild();
                leftChild.temporary()[i] = smootherStep.eval(node.phasefield()[i])*weights[i];
            }
            if(node.hasRightChild()) {
                Node rightChild = node.rightChild();
                rightChild.temporary()[i] = smootherStep.eval(-node.phasefield()[i])*weights[i];
            }
        }
    }

    size_t i = 0;
    for(Node node : tree.nodesOnLevel(levelToOptimize)) {
        for(Functional const& f : objectives) {
            f(phasesOnLevel[i].asContiguous(), node.temporary(), cost, gradientsOnLevel[i].asContiguous(), nullptr);
        }
        ++i;
    }
}

}
