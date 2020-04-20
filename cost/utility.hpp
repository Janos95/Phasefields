//
// Created by janos on 2/23/20.
//

#pragma once


#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/StlMath.h>

#include <algorithm>

class F
{
public:
    F(const double a, const double b):
            m_a(a), m_b(b), m_c1(1./std::pow(1.+a, 2)), m_c2(1./std::pow(b-1., 2))
    {
    }

    template<class T>
    T operator()(const T& x) const {
        if(x < m_a)
            return pow(x - m_a, 2) * m_c1;
        if(x < m_b)
            return 0.;

        return pow(m_b - x, 2) * m_c2;
    }


private:
    const double m_a,m_b,m_c1,m_c2;
};


class W
{
public:

    W(const double a, const double b): m_a(a), m_b(b), m_c3(-30. / std::pow(a - b, 5))
    {
    }

    template<class T>
    T operator()(const T& x) const{
        if(x < m_a)
            return 0.;
        if(x < m_b)
            return pow(x-m_a, 2) * pow(x - m_b, 2) * m_c3;

        return 0.;
    }
private:
    const double m_a,m_b,m_c3;
};


class FGrad
{
public:
    FGrad(const double a, const double b):
            m_a(a), m_b(b), m_c1(1./std::pow(1.+a, 2)), m_c2(1./std::pow(b-1., 2))
    {
    }

    inline double operator()(const double x) const {
        if(x < m_a)
            return 2.*(x-m_a) * m_c1;
        if(x < m_b)
            return 0;

        return 2.*(x - m_b) * m_c2;
    }


private:
    const double m_a,m_b,m_c1,m_c2;
};


class WGrad
{
public:

    WGrad(const double a, const double b): m_a(a), m_b(b), m_c3(-30. / std::pow(a - b, 5))
    {
    }

    inline double operator()(const double x) const{
        if(x < m_a)
            return 0.;
        if(x < m_b)
            return 2. * ((x - m_a) * pow(x - m_b, 2)  + pow(x - m_a, 2) * (x - m_b)) * m_c3;
        return 0.;
    }
private:
    const double m_a,m_b,m_c3;
};


auto computeAreas(const Containers::ArrayView<const Vector3ui>& indices, const Containers::StridedArrayView1D<const Vector3>& data){
    Containers::Array<Float> areas(Containers::NoInit, indices.size());
    for (int j = 0; j < indices.size(); ++j) {
        auto& triangle = indices[j];
        auto a = data[triangle[0]] - data[triangle[1]];
        auto b = data[triangle[0]] - data[triangle[2]];
        areas[j] = Math::cross(a, b).length() * .5f;
    }
}

template<class R, class T>
void updateWeight(const int target, const T& w, R& neighbors){
    auto it = std::find_if(neighbors.begin(), neighbors.end(), [target](const auto& n){ return n.vertex == target; });
    assert(it != neighbors.end());
    it->weight = w;
}

class StoppingCriteria
{
public:

    StoppingCriteria() = default;

    StoppingCriteria(const int source, const int numComponents, const Containers::Array<int>& components):
        m_startComponent(components[source]),
        m_numComponentsToFind(numComponents - m_startComponent - 1),
        m_components(&components),
        m_found(Containers::DirectInit, m_numComponentsToFind + 1, false),
        m_targetVertices(Containers::DirectInit, m_numComponentsToFind + 1, -1)
    {
    }

    bool operator()(int target)
    {
        auto comp = (*m_components)[target];
        auto idx = comp - m_startComponent;
        //if target is not in the interface then comp is -1.
        if(comp <= m_startComponent || m_found[idx])
            return false;

        m_found[idx] = true;
        m_targetVertices[idx] = target;
        auto stop = ++m_numComponentsFound == m_numComponentsToFind;
        return stop;
    }

    bool foundAll() const {
        return std::all_of(m_found.begin()+1, m_found.end(), [](const auto& x){ return x; });
    }

    int target(int i) const{
        CORRADE_INTERNAL_ASSERT(i > m_startComponent);
        auto target = m_targetVertices[i - m_startComponent];
        CORRADE_INTERNAL_ASSERT(target >= 0);
        return target;
    }

private:
    int m_startComponent;
    int m_numComponentsToFind;
    const Containers::Array<int>* m_components;

    int m_numComponentsFound = 0;
    Corrade::Containers::Array<bool> m_found;
    Corrade::Containers::Array<int> m_targetVertices;
};



