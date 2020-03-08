//
// Created by janos on 27.02.20.
//
#include "color_callback.hpp"
#include "colormaps.hpp"

#include "object.hpp"
#include "scene.hpp"

#include <Eigen/Core>

#include <ceres/iteration_callback.h>

#include <Corrade/Containers/Reference.h>
#include <Corrade/Containers/StridedArrayView.h>

#include <Corrade/Utility/Algorithms.h>
#include <Magnum/MeshTools/Interleave.h>

#include <mutex>
#include <algorithm>

using namespace Magnum;
using namespace Corrade;


ColorCallback::ColorCallback(Eigen::VectorXd& U, Containers::Array<char>&& data):
    m_U(U),
    m_vertexData(std::move(data))
{
}

ceres::CallbackReturnType ColorCallback::operator()(const ceres::IterationSummary& summary)
{
    std::lock_guard lock(m_mutex);
    if(!m_current){
        Containers::StridedArrayView1D<void> erasedView(
                m_vertexData,
                m_vertexData.begin() + 6 * sizeof(Float) /*start*/,
                m_U->size() /*size */ ,
                10 * sizeof(Float) /* stride */);
        auto view = arrayCast<Color4>(erasedView);
        std::transform(m_U->begin(), m_U->end(), view.begin(), jet_colormap);
        m_current = false;
    }
    return ceres::SOLVER_CONTINUE;
}

void ColorCallback::operator()(Scene& scene){
    if(!m_current){
        auto& vertices = scene.getObject("mesh")->vertices;
        auto data = vertices.map(
                0,
                m_vertexData.size(),
                GL::Buffer::MapFlag::Write|GL::Buffer::MapFlag::InvalidateBuffer);
                /* GL::Buffer::MapFlag::Read|GL::Buffer::MapFlag::Write)); */
        CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
        {
            std::lock_guard lock(m_mutex);
            Utility::copy(m_vertexData, data);
            m_current = true;
        }
        vertices.unmap();
    }
}