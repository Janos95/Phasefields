//
// Created by janos on 27.11.19.
//

#pragma once

#include "scene.hpp"

#include <Eigen/Core>

#include <ceres/iteration_callback.h>

#include <Corrade/Containers/Reference.h>
#include <Corrade/Containers/StridedArrayView.h>

#include <mutex>


enum class LayoutFlag: Magnum::UnsignedByte {
    SmoothNormals = 1 << 1,
    FlatNormals = 1 << 2,
};

using LayoutFlags = Corrade::Containers::EnumSet<LayoutFlag>;

CORRADE_ENUMSET_OPERATORS(LayoutFlags)


class ColorCallback : public ceres::IterationCallback
{
public:
    explicit ColorCallback(
            Eigen::VectorXd& U,
            Corrade::Containers::Array<char>&& data);

    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override;

    void operator()(Scene& scene);

private:

    std::atomic_bool m_current = false;
    LayoutFlags m_layout;
    Magnum::UnsignedInt m_stride;

    Corrade::Containers::Reference<Eigen::VectorXd> m_U;
    Corrade::Containers::Array<char> m_vertexData;

    mutable std::mutex m_mutex;
};