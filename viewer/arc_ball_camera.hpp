//
// Created by janos on 2/16/20.
//


#pragma once

#include "arc_ball.hpp"

#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/AbstractTranslationRotation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>

/* Arcball camera implementation integrated into the SceneGraph */
class ArcBallCamera: public ArcBall {
public:
    template<class Transformation> ArcBallCamera(
            SceneGraph::Scene<Transformation>& scene,
            const Vector3& cameraPosition, const Vector3& viewCenter,
            const Vector3& upDir, Deg fov, const Vector2i& windowSize,
            const Vector2i& viewportSize):
            ArcBall{cameraPosition, viewCenter, upDir, fov, windowSize}
    {
        /* Create a camera object of a concrete type */
        auto* cameraObject = new SceneGraph::Object<Transformation>{&scene};
        (*(_camera = new SceneGraph::Camera3D{*cameraObject}))
                .setAspectRatioPolicy(SceneGraph::AspectRatioPolicy::Extend)
                .setProjectionMatrix(Matrix4::perspectiveProjection(
                        fov, Vector2{windowSize}.aspectRatio(), 0.01f, 100.0f))
                .setViewport(viewportSize);

        /* Save the abstract transformation interface and initialize the
           camera position through that */
        (*(_cameraObject = cameraObject))
                .rotate(transformation().rotation())
                .translate(transformation().translation());
    }

    /* Update screen and viewport size after the window has been resized */
    void reshape(const Vector2i& windowSize, const Vector2i& viewportSize) {
        _windowSize = windowSize;
        _camera->setViewport(viewportSize);
    }

    /* Update the SceneGraph camera if arcball has been changed */
    bool update() {
        /* call the internal update */
        if(!updateTransformation()) return false;

        (*_cameraObject)
                .resetTransformation()
                .rotate(transformation().rotation())
                .translate(transformation().translation());
        return true;
    }

    /* Draw objects using the internal scenegraph camera */
    void draw(SceneGraph::DrawableGroup3D& drawables) {
        _camera->draw(drawables);
    }

    Vector2i viewport() const {
        return _camera->viewport();
    }

    auto projectionMatrix() const {
        return _camera->projectionMatrix();
    }

private:
    SceneGraph::AbstractTranslationRotation3D* _cameraObject{};
    SceneGraph::Camera3D* _camera{};
};
