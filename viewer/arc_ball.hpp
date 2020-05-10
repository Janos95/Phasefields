//
// Created by janos on 2/16/20.
//

#pragma once

#include <Magnum/Math/Quaternion.h>
#include <Magnum/Math/DualQuaternion.h>
#include <Magnum/Magnum.h>

namespace Mg = Magnum;
/* Implementation of Ken Shoemake's arcball camera with smooth navigation
   feature: https://www.talisman.org/~erlkonig/misc/shoemake92-arcball.pdf */
class ArcBall {
public:
    ArcBall(const Mg::Vector3& cameraPosition, const Mg::Vector3& viewCenter,
            const Mg::Vector3& upDir, Mg::Deg fov, const Mg::Vector2i& windowSize);

    /* Set the camera view parameters: eye position, view center, up
       direction */
    void setViewParameters(const Mg::Vector3& eye, const Mg::Vector3& viewCenter,
                           const Mg::Vector3& upDir);

    /* Reset the camera to its initial position, view center, and up dir */
    void reset();

    /* Update screen size after the window has been resized */
    void reshape(const Mg::Vector2i& windowSize) { _windowSize = windowSize; }

    /* Update any unfinished transformation due to lagging, return true if
       the camera matrices have changed */
    bool updateTransformation();

    /* Get/set the amount of lagging such that the camera will (slowly)
       smoothly navigate. Lagging must be in [0, 1) */
    Mg::Float lagging() const { return _lagging; }
    void setLagging(Mg::Float lagging);

    /* Initialize the first (screen) mouse position for camera
       transformation. This should be called in mouse pressed event. */
    void initTransformation(const Mg::Vector2i& mousePos);

    /* Rotate the camera from the previous (screen) mouse position to the
       current (screen) position */
    void rotate(const Mg::Vector2i& mousePos);

    /* Translate the camera from the previous (screen) mouse position to
       the current (screen) mouse position */
    void translate(const Mg::Vector2i& mousePos);

    /* Translate the camera by the delta amount of (NDC) mouse position.
       Note that NDC position must be in [-1, -1] to [1, 1]. */
    void translateDelta(const Mg::Vector2& translationNDC);

    /* Zoom the camera (positive delta = zoom in, negative = zoom out) */
    void zoom(Mg::Float delta);

    /* Get the camera's view transformation as a qual quaternion */
    const Mg::DualQuaternion& view() const { return _view; }

    /* Get the camera's view transformation as a matrix */
    Mg::Matrix4 viewMatrix() const { return _view.toMatrix(); }

    /* Get the camera's inverse view matrix (which also produces
       transformation of the camera) */
    Mg::Matrix4 inverseViewMatrix() const { return _inverseView.toMatrix(); }

    /* Get the camera's transformation as a dual quaternion */
    const Mg::DualQuaternion& transformation() const { return _inverseView; }

    /* Get the camera's transformation matrix */
    Mg::Matrix4 transformationMatrix() const { return _inverseView.toMatrix(); }

    /* Return the distance from the camera position to the center view */
    Mg::Float viewDistance() const { return Mg::Math::abs(_targetZooming); }

protected:
    /* Update the camera transformations */
    void updateInternalTransformations();

    /* Transform from screen coordinate to NDC - normalized device
       coordinate. The top-left of the screen corresponds to [-1, 1] NDC,
       and the bottom right is [1, -1] NDC. */
    Mg::Vector2 screenCoordToNDC(const Mg::Vector2i& mousePos) const;

    Mg::Deg _fov;
    Mg::Vector2i _windowSize;

    Mg::Vector2 _prevMousePosNDC;
    Mg::Float _lagging{};

    Mg::Vector3 _targetPosition, _currentPosition, _positionT0;
    Mg::Quaternion _targetQRotation, _currentQRotation, _qRotationT0;
    Mg::Float _targetZooming, _currentZooming, _zoomingT0;
    Mg::DualQuaternion _view, _inverseView;
};

