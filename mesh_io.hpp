//
// Created by janos on 03.04.20.
//

#pragma once

#include "viewer.hpp"

#include <string>
#include <mutex>

class MeshIO : public Viewer::AbstractEventHandler {
public:

    void drawEvent() override;
    bool loadMesh(std::string const& path);
    bool loadBlob(std::string const& path);
    bool saveMesh(std::string const& path);
    bool saveBlob(std::string const& path);

    std::vector<char> m_inputMesh;
    std::vector<char> m_outputMesh;

    std::vector<char> m_inputBlob;
    std::vector<char> m_outputBlob;

    Scene& m_scene;
};



