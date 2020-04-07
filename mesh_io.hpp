//
// Created by janos on 03.04.20.
//

#pragma once

#include "viewer.hpp"
#include "phasefield_data.hpp"

#include <string>
#include <mutex>

class MeshIO : public Viewer::AbstractEventHandler {
public:

    explicit MeshIO(PhasefieldData& data):
        m_phasefieldData(data),
        m_inputMesh(100,0),
        m_outputMesh(100,0),
        m_inputBlob(100,0),
        m_outputBlob(100,0)
        {
        }

    void drawImGui() override;
    bool loadMesh(std::string const& path);
    bool loadBlob(std::string const& path);
    bool saveMesh(std::string const& path);
    bool saveBlob(std::string const& path);

private:
    std::vector<char> m_inputMesh;
    std::vector<char> m_outputMesh;

    std::vector<char> m_inputBlob;
    std::vector<char> m_outputBlob;

    PhasefieldData& m_phasefieldData;
};



