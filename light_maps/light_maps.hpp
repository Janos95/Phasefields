// IBL light map
// Author: Max Schwarz <max.schwarz@ais.uni-bonn.de>


#include <Magnum/GL/CubeMapTexture.h>
#include <Magnum/GL/Texture.h>

#include <memory>
#include <string>

class LightMap
{
public:
    LightMap();
    explicit LightMap(const std::string& path);

    bool load(const std::string& path);

    constexpr const std::string& path() const
    { return m_path; }

    inline Magnum::GL::CubeMapTexture& irradianceMap()
    { return m_irradiance; }

    inline Magnum::GL::CubeMapTexture& prefilterMap()
    { return m_prefilter; }

    inline Magnum::GL::Texture2D& brdfLUT()
    { return m_brdfLUT; }

    inline Magnum::GL::CubeMapTexture& cubeMap()
    { return m_cubeMap; }

private:
    std::string m_path;

    Magnum::GL::CubeMapTexture m_cubeMap;
    Magnum::GL::CubeMapTexture m_irradiance;
    Magnum::GL::CubeMapTexture m_prefilter;
    Magnum::GL::Texture2D m_brdfLUT;
};


