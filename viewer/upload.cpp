//
// Created by janos on 13.03.20.
//

#include "upload.hpp"

#include <Magnum/GL/TextureFormat.h>
#include <Magnum/GL/Sampler.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/Duplicate.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/GenerateFlatNormals.h>
#include <Magnum/MeshTools/GenerateNormals.h>
#include <Magnum/Shaders/Generic.h>

using namespace Magnum;
using namespace Corrade;



void upload(DrawableData& drawableData) {
    auto& [vertices, indices, mesh, meshData] = drawableData;

    vertices.setData(meshData.vertexData());
    indices.setData(meshData.indexData());

    CORRADE_ASSERT((!meshData.isIndexed() || indices.id()) && vertices.id(),
                   "upload: invalid external buffer(s)",);

    /* Basics */
    mesh.setPrimitive(meshData.primitive());

    /* Vertex data */
    GL::Buffer verticesRef = GL::Buffer::wrap(vertices.id(), GL::Buffer::TargetHint::Array);
    GL::Buffer indicesRef = GL::Buffer::wrap(indices.id(), GL::Buffer::TargetHint::Array);
    for(UnsignedInt i = 0; i != meshData.attributeCount(); ++i) {
        Containers::Optional<GL::DynamicAttribute> attribute;

        /* Ignore implementation-specific formats because GL needs three
           separate values to describe them so there's no way to put them in a
           single 32-bit value :( */
        const VertexFormat format = meshData.attributeFormat(i);
        if(isVertexFormatImplementationSpecific(format)) {
            Warning{} << "upload: ignoring attribute" << meshData.attributeName(i) << "with an implementation-specific format" << reinterpret_cast<void*>(vertexFormatUnwrap(format));
            continue;
        }

        switch(meshData.attributeName(i)) {
            case Trade::MeshAttribute::Position:
                /* Pick 3D position always, the format will properly reduce it
                   to a 2-component version if needed */
                attribute.emplace(Shaders::Generic3D::Position{}, format);
                break;
            case Trade::MeshAttribute::TextureCoordinates:
                /** @todo have Generic2D derived from Generic that has all
                    attribute definitions common for 2D and 3D */
                attribute.emplace(Shaders::Generic2D::TextureCoordinates{}, format);
                break;
            case Trade::MeshAttribute::Color:
                /** @todo have Generic2D derived from Generic that has all
                    attribute definitions common for 2D and 3D */
                /* Pick Color4 always, the format will properly reduce it to a
                   3-component version if needed */
                attribute.emplace(Shaders::Generic2D::Color4{}, format);
                break;
            case Trade::MeshAttribute::Tangent:
                /* Pick Tangent4 always, the format will properly reduce it to
                   a 3-component version if needed */
                attribute.emplace(Shaders::Generic3D::Tangent4{}, format);
                break;
            case Trade::MeshAttribute::Bitangent:
                attribute.emplace(Shaders::Generic3D::Bitangent{}, format);
                break;
            case Trade::MeshAttribute::Normal:
                attribute.emplace(Shaders::Generic3D::Normal{}, format);
                break;

                /* So it doesn't yell that we didn't handle a known attribute */
            case Trade::MeshAttribute::Custom: break; /* LCOV_EXCL_LINE */
        }

        if(!attribute) {
            Warning{} << "upload: ignoring unknown attribute" << meshData.attributeName(i);
            continue;
        }

        mesh.addVertexBuffer(verticesRef,
                             meshData.attributeOffset(i), meshData.attributeStride(i),
                             *attribute);

    }

    if(meshData.isIndexed()) {
        mesh.setIndexBuffer(indicesRef, 0, meshData.indexType())
                .setCount(meshData.indexCount());
    } else mesh.setCount(meshData.vertexCount());

}


Trade::MeshData preprocess(Trade::MeshData& meshData, CompileFlag flags) {

    CORRADE_INTERNAL_ASSERT(meshData.primitive() == MeshPrimitive::Triangles);
    CORRADE_INTERNAL_ASSERT(meshData.attributeCount(Trade::MeshAttribute::Position));
    /* This could fire if we have 2D positions or for packed formats */
    CORRADE_INTERNAL_ASSERT(meshData.attributeFormat(Trade::MeshAttribute::Position) == VertexFormat::Vector3);

    /* If the data already have a normal array, reuse its location,
       otherwise mix in an extra one */
    Trade::MeshAttributeData normalAttribute;
    Containers::ArrayView<const Trade::MeshAttributeData> extra;
    if (!meshData.hasAttribute(Trade::MeshAttribute::Normal)) {
        normalAttribute = Trade::MeshAttributeData{
                Trade::MeshAttribute::Normal, VertexFormat::Vector3,
                nullptr};
        extra = {&normalAttribute, 1};
        /* If we reuse a normal location, expect correct type */
    } else
        CORRADE_INTERNAL_ASSERT(meshData.attributeFormat(Trade::MeshAttribute::Normal) == VertexFormat::Vector3);

    /* If we want flat normals, we need to first duplicate everything using
       the index buffer. Otherwise just interleave the potential extra
       normal attribute in. */
    Trade::MeshData generated{MeshPrimitive::Points, 0};
    if (flags & CompileFlag::GenerateFlatNormals && meshData.isIndexed())
        generated = MeshTools::duplicate(meshData, extra);
    else
        generated = MeshTools::interleave(meshData, extra);

    /* Generate the normals. If we don't have the index buffer, we can only
       generate flat ones. */
    if (flags & CompileFlag::GenerateFlatNormals || !meshData.isIndexed())
        MeshTools::generateFlatNormalsInto(
                generated.attribute<Vector3>(Trade::MeshAttribute::Position),
                generated.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));
    else
        MeshTools::generateSmoothNormalsInto(generated.indices(),
                                             generated.attribute<Vector3>(Trade::MeshAttribute::Position),
                                             generated.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));

    return generated;
}
