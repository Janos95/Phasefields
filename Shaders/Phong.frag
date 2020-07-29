//shamelessly stolen from https://github.com/mosra/magnum

#if defined(OBJECT_ID) && !defined(GL_ES) && !defined(NEW_GLSL)
#extension GL_EXT_gpu_shader4: require
#endif

#ifndef NEW_GLSL
#define in varying
#define fragmentColor gl_FragColor
#define texture texture2D
#endif

#ifndef RUNTIME_CONST
#define const
#endif

#ifdef AMBIENT_TEXTURE
#ifdef EXPLICIT_TEXTURE_LAYER
layout(binding = 0)
#endif
uniform lowp sampler2D ambientTexture;
#endif

#ifdef EXPLICIT_UNIFORM_LOCATION
layout(location = 4)
#endif
uniform lowp vec4 ambientColor
    #ifndef GL_ES
    #ifndef AMBIENT_TEXTURE
    = vec4(0.0)
    #else
    = vec4(1.0)
    #endif
    #endif
    ;

#if LIGHT_COUNT
#ifdef DIFFUSE_TEXTURE
#ifdef EXPLICIT_TEXTURE_LAYER
layout(binding = 1)
#endif
uniform lowp sampler2D diffuseTexture;
#endif

#ifdef EXPLICIT_UNIFORM_LOCATION
layout(location = 5)
#endif
uniform lowp vec4 diffuseColor
    #ifndef GL_ES
    = vec4(1.0)
    #endif
    ;

#ifdef SPECULAR_TEXTURE
#ifdef EXPLICIT_TEXTURE_LAYER
layout(binding = 2)
#endif
uniform lowp sampler2D specularTexture;
#endif

#ifdef NORMAL_TEXTURE
#ifdef EXPLICIT_TEXTURE_LAYER
layout(binding = 3)
#endif
uniform lowp sampler2D normalTexture;
#endif

#ifdef EXPLICIT_UNIFORM_LOCATION
layout(location = 6)
#endif
uniform lowp vec4 specularColor
    #ifndef GL_ES
    = vec4(1.0)
    #endif
    ;

#ifdef EXPLICIT_UNIFORM_LOCATION
layout(location = 7)
#endif
uniform mediump float shininess
    #ifndef GL_ES
    = 80.0
    #endif
    ;
#endif

#if LIGHT_COUNT
/* Needs to be last because it uses locations 10 + LIGHT_COUNT to
   10 + 2*LIGHT_COUNT - 1. Location 10 is lightPositions. Also it can't be
   specified as 10 + LIGHT_COUNT because that requires ARB_enhanced_layouts. */
#ifdef EXPLICIT_UNIFORM_LOCATION
layout(location = LIGHT_COLORS_LOCATION) /* I fear this will blow up some drivers */
#endif
uniform lowp vec4 lightColors[LIGHT_COUNT]
    #ifndef GL_ES
    = vec4[](LIGHT_COLOR_INITIALIZER)
    #endif
    ;
#endif

#if LIGHT_COUNT
in mediump vec3 transformedNormal;
#ifdef NORMAL_TEXTURE
in mediump vec3 transformedTangent;
#endif
in highp vec3 lightDirections[LIGHT_COUNT];
in highp vec3 cameraDirection;
#endif

#if defined(AMBIENT_TEXTURE) || defined(DIFFUSE_TEXTURE) || defined(SPECULAR_TEXTURE) || defined(NORMAL_TEXTURE)
in mediump vec2 interpolatedTextureCoordinates;
#endif

#ifdef VERTEX_COLOR
in lowp vec4 interpolatedVertexColor;
#endif

#ifdef NEW_GLSL
#ifdef EXPLICIT_ATTRIB_LOCATION
layout(location = COLOR_OUTPUT_ATTRIBUTE_LOCATION)
#endif
out lowp vec4 fragmentColor;
#endif

void main() {
    lowp const vec4 finalAmbientColor =
        #ifdef AMBIENT_TEXTURE
        texture(ambientTexture, interpolatedTextureCoordinates)*
        #endif
        #ifdef VERTEX_COLOR
        interpolatedVertexColor*
        #endif
        ambientColor;
    #if LIGHT_COUNT
    lowp const vec4 finalDiffuseColor =
        #ifdef DIFFUSE_TEXTURE
        texture(diffuseTexture, interpolatedTextureCoordinates)*
        #endif
        #ifdef VERTEX_COLOR
        interpolatedVertexColor*
        #endif
        diffuseColor;
    lowp const vec4 finalSpecularColor =
        #ifdef SPECULAR_TEXTURE
        texture(specularTexture, interpolatedTextureCoordinates)*
        #endif
        specularColor;
    #endif

    /* Ambient color */
    fragmentColor = finalAmbientColor;

    #if LIGHT_COUNT
    /* Normal */
    mediump vec3 normalizedTransformedNormal = normalize(transformedNormal);
    #ifdef NORMAL_TEXTURE
    mediump vec3 normalizedTransformedTangent = normalize(transformedTangent);
    mediump mat3 tbn = mat3(
        normalizedTransformedTangent,
        normalize(cross(normalizedTransformedNormal,
                        normalizedTransformedTangent)),
        normalizedTransformedNormal
    );
    normalizedTransformedNormal = tbn*(texture(normalTexture, interpolatedTextureCoordinates).rgb*2.0 - vec3(1.0));
    #endif

    /* Add diffuse color for each light */
    for(int i = 0; i < LIGHT_COUNT; ++i) {
        highp vec3 normalizedLightDirection = normalize(lightDirections[i]);
        lowp float intensity = max(0.0, dot(normalizedTransformedNormal, normalizedLightDirection));
        fragmentColor += vec4(finalDiffuseColor.rgb*lightColors[i].rgb*intensity, lightColors[i].a*finalDiffuseColor.a/float(LIGHT_COUNT));

        /* Add specular color, if needed */
        if(intensity > 0.001) {
            highp vec3 reflection = reflect(-normalizedLightDirection, normalizedTransformedNormal);
            mediump float specularity = clamp(pow(max(0.0, dot(normalize(cameraDirection), reflection)), shininess), 0.0, 1.0);
            fragmentColor += vec4(finalSpecularColor.rgb*specularity, finalSpecularColor.a);
        }
    }
    #endif
}
