uniform highp sampler2D depthTexture;

in mediump vec2 textureCoordinates;

out highp vec3 reinterpretDepth;

void main() {
    /* Convert from range 0.0 - 1.0 to 0 - 0xffffff (we have a 24bit depth and
       floats have 24bit mantissa so this should preserve everything), then
       separate that into three 8bit values and then unpack each 8bit value
       back to 0.0 - 1.0 again in order to make the WebGL RGBA8 pipeline happy.
       All the fancy packing algos from https://stackoverflow.com/q/9882716 and
       elsewhere were not treating depth = 1.0 correctly, so I'm doing my own
       thing here. */
    highp float depth = texture(depthTexture, textureCoordinates).r;
    highp uint depthI = uint(depth*float(0xffffffu));
    highp uvec3 depthIV = uvec3((depthI >> 16), (depthI >> 8), depthI) & 0xffu;
    reinterpretDepth = vec3(depthIV)/255.0;
}