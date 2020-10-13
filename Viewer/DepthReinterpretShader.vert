out mediump vec2 textureCoordinates;

void main() {
    /*  -1   0   1       3
       1 0-------+-------2
         |       |     /
       0 |       |   /
         |       | /
      -1 +-------+
         |     /
         |   /
         | /
      -3 1              */
    gl_Position = vec4(gl_VertexID == 2 ?  3.0 : -1.0,
    gl_VertexID == 1 ? -3.0 :  1.0, 0.0, 1.0);

    /*   0  0.5  1       2
       1 0-------+-------2
         |       |     /
     0.5 |       |   /
         |       | /
       0 +-------+
         |     /
         |   /
         | /
      -1 1              */
    textureCoordinates = vec2(gl_VertexID == 2 ?  2.0 : 0.0,
    gl_VertexID == 1 ? -1.0 : 1.0);
}