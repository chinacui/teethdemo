#version 120
attribute highp vec4 a_pos;
attribute vec3 a_texcoords;
varying vec3 v_texCoord;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
void main() {
    gl_Position = mvp_matrix*a_pos;
    v_texCoord=a_texcoords;
 

}