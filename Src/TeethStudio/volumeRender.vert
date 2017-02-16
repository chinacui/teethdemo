#version 120

attribute highp vec4 a_pos;
varying vec4 v_pos;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;


void main() {

    gl_Position = mvp_matrix*a_pos;
    v_pos=gl_Position;

}