#version 120
attribute highp vec4 a_pos;
attribute highp vec3 a_color;
varying vec3 v_color;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;

void main(void)
{

	gl_Position = mvp_matrix * a_pos;
	v_color=a_color;
 }