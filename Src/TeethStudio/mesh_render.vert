#version 120
attribute highp vec4 a_pos;
attribute highp vec3 a_normal;
attribute highp vec3 a_color;
varying vec3 v_color;
varying vec2 v_texcoords;
attribute vec2 a_texcoords;
uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
varying highp vec4 fp;
varying highp vec3 fn;
void main(void)
{
	fp = mv_matrix * a_pos;
	fn = mat3(mv_matrix)* a_normal;
	gl_Position = mvp_matrix * a_pos;
	v_texcoords=a_texcoords;
	v_color=a_color;
 }