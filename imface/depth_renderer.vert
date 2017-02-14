#version 120
varying vec4 v_pos;
void main() {
  gl_Position = projectionMatrix *modelViewMatrix*vec4(position,1.0);
  v_pos=gl_Position;

}