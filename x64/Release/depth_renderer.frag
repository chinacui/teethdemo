#version 120
varying vec4 v_pos;
void main()
{


  gl_FragColor.r=gl_FragCoord.z/gl_FragCoord.w;

}