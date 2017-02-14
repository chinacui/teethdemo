#version 120
varying vec3 v_texCoord;

void main()
{

  gl_FragColor.rgb=v_texCoord;
  gl_FragColor.a=gl_FragCoord.z/gl_FragCoord.w;

}