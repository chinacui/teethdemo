#version 120

varying vec3 v_color;



void main(void) 
{
 

	vec4 vcolor=vec4(v_color.x,v_color.y,v_color.z,1.0);;
	gl_FragColor = vcolor;

	
}
