#version 120
varying highp vec4 fp;
varying highp vec3 fn;
varying vec3 v_color;


uniform highp vec4 u_light_pos;
uniform highp vec4 u_light_diff;
uniform highp vec4 u_light_spec;
uniform highp vec4 u_light_amb;
uniform float u_spec_power;

uniform sampler2D u_sampler_exture;
varying vec2 v_texcoords;

void main(void) 
{
 
	vec3 L = u_light_pos.xyz - fp.xyz;
	vec3 V = -fp.xyz;
	vec3 N;
	if(fn == vec3(0.0,0.0,0.0))
	    N = vec3(0.0,0.0,0.0);
	else
		N = normalize(fn);

	L = normalize(L);
	V = normalize(V);
	vec4 vcolor;
	vcolor=vec4(vcolor.x,vcolor.y,vcolor.z,1.0);
	vec3 R = reflect(-L, N);
	vec4 diffuse;
	diffuse = max(abs(dot(N,L)),0) * u_light_diff*vcolor;
	vec4 specular = pow(max(dot(R,V), 0.0), u_spec_power) * u_light_spec;
		
	gl_FragColor = vcolor*u_light_amb + diffuse + specular;
		
	
}
