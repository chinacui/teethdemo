#version 120
varying highp vec4 fp;
varying highp vec3 fn;
varying vec3 v_color;


uniform highp vec4 u_light_pos;
uniform highp vec4 u_light_diff;
uniform highp vec4 u_light_spec;
uniform highp vec4 u_light_amb;
uniform float u_spec_power;
uniform int use_texture;
uniform int is_shinning;
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
	if (use_texture == 1)
		vcolor = texture2D(u_sampler_exture, v_texcoords);
	else
		vcolor=vec4(v_color.x,v_color.y,v_color.z,1.0);
	vec3 R = reflect(-L, N);
	vec4 diffuse;
	diffuse = max(abs(dot(N,L)),0) * u_light_diff*vcolor;
	vec4 specular = pow(max(dot(R,V), 0.0), u_spec_power) * u_light_spec;
		
	gl_FragColor = vcolor*u_light_amb + diffuse ;//+ specular;
	if(is_shinning==1)
	gl_FragColor=vec4(gl_FragColor.x*1.3,gl_FragColor.y*1.0,gl_FragColor.z*1.1,gl_FragColor.a);

	
}
