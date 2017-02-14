#version 120


const float steps =150.0;
varying vec4 v_pos;
uniform sampler2D u_front_texcoord_map;
uniform sampler2D u_back_texcoord_map;
uniform sampler3D u_volume_data;
uniform sampler2D u_environment_depth;
uniform sampler2D u_transfer_function;

uniform int u_render_as_image;
uniform int u_with_env_depth;



 void main() {


    vec2 screenTexCoord = v_pos.xy/v_pos.w;
    screenTexCoord.x = 0.5*screenTexCoord.x + 0.5;
    screenTexCoord.y = 0.5*screenTexCoord.y + 0.5;



    vec4 frontTexCoord = texture2D(u_front_texcoord_map,screenTexCoord);
	vec4 backTexCoord = texture2D(u_back_texcoord_map,screenTexCoord);
	float minimal=0.001;

	if(u_render_as_image==1)
	{	
		gl_FragColor=texture3D(u_volume_data,frontTexCoord.xyz);
		
	}
	else
	{
		
      vec4 currentPos = frontTexCoord;

      float frontDepth=frontTexCoord.a;
      float backDepth=backTexCoord.a;
	  float envDepth=999999999.0;

	 if(u_with_env_depth==1&&texture2D(u_environment_depth,screenTexCoord).r>0.0)
	 {
	  envDepth=texture2D(u_environment_depth,screenTexCoord).r;
	 }
		 
        

     float depthStep=(backDepth-frontDepth)/steps;

     float currentDepth=frontDepth;
     

     vec3 dir = backTexCoord.xyz - currentPos.xyz;
     vec3 Step = dir/steps;
 

     vec4 accum = vec4(0, 0, 0, 0);
	 float accAlpha = 0.0;


     for(float i = 0.0; i<=steps; i+=1.0)
     {

   
      float value=texture3D(u_volume_data,currentPos.xyz).r;
	  vec4 srcValue=vec4(value);
	  srcValue.a= texture2D(u_transfer_function,vec2(value,0.1)).r;
	   

	  srcValue.rgb *= srcValue.a;
      
	  accAlpha += (1.0 - pow((1.0 - srcValue.a), depthStep));
      accum += ((1.0 - accum.a) * srcValue);

      currentPos.xyz += Step;

      currentDepth+=depthStep;

     if(currentPos.x< 0.0 || currentPos.y < 0.0 || currentPos.z<0.0||currentPos.x> 1.0 || currentPos.y >1.0 || currentPos.z>1.0 ||accAlpha>=1.0||currentDepth>=envDepth)
            break;

     }
          
     gl_FragColor=accum;
			

  }
}