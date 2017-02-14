#version 120


const float steps =250.0;
varying vec4 v_pos;
uniform sampler2D u_front_texcoord_map;
uniform sampler2D u_back_texcoord_map;
uniform sampler2D u_volume_data;
uniform sampler2D u_environment_depth;
uniform sampler2D u_transfer_function;

uniform int u_render_as_image;
uniform int u_with_env_depth;


uniform float u_number_of_slices;
uniform float u_slices_over_x;
uniform float u_slices_over_y;



vec4 texture3D(sampler2D tex,vec3 volpos)
{

 		float s1,s2;
        float dx1,dy1;
        float dx2,dy2;

        vec2 texpos1,texpos2;

        s1 = floor(volpos.z*(u_number_of_slices-1));
        s2 = ceil(volpos.z*(u_number_of_slices-1));
        //s2 = s1+1.0;
        //if(s2>u_number_of_slices)
        //	s2=u_number_of_slices;

        dx1 = fract(s1/u_slices_over_x);
        dy1 = floor(s1/u_slices_over_x)/u_slices_over_y;

        dx2 = fract(s2/u_slices_over_x);
        dy2 = floor(s2/u_slices_over_x)/u_slices_over_y;
        
			//	if(volpos.x==1.0)
			//	volpos.x=0.999;
			//	if(volpos.y==1.0)
			//	volpos.y=0.999;
				
			//	if(volpos.x<=0.0)
			//	volpos.x=0.001;
			//	if(volpos.y<=0.0)
			//	volpos.y=0.001;
				
        texpos1.x = dx1+(volpos.x/u_slices_over_x);
        texpos1.y = dy1+(volpos.y/u_slices_over_y);
        
				//if(texpos1.x==1.0)
				//texpos1.x=0.999;
				//if(texpos1.y==1.0)
				//texpos1.y=0.999;
				
        texpos2.x = dx2+(volpos.x/u_slices_over_x);
        texpos2.y = dy2+(volpos.y/u_slices_over_y);
				//if(texpos2.x==1.0)
			//	texpos2.x=0.999;
			//	if(texpos2.y==1.0)
			//	texpos2.y=0.999;
        vec4 p1=texture2D(tex,texpos1);
        vec4 p2=texture2D(tex,texpos2);

       
		return mix( p1, p2, (volpos.z*u_number_of_slices)-s1);
}


 void main() {


    vec2 screenTexCoord = v_pos.xy/v_pos.w;
    screenTexCoord.x = 0.5*screenTexCoord.x + 0.5;
    screenTexCoord.y = 0.5*screenTexCoord.y + 0.5;



    vec4 frontTexCoord = texture2D(u_front_texcoord_map,screenTexCoord);
	vec4 backTexCoord = texture2D(u_back_texcoord_map,screenTexCoord);
	//gl_FragColor=vec4(backTexCoord.xyz,1.0);
		//gl_FragColor=vec4(texture2D(u_volume_data,screenTexCoord).xyz,1.0);
	//return;
	float minimal=0.001;

	if(u_render_as_image==1)
	{	
		//gl_FragColor=vec4(texture2D(u_volume_data,screenTexCoord).xyz,1.0);
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
		 
        
// sample for rendeirng
     float depthStep=(backDepth-frontDepth)/steps;

     float currentDepth=frontDepth;
     

     vec3 dir = backTexCoord.xyz - currentPos.xyz;
     vec3 Step = dir/steps;
 

     vec4 accum = vec4(0, 0, 0, 0);
	 float accAlpha = 0.0;


     for(float i = 0.0; i<=steps; i+=1.0)
     {

	vec4 volume_value=texture3D(u_volume_data,currentPos.xyz);
     
	
	  vec4 srcValue=volume_value;
	  
	  //srcValue.a= texture2D(u_transfer_function,vec2(value,0.1)).r;
	 //  srcValue.a=1.0;
	  //srcValue.a=srcValue.r;
	 //if(srcValue.r<0.1)
	 //srcValue.a=0;
	 
	  
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