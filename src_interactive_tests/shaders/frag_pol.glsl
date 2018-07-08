#version 330 core

in vec3 vertexWCS;
in vec3 normalECS;
in vec3 vectorVertexToEyeECS;
in vec3 vectorVertexToLightECS;
in vec4 fragmentRGBA;

out vec4 color;

uniform mat4 MV;
uniform vec3 lightWCS;

void main() {

    vec3  lightColor       = vec3(1,1,1);
    float lightPower       = 100.0f;
     
    vec3  diffuseColor     = fragmentRGBA.xyz;
    vec3  ambientColor     = vec3(0.4, 0.4, 0.4) * diffuseColor;
    vec3  specularColor    = vec3(0.3, 0.3, 0.3);

    float distLightVertex  = length( lightWCS - vertexWCS );

    vec3  nECS             = normalize( normalECS );

    vec3  vertexToLightECS = normalize( vectorVertexToLightECS );

    vec3  vertexToEyeECS   = normalize( vectorVertexToEyeECS );

    float cosTheta         = clamp(dot(nECS, vertexToLightECS ), 0.0, 1.0);
     
    vec3  refECS           = reflect(-1.0* vertexToLightECS, nECS);

    float cosAlpha         = clamp(dot(vertexToEyeECS,refECS), 0.0,1.0 );
     
    vec3  coeff            = (lightColor * lightPower) /
                                           (distLightVertex*distLightVertex);
    vec3  rgb              = ambientColor  + 
                             diffuseColor  * coeff * cosTheta +
                             specularColor * coeff * pow( cosAlpha, 5.0);
    color                  = vec4(rgb, fragmentRGBA.w);
//    color                  = vec4(rgb, 1.0);

}
