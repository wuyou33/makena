#version 330 core

in vec3 vertexWCS;
in vec3 fragmentRGB;

out vec3 color;

uniform vec3 lightWCS;

void main() {

    vec3  lightColor       = vec3(1,1,1);
    float lightPower       = 50.0f;
     
    vec3  diffuseColor     = fragmentRGB;
    vec3  ambientColor     = vec3(0.3, 0.3, 0.3) * diffuseColor;

    float distLightVertex  = length( lightWCS - vertexWCS );

    vec3  coeff            = (lightColor * lightPower) /
                                             (distLightVertex*distLightVertex);

    color                  = ambientColor + diffuseColor * coeff;

}
