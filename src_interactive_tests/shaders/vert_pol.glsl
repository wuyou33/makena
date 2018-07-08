#version 330 core

layout(location = 0) in vec3 vertexLCS;
layout(location = 1) in vec4 vertexRGBA;
layout(location = 2) in vec3 normalLCS;

out vec3 vertexWCS;
out vec3 normalECS;
out vec3 vectorVertexToEyeECS;
out vec3 vectorVertexToLightECS;
out vec4 fragmentRGBA;

uniform mat4 MVP;
uniform mat4 M;
uniform mat4 V;
uniform vec3 lightWCS;

void main() {

    gl_Position            = MVP * vec4(vertexLCS, 1.0);

    vertexWCS              = (M * vec4(vertexLCS, 1.0)).xyz;

    vec3 vertexECS         = ( V * M * vec4(vertexLCS, 1.0) ).xyz;

    normalECS              = ( V * M * vec4(normalLCS,0.0) ).xyz;

    vectorVertexToEyeECS   = vec3(0,0,0) - vertexECS;

    vec3 lightECS          = ( V * vec4(lightWCS,1.0) ).xyz;

    vectorVertexToLightECS = lightECS - vertexECS;

    fragmentRGBA           = vertexRGBA;

}

