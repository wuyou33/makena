#version 330 core

layout(location = 0) in vec3 vertexLCS;
layout(location = 1) in vec3 vertexRGB;

out vec3 vertexWCS;
out vec3 fragmentRGB;

uniform mat4 MVP;
uniform mat4 M;

void main() {

    gl_Position            = MVP * vec4(vertexLCS, 1.0);

    vertexWCS              = (M * vec4(vertexLCS, 1.0)).xyz;

    fragmentRGB            = vertexRGB;

}

