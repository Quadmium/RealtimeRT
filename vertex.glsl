#version 450 core

in vec2 in_position;
out vec2 position;


void main() {
    gl_Position = vec4(in_position.x, in_position.y, 0.0, 1.0);
    position = in_position;
}