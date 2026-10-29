#version 450

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec4 in_color;

layout(binding = 0) uniform SceneUbo {
    mat4 view;
    mat4 proj;
    vec4 camera_pos;
    vec4 light_dir;
    vec2 viewport_size;
    float fov;
    float padding;
} ubo;

layout(location = 0) out vec4 out_color;

void main() {
    gl_Position = ubo.proj * (ubo.view * vec4(in_position, 1.0));
    out_color = in_color;
}
