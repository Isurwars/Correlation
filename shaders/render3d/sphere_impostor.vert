#version 450

// Quad corner attributes (-1.0 to 1.0)
layout(location = 0) in vec2 in_corner;

// Per-instance atom attributes
layout(location = 1) in vec3 in_center;
layout(location = 2) in float in_radius;
layout(location = 3) in vec4 in_color;
layout(location = 4) in vec4 in_halo_color;
layout(location = 5) in vec2 in_style; // x = halo_radius_ratio, y = dim_factor

layout(binding = 0) uniform SceneUbo {
    mat4 view;
    mat4 proj;
    vec4 camera_pos;
    vec4 light_dir;
    vec2 viewport_size;
    float fov;
    float padding;
} ubo;

layout(location = 0) out vec2 out_corner;
layout(location = 1) out vec3 out_view_center;
layout(location = 2) out float out_radius;
layout(location = 3) out vec4 out_color;
layout(location = 4) out vec4 out_halo_color;
layout(location = 5) out vec2 out_style;

void main() {
    vec4 view_center = ubo.view * vec4(in_center, 1.0);
    float halo_multiplier = (in_halo_color.a > 0.0) ? max(in_style.x, 1.0) : 1.0;
    float bounding_radius = in_radius * halo_multiplier * 1.15;

    vec3 view_pos = view_center.xyz + vec3(in_corner * bounding_radius, 0.0);
    gl_Position = ubo.proj * vec4(view_pos, 1.0);

    out_corner = in_corner * bounding_radius;
    out_view_center = view_center.xyz;
    out_radius = in_radius;
    out_color = in_color;
    out_halo_color = in_halo_color;
    out_style = in_style;
}
