#version 450

// Quad corner attributes: x in [-1, 1], y in [0, 1]
layout(location = 0) in vec2 in_corner;

// Per-instance bond attributes
layout(location = 1) in vec3 in_start;
layout(location = 2) in vec3 in_end;
layout(location = 3) in float in_radius;
layout(location = 4) in vec4 in_color;

layout(binding = 0) uniform SceneUbo {
    mat4 view;
    mat4 proj;
    vec4 camera_pos;
    vec4 light_dir;
    vec2 viewport_size;
    float fov;
    float padding;
} ubo;

layout(location = 0) out vec3 out_view_start;
layout(location = 1) out vec3 out_view_end;
layout(location = 2) out float out_radius;
layout(location = 3) out vec4 out_color;
layout(location = 4) out vec3 out_fragment_view_pos;

void main() {
    vec3 view_a = (ubo.view * vec4(in_start, 1.0)).xyz;
    vec3 view_b = (ubo.view * vec4(in_end, 1.0)).xyz;

    vec3 axis = view_b - view_a;
    float len = length(axis);
    vec3 axis_dir = (len > 1e-4) ? (axis / len) : vec3(0.0, 1.0, 0.0);

    vec3 mid = 0.5 * (view_a + view_b);
    vec3 side = normalize(cross(axis_dir, mid));
    if (length(side) < 1e-4) {
        side = vec3(1.0, 0.0, 0.0);
    }

    float r_expand = in_radius * 1.25;
    vec3 pos_along_axis = mix(view_a - axis_dir * (r_expand * 0.2),
                              view_b + axis_dir * (r_expand * 0.2),
                              in_corner.y);
    vec3 final_view_pos = pos_along_axis + side * (in_corner.x * r_expand);

    gl_Position = ubo.proj * vec4(final_view_pos, 1.0);

    out_view_start = view_a;
    out_view_end = view_b;
    out_radius = in_radius;
    out_color = in_color;
    out_fragment_view_pos = final_view_pos;
}
