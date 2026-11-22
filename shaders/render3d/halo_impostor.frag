#version 450

layout(location = 0) in vec2 in_corner;
layout(location = 1) in vec3 in_view_center;
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
    vec4 fog_params;
    vec4 bg_color;
} ubo;

layout(location = 0) out vec4 out_frag_color;

void main() {
    if (in_halo_color.a <= 0.0) {
        discard;
    }

    float r = in_radius;
    float halo_radius = r * max(in_style.x, 1.0);
    float dist_center = length(in_corner);
    if (dist_center > halo_radius) {
        discard;
    }

    // Smooth radial falloff from core boundary to halo edge
    float inner_edge = (in_style.y > 0.0) ? (r * 0.95) : (r * 0.25);
    float halo_alpha = smoothstep(halo_radius, inner_edge, dist_center);
    float out_alpha = in_halo_color.a * halo_alpha * 0.85;
    if (out_alpha <= 0.002) {
        discard;
    }

    out_frag_color = vec4(in_halo_color.rgb, out_alpha);

    vec4 clip_center = ubo.proj * vec4(in_view_center, 1.0);
    gl_FragDepth = clip_center.z / clip_center.w;
}
