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
} ubo;

layout(location = 0) out vec4 out_frag_color;

void main() {
    vec3 ray_target = in_view_center + vec3(in_corner, 0.0);
    vec3 ray_dir = normalize(ray_target);

    vec3 center = in_view_center;
    float r = in_radius;

    float b = dot(ray_dir, center);
    float c = dot(center, center) - r * r;
    float discriminant = b * b - c;

    if (discriminant < 0.0) {
        if (in_halo_color.a > 0.0) {
            float dist_center = length(in_corner);
            float halo_radius = r * max(in_style.x, 1.0);
            if (dist_center <= halo_radius) {
                float halo_alpha = smoothstep(halo_radius, r * 0.95, dist_center);
                out_frag_color = vec4(in_halo_color.rgb, in_halo_color.a * halo_alpha * 0.85);
                vec4 clip_center = ubo.proj * vec4(center, 1.0);
                gl_FragDepth = clip_center.z / clip_center.w;
                return;
            }
        }
        discard;
    }

    float t = b - sqrt(discriminant);
    vec3 hit_point = ray_dir * t;
    vec3 normal = normalize(hit_point - center);

    vec3 light_dir_view = normalize((ubo.view * vec4(ubo.light_dir.xyz, 0.0)).xyz);
    vec3 view_dir = -ray_dir;
    vec3 halfway = normalize(light_dir_view + view_dir);

    float diff = max(dot(normal, light_dir_view), 0.0);
    float spec = pow(max(dot(normal, halfway), 0.0), 32.0);

    float ambient = 0.35;
    float dim = in_style.y;
    vec3 shaded = (ambient + (1.0 - ambient) * diff) * in_color.rgb * dim + 0.3 * spec * vec3(1.0);

    float dist_perp = sqrt(max(dot(center, center) - b * b, 0.0));
    float alpha = in_color.a * (1.0 - smoothstep(r * 0.985, r, dist_perp));

    out_frag_color = vec4(shaded, alpha);

    vec4 clip_hit = ubo.proj * vec4(hit_point, 1.0);
    gl_FragDepth = clip_hit.z / clip_hit.w;
}
