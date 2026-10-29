#version 450

layout(location = 0) in vec3 in_view_start;
layout(location = 1) in vec3 in_view_end;
layout(location = 2) in float in_radius;
layout(location = 3) in vec4 in_color;
layout(location = 4) in vec3 in_fragment_view_pos;

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
    vec3 ray_dir = normalize(in_fragment_view_pos);

    vec3 a = in_view_start;
    vec3 b = in_view_end;
    vec3 axis = b - a;
    float axis_len = length(axis);
    if (axis_len < 1e-4) {
        discard;
    }
    vec3 u = axis / axis_len;

    vec3 m = -a; // ray origin is (0,0,0) in view space

    vec3 d_perp = ray_dir - dot(ray_dir, u) * u;
    vec3 m_perp = m - dot(m, u) * u;

    float r = in_radius;
    float quad_a = dot(d_perp, d_perp);
    if (quad_a < 1e-6) {
        discard;
    }

    float quad_b = dot(d_perp, m_perp);
    float quad_c = dot(m_perp, m_perp) - r * r;
    float discriminant = quad_b * quad_b - quad_a * quad_c;

    if (discriminant < 0.0) {
        discard;
    }

    float t = (-quad_b - sqrt(discriminant)) / quad_a;
    if (t < 0.0) {
        discard;
    }

    vec3 hit_point = ray_dir * t;
    float h = dot(hit_point - a, u);

    if (h < 0.0 || h > axis_len) {
        discard;
    }

    vec3 normal = normalize((hit_point - a) - h * u);

    vec3 light_dir_view = normalize((ubo.view * vec4(ubo.light_dir.xyz, 0.0)).xyz);
    vec3 view_dir = -ray_dir;
    vec3 halfway = normalize(light_dir_view + view_dir);

    float diff = max(dot(normal, light_dir_view), 0.0);
    float spec = pow(max(dot(normal, halfway), 0.0), 24.0);

    float ambient = 0.35;
    vec3 shaded = (ambient + (1.0 - ambient) * diff) * in_color.rgb + 0.25 * spec * vec3(1.0);

    float dist_axis = length((hit_point - a) - h * u);
    float alpha = in_color.a * (1.0 - smoothstep(r * 0.95, r, dist_axis));

    out_frag_color = vec4(shaded, alpha);

    vec4 clip_hit = ubo.proj * vec4(hit_point, 1.0);
    gl_FragDepth = clip_hit.z / clip_hit.w;
}
