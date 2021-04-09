#version 450 core
in vec2 position;
out vec4 frag_color;

uniform int scene;
uniform float g_seed;
uniform float blend;

uniform vec3 camera_forward;
uniform vec3 camera_up;
uniform vec3 camera_right;
uniform vec3 camera_pos;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 1000.5,seed += 1000.5)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 1000.5,seed += 1000.5)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 1000.5, seed += 1000.5)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

vec3 random_in_unit_sphere(inout float seed)
{
  vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
  float phi = h.y;
  float r = pow(h.z, 1.0/3.0);
  return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 random_unit_vector(inout float seed) {
  return normalize(random_in_unit_sphere(seed));
}

bool near_zero(vec3 e) {
  const float s = 0.0000001;
  return abs(e.x) < s && abs(e.y) < s && abs(e.z) < s;
}

float reflectance(float c, float refract_ratio) {
  float r0 = (1 - refract_ratio) / (1 + refract_ratio);
  r0 = r0 * r0;
  return r0 + (1 - r0) * pow((1 - c), 5);
}

// Pack two vec3's into one (losing half accuracy)
vec3 pack_two(vec3 a, vec3 b) {
  uint ux = packHalf2x16(vec2(a.x, b.x));
  uint uy = packHalf2x16(vec2(a.y, b.y));
  uint uz = packHalf2x16(vec2(a.z, b.z));
  return vec3(uintBitsToFloat(ux), uintBitsToFloat(uy), uintBitsToFloat(uz));
}

void unpack_two(vec3 p, out vec3 a, out vec3 b) {
  uint ux = floatBitsToUint(p.x);
  uint uy = floatBitsToUint(p.y);
  uint uz = floatBitsToUint(p.z);
  vec2 unpack_x = unpackHalf2x16(ux);
  vec2 unpack_y = unpackHalf2x16(uy);
  vec2 unpack_z = unpackHalf2x16(uz);
  a = vec3(unpack_x.x, unpack_y.x, unpack_z.x);
  b = vec3(unpack_x.y, unpack_y.y, unpack_z.y);
}

struct Triangle
{
	vec3 vertex0;
	vec3 vertex1;
	vec3 vertex2;
};

struct Ray
{
	vec3 origin;
	vec3 direction;
};
vec3 Ray_at(Ray self, float t) {
  return self.origin + t * self.direction;
}

struct Material
{
  int type;
  vec3 albedo;
};
const int Material_Lambertian = 0;
const int Material_Metal = 1;
const int Material_Dielectric = 2;
const int Material_Checkers = 3;
const int Material_Light = 4;

struct HitResult
{
  bool hit;
  vec3 point;
  vec3 normal;
  float t;
  bool front;
  Material material;
};
void HitResult_set_normal(inout HitResult h, Ray r, vec3 outward_normal) {
  h.front = dot(r.direction, outward_normal) < 0;
  h.normal = h.front ? outward_normal : -outward_normal;
}

struct Sphere
{
  vec3 center;
  float radius;
  Material material;
};
HitResult Sphere_hit(Sphere self, Ray r, float t_min, float t_max) {
  HitResult res;
  res.hit = false;
  
  vec3 d_unit = normalize(r.direction);
  vec3 f = r.origin - self.center;
  float a = dot(r.direction, r.direction);
  float b = 2 * dot(f, r.direction);
  float c = dot(f, f) - self.radius * self.radius;

  vec3 tmp = (f - dot(f, d_unit) * d_unit);
  float b2_minus_4ac = 4 * a * (self.radius * self.radius - dot(tmp, tmp));

  // Return early for invaid determinant
  if (b2_minus_4ac < 0) {
    return res;
  }

  float q = -0.5 * (b + (b >= 0 ? 1 : -1) * sqrt(b2_minus_4ac));

  // Calculate two solutions
  float t0 = c / q;
  float t1 = q / a;

  bool t0_valid = t0 >= t_min && t0 <= t_max;
  bool t1_valid = t1 >= t_min && t1 <= t_max;

  if (t0_valid && t1_valid) {
    res.t = min(t0, t1);
  } else if (t0_valid) {
    res.t = t0;
  } else if (t1_valid) {
    res.t = t1;
  } else {
    return res;
  }

  res.hit = true;
  res.point = Ray_at(r, res.t);
  HitResult_set_normal(res, r, (res.point - self.center) / self.radius);
  res.material = self.material;

  return res;
}

struct ScatterResult {
  int effect;
  vec3 attenuation;
  Ray scattered;
};
const int ScatterResult_Effect_Absorbed = 0;
const int ScatterResult_Effect_Scattered = 1;
const int ScatterResult_Effect_Emitted = 2;
ScatterResult Material_scatter(Ray r, HitResult hit, inout float seed) {
  Material material = hit.material;
  if (material.type == Material_Lambertian) {
    vec3 scatter_dir = hit.normal + random_unit_vector(seed);

    if (near_zero(scatter_dir)) {
      scatter_dir = hit.normal;
    }

    return ScatterResult(ScatterResult_Effect_Scattered, material.albedo, Ray(hit.point, normalize(scatter_dir)));
  } else if (material.type == Material_Metal) {
    vec3 reflected = reflect(normalize(r.direction), hit.normal);
    Ray scattered = Ray(hit.point, reflected);
    return ScatterResult(dot(scattered.direction, hit.normal) > 0 ? ScatterResult_Effect_Scattered : ScatterResult_Effect_Absorbed, material.albedo, scattered);
  } else if (material.type == Material_Dielectric) {
    float ir = 1.5;
    float refract_ratio = hit.front ? (1.0 / ir) : ir;

    vec3 unit_dir = normalize(r.direction);
    float c = min(dot(-unit_dir, hit.normal), 1.0);
    float s = sqrt(1.0 - c * c);

    bool no_refract = refract_ratio * s > 1.0;
    vec3 dir;

    if (no_refract || reflectance(c, refract_ratio) > hash1(seed)) {
      dir = reflect(unit_dir, hit.normal);
    } else {
      dir = refract(unit_dir, hit.normal, refract_ratio);
    }

    return ScatterResult(ScatterResult_Effect_Scattered, material.albedo, Ray(hit.point, dir));
  } else if (material.type == Material_Checkers) {
    vec3 albedo1;
    vec3 albedo2;
    unpack_two(material.albedo, albedo1, albedo2);

    vec3 scatter_dir = hit.normal + random_unit_vector(seed);

    if (near_zero(scatter_dir)) {
      scatter_dir = hit.normal;
    }

    float scale = 2;
    int x = int(10000 + hit.point.x * scale);
    int z = int(10000 + hit.point.z * scale);

    vec3 final_albedo = albedo1;

    if (x % 2 == z % 2) {
      final_albedo = albedo2;
    }

    return ScatterResult(ScatterResult_Effect_Scattered, final_albedo, Ray(hit.point, normalize(scatter_dir)));
  } else if (material.type == Material_Light) {
    return ScatterResult(ScatterResult_Effect_Emitted, material.albedo, Ray(vec3(0, 0, 0), vec3(0, 0, 0)));
  }
}

void main()
{
  float seed = g_seed + baseHash(floatBitsToUint(position));

  float focal = 1.0;

  const int max_spheres = 10;
  int n_spheres = 0;
  Sphere all_spheres[max_spheres];
  int n_rays = 4;
  bool have_sky = true;

  if (scene == 1) {
    have_sky = true;
    n_rays = 4;

    Sphere s;
    s.center = vec3(0, 0, -2);
    s.radius = 0.5;
    s.material = Material(Material_Lambertian, vec3(0.5, 0.5, 0.5));

    Sphere ground;
    ground.center = vec3(0.0, -100.5, -1.0);
    ground.radius = 100;
    ground.material = Material(Material_Lambertian, vec3(0.5, 0.5, 0.5));

    n_spheres = 2;
    all_spheres[0] = s;
    all_spheres[1] = ground;
  } else if (scene == 2) {
    have_sky = true;
    n_rays = 4;

    Sphere ground;
    ground.center = vec3(0.0, -100.5, -1.0);
    ground.radius = 100;
    ground.material = Material(Material_Checkers, pack_two(vec3(1, 0, 0), vec3(1, 1, 0)));

    Sphere s1;
    s1.center = vec3(1, 0, -2);
    s1.radius = 0.5;
    s1.material = Material(Material_Metal, vec3(0.5, 0.5, 0.5));
    Sphere s2;
    s2.center = vec3(-1, 0, -2);
    s2.radius = 0.5;
    s2.material = Material(Material_Lambertian, vec3(0, 0, 1));
    Sphere s3;
    s3.center = vec3(0, 0.5, -1);
    s3.radius = 0.5;
    s3.material = Material(Material_Dielectric, vec3(1, 1, 1));
    Sphere s4;
    s4.center = vec3(0, 0.5, -1);
    s4.radius = -0.45;
    s4.material = Material(Material_Dielectric, vec3(1, 1, 1));

    n_spheres = 5;
    all_spheres[0] = ground;
    all_spheres[1] = s1;
    all_spheres[2] = s2;
    all_spheres[3] = s3;
    all_spheres[4] = s4;
  } else if (scene == 3) {
    have_sky = false;
    n_rays = 20;

    Sphere s;
    s.center = vec3(0, 0.5, -4);
    s.radius = 1;
    s.material = Material(Material_Light, vec3(1, 1, 1));

    Sphere s2;
    s2.center = vec3(-1, 0, -2);
    s2.radius = 0.5;
    s2.material = Material(Material_Lambertian, vec3(1, 0, 0));
    Sphere s3;
    s3.center = vec3(1, 0, -2);
    s3.radius = 0.5;
    s3.material = Material(Material_Lambertian, vec3(0, 0, 1));

    Sphere ground;
    ground.center = vec3(0.0, -100.5, -1.0);
    ground.radius = 100;
    ground.material = Material(Material_Lambertian, vec3(1, 1, 1));

    n_spheres = 4;
    all_spheres[0] = s;
    all_spheres[1] = ground;
    all_spheres[2] = s2;
    all_spheres[3] = s3;
  } else if (scene == 4) {
    have_sky = true;
    n_rays = 4;

    Sphere ground;
    ground.center = vec3(0.0, -100.5, -1.0);
    ground.radius = 100;
    ground.material = Material(Material_Metal, vec3(0.5, 0.5, 0.5));

    Sphere s1;
    s1.center = vec3(0.5, -0.325, -3);
    s1.radius = 0.2;
    s1.material = Material(Material_Dielectric, vec3(1, 1, 1));
    Sphere s2;
    s2.center = vec3(0, -0.325, -3);
    s2.radius = 0.2;
    s2.material = Material(Material_Dielectric, vec3(1, 1, 1));
    Sphere s3;
    s3.center = vec3(-0.5, -0.325, -3);
    s3.radius = 0.2;
    s3.material = Material(Material_Dielectric, vec3(1, 1, 1));
    Sphere s4;
    s4.center = vec3(0, -0.22, -6);
    s4.radius = 0.4;
    s4.material = Material(Material_Dielectric, vec3(1, 1, 1));
    Sphere s5;
    s5.center = vec3(1, -0.14, -3.5);
    s5.radius = 0.4;
    s5.material = Material(Material_Checkers, pack_two(vec3(1, 1, 1), vec3(0, 0, 1)));
    Sphere s6;
    s6.center = vec3(-1, -0.14, -3.5);
    s6.radius = 0.4;
    s6.material = Material(Material_Lambertian, vec3(0, 0, 1));

    n_spheres = 7;
    all_spheres[0] = ground;
    all_spheres[1] = s1;
    all_spheres[2] = s2;
    all_spheres[3] = s3;
    all_spheres[4] = s4;
    all_spheres[5] = s5;
    all_spheres[6] = s6;
  }

  vec3 color_result = vec3(0, 0, 0);

  for (int ray_ctr = 0; ray_ctr < n_rays; ++ray_ctr) {
    vec2 rand_offset = hash2(seed) / 500;
    Ray ray;
    ray.origin = camera_pos;
    ray.direction = normalize(focal * camera_forward + (position.x/2+rand_offset.x) * camera_right + (position.y/2+rand_offset.y) * camera_up);

    vec3 current_ray_color = vec3(0, 0, 0);
    vec3 attenuation = vec3(1.0, 1.0, 1.0);
    int depth = 10;

    while (true) {
      if (depth <= 0) {
        break;
      }

      HitResult closest;
      closest.hit = false;
      for (int i = 0; i < n_spheres; ++i) {
        HitResult res = Sphere_hit(all_spheres[i], ray, 0.001, 10000);
        if (res.hit) {
          if (!closest.hit) {
            closest = res;
          } else if (res.t < closest.t) {
            closest = res;
          }
        }
      }

      if (!closest.hit) {
        if (have_sky) {
          vec3 unit_direction = normalize(ray.direction);
          float t = 0.5*(unit_direction.y + 1.0);
          vec3 color = (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
          color_result += attenuation * color;
        }
        break;
      }

      ScatterResult scatter = Material_scatter(ray, closest, seed);
      if (scatter.effect == ScatterResult_Effect_Absorbed) {
        break;
      }

      if (scatter.effect == ScatterResult_Effect_Emitted) {
        color_result += attenuation * scatter.attenuation;
        break;
      }
   
      attenuation *= scatter.attenuation;
      ray = scatter.scattered;
      depth -= 1;
    }

    color_result += current_ray_color;
  }

  color_result /= n_rays;

  frag_color = vec4(pow(color_result, vec3(0.5, 0.5, 0.5)), blend);
}