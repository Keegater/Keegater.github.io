#version 300 es
precision highp float;

// time for animation
uniform float u_time;
uniform vec2 u_resolution;
uniform int u_maxBounces;
uniform float u_ambientStrength;

// Camera params
uniform vec3 u_camPos;
uniform vec3 u_camForward;
uniform vec3 u_camRight;
uniform vec3 u_camUp;

// Element controls
uniform vec3  u_centerLightColor;
uniform float u_centerLightBrightness;
uniform float u_leftPaddleY;
uniform float u_rightPaddleY;
uniform vec2  u_ballPos;
uniform float u_ballRadius;
uniform sampler2D u_iceHeightMap;
uniform sampler2D u_iceNormalMap;
uniform sampler2D u_lavaTexture;
uniform float    u_iceBumpStrength;
uniform float    u_iceNormalBlend;

out vec4 fragColor;

// Material types
const int MATERIAL_DIFFUSE = 0;
const int MATERIAL_REFLECTIVE = 1;
const int MATERIAL_REFRACTIVE = 2;
const int MATERIAL_LIGHT = 3;

// Number of light spheres on wall
const int NUM_LIGHT_SPHERES = 7;

// Ray structure
struct Ray {
    vec3 origin;
    vec3 direction;
};

// Sphere structure
struct Sphere {
    vec3 center;
    float radius;
    vec3 color;
    int material;
    float reflectivity;
    float refractiveIndex;
};

// Box
struct Box {
    vec3 center;
    vec3 halfSize;
    vec3 color;
    int material;
    float reflectivity;
    bool isLava;
};

// Hit record structure
struct HitData {
    bool hit;
    float t;
    vec3 point;
    vec3 normal;
    vec3 color;
    int material;
    float reflectivity;
    float refractiveIndex;
    bool frontFace;
    bool isLava;
};

// Scene objects
Sphere spheres[1];                 // ball sphere
Sphere lightSpheres[NUM_LIGHT_SPHERES]; // light spheres
Box boxes[4];                      // left paddle, right paddle, top/bottom bounds
Box iceFloor;                      // reflective ice floor

// Initialize scene objects
void initScene(float time) {

    // Ball
    spheres[0] = Sphere(
        vec3(u_ballPos, 0.0),
        u_ballRadius,
        vec3(1.0, 0.3, 0.3),
        MATERIAL_DIFFUSE,
        0.0,
        1.0
    );

    // Left paddle (AI Opponent)
    boxes[0] = Box(
        vec3(-4.5, u_leftPaddleY, 0.0),
        vec3(0.05, 0.5, 0.2),
        vec3(0.8, 0.8, 0.8),
        MATERIAL_DIFFUSE,
        0.0,
        false
    );

    // Right paddle (player)
    boxes[1] = Box(
        vec3(4.5, u_rightPaddleY, 0.0),
        vec3(0.05, 0.5, 0.2),
        vec3(0.8, 0.8, 0.8),
        MATERIAL_DIFFUSE,
        0.0,
        false
    );

    // Bound boards top and bottom
    float boardTop = 3.25;
    float boardBottom = -3.25;
    float boardHalfWidth = 4.5;
    float borderHalfThickness = 0.25;

    // Set up light spheres
    vec3 lightCol = clamp(u_centerLightColor, 0.0, 1.0) * u_centerLightBrightness;
    float lightsY = boardTop - 0.01;
    float lightsZ = 1.5;
    for (int i = 0; i < NUM_LIGHT_SPHERES; i++) {
        float t = (float(i) + 0.5) / float(NUM_LIGHT_SPHERES);
        float x = mix(-boardHalfWidth - 1.5, boardHalfWidth + 1.5, t);
        lightSpheres[i] = Sphere(
            vec3(x, lightsY, lightsZ),
            0.08,
            lightCol,
            MATERIAL_LIGHT,
            0.0,
            1.0
        );
    }

    // Top border
    boxes[2] = Box(
        vec3(0.0, boardTop + borderHalfThickness, 0.0),
        vec3(boardHalfWidth + 2.0, borderHalfThickness, 4),
        vec3(0.01, 0.01, 0.01),
        MATERIAL_DIFFUSE,
        0.0,
        true
    );

    // Bottom border
    boxes[3] = Box(
        vec3(0.0, boardBottom - borderHalfThickness, 0.0),
        vec3(boardHalfWidth + 2.0, borderHalfThickness, 0.25),
        vec3(0.01, 0.01, 0.01),
        MATERIAL_DIFFUSE,
        0.0,
        false
    );

    // Ice floor
    iceFloor = Box(
        vec3(0.0, 0.0, -0.8),
        vec3(boardHalfWidth + 2.0, (boardTop - boardBottom) * 0.55, 0.6),
        vec3(0.1, 0.45, 0.8),
        MATERIAL_REFLECTIVE,
        0.6,
        false
    );

}

// Sphere intersection
bool intersectSphere(Ray ray, Sphere sphere, out float t) {
    vec3 oc = ray.origin - sphere.center;
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        return false;
    }

    float sqrtD = sqrt(discriminant);
    float t0 = (-b - sqrtD) / (2.0 * a);
    float t1 = (-b + sqrtD) / (2.0 * a);

    t = t0;
    if (t < 0.001) {
        t = t1;
        if (t < 0.001) {
            return false;
        }
    }
    return true;
}

// Box intersection
bool intersectBox(Ray ray, Box box, out float tOut) {
    // transform ray into box local coordinates centered at box.center
    vec3 rayOrigin = ray.origin - box.center;
    vec3 rayDirection = ray.direction;

    vec3 invD = 1.0 / rayDirection;

    vec3 t0s = (-box.halfSize - rayOrigin) * invD;
    vec3 t1s = ( box.halfSize - rayOrigin) * invD;

    vec3 tsmaller = min(t0s, t1s);
    vec3 tbigger  = max(t0s, t1s);

    float tmin = max(max(tsmaller.x, tsmaller.y), tsmaller.z);
    float tmax = min(min(tbigger.x, tbigger.y), tbigger.z);

    if (tmax < 0.0) {
        return false;
    }
    if (tmin > tmax) {
        return false;
    }

    float tHit = tmin;
    if (tHit < 0.001) {
        tHit = tmax;
        if (tHit < 0.001) {
            return false;
        }
    }

    tOut = tHit;
    return true;
}

// compute box normal using hit point
vec3 boxNormalAtPoint(vec3 p, Box box) {
    vec3 local = p - box.center;
    vec3 absLocal = abs(local);
    vec3 n = vec3(0.0);

    // compare normalized distances to half-size to find which face was hit
    float bx = absLocal.x / box.halfSize.x;
    float by = absLocal.y / box.halfSize.y;
    float bz = absLocal.z / box.halfSize.z;

    if (bx > by && bx > bz) {
        n.x = sign(local.x);
    } else if (by > bx && by > bz) {
        n.y = sign(local.y);
    } else {
        n.z = sign(local.z);
    }
    return n;
}

// Ice floor displacement texture
float iceHeight(vec2 p) {
    // Tile across floor
    vec2 uv = fract(p * 0.25);
    float hTex = texture(u_iceHeightMap, uv).r; // 0,1
    return hTex * 2.0 - 1.0; // remap to -1,1
}

vec3 iceBumpNormal(vec3 worldPos) {
    vec2 p = worldPos.xy * 0.6;
    float eps = 0.05;
    float hC = iceHeight(p);
    float hX1 = iceHeight(p + vec2(eps, 0.0));
    float hX2 = iceHeight(p - vec2(eps, 0.0));
    float hY1 = iceHeight(p + vec2(0.0, eps));
    float hY2 = iceHeight(p - vec2(0.0, eps));

    float dhdx = (hX1 - hX2) / (2.0 * eps);
    float dhdy = (hY1 - hY2) / (2.0 * eps);

    float bumpStrength = u_iceBumpStrength;
    vec3 n = normalize(vec3(-dhdx * bumpStrength, -dhdy * bumpStrength, 1.0));

    vec2 uv = fract(worldPos.xy * 0.25);
    vec3 nTex = texture(u_iceNormalMap, uv).xyz * 2.0 - 1.0;
    nTex = normalize(nTex);

    // Blend height normal and normal map normal
    n = normalize(mix(n, nTex, u_iceNormalBlend));
    return n;
}

vec3 lavaColor(vec3 worldPos) {
    vec2 uv = vec2(worldPos.x * 0.12 + u_time * 0.08,
                   worldPos.z * 0.18 + u_time * 0.03);
    uv = fract(uv);
    return texture(u_lavaTexture, uv).rgb;
}

HitData traceScene(Ray ray, bool includeLightSphere) {
    HitData closest;
    closest.hit = false;
    closest.t = 1e10;
    closest.isLava = false;

    // Light sphere intersection
    if (includeLightSphere) {
        for (int i = 0; i < NUM_LIGHT_SPHERES; i++) {
            float t;
            if (intersectSphere(ray, lightSpheres[i], t)) {
                if (t < closest.t) {
                    closest.hit = true;
                    closest.t = t;
                    closest.point = ray.origin + t * ray.direction;
                    vec3 n = normalize(closest.point - lightSpheres[i].center);
                    closest.frontFace = dot(ray.direction, n) < 0.0;
                    closest.normal = closest.frontFace ? n : -n;
                    closest.color = lightSpheres[i].color;
                    closest.material = lightSpheres[i].material;
                    closest.reflectivity = lightSpheres[i].reflectivity;
                    closest.refractiveIndex = lightSpheres[i].refractiveIndex;
                    closest.isLava = false;
                }
            }
        }
    }

    // Box intersections
    for (int i = 0; i < 4; i++) {
        float t;
        if (intersectBox(ray, boxes[i], t)) {
            if (t < closest.t) {
                closest.hit = true;
                closest.t = t;
                closest.point = ray.origin + t * ray.direction;
                vec3 n = boxNormalAtPoint(closest.point, boxes[i]);
                closest.frontFace = dot(ray.direction, n) < 0.0;
                closest.normal = closest.frontFace ? n : -n;
                closest.color = boxes[i].color;
                closest.material = boxes[i].material;
                closest.reflectivity = boxes[i].reflectivity;
                closest.refractiveIndex = 1.0;
                closest.isLava = boxes[i].isLava;
            }
        }
    }

    // Ball sphere intersection
    float tBall;
    if (intersectSphere(ray, spheres[0], tBall)) {
        if (tBall < closest.t) {
            closest.hit = true;
            closest.t = tBall;
            closest.point = ray.origin + tBall * ray.direction;
            vec3 n = normalize(closest.point - spheres[0].center);
            closest.frontFace = dot(ray.direction, n) < 0.0;
            closest.normal = closest.frontFace ? n : -n;
            closest.color = spheres[0].color;
            closest.material = spheres[0].material;
            closest.reflectivity = spheres[0].reflectivity;
            closest.refractiveIndex = spheres[0].refractiveIndex;
            closest.isLava = false;
        }
    }

    // Ice floor intersection
    float tFloor;
    if (intersectBox(ray, iceFloor, tFloor)) {
        if (tFloor < closest.t) {
            closest.hit = true;
            closest.t = tFloor;
            closest.point = ray.origin + tFloor * ray.direction;
            vec3 n = iceBumpNormal(closest.point);
            closest.frontFace = dot(ray.direction, n) < 0.0;
            closest.normal = closest.frontFace ? n : -n;
            closest.color = iceFloor.color;
            closest.material = iceFloor.material;
            closest.reflectivity = iceFloor.reflectivity;
            closest.refractiveIndex = 1.0;
            closest.isLava = false;
        }
    }

    return closest;
}

// refraction
vec3 refraction(vec3 I, vec3 N, float eta) {
    return refract(I, N, eta);
}

// Schlick
float schlick(float cosine, float ref_idx) {
    float r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
}

vec3 trace(Ray ray, int maxDepth) {
    vec3 color = vec3(0.0);
    vec3 attenuation = vec3(1.0);

    for (int depth = 0; depth < 8; depth++) {
        if (depth >= maxDepth) break;

        // Check if ray hits a light sphere
        HitData lightHit = traceScene(ray, true);
        if (lightHit.hit && lightHit.material == MATERIAL_LIGHT) {
            color += attenuation * lightHit.color;
            break;
        }

        // Trace only non-light geometry
        HitData hit = traceScene(ray, false);

        if (!hit.hit) {
            // No hit, gradient sky
            vec3 unitDir = normalize(ray.direction);
            float tEnv = 0.5 * (unitDir.z + 1.0);
            tEnv = clamp(tEnv, 0.0, 1.0);
            vec3 skyColor = mix(vec3(0.005, 0.005, 0.01), vec3(0.04, 0.07, 0.10), tEnv);

            color += attenuation * skyColor;
            break;
        }

        // Direct lighting from the light spheres
        vec3 surfaceColor = hit.color;
        if (hit.isLava) {
            surfaceColor = lavaColor(hit.point);
        }

        vec3 viewDir = normalize(-ray.direction);
        vec3 directLight = vec3(0.0);
        vec3 ambientLight = u_ambientStrength * surfaceColor;

        for (int i = 0; i < NUM_LIGHT_SPHERES; i++) {
            vec3 lightPos = lightSpheres[i].center;
            vec3 lightColor = lightSpheres[i].color;

            vec3 toLight = lightPos - hit.point;
            float lightDist = length(toLight);
            vec3 lightDir = toLight / max(lightDist, 0.0001);

            // Shadow ray - block light if an object is between the point and the light
            Ray shadowRay = Ray(hit.point + lightDir * 0.001, lightDir);
            HitData shadowHit = traceScene(shadowRay, false);
            bool inShadow = shadowHit.hit && shadowHit.t < lightDist;

            if (!inShadow) {
                float diff = max(dot(hit.normal, lightDir), 0.0);
                vec3 halfDir = normalize(lightDir + viewDir);
                float spec = pow(max(dot(hit.normal, halfDir), 0.0), 64.0);
                directLight += surfaceColor * diff * lightColor + spec * lightColor;
            }
        }

        if (hit.material == MATERIAL_DIFFUSE) {
            // Local illumination only (no more bounces)
            vec3 lighting = directLight + ambientLight;
            if (hit.isLava) {
                lighting += surfaceColor * 0.2;
            }
            color += attenuation * lighting;
            break;

        } else if (hit.material == MATERIAL_REFRACTIVE) {
            // Refract or reflect
            float etaRatio = hit.frontFace ? (1.0 / hit.refractiveIndex) : hit.refractiveIndex;
            float cosTheta = min(dot(-ray.direction, hit.normal), 1.0);
            float sin2Theta = 1.0 - cosTheta * cosTheta;
            bool total_internal_reflection = etaRatio * etaRatio * sin2Theta > 1.0;
            float schlick_value = schlick(cosTheta, hit.refractiveIndex);
            vec3 direction;
            if (total_internal_reflection || schlick_value > 0.5) {
                direction = reflect(ray.direction, hit.normal);
            } else {
                direction = refraction(ray.direction, hit.normal, etaRatio);
            }
            ray = Ray(hit.point + direction * 0.001, direction);
            attenuation *= hit.color;

        } else if (hit.material == MATERIAL_REFLECTIVE) {
            // Ice Floor
            color += attenuation * (0.02 * directLight + 0.02 * ambientLight);
            attenuation *= hit.reflectivity;
            ray = Ray(hit.point + hit.normal * 0.001, reflect(ray.direction, hit.normal));

        } else {
            // Fallback for other refractive materials
            color += attenuation * (0.02 * directLight + 0.02 * ambientLight);
            attenuation *= hit.reflectivity;
            ray = Ray(hit.point + hit.normal * 0.001, reflect(ray.direction, hit.normal));
        }

        if (length(attenuation) < 0.01) break;
    }

    return color;
}

void main() {
    initScene(u_time);

    vec2 uv = (gl_FragCoord.xy - 0.5 * u_resolution) / u_resolution.y;

    // Build ray from camera
    float focalDist = 1.5;
    vec3 direction = normalize(
        u_camForward * focalDist +
        u_camRight * uv.x +
        u_camUp * uv.y
    );

    Ray ray = Ray(u_camPos, direction);
    vec3 color = trace(ray, u_maxBounces);

    color = pow(color, vec3(1.0 / 2.2));

    fragColor = vec4(color, 1.0);
}