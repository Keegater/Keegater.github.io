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
uniform float    u_iceBumpStrength;
uniform float    u_iceNormalBlend;

out vec4 fragColor;

// Material types
const int MATERIAL_DIFFUSE = 0;
const int MATERIAL_REFLECTIVE = 1;
const int MATERIAL_REFRACTIVE = 2;
const int MATERIAL_LIGHT = 3;

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
};

// Scene objects
Sphere spheres[1];          // ball sphere
Sphere lightSphere;         // not using rn
Box boxes[4];               // left paddle, right paddle, top/bottom bounds
Box iceFloor;               // reflective ice floor

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
        0.0
    );

    // Right paddle (player)
    boxes[1] = Box(
        vec3(4.5, u_rightPaddleY, 0.0),
        vec3(0.05, 0.5, 0.2),
        vec3(0.8, 0.8, 0.8),
        MATERIAL_DIFFUSE,
        0.0
    );

    // Bound boards top and bottom
    float boardTop = 3.25;
    float boardBottom = -3.25;
    float boardHalfWidth = 4.5;
    float borderHalfThickness = 0.25;

    // Top border
    boxes[2] = Box(
        vec3(0.0, boardTop + borderHalfThickness, 0.0),
        vec3(boardHalfWidth + 1.0, borderHalfThickness, 0.25),
        vec3(0.0, 0.0, 0.0),
        MATERIAL_DIFFUSE,
        0.0
    );

    // Bottom border
    boxes[3] = Box(
        vec3(0.0, boardBottom - borderHalfThickness, 0.0),
        vec3(boardHalfWidth + 1.0, borderHalfThickness, 0.25),
        vec3(0.0, 0.0, 0.0),
        MATERIAL_DIFFUSE,
        0.0
    );

    // Ice floor
    iceFloor = Box(
        vec3(0.0, 0.0, -0.8),
        vec3(boardHalfWidth + 1.0, (boardTop - boardBottom) * 0.55, 0.6),
        vec3(0.1, 0.45, 0.8),
        MATERIAL_REFLECTIVE,
        0.6
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

HitData traceScene(Ray ray, bool includeLightSphere) {
    HitData closest;
    closest.hit = false;
    closest.t = 1e10;
/*
    // Light sphere intersection
    if (includeLightSphere) {
        float t;
        if (intersectSphere(ray, lightSphere, t)) {
            if (t < closest.t) {
                closest.hit = true;
                closest.t = t;
                closest.point = ray.origin + t * ray.direction;
                vec3 n = normalize(closest.point - lightSphere.center);
                closest.frontFace = dot(ray.direction, n) < 0.0;
                closest.normal = closest.frontFace ? n : -n;
                closest.color = lightSphere.color;
                closest.material = lightSphere.material;
                closest.reflectivity = lightSphere.reflectivity;
                closest.refractiveIndex = lightSphere.refractiveIndex;
            }
        }
    }
    */

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

/*
        HitData lightHit = traceScene(ray, true);
        if (lightHit.hit && lightHit.material == MATERIAL_LIGHT) {
            color += attenuation * lightHit.color;
            break;
        }
*/

        HitData hit = traceScene(ray, false);

        if (!hit.hit) {
            // No hit
            vec3 skyColor = vec3(0.047, 0.047, 0.047);
            color += attenuation * skyColor;
            break;
        }
/*
        // Lighting from central light sphere
        vec3 viewDir = normalize(-ray.direction);
        vec3 lightDir1 = normalize(lightSphere.center - hit.point);
        float lightDist1 = length(lightSphere.center - hit.point);

        Ray shadowRay1 = Ray(hit.point + lightDir1 * 0.001, lightDir1);
        HitData shadowHit1 = traceScene(shadowRay1, false);
        bool inShadow1 = shadowHit1.hit && shadowHit1.t < lightDist1 - 0.001;

        float diff1 = max(dot(hit.normal, lightDir1), 0.0);
        vec3 halfDir1 = normalize(lightDir1 + viewDir);
        float spec1 = pow(max(dot(hit.normal, halfDir1), 0.0), 64.0);
        vec3 lightColor1 = lightSphere.color;

        vec3 contrib1 = vec3(0.0);
        if (!inShadow1) {
            contrib1 = hit.color * diff1 * lightColor1 + spec1 * lightColor1;
        }

        vec3 directLight = contrib1;
*/
        vec3 lightDir = normalize(vec3(1.0, 1.0, 0.8));
        vec3 lightColor = vec3(1.0, 1.0, 1.0);

        float diff = max(dot(hit.normal, lightDir), 0.0);
        vec3 directLight = diff * lightColor;

        vec3 ambientLight = u_ambientStrength * hit.color;

        if (hit.material == MATERIAL_DIFFUSE) {
            vec3 lighting = directLight + ambientLight;
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
            // Ice floor
            vec3 glowColor = u_centerLightColor * u_centerLightBrightness * 0.15;
            color += attenuation * (0.02 * directLight + 0.02 * ambientLight + glowColor);
            attenuation *= hit.reflectivity;
            ray = Ray(hit.point + hit.normal * 0.001, reflect(ray.direction, hit.normal));

        } else {
            // Fallback for other reflective materials
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