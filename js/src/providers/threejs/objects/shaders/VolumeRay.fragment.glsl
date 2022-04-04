#include <common>
#include <clipping_planes_pars_fragment>
#include <lights_pars_begin>

precision highp sampler3D;

uniform vec3 lightMapSize;
uniform vec2 lightMapRenderTargetSize;
uniform sampler2D shadowTexture;

uniform mat4 transform;
uniform sampler3D childTextureLowerBits;
uniform sampler3D childTextureUpperBits;
uniform sampler3D volumeTexture;
uniform int Npack;

uniform sampler2D colormap;
uniform sampler2D jitterTexture;
uniform float focal_length;
uniform float focal_plane;
uniform float low;
uniform float high;
uniform mat4 modelViewMatrix;
uniform float samples;
uniform float alpha_coef;
uniform float gradient_step;

uniform vec4 scale;
uniform vec4 translation;

varying vec3 localPosition;
varying vec3 transformedCameraPosition;
varying vec3 transformedWorldPosition;

float inv_range;

struct Ray {
    vec3 origin;
    vec3 direction;
    vec3 inv_direction;
    uint a;
    int sign[3];
};

vec3 aabb[2] = vec3[2](
    vec3(-0.5, -0.5, -0.5),
    vec3(0.5, 0.5, 0.5)
);


Ray makeRay(vec3 origin, vec3 direction) {
    vec3 inv_direction = vec3(1.0) / direction;
    uint a = (
        ((inv_direction.x < 0.0) ? uint(4) : uint(0)) +
        ((inv_direction.y < 0.0) ? uint(2) : uint(0)) +
        ((inv_direction.z < 0.0) ? uint(1) : uint(0))
    );
    return Ray(
        origin,
        direction,
        inv_direction,
        a,
        int[3](
            ((inv_direction.x < 0.0) ? 1 : 0),
            ((inv_direction.y < 0.0) ? 1 : 0),
            ((inv_direction.z < 0.0) ? 1 : 0)
        )
    );
}

/*
    From: https://github.com/hpicgs/cgsee/wiki/Ray-Box-Intersection-on-the-GPU
*/
void intersect(
    in Ray ray, in vec3 aabb[2],
    out float txmin, out float txmax,
    out float tymin, out float tymax,
    out float tzmin, out float tzmax,
    out float tmin, out float tmax
){
    txmin = (aabb[ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
    txmax = (aabb[1-ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
    tymin = (aabb[ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
    tymax = (aabb[1-ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
    tzmin = (aabb[ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
    tzmax = (aabb[1-ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
    tmin = max(max(txmin, tymin), tzmin);
    tmax = min(min(txmax, tymax), tzmax);
}


vec3 index2coords(const in int index, const in int N) {
    int i, j, k;
    int N2 = N * N;
    i = ((index) / N2);
    j = ((index) - i * N2) / N;
    k = ((index) - i * N2 - j * N);

    return vec3(i, j, k);
}

int getChildren(const in int inode, out int children[8]) {
    int Nchild = 0;
    vec3 coords;
    for (int ind_child = 0; ind_child < 8; ind_child++) {
        if ((inode + ind_child < 0) && (inode + ind_child >= 125)) {
            children[ind_child] = -1;
            continue;
        }

        // TODO: do not hardcode shape of 5x5x5
        coords = index2coords(inode + ind_child, 5);

        children[ind_child] = (
            int(texture(childTextureLowerBits, coords)) +
            int(texture(childTextureUpperBits, coords)) * 65536
        );

        if (children[ind_child] >= 0) {
            Nchild += 1;
        }
    }

    return Nchild;
}

void integrateData(const in int inode, inout vec4 value) {
    float px;
    vec4 pxColor = vec4(0.0);

    // FIXME: inode should be 0-indexed
    px = texture(volumeTexture, index2coords(inode - 1, 2)).x;
    float scaled_px = px ; // TOOD: (px - low) * inv_range;
    float step = 1.0;

    if (scaled_px > 0.0) {
        scaled_px = min(scaled_px, 0.99);

        pxColor = texture(colormap, vec2(scaled_px, 0.5));

        // FIXME: compute step size from node depth and tin/tout
        pxColor.a = 1.0 - pow(1.0 - pxColor.a, step * alpha_coef);
        pxColor.a *= (1.0 - value.a);

        pxColor.rgb *= pxColor.a;

        value += pxColor;
    }
}

uint first_node(
    const in float tx0, const in float ty0, const in float tz0,
    const in float txm, const in float tym, const in float tzm
) {
    uint index = uint(0);

    if (tx0 >= max(ty0, tz0)) {              // enters YZ plane
        if (tym < tx0) index |= uint(2);
        if (tzm < tx0) index |= uint(1);
    } else if (ty0 >= max(tx0, tz0))  {      // enters XZ plane
        if (txm < ty0) index |= uint(4);
        if (tzm < ty0) index |= uint(1);
    } else {                                 // enters XY plane
        if (txm < tz0) index |= uint(4);
        if (tym < tz0) index |= uint(2);
    }
    return index;
}

uint next_node(
    const in float tx, const in float ty, const in float tz,
    const in uint ix, const in uint iy, const in uint iz
) {
    if (tx < min(ty, tz)) {         // YZ plane
        return ix;
    } else if (ty < min(tx, tz)) {  // XZ plane
        return iy;
    } else {                        // XY plane
        return iz;
    }
}

void proc_subtree (
    const in float tx0, const in float ty0, const in float tz0,
    const in float tx1, const in float ty1, const in float tz1,
    const in int inode, const in uint a,
    inout vec4 value
) {
    if (tx1 < 0. || ty1 < 0. || tz1 < 0.) return;

    // Compute children
    int children[8];
    int Nchild;
    Nchild = getChildren(inode, children);

    // Process leaf node
    if (Nchild == 0) {
        integrateData(inode, value);
    }

    float txm, tym, tzm;
    txm = (tx0 + tx1) * 0.5;
    tym = (ty0 + ty1) * 0.5;
    tzm = (tz0 + tz1) * 0.5;

    uint ind_child = first_node(tx0, ty0, tz0, txm, tym, tzm);

    do {
        switch (ind_child) {
        case uint(0):
            proc_subtree(tx0, ty0, tz0, txm, tym, tzm, children[a], a, value);
            ind_child = next_node(txm, tym, tzm, uint(4), uint(2), uint(1));
            break;
        case uint(1):
            proc_subtree(tx0, ty0, tzm, txm, tym, tz1, children[uint(1)^a], a, value);
            ind_child = next_node(txm, tym, tz1, uint(5), uint(3), uint(8));
            break;
        case uint(2):
            proc_subtree(tx0, tym, tz0, txm, ty1, tzm, children[uint(2)^a], a, value);
            ind_child = next_node(txm, ty1, tzm, uint(6), uint(8), uint(3));
            break;
        case uint(3):
            proc_subtree(tx0, tym, tzm, txm, ty1, tz1, children[uint(3)^a], a, value);
            ind_child = next_node(txm, ty1, tz1, uint(7), uint(8), uint(8));
            break;
        case uint(4):
            proc_subtree(txm, ty0, tz0, tx1, tym, tzm, children[uint(4)^a], a, value);
            ind_child = next_node(tx1, tym, tzm, uint(8), uint(6), uint(5));
            break;
        case uint(5):
            proc_subtree(txm, ty0, tzm, tx1, tym, tz1, children[uint(5)^a], a, value);
            ind_child = next_node(tx1, tym, tz1, uint(8), uint(7), uint(8));
            break;
        case uint(6):
            proc_subtree(txm, tym, tz0, tx1, ty1, tzm, children[uint(6)^a], a, value);
            ind_child = next_node(tx1, ty1, tzm, uint(8), uint(8), uint(7));
            break;
        case uint(7):
            proc_subtree(txm, tym, tzm, tx1, ty1, tz1, children[uint(7)^a], a, value);
            ind_child = uint(8);
            break;
        };
    } while (ind_child < uint(8));
}

void main() {
    float tmin = 0.0;
    float txmin = 0.0;
    float tymin = 0.0;
    float tzmin = 0.0;
    float tmax = 0.0;
    float txmax = 0.0;
    float tymax = 0.0;
    float tzmax = 0.0;
    float px = 0.0;
    float shadow = 0.0;
    vec4 pxColor = vec4(0.0, 0.0, 0.0, 0.0);

    inv_range = 1.0 / (high - low);
    aabb[0] = aabb[0] * scale.xyz + translation.xyz;
    aabb[1] = aabb[1] * scale.xyz + translation.xyz;

    vec4 accuColor = vec4(0.0, 0.0, 0.0, 0.0);

    vec4 value = vec4(0.0, 0.0, 0.0, 0.0);
    vec3 direction = normalize(transformedWorldPosition - transformedCameraPosition);

    // Find entry point
    Ray ray = makeRay(transformedCameraPosition, direction);
    intersect(ray, aabb, txmin, txmax, tymin, tymax, tzmin, tzmax, tmin, tmax);

    // If entry point is smaller than exit point, find path in octree
    if (max(max(txmin, tymin), tzmin) < min(min(txmax, tymax), tzmax)) {
        proc_subtree(
            txmin, tymin, tzmin,
            txmax, tymax, tzmax,
            0,
            ray.a,
            value
        );
    }

    gl_FragColor = value;
}