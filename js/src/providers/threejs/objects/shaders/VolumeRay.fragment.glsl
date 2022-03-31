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

void getChildren(in int inode, out int[8] children) {

}

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

vec3 vec3tovec2(vec3 v) {
	return v;
}

vec3 worldGetNormal(in float px, in vec3 pos)
{
	return normalize(
		vec3(px -  texture(volumeTexture, vec3tovec2(pos + vec3(gradient_step, 0, 0))).x,
			 px -  texture(volumeTexture, vec3tovec2(pos + vec3(0, gradient_step, 0))).x,
			 px -  texture(volumeTexture, vec3tovec2(pos + vec3(0, 0, gradient_step))).x
		 )
	);
}

float getShadow(vec3 textcoord, vec2 sliceCount)
{
	float zidx1 = floor(textcoord.z * lightMapSize.z);
	float zidx2 = ceil(textcoord.z * lightMapSize.z);

	float shadow1 = texture2D(shadowTexture,
		vec2(
			floor(mod(zidx1, sliceCount.x)) * lightMapSize.x / lightMapRenderTargetSize.x,
			floor(zidx1 / sliceCount.x)* lightMapSize.y / lightMapRenderTargetSize.y
		)
		+ vec2(textcoord.x / sliceCount.x, textcoord.y / sliceCount.y)
	).r;

	float shadow2 = texture2D(shadowTexture,
		vec2(
			floor(mod(zidx2, sliceCount.x)) * lightMapSize.x / lightMapRenderTargetSize.x,
			floor(zidx2 / sliceCount.x)* lightMapSize.y / lightMapRenderTargetSize.y
		)
		+ vec2(textcoord.x / sliceCount.x, textcoord.y / sliceCount.y)
	).r;

	return mix(shadow1, shadow2, textcoord.z * lightMapSize.z - zidx1);
}

// void proc_subtree (
//     in float tx0, in float tx1,
//     in float ty0, in float ty1,
//     in float tz0, in float tz1,
//     in uint a, in int inode,
// ) {
//     if (tx1 < 0 || ty1 < 0 || tz1 < 0) return;

// }

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

	vec4 value = vec4(0.0, 0.0, 0.0, 1.0);
	vec3 direction = normalize(transformedWorldPosition - transformedCameraPosition);

	// Find entry point
	Ray ray = makeRay(transformedCameraPosition, direction);
	intersect(ray, aabb, txmin, txmax, tymin, tymax, tzmin, tzmax, tmin, tmax);

	value = texture(volumeTexture, vec3(1, 0, 0));
	value.a = 1.0;

	gl_FragColor = value;
}