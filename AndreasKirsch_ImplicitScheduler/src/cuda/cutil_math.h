////////////////////////////////////////////////////////////////////////////////
// Essential functions on CUDA built in type
// Linh Ha lha@cs.utah.edu
// The vector / vector and vector * vector is point-wise
// Use it at your own risk
////////////////////////////////////////////////////////////////////////////////

#ifndef CUTIL_MATH_H
#define CUTIL_MATH_H

#include "cuda_runtime.h"

typedef unsigned int uint;
typedef unsigned short ushort;

#ifndef __CUDACC__
#include <stdlib.h>
#include <math.h>

////////////////////////////////////////////////////////////////////////////////
// host implementations of CUDA functions
////////////////////////////////////////////////////////////////////////////////

inline float fminf(float a, float b){
    return a < b ? a : b;
}

inline float fmaxf(float a, float b){
    return a > b ? a : b;
}

inline int max(int a, int b){
    return a > b ? a : b;
}

inline int min(int a, int b){
    return a < b ? a : b;
}

inline uint max(uint a, uint b){
    return a > b ? a : b;
}

inline uint min(uint a, uint b){
    return a < b ? a : b;
}

inline float rsqrtf(float x){
    return 1.0f / sqrtf(x);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Base functions
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ float clamp(float f, float a, float b){
    return fmaxf(a, fminf(f, b));
}
inline __device__ __host__ int clamp(int f, int a, int b){
    return max(a, min(f, b));
}
inline __device__ __host__ uint clamp(uint f, uint a, uint b){
    return max(a, min(f, b));
}

// smoothstep
inline __device__ __host__ float smoothstep(float a, float b, float x)
{
	float y = clamp((x - a) / (b - a), 0.0f, 1.0f);
	return (y*y*(3.0f - (2.0f*y)));
}
////////////////////////////////////////////////////////////////////////////////
// additional constructors
////////////////////////////////////////////////////////////////////////////////
__inline__ __host__ __device__ float2 make_float2(float s){
    return make_float2(s, s);
}
inline __host__ __device__ float2 make_float2(int2 a){
    return make_float2(float(a.x), float(a.y));
}
inline __host__ __device__ float2 make_float2(uint2 a){
    return make_float2(float(a.x), float(a.y));
}

inline __host__ __device__ int2 make_int2(int s){
    return make_int2(s, s);
}
inline __host__ __device__ int2 make_int2(uint2 a){
    return make_int2(int(a.x), int(a.y));
}

inline __host__ __device__ uint2 make_uint2(uint s){
    return make_uint2(s, s);
}

// discards z
inline __host__ __device__ float2 make_float2(float3 a){
    return make_float2(a.x, a.y);
}
inline __host__ __device__ int2 make_int2(int3 a){
    return make_int2(a.x, a.y);
}
inline __host__ __device__ uint2 make_uint2(uint3 a){
    return make_uint2(a.x, a.y);
}

// additional constructors
inline __host__ __device__ float3 make_float3(float s){
    return make_float3(s, s, s);
}
inline __host__ __device__ float3 make_float3(float2 a){
    return make_float3(a.x, a.y, 0.0f);
}
inline __host__ __device__ float3 make_float3(float2 a, float s){
    return make_float3(a.x, a.y, s);
}

inline __host__ __device__ int3 make_int3(int s){
    return make_int3(s, s, s);
}
inline __host__ __device__ int3 make_int3(int2 a){
    return make_int3(a.x, a.y, 0);
}
inline __host__ __device__ int3 make_int3(int2 a, int s){
    return make_int3(a.x, a.y, s);
}

inline __host__ __device__ uint3 make_uint3(uint s){
    return make_uint3(s, s, s);
}
inline __host__ __device__ uint3 make_uint3(uint2 a){
    return make_uint3(a.x, a.y, 0);
}
inline __host__ __device__ uint3 make_uint3(uint2 a, uint s){
    return make_uint3(a.x, a.y, s);
}

// change type
inline __host__ __device__ float3 make_float3(int3 a){
    return make_float3(float(a.x), float(a.y), float(a.z));
}
inline __host__ __device__ float3 make_float3(uint3 a){
    return make_float3(float(a.x), float(a.y), float(a.z));
}
inline __host__ __device__ int3 make_int3(uint3 a){
    return make_int3(int(a.x), int(a.y), int(a.z));
}

// discards w
inline __host__ __device__ float3 make_float3(float4 a){
    return make_float3(a.x, a.y, a.z);
}
inline __host__ __device__ int3 make_int3(int4 a){
    return make_int3(a.x, a.y, a.z);
}
inline __host__ __device__ uint3 make_uint3(uint4 a){
    return make_uint3(a.x, a.y, a.z);
}

// additional constructors float4
inline __host__ __device__ float4 make_float4(float s){
    return make_float4(s, s, s, s);
}
inline __host__ __device__ float4 make_float4(float3 a){
    return make_float4(a.x, a.y, a.z, 0.0f);
}
inline __host__ __device__ float4 make_float4(float3 a, float w){
    return make_float4(a.x, a.y, a.z, w);
}

inline __host__ __device__ int4 make_int4(int s){
    return make_int4(s, s, s, s);
}
inline __host__ __device__ int4 make_int4(int3 a){
    return make_int4(a.x, a.y, a.z, 0);
}
inline __host__ __device__ int4 make_int4(int3 a, int w){
    return make_int4(a.x, a.y, a.z, w);
}

inline __host__ __device__ uint4 make_uint4(uint s){
    return make_uint4(s, s, s, s);
}
inline __host__ __device__ uint4 make_uint4(uint3 a){
    return make_uint4(a.x, a.y, a.z, 0);
}
inline __host__ __device__ uint4 make_uint4(uint3 a, uint w){
    return make_uint4(a.x, a.y, a.z, w);
}

inline __host__ __device__ float4 make_float4(int4 a){
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}
inline __host__ __device__ float4 make_float4(uint4 a){
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}
inline __host__ __device__ int4 make_int4(uint4 a){
    return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}

////////////////////////////////////////////////////////////////////////////////
// negate
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float2 operator-(float2 a){
    return make_float2(-a.x, -a.y);
}

inline __host__ __device__ int2 operator-(int2 a){
    return make_int2(-a.x, -a.y);
}

inline __host__ __device__ float3 operator-(float3 a){
    return make_float3(-a.x, -a.y, -a.z);
}
inline __host__ __device__ int3 operator-(int3 a){
    return make_int3(-a.x, -a.y, -a.z);
}

inline __host__ __device__ float4 operator-(float4 a){
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}

inline __host__ __device__ int4 operator-(int4 a){
    return make_int4(-a.x, -a.y, -a.z, -a.w);
}
////////////////////////////////////////////////////////////////////////////////
// addition
////////////////////////////////////////////////////////////////////////////////
// addition 2
inline __host__ __device__ float2 operator+(float2 a, float2 b){
    return make_float2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(float2 &a, float2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ float2 operator+(float2 a, float b){
    return make_float2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(float2 &a, float b){
    a.x += b; a.y += b;
}

inline __host__ __device__ int2 operator+(int2 a, int2 b){
    return make_int2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(int2 &a, int2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ int2 operator+(int2 a, int b){
    return make_int2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(int2 &a, int b){
    a.x += b; a.y += b;
}

inline __host__ __device__ uint2 operator+(uint2 a, uint2 b){
    return make_uint2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(uint2 &a, uint2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ uint2 operator+(uint2 a, uint b){
    return make_uint2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(uint2 &a, uint b){
    a.x += b; a.y += b;
}

inline __host__ __device__ float2 operator+(float b, float2 a){
    return make_float2(a.x + b, a.y + b);
}
inline __host__ __device__ uint2 operator+(uint b, uint2 a){
    return make_uint2(a.x + b, a.y + b);
}
inline __host__ __device__ int2 operator+(int b, int2 a){
    return make_int2(a.x + b, a.y + b);
}

// addition 3
inline __host__ __device__ float3 operator+(float3 a, float3 b){
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(float3 &a, float3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ float3 operator+(float3 a, float b){
    return make_float3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(float3 &a, float b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ int3 operator+(int3 a, int3 b){
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(int3 &a, int3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ int3 operator+(int3 a, int b){
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(int3 &a, int b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ uint3 operator+(uint3 a, uint3 b){
    return make_uint3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(uint3 &a, uint3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ uint3 operator+(uint3 a, uint b){
    return make_uint3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(uint3 &a, uint b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ int3 operator+(int b, int3 a){
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ uint3 operator+(uint b, uint3 a){
    return make_uint3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ float3 operator+(float b, float3 a){
    return make_float3(a.x + b, a.y + b, a.z + b);
}

// addition 4
inline __host__ __device__ float4 operator+(float4 a, float4 b){
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(float4 &a, float4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ float4 operator+(float4 a, float b){
    return make_float4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(float4 &a, float b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline __host__ __device__ int4 operator+(int4 a, int4 b){
    return make_int4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(int4 &a, int4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ int4 operator+(int4 a, int b){
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(int4 &a, int b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline __host__ __device__ uint4 operator+(uint4 a, uint4 b){
    return make_uint4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(uint4 &a, uint4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ uint4 operator+(uint4 a, uint b){
    return make_uint4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(uint4 &a, uint b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline __host__ __device__ float4 operator+(float b, float4 a){
    return make_float4(a.x + b, a.y + b, a.z + b,  a.w + b);
}

inline __host__ __device__ int4 operator+(int b, int4 a){
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}

inline __host__ __device__ uint4 operator+(uint b, uint4 a){
    return make_uint4(a.x + b, a.y + b, a.z + b,  a.w + b);
}

////////////////////////////////////////////////////////////////////////////////
// subtract
////////////////////////////////////////////////////////////////////////////////
// subtract 2
inline __host__ __device__ float2 operator-(float2 a, float2 b){
    return make_float2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(float2 &a, float2 b){
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ float2 operator-(float2 a, float b){
    return make_float2(a.x - b, a.y - b);
}
inline __host__ __device__ void operator-=(float2 &a, float b){
    a.x -= b; a.y -= b;
}

inline __host__ __device__ int2 operator-(int2 a, int2 b){
    return make_int2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(int2 &a, int2 b){
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ int2 operator-(int2 a, int b){
    return make_int2(a.x - b, a.y - b);
}
inline __host__ __device__ void operator-=(int2 &a, int b){
    a.x -= b; a.y -= b;
}

inline __host__ __device__ uint2 operator-(uint2 a, uint2 b){
    return make_uint2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(uint2 &a, uint2 b){
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ uint2 operator-(uint2 a, uint b){
    return make_uint2(a.x - b, a.y - b);
}
inline __host__ __device__ void operator-=(uint2 &a, uint b){
    a.x -= b; a.y -= b;
}

inline __host__ __device__ float2 operator-(float b, float2 a){
    return make_float2(b - a.x, b - a.y);
}

inline __host__ __device__ int2 operator-(int b, int2 a){
    return make_int2(b - a.x, b - a.y);
}

inline __host__ __device__ uint2 operator-(uint b, uint2 a){
    return make_uint2(b - a.x, b - a.y);
}

// subtract 3
inline __host__ __device__ float3 operator-(float3 a, float3 b){
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(float3 &a, float3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ float3 operator-(float3 a, float b){
    return make_float3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(float3 &a, float b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ int3 operator-(int3 a, int3 b){
    return make_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(int3 &a, int3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ int3 operator-(int3 a, int b){
    return make_int3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(int3 &a, int b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ uint3 operator-(uint3 a, uint3 b){
    return make_uint3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(uint3 &a, uint3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ uint3 operator-(uint3 a, uint b){
    return make_uint3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(uint3 &a, uint b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ float3 operator-(float b, float3 a){
    return make_float3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ int3 operator-(int b, int3 a){
    return make_int3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ uint3 operator-(uint b, uint3 a){
    return make_uint3(b - a.x, b - a.y, b - a.z);
}

// subtract 4
inline __host__ __device__ float4 operator-(float4 a, float4 b){
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(float4 &a, float4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ float4 operator-(float4 a, float b){
    return make_float4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ void operator-=(float4 &a, float b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline __host__ __device__ int4 operator-(int4 a, int4 b){
    return make_int4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(int4 &a, int4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ int4 operator-(int4 a, int b){
    return make_int4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ void operator-=(int4 &a, int b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline __host__ __device__ uint4 operator-(uint4 a, uint4 b){
    return make_uint4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(uint4 &a, uint4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ uint4 operator-(uint4 a, uint b){
    return make_uint4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ void operator-=(uint4 &a, uint b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline __host__ __device__ float4 operator-(float b, float4 a){
    return make_float4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline __host__ __device__ uint4 operator-(uint b, uint4 a){
    return make_uint4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline __host__ __device__ int4 operator-(int b, int4 a){
    return make_int4(b - a.x, b - a.y, b - a.z, b - a.w);
}

////////////////////////////////////////////////////////////////////////////////
// Multiply
////////////////////////////////////////////////////////////////////////////////
// Multiply 2
inline __host__ __device__ float2 operator*(float2 a, float2 b){
    return make_float2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(float2 &a, float2 b){
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ float2 operator*(float2 a, float b){
    return make_float2(a.x * b, a.y * b);
}
inline __host__ __device__ void operator*=(float2 &a, float b){
    a.x *= b; a.y *= b;
}

inline __host__ __device__ int2 operator*(int2 a, int2 b){
    return make_int2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(int2 &a, int2 b){
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ int2 operator*(int2 a, int b){
    return make_int2(a.x * b, a.y * b);
}
inline __host__ __device__ void operator*=(int2 &a, int b){
    a.x *= b; a.y *= b;
}

inline __host__ __device__ uint2 operator*(uint2 a, uint2 b){
    return make_uint2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(uint2 &a, uint2 b){
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ uint2 operator*(uint2 a, uint b){
    return make_uint2(a.x * b, a.y * b);
}
inline __host__ __device__ void operator*=(uint2 &a, uint b){
    a.x *= b; a.y *= b;
}

inline __host__ __device__ float2 operator*(float b, float2 a){
    return make_float2(b * a.x, b * a.y);
}
inline __host__ __device__ int2 operator*(int b, int2 a){
    return make_int2(b * a.x, b * a.y);
}
inline __host__ __device__ uint2 operator*(uint b, uint2 a){
    return make_uint2(b * a.x, b * a.y);
}

// Multiply 3
inline __host__ __device__ float3 operator*(float3 a, float3 b){
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(float3 &a, float3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ float3 operator*(float3 a, float b){
    return make_float3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ void operator*=(float3 &a, float b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ int3 operator*(int3 a, int3 b){
    return make_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(int3 &a, int3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ int3 operator*(int3 a, int b){
    return make_int3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ void operator*=(int3 &a, int b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ uint3 operator*(uint3 a, uint3 b){
    return make_uint3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(uint3 &a, uint3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ uint3 operator*(uint3 a, uint b){
    return make_uint3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ void operator*=(uint3 &a, uint b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ float3 operator*(float b, float3 a){
    return make_float3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ int3 operator*(int b, int3 a){
    return make_int3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ uint3 operator*(uint b, uint3 a){
    return make_uint3(b * a.x, b * a.y, b * a.z);
}

//multiply 4
inline __host__ __device__ float4 operator*(float4 a, float4 b){
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(float4 &a, float4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ float4 operator*(float4 a, float b){
    return make_float4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ void operator*=(float4 &a, float b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __host__ __device__ int4 operator*(int4 a, int4 b){
    return make_int4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(int4 &a, int4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ int4 operator*(int4 a, int b){
    return make_int4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ void operator*=(int4 &a, int b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __host__ __device__ uint4 operator*(uint4 a, uint4 b){
    return make_uint4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(uint4 &a, uint4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ uint4 operator*(uint4 a, uint b){
    return make_uint4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ void operator*=(uint4 &a, uint b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __host__ __device__ float4 operator*(float b, float4 a){
    return make_float4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ uint4 operator*(uint b, uint4 a){
    return make_uint4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ int4 operator*(int b, int4 a){
    return make_int4(b * a.x, b * a.y, b * a.z, b * a.w);
}

////////////////////////////////////////////////////////////////////////////////
// Divide
////////////////////////////////////////////////////////////////////////////////
// divide 2
inline __host__ __device__ float2 operator/(float2 a, float2 b){
    return make_float2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ void operator/=(float2 &a, float2 b){
    a.x /= b.x; a.y /= b.y;
}
inline __host__ __device__ float2 operator/(float2 a, float b){
    float inv = 1.f / b;
    return make_float2(a.x * inv, a.y *inv);
}
inline __host__ __device__ void operator/=(float2 &a, float b){
    float inv = 1.f / b;
    a.x *= inv; a.y *= inv;
}
inline __host__ __device__ float2 operator/(float b, float2 a){
    return make_float2(b /a.x, b /a.y);
}

// divide 3
inline __host__ __device__ float3 operator/(float3 a, float3 b){
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ void operator/=(float3 &a, float3 b){
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}
inline __host__ __device__ float3 operator/(float3 a, float b){
    float inv = 1.f / b;
    return make_float3(a.x * inv, a.y * inv, a.z * inv);
}
inline __host__ __device__ void operator/=(float3 &a, float b){
    float inv = 1.f / b;
    a.x *= inv; a.y *= inv; a.z *= inv;
}
inline __host__ __device__ float3 operator/(float b, float3 a){
    return make_float3(b / a.x, b / a.y, b / a.z);
}

// divide 4
inline __host__ __device__ float4 operator/(float4 a, float4 b){
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z,  a.w / b.w);
}
inline __host__ __device__ void operator/=(float4 &a, float4 b){
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}
inline __host__ __device__ float4 operator/(float4 a, float b){
    float inv = 1.f / b;
    return make_float4(a.x * inv, a.y * inv, a.z * inv,  a.w * inv);
}
inline __host__ __device__ void operator/=(float4 &a, float b){
    float inv = 1.f / b;
    a.x *= inv; a.y *= inv; a.z *= inv; a.w *= inv;
}
inline __host__ __device__ float4 operator/(float b, float4 a){
    return make_float4(b / a.x, b / a.y, b / a.z, b / a.w);
}

////////////////////////////////////////////////////////////////////////////////
// min
////////////////////////////////////////////////////////////////////////////////
inline  __host__ __device__ float2 fminf(float2 a, float2 b){
	return make_float2(fminf(a.x,b.x), fminf(a.y,b.y));
}
inline __host__ __device__ float3 fminf(float3 a, float3 b){
	return make_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}
inline  __host__ __device__ float4 fminf(float4 a, float4 b){
	return make_float4(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z), fminf(a.w,b.w));
}

inline __host__ __device__ int2 min(int2 a, int2 b){
    return make_int2(min(a.x,b.x), min(a.y,b.y));
}
inline __host__ __device__ int3 min(int3 a, int3 b){
    return make_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline __host__ __device__ int4 min(int4 a, int4 b){
    return make_int4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}

inline __host__ __device__ uint2 min(uint2 a, uint2 b){
    return make_uint2(min(a.x,b.x), min(a.y,b.y));
}
inline __host__ __device__ uint3 min(uint3 a, uint3 b){
    return make_uint3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline __host__ __device__ uint4 min(uint4 a, uint4 b){
    return make_uint4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}

////////////////////////////////////////////////////////////////////////////////
// max
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float2 fmaxf(float2 a, float2 b){
	return make_float2(fmaxf(a.x,b.x), fmaxf(a.y,b.y));
}
inline __host__ __device__ float3 fmaxf(float3 a, float3 b){
	return make_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}
inline __host__ __device__ float4 fmaxf(float4 a, float4 b){
	return make_float4(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z), fmaxf(a.w,b.w));
}

inline __host__ __device__ int2 max(int2 a, int2 b){
    return make_int2(max(a.x,b.x), max(a.y,b.y));
}
inline __host__ __device__ int3 max(int3 a, int3 b){
    return make_int3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}
inline __host__ __device__ int4 max(int4 a, int4 b){
    return make_int4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}

inline __host__ __device__ uint2 max(uint2 a, uint2 b){
    return make_uint2(max(a.x,b.x), max(a.y,b.y));
}
inline __host__ __device__ uint3 max(uint3 a, uint3 b){
    return make_uint3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}
inline __host__ __device__ uint4 max(uint4 a, uint4 b){
    return make_uint4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}

////////////////////////////////////////////////////////////////////////////////
// lerp
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ float2 lerp(float2 a, float2 b, float t){
    return a + t*(b-a);
}
inline __device__ __host__ float3 lerp(float3 a, float3 b, float t){
    return a + t*(b-a);
}
inline __device__ __host__ float4 lerp(float4 a, float4 b, float t){
    return a + t*(b-a);
}

////////////////////////////////////////////////////////////////////////////////
// clamp
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ float2 clamp(float2 v, float a, float b){
    return make_float2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ float2 clamp(float2 v, float2 a, float2 b){
    return make_float2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ float3 clamp(float3 v, float a, float b){
    return make_float3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ float3 clamp(float3 v, float3 a, float3 b){
    return make_float3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ float4 clamp(float4 v, float a, float b){
    return make_float4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ float4 clamp(float4 v, float4 a, float4 b){
    return make_float4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}


inline __device__ __host__ int2 clamp(int2 v, int a, int b){
    return make_int2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ int2 clamp(int2 v, int2 a, int2 b){
    return make_int2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ int3 clamp(int3 v, int a, int b){
    return make_int3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ int3 clamp(int3 v, int3 a, int3 b){
    return make_int3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ int4 clamp(int4 v, int a, int b){
    return make_int4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ int4 clamp(int4 v, int4 a, int4 b){
    return make_int4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

inline __device__ __host__ uint2 clamp(uint2 v, uint a, uint b){
    return make_uint2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ uint2 clamp(uint2 v, uint2 a, uint2 b){
    return make_uint2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ uint3 clamp(uint3 v, uint a, uint b){
    return make_uint3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ uint3 clamp(uint3 v, uint3 a, uint3 b){
    return make_uint3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ uint4 clamp(uint4 v, uint a, uint b){
    return make_uint4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ uint4 clamp(uint4 v, uint4 a, uint4 b){
    return make_uint4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

////////////////////////////////////////////////////////////////////////////////
// dot product
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float dot(float2 a, float2 b){
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ float dot(float3 a, float3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ float dot(float4 a, float4 b){
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline __host__ __device__ int dot(int2 a, int2 b){
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ int dot(int3 a, int3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ int dot(int4 a, int4 b){
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline __host__ __device__ uint dot(uint2 a, uint2 b){
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ uint dot(uint3 a, uint3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ uint dot(uint4 a, uint4 b){
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

////////////////////////////////////////////////////////////////////////////////
// length
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float length(float2 v){
    return sqrtf(dot(v, v));
}
inline __host__ __device__ float length(float3 v){
    return sqrtf(dot(v, v));
}
inline __host__ __device__ float length(float4 r){
    return sqrtf(dot(r, r));
}

////////////////////////////////////////////////////////////////////////////////
// normalize
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float2 normalize(float2 v){
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ float3 normalize(float3 v){
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ float4 normalize(float4 v){
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}

////////////////////////////////////////////////////////////////////////////////
// floor
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float2 floor(const float2 v){
    return make_float2(floor(v.x), floor(v.y));
}
inline __host__ __device__ float3 floor(const float3 v){
    return make_float3(floor(v.x), floor(v.y), floor(v.z));
}
inline __host__ __device__ float4 floor(const float4 v){
    return make_float4(floor(v.x), floor(v.y), floor(v.z), floor(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// reflect
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float2 reflect(float2 i, float2 n){
	return i - 2.0f * n * dot(n,i);
}
inline __host__ __device__ float3 reflect(float3 i, float3 n){
	return i - 2.0f * n * dot(n,i);
}
inline __host__ __device__ float4 reflect(float4 i, float4 n){
	return i - 2.0f * n * dot(n,i);
}

////////////////////////////////////////////////////////////////////////////////
// absolute value
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float2 fabsf(float2 v){
	return make_float2(fabsf(v.x), fabsf(v.y));
}
inline __host__ __device__ float3 fabsf(float3 v){
	return make_float3(fabsf(v.x), fabsf(v.y), fabsf(v.z));
}
inline __host__ __device__ float4 fabsf(float4 v)
{
	return make_float4(fabsf(v.x), fabsf(v.y), fabsf(v.z), fabsf(v.w));
}

inline __host__ __device__ int2 abs(int2 v){
	return make_int2(abs(v.x), abs(v.y));
}
inline __host__ __device__ int3 abs(int3 v){
	return make_int3(abs(v.x), abs(v.y), abs(v.z));
}
inline __host__ __device__ int4 abs(int4 v)
{
	return make_int4(abs(v.x), abs(v.y), abs(v.z), abs(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// cross product
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ float3 cross(float3 a, float3 b)
{
    return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

////////////////////////////////////////////////////////////////////////////////
// Double2
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 make_double2(double s){
    return make_double2(s, s);
}

inline __host__ __device__ double2 make_double2(int2 a){
    return make_double2(double(a.x), double(a.y));
}

// negate
inline __host__ __device__ double2 operator-(double2 &a){
    return make_double2(-a.x, -a.y);
}

// addition
inline __host__ __device__ double2 operator+(double2 a, double2 b){
    return make_double2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(double2 &a, double2 b){
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ double2 operator+(double2 a, double b){
    return make_double2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(double2 &a, double b){
    a.x += b; a.y += b;
}
// subtract
inline __host__ __device__ double2 operator-(double2 a, double2 b)
{
    return make_double2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(double2 &a, double2 b)
{
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ double2 operator-(double2 a, double b)
{
    return make_double2(a.x - b, a.y - b);
}
inline __host__ __device__ void operator-=(double2 &a, double b)
{
    a.x -= b; a.y -= b;
}

// multiply
inline __host__ __device__ double2 operator*(double2 a, double2 b)
{
    return make_double2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(double2 &a, double2 b)
{
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ double2 operator*(double2 a, double s)
{
    return make_double2(a.x * s, a.y * s);
}
inline __host__ __device__ void operator*=(double2 &a, double s)
{
    a.x *= s; a.y *= s;
}
inline __host__ __device__ double2 operator*(double s, double2 a)
{
    return make_double2(a.x * s, a.y * s);
}

// divide
inline __host__ __device__ double2 operator/(double2 a, double2 b){
    return make_double2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ double2 operator/(double2 a, double s){
    double inv = 1.0 / s;
    return a * inv;
}
inline __host__ __device__ void operator/=(double2 &a, double s){
    double inv = 1.0 / s;
    a.x *= inv; a.y *= inv;
}
inline __host__ __device__ double2 operator/(double s, double2 a){
    return make_double2(s / a.x, s / a.y );
}

//lerp
inline __device__ __host__ double2 lerp(double2 a, double2 b, double t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ double clamp(double f, double a, double b){
    return fmax(a, fmin(f, b));
}
inline __device__ __host__ double2 clamp(double2 v, double a, double b){
    return make_double2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ double2 clamp(double2 v, double2 a, double2 b){
    return make_double2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}

// dot product
inline __host__ __device__ double dot(double2 a, double2 b)
{
    return a.x * b.x + a.y * b.y;
}

// length
inline __host__ __device__ double length(double2 v)
{
    return sqrt(dot(v, v));
}

// normalize
inline __host__ __device__ double2 normalize(double2 v)
{
    double invLen = 1.0 / sqrt(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ double2 floor(const double2 v)
{
    return make_double2(floor(v.x), floor(v.y));
}

// reflect
inline __host__ __device__ double2 reflect(double2 i, double2 n)
{
    return i - 2.0 * n * dot(n,i);
}

// absolute value
inline __host__ __device__ double2 fabs(double2 v){
	return make_double2(fabs(v.x), fabs(v.y));
}

#endif
