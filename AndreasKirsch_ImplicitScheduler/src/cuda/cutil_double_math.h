/*
 * cutil_d3math.h
 *
 *  Created on: Jun 24, 2011
 *      Author: andreas
 */
////////////////////////////////////////////////////////////////////////////////
// Essential functions on CUDA built in type
// Linh Ha lha@cs.utah.edu
// The vector / vector and vector * vector is point-wise
// Use it at your own risk
////////////////////////////////////////////////////////////////////////////////

#ifndef CUTIL_DOUBLE_MATH_H
#define CUTIL_DOUBLE_MATH_H

#include "cutil_math.h"

#ifndef __CUDACC__
#include <math.h>

inline double fmin(double a, double b){
    return a < b ? a : b;
}

inline double fmax(double a, double b){
    return a > b ? a : b;
}

inline double rsqrt(double x){
    return 1.0f / sqrt(x);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// additional constructors
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ double3 make_double3(double s){
    return make_double3(s, s, s);
}
inline __host__ __device__ double3 make_double3(float2 a){
    return make_double3(a.x, a.y, 0.0f);
}
inline __host__ __device__ double3 make_double3(float2 a, double s){
    return make_double3(a.x, a.y, s);
}

// change type
inline __host__ __device__ double3 make_double3(float3 a){
    return make_double3(double(a.x), double(a.y), double(a.z));
}
inline __host__ __device__ double3 make_double3(int3 a){
    return make_double3(double(a.x), double(a.y), double(a.z));
}
inline __host__ __device__ double3 make_double3(uint3 a){
    return make_double3(double(a.x), double(a.y), double(a.z));
}
// discards w
inline __host__ __device__ double3 make_double3(double4 a){
    return make_double3(a.x, a.y, a.z);
}

// additional constructors double4
inline __host__ __device__ double4 make_double4(double s){
    return make_double4(s, s, s, s);
}
inline __host__ __device__ double4 make_double4(double3 a){
    return make_double4(a.x, a.y, a.z, 0.0f);
}
inline __host__ __device__ double4 make_double4(double3 a, double w){
    return make_double4(a.x, a.y, a.z, w);
}

inline __host__ __device__ double4 make_double4(float4 a){
    return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));
}
inline __host__ __device__ double4 make_double4(int4 a){
    return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));
}
inline __host__ __device__ double4 make_double4(uint4 a){
    return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));
}

////////////////////////////////////////////////////////////////////////////////
// negate
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ double3 operator-(const double3 &a){
    return make_double3(-a.x, -a.y, -a.z);
}

inline __host__ __device__ double4 operator-(const double4 &a){
    return make_double4(-a.x, -a.y, -a.z, -a.w);
}

////////////////////////////////////////////////////////////////////////////////
// addition
////////////////////////////////////////////////////////////////////////////////

// addition 3
inline __host__ __device__ double3 operator+(double3 a, double3 b){
    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(double3 &a, double3 b){
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ double3 operator+(double3 a, double b){
    return make_double3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(double3 &a, double b){
    a.x += b; a.y += b; a.z += b;
}

inline __host__ __device__ double3 operator+(double b, double3 a){
    return make_double3(a.x + b, a.y + b, a.z + b);
}

// addition 4
inline __host__ __device__ double4 operator+(double4 a, double4 b){
    return make_double4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(double4 &a, double4 b){
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ double4 operator+(double4 a, double b){
    return make_double4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(double4 &a, double b){
    a.x += b; a.y += b; a.z += b; a.w += b;
}

inline __host__ __device__ double4 operator+(double b, double4 a){
    return make_double4(a.x + b, a.y + b, a.z + b,  a.w + b);
}

////////////////////////////////////////////////////////////////////////////////
// subtract
////////////////////////////////////////////////////////////////////////////////
// subtract 3
inline __host__ __device__ double3 operator-(double3 a, double3 b){
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(double3 &a, double3 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ double3 operator-(double3 a, double b){
    return make_double3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(double3 &a, double b){
    a.x -= b; a.y -= b; a.z -= b;
}

inline __host__ __device__ double3 operator-(double b, double3 a){
    return make_double3(b - a.x, b - a.y, b - a.z);
}

// subtract 4
inline __host__ __device__ double4 operator-(double4 a, double4 b){
    return make_double4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(double4 &a, double4 b){
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ double4 operator-(double4 a, double b){
    return make_double4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ void operator-=(double4 &a, double b){
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}

inline __host__ __device__ double4 operator-(double b, double4 a){
    return make_double4(b - a.x, b - a.y, b - a.z, b - a.w);
}

////////////////////////////////////////////////////////////////////////////////
// Multiply
////////////////////////////////////////////////////////////////////////////////
// Multiply 3
inline __host__ __device__ double3 operator*(double3 a, double3 b){
    return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(double3 &a, double3 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ double3 operator*(double3 a, double b){
    return make_double3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ void operator*=(double3 &a, double b){
    a.x *= b; a.y *= b; a.z *= b;
}

inline __host__ __device__ double3 operator*(double b, double3 a){
    return make_double3(b * a.x, b * a.y, b * a.z);
}

//multiply 4
inline __host__ __device__ double4 operator*(double4 a, double4 b){
    return make_double4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(double4 &a, double4 b){
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ double4 operator*(double4 a, double b){
    return make_double4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ void operator*=(double4 &a, double b){
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __host__ __device__ double4 operator*(double b, double4 a){
    return make_double4(b * a.x, b * a.y, b * a.z, b * a.w);
}

////////////////////////////////////////////////////////////////////////////////
// Divide
////////////////////////////////////////////////////////////////////////////////
// divide 3
inline __host__ __device__ double3 operator/(double3 a, double3 b){
    return make_double3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ void operator/=(double3 &a, double3 b){
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}
inline __host__ __device__ double3 operator/(double3 a, double b){
    double inv = 1.f / b;
    return make_double3(a.x * inv, a.y * inv, a.z * inv);
}
inline __host__ __device__ void operator/=(double3 &a, double b){
    double inv = 1.f / b;
    a.x *= inv; a.y *= inv; a.z *= inv;
}
inline __host__ __device__ double3 operator/(double b, double3 a){
    return make_double3(b / a.x, b / a.y, b / a.z);
}

// divide 4
inline __host__ __device__ double4 operator/(double4 a, double4 b){
    return make_double4(a.x / b.x, a.y / b.y, a.z / b.z,  a.w / b.w);
}
inline __host__ __device__ void operator/=(double4 &a, double4 b){
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}
inline __host__ __device__ double4 operator/(double4 a, double b){
    double inv = 1.f / b;
    return make_double4(a.x * inv, a.y * inv, a.z * inv,  a.w * inv);
}
inline __host__ __device__ void operator/=(double4 &a, double b){
    double inv = 1.f / b;
    a.x *= inv; a.y *= inv; a.z *= inv; a.w *= inv;
}
inline __host__ __device__ double4 operator/(double b, double4 a){
    return make_double4(b / a.x, b / a.y, b / a.z, b / a.w);
}

////////////////////////////////////////////////////////////////////////////////
// min
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 fmin(double3 a, double3 b){
	return make_double3(fmin(a.x,b.x), fmin(a.y,b.y), fmin(a.z,b.z));
}
inline  __host__ __device__ double4 fmin(double4 a, double4 b){
	return make_double4(fmin(a.x,b.x), fmin(a.y,b.y), fmin(a.z,b.z), fmin(a.w,b.w));
}

////////////////////////////////////////////////////////////////////////////////
// max
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 fmax(double3 a, double3 b){
	return make_double3(fmax(a.x,b.x), fmax(a.y,b.y), fmax(a.z,b.z));
}
inline __host__ __device__ double4 fmax(double4 a, double4 b){
	return make_double4(fmax(a.x,b.x), fmax(a.y,b.y), fmax(a.z,b.z), fmax(a.w,b.w));
}

////////////////////////////////////////////////////////////////////////////////
// lerp
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ double3 lerp(double3 a, double3 b, double t){
    return a + t*(b-a);
}
inline __device__ __host__ double4 lerp(double4 a, double4 b, double t){
    return a + t*(b-a);
}

////////////////////////////////////////////////////////////////////////////////
// clamp
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ double3 clamp(double3 v, double a, double b){
    return make_double3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ double3 clamp(double3 v, double3 a, double3 b){
    return make_double3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ double4 clamp(double4 v, double a, double b){
    return make_double4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ double4 clamp(double4 v, double4 a, double4 b){
    return make_double4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}
////////////////////////////////////////////////////////////////////////////////
// dot product
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double dot(double3 a, double3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ double dot(double4 a, double4 b){
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

////////////////////////////////////////////////////////////////////////////////
// length
////////////////////////////////////////////////////////////////////////////////

inline __host__ __device__ double length(double3 v){
    return sqrt(dot(v, v));
}
inline __host__ __device__ double length(double4 r){
    return sqrt(dot(r, r));
}

////////////////////////////////////////////////////////////////////////////////
// normalize
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 normalize(double3 v){
    double invLen = rsqrt(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ double4 normalize(double4 v){
    double invLen = rsqrt(dot(v, v));
    return v * invLen;
}

////////////////////////////////////////////////////////////////////////////////
// floor
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 floor(const double3 v){
    return make_double3(floor(v.x), floor(v.y), floor(v.z));
}
inline __host__ __device__ double4 floor(const double4 v){
    return make_double4(floor(v.x), floor(v.y), floor(v.z), floor(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// reflect
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 reflect(double3 i, double3 n){
	return i - 2.0f * n * dot(n,i);
}
inline __host__ __device__ double4 reflect(double4 i, double4 n){
	return i - 2.0f * n * dot(n,i);
}

////////////////////////////////////////////////////////////////////////////////
// absolute value
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 fabs(double3 v){
	return make_double3(fabs(v.x), fabs(v.y), fabs(v.z));
}
inline __host__ __device__ double4 fabs(double4 v)
{
	return make_double4(fabs(v.x), fabs(v.y), fabs(v.z), fabs(v.w));
}

////////////////////////////////////////////////////////////////////////////////
// cross product
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 cross(double3 a, double3 b)
{
    return make_double3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

#endif /* CUTIL_DOUBLE_MATH_H */
