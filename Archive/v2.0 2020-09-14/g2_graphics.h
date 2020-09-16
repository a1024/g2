//g2_graphics.h - Include for G2 graphics API.
//Copyright (C) 2012-2020  Ayman Wagih Mohsen
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include		"g2_common.h"
#include		"g2_sse2.h"//lighten_sse2() for dim_screen()
#include		<Windows.h>
#include		<gl/gl.h>
#include		<gl/glu.h>
#include		<map>
#include		<list>
extern int		w, h, X0, Y0;
extern HDC		ghDC;
extern bool		usingOpenGL;
extern long long broken;

extern int		fontH, *fontW;
extern bool		SSE4_1;
extern int		simd_method;

//vector algebra:
//extern const float _pi, _2pi, pi_2, inv_2pi, sqrt2, torad, infinity, inv255, inv256, inv128;
template<typename Type>struct	tvec2
{
	Type x, y;
	tvec2():x(0), y(0){}
	tvec2(Type x, Type y):x(x), y(y){}
	void set(Type x, Type y){this->x=x, this->y=y;}
	template<typename T2>tvec2<Type>& operator=(tvec2<T2> const &v){x=(Type)v.x, y=(Type)v.y;return *this;}
	tvec2& operator+=(tvec2 const &b){x+=b.x, y+=b.y; return *this;}
	tvec2& operator-=(tvec2 const &b){x-=b.x, y-=b.y; return *this;}
	tvec2& operator+=(Type x){this->x+=x, y+=x; return *this;}
	tvec2& operator-=(Type x){this->x-=x, y-=x; return *this;}
	tvec2& operator*=(Type x){this->x*=x, y*=x; return *this;}
	tvec2& operator/=(Type x){this->x/=x, y/=x; return *this;}
	Type dot(tvec2 const &other)const{return x*other.x+y*other.y;}
	Type cross(tvec2 const &other)const{return x*other.y-y*other.x;}
	Type magnitude()const{return sqrt(x*x+y*y);}
	Type mag_sq()const{return x*x+y*y;}
	Type angle()const{return atan(y/x);}
	Type angle2()const{return atan2(y, x);}
};
template<typename Type>inline tvec2<Type>	operator*(tvec2<Type> const &p, Type x){return tvec2<Type>(p.x*x, p.y*x);}
template<typename Type>inline tvec2<Type>	operator*(Type x, tvec2<Type> const &p){return tvec2<Type>(p.x*x, p.y*x);}
template<typename Type>inline tvec2<Type>	operator/(tvec2<Type> const &p, Type x){return tvec2<Type>(p.x/x, p.y/x);}
template<typename Type>inline tvec2<Type>	operator+(tvec2<Type> const &a, tvec2<Type> const &b){return tvec2<Type>(a.x+b.x, a.y+b.y);}
template<typename Type>inline tvec2<Type>	operator-(tvec2<Type> const &a, tvec2<Type> const &b){return tvec2<Type>(a.x-b.x, a.y-b.y);}
template<typename Type>inline tvec2<Type>	operator+(tvec2<Type> const &p, Type x){return tvec2<Type>(p.x+x, p.y+x);}
template<typename Type>inline tvec2<Type>	operator-(tvec2<Type> const &p, Type x){return tvec2<Type>(p.x-x, p.y-x);}
template<typename Type>inline bool			operator==(tvec2<Type> const &a, tvec2<Type> const &b){return a.x==b.x&&a.y==b.y;}
template<typename Type>inline bool			operator!=(tvec2<Type> const &a, tvec2<Type> const &b){return a.x!=b.x||a.y!=b.y;}
template<typename Type>inline tvec2<Type>	operator-(tvec2<Type> const &p){return tvec2(-p.x, -p.y);}
template<typename Type>struct	tmat2
{
	Type a, b,//row 0
		c, d;//row 1
	tmat2():a(0), b(0), c(0), d(0){}
	tmat2(Type a, Type b, Type c, Type d):a(a), b(b), c(c), d(d){}
	void set(Type a, Type b, Type c, Type d){this->a=a, this->b=b, this->c=c, this->d=d;}
	void set_v(tvec2<Type> &col1, tvec2<Type> &col2)
	{
		a=col1.x, b=col2.x;
		c=col1.y, d=col2.y;
	}
};
template<typename Type>inline tvec2<Type>	operator*(tmat2<Type> const &m, tvec2<Type> const &v){return tvec2<Type>(m.a*v.x+m.b*v.y, m.c*v.x+m.d*v.y);}
template<typename Type>struct	tvec3
{
	Type x, y, z;
	tvec3():x(0), y(0), z(0){}
	tvec3(Type x, Type y, Type z):x(x), y(y), z(z){}
	tvec3(Type gain):x(gain), y(gain), z(gain){}
	tvec3(Type *p):x(p[0]), y(p[1]), z(p[2]){}
	template<typename T2>tvec3(tvec3<T2> const &v):x((Type)v.x), y((Type)v.y), z((Type)v.z){}
	void set(Type x, Type y, Type z){this->x=x, this->y=y, this->z=z;}
	Type& operator[](int idx){return (&x)[idx];}
	Type operator[](int idx)const{return (&x)[idx];}
	template<typename T2>tvec3<Type>& operator=(tvec3<T2> const &v){x=(Type)v.x, y=(Type)v.y, z=(Type)v.z;return *this;}
	tvec3& operator+=(tvec3 const &b){x+=b.x, y+=b.y, z+=b.z; return *this;}
	tvec3& operator-=(tvec3 const &b){x-=b.x, y-=b.y, z-=b.z; return *this;}
	tvec3& operator+=(Type x){this->x+=x, y+=x, z+=x; return *this;}
	tvec3& operator-=(Type x){this->x-=x, y-=x, z-=x; return *this;}
	tvec3& operator*=(Type x){this->x*=x, y*=x, z*=x; return *this;}
	tvec3& operator/=(Type x){this->x/=x, y/=x, z/=x; return *this;}
	template<typename T2>operator tvec3<T2>()
	{
		return tvec3<T2>((T2)x, (T2)y, (T2)z);
	}
	Type dot(tvec3 const &other)const{return x*other.x+y*other.y+z*other.z;}
	tvec3 cross(tvec3 const &other)const{return tvec3<Type>(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x);}
	tvec3 triple_product(tvec3 const &b, tvec3 const &c)const;
	Type magnitude()const{return sqrt(x*x+y*y+z*z);}
	Type mag_sq()const{return x*x+y*y+z*z;}
	Type theta()const{return atan(z/sqrt(x*x+y*y));}//vertical angle
	//Type theta2()const{return atan2(y, x);}
	Type phi(){return atan(y/x);}//horizontal angle
	Type phi2(){return atan2(y, x);}
	bool isnan(){return x!=x||y!=y||z!=z;}
	bool isnan_or_inf(){return x!=x||y!=y||z!=z||abs(x)==infinity||abs(y)==infinity||abs(z)==infinity;}
};
template<typename Type>inline tvec3<Type>	operator*(tvec3<Type> const &p, Type x){return tvec3<Type>(p.x*x, p.y*x, p.z*x);}
template<typename Type>inline tvec3<Type>	operator*(Type x, tvec3<Type> const &p){return tvec3<Type>(p.x*x, p.y*x, p.z*x);}
template<typename Type>inline tvec3<Type>	operator/(tvec3<Type> const &p, Type x){return tvec3<Type>(p.x/x, p.y/x, p.z/x);}
template<typename Type>inline tvec3<Type>	operator+(tvec3<Type> const &a, tvec3<Type> const &b){return tvec3<Type>(a.x+b.x, a.y+b.y, a.z+b.z);}
template<typename T1, typename T2>inline tvec3<T1> operator-(tvec3<T1> const &a, tvec3<T2> const &b){return tvec3<T1>(a.x-(T2)b.x, a.y-(T2)b.y, a.z-(T2)b.z);}
template<typename Type>inline tvec3<Type>	operator+(tvec3<Type> const &p, Type x){return tvec3<Type>(p.x+x, p.y+x, p.z+x);}
template<typename Type>inline tvec3<Type>	operator-(tvec3<Type> const &p, Type x){return tvec3<Type>(p.x-x, p.y-x, p.z-x);}
template<typename Type>inline bool			operator==(tvec3<Type> const &a, tvec3<Type> const &b)
{
	Type tolerance=(Type)1e-10;
	//Type sd=abs(a.x-b.x)+abs(a.y-b.y)+abs(a.z-b.z);
	//if(sd<(Type)1e-6&&sd>tolerance)
	//	int LOL_1=0;
	return abs(a.x-b.x)+abs(a.y-b.y)+abs(a.z-b.z)<tolerance;
//	return abs(a.x-b.x)<tolerance&&abs(a.y-b.y)<tolerance&&abs(a.z-b.z)<tolerance;
//	return a.x==b.x&&a.y==b.y&&a.z==b.z;
}
template<typename Type>inline bool			operator!=(tvec3<Type> const &a, tvec3<Type> const &b){return a.x!=b.x||a.y!=b.y||a.z!=b.z;}
template<typename Type>inline tvec3<Type>	operator-(tvec3<Type> const &p){return tvec3<Type>(-p.x, -p.y, -p.z);}
template<typename Type>inline tvec3<Type>	tvec3<Type>::triple_product(tvec3<Type> const &b, tvec3<Type> const &c)const{return this->dot(c)*b-this->dot(b)*c;}
template<typename Type>inline tvec3<Type>	normalize(tvec3<Type> const &v){Type invm=1/v.magnitude(); return invm*v;}
typedef tvec2<float> vec2;
typedef tvec2<double> dvec2;
typedef tmat2<float> mat2;
typedef tmat2<double> dmat2;
typedef tvec3<float> vec3;
typedef tvec3<double> dvec3;

inline bool is_nan_or_inf(double x){return (((int*)&x)[1]&0x7FF00000)==0x7FF00000;}
inline int count_nan_inf(dvec3 const &v){return is_nan_or_inf(v.x)+is_nan_or_inf(v.y)+is_nan_or_inf(v.z);}

//vector algebra for OpenGL:
struct			mat3
{
	vec3	c[3];
	mat3(){}
	mat3(float gain){memset(c, 0, 9<<2), c[0][0]=c[1][1]=c[2][2]=gain;}
	mat3(vec3 const &c0, vec3 const &c1, vec3 const &c2){c[0]=c0, c[1]=c1, c[2]=c2;}
	float* data(){return &c[0][0];}
};
struct			vec4
{
	union
	{
		struct{float x, y, z, w;};
		struct{float r, g, b, a;};
	};
	vec4(){memset(&x, 0, 4<<2);}
//	vec4(float gain):x(gain), y(gain), z(gain), w(gain){}
	vec4(float x, float y, float z, float w):x(x), y(y), z(z), w(w){}
	vec4(vec3 const &v, float w):x(v.x), y(v.y), z(v.z), w(w){}
	vec4(__m128 const &v){_mm_storeu_ps(&x, v);}
	float& operator[](int idx){return (&x)[idx];}
	float operator[](int idx)const{return (&x)[idx];}
	operator __m128(){return _mm_loadu_ps(&x);}
	operator vec3(){return vec3(x, y, z);}
	void set(float x, float y, float z, float w){this->x=x, this->y=y, this->z=z, this->w=w;}
	void setzero(){_mm_storeu_ps(&x, _mm_setzero_ps());}
	float dot(vec4 const &other)
	{
		__m128 a=_mm_loadu_ps(&x), b=_mm_loadu_ps(&other.x);
		a=_mm_mul_ps(a, b);
		a=_mm_hadd_ps(a, a);
		a=_mm_hadd_ps(a, a);
		float result;
		_mm_store_ss(&result, a);//does not need to be aligned
		return result;
	}
};
inline vec4		operator+(vec4 const &a, vec4 const &b)
{
	__m128 va=_mm_loadu_ps(&a.x), vb=_mm_loadu_ps(&b.x);
	va=_mm_add_ps(va, vb);
	vec4 r;
	_mm_storeu_ps(&r.x, va);
	return r;
}
inline vec4		operator*(vec4 const &v, float s)
{
	__m128 vv=_mm_loadu_ps(&v.x), vs=_mm_set1_ps(s);
	vv=_mm_mul_ps(vv, vs);
	vec4 r;
	_mm_storeu_ps(&r.x, vv);
	return r;
}
inline vec4		operator*(float s, vec4 const &v){return v*s;}
struct			mat4//column major
{
	vec4	c[4];
	mat4(){}
	mat4(const float *v){memcpy(c, v, 16<<2);}
	mat4(float gain){memset(c, 0, 16<<2), c[0][0]=c[1][1]=c[2][2]=c[3][3]=gain;}
	mat4(vec4 const &c0, vec4 const &c1, vec4 const &c2, vec4 const &c3){c[0]=c0, c[1]=c1, c[2]=c2, c[3]=c3;}
	mat4(
		float a11, float a12, float a13, float a14,
		float a21, float a22, float a23, float a24,
		float a31, float a32, float a33, float a34,
		float a41, float a42, float a43, float a44)
	{
		c[0].set(a11, a21, a31, a41);
		c[1].set(a12, a22, a32, a42);
		c[2].set(a13, a23, a33, a43);
		c[3].set(a14, a24, a34, a44);
	}
	void set(
		float a11, float a12, float a13, float a14,
		float a21, float a22, float a23, float a24,
		float a31, float a32, float a33, float a34,
		float a41, float a42, float a43, float a44)
	{
		c[0].set(a11, a21, a31, a41);
		c[1].set(a12, a22, a32, a42);
		c[2].set(a13, a23, a33, a43);
		c[3].set(a14, a24, a34, a44);
	}
	operator mat3(){return mat3((vec3)c[0], (vec3)c[1], (vec3)c[2]);}
	float* data(){return &c[0][0];}
	void setzero(){memset(c, 0, 16<<2);}
	//void setrow(int idx, vec4 const &v){float *p=&c[0][0]+idx; p[0]=v.x, p[4]=v.y, p[8]=v.z, p[12]=v.w;}
	vec4& operator[](int idx){return c[idx];}
	vec4 operator[](int idx)const{return c[idx];}
};
inline mat4		transpose(mat4 const &m)
{
	__m128 r0=_mm_loadu_ps(&m[0][0]), r1=_mm_loadu_ps(&m[1][0]), r2=_mm_loadu_ps(&m[2][0]), r3=_mm_loadu_ps(&m[3][0]);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	mat4 m2;
	_mm_storeu_ps(&m2[0][0], r0);
	_mm_storeu_ps(&m2[1][0], r1);
	_mm_storeu_ps(&m2[2][0], r2);
	_mm_storeu_ps(&m2[3][0], r3);
	return m2;
}
inline vec4		operator*(mat4 const &m, vec4 const &v)
{
	__m128 c0=_mm_loadu_ps(&m[0][0]), c1=_mm_loadu_ps(&m[1][0]), c2=_mm_loadu_ps(&m[2][0]), c3=_mm_loadu_ps(&m[3][0]);
	//	vv=_mm_loadu_ps(&v.x);
	c0=_mm_mul_ps(c0, _mm_set1_ps(v.x));
	c1=_mm_mul_ps(c1, _mm_set1_ps(v.y));
	c2=_mm_mul_ps(c2, _mm_set1_ps(v.z));
	c3=_mm_mul_ps(c3, _mm_set1_ps(v.w));
	c0=_mm_add_ps(c0, c1);
	c0=_mm_add_ps(c0, c2);
	c0=_mm_add_ps(c0, c3);
	vec4 r;
	_mm_storeu_ps(&r[0], c0);
	return r;
}
inline mat4		operator*(mat4 const &a, mat4 const &b)
{
	mat4 c=transpose(a);
	float d[]=
	{
		c[0].dot(b[0]), c[0].dot(b[1]), c[0].dot(b[2]), c[0].dot(b[3]),
		c[1].dot(b[0]), c[1].dot(b[1]), c[1].dot(b[2]), c[1].dot(b[3]),
		c[2].dot(b[0]), c[2].dot(b[1]), c[2].dot(b[2]), c[2].dot(b[3]),
		c[3].dot(b[0]), c[3].dot(b[1]), c[3].dot(b[2]), c[3].dot(b[3])
	};
	return mat4(
		vec4(c[0].dot(b[0]), c[1].dot(b[0]), c[2].dot(b[0]), c[3].dot(b[0])),
		vec4(c[0].dot(b[1]), c[1].dot(b[1]), c[2].dot(b[1]), c[3].dot(b[1])),
		vec4(c[0].dot(b[2]), c[1].dot(b[2]), c[2].dot(b[2]), c[3].dot(b[2])),
		vec4(c[0].dot(b[3]), c[1].dot(b[3]), c[2].dot(b[3]), c[3].dot(b[3])));
}
inline mat4		translate(mat4 const &m, vec3 const &delta)//from glm
{
	vec4 v2(delta, 1);
	mat4 r=m;
	r[3]=m[0]*v2[0]+m[1]*v2[1]+m[2]*v2[2]+m[3];
	return r;
}
inline mat4		rotate(mat4 const &m, float angle, vec3 const &dir)//from glm
{
	float ca=cos(angle), sa=sin(angle);
	vec3 axis=normalize(dir), temp=(1-ca)*axis;
	mat4 rotate(
		vec4(ca+temp[0]*axis[0],				temp[0]*axis[1]+sa*axis[2],		temp[0]*axis[2]-sa*axis[1],	0),//col 1
		vec4(	temp[1]*axis[0]-sa*axis[2],	ca+	temp[1]*axis[1],				temp[1]*axis[2]+sa*axis[0],	0),
		vec4(	temp[2]*axis[0]+sa*axis[1],		temp[2]*axis[1]-sa*axis[0],	ca+	temp[2]*axis[2],			0),
		vec4(0, 0, 0, 1));//col 4
	return m*rotate;
}
inline mat4		scale(mat4 const &m, vec3 const &ammount)
{
	mat4 r(m[0]*ammount.x, m[1]*ammount.y, m[2]*ammount.z, m[3]);
	return r;
}
inline mat4		lookAt(vec3 const &cam, vec3 const &center, vec3 const &up)//from glm
{
	vec3 f=normalize(center-cam),
		u=normalize(up), s=normalize(f.cross(u));
	u=s.cross(f);
	mat4 r(
		vec4(s, -s.dot(cam)),
		vec4(u, -u.dot(cam)),
		vec4(-f, f.dot(cam)),
		vec4(0, 0, 0, 1));
	r=transpose(r);
	return r;
}
inline mat4		matrixFPSViewRH(vec3 const &_cam, float pitch, float yaw)
{//https://www.3dgep.com/understanding-the-view-matrix/
	vec3 cam(_cam.y, _cam.z, _cam.x);//why yzx?
	float cos_p=cos(pitch), sin_p=sin(pitch), cos_y=cos(yaw), sin_y=sin(yaw);
	vec3
			xaxis(cos_y, 0, -sin_y),
			yaxis(sin_y*sin_p, cos_p, cos_y*sin_p),
			zaxis(sin_y*cos_p, -sin_p, cos_p*cos_y);
	
	return mat4(
		vec4(xaxis.z, yaxis.z, zaxis.z, 0),//why zxy?
		vec4(xaxis.x, yaxis.x, zaxis.x, 0),
		vec4(xaxis.y, yaxis.y, zaxis.y, 0),
		vec4(-xaxis.dot(cam), -yaxis.dot(cam), -zaxis.dot(cam), 1));
}
inline mat4		perspective(float tanfov, float ar, float znear, float zfar)
{
	return mat4(
		vec4(1/tanfov, 0, 0, 0),
		vec4(0, ar/tanfov, 0, 0),
		vec4(0, 0, -(zfar+znear)/(zfar-znear), -1),
		vec4(0, 0, -2*zfar*znear/(zfar-znear), 0));
}
struct			ivec4
{
	union
	{
		struct{int x, y, z, w;};
		struct{int x1, y1, dx, dy;};
	};
	ivec4(){_mm_storeu_si128((__m128i*)&x, _mm_setzero_si128());}
	ivec4(int x, int y, int z, int w):x(x), y(y), z(z), w(w){}
	void set(int x, int y, int z, int w){this->x=x, this->y=y, this->z=z, this->w=w;}
};
inline mat4		GetTransformInverseNoScale(const mat4 &inM)// Requires this matrix to be transform matrix, NoScale version requires this matrix be of scale 1
{
#define MakeShuffleMask(x,y,z,w)           (x | (y<<2) | (z<<4) | (w<<6))

// vec(0, 1, 2, 3) -> (vec[x], vec[y], vec[z], vec[w])
#define VecSwizzleMask(vec, mask)          _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(vec), mask))
#define VecSwizzle(vec, x, y, z, w)        VecSwizzleMask(vec, MakeShuffleMask(x,y,z,w))
#define VecSwizzle1(vec, x)                VecSwizzleMask(vec, MakeShuffleMask(x,x,x,x))
// special swizzle
#define VecSwizzle_0022(vec)               _mm_moveldup_ps(vec)
#define VecSwizzle_1133(vec)               _mm_movehdup_ps(vec)

// return (vec1[x], vec1[y], vec2[z], vec2[w])
#define VecShuffle(vec1, vec2, x,y,z,w)    _mm_shuffle_ps(vec1, vec2, MakeShuffleMask(x,y,z,w))
// special shuffle
#define VecShuffle_0101(vec1, vec2)        _mm_movelh_ps(vec1, vec2)
#define VecShuffle_2323(vec1, vec2)        _mm_movehl_ps(vec2, vec1)
	mat4 r;

	// transpose 3x3, we know m03 = m13 = m23 = 0
	__m128 t0 = VecShuffle_0101(inM[0], inM[1]); // 00, 01, 10, 11
	__m128 t1 = VecShuffle_2323(inM[0], inM[1]); // 02, 03, 12, 13
	r[0] = VecShuffle(t0, inM[2], 0,2,0,3); // 00, 10, 20, 23(=0)
	r[1] = VecShuffle(t0, inM[2], 1,3,1,3); // 01, 11, 21, 23(=0)
	r[2] = VecShuffle(t1, inM[2], 0,2,2,3); // 02, 12, 22, 23(=0)

	// last line
	r[3] =                  _mm_mul_ps(r[0], VecSwizzle1(inM[3], 0));
	r[3] = _mm_add_ps(r[3], _mm_mul_ps(r[1], VecSwizzle1(inM[3], 1)));
	r[3] = _mm_add_ps(r[3], _mm_mul_ps(r[2], VecSwizzle1(inM[3], 2)));
	r[3] = _mm_sub_ps(_mm_setr_ps(0.f, 0.f, 0.f, 1.f), r[3]);

#undef MakeShuffleMask
#undef VecSwizzleMask
#undef VecSwizzle
#undef VecSwizzle1
#undef VecSwizzle_0022
#undef VecSwizzle_1133
#undef VecShuffle
#undef VecShuffle_0101
#undef VecShuffle_2323

	return r;
}
inline mat3		normalMatrix(mat4 const &m)//inverse transpose of top left 3x3 submatrix
{
	mat4 r=GetTransformInverseNoScale(m);
	return (mat3)transpose(r);
}

struct			Camera
{
	dvec3 p0;
	dvec2 a0;
	double tanfov0, dcam0, da0, mouse_sensitivity;
	//const dvec3 p0;
	//const dvec2 a0;
	//const double tanfov0, dcam0, da0, mouse_sensitivity;
	dvec3		p;//position
	dvec2		a;//ax: yaw, ay: pitch
	double		tanfov, dcam, da, da_tfov, cax, sax, cay, say;
	Camera():p0(5, 5, 5), a0(225*G2::_pi_180, 324.7356103172454f*G2::_pi_180), tanfov0(1), dcam0(0.04), da0(2*G2::_pi_180), mouse_sensitivity(0.003),
		p(p0), a(a0), dcam(dcam0), tanfov(tanfov0), da(da0), da_tfov(tanfov), cax(cos(a.x)), sax(sin(a.x)), cay(cos(a.y)), say(sin(a.y)){}
	Camera(double x, double y, double z, double ax, double ay, double tanfov):
		p0(x, y, z), a0(ax, ay), tanfov0(tanfov), dcam0(0.04), da0(2*G2::_pi_180), mouse_sensitivity(0.003),
		p(x, y, z), a(ax, ay), tanfov(tanfov), dcam(dcam0), da(da0), da_tfov(tanfov), cax(cos(a.x)), sax(sin(a.x)), cay(cos(a.y)), say(sin(a.y)){}
	void set(double x, double y, double z, double ax, double ay, double tanfov)
	{
		p0.set(x, y, z), a0.set(ax, ay), tanfov0=tanfov, dcam0=0.04, da0=2*G2::_pi_180, mouse_sensitivity=0.003;
		p.set(x, y, z), a.set(ax, ay), this->tanfov=tanfov, dcam=dcam0, da=da0, da_tfov=tanfov, cax=cos(a.x), sax=sin(a.x), cay=cos(a.y), say=sin(a.y);
	}
	void moveFastForward(){p.x+=10*dcam*cax*cay,	p.y+=10*dcam*sax*cay,	p.z+=10*dcam*say;}
	void moveFastLeft	(){p.x-=10*dcam*sax,		p.y+=10*dcam*cax;}
	void moveFastBack	(){p.x-=10*dcam*cax*cay,	p.y-=10*dcam*sax*cay,	p.z-=10*dcam*say;}
	void moveFastRight	(){p.x+=10*dcam*sax,		p.y-=10*dcam*cax;}
	void moveFastUp		(){p.z+=10*dcam;}
	void moveFastDown	(){p.z-=10*dcam;}
	void moveForward	(){p.x+=dcam*cax*cay,	p.y+=dcam*sax*cay,	p.z+=dcam*say;}
	void moveLeft		(){p.x-=dcam*sax,		p.y+=dcam*cax;}
	void moveBack		(){p.x-=dcam*cax*cay,	p.y-=dcam*sax*cay,	p.z-=dcam*say;}
	void moveRight		(){p.x+=dcam*sax,		p.y-=dcam*cax;}
	void moveUp			(){p.z+=dcam;}
	void moveDown		(){p.z-=dcam;}
	void turnUp			(){update_angle(a.y+=da_tfov*da, cay, say);}
	void turnDown		(){update_angle(a.y-=da_tfov*da, cay, say);}
	void turnLeft		(){update_angle(a.x+=da_tfov*da, cax, sax);}
	void turnRight		(){update_angle(a.x-=da_tfov*da, cax, sax);}
	//void turnMouse		(unsigned long lParam)
	//{
	//	update_angle(ax+=mouse_sensitivity*(X0-((short*)&lParam)[0]), cax, sax);
	//	update_angle(ay+=mouse_sensitivity*(Y0-((short*)&lParam)[1]), cay, say);
	//}
	void turnMouse		(long lParam)
	{
		short dx=(short&)lParam, dy=((short*)&lParam)[1];
		a.x+=mouse_sensitivity*da_tfov*(X0-dx), update_angle(a.x, cax, sax);
		a.y+=mouse_sensitivity*da_tfov*(Y0-dy), update_angle(a.y, cay, say);
	}
	void zoomIn			(){tanfov/=1.1f, da_tfov=tanfov>1?1:tanfov;}
	void zoomOut		(){tanfov*=1.1f, da_tfov=tanfov>1?1:tanfov;}
	void reset(){p=p0, a=a0, tanfov=tanfov0,	dcam=dcam0, da_tfov=tanfov;}
	void teleport(dvec3 const &p, double ax, double ay, double tanfov)
	{
		this->p=p;
		update_angle(a.x=ax, cax, sax);
		update_angle(a.y=ay, cay, say);
		this->tanfov=tanfov;
		da_tfov=tanfov>1?1:tanfov;
	}
	void faster(double A){dcam*=A;}
	void faster(){dcam*=2;}
	void slower(){dcam*=0.5;}

	//vertex conversions - double precision
	void relworld2camscreen(dvec3 const &d, dvec3 &cp, dvec2 &s)const
	{
		double cpt=d.x*cax+d.y*sax;
		cp.x=d.x*sax-d.y*cax, cp.y=cpt*say-d.z*cay, cp.z=cpt*cay+d.z*say;
		cpt=X0/(cp.z*tanfov), s.x=X0+cp.x*cpt, s.y=Y0+cp.y*cpt;
	}
	void world2camscreen(dvec3 const &p_world, dvec3 &cp, dvec2 &s)const{relworld2camscreen(p_world-p, cp, s);}
	void relworld2cam(dvec3 const &d, dvec3 &cp)const
	{
		double temp=d.x*cax+d.y*sax;
		cp.x=d.x*sax-d.y*cax, cp.y=temp*say-d.z*cay, cp.z=temp*cay+d.z*say;
	}
	void world2cam(dvec3 const &p, dvec3 &cp)const{relworld2cam(p-this->p, cp);}
	void cam2screen(dvec3 const &cp, dvec2 &s)const
	{
		double temp=X0/(cp.z*tanfov);
		s.set(X0+cp.x*temp, Y0+cp.y*temp);
	}

	//vertex conversions - single precision
	void world2camscreen(dvec3 const &p_world, vec3 &cp, vec2 &s)const
	{
		dvec3 dcp;
		dvec2 ds;
		relworld2camscreen(p_world-p, dcp, ds);
		cp=dcp, s=ds;
	}

	void scaleXabout(double const &VX, double A){p.x=VX+(p.x-VX)*A;}
	void scaleYabout(double const &VY, double A){p.y=VY+(p.y-VY)*A;}
	void scaleZabout(double const &VZ, double A){p.z=VZ+(p.z-VZ)*A;}
	void scaleXYZabout(double const &VX, double const &VY, double const &VZ, double A)
	{
		scaleXabout(VX, A);
		scaleYabout(VY, A);
		scaleZabout(VZ, A);
	}
};

struct			Bitmap
{
	int			w, h, *rgb;
	HBITMAP		hBitmap;
	Bitmap():w(0), h(0), rgb(0), hBitmap(0){}
	void		set(int w, int h);
//	void		set(int w, int h, int *&rgb);
	void		resize(int w, int h);
	void		finish();
	void		use();
	void		drop();
};
extern Bitmap	gBitmap;

struct			Region
{
	static const Region *current;
	HRGN__		*hRgn;
//	ivec4		p;
	int bx1, bx2, by1, by2, bw, bh, X0, Y0;
//	Region(int x1, int y1, int x2, int y2);
//	~Region();
	void		create(int x1, int y1, int x2, int y2);
	void		destroy();
	void		use();
	void		drop();
};

void			capture_window(int *rgb);
void			display_texture_fullwindow(int *rgb);

//text API:
double			changeFont(unsigned wParam);
double			largerFont();
double			smallerFont();
double			setFont(int newFont);
void			selectFont();
void			deselectFont();
int				getTextWidth(const char *a, int length);
int				getTextWidth(const char *a, int i, int f);

int				getBkMode();
int				setBkMode(int mode);//TODO: merge mode & color
int				getBkColor();
int				setBkColor(int color);
int				getTextColor();
int				setTextColor(int color);

extern const int g_buf_size;
extern char		g_buf[65536];
//extern char		g_buf[1024];
void			emergencyPrint(int x, int y, const char* format, ...);
void			GUIPrint(int x, int y, const char* format, ...);
void			GUIPrint(int x, int y, int value);
int				print(int x, int y, int tab_origin, const char *a, int length);
int				print(int x, int y, const char *a, int length);

//2D API:
struct			PenBrush
{
	HPEN hPen;
	HBRUSH hBrush;
	int pen_color, brush_color;
	PenBrush():hPen(0), hBrush(0), pen_color(0), brush_color(0){}
	PenBrush(int color);
	~PenBrush();
	void use();
	void drop();
};
struct			Pen
{
	HPEN hPen;
	int color;
	Pen():hPen(0), color(0){}
	Pen(int color);
	~Pen();
	void set(int color);
	void destroy();
	void use();
	void drop();
};
void			line(int x1, int y1, int x2, int y2);
void			moveTo(int x, int y);
void			lineTo(int x, int y);
void			rectangle(int x1, int y1, int x2, int y2);
void			setPixel(int x, int y, int color);

void			vertical_line_mul(int x, int y1, int y2, double Mr, double Mg, double Mb);
void			vertical_line_equation(int x, int y1, int y2, double Ar, double Ag, double Ab, double a);
void			vertical_line_equation(int x, int y1, int y2, double a);
void			vertical_line_antiequation(int x, int y1, int y2, double Ar, double Ag, double Ab, double a, double Br, double Bg, double Bb);
void			vertical_line_antiequation(int x, int y1, int y2, double a);
void			vertical_line_logic_inequality(int x, int y1, int y2, double Ar, double Ag, double Ab);
void			vertical_line_logic_inequality(int x, int y1, int y2);
void			dim_screen();

//3D API:
void			newframe();
void			clear_depth_buffer();
void			resize_3D(int w, int h);

void			gl_line_3D(dvec3 const &p1, dvec3 const &p2, int color);
//void			line_3D(Camera const &cam, Region const &r, dvec2 const &s1, dvec3 const &cp1, dvec2 const &s2, dvec3 const &cp2, int lineColor);
void			line_3D(Camera const &cam, dvec2 const &s1, dvec3 const &cp1, dvec2 const &s2, dvec3 const &cp2, int lineColor);
void			point_3D(Camera const &cam, dvec3 &p_world, int lineColor);
void			point_3D_cam(Camera const &cam, dvec3 &cp, int lineColor);//software line API
void			point_3D_2x2(Camera const &cam, dvec3 &p_world, int Rcolor);//C3D
void			point_3D_2x2(Camera const &cam, dvec3 &p_world, int Rcolor, int Icolor, int Jcolor, int Kcolor);//C3D

extern bool		lsw_transparency_multiply;
void			render_solid_transparent(Camera const &cam, dvec3 const &p1, dvec3 const &p2, dvec3 const &p3, int color);
void			render_solid_transparent(Camera const &cam, vec3 const &p1, vec3 const &p2, vec3 const &p3, int color);//C3D contour
void			render_textured_transparent(Camera const &cam, dvec3 const &p1, dvec3 const &p2, dvec3 const &p3, int *texture, int txw, int txh, vec2 const &tx1, mat2 const &txm);//C3D planes

//basics:
void			initiate();
void			finish();
void			show();
void			resize_2D();

//OpenGL API prerequisites
void			copy_to_clipboard(const char *a, int size);
int				floor_log2(unsigned n);
void			prof_add(const char *label, int divisor);
void			vbo_to_clipboard(vec3 const *vn, int vcount, int const *indices, int icount, int ishr);

//
//OpenGL API
//
inline void		check(int line_number)
{
	int code=glGetError();
	if(code&&!broken)
		((int*)&broken)[1]|=line_number, *(int*)&broken=code;
}
inline void		error(int line_number)
{
	if(!broken)
		((int*)&broken)[1]|=line_number, *(int*)&broken|=glGetError();
}
void			gl_initiate(HDC hDC, int w, int h);
void			gl_finish();
inline void		gl_resize(int w, int h){glViewport(0, 0, w, h);}
inline void		gl_newframe(){glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);}
inline void		gl_disabledepthtest(){glDisable(GL_DEPTH_TEST);}
inline void		gl_enabledepthtest(){glEnable(GL_DEPTH_TEST);}
struct			GPUBuffer
{
	unsigned	VBO, EBO;
	int			vertices_stride, vertices_start, normals_stride, normals_start, n_vertices;
	GPUBuffer():VBO(0), EBO(0), n_vertices(0){}
	void create_VN_I(float *VVVNNN, int n_floats, int *indices, int n_ints);
};
int				gl_print(int x, int y, int tab_origin, const char *format, ...);
namespace		GL2_2D
{
	extern bool	continuous;
	void		curve_begin();
	void		curve_point(float x, float y);
	void		draw_curve();
}
namespace		GL2_L3D
{
	void		begin();
	void		push_surface(vec3 const *vn, int vcount_x2, int *idx, int icount, int color);
	void		end();
	void		draw(Camera const &cam, vec3 const &lightpos);
	void		draw_buffer(Camera const &cam, GPUBuffer const &buffer, vec3 const &modelpos, vec3 const &lightpos);
}
namespace		GL2_3D
{
	void		begin();
	void		begin_transparent();
	void		push_square(float x1, float x2, float y1, float y2, float z, int *tx, int txw, int txh);
	//void		push_triangle(vec3 const &a, vec3 const &b, vec3 const &c, int color, int *tx=0, int txw=0, int txh=0);
	void		push_triangle(vec3 const &a, vec3 const &b, vec3 const &c, int color);
	void		push_line_segment(vec3 const &p1, vec3 const &p2, int color);
	void		push_point(vec3 const &p, int color);
	void		end();
	void		draw(Camera const &cam);
}