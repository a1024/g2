//best viewed with tab size of 4 spaces
//g2_avx.cpp - All AVX functions are implemented here.
//All AVX functions are compiled in a separate project as a static library.
//Copyright (C) 2012-2020  Ayman Wagih Mohsen, unless source link provided.
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

//#pragma warning(push)
//#pragma warning(disable:4752)
#ifdef _DEBUG
#pragma warning(push)
#pragma warning(disable:4996)
#endif
#include"g2_common.h"
#include"g2_avx.h"

#ifndef __AVX__
//	#define	__AVX2__
	#define	__AVX__//defined by /arch:AVX, define in VS2010
#endif

	#define	MAX_VECTOR_SIZE 256
#include	"vectorclass.h"
#include	"vectormath_exp.h"
#include	"vectormath_trig.h"
#include	"vectormath_hyp.h"
#include	"vectormath_common.h"
//#include	<immintrin.h>
#ifdef _DEBUG
#pragma warning(pop)
#endif
	#define	VLEAVE	_mm256_zeroupper()
//	#define	VLEAVE
bool		avx_supported=false;
static __m256d minus_one, sign_mask;
struct Vect4d
{
	__declspec(align(32)) __m256d v;
	Vect4d(){}
//	Vect4d():v(_mm256_setzero_pd()){}//CRASH on non-AVX CPU
#ifdef ALIGNED_INTRINSICS
	Vect4d(VectP const &x):v(_mm256_load_pd(x.r)){}
	Vect4d(double const *x):v(_mm256_load_pd(x)){}
#else
	Vect4d(VectP const &x):v(_mm256_loadu_pd(x.r)){}
	Vect4d(double const *x):v(_mm256_loadu_pd(x)){}
//	Vect4d(double *a):v(_mm256_loadu_pd(a)){}
#endif
	Vect4d(__m256d const &v)
	{
		if(avx_supported)
			this->v=v;
	}
	Vect4d(double x):v(_mm256_set1_pd(x)){}
	Vect4d(double v0, double v1, double v2, double v3):v(_mm256_set_pd(v3, v2, v1, v0)){}
	operator __m256d()const{return v;}
	operator Vec4d()const{return v;}
	void set(double v0, double v1, double v2, double v3){v=_mm256_set_pd(v3, v2, v1, v0);}
	void setzero(){v=_mm256_setzero_pd();}
	double& v0(){return v.m256d_f64[0];}
	double& v1(){return v.m256d_f64[1];}
	double& v2(){return v.m256d_f64[2];}
	double& v3(){return v.m256d_f64[3];}
	double v0()const{return v.m256d_f64[0];}
	double v1()const{return v.m256d_f64[1];}
	double v2()const{return v.m256d_f64[2];}
	double v3()const{return v.m256d_f64[3];}
	double& get(int component){return v.m256d_f64[component];}
	double get(int component)const{return v.m256d_f64[component];}
	Vect4d r_is_false()const{return _mm256_cmp_pd(v, _mm256_setzero_pd(), _CMP_EQ_OQ);}
	Vect4d r_is_true()const{return _mm256_xor_pd(_mm256_cmp_pd(v, _mm256_setzero_pd(), _CMP_EQ_OQ), minus_one);}
	Vect4d& operator=(__m256d v){this->v=v; return *this;}
	Vect4d& operator+=(Vect4d const &other){v=_mm256_add_pd(v, other); return *this;}
	Vect4d& operator-=(Vect4d const &other){v=_mm256_sub_pd(v, other); return *this;}
	Vect4d& operator*=(Vect4d const &other){v=_mm256_mul_pd(v, other); return *this;}
	Vect4d& operator/=(Vect4d const &other){v=_mm256_div_pd(v, other); return *this;}
	Vect4d& operator%=(Vect4d const &q)//x%q = x-floor(x/q)*q	SSE4.1
	{
		__m256d v_q=_mm256_div_pd(v, q);
		v_q=_mm256_round_pd(v_q, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);
		v_q=_mm256_mul_pd(v_q, q);
		v=_mm256_sub_pd(v, v_q);
		return *this;
	}

	Vect4d& operator&=(Vect4d const &other){v=_mm256_and_pd(v, other); return *this;}
	Vect4d& operator|=(Vect4d const &other){v=_mm256_or_pd(v, other); return *this;}
	Vect4d& operator^=(Vect4d const &other){v=_mm256_xor_pd(v, other); return *this;}

	Vect4d& ceil_this(){v=_mm256_ceil_pd(v);return *this;}
	Vect4d ceil()const{return _mm256_ceil_pd(v);}
	Vect4d& floor_this(){v=_mm256_floor_pd(v);return *this;}
	Vect4d floor()const{return _mm256_floor_pd(v);}
	Vect4d& round_this(){v=_mm256_round_pd(v, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);return *this;}
	Vect4d round()const{return _mm256_round_pd(v, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);}
	Vect4d complement()const{return _mm256_xor_pd(v, minus_one);}
	Vect4d& abs_this(){v=_mm256_and_pd(v, sign_mask); return *this;}
	Vect4d abs()const{return _mm256_and_pd(v, sign_mask);}
};
Vect4d m_zero, m_one, m_two, m_ones_mask,
	m_sign_mask, m_sign_mask_complement,
	m_pi, m_pi_2,
	m_half, m_third, m_ln2, m_ln10, m_1_ln10,
	m_inf, m_qnan,
	m_one_percent, m_phi, m_sqrt5,
	m_inv_ln10,
	m_ln_phi, m_inv_sqrt5;
__forceinline Vect4d operator~(Vect4d const &x){return _mm256_xor_pd(x, m_ones_mask);}
__forceinline Vect4d operator*(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_mul_pd(a, b));}
__forceinline Vect4d operator/(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_div_pd(a, b));}
__forceinline Vect4d operator%(Vect4d const &x, Vect4d const &q)
{
	__m256d mask=_mm256_xor_pd(_mm256_cmp_pd(q, _mm256_setzero_pd(), _CMP_EQ_OQ), minus_one);//zero: 0, nonzero: FFFF
	__m256d x_q=_mm256_div_pd(x, q);
	x_q=_mm256_floor_pd(x_q);
	x_q=_mm256_mul_pd(x_q, q);
	x_q=_mm256_sub_pd(x, x_q);
	return _mm256_and_pd(x_q, mask);
}
__forceinline Vect4d operator-(Vect4d const &x){return Vect4d(_mm256_xor_pd(x, m_sign_mask_complement));}
__forceinline Vect4d operator+(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_add_pd(a, b));}
__forceinline Vect4d operator-(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_sub_pd(a, b));}
__forceinline Vect4d operator>(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_cmp_pd(a, b, _CMP_GT_OQ));}
__forceinline Vect4d operator<(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_cmp_pd(a, b, _CMP_LT_OQ));}
__forceinline Vect4d operator>=(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_cmp_pd(a, b, _CMP_GE_OQ));}
__forceinline Vect4d operator<=(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_cmp_pd(a, b, _CMP_LE_OQ));}
__forceinline Vect4d operator==(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_cmp_pd(a, b, _CMP_EQ_OQ));}
__forceinline Vect4d operator!=(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_xor_pd(_mm256_cmp_pd(a, b, _CMP_EQ_OQ), m_ones_mask));}
__forceinline Vect4d operator&(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_and_pd(a, b));}
__forceinline Vect4d operator|(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_or_pd(a, b));}
__forceinline Vect4d operator^(Vect4d const &a, Vect4d const &b){return Vect4d(_mm256_xor_pd(a, b));}
__forceinline Vect4d and(Vect4d const &a, Vect4d const &b){return a&b;}
__forceinline Vect4d sqrt(Vect4d const &x){return Vect4d(_mm256_sqrt_pd(x));}
__forceinline Vect4d is_not_nan(Vect4d const &x){return _mm256_cmp_pd(x, x, _CMP_EQ_OQ);}
__forceinline Vect4d isinf(Vect4d const &x){return _mm256_cmp_pd(x.abs(), m_inf, _CMP_EQ_OQ);}
#pragma warning(push)
#pragma warning(disable:4172)
struct Quat4d;
struct Comp4d
{
	Vect4d r, i;
	Comp4d(){}
	Comp4d(CompP const &x):r(x.r), i(x.i){}
	Comp4d(__m256d const &r, __m256d const &i):r(r), i(i){}
	Comp4d(Vect4d const &r, Vect4d const &i):r(r), i(i){}
	Vect4d r_is_false()const{return r.r_is_false();}
	Vect4d c_is_false()const{return r.r_is_false()&i.r_is_false();}
	Vect4d r_is_true()const{return r!=m_zero;}
	Vect4d c_is_true()const{return ((r==m_zero)&(i==m_zero)).complement();}
	Comp4d floor()const{return Comp4d(r.floor(), i.floor());}
	Comp4d ceil()const{return Comp4d(r.ceil(), i.ceil());}
	Comp4d round()const{return Comp4d(r.round(), i.round());}
	Vect4d abs()const{return sqrt(r*r+i*i);}
	Vect4d arg()const{return Vect4d(::atan2(i, r));}
	Comp4d& operator+=(Comp4d const &b){r+=b.r, i+=b.i; return *this;}
	Comp4d& operator-=(Comp4d const &b){r-=b.r, i-=b.i; return *this;}
	Comp4d& operator*=(Comp4d const &b)
	{
		Vect4d rr=r*b.r-i*b.i, ri=r*b.i+i*b.r;
		r=rr, i=ri;
		return *this;
	}
	Comp4d& operator*=(Vect4d const &b){r*=b, i*=b; return *this;}
	Comp4d& operator/=(Comp4d const &b)
	{
		Vect4d _1_mag_b=m_one/sqrt(b.r*b.r+b.i*b.i);
		Vect4d
			rr=(b.r*r+b.i*i)*_1_mag_b,
			ri=(b.r*i-b.i*r)*_1_mag_b;
		r=rr, i=ri;
		return *this;
	}
	Comp4d& operator/=(Vect4d const &br){r/=br, i/=br; return *this;}
	Quat4d& operator/=(Quat4d const &b);
	Comp4d& operator^=(Comp4d const &b)
	{
		Comp4d t(::log(sqrt(r*r+i*i)), atan2(i, r));
		t*=b;
		Vect4d r0=::exp(t.r);
		Vec4d sin_ti, cos_ti;
		sin_ti=sincos(&cos_ti, t.i);
		r=r0*Vect4d(cos_ti), i=r0*Vect4d(sin_ti);
		return *this;
	}
	Comp4d& operator^=(Vect4d const &br)
	{
		Comp4d t(::log(sqrt(r*r+i*i)), atan2(i, r));
		t*=br;
		Vect4d r0=::exp(t.r);
		Vec4d sin_ti, cos_ti;
		sin_ti=sincos(&cos_ti, t.i);
		r=r0*Vect4d(cos_ti), i=r0*Vect4d(sin_ti);
		return *this;
	}
};
__forceinline Comp4d operator*(Comp4d const &a, Comp4d const &b){return Comp4d(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r);}
__forceinline Comp4d operator*(Comp4d const &a, Vect4d const &br){return Comp4d(a.r*br, a.i*br);}
__forceinline Comp4d operator*(Vect4d const &ar, Comp4d const &b){return Comp4d(ar*b.r, ar*b.i);}
__forceinline Comp4d operator/(Comp4d const &a, Comp4d const &b)
{
	Vect4d _1_mag_b=m_one/(b.r*b.r+b.i*b.i);
	return Comp4d((b.r*a.r+b.i*a.i)*_1_mag_b, (b.r*a.i-b.i*a.r)*_1_mag_b);
}
__forceinline Comp4d operator/(Comp4d const &a, Vect4d const &br){return Comp4d(a.r/br, a.i/br);}
__forceinline Comp4d operator/(Vect4d const &ar, Comp4d const &b)
{
	Vect4d _ar_mag_b=ar/(b.r*b.r+b.i*b.i);
	return Comp4d(b.r*_ar_mag_b, -b.i*_ar_mag_b);
}
__forceinline Comp4d operator+(Comp4d const &a, Comp4d const &b){return Comp4d(a.r+b.r, a.i+b.i);}
__forceinline Comp4d operator+(Comp4d const &a, Vect4d const &b){return Comp4d(a.r+b, a.i);}
__forceinline Comp4d operator+(Vect4d const &a, Comp4d const &b){return Comp4d(a+b.r, b.i);}
__forceinline Comp4d operator-(Comp4d const &a, Comp4d const &b){return Comp4d(a.r-b.r, a.i-b.i);}
__forceinline Comp4d operator-(Comp4d const &a, Vect4d const &b){return Comp4d(a.r-b, a.i);}
__forceinline Comp4d operator-(Vect4d const &a, Comp4d const &b){return Comp4d(a-b.r, -b.i);}
__forceinline Comp4d operator-(Comp4d const &a){return Comp4d(-a.r, -a.i);}
__forceinline Vect4d operator==(Comp4d const &a, Comp4d const &b){return (a.r==b.r)&(a.i==b.i);}
__forceinline Vect4d operator==(Comp4d const &a, Vect4d const &b){return (a.r==b)&(a.i==m_zero);}
__forceinline Vect4d operator==(Vect4d const &a, Comp4d const &b){return (a==b.r)&(m_zero==b.i);}
__forceinline Vect4d operator!=(Comp4d const &a, Comp4d const &b){return (a.r!=b.r)|(a.i!=b.i);}
__forceinline Vect4d operator!=(Comp4d const &a, Vect4d const &b){return (a.r!=b)|(a.i!=m_zero);}
__forceinline Vect4d operator!=(Vect4d const &a, Comp4d const &b){return (a!=b.r)|(m_zero!=b.i);}
__forceinline Comp4d operator|(Comp4d const &a, Comp4d const &b){return Comp4d(a.r|b.r, a.i|b.i);}
__forceinline Comp4d operator|(Comp4d const &a, Vect4d const &b){return Comp4d(a.r|b, a.i);}
__forceinline Comp4d operator|(Vect4d const &a, Comp4d const &b){return Comp4d(a|b.r, b.i);}
__forceinline Comp4d and(Comp4d const &a, Vect4d const &b){return Comp4d(a.r&b, a.i&b);}
__forceinline Comp4d and(Vect4d const &a, Comp4d const &b){return Comp4d(a&b.r, a&b.i);}
__forceinline Comp4d operator%(Comp4d const &a, Comp4d const &b){return and((a-(a/b).floor()*b), a.c_is_true());}
__forceinline Comp4d operator%(Comp4d const &a, Vect4d const &b){return and((a-(a/b).floor()*b), a.c_is_true());}
__forceinline Comp4d operator%(Vect4d const &a, Comp4d const &b){return and((a-(a/b).floor()*b), a.r_is_true());}
__forceinline Comp4d log(Comp4d const &x){return Comp4d(::log(sqrt(x.r*x.r+x.i*x.i)), ::atan2(x.i, x.r));}
__forceinline Comp4d exp(Comp4d const &x)
{
	Vec4d sin_xi, cos_xi;
	sin_xi=sincos(&cos_xi, x.i);
	Vect4d exp_xr=::exp(x.r);
	return Comp4d(exp_xr*Vect4d(cos_xi), exp_xr*Vect4d(sin_xi));
}
__forceinline Comp4d sqrt(Comp4d const &x)
{
	auto s=sqrt(m_two*(x.r+x.abs()));//3 sqrts
	auto i=sqrt(-x.r), mask=s==m_zero;
	return Comp4d(s*m_half, x.i/s&mask.complement()|i&mask);
}
__forceinline Comp4d operator^(Comp4d const &a, Comp4d const &b)
{
	Vect4d mask=a.r==m_zero&a.i==m_zero&b.r==m_zero&b.i==m_zero, mask_c=mask.complement();
	Comp4d t(::log(sqrt(a.r*a.r+a.i*a.i)), atan2(a.i, a.r));
	t*=b;
	Vect4d r0=::exp(t.r);
	Vec4d sin_ti, cos_ti;
	sin_ti=sincos(&cos_ti, t.i);
	return Comp4d(r0*Vect4d(cos_ti)&mask_c|m_one&mask, r0*Vect4d(sin_ti)&mask_c);
}
__forceinline Comp4d operator^(Comp4d const &a, Vect4d const &br)
{
	Vect4d mask=a.r==m_zero&a.i==m_zero&br==m_zero, mask_c=mask.complement();
	Comp4d t(::log(sqrt(a.r*a.r+a.i*a.i)), atan2(a.i, a.r));
	t*=br;
	Vect4d r0=::exp(t.r);
	Vec4d sin_ti, cos_ti;
	sin_ti=sincos(&cos_ti, t.i);
	return Comp4d(r0*Vect4d(cos_ti)&mask_c|m_one&mask, r0*Vect4d(sin_ti)&mask_c);
}
__forceinline Comp4d operator^(Vect4d const &ar, Comp4d const &b)
{
	Vect4d mask=ar==m_zero&b.r==m_zero&b.i==m_zero, mask_c=mask.complement();
	Comp4d t(::log(ar.abs()), atan2(_mm256_setzero_pd(), ar));
	t*=b;
	Vect4d r0=::exp(t.r);
	Vec4d sin_ti, cos_ti;
	sin_ti=sincos(&cos_ti, t.i);
	return Comp4d(r0*Vect4d(cos_ti)&mask_c|m_one&mask, r0*Vect4d(sin_ti)&mask_c);
}
struct Quat4d
{
	Vect4d r, i, j, k;
	Quat4d(){}
	Quat4d(QuatP const &x):r(x.r), i(x.i), j(x.j), k(x.k){}
	Quat4d(__m256d const &r, __m256d const &i, __m256d const &j, __m256d const &k):r(r), i(i), j(j), k(k){}
	Quat4d(Vect4d const &r, Vect4d const &i, Vect4d const &j, Vect4d const &k):r(r), i(i), j(j), k(k){}
	Quat4d(Comp4d const &x):r(x.r), i(x.i){}
	void setzero(){r.setzero(), i.setzero(), j.setzero(), k.setzero();}
	Vect4d r_is_false()const{return r.r_is_false();}
	Vect4d c_is_false()const{return r.r_is_false()&i.r_is_false();}
	Vect4d q_is_false()const{return r.r_is_false()&i.r_is_false()&j.r_is_false()&k.r_is_false();}
	Vect4d r_is_true()const{return r!=m_zero;}
	Vect4d c_is_true()const{return ((r==m_zero)&(i==m_zero)).complement();}
	Vect4d q_is_true()const{return ((r==m_zero)&(i==m_zero)&(j==m_zero)&(k==m_zero)).complement();}
	Quat4d floor()const{return Quat4d(r.floor(), i.floor(), j.floor(), k.floor());}
	Quat4d ceil()const{return Quat4d(r.ceil(), i.ceil(), j.ceil(), k.ceil());}
	Quat4d round()const{return Quat4d(r.round(), i.round(), j.round(), k.round());}
	Vect4d abs()const{return sqrt(r*r+i*i+j*j+k*k);}
	Quat4d& operator+=(Quat4d const &b){r+=b.r, i+=b.i, j+=b.j, k+=b.k;return *this;}
	Quat4d& operator+=(Comp4d const &b){r+=b.r, i+=b.i;return *this;}
	Quat4d& operator+=(Vect4d const &br){r+=br;return *this;}
	Quat4d& operator-=(Quat4d const &b){r-=b.r, i-=b.i, j-=b.j, k-=b.k;return *this;}
	Quat4d& operator-=(Comp4d const &b){r-=b.r, i-=b.i;return *this;}
	Quat4d& operator-=(Vect4d const &br){r-=br;return *this;}
	Quat4d& operator*=(Quat4d const &b)
	{
		Vect4d
			rr=r*b.r+i*b.i+j*b.j+k*b.k,
			ri=r*b.i+i*b.r+j*b.k-k*b.j,
			rj=r*b.j-i*b.k+j*b.r+k*b.i,
			rk=r*b.k+i*b.j+j*b.i+k*b.r;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat4d& operator*=(Comp4d const &b)
	{
		Vect4d
			rr=r*b.r+i*b.i,
			ri=r*b.i+i*b.r,
			rj=j*b.r+k*b.i,
			rk=j*b.i+k*b.r;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat4d& operator*=(Vect4d const &b){r*=b, i*=b, j*=b, k*=b;return *this;}
	Quat4d& operator/=(Quat4d const &b)
	{
		Vect4d _1_mag_y=m_one/sqrt(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
		Vect4d
			rr=(b.r*r+b.i*i+b.j*j+b.k*k)*_1_mag_y,
			ri=(b.r*i-b.i*r-b.j*k+b.k*j)*_1_mag_y,
			rj=(b.r*j+b.i*k-b.j*r-b.k*i)*_1_mag_y,
			rk=(b.r*k-b.i*j+b.j*i-b.k*r)*_1_mag_y;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat4d& operator/=(Comp4d const &b)
	{
		Vect4d _1_mag_y=m_one/sqrt(b.r*b.r+b.i*b.i);
		Vect4d
			rr=(b.r*r+b.i*i)*_1_mag_y,
			ri=(b.r*i-b.i*r)*_1_mag_y,
			rj=(b.r*j+b.i*k)*_1_mag_y,
			rk=(b.r*k-b.i*j)*_1_mag_y;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat4d& operator/=(Vect4d const &b){r/=b, i/=b, j/=b, k/=b;return *this;}
	Quat4d& operator^=(Quat4d const &b)
	{
		Vect4d mag_v=i*i+j*j+k*k;
		Vect4d ln_mag_a=::log(sqrt(r*r+mag_v));
		mag_v=sqrt(mag_v);
		Vect4d v_mul=acos((r/ln_mag_a));
		v_mul/=mag_v;
		Quat4d t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
		t=Quat4d(
			b.r*t.r+b.i*t.i+b.j*t.j+b.k*t.k,
			b.r*t.i+b.i*t.r+b.j*t.k-b.k*t.j,
			b.r*t.j-b.i*t.k+b.j*t.r+b.k*t.i,
			b.r*t.k+b.i*t.j+b.j*t.i+b.k*t.r);
		Vect4d mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
		Vec4d sin_mu, cos_mu;
		sin_mu=sincos(&cos_mu, mag_u);
		Vect4d exp_tr=::exp(t.r);
		v_mul=exp_tr*Vect4d(sin_mu)/mag_u;
		r=exp_tr*Vect4d(cos_mu), i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
		return *this;
	}
	Quat4d& operator^=(Comp4d const &b)
	{
		Vect4d mag_v=i*i+j*j+k*k;
		Vect4d ln_mag_a=::log(sqrt(r*r+mag_v));
		mag_v=sqrt(mag_v);
		Vect4d v_mul=acos((r/ln_mag_a));
		v_mul/=mag_v;
		Quat4d t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
		t=Quat4d(
			b.r*t.r+b.i*t.i,
			b.r*t.i+b.i*t.r,
			b.r*t.j-b.i*t.k,
			b.r*t.k+b.i*t.j);
		Vect4d mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
		Vec4d sin_mu, cos_mu;
		sin_mu=sincos(&cos_mu, mag_u);
		Vect4d exp_tr=::exp(t.r);
		v_mul=exp_tr*Vect4d(sin_mu)/mag_u;
		r=exp_tr*Vect4d(cos_mu), i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
		return *this;
	}
	Quat4d& operator^=(Vect4d const &br)
	{
		Vect4d mag_v=i*i+j*j+k*k;
		Vect4d ln_mag_a=::log(sqrt(r*r+mag_v));
		mag_v=sqrt(mag_v);
		Vect4d v_mul=acos((r/ln_mag_a));
		v_mul/=mag_v;
		Quat4d t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
		t=Quat4d(br*t.r, br*t.i, br*t.j, br*t.k);
		Vect4d mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
		Vec4d sin_mu, cos_mu;
		sin_mu=sincos(&cos_mu, mag_u);
		Vect4d exp_tr=::exp(t.r);
		v_mul=exp_tr*Vect4d(sin_mu)/mag_u;
		r=exp_tr*Vect4d(cos_mu), i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
		return *this;
	}
};
__forceinline Vect4d operator==(Quat4d const &a, Quat4d const &b){return (a.r==b.r)&(a.i==b.i)&(a.j==b.j)&(a.k==b.k);}
__forceinline Vect4d operator==(Quat4d const &a, Comp4d const &b){return (a.r==b.r)&(a.i==b.i)&(a.j==m_zero)&(a.k==m_zero);}
__forceinline Vect4d operator==(Quat4d const &a, Vect4d const &b){return (a.r==b)&(a.i==m_zero)&(a.j==m_zero)&(a.k==m_zero);}
__forceinline Vect4d operator==(Comp4d const &a, Quat4d const &b){return (a.r==b.r)&(a.i==b.i)&(m_zero==b.j)&(m_zero==b.k);}
__forceinline Vect4d operator==(Vect4d const &a, Quat4d const &b){return (a==b.r)&(m_zero==b.i)&(m_zero==b.j)&(m_zero==b.k);}

__forceinline Vect4d operator!=(Quat4d const &a, Quat4d const &b){return (a.r!=b.r)|(a.i!=b.i)|(a.j!=b.j)|(a.k!=b.k);}
__forceinline Vect4d operator!=(Quat4d const &a, Comp4d const &b){return (a.r!=b.r)|(a.i!=b.i)|(a.j!=m_zero)|(a.k!=m_zero);}
__forceinline Vect4d operator!=(Quat4d const &a, Vect4d const &b){return (a.r!=b)|(a.i!=m_zero)|(a.j!=m_zero)|(a.k!=m_zero);}
__forceinline Vect4d operator!=(Comp4d const &a, Quat4d const &b){return (a.r!=b.r)|(a.i!=b.i)|(m_zero!=b.j)|(m_zero!=b.k);}
__forceinline Vect4d operator!=(Vect4d const &a, Quat4d const &b){return (a!=b.r)|(m_zero!=b.i)|(m_zero!=b.j)|(m_zero!=b.k);}

__forceinline Quat4d operator+(Quat4d const &a , Quat4d const &b ){return Quat4d(a.r+b.r, a.i+b.i, a.j+b.j, a.k+b.k);}
__forceinline Quat4d operator+(Quat4d const &a , Comp4d const &b ){return Quat4d(a.r+b.r, a.i+b.i, a.j, a.k);}
__forceinline Quat4d operator+(Quat4d const &a , Vect4d const &br){return Quat4d(a.r+br, a.i, a.j, a.k);}
__forceinline Quat4d operator+(Comp4d const &a , Quat4d const &b ){return Quat4d(a.r+b.r, a.i+b.i, b.j, b.k);}
__forceinline Quat4d operator+(Vect4d const &ar, Quat4d const &b ){return Quat4d(ar+b.r, b.i, b.j, b.k);}

__forceinline Quat4d operator-(Quat4d const &a , Quat4d const &b ){return Quat4d(a.r-b.r, a.i-b.i, a.j-b.j, a.k-b.k);}
__forceinline Quat4d operator-(Quat4d const &a , Comp4d const &b ){return Quat4d(a.r-b.r, a.i-b.i, a.j, a.k);}
__forceinline Quat4d operator-(Quat4d const &a , Vect4d const &br){return Quat4d(a.r-br, a.i, a.j, a.k);}
__forceinline Quat4d operator-(Comp4d const &a , Quat4d const &b ){return Quat4d(a.r-b.r, a.i-b.i, -b.j, -b.k);}
__forceinline Quat4d operator-(Vect4d const &ar, Quat4d const &b ){return Quat4d(ar-b.r, -b.i, -b.j, -b.k);}
__forceinline Quat4d operator-(Quat4d const &a){return Quat4d(-a.r, -a.i, -a.j, -a.k);}
__forceinline Quat4d operator*(Quat4d const &a, Quat4d const &b)
{
	return Quat4d(
		a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,
		a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,
		a.r*b.j-a.i*b.k+a.j*b.r+a.k*b.i,
		a.r*b.k+a.i*b.j-a.j*b.i+a.k*b.r);
}
__forceinline Quat4d operator*(Quat4d const &a, Comp4d const &b)
{
	return Quat4d(
		a.r*b.r-a.i*b.i,
		a.r*b.i+a.i*b.r,
		a.j*b.r+a.k*b.i,
		-a.j*b.i+a.k*b.r);
}
__forceinline Quat4d operator*(Quat4d const &a, Vect4d const &br){return Quat4d(a.r*br, a.i*br, a.j*br, a.k*br);}
__forceinline Quat4d operator*(Comp4d const &a, Quat4d const &b)
{
	return Quat4d(
		a.r*b.r-a.i*b.i,
		a.r*b.i+a.i*b.r,
		a.r*b.j-a.i*b.k,
		a.r*b.k+a.i*b.j);
}
__forceinline Quat4d operator*(Vect4d const &ar, Quat4d const &b){return Quat4d(ar*b.r, ar*b.i, ar*b.j, ar*b.k);}
__forceinline Quat4d operator/(Quat4d const &a, Quat4d const &b)
{
	Vect4d _1_mag_y=m_one/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat4d(
		(b.r*a.r+b.i*a.i+b.j*a.j+b.k*a.k)*_1_mag_y,
		(b.r*a.i-b.i*a.r-b.j*a.k+b.k*a.j)*_1_mag_y,
		(b.r*a.j+b.i*a.k-b.j*a.r-b.k*a.i)*_1_mag_y,
		(b.r*a.k-b.i*a.j+b.j*a.i-b.k*a.r)*_1_mag_y);
}
__forceinline Quat4d& Comp4d::operator/=(Quat4d const &b)
{
	Vect4d _1_mag_y=m_one/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat4d(
		(b.r*r+b.i*i)*_1_mag_y,
		(b.r*i-b.i*r)*_1_mag_y,
		(-b.j*r-b.k*i)*_1_mag_y,
		(b.j*i-b.k*r)*_1_mag_y);
}
__forceinline Quat4d operator/(Quat4d const &a, Comp4d const &b)
{
	Vect4d _1_mag_y=m_one/(b.r*b.r+b.i*b.i);
	return Quat4d(
		(b.r*a.r+b.i*a.i)*_1_mag_y,
		(b.r*a.i-b.i*a.r)*_1_mag_y,
		(b.r*a.j+b.i*a.k)*_1_mag_y,
		(b.r*a.k-b.i*a.j)*_1_mag_y);
}
__forceinline Quat4d operator/(Quat4d const &a, Vect4d const &br)
{
	Vect4d _1_mag_y=m_one/br;
	return Quat4d(a.r*_1_mag_y, a.i*_1_mag_y, a.j*_1_mag_y, a.k*_1_mag_y);
}
__forceinline Quat4d operator/(Comp4d const &a, Quat4d const &b)
{
	Vect4d _1_mag_y=m_one/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat4d(
		(b.r*a.r+b.i*a.i)*_1_mag_y,
		(b.r*a.i-b.i*a.r)*_1_mag_y,
		(-b.j*a.r-b.k*a.i)*_1_mag_y,
		(b.j*a.i-b.k*a.r)*_1_mag_y);
}
__forceinline Quat4d operator/(Vect4d const &ar, Quat4d const &b)
{
	Vect4d _ar_mag_y=ar/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat4d(b.r*_ar_mag_y, -b.i*_ar_mag_y, -b.j*_ar_mag_y, -b.k*_ar_mag_y);
}
__forceinline Quat4d log(Quat4d const &x)
{
	Vect4d mag_v=x.i*x.i+x.j*x.j+x.k*x.k;

	Vect4d real=mag_v==m_zero, real_c=real.complement(), rr=log(x.r.abs().v);

	Vect4d mag_x=sqrt(x.r*x.r+mag_v);
	mag_v=sqrt(mag_v);
	Vect4d u_mul=Vect4d(::acos((x.r/mag_x).v))/mag_v&real_c;
	Vect4d ln_mag_x=::log(mag_x.v);
	return Quat4d(rr&real|ln_mag_x&real_c, m_pi&real|x.i*u_mul, x.j*u_mul, x.k*u_mul);
}
__forceinline Quat4d exp(Quat4d const &x)
{
	Vect4d exp_r=::exp(x.r);
	Vect4d mag_v=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
	Vec4d sin_v, cos_v;
	sin_v=sincos(&cos_v, mag_v);
	Vect4d v_mul=exp_r*Vect4d(sin_v)/mag_v;
	return Quat4d(exp_r*Vect4d(cos_v), x.i*v_mul, x.j*v_mul, x.k*v_mul);
}
__forceinline Quat4d sqrt(Quat4d const &x)
{
	auto s=sqrt(m_two*(x.r+x.abs()));
	auto i=sqrt(-x.r), mask=s==m_zero, mask_c=mask.complement();
	auto _1_s=m_one/s&mask_c;
	return Quat4d(s*m_half, x.i*_1_s|i&mask, x.j*_1_s, x.k*_1_s);
//	return exp(m_half*log(x));
}
__forceinline Quat4d operator^(Quat4d const &a, Quat4d const &b)
{
	Vect4d mask=a.r==m_zero&a.i==m_zero&a.j==m_zero&a.k==m_zero&b.r==m_zero&b.i==m_zero&b.j==m_zero&b.k==m_zero, mask_c=mask.complement();
	Quat4d q=exp(b*log(a));
	return Quat4d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat4d operator^(Quat4d const &a, Comp4d const &b)
{
	Vect4d mask=a.r==m_zero&a.i==m_zero&a.j==m_zero&a.k==m_zero&b.r==m_zero&b.i==m_zero, mask_c=mask.complement();
	Quat4d q=exp(b*log(a));
	return Quat4d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat4d operator^(Quat4d const &a, Vect4d const &b)
{
	Vect4d mask=a.r==m_zero&a.i==m_zero&a.j==m_zero&a.k==m_zero&b==m_zero, mask_c=mask.complement();
	Quat4d q=exp(b*log(a));
	return Quat4d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat4d operator^(Comp4d const &a, Quat4d const &b)
{
	Vect4d mask=a.r==m_zero&a.i==m_zero&b.r==m_zero&b.i==m_zero&b.j==m_zero&b.k==m_zero, mask_c=mask.complement();
	Quat4d q=exp(b*log(a));
	return Quat4d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat4d operator^(Vect4d const &a, Quat4d const &b)
{
	Vect4d mask=a==m_zero&b.r==m_zero&b.i==m_zero&b.j==m_zero&b.k==m_zero, mask_c=mask.complement();
	Quat4d q=exp(b*Vect4d(::log(a)));
	return Quat4d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat4d operator|(Quat4d const &a, Quat4d const &b){return Quat4d(a.r|b.r, a.i|b.i, a.j|b.j, a.k|b.k);}
__forceinline Quat4d operator|(Quat4d const &a, Comp4d const &b){return Quat4d(a.r|b.r, a.i|b.i, a.j, a.k);}
__forceinline Quat4d operator|(Quat4d const &a, Vect4d const &b){return Quat4d(a.r|b, a.i, a.j, a.k);}
__forceinline Quat4d operator|(Comp4d const &a, Quat4d const &b){return Quat4d(a.r|b.r, a.i|b.i, b.j, b.k);}
__forceinline Quat4d operator|(Vect4d const &a, Quat4d const &b){return Quat4d(a|b.r, b.i, b.j, b.k);}
__forceinline Quat4d and(Quat4d const &a, Vect4d const &b){return Quat4d(a.r&b, a.i&b, a.j&b, a.k&b);}
__forceinline Quat4d and(Vect4d const &a, Quat4d const &b){return Quat4d(a&b.r, a&b.i, a&b.j, a&b.k);}
__forceinline Quat4d operator%(Quat4d const &a, Quat4d const &b){return and((a-(a/b).floor()*b), a.q_is_true());}
__forceinline Quat4d operator%(Quat4d const &a, Comp4d const &b){return and((a-(a/b).floor()*b), a.q_is_true());}
__forceinline Quat4d operator%(Quat4d const &a, Vect4d const &b){return and((a-(a/b).floor()*b), a.q_is_true());}
__forceinline Quat4d operator%(Comp4d const &a, Quat4d const &b){return and((a-(a/b).floor()*b), a.c_is_true());}
__forceinline Quat4d operator%(Vect4d const &a, Quat4d const &b){return and((a-(a/b).floor()*b), a.r_is_true());}
Comp4d m_i;
void avx_initialize()
{
	minus_one=_mm256_castsi256_pd(_mm256_set1_epi32(-1)), sign_mask=_mm256_castsi256_pd(_mm256_set_epi32(0x7FFFFFFF, -1, 0x7FFFFFFF, -1, 0x7FFFFFFF, -1, 0x7FFFFFFF, -1));

	m_zero=0., m_one=1, m_two=2, m_ones_mask=_mm256_castsi256_pd(_mm256_set1_epi32(-1)),
	m_sign_mask=sign_mask, m_sign_mask_complement=m_sign_mask.complement(),
	m_pi=G2::_pi, m_pi_2=G2::_pi_2,
	m_half=0.5, m_third=1./3, m_ln2=G2::_ln2, m_ln10=G2::_ln10, m_1_ln10=1/G2::_ln10,
	m_inf=_HUGE, m_qnan=G2::_qnan,
	m_one_percent=0.01, m_phi=G2::_phi, m_sqrt5=G2::_sqrt5;

	m_inv_ln10=m_one/m_ln10;
	m_ln_phi=::log(G2::_phi), m_inv_sqrt5=1/G2::_sqrt5;

	m_i.r=0., m_i.i=1;
}
#pragma warning(pop)
void minmax_avx(double *a, int size, double *lo_hi)//size multiple of 4, 32 byte aligned
{
	//if((int)a&31)
	//	return;

	__m256d min=_mm256_load_pd(a), max=min;
	for(int k=4;k+3<size;k+=4)
	{
		__m256d vk=_mm256_load_pd(a+k);		//1.576 c/v on R7 2700
		min=_mm256_min_pd(min, vk);
		max=_mm256_max_pd(max, vk);
	}
	__m256d
		t0=_mm256_shuffle_pd(min, max, 0),	//{max2, min2, max0, min0}	//{qword3, qword2, qword1, qword0}
		t1=_mm256_shuffle_pd(min, max, 0xF),//{max3, min3, max1, min1}
		t2=_mm256_min_pd(t0, t1),	//{~,		minb2,	~,		minb0}
		t3=_mm256_max_pd(t0, t1),	//{maxb3,	~,		maxb1,	~}
		t5=_mm256_shuffle_pd(t2, t3, 0xA);//{maxb3, minb2, maxb1, minb0}
	__m128d u0=_mm256_castpd256_pd128(t5);//{maxb1, minb0}
	t5=_mm256_permute2f128_pd(t5, t5, 1);
	__m128d u1=_mm256_castpd256_pd128(t5);//{maxb3, minb2}
	__m128d u2=_mm_min_pd(u0, u1), u3=_mm_max_pd(u0, u1),//{~, minc0}	{maxc1, ~}
		u_minmax=_mm_shuffle_pd(u2, u3, 2);//{maxc1, minc0}
	_mm_store_pd(lo_hi, u_minmax);
	//lo=min.m256d_f64[0];
	//if(lo>min.m256d_f64[1])
	//	lo=min.m256d_f64[1];
	//if(lo>min.m256d_f64[2])
	//	lo=min.m256d_f64[2];
	//if(lo>min.m256d_f64[3])
	//	lo=min.m256d_f64[3];
	//hi=max.m256d_f64[0];
	//if(hi<max.m256d_f64[1])
	//	hi=max.m256d_f64[1];
	//if(hi<max.m256d_f64[2])
	//	hi=max.m256d_f64[2];
	//if(hi<max.m256d_f64[3])
	//	hi=max.m256d_f64[3];
	VLEAVE;
}
namespace	G2
{
	namespace avx
	{
		void zeroall(){_mm256_zeroall();}
		void zeroupper(){_mm256_zeroupper();}
		__forceinline void assign(VectP &p, Vect4d const &v){_mm256_storeu_pd(p.r, v);}
		//__forceinline void assign(double *p, Vect4d const &v){_mm256_storeu_pd(p, v.v);}
		__forceinline void assign(CompP &p, Comp4d const &v){_mm256_storeu_pd(p.r, v.r), _mm256_storeu_pd(p.i, v.i);}
		__forceinline void assign(QuatP &p, Quat4d const &v){_mm256_storeu_pd(p.r, v.r), _mm256_storeu_pd(p.i, v.i), _mm256_storeu_pd(p.j, v.j), _mm256_storeu_pd(p.k, v.k);}

		void r_r_setzero				(VectP &r, VectP const&)					{assign(r, _mm256_setzero_pd());VLEAVE;}
		void c_c_setzero				(CompP &r, CompP const&)					{assign(r, Comp4d(_mm256_setzero_pd(), _mm256_setzero_pd()));VLEAVE;}
		void q_q_setzero				(QuatP &r, QuatP const&)					{assign(r, Quat4d(_mm256_setzero_pd(), _mm256_setzero_pd(), _mm256_setzero_pd(), _mm256_setzero_pd()));VLEAVE;}

		void  r_r_ceil					(VectP &r, VectP const &x)					{assign(r, Vect4d(x).ceil());VLEAVE;}
		void  c_c_ceil					(CompP &r, CompP const &x)					{assign(r, Comp4d(x).ceil());VLEAVE;}
		void  q_q_ceil					(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x).ceil());VLEAVE;}

		void  r_r_floor					(VectP &r, VectP const &x)					{assign(r, Vect4d(x).floor());VLEAVE;}
		void  c_c_floor					(CompP &r, CompP const &x)					{assign(r, Comp4d(x).floor());VLEAVE;}
		void  q_q_floor					(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x).floor());VLEAVE;}

		void  r_r_round					(VectP &r, VectP const &x)					{assign(r, Vect4d(x).round());VLEAVE;}
		void  c_c_round					(CompP &r, CompP const &x)					{assign(r, Comp4d(x).round());VLEAVE;}
		void  q_q_round					(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x).round());VLEAVE;}

		void  r_r_abs					(VectP &r, VectP const &x)					{assign(r, Vect4d(x).abs());VLEAVE;}
		void  r_c_abs					(VectP &r, CompP const &x)					{assign(r, Comp4d(x).abs());VLEAVE;}
		void  r_q_abs					(VectP &r, QuatP const &x)					{assign(r, Quat4d(x).abs());VLEAVE;}

		void  r_r_arg					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, m_pi&(rx<m_zero)|m_qnan&(rx==m_zero));VLEAVE;}
		void  r_c_arg					(VectP &r, CompP const &x)
		{
			Comp4d cx=x;
			Vect4d mask=(cx.r==m_zero)&(cx.i==m_zero);//zero: FFFF, nonzero: 0
			assign(r, Vect4d(::atan2(cx.i.v, cx.r.v))&mask.complement()|m_qnan&mask);
			VLEAVE;
		}
		void  r_q_arg					(VectP &r, QuatP const &x)					{Quat4d qx=x; assign(r, Vect4d(::acos((qx.r/qx.abs()).v)));VLEAVE;}

		void  r_c_real					(VectP &r, CompP const &x)					{assign(r, Vect4d(x.r));VLEAVE;}

		void  r_c_imag					(VectP &r, CompP const &x)					{assign(r, Vect4d(x.i));VLEAVE;}

		//r_conjugate: assign
		void c_c_conjugate				(CompP &r, CompP const &x)					{assign(r, Comp4d(Vect4d(x.r), -Vect4d(x.i)));VLEAVE;}
		void q_q_conjugate				(QuatP &r, QuatP const &x)					{assign(r, Quat4d(Vect4d(x.r), -Vect4d(x.i), -Vect4d(x.j), -Vect4d(x.k)));VLEAVE;}

		void  c_r_polar					(CompP &r, VectP const &x)					{Vect4d rx=x; assign(r, Comp4d(rx&m_sign_mask, m_pi&(rx<m_zero)|m_qnan&(rx==m_zero)));VLEAVE;}
		void  c_c_polar					(CompP &r, CompP const &x)
		{
			Comp4d cx=x;
			Vect4d mag=cx.abs(), z_mask=mag==m_zero, arg=::atan2(cx.i.v, cx.r.v);
			assign(r, Comp4d(mag, arg&z_mask|m_qnan&z_mask.complement()));
			VLEAVE;
		}
		void  c_q_polar					(CompP &r, QuatP const &x)
		{
			Quat4d qx=x;
			Vect4d mag=qx.abs();
			assign(r, Comp4d(mag, Vect4d(::acos((qx.r/mag).v))));
			VLEAVE;
		}

		//r_cartesian	assign
		void  c_c_cartesian				(CompP &r, CompP const &x)
		{
			Comp4d cx=x;
			Vec4d sin_i, cos_i;
			sin_i=sincos(&cos_i, cx.i.v);
			assign(r, Comp4d(cx.r*Vect4d(cos_i), cx.r*Vect4d(sin_i)));
			VLEAVE;
		}
		void  q_q_cartesian				(QuatP &r, QuatP const &x)
		{
			Quat4d qx=x;
			Vec4d sin_i, cos_i, sin_j, cos_j, sin_k, cos_k;
			sin_i=sincos(&cos_i, qx.i.v);
			sin_j=sincos(&cos_j, qx.j.v);
			sin_k=sincos(&cos_k, qx.k.v);
			cos_k*=qx.r.v;
			assign(r, Quat4d(
				Vect4d(cos_i)*Vect4d(cos_j)*Vect4d(cos_k),
				Vect4d(sin_i)*Vect4d(cos_j)*Vect4d(cos_k),
				Vect4d(sin_j)*Vect4d(cos_k),
				qx.r*Vect4d(sin_k)));
			VLEAVE;
		}

		void r_rr_plus					(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)+Vect4d(y));VLEAVE;}
		void c_rc_plus					(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)+Comp4d(y));VLEAVE;}
		void q_rq_plus					(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)+Quat4d(y));VLEAVE;}
		void c_cr_plus					(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)+Vect4d(y));VLEAVE;}
		void c_cc_plus					(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)+Comp4d(y));VLEAVE;}
		void q_cq_plus					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)+Quat4d(y));VLEAVE;}
		void q_qr_plus					(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)+Vect4d(y));VLEAVE;}
		void q_qc_plus					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)+Comp4d(y));VLEAVE;}
		void q_qq_plus					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)+Quat4d(y));VLEAVE;}

		void  r_r_minus					(VectP &r, VectP const &x)					{assign(r, -Vect4d(x));VLEAVE;}
		void  c_c_minus					(CompP &r, CompP const &x)					{assign(r, -Comp4d(x));VLEAVE;}
		void  q_q_minus					(QuatP &r, QuatP const &x)					{assign(r, -Quat4d(x));VLEAVE;}
		void r_rr_minus					(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)-Vect4d(y));VLEAVE;}
		void c_rc_minus					(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)-Comp4d(y));VLEAVE;}
		void q_rq_minus					(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)-Quat4d(y));VLEAVE;}
		void c_cr_minus					(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)-Vect4d(y));VLEAVE;}
		void c_cc_minus					(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)-Comp4d(y));VLEAVE;}
		void q_cq_minus					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)-Quat4d(y));VLEAVE;}
		void q_qr_minus					(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)-Vect4d(y));VLEAVE;}
		void q_qc_minus					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)-Comp4d(y));VLEAVE;}
		void q_qq_minus					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)-Quat4d(y));VLEAVE;}

		void r_rr_multiply				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)*Vect4d(y));VLEAVE;}
		void c_rc_multiply				(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)*Comp4d(y));VLEAVE;}
		void q_rq_multiply				(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)*Quat4d(y));VLEAVE;}
		void c_cr_multiply				(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)*Vect4d(y));VLEAVE;}
		void c_cc_multiply				(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)*Comp4d(y));VLEAVE;}//(xr+i*xi)(yr+i*yi) = xr*yr-xi*yi+i(xr*yi+xi*yr)
		void q_cq_multiply				(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)*Quat4d(y));VLEAVE;}
		void q_qr_multiply				(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)*Vect4d(y));VLEAVE;}
		void q_qc_multiply				(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)*Comp4d(y));VLEAVE;}
		void q_qq_multiply				(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)*Quat4d(y));VLEAVE;}

		void  r_r_divide				(VectP &r, VectP const &y)					{assign(r, m_one/Vect4d(y));VLEAVE;}
		void  c_c_divide				(CompP &r, CompP const &y)					{assign(r, m_one/Comp4d(y));VLEAVE;}
		void  q_q_divide				(QuatP &r, QuatP const &y)					{assign(r, m_one/Quat4d(y));VLEAVE;}
		void r_rr_divide				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)/Vect4d(y));VLEAVE;}
		void c_rc_divide				(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)/Comp4d(y));VLEAVE;}
		void q_rq_divide				(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)/Quat4d(y));VLEAVE;}
		void c_cr_divide				(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)/Vect4d(y));VLEAVE;}
		void c_cc_divide				(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)/Comp4d(y));VLEAVE;}
		void q_cq_divide				(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)/Quat4d(y));VLEAVE;}
		void q_qr_divide				(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)/Vect4d(y));VLEAVE;}
		void q_qc_divide				(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)/Comp4d(y));VLEAVE;}
		void q_qq_divide				(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)/Quat4d(y));VLEAVE;}

		void r_rr_logic_divides			(VectP &r, VectP const &y, VectP const &x)	{auto t=Vect4d(x)/Vect4d(y); assign(r, (t==t.floor())&m_one);VLEAVE;}//rc_divides, rq_divides: applied to each component
		void r_rc_logic_divides			(VectP &r, VectP const &y, CompP const &x)	{auto t=Comp4d(x)/Vect4d(y); assign(r, (t==t.floor())&m_one);VLEAVE;}
		void r_rq_logic_divides			(VectP &r, VectP const &y, QuatP const &x)	{auto t=Quat4d(x)/Vect4d(y); assign(r, (t==t.floor())&m_one);VLEAVE;}
		void r_cr_logic_divides			(VectP &r, CompP const &y, VectP const &x)	{auto t=Vect4d(x)/Comp4d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&m_one);VLEAVE;}
		void r_cc_logic_divides			(VectP &r, CompP const &y, CompP const &x)	{auto t=Comp4d(x)/Comp4d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&m_one);VLEAVE;}
		void r_cq_logic_divides			(VectP &r, CompP const &y, QuatP const &x)	{auto t=Quat4d(x)/Comp4d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);VLEAVE;}
		void r_qr_logic_divides			(VectP &r, QuatP const &y, VectP const &x)	{auto t=Vect4d(x)/Quat4d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);VLEAVE;}
		void r_qc_logic_divides			(VectP &r, QuatP const &y, CompP const &x)	{auto t=Comp4d(x)/Quat4d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);VLEAVE;}
		void r_qq_logic_divides			(VectP &r, QuatP const &y, QuatP const &x)	{auto t=Quat4d(x)/Quat4d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);VLEAVE;}

		void c_cr_pow					(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)^Vect4d(y));VLEAVE;}
		void c_cc_pow					(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)^Comp4d(y));VLEAVE;}
		void q_cq_pow					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)^Quat4d(y));VLEAVE;}
		void q_qr_pow					(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)^Vect4d(y));VLEAVE;}
		void q_qc_pow					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)^Comp4d(y));VLEAVE;}
		void q_qq_pow					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)^Quat4d(y));VLEAVE;}

		void  c_c_ln					(CompP &r, CompP const &x)					{assign(r, log(Comp4d(x)));VLEAVE;}
		void  q_q_ln					(QuatP &r, QuatP const &x)					{assign(r, log(Quat4d(x)));VLEAVE;}
	
		void  c_c_log					(CompP &r, CompP const &x)					{assign(r, log(Comp4d(x))*m_inv_ln10);VLEAVE;}
		void  q_q_log					(QuatP &r, QuatP const &x)					{assign(r, log(Quat4d(x))*m_inv_ln10);VLEAVE;}
		void c_cr_log					(CompP &r, CompP const &x, VectP const &y)	{assign(r, log(Comp4d(x))/log(Comp4d(Vect4d(y), m_zero)));VLEAVE;}
		void c_cc_log					(CompP &r, CompP const &x, CompP const &y)	{assign(r, log(Comp4d(x))/log(Comp4d(y)));VLEAVE;}
		void q_cq_log					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, log(Comp4d(x))/log(Quat4d(y)));VLEAVE;}
		void q_qc_log					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, log(Quat4d(x))/log(Comp4d(y)));VLEAVE;}
		void q_qq_log					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, log(Quat4d(x))/log(Quat4d(y)));VLEAVE;}
	
	/*	typedef std::complex<double> Complex1d;
		typedef boost::math::quaternion<double> Quaternion1d;
		Complex1d tetrate(double x, double y)
		{
			double log_x=std::log((double)x), ry=y;
			int h=int(::floor(ry-.5)+1);//rounded real part
			Complex1d t(ry-(double)h);//real-round(real)+i*imag
			{
				auto q=std::sqrt(log_x);
				t=(t+1.)*(
					(
						log_x>_ln2?
							(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+_ln2)/log_x-1.)
							+t*	(
									(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+_ln2)/log_x-1.)
									+t*	(
											(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+_ln2)/log_x-1.)
											+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+_ln2)/log_x-1.)
										)
								)
						:
							(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+_ln2)/log_x)
							+t*	(
									(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+_ln2)/log_x)
									+t*	(
											(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+_ln2)/log_x)
											+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+_ln2)/log_x)
										)
								)
					)*t+1.)+std::log(t+2.)/log_x-_ln2/log_x*(1.+t);
			}
			for(;h>0;--h)
				t=std::exp(log_x*t);
			for(;h<0;++h)
				t=std::log(t)/log_x;
			return t;
		}
		Complex1d tetrate(Complex1d &x, double y)
		{
			Complex1d log_x=std::log((Complex1d)x);
			double ry=y;
			int h=int(::floor(ry-.5)+1);//rounded real part
			Complex1d t(ry-(double)h);//real-round(real)+i*imag
			{
				auto q=std::sqrt(log_x);
				t=(t+1.)*(
					(
						log_x.real()>_ln2?
							(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+_ln2)/log_x-1.)
							+t*	(
									(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+_ln2)/log_x-1.)
									+t*	(
											(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+_ln2)/log_x-1.)
											+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+_ln2)/log_x-1.)
										)
								)
						:
							(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+_ln2)/log_x)
							+t*	(
									(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+_ln2)/log_x)
									+t*	(
											(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+_ln2)/log_x)
											+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+_ln2)/log_x)
										)
								)
					)*t+1.)+std::log(t+2.)/log_x-_ln2/log_x*(1.+t);
			}
			for(;h>0;--h)
				t=std::exp(log_x*t);
			for(;h<0;++h)
				t=std::log(t)/log_x;
			return t;
		}
		Quaternion1d tetrate(Quaternion1d &x, double y)
		{
			Quaternion1d qx=x;
			double ry=y;
			if(ry<-1)
			{
				int steps=int(abs(ry));
				Quaternion1d t(ry-::floor(ry)), lrx=log(qx);
				for(int k=0;k<steps;++k)
					t=log(t)/lrx;
				return t;
			}
			else if(ry<=0)
				return Quaternion1d(1+ry);
			else
			{
				int h=int(ry)+1;
				Quaternion1d t(ry-::floor(ry));
				for(int k=0;k<h;++k)
					t=pow(qx, t);
				return t;
			}
		}
		Complex1d tetrate(double x, Complex1d const &y)
		{
			double log_x=std::log((double)x);
			std::complex<double> cy=y;
			int h=int(::floor(cy.real()-.5)+1);//rounded real part
			std::complex<double> t(cy-(double)h);//real-round(real)+i*imag
			{
				auto q=std::sqrt(log_x);
				t=(t+1.)*(
					(
						log_x>_ln2?
							(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+_ln2)/log_x-1.)
							+t*	(
									(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+_ln2)/log_x-1.)
									+t*	(
											(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+_ln2)/log_x-1.)
											+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+_ln2)/log_x-1.)
										)
								)
						:
							(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+_ln2)/log_x)
							+t*	(
									(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+_ln2)/log_x)
									+t*	(
											(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+_ln2)/log_x)
											+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+_ln2)/log_x)
										)
								)
					)*t+1.)+std::log(t+2.)/log_x-_ln2/log_x*(1.+t);
			}
			for(;h>0;--h)
				t=std::exp(log_x*t);
			for(;h<0;++h)
				t=std::log(t)/log_x;
			return t;
		}
		Complex1d tetrate(Complex1d const &x, Complex1d const &y)
		{
			std::complex<double> log_x=std::log((std::complex<double>)x), cy=y;
		//	if(log_x.real()<.03)//abs(log_x)<1.03045453395352
		//	{
		//		if(cy.real()<-1.)
		//			return -30.;
		//		return 1.;
		//	}
			int h=int(::floor(cy.real()-.5)+1);
			std::complex<double> t(cy-(double)h);//real-round(real)+i*imag
			{
			//	bool unassigned=true;
			//	if(log_x.real()<.001)//abs(log_x)<1.00100050016671
			//	{
			//		if(t.real()>-1)//real-round(real)>-1
			//			unassigned=false, t=1.;
			//		else if(t.real()<-1)
			//			unassigned=false, t=-990.;
			//	}
			//	if(unassigned)
			//	{
					std::complex<double> q=std::sqrt(log_x);
					t=(t+1.)*(
						(
							log_x.real()>_ln2?
								(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+_ln2)/log_x-1.)
								+t*	(
										(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+_ln2)/log_x-1.)
										+t*	(
												(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+_ln2)/log_x-1.)
												+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+_ln2)/log_x-1.)
											)
									)
							:
								(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+_ln2)/log_x)
								+t*	(
										(1.1-2.608785958462561*(1.-0.6663562294911147*sqrt(log_x))*sqrt(log_x)-(-0.625+_ln2)/log_x)
										+t*	(
												(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+_ln2)/log_x)
												+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+_ln2)/log_x)
											)
									)
						)*t+1.)+log(t+2.)/log_x-_ln2/log_x*(1.+t);
			//	}
			}
			for(;h>0;--h)
				t=exp(log_x*t);
			for(;h<0;++h)
				t=std::log(t)/log_x;
			return t;
		}
		Comp4d c_rr_tetrate					(Vect4d const &x, Vect4d const &y)
		{
			Complex1d
				lo=tetrate(x.lo(), y.lo()),
				hi=tetrate(x.hi(), y.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_rc_tetrate					(Vect4d const &xr, Comp4d const &y)
		{
			Complex1d
				lo=tetrate(xr.lo(), Complex1d(y.r.lo(), y.i.lo())),
				hi=tetrate(xr.hi(), Complex1d(y.r.hi(), y.i.hi()));
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_cr_tetrate					(Comp4d const &x, Vect4d const &y)
		{
			Complex1d
				lo=tetrate(Complex1d(x.r.lo(), x.i.lo()), y.lo()),
				hi=tetrate(Complex1d(x.r.hi(), x.i.hi()), y.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_cc_tetrate					(Comp4d const &x, Comp4d const &y)
		{
			Complex1d
				lo=tetrate(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.r.lo(), y.i.lo())),
				hi=tetrate(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.r.hi(), y.i.hi()));
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Quat4d q_qr_tetrate					(Quat4d const &x, Vect4d const &y)
		{
			Quaternion1d
				lo=tetrate(Quaternion1d(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo()), y.lo()),
				hi=tetrate(Quaternion1d(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi()), y.hi());
			return Quat4d(Vect4d(lo.R_component_1(), hi.R_component_1()), Vect4d(lo.R_component_2(), hi.R_component_2()),
				Vect4d(lo.R_component_3(), hi.R_component_3()), Vect4d(lo.R_component_4(), hi.R_component_4()));
		}
	
		Complex1d pentate(double x, double y)
		{
			long long h=convert_d2ll_zero(y);
		//	long long h=y.r!=y.r||y.r<-ll_max||y.r>ll_max?0:long long(y);
		//	long long h=std::isnan(y.r)||std::isinf(y.r)?0:long long(y);
			if(h<-2)	return _HUGE;//1/::sin(0);
			if(h==-2)	return -1;
			if(h==-1)	return 0;
			if(h==0)	return 1;
			if(h==1)	return x;
			double rx=x;
			Complex1d result(rx);
			for(int k=0;k<h;++k)
				result=Tetrate::fit1(rx, result);
			return result;
		}
		Complex1d pentate(Complex1d const &x, double y)
		{
			long long h=convert_d2ll_zero(y);
			if(h<-2)	return _HUGE;//1/::sin(0);
			if(h==-2)	return -1;
			if(h==-1)	return 0;
			if(h==0)	return 1;
			if(h==1)	return x;
			Complex1d cx=x, result(cx);
			for(int k=0;k<h;++k)
				result=Tetrate::fit1(cx, result);
			return result;
		}
		Comp4d c_rr_pentate					(Vect4d const &x, Vect4d const &y)
		{
			Complex1d
				lo=pentate(x.lo(), y.lo()),
				hi=pentate(x.hi(), y.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_cr_pentate					(Comp4d const &x, Vect4d const &y)
		{
			Complex1d
				lo=pentate(Complex1d(x.r.lo(), x.i.lo()), y.lo()),
				hi=pentate(Complex1d(x.r.hi(), x.i.hi()), y.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		bool disc_rr_pentate_i			(Value const &x0, Value const &y0, Value const &x1, Value const &y1){return false;}//
		bool disc_cr_pentate_i			(Value const &x0, Value const &y0, Value const &x1, Value const &y1){return false;}////*/

		void  r_r_bitwise_shift_left_l	(VectP &r, VectP const &x)					{assign(r, Vect4d(::exp((Vect4d(x).floor()*m_ln2).v)));VLEAVE;}//<<x = 2^x
		void  c_c_bitwise_shift_left_l	(CompP &r, CompP const &x)					{assign(r, exp(Comp4d(x).floor()*m_ln2));VLEAVE;}
		void  q_q_bitwise_shift_left_l	(QuatP &r, QuatP const &x)					{assign(r, exp(Quat4d(x).floor()*m_ln2));VLEAVE;}
		void  r_r_bitwise_shift_left_r	(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, rx+rx);VLEAVE;}//x<< = x*2
		void  c_c_bitwise_shift_left_r	(CompP &r, CompP const &x)					{Comp4d cx=x; assign(r, cx+cx);VLEAVE;}
		void  q_q_bitwise_shift_left_r	(QuatP &r, QuatP const &x)					{Quat4d qx=x; assign(r, qx+qx);VLEAVE;}
		void r_rr_bitwise_shift_left	(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)*Vect4d(::exp((Vect4d(y).floor()*m_ln2).v)));VLEAVE;}//x<<y = x*2^y
		void c_rc_bitwise_shift_left	(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)*exp(Comp4d(y).floor()*m_ln2));VLEAVE;}
		void q_rq_bitwise_shift_left	(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)*exp(Quat4d(y).floor()*m_ln2));VLEAVE;}
		void c_cr_bitwise_shift_left	(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)*Vect4d(::exp((Vect4d(y).floor()*m_ln2).v)));VLEAVE;}
		void c_cc_bitwise_shift_left	(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)*exp(Comp4d(y).floor()*m_ln2));VLEAVE;}
		void q_cq_bitwise_shift_left	(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)*exp(Quat4d(y).floor()*m_ln2));VLEAVE;}
		void q_qr_bitwise_shift_left	(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)*Vect4d(::exp((Vect4d(y).floor()*m_ln2).v)));VLEAVE;}
		void q_qc_bitwise_shift_left	(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)*exp(Comp4d(y).floor()*m_ln2));VLEAVE;}
		void q_qq_bitwise_shift_left	(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)*exp(Quat4d(y).floor()*m_ln2));VLEAVE;}

		void  r_r_bitwise_shift_right_l	(VectP &r, VectP const &x)					{assign(r, Vect4d(::exp((-Vect4d(x).floor()*m_ln2).v)));VLEAVE;}//>>x = 2^-x = exp(-x*ln2)
		void  c_c_bitwise_shift_right_l	(CompP &r, CompP const &x)					{assign(r, exp(-Comp4d(x).floor()*m_ln2));VLEAVE;}
		void  q_q_bitwise_shift_right_l	(QuatP &r, QuatP const &x)					{assign(r, exp(-Quat4d(x).floor()*m_ln2));VLEAVE;}
		void  r_r_bitwise_shift_right_r	(VectP &r, VectP const &x)					{assign(r, Vect4d(x)*m_half);VLEAVE;}//x>> = x/2
		void  c_c_bitwise_shift_right_r	(CompP &r, CompP const &x)					{assign(r, Comp4d(x)*m_half);VLEAVE;}
		void  q_q_bitwise_shift_right_r	(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x)*m_half);VLEAVE;}
		void r_rr_bitwise_shift_right	(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)*Vect4d(::exp((-Vect4d(y).floor()*m_ln2).v)));VLEAVE;}
		void c_rc_bitwise_shift_right	(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)*exp(-Comp4d(y).floor()*m_ln2));VLEAVE;}
		void q_rq_bitwise_shift_right	(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)*exp(-Quat4d(y).floor()*m_ln2));VLEAVE;}
		void c_cr_bitwise_shift_right	(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)*Vect4d(::exp((-Vect4d(y).floor()*m_ln2).v)));VLEAVE;}
		void c_cc_bitwise_shift_right	(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)*exp(-Comp4d(y).floor()*m_ln2));VLEAVE;}
		void q_cq_bitwise_shift_right	(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)*exp(-Quat4d(y).floor()*m_ln2));VLEAVE;}
		void q_qr_bitwise_shift_right	(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)*Vect4d(::exp((-Vect4d(y).floor()*m_ln2).v)));VLEAVE;}
		void q_qc_bitwise_shift_right	(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)*exp(-Comp4d(y).floor()*m_ln2));VLEAVE;}
		void q_qq_bitwise_shift_right	(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)*exp(-Quat4d(y).floor()*m_ln2));VLEAVE;}

	/*	__forceinline Vect4d bitwise_not(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			Vect4d xc((double)~convert_d2ll(x.lo()), (double)~convert_d2ll(x.hi()));
			return xc&mask;
		}
		Vect4d  r_r_bitwise_not				(Vect4d const &x)					{return bitwise_not(x);}
		Comp4d  c_c_bitwise_not				(Comp4d const &x)					{return Comp4d(bitwise_not(x.r), bitwise_not(x.i));}
		Quat4d  q_q_bitwise_not				(Quat4d const &x)					{return Quat4d(bitwise_not(x.r), bitwise_not(x.i), bitwise_not(x.j), bitwise_not(x.k));}

		__forceinline Vect4d bitwise_and(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			return Vect4d((double)!~convert_d2ll(x.lo()), (double)!~convert_d2ll(x.hi()))&mask;
		}
		__forceinline Vect4d bitwise_and(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect4d(double(convert_d2ll(x.lo())&convert_d2ll(y.lo())), double(convert_d2ll(x.hi())&convert_d2ll(y.hi())))&mask;
		}
		__forceinline Vect4d bitwise_and_ll(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect4d v;
			(long long&)v.lo()=convert_d2ll(x.lo())&convert_d2ll(y.lo());
			(long long&)v.hi()=convert_d2ll(x.hi())&convert_d2ll(y.hi());
			return v&mask;
		}
		__forceinline Vect4d convert_ll2d(Vect4d const &x){Vect4d x2=x; return Vect4d((double)(long long&)x2.lo(), (double)(long long&)x2.hi());}
		Vect4d  r_r_bitwise_and				(Vect4d const &x)					{return bitwise_and(x);}
		Comp4d  c_c_bitwise_and				(Comp4d const &x)					{return Comp4d(bitwise_and(x.r), bitwise_and(x.i));}
		Quat4d  q_q_bitwise_and				(Quat4d const &x)					{return Quat4d(bitwise_and(x.r), bitwise_and(x.i), bitwise_and(x.j), bitwise_and(x.k));}
		Vect4d r_rr_bitwise_and				(Vect4d const &x, Vect4d const &y)	{return bitwise_and(x, y);}
		Comp4d c_rc_bitwise_and				(Vect4d const &x, Comp4d const &y)	{return Comp4d(bitwise_and(x, y.r), bitwise_and(x, y.i));}
		Quat4d q_rq_bitwise_and				(Vect4d const &x, Quat4d const &y)	{return Quat4d(bitwise_and(x, y.r), bitwise_and(x, y.i), bitwise_and(x, y.j), bitwise_and(x, y.k));}
		Comp4d c_cr_bitwise_and				(Comp4d const &x, Vect4d const &y)	{return Comp4d(bitwise_and(x.r, y), bitwise_and(x.i, y));}
		Comp4d c_cc_bitwise_and				(Comp4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i);
			return Comp4d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr));
		}
		Quat4d q_cq_bitwise_and				(Comp4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i),
				xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j),
				xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k);
			return Quat4d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xr_yj-xi_yk), convert_ll2d(xr_yk+xi_yj));
		}
		Quat4d q_qr_bitwise_and				(Quat4d const &x, Vect4d const &y)	{return Quat4d(bitwise_and(x.r, y), bitwise_and(x.i, y), bitwise_and(x.j, y), bitwise_and(x.k, y));}
		Quat4d q_qc_bitwise_and				(Quat4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i);
			return Quat4d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xj_yr+xk_yi), convert_ll2d(-xj_yi+xk_yr));
		}
		Quat4d q_qq_bitwise_and				(Quat4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i),
				xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j), xj_yj=bitwise_and_ll(x.j, y.j), xk_yj=bitwise_and_ll(x.k, y.j),
				xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k), xj_yk=bitwise_and_ll(x.j, y.k), xk_yk=bitwise_and_ll(x.k, y.k);
			return Quat4d(convert_ll2d(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d(xr_yk+xi_yj-xj_yi+xk_yr));
		}

		__forceinline Vect4d bitwise_nand(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			return Vect4d((double)(~convert_d2ll(x.lo())!=0), (double)(~convert_d2ll(x.hi())!=0))&mask;
		}
		__forceinline Vect4d bitwise_nand(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect4d(double(~(convert_d2ll(x.lo())&convert_d2ll(y.lo()))), double(~(convert_d2ll(x.hi())&convert_d2ll(y.hi()))))&mask;
		}
		__forceinline Vect4d bitwise_nand_ll(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect4d v;
			(long long&)v.lo()=~(convert_d2ll(x.lo())&convert_d2ll(y.lo()));
			(long long&)v.hi()=~(convert_d2ll(x.hi())&convert_d2ll(y.hi()));
			return v&mask;
		}
		Vect4d  r_r_bitwise_nand			(Vect4d const &x)					{return bitwise_nand(x);}
		Comp4d  c_c_bitwise_nand			(Comp4d const &x)					{return Comp4d(bitwise_nand(x.r), bitwise_nand(x.i));}
		Quat4d  q_q_bitwise_nand			(Quat4d const &x)					{return Quat4d(bitwise_nand(x.r), bitwise_nand(x.i), bitwise_nand(x.j), bitwise_nand(x.k));}
		Vect4d r_rr_bitwise_nand			(Vect4d const &x, Vect4d const &y)	{return bitwise_nand(x, y);}
		Comp4d c_rc_bitwise_nand			(Vect4d const &x, Comp4d const &y)	{return Comp4d(bitwise_nand(x, y.r), bitwise_nand(x, y.i));}
		Quat4d q_rq_bitwise_nand			(Vect4d const &x, Quat4d const &y)	{return Quat4d(bitwise_nand(x, y.r), bitwise_nand(x, y.i), bitwise_nand(x, y.j), bitwise_nand(x, y.k));}
		Comp4d c_cr_bitwise_nand			(Comp4d const &x, Vect4d const &y)	{return Comp4d(bitwise_nand(x.r, y), bitwise_nand(x.i, y));}
		Comp4d c_cc_bitwise_nand			(Comp4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i);
			return Comp4d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr));
		}
		Quat4d q_cq_bitwise_nand			(Comp4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i),
				xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j),
				xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k);
			return Quat4d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xr_yj-xi_yk), convert_ll2d(xr_yk+xi_yj));
		}
		Quat4d q_qr_bitwise_nand			(Quat4d const &x, Vect4d const &y)	{return Quat4d(bitwise_nand(x.r, y), bitwise_nand(x.i, y), bitwise_nand(x.j, y), bitwise_nand(x.k, y));}
		Quat4d q_qc_bitwise_nand			(Quat4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i);
			return Quat4d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xj_yr+xk_yi), convert_ll2d(-xj_yi+xk_yr));
		}
		Quat4d q_qq_bitwise_nand			(Quat4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i),
				xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j), xj_yj=bitwise_nand_ll(x.j, y.j), xk_yj=bitwise_nand_ll(x.k, y.j),
				xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k), xj_yk=bitwise_nand_ll(x.j, y.k), xk_yk=bitwise_nand_ll(x.k, y.k);
			return Quat4d(convert_ll2d(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d(xr_yk+xi_yj-xj_yi+xk_yr));
		}
	
		__forceinline Vect4d bitwise_or(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			return Vect4d((double)(convert_d2ll(x.lo())!=0), (double)(convert_d2ll(x.hi())!=0))&mask;
		}
		__forceinline Vect4d bitwise_or(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect4d(double(convert_d2ll(x.lo())|convert_d2ll(y.lo())), double(convert_d2ll(x.hi())|convert_d2ll(y.hi())))&mask;
		}
		__forceinline Vect4d bitwise_or_ll_c(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect4d v;
			(long long&)v.lo()=~convert_d2ll(x.lo())&~convert_d2ll(y.lo());
			(long long&)v.hi()=~convert_d2ll(x.hi())&~convert_d2ll(y.hi());
			return v&mask;
		}
		__forceinline Vect4d convert_ll2d_c(Vect4d const &x){Vect4d x2=x; return Vect4d((double)~(long long)(long long&)x2.lo(), (double)~(long long)(long long&)x2.hi());}
		Vect4d  r_r_bitwise_or				(Vect4d const &x)					{return bitwise_or(x);}
		Comp4d  c_c_bitwise_or				(Comp4d const &x)					{return Comp4d(bitwise_or(x.r), bitwise_or(x.i));}
		Quat4d  q_q_bitwise_or				(Quat4d const &x)					{return Quat4d(bitwise_or(x.r), bitwise_or(x.i), bitwise_or(x.j), bitwise_or(x.k));}
		Vect4d r_rr_bitwise_or				(Vect4d const &x, Vect4d const &y)	{return bitwise_or(x, y);}
		Comp4d c_rc_bitwise_or				(Vect4d const &x, Comp4d const &y)	{return Comp4d(bitwise_or(x, y.r), bitwise_or(x, y.i));}
		Quat4d q_rq_bitwise_or				(Vect4d const &x, Quat4d const &y)	{return Quat4d(bitwise_or(x, y.r), bitwise_or(x, y.i), bitwise_or(x, y.j), bitwise_or(x, y.k));}
		Comp4d c_cr_bitwise_or				(Comp4d const &x, Vect4d const &y)	{return Comp4d(bitwise_or(x.r, y), bitwise_or(x.i, y));}
		Comp4d c_cc_bitwise_or				(Comp4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i);
			return Comp4d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr));
		}
		Quat4d q_cq_bitwise_or				(Comp4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i),
				xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j),
				xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k);
			return Quat4d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xr_yj-xi_yk), convert_ll2d_c(xr_yk+xi_yj));
		}
		Quat4d q_qr_bitwise_or				(Quat4d const &x, Vect4d const &y)	{return Quat4d(bitwise_or(x.r, y), bitwise_or(x.i, y), bitwise_or(x.j, y), bitwise_or(x.k, y));}
		Quat4d q_qc_bitwise_or				(Quat4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i);
			return Quat4d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xj_yr+xk_yi), convert_ll2d_c(-xj_yi+xk_yr));
		}
		Quat4d q_qq_bitwise_or				(Quat4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i),
				xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j), xj_yj=bitwise_or_ll_c(x.j, y.j), xk_yj=bitwise_or_ll_c(x.k, y.j),
				xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k), xj_yk=bitwise_or_ll_c(x.j, y.k), xk_yk=bitwise_or_ll_c(x.k, y.k);
			return Quat4d(convert_ll2d_c(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d_c(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d_c(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d_c(xr_yk+xi_yj-xj_yi+xk_yr));
		}
	
		__forceinline Vect4d bitwise_nor(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			return Vect4d((double)!convert_d2ll(x.lo()), (double)!convert_d2ll(x.hi()))&mask;
		}
		__forceinline Vect4d bitwise_nor(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect4d(double(~(convert_d2ll(x.lo())|convert_d2ll(y.lo()))), double(~(convert_d2ll(x.hi())|convert_d2ll(y.hi()))))&mask;
		}
		__forceinline Vect4d bitwise_nor_ll_c(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect4d v;
			(long long&)v.lo()=~convert_d2ll(x.lo())|~convert_d2ll(y.lo());
			(long long&)v.hi()=~convert_d2ll(x.hi())|~convert_d2ll(y.hi());
			return v&mask;
		}
		Vect4d  r_r_bitwise_nor				(Vect4d const &x)					{return bitwise_nor(x);}
		Comp4d  c_c_bitwise_nor				(Comp4d const &x)					{return Comp4d(bitwise_nor(x.r), bitwise_nor(x.i));}
		Quat4d  q_q_bitwise_nor				(Quat4d const &x)					{return Quat4d(bitwise_nor(x.r), bitwise_nor(x.i), bitwise_nor(x.j), bitwise_nor(x.k));}
		Vect4d r_rr_bitwise_nor				(Vect4d const &x, Vect4d const &y)	{return bitwise_nor(x, y);}
		Comp4d c_rc_bitwise_nor				(Vect4d const &x, Comp4d const &y)	{return Comp4d(bitwise_nor(x, y.r), bitwise_nor(x, y.i));}
		Quat4d q_rq_bitwise_nor				(Vect4d const &x, Quat4d const &y)	{return Quat4d(bitwise_nor(x, y.r), bitwise_nor(x, y.i), bitwise_nor(x, y.j), bitwise_nor(x, y.k));}
		Comp4d c_cr_bitwise_nor				(Comp4d const &x, Vect4d const &y)	{return Comp4d(bitwise_nor(x.r, y), bitwise_nor(x.i, y));}
		Comp4d c_cc_bitwise_nor				(Comp4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i);
			return Comp4d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr));
		}
		Quat4d q_cq_bitwise_nor				(Comp4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i),
				xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j),
				xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k);
			return Quat4d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xr_yj-xi_yk), convert_ll2d_c(xr_yk+xi_yj));
		}
		Quat4d q_qr_bitwise_nor				(Quat4d const &x, Vect4d const &y)	{return Quat4d(bitwise_nor(x.r, y), bitwise_nor(x.i, y), bitwise_nor(x.j, y), bitwise_nor(x.k, y));}
		Quat4d q_qc_bitwise_nor				(Quat4d const &x, Comp4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i);
			return Quat4d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xj_yr+xk_yi), convert_ll2d_c(-xj_yi+xk_yr));
		}
		Quat4d q_qq_bitwise_nor				(Quat4d const &x, Quat4d const &y)
		{
			Vect4d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i),
				xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j), xj_yj=bitwise_nor_ll_c(x.j, y.j), xk_yj=bitwise_nor_ll_c(x.k, y.j),
				xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k), xj_yk=bitwise_nor_ll_c(x.j, y.k), xk_yk=bitwise_nor_ll_c(x.k, y.k);
			return Quat4d(convert_ll2d_c(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d_c(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d_c(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d_c(xr_yk+xi_yj-xj_yi+xk_yr));
		}
	
		__forceinline Vect4d bitwise_xor(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			return Vect4d((double)G2::bitwise_xor(x.lo()), (double)bitwise_xor(x.hi()))&mask;
		}
		__forceinline Vect4d bitwise_xor(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()&is_not_nan(y)&isinf(y).complement();
			return Vect4d(double(convert_d2ll(x.lo())^convert_d2ll(y.lo())), double(convert_d2ll(x.hi())^convert_d2ll(y.hi())))&mask;
		}
		Vect4d  r_r_bitwise_xor				(Vect4d const &x)					{return bitwise_xor(x);}
		Comp4d  c_c_bitwise_xor				(Comp4d const &x)					{return Comp4d(bitwise_xor(x.r), bitwise_xor(x.i));}
		Quat4d  q_q_bitwise_xor				(Quat4d const &x)					{return Quat4d(bitwise_xor(x.r), bitwise_xor(x.i), bitwise_xor(x.j), bitwise_xor(x.k));}
		Vect4d r_rr_bitwise_xor				(Vect4d const &x, Vect4d const &y)	{return bitwise_xor(x, y);}
		Comp4d c_rc_bitwise_xor				(Vect4d const &x, Comp4d const &y)	{return Comp4d(bitwise_xor(x, y.r), y.i);}
		Quat4d q_rq_bitwise_xor				(Vect4d const &x, Quat4d const &y)	{return Quat4d(bitwise_xor(x, y.r), y.i, y.j, y.k);}
		Comp4d c_cr_bitwise_xor				(Comp4d const &x, Vect4d const &y)	{return Comp4d(bitwise_xor(x.r, y), x.i);}
		Comp4d c_cc_bitwise_xor				(Comp4d const &x, Comp4d const &y)	{return Comp4d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i));}
		Quat4d q_cq_bitwise_xor				(Comp4d const &x, Quat4d const &y)	{return Quat4d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), y.j, y.k);}
		Quat4d q_qr_bitwise_xor				(Quat4d const &x, Vect4d const &y)	{return Quat4d(bitwise_xor(x.r, y), x.i, x.j, x.k);}
		Quat4d q_qc_bitwise_xor				(Quat4d const &x, Comp4d const &y)	{return Quat4d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), x.j, x.k);}
		Quat4d q_qq_bitwise_xor				(Quat4d const &x, Quat4d const &y)	{return Quat4d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), bitwise_xor(x.j, y.j), bitwise_xor(x.k, y.k));}
	
		__forceinline Vect4d bitwise_xnor(Vect4d const &x)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement();
			return Vect4d((double)!bitwise_xor(x.lo()), (double)!bitwise_xor(x.hi()))&mask;
		}
		__forceinline Vect4d bitwise_xnor(Vect4d const &x, Vect4d const &y)
		{
			Vect4d mask=is_not_nan(x)&isinf(x).complement()&is_not_nan(y)&isinf(y).complement();
			return Vect4d((double)~(convert_d2ll(x.lo())^convert_d2ll(y.lo())), (double)~(convert_d2ll(x.hi())^convert_d2ll(y.hi())))&mask;
		}
		Vect4d  r_r_bitwise_xnor			(Vect4d const &x)					{return bitwise_xnor(x);}
		Comp4d  c_c_bitwise_xnor			(Comp4d const &x)					{return Comp4d(bitwise_xnor(x.r), bitwise_xnor(x.i));}
		Quat4d  q_q_bitwise_xnor			(Quat4d const &x)					{return Quat4d(bitwise_xnor(x.r), bitwise_xnor(x.i), bitwise_xnor(x.j), bitwise_xnor(x.k));}
		Vect4d r_rr_bitwise_xnor			(Vect4d const &x, Vect4d const &y)	{return bitwise_xnor(x, y);}
		Comp4d c_rc_bitwise_xnor			(Vect4d const &x, Comp4d const &y)	{return Comp4d(bitwise_xnor(x, y.r), y.i);}
		Quat4d q_rq_bitwise_xnor			(Vect4d const &x, Quat4d const &y)	{return Quat4d(bitwise_xnor(x, y.r), y.i, y.j, y.k);}
		Comp4d c_cr_bitwise_xnor			(Comp4d const &x, Vect4d const &y)	{return Comp4d(bitwise_xnor(x.r, y), x.i);}
		Comp4d c_cc_bitwise_xnor			(Comp4d const &x, Comp4d const &y)	{return Comp4d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i));}
		Quat4d q_cq_bitwise_xnor			(Comp4d const &x, Quat4d const &y)	{return Quat4d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), y.j, y.k);}
		Quat4d q_qr_bitwise_xnor			(Quat4d const &x, Vect4d const &y)	{return Quat4d(bitwise_xnor(x.r, y), x.i, x.j, x.k);}
		Quat4d q_qc_bitwise_xnor			(Quat4d const &x, Comp4d const &y)	{return Quat4d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), x.j, x.k);}
		Quat4d q_qq_bitwise_xnor			(Quat4d const &x, Quat4d const &y)	{return Quat4d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), bitwise_xnor(x.j, y.j), bitwise_xnor(x.k, y.k));}//*/
	
		void  r_r_logic_equal			(VectP &r, VectP const &x)					{assign(r, (Vect4d(x)==m_zero)&m_one);VLEAVE;}
		void  r_c_logic_equal			(VectP &r, CompP const &x)					{assign(r, Comp4d(x).c_is_true().complement()&m_one);VLEAVE;}
		void  r_q_logic_equal			(VectP &r, QuatP const &x)					{assign(r, Quat4d(x).q_is_true().complement()&m_one);VLEAVE;}
		void r_rr_logic_equal			(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x)==Vect4d(y))&m_one);VLEAVE;}
		void r_rc_logic_equal			(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x)==Comp4d(y))&m_one);VLEAVE;}
		void r_rq_logic_equal			(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x)==Quat4d(y))&m_one);VLEAVE;}
		void r_cr_logic_equal			(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp4d(x)==Vect4d(y))&m_one);VLEAVE;}
		void r_cc_logic_equal			(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp4d(x)==Comp4d(y))&m_one);VLEAVE;}
		void r_cq_logic_equal			(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp4d(x)==Quat4d(y))&m_one);VLEAVE;}
		void r_qr_logic_equal			(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat4d(x)==Vect4d(y))&m_one);VLEAVE;}
		void r_qc_logic_equal			(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat4d(x)==Comp4d(y))&m_one);VLEAVE;}
		void r_qq_logic_equal			(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat4d(x)==Quat4d(y))&m_one);VLEAVE;}
	
		void  r_r_logic_not_equal		(VectP &r, VectP const &x)					{assign(r, Vect4d(x).r_is_true()&m_one);VLEAVE;}
		void  r_c_logic_not_equal		(VectP &r, CompP const &x)					{assign(r, Comp4d(x).c_is_true()&m_one);VLEAVE;}
		void  r_q_logic_not_equal		(VectP &r, QuatP const &x)					{assign(r, Quat4d(x).q_is_true()&m_one);VLEAVE;}
		void r_rr_logic_not_equal		(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x)!=Vect4d(y))&m_one);VLEAVE;}
		void r_rc_logic_not_equal		(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x)!=Comp4d(y))&m_one);VLEAVE;}
		void r_rq_logic_not_equal		(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x)!=Quat4d(y))&m_one);VLEAVE;}
		void r_cr_logic_not_equal		(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp4d(x)!=Vect4d(y))&m_one);VLEAVE;}
		void r_cc_logic_not_equal		(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp4d(x)!=Comp4d(y))&m_one);VLEAVE;}
		void r_cq_logic_not_equal		(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp4d(x)!=Quat4d(y))&m_one);VLEAVE;}
		void r_qr_logic_not_equal		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat4d(x)!=Vect4d(y))&m_one);VLEAVE;}
		void r_qc_logic_not_equal		(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat4d(x)!=Comp4d(y))&m_one);VLEAVE;}
		void r_qq_logic_not_equal		(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat4d(x)!=Quat4d(y))&m_one);VLEAVE;}
	
		void  r_r_logic_less_l			(VectP &r, VectP const &x)					{assign(r, (m_zero<Vect4d(x))&m_one);VLEAVE;}
		void  r_c_logic_less_l			(VectP &r, CompP const &x)					{assign(r, (m_zero<Vect4d(x.r))&m_one);VLEAVE;}
		void  r_q_logic_less_l			(VectP &r, QuatP const &x)					{assign(r, (m_zero<Vect4d(x.r))&m_one);VLEAVE;}
		void  r_r_logic_less_r			(VectP &r, VectP const &x)					{assign(r, (Vect4d(x)<m_zero)&m_one);VLEAVE;}
		void  r_c_logic_less_r			(VectP &r, CompP const &x)					{assign(r, (Vect4d(x.r)<m_zero)&m_one);VLEAVE;}
		void  r_q_logic_less_r			(VectP &r, QuatP const &x)					{assign(r, (Vect4d(x.r)<m_zero)&m_one);VLEAVE;}
		void r_rr_logic_less			(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x)<Vect4d(y))&m_one);VLEAVE;}
		void r_rc_logic_less			(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x)<Vect4d(y.r))&m_one);VLEAVE;}
		void r_rq_logic_less			(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x)<Vect4d(y.r))&m_one);VLEAVE;}
		void r_cr_logic_less			(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)<Vect4d(y))&m_one);VLEAVE;}
		void r_cc_logic_less			(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)<Vect4d(y.r))&m_one);VLEAVE;}
		void r_cq_logic_less			(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)<Vect4d(y.r))&m_one);VLEAVE;}
		void r_qr_logic_less			(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)<Vect4d(y))&m_one);VLEAVE;}
		void r_qc_logic_less			(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)<Vect4d(y.r))&m_one);VLEAVE;}
		void r_qq_logic_less			(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)<Vect4d(y.r))&m_one);VLEAVE;}
	
		void  r_r_logic_less_equal_l	(VectP &r, VectP const &x)					{assign(r, (m_zero<=Vect4d(x))&m_one);VLEAVE;}
		void  r_c_logic_less_equal_l	(VectP &r, CompP const &x)					{assign(r, (m_zero<=Vect4d(x.r))&m_one);VLEAVE;}
		void  r_q_logic_less_equal_l	(VectP &r, QuatP const &x)					{assign(r, (m_zero<=Vect4d(x.r))&m_one);VLEAVE;}
		void  r_r_logic_less_equal_r	(VectP &r, VectP const &x)					{assign(r, (Vect4d(x)<=m_zero)&m_one);VLEAVE;}
		void  r_c_logic_less_equal_r	(VectP &r, CompP const &x)					{assign(r, (Vect4d(x.r)<=m_zero)&m_one);VLEAVE;}
		void  r_q_logic_less_equal_r	(VectP &r, QuatP const &x)					{assign(r, (Vect4d(x.r)<=m_zero)&m_one);VLEAVE;}
		void r_rr_logic_less_equal		(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x)<=Vect4d(y))&m_one);VLEAVE;}
		void r_rc_logic_less_equal		(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x)<=Vect4d(y.r))&m_one);VLEAVE;}
		void r_rq_logic_less_equal		(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x)<=Vect4d(y.r))&m_one);VLEAVE;}
		void r_cr_logic_less_equal		(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)<=Vect4d(y))&m_one);VLEAVE;}
		void r_cc_logic_less_equal		(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)<=Vect4d(y.r))&m_one);VLEAVE;}
		void r_cq_logic_less_equal		(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)<=Vect4d(y.r))&m_one);VLEAVE;}
		void r_qr_logic_less_equal		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)<=Vect4d(y))&m_one);VLEAVE;}
		void r_qc_logic_less_equal		(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)<=Vect4d(y.r))&m_one);VLEAVE;}
		void r_qq_logic_less_equal		(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)<=Vect4d(y.r))&m_one);VLEAVE;}
	
		void  r_r_logic_greater_l		(VectP &r, VectP const &x)					{assign(r, (m_zero>Vect4d(x))&m_one);VLEAVE;}
		void  r_c_logic_greater_l		(VectP &r, CompP const &x)					{assign(r, (m_zero>Vect4d(x.r))&m_one);VLEAVE;}
		void  r_q_logic_greater_l		(VectP &r, QuatP const &x)					{assign(r, (m_zero>Vect4d(x.r))&m_one);VLEAVE;}
		void  r_r_logic_greater_r		(VectP &r, VectP const &x)					{assign(r, (Vect4d(x)>m_zero)&m_one);VLEAVE;}
		void  r_c_logic_greater_r		(VectP &r, CompP const &x)					{assign(r, (Vect4d(x.r)>m_zero)&m_one);VLEAVE;}
		void  r_q_logic_greater_r		(VectP &r, QuatP const &x)					{assign(r, (Vect4d(x.r)>m_zero)&m_one);VLEAVE;}
		void r_rr_logic_greater			(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x)>Vect4d(y))&m_one);VLEAVE;}
		void r_rc_logic_greater			(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x)>Vect4d(y.r))&m_one);VLEAVE;}
		void r_rq_logic_greater			(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x)>Vect4d(y.r))&m_one);VLEAVE;}
		void r_cr_logic_greater			(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)>Vect4d(y))&m_one);VLEAVE;}
		void r_cc_logic_greater			(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)>Vect4d(y.r))&m_one);VLEAVE;}
		void r_cq_logic_greater			(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)>Vect4d(y.r))&m_one);VLEAVE;}
		void r_qr_logic_greater			(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)>Vect4d(y))&m_one);VLEAVE;}
		void r_qc_logic_greater			(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)>Vect4d(y.r))&m_one);VLEAVE;}
		void r_qq_logic_greater			(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)>Vect4d(y.r))&m_one);VLEAVE;}
	
		void  r_r_logic_greater_equal_l	(VectP &r, VectP const &x)					{assign(r, (m_zero>=Vect4d(x))&m_one);VLEAVE;}
		void  r_c_logic_greater_equal_l	(VectP &r, CompP const &x)					{assign(r, (m_zero>=Vect4d(x.r))&m_one);VLEAVE;}
		void  r_q_logic_greater_equal_l	(VectP &r, QuatP const &x)					{assign(r, (m_zero>=Vect4d(x.r))&m_one);VLEAVE;}
		void  r_r_logic_greater_equal_r	(VectP &r, VectP const &x)					{assign(r, (Vect4d(x)>=m_zero)&m_one);VLEAVE;}
		void  r_c_logic_greater_equal_r	(VectP &r, CompP const &x)					{assign(r, (Vect4d(x.r)>=m_zero)&m_one);VLEAVE;}
		void  r_q_logic_greater_equal_r	(VectP &r, QuatP const &x)					{assign(r, (Vect4d(x.r)>=m_zero)&m_one);VLEAVE;}
		void r_rr_logic_greater_equal	(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x)>=Vect4d(y))&m_one);VLEAVE;}
		void r_rc_logic_greater_equal	(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x)>=Vect4d(y.r))&m_one);VLEAVE;}
		void r_rq_logic_greater_equal	(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x)>=Vect4d(y.r))&m_one);VLEAVE;}
		void r_cr_logic_greater_equal	(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)>=Vect4d(y))&m_one);VLEAVE;}
		void r_cc_logic_greater_equal	(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)>=Vect4d(y.r))&m_one);VLEAVE;}
		void r_cq_logic_greater_equal	(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)>=Vect4d(y.r))&m_one);VLEAVE;}
		void r_qr_logic_greater_equal	(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect4d(x.r)>=Vect4d(y))&m_one);VLEAVE;}
		void r_qc_logic_greater_equal	(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect4d(x.r)>=Vect4d(y.r))&m_one);VLEAVE;}
		void r_qq_logic_greater_equal	(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect4d(x.r)>=Vect4d(y.r))&m_one);VLEAVE;}
	
		void  r_r_logic_not				(VectP &r, VectP const &x)					{assign(r, (Vect4d(x)==m_zero)&m_one);VLEAVE;}
		void  r_c_logic_not				(VectP &r, CompP const &x)					{Comp4d cx=x; assign(r, (cx.r==m_zero)&(cx.i==m_zero)&m_one);VLEAVE;}
		void  r_q_logic_not				(VectP &r, QuatP const &x)					{Quat4d qx=x; assign(r, (qx.r==m_zero)&(qx.i==m_zero)&(qx.j==m_zero)&(qx.k==m_zero)&m_one);VLEAVE;}
	
		void r_rr_logic_and				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x).r_is_true()&Vect4d(y).r_is_true()&m_one);VLEAVE;}
		void r_rc_logic_and				(VectP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x).r_is_true()&Comp4d(y).c_is_true()&m_one);VLEAVE;}
		void r_rq_logic_and				(VectP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x).r_is_true()&Quat4d(y).q_is_true()&m_one);VLEAVE;}
		void r_cr_logic_and				(VectP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x).c_is_true()&Vect4d(y).r_is_true()&m_one);VLEAVE;}
		void r_cc_logic_and				(VectP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x).c_is_true()&Comp4d(y).c_is_true()&m_one);VLEAVE;}
		void r_cq_logic_and				(VectP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x).c_is_true()&Quat4d(y).q_is_true()&m_one);VLEAVE;}
		void r_qr_logic_and				(VectP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x).q_is_true()&Vect4d(y).r_is_true()&m_one);VLEAVE;}
		void r_qc_logic_and				(VectP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x).q_is_true()&Comp4d(y).c_is_true()&m_one);VLEAVE;}
		void r_qq_logic_and				(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x).q_is_true()&Quat4d(y).q_is_true()&m_one);VLEAVE;}
	
		void r_rr_logic_or				(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x).r_is_true()|Vect4d(y).r_is_true())&m_one);VLEAVE;}
		void r_rc_logic_or				(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x).r_is_true()|Comp4d(y).c_is_true())&m_one);VLEAVE;}
		void r_rq_logic_or				(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x).r_is_true()|Quat4d(y).q_is_true())&m_one);VLEAVE;}
		void r_cr_logic_or				(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp4d(x).c_is_true()|Vect4d(y).r_is_true())&m_one);VLEAVE;}
		void r_cc_logic_or				(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp4d(x).c_is_true()|Comp4d(y).c_is_true())&m_one);VLEAVE;}
		void r_cq_logic_or				(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp4d(x).c_is_true()|Quat4d(y).q_is_true())&m_one);VLEAVE;}
		void r_qr_logic_or				(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat4d(x).q_is_true()|Vect4d(y).r_is_true())&m_one);VLEAVE;}
		void r_qc_logic_or				(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat4d(x).q_is_true()|Comp4d(y).c_is_true())&m_one);VLEAVE;}
		void r_qq_logic_or				(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat4d(x).q_is_true()|Quat4d(y).q_is_true())&m_one);VLEAVE;}
	
		void r_rr_logic_xor				(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect4d(x).r_is_true()^Vect4d(y).r_is_true())&m_one);VLEAVE;}
		void r_rc_logic_xor				(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect4d(x).r_is_true()^Comp4d(y).c_is_true())&m_one);VLEAVE;}
		void r_rq_logic_xor				(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect4d(x).r_is_true()^Quat4d(y).q_is_true())&m_one);VLEAVE;}
		void r_cr_logic_xor				(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp4d(x).c_is_true()^Vect4d(y).r_is_true())&m_one);VLEAVE;}
		void r_cc_logic_xor				(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp4d(x).c_is_true()^Comp4d(y).c_is_true())&m_one);VLEAVE;}
		void r_cq_logic_xor				(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp4d(x).c_is_true()^Quat4d(y).q_is_true())&m_one);VLEAVE;}
		void r_qr_logic_xor				(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat4d(x).q_is_true()^Vect4d(y).r_is_true())&m_one);VLEAVE;}
		void r_qc_logic_xor				(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat4d(x).q_is_true()^Comp4d(y).c_is_true())&m_one);VLEAVE;}
		void r_qq_logic_xor				(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat4d(x).q_is_true()^Quat4d(y).q_is_true())&m_one);VLEAVE;}
	
		void r_rr_condition_zero	(VectP &r, VectP const &x, VectP const &y)	{Vect4d rx=x, mask=rx.r_is_true();	assign(r, rx&mask|Vect4d(y)&mask.complement());VLEAVE;}
		void c_rc_condition_zero	(CompP &r, VectP const &x, CompP const &y)	{Vect4d rx=x, mask=rx.r_is_true();	assign(r, rx&mask|and(Comp4d(y), mask.complement()));VLEAVE;}
		void q_rq_condition_zero	(QuatP &r, VectP const &x, QuatP const &y)	{Vect4d rx=x, mask=rx.r_is_true();	assign(r, rx&mask|and(Quat4d(y), mask.complement()));VLEAVE;}
		void c_cr_condition_zero	(CompP &r, CompP const &x, VectP const &y)	{Comp4d cx=x; Vect4d mask=cx.c_is_true();	assign(r, and(cx, mask)|Vect4d(y)&mask.complement());VLEAVE;}
		void c_cc_condition_zero	(CompP &r, CompP const &x, CompP const &y)	{Comp4d cx=x; Vect4d mask=cx.c_is_true();	assign(r, and(cx, mask)|and(Comp4d(y), mask.complement()));VLEAVE;}
		void q_cq_condition_zero	(QuatP &r, CompP const &x, QuatP const &y)	{Comp4d cx=x; Vect4d mask=cx.c_is_true();	assign(r, and(cx, mask)|and(Quat4d(y), mask.complement()));VLEAVE;}
		void q_qr_condition_zero	(QuatP &r, QuatP const &x, VectP const &y)	{Quat4d qx=x; Vect4d mask=qx.q_is_true();	assign(r, and(qx, mask)|Vect4d(y)&mask.complement());VLEAVE;}
		void q_qc_condition_zero	(QuatP &r, QuatP const &x, CompP const &y)	{Quat4d qx=x; Vect4d mask=qx.q_is_true();	assign(r, and(qx, mask)|and(Comp4d(y), mask.complement()));VLEAVE;}
		void q_qq_condition_zero	(QuatP &r, QuatP const &x, QuatP const &y)	{Quat4d qx=x; Vect4d mask=qx.q_is_true();	assign(r, and(qx, mask)|and(Quat4d(y), mask.complement()));VLEAVE;}

		void  r_r_percent				(VectP &r, VectP const &x)					{assign(r, Vect4d(x)*m_one_percent);VLEAVE;}
		void  c_c_percent				(CompP &r, CompP const &x)					{assign(r, Comp4d(x)*m_one_percent);VLEAVE;}
		void  q_q_percent				(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x)*m_one_percent);VLEAVE;}
	
		void r_rr_modulo				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(x)%Vect4d(y));VLEAVE;}
		void c_rc_modulo				(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect4d(x)%Comp4d(y));VLEAVE;}
		void q_rq_modulo				(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect4d(x)%Quat4d(y));VLEAVE;}
		void c_cr_modulo				(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp4d(x)%Vect4d(y));VLEAVE;}
		void c_cc_modulo				(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp4d(x)%Comp4d(y));VLEAVE;}
		void q_cq_modulo				(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp4d(x)%Quat4d(y));VLEAVE;}
		void q_qr_modulo				(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat4d(x)%Vect4d(y));VLEAVE;}
		void q_qc_modulo				(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat4d(x)%Comp4d(y));VLEAVE;}
		void q_qq_modulo				(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat4d(x)%Quat4d(y));VLEAVE;}

		void  r_r_sgn					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, ((rx>m_zero)&m_one)-((rx<m_zero)&m_one));VLEAVE;}
		void  c_c_sgn					(CompP &r, CompP const &x)					{Comp4d cx=x; Vect4d mag=cx.abs(); assign(r, and((cx/mag), mag.r_is_true()));VLEAVE;}
		void  q_q_sgn					(QuatP &r, QuatP const &x)					{Quat4d qx=x; Vect4d mag=qx.abs(); assign(r, and((qx/mag), mag.r_is_true()));VLEAVE;}
		
		__forceinline Comp4d sq(Comp4d const &x){Vect4d ri=x.r*x.i; return Comp4d(x.r*x.r-x.i*x.i, ri+ri);}
		//__forceinline Comp4d sq(Comp4d const &x){return Comp4d(x.r*x.r-x.i*x.i, m_two*x.r*x.i);}
		__forceinline Quat4d sq(Quat4d const &x)
		{
			auto _2r=x.r+x.r;
			return Quat4d(x.r*x.r-x.i*x.i-x.j*x.j-x.k*x.k, x.i*_2r, x.j*_2r, x.k*_2r);
		}
		void  r_r_sq					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, rx*rx);VLEAVE;}
		void  c_c_sq					(CompP &r, CompP const &x)					{assign(r, sq(Comp4d(x)));VLEAVE;}
		void  q_q_sq					(QuatP &r, QuatP const &x)					{assign(r, sq(Quat4d(x)));VLEAVE;}

		void  c_c_sqrt					(CompP &r, CompP const &x)					{assign(r, sqrt(Comp4d(x)));VLEAVE;}
		void  q_q_sqrt					(QuatP &r, QuatP const &x)					{assign(r, sqrt(Quat4d(x)));VLEAVE;}

		void  r_r_invsqrt				(VectP &r, VectP const &x)
		{
			Vect4d rx=x;
			const Vect4d m_one_and_half=_mm256_set1_pd(1.5);
#ifdef __AVX2__
			__m256i mask=_mm256_set_epi32(
				0x5FE6EC85, 0xE7DE30DA,
				0x5FE6EC85, 0xE7DE30DA,
				0x5FE6EC85, 0xE7DE30DA,
				0x5FE6EC85, 0xE7DE30DA);
			Vect4d t0=_mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(rx.v), 1));
			t0.v=_mm256_castsi256_pd(_mm256_sub_epi64(mask, _mm256_castpd_si256(t0.v)));
#else
			const long long mask=0x5FE6EC85E7DE30DA;
			double t[4];
			(long long&)t[0]=mask-((long long&)x.r[0]>>1);
			(long long&)t[1]=mask-((long long&)x.r[1]>>1);
			(long long&)t[2]=mask-((long long&)x.r[2]>>1);
			(long long&)t[3]=mask-((long long&)x.r[3]>>1);
			Vect4d t0=t;
#endif
			assign(r, t0*(m_one_and_half-m_half*rx*t0*t0));
			VLEAVE;
		}

		void  r_r_cbrt					(VectP &r, VectP const &x)					{assign(r, Vect4d(::cbrt(Vect4d(x).v)));VLEAVE;}
		void  c_c_cbrt					(CompP &r, CompP const &x)					{assign(r, exp(m_third*log(Comp4d(x))));VLEAVE;}//optimize
		void  q_q_cbrt					(QuatP &r, QuatP const &x)
		{
			Quat4d qx=x;
			Vect4d mag_v=qx.i*qx.i+qx.j*qx.j+qx.k*qx.k;

			Vect4d real=mag_v==m_zero, real_c=real.complement(), rr=cbrt(qx.r.v);

			Vect4d mag_x=sqrt(qx.r*qx.r+mag_v);
			mag_v=sqrt(mag_v);
			Vect4d u_mul=Vect4d(::acos((qx.r/mag_x).v))/(mag_v*m_third);
			Vect4d ln_mag_x=::log(mag_x.v);
			Quat4d result(ln_mag_x*m_third, m_pi&real|qx.i*u_mul, qx.j*u_mul, qx.k*u_mul);
			result=exp(result);
			assign(r, Quat4d(rr&real|result.r&real_c, result.i&real_c, result.j&real_c, result.k&real_c));

		//	assign(r, exp(m_third*log(Quat4d(x))));
			VLEAVE;
		}

		void  r_r_gauss					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, Vect4d(::exp((-rx*rx).v)));VLEAVE;}
		void  c_c_gauss					(CompP &r, CompP const &x)					{assign(r, exp(-sq(Comp4d(x))));VLEAVE;}
		void  q_q_gauss					(QuatP &r, QuatP const &x)					{assign(r, exp(-sq(Quat4d(x))));VLEAVE;}

	/*	void  r_r_erf					(VectP &r, VectP const &x)					{assign(r, Vect4d(boost::math::erf(x.lo()), boost::math::erf(x.hi()));}

		__forceinline double zeta(double x)
		{
			try
			{
				return boost::math::zeta(x);
			}
			catch(...)
			{
				return _qnan;
			}
		}
		Vect4d  r_r_zeta					(Vect4d const &x)					{return Vect4d(zeta(x.lo()), zeta(x.hi()));}

		Vect4d  r_r_tgamma					(Vect4d const &x)					{return Vect4d(tgamma(x.lo()), tgamma(x.hi()));}
		Comp4d  c_c_tgamma					(Comp4d const &x)
		{
			Complex1d lo(x.r.lo(), x.i.lo());	lo=tgamma(lo);
			Complex1d hi(x.r.hi(), x.i.hi());	hi=tgamma(hi);
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Quat4d  q_q_tgamma					(Quat4d const &x)
		{
			Quaternion1d lo(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo());	lo=tgamma(lo);
			Quaternion1d hi(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi());	hi=tgamma(hi);
			return Quat4d(Vect4d(lo.R_component_1(), hi.R_component_1()), Vect4d(lo.R_component_2(), hi.R_component_2()),
				Vect4d(lo.R_component_3(), hi.R_component_3()), Vect4d(lo.R_component_4(), hi.R_component_4()));
		}
		double tgamma(double x, double y)//upper incomplete gamma, G(x, 0) = G(x)
		{
			try
			{
				if(x>0)
				{
					if(y>=0)
						return boost::math::tgamma((double)x, (double)y);
				}
				else if(x==0)
					return _HUGE;
				long long cInf=0x7FF8000000000010;//x.r>0&&y.r<0||x.r<0
				return (double&)cInf;
				if(x.r<0)
				//{
				//	long long cInf=0x7FF8000000000010;
				//	return (double&)cInf;
				//}
				//else if(x.r==0)
				//	return _HUGE;
				//else if(y.r<0)
				//{
				//	long long cInf=0x7FF8000000000010;
				//	return (double&)cInf;
				//}
				//return boost::math::tgamma((double)x, (double)y);
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return _HUGE;
				if(valueBits==0x8000000000000000)
					return -_HUGE;
				long long divergentBits=0x7FF8000000000010;
				return (double&)divergentBits;
			}
			catch(std::overflow_error&)
			{
				return _HUGE;
			}
		}
		Vect4d r_rr_tgamma					(Vect4d const &x, Vect4d const &y)	{return Vect4d(tgamma(x.lo(), y.lo()), tgamma(x.hi(), y.hi()));}
		bool disc_r_tgamma_i			(Value const &x0, Value const &x1){return (x0.r<=0||x1.r<=0)&&_1d_int_in_range(x0.r, x1.r);}
		bool disc_c_tgamma_i			(Value const &x0, Value const &x1)
		{
			if(x0.r==x1.r)	return false;
			if(x0.i==x1.i)	return x0.i==0&&_1d_int_in_range(x0.r, x1.r);
			if(std::signbit(x0.i)!=std::signbit(x1.i))
			{
				double t=_1d_zero_crossing(x0.r, x0.i, x1.r, x1.i);
				return t<=0&&t==std::floor(t);
			}
			return false;
		}
		bool disc_q_tgamma_i			(Value const &x0, Value const &x1){return false;}//
		bool disc_rr_tgamma_i			(Value const &x0, Value const &y0, Value const &x1, Value const &y1){return false;}////*/

		void  r_r_loggamma				(VectP &r, VectP const &x)
		{
			__declspec(align(16)) const Vect4d coef[6]={_mm256_set1_pd(76.18009172947146), _mm256_set1_pd(-86.50532032941677), _mm256_set1_pd(24.01409824083091),//https://jamesmccaffrey.wordpress.com/2013/06/19/the-log-gamma-function-with-c/
				_mm256_set1_pd(-1.231739572450155), _mm256_set1_pd(0.1208650973866179e-2), _mm256_set1_pd(-0.5395239384953e-5)};
			Vect4d rx=x, denom=rx+m_one, y=rx+Vect4d(_mm256_set1_pd(5.5)), series=Vect4d(_mm256_set1_pd(1.000000000190015));
			for(int i=0;i<6;++i)
			{
				series+=coef[i]/denom;
				denom+=m_one;
			}
			assign(r, Vect4d(_mm256_set1_pd(_ln_sqrt_2pi))+(rx+m_half)*Vect4d(::log(y.v))-y+Vect4d(::log((series/rx).v)));
			VLEAVE;
		}

	/*	Vect4d  r_r_factorial				(Vect4d const &x)					{return Vect4d(tgamma(x.lo()+1), tgamma(x.hi()+1));}
		Comp4d  c_c_factorial				(Comp4d const &x)
		{
			Complex1d lo(x.r.lo()+1, x.i.lo());	lo=tgamma(lo);
			Complex1d hi(x.r.hi()+1, x.i.hi());	hi=tgamma(hi);
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Quat4d  q_q_factorial				(Quat4d const &x)
		{
			Quaternion1d lo(x.r.lo()+1, x.i.lo(), x.j.lo(), x.k.lo());	lo=tgamma(lo);
			Quaternion1d hi(x.r.hi()+1, x.i.hi(), x.j.hi(), x.k.hi());	hi=tgamma(hi);
			return Quat4d(Vect4d(lo.R_component_1(), hi.R_component_1()), Vect4d(lo.R_component_2(), hi.R_component_2()),
				Vect4d(lo.R_component_3(), hi.R_component_3()), Vect4d(lo.R_component_4(), hi.R_component_4()));
		}
		bool disc_r_factorial_i			(Value const &x0, Value const &x1)
		{
			Value _x0=x0.r+1, _x1=x1.r+1;
			return disc_r_tgamma_i(_x0, _x1);
		}
		bool disc_c_factorial_i			(Value const &x0, Value const &x1)
		{
			Value _x0=x0, _x1=x1.r;
			_x0.r+=1, _x1.r+=1;
			return disc_c_tgamma_i(_x0, _x1);
		}
		bool disc_q_factorial_i			(Value const &x0, Value const &x1)
		{
			Value _x0=x0, _x1=x1.r;
			_x0.r+=1, _x1.r+=1;
			return disc_q_tgamma_i(_x0, _x1);
		}

		__forceinline double permutation(double x, double y)
		{
			try
			{
				double rx=x, ry=y;
				return boost::math::tgamma(rx+1)/boost::math::tgamma(rx-ry+1);
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return _HUGE;
				if(valueBits==0x8000000000000000)
					return -_HUGE;
				long long divergentBits=0x7FF8000000000010;
				return *(double*)&divergentBits;
			}
			catch(std::overflow_error&)
			{
				return _HUGE;
			}
		}
		Complex1d permutation(Complex1d const &x, Complex1d const &y)
		{
			try
			{
				std::complex<double> cx=x, cy=y;
				return tgamma(cx+1.)/tgamma(cx-cy+1.);
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return _HUGE;
				if(valueBits==0x8000000000000000)
					return -_HUGE;
				long long divergentBits=0x7FF8000000000010;
				return *(double*)&divergentBits;
			}
			catch(std::overflow_error&)
			{
				return _HUGE;
			}
		}
		Quaternion1d permutation(Quaternion1d const &x, Quaternion1d const &y)
		{
			try
			{
				boost::math::quaternion<double> qx=x, qy=y;
				return tgamma(qx+1.)/tgamma(qx-qy+1.);
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return Quaternion1d(_HUGE);
				if(valueBits==0x8000000000000000)
					return Quaternion1d(-_HUGE);
				long long divergentBits=0x7FF8000000000010;
				return Quaternion1d((double&)divergentBits);
			}
			catch(std::overflow_error&)
			{
				return Quaternion1d(_HUGE);
			}
		}
		Vect4d  r_r_permutation				(Vect4d const &x)					{return m_one;}
		Comp4d  c_c_permutation				(Comp4d const &x)					{return Comp4d(m_one, m_zero);}
		Quat4d  q_q_permutation				(Quat4d const &x)					{return Quat4d(m_one, m_zero, m_zero, m_zero);}
		Vect4d r_rr_permutation				(Vect4d const &x, Vect4d const &y)	{return Vect4d(permutation(x.lo(), y.lo()), permutation(x.hi(), y.hi()));}
		Comp4d c_cr_permutation				(Comp4d const &x, Vect4d const &y)
		{
			Complex1d lo=permutation(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.lo(), 0));
			Complex1d hi=permutation(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.hi(), 0));
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_cc_permutation				(Comp4d const &x, Comp4d const &y)
		{
			Complex1d lo=permutation(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.r.lo(), y.i.lo()));
			Complex1d hi=permutation(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.r.hi(), y.i.hi()));
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Quat4d q_qq_permutation				(Quat4d const &x, Quat4d const &y)
		{
			Quaternion1d lo=permutation(Quaternion1d(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo()), Quaternion1d(y.r.lo(), y.i.lo(), y.j.lo(), y.k.lo()));
			Quaternion1d hi=permutation(Quaternion1d(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi()), Quaternion1d(y.r.hi(), y.i.hi(), y.j.hi(), y.k.hi()));
			return Quat4d(Vect4d(lo.R_component_1(), hi.R_component_1()), Vect4d(lo.R_component_2(), hi.R_component_2()),
				Vect4d(lo.R_component_3(), hi.R_component_3()), Vect4d(lo.R_component_4(), hi.R_component_4()));
		}
	
		double combination(double x, double y)
		{
			try
			{
				double rx=x, ry=y;
				double r=boost::math::tgamma(rx+1)/(boost::math::tgamma(rx-ry+1)*boost::math::tgamma(ry+1));
				return r;
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return _HUGE;
				if(valueBits==0x8000000000000000)
					return -_HUGE;
				long long divergentBits=0x7FF8000000000010;
				return *(double*)&divergentBits;
			}
			catch(std::overflow_error&)
			{
				return _HUGE;
			}
		}
		Complex1d combination(Complex1d const &x, Complex1d const &y)
		{
			try
			{
				std::complex<double> cx=x, cy=y;
				return tgamma(cx+1.)/(tgamma(cx-cy+1.)*tgamma(cy+1.));
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return _HUGE;
				if(valueBits==0x8000000000000000)
					return -_HUGE;
				long long divergentBits=0x7FF8000000000010;
				return *(double*)&divergentBits;
			}
			catch(std::overflow_error&)
			{
				return _HUGE;
			}
		}
		Quaternion1d combination(Quaternion1d const &x, Quaternion1d const &y)
		{
			try
			{
				boost::math::quaternion<double> qx=x, qy=y;
				return tgamma(qx+1.)/(tgamma(qx-qy+1.)*tgamma(qy+1.));
			}
			catch(std::domain_error&)
			{
				long long valueBits=(long long&)x;
				if(!valueBits)
					return Quaternion1d(_HUGE);
				if(valueBits==0x8000000000000000)
					return Quaternion1d(-_HUGE);
				long long divergentBits=0x7FF8000000000010;
				return Quaternion1d((double&)divergentBits);
			}
			catch(std::overflow_error&)
			{
				return Quaternion1d(_HUGE);
			}
		}
		Vect4d  r_r_combination				(Vect4d const &x)					{return m_one;}
		Comp4d  c_c_combination				(Comp4d const &x)					{return Comp4d(m_one, m_zero);}
		Quat4d  q_q_combination				(Quat4d const &x)					{return Quat4d(m_one, m_zero, m_zero, m_zero);}
		Vect4d r_rr_combination				(Vect4d const &x, Vect4d const &y)	{return Vect4d(combination(x.lo(), y.lo()), combination(x.hi(), y.hi()));}
		Comp4d c_cr_combination				(Comp4d const &x, Vect4d const &y)
		{
			Complex1d lo=combination(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.lo(), 0));
			Complex1d hi=combination(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.hi(), 0));
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_cc_combination				(Comp4d const &x, Comp4d const &y)
		{
			Complex1d lo=combination(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.r.lo(), y.i.lo()));
			Complex1d hi=combination(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.r.hi(), y.i.hi()));
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Quat4d q_qq_combination				(Quat4d const &x, Quat4d const &y)
		{
			Quaternion1d lo=combination(Quaternion1d(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo()), Quaternion1d(y.r.lo(), y.i.lo(), y.j.lo(), y.k.lo()));
			Quaternion1d hi=combination(Quaternion1d(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi()), Quaternion1d(y.r.hi(), y.i.hi(), y.j.hi(), y.k.hi()));
			return Quat4d(Vect4d(lo.R_component_1(), hi.R_component_1()), Vect4d(lo.R_component_2(), hi.R_component_2()),
				Vect4d(lo.R_component_3(), hi.R_component_3()), Vect4d(lo.R_component_4(), hi.R_component_4()));
		}//*/

		__forceinline Quat4d sgnu(Quat4d const &x){return Quat4d(m_zero, x.i, x.j, x.k);}
		__forceinline Quat4d acosh(Quat4d const &x){return log(x+sqrt(sq(x)-m_one));}
		__forceinline Comp4d asinh(Comp4d const &x){return log(x+sqrt(sq(x)+m_one));}
		__forceinline Quat4d asinh(Quat4d const &x){return log(x+sqrt(sq(x)+m_one));}
		__forceinline Vect4d sinhc(Vect4d const &x){return Vect4d(::sinh(x.v))/x;}
		__forceinline Comp4d cos(Comp4d const &x)
		{
			Comp4d exp_ix=exp(m_i*x);
			return m_half*exp_ix+m_half/exp_ix;
		//	return (exp(m_i*x)+exp(-m_i*x))*m_half;
		}
		__forceinline Quat4d cos(Quat4d const &x)//boost::math
		{
			Vect4d z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
			Vec4d sin_xr, cos_xr;
			sin_xr=sincos(&cos_xr, x.r.v);
			Vect4d w=-Vect4d(sin_xr)*sinhc(z);
			return Quat4d(Vect4d(cos_xr)*Vect4d(::cosh(z.v)), w*x.i, w*x.j, w*x.k);
		}
		__forceinline Comp4d sin(Comp4d const &x)
		{
			Comp4d exp_ix=exp(m_i*x);
			return m_half*exp_ix-m_half/exp_ix;
		//	return (exp(m_i*x)-exp(-m_i*x))*(-m_half*m_i);
		}
		__forceinline Quat4d sin(Quat4d const &x)//boost::math
		{
			Vect4d z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
			Vec4d sin_xr, cos_xr;
			sin_xr=sincos(&cos_xr, x.r.v);
			Vect4d w=-Vect4d(cos_xr)*sinhc(z);
			return Quat4d(Vect4d(sin_xr)*Vect4d(::cosh(z.v)), w*x.i, w*x.j, w*x.k);
		}
		void  r_r_cos					(VectP &r, VectP const &x)					{assign(r, Vect4d(::cos(Vect4d(x).v)));VLEAVE;}
		void  c_c_cos					(CompP &r, CompP const &x)					{assign(r, cos(Comp4d(x)));VLEAVE;}
		void  q_q_cos					(QuatP &r, QuatP const &x)					{assign(r, cos(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d acos(Comp4d const &x){return -m_i*log(x+sqrt(sq(x)-m_one));}
		__forceinline Quat4d acos(Quat4d const &x){return -sgnu(x)*acosh(x);}
		void  c_c_acos					(CompP &r, CompP const &x)					{assign(r, acos(Comp4d(x)));VLEAVE;}
		void  q_q_acos					(QuatP &r, QuatP const &x)					{assign(r, acos(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d cosh(Comp4d const &x)
		{
			Comp4d exp_x=exp(x);
			return m_half*exp_x+m_half/exp_x;
		//	return (exp(x)+exp(-x))*m_half;
		}
		__forceinline Quat4d cosh(Quat4d const &x)
		{
			Quat4d exp_x=exp(x);
			return m_half*exp_x+m_half/exp_x;
		//	return (exp(x)+exp(-x))*m_half;
		}
		void  r_r_cosh					(VectP &r, VectP const &x)					{assign(r, Vect4d(::cosh(Vect4d(x).v)));VLEAVE;}
		void  c_c_cosh					(CompP &r, CompP const &x)					{assign(r, cosh(Comp4d(x)));VLEAVE;}
		void  q_q_cosh					(QuatP &r, QuatP const &x)					{assign(r, cosh(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d acosh(Comp4d const &x){return log(x+sqrt(sq(x)-m_one));}
		void  c_c_acosh					(CompP &r, CompP const &x)					{assign(r, acosh(Comp4d(x)));VLEAVE;}
		void  q_q_acosh					(QuatP &r, QuatP const &x)					{assign(r, acosh(Quat4d(x)));VLEAVE;}

		void  r_r_cosc					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, Vect4d(::cos(rx.v))/rx);VLEAVE;}
		void  c_c_cosc					(CompP &r, CompP const &x)					{Comp4d cx=x; assign(r, cos(cx)/cx);VLEAVE;}
		void  q_q_cosc					(QuatP &r, QuatP const &x)					{Quat4d qx=x; assign(r, cos(qx)/qx);VLEAVE;}

		void  r_r_sec					(VectP &r, VectP const &x)					{assign(r, m_one/Vect4d(::cos(Vect4d(x).v)));VLEAVE;}
		void  c_c_sec					(CompP &r, CompP const &x)					{assign(r, m_one/cos(Comp4d(x)));VLEAVE;}
		void  q_q_sec					(QuatP &r, QuatP const &x)					{assign(r, m_one/cos(Quat4d(x)));VLEAVE;}

		void  c_c_asec					(CompP &r, CompP const &x)					{assign(r, acos(m_one/Comp4d(x)));VLEAVE;}
		void  q_q_asec					(QuatP &r, QuatP const &x)					{assign(r, acos(m_one/Quat4d(x)));VLEAVE;}

		void  r_r_sech					(VectP &r, VectP const &x)					{assign(r, m_one/Vect4d(::cosh(Vect4d(x).v)));VLEAVE;}
		void  c_c_sech					(CompP &r, CompP const &x)					{assign(r, m_one/cosh(Comp4d(x)));VLEAVE;}
		void  q_q_sech					(QuatP &r, QuatP const &x)					{assign(r, m_one/cosh(Quat4d(x)));VLEAVE;}

		void  c_c_asech					(CompP &r, CompP const &x)					{assign(r, acosh(m_one/Comp4d(x)));VLEAVE;}
		void  q_q_asech					(QuatP &r, QuatP const &x)					{assign(r, acosh(m_one/Quat4d(x)));VLEAVE;}

		void  r_r_sin					(VectP &r, VectP const &x)					{assign(r, Vect4d(::sin(Vect4d(x).v)));VLEAVE;}
		void  c_c_sin					(CompP &r, CompP const &x)					{assign(r, sin(Comp4d(x)));VLEAVE;}
		void  q_q_sin					(QuatP &r, QuatP const &x)					{assign(r, sin(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d asin(Comp4d const &x){return -m_i*log(m_i*x+sqrt(m_one-sq(x)));}
		__forceinline Quat4d asin(Quat4d const &x){Quat4d t=sgnu(x); return -t*asinh(x*t);}
		void  c_c_asin					(CompP &r, CompP const &x)					{assign(r, asin(Comp4d(x)));VLEAVE;}
		void  q_q_asin					(QuatP &r, QuatP const &x)					{assign(r, asin(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d sinh(Comp4d const &x)
		{
			Comp4d exp_x=exp(x);
			return m_half*exp_x-m_half/exp_x;
		//	return (exp(x)-exp(-x))*m_half;
		}
		__forceinline Quat4d sinh(Quat4d const &x)
		{
			Quat4d exp_x=exp(x);
			return m_half*exp_x-m_half/exp_x;
		//	return (exp(x)-exp(-x))*m_half;
		}
		void  r_r_sinh					(VectP &r, VectP const &x)					{assign(r, Vect4d(::sinh(Vect4d(x).v)));VLEAVE;}
		void  c_c_sinh					(CompP &r, CompP const &x)					{assign(r, sinh(Comp4d(x)));VLEAVE;}
		void  q_q_sinh					(QuatP &r, QuatP const &x)					{assign(r, sinh(Quat4d(x)));VLEAVE;}

		void  r_r_asinh					(VectP &r, VectP const &x)					{assign(r, Vect4d(::asinh(Vect4d(x).v)));VLEAVE;}
		void  c_c_asinh					(CompP &r, CompP const &x)					{assign(r, asinh(Comp4d(x)));VLEAVE;}
		void  q_q_asinh					(QuatP &r, QuatP const &x)					{assign(r, asinh(Quat4d(x)));VLEAVE;}

		void  r_r_sinc					(VectP &r, VectP const &x)					{Vect4d rx=x, mask=rx==m_zero; assign(r, Vect4d(::sin(rx.v))/rx&mask.complement()|m_one&mask);VLEAVE;}
		void  c_c_sinc					(CompP &r, CompP const &x)					{Comp4d cx=x; Vect4d mask=cx.r==m_zero&cx.i==m_zero; assign(r, and(sin(cx)/cx, mask.complement())|m_one&mask);VLEAVE;}
		void  q_q_sinc					(QuatP &r, QuatP const &x)					{Quat4d qx=x; Vect4d mask=qx.r==m_zero&qx.i==m_zero&qx.j==m_zero&qx.k==m_zero; assign(r, and(sin(qx)/qx, mask.complement())|m_one&mask);VLEAVE;}

		void  r_r_sinhc					(VectP &r, VectP const &x)					{Vect4d rx=x, mask=rx==m_zero; assign(r, Vect4d(::sinh(rx.v))/rx&mask.complement()|m_one&mask);VLEAVE;}
		void  c_c_sinhc					(CompP &r, CompP const &x)					{Comp4d cx=x; Vect4d mask=cx.r==m_zero&cx.i==m_zero; assign(r, and(sinh(cx)/cx, mask.complement())|m_one&mask);VLEAVE;}
		void  q_q_sinhc					(QuatP &r, QuatP const &x)					{Quat4d qx=x; Vect4d mask=qx.r==m_zero&qx.i==m_zero&qx.j==m_zero&qx.k==m_zero; assign(r, and(sinh(qx)/qx, mask.complement())|m_one&mask);VLEAVE;}

		void  r_r_csc					(VectP &r, VectP const &x)					{assign(r, m_one/Vect4d(::sin(Vect4d(x).v)));VLEAVE;}
		void  c_c_csc					(CompP &r, CompP const &x)					{assign(r, m_one/sin(Comp4d(x)));VLEAVE;}
		void  q_q_csc					(QuatP &r, QuatP const &x)					{assign(r, m_one/sin(Quat4d(x)));VLEAVE;}

		void  c_c_acsc					(CompP &r, CompP const &x)					{assign(r, asin(m_one/Comp4d(x)));VLEAVE;}
		void  q_q_acsc					(QuatP &r, QuatP const &x)					{assign(r, asin(m_one/Quat4d(x)));VLEAVE;}

		void  r_r_csch					(VectP &r, VectP const &x)					{assign(r, m_one/Vect4d(::sinh(Vect4d(x).v)));VLEAVE;}
		void  c_c_csch					(CompP &r, CompP const &x)					{assign(r, m_one/sinh(Comp4d(x)));VLEAVE;}
		void  q_q_csch					(QuatP &r, QuatP const &x)					{assign(r, m_one/sinh(Quat4d(x)));VLEAVE;}

		void  r_r_acsch					(VectP &r, VectP const &x)					{assign(r, Vect4d(::asinh((m_one/Vect4d(x)).v)));VLEAVE;}
		void  c_c_acsch					(CompP &r, CompP const &x)					{assign(r, asinh(m_one/Comp4d(x)));VLEAVE;}
		void  q_q_acsch					(QuatP &r, QuatP const &x)					{assign(r, asinh(m_one/Quat4d(x)));VLEAVE;}

		__forceinline Comp4d tan(Comp4d const &x)
		{
			const Comp4d m_two_i=m_i+m_i;
			Comp4d exp_2ix=exp(m_two_i*x);
			return (exp_2ix-m_one)/((exp_2ix+m_one)*m_i);
		}
		__forceinline Quat4d tan(Quat4d const &x)
		{
			const Quat4d m_two_i=m_i+m_i;
			Quat4d exp_2ix=exp(m_two_i*x);
			return (exp_2ix-m_one)/((exp_2ix+m_one)*m_i);
		}
		void  r_r_tan					(VectP &r, VectP const &x)					{assign(r, Vect4d(::tan(Vect4d(x).v)));VLEAVE;}
		void  c_c_tan					(CompP &r, CompP const &x)					{assign(r, tan(Comp4d(x)));VLEAVE;}
		void  q_q_tan					(QuatP &r, QuatP const &x)					{assign(r, tan(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d atan(Comp4d const &x){return (m_i*m_half)*log((m_i+x)/(m_i-x));}
		__forceinline Quat4d atan(Quat4d const &x){return (m_i*m_half)*log((m_i+x)/(m_i-x));}
		void  r_r_atan					(VectP &r, VectP const &x)					{assign(r, Vect4d(::atan(Vect4d(x).v)));VLEAVE;}
		void  c_c_atan					(CompP &r, CompP const &x)					{assign(r, atan(Comp4d(x)));VLEAVE;}
		void  q_q_atan					(QuatP &r, QuatP const &x)					{assign(r, atan(Quat4d(x)));VLEAVE;}
		void r_rr_atan					(VectP &r, VectP const &y, VectP const &x)	{assign(r, Vect4d(::atan2(Vect4d(y).v, Vect4d(x).v)));VLEAVE;}
		void c_rc_atan					(CompP &r, VectP const &y, CompP const &x)
		{
			Comp4d cx=x;	Vect4d ry=y;
			Comp4d tr=atan(ry/cx);
			Vect4d mask_y=ry<m_zero, mask_x=cx.r<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void q_rq_atan					(QuatP &r, VectP const &y, QuatP const &x)
		{
			Quat4d qx=x;	Vect4d ry=y;
			Quat4d tr=atan(ry/qx);
			Vect4d mask_y=ry<m_zero, mask_x=qx.r<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void c_cr_atan					(CompP &r, CompP const &y, VectP const &x)
		{
			Vect4d rx=x;	Comp4d cy=y;
			Comp4d tr=atan(cy/rx);
			Vect4d mask_y=cy.r<m_zero, mask_x=rx<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void c_cc_atan					(CompP &r, CompP const &y, CompP const &x)
		{
			Comp4d cx=x;	Comp4d cy=y;
			Comp4d tr=atan(cy/cx);
			Vect4d mask_y=cy.r<m_zero, mask_x=cx.r<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void q_cq_atan					(QuatP &r, CompP const &y, QuatP const &x)
		{
			Quat4d qx=x;	Comp4d cy=y;
			auto tr=atan(cy/qx);
			Vect4d mask_y=cy.r<m_zero, mask_x=qx.r<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void q_qr_atan					(QuatP &r, QuatP const &y, VectP const &x)
		{
			Vect4d rx=x;	Quat4d qy=y;
			auto tr=atan(qy/rx);
			Vect4d mask_y=qy.r<m_zero, mask_x=rx<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void q_qc_atan					(QuatP &r, QuatP const &y, CompP const &x)
		{
			Comp4d cx=x;	Quat4d qy=y;
			auto tr=atan(qy/cx);
			Vect4d mask_y=qy.r<m_zero, mask_x=cx.r<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}
		void q_qq_atan					(QuatP &r, QuatP const &y, QuatP const &x)
		{
			Quat4d qx=x;	Quat4d qy=y;
			auto tr=atan(qy/qx);
			Vect4d mask_y=y.r<m_zero, mask_x=qx.r<m_zero;
			Vect4d addition=m_pi&mask_y;
			Vect4d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
			VLEAVE;
		}

		__forceinline Comp4d tanh(Comp4d const &x){Comp4d e2x=exp(x+x); return (e2x-m_one)/(e2x+m_one);}
		__forceinline Quat4d tanh(Quat4d const &x){Quat4d e2x=exp(x+x); return (e2x-m_one)/(e2x+m_one);}
		void  r_r_tanh					(VectP &r, VectP const &x)					{assign(r, Vect4d(::tanh(Vect4d(x).v)));}
		void  c_c_tanh					(CompP &r, CompP const &x)					{assign(r, tanh(Comp4d(x)));VLEAVE;}
		void  q_q_tanh					(QuatP &r, QuatP const &x)					{assign(r, tanh(Quat4d(x)));VLEAVE;}

		__forceinline Comp4d atanh(Comp4d const &x){return m_half*log((m_one+x)/(m_one-x));}
		__forceinline Quat4d atanh(Quat4d const &x){return m_half*log((m_one+x)/(m_one-x));}
		void  c_c_atanh					(CompP &r, CompP const &x)					{assign(r, atanh(Comp4d(x)));VLEAVE;}
		void  q_q_atanh					(QuatP &r, QuatP const &x)					{assign(r, atanh(Quat4d(x)));VLEAVE;}

		void  r_r_tanc					(VectP &r, VectP const &x)					{Vect4d rx=x, mask=rx==m_zero; assign(r, Vect4d(::tan(rx.v))/rx&mask.complement()|m_one&mask);VLEAVE;}
		void  c_c_tanc					(CompP &r, CompP const &x)					{Comp4d cx=x; Vect4d mask=cx.r==m_zero&cx.i==m_zero; assign(r, and(tan(cx)/cx, mask.complement())|m_one&mask);VLEAVE;}
		void  q_q_tanc					(QuatP &r, QuatP const &x)					{Quat4d qx=x; Vect4d mask=qx.r==m_zero&qx.i==m_zero&qx.j==m_zero&qx.k==m_zero; assign(r, and(tan(qx)/qx, mask.complement())|m_one&mask);VLEAVE;}

		void  r_r_cot					(VectP &r, VectP const &x)					{assign(r, m_one/Vect4d(::tan(Vect4d(x).v)));VLEAVE;}
		void  c_c_cot					(CompP &r, CompP const &x)					{assign(r, m_one/tan(Comp4d(x)));VLEAVE;}
		void  q_q_cot					(QuatP &r, QuatP const &x)					{assign(r, m_one/tan(Quat4d(x)));VLEAVE;}

		void  r_r_acot					(VectP &r, VectP const &x)
		{
			Vect4d rx=x, mask=rx==m_zero;
			assign(r, Vect4d(::atan((m_one/rx).v))&mask.complement()|m_pi_2&mask);
			VLEAVE;
		}
		void  c_c_acot					(CompP &r, CompP const &x)
		{
			Comp4d cx=x;
			Vect4d mask=(cx.r==m_zero)&(cx.i==m_zero);
			assign(r, and(atan(m_one/cx), mask.complement())|m_pi_2&mask);
			VLEAVE;
		}
		void  q_q_acot					(QuatP &r, QuatP const &x)
		{
			Quat4d qx=x;
			Vect4d mask=(qx.r==m_zero)&(qx.i==m_zero)&(qx.j==m_zero)&(qx.k==m_zero);
			assign(r, and(atan(m_one/qx), mask.complement())|m_pi_2&mask);
			VLEAVE;
		}

		void  r_r_coth					(VectP &r, VectP const &x)					{assign(r, m_one/Vect4d(::tanh(Vect4d(x).v)));VLEAVE;}
		void  c_c_coth					(CompP &r, CompP const &x)					{assign(r, m_one/tanh(Comp4d(x)));VLEAVE;}
		void  q_q_coth					(QuatP &r, QuatP const &x)					{assign(r, m_one/tanh(Quat4d(x)));VLEAVE;}

		void  c_c_acoth					(CompP &r, CompP const &x)					{assign(r, atanh(m_one/Comp4d(x)));VLEAVE;}
		void  q_q_acoth					(QuatP &r, QuatP const &x)					{assign(r, atanh(m_one/Quat4d(x)));VLEAVE;}

		void  r_r_exp					(VectP &r, VectP const &x)					{assign(r, Vect4d(::exp(Vect4d(x).v)));VLEAVE;}
		void  c_c_exp					(CompP &r, CompP const &x)					{assign(r, exp(Comp4d(x)));VLEAVE;}
		void  q_q_exp					(QuatP &r, QuatP const &x)					{assign(r, exp(Quat4d(x)));VLEAVE;}
	
		__forceinline Vect4d pow(Vect4d const &x, Vect4d const &y){return Vect4d(::exp((y*Vect4d(::log(x.v))).v));}
		void  r_r_fib					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, (Vect4d(::exp((rx*m_ln_phi).v))-Vect4d(::cos((m_pi*rx).v))*Vect4d(::exp((-rx*m_ln_phi).v)))*m_inv_sqrt5);VLEAVE;}
	//	void  r_r_fib						(Vect4d const &x)					{return (pow(m_phi, x)-Vect4d(::cos((m_pi*x).v))*pow(m_phi, -x))/m_sqrt5;}
	//	void  r_r_fib						(Vect4d const &x)					{return (Vect4d(::pow(m_phi.v, x.v))-Vect4d(::cos((m_pi*x).v))*Vect4d(::pow(m_phi.v, (-x).v)))/m_sqrt5;}
		void  c_c_fib					(CompP &r, CompP const &x)					{Comp4d cx=x; assign(r, (exp(cx*m_ln_phi)-cos(m_pi*cx)*exp(-cx*m_ln_phi))*m_inv_sqrt5);VLEAVE;}
	//	void  c_c_fib					(CompP &r, CompP const &x)					{Comp4d cx=x; assign(r, ((m_phi^cx)-cos(m_pi*cx)*(m_phi^-cx))/m_sqrt5;}
		void  q_q_fib					(QuatP &r, QuatP const &x)					{Quat4d qx=x; assign(r, (exp(qx*m_ln_phi)-cos(m_pi*qx)*exp(-qx*m_ln_phi))*m_inv_sqrt5);VLEAVE;}
	
	/*	const double rand_norm=9.31322574615479e-010;
		Vect4d my_rand(){return Vect4d((rand()<<15|rand())*rand_norm, (rand()<<15|rand())*rand_norm);}
		Vect4d  r_r_random					(Vect4d const &x)					{return my_rand();}
		Comp4d  c_c_random					(Comp4d const &x)					{return Comp4d(my_rand(), my_rand());}
		Quat4d  q_q_random					(Quat4d const &x)					{return Quat4d(my_rand(), my_rand(), my_rand(), my_rand());}
		Vect4d r_rr_random					(Vect4d const &x, Vect4d const &y)	{return my_rand();}
		Comp4d c_cr_random					(Comp4d const &x, Vect4d const &y)	{return Comp4d(my_rand(), my_rand());}
		Comp4d c_cc_random					(Comp4d const &x, Comp4d const &y)	{return Comp4d(my_rand(), my_rand());}
		Quat4d q_qq_random					(Quat4d const &x, Quat4d const &y)	{return Quat4d(my_rand(), my_rand(), my_rand(), my_rand());}

		__forceinline double beta(double x, double y)
		{
			try
			{
				return boost::math::beta(x, y);
			}
			catch(...){return _qnan;}
		}
		Vect4d  r_r_beta					(Vect4d const &x)					{return Vect4d(beta(x.lo(), x.lo()), beta(x.hi(), x.hi()));}
		Vect4d r_rr_beta					(Vect4d const &x, Vect4d const &y)	{return Vect4d(beta(x.lo(), y.lo()), beta(x.hi(), y.hi()));}

		__forceinline double cyl_bessel_j(double x, double y)
		{
			try
			{
				return boost::math::cyl_bessel_j(x, y);
			}
			catch(...){return _qnan;}
		}
		Vect4d  r_r_cyl_bessel_j			(Vect4d const &x)					{return Vect4d(cyl_bessel_j(0, x.lo()), cyl_bessel_j(0, x.hi()));}
		Vect4d r_rr_cyl_bessel_j			(Vect4d const &x, Vect4d const &y)	{return Vect4d(cyl_bessel_j(x.lo(), y.lo()), cyl_bessel_j(x.hi(), y.hi()));}

		__forceinline double cyl_neumann(double x, double y)
		{
			if(x>0)
				int LOL_1=0;
			if(y!=y||y<0)
				return _qnan;
			try
			{
				return boost::math::cyl_neumann(x, y);
			}
			catch(std::domain_error&){return _qnan;}
			catch(std::overflow_error&){return -_HUGE;}
			catch(...){return _qnan;}
		}
		Vect4d  r_r_cyl_neumann				(Vect4d const &x)					{return Vect4d(cyl_neumann(0, x.lo()), cyl_neumann(0, x.hi()));}
		Vect4d r_rr_cyl_neumann				(Vect4d const &x, Vect4d const &y)	{return Vect4d(cyl_neumann(x.lo(), y.lo()), cyl_neumann(x.hi(), y.hi()));}

		__forceinline Complex1d r_hankel1(double x, double y)
		{
			try
			{
				return boost::math::cyl_bessel_j(x, y)+Complex1d(0, 1)*boost::math::cyl_neumann(x, y);
			}
			catch(std::domain_error&){return _qnan;}
			catch(std::overflow_error&){return -_HUGE;}
			catch(...){return _qnan;}
		}
		Comp4d  c_r_hankel1					(Vect4d const &x)
		{
			Complex1d lo=r_hankel1(x.lo(), x.lo()), hi=r_hankel1(x.hi(), x.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d  c_c_hankel1					(Comp4d const &x)
		{
			Complex1d lo=r_hankel1(x.r.lo(), x.i.lo()), hi=r_hankel1(x.r.hi(), x.i.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}
		Comp4d c_rr_hankel1					(Vect4d const &x, Vect4d const &y)
		{
			Complex1d lo=r_hankel1(x.lo(), y.lo()), hi=r_hankel1(x.hi(), y.hi());
			return Comp4d(Vect4d(lo.real(), hi.real()), Vect4d(lo.imag(), hi.imag()));
		}//*/

		__forceinline Vect4d sgn(Vect4d const &x){return ((x>m_zero)&m_one)-((x<m_zero)&m_one);}
		__forceinline Comp4d sgn(Comp4d const &x){Vect4d mag=x.abs(); return and(x/mag, mag.r_is_true());}
		__forceinline Quat4d sgn(Quat4d const &x){Vect4d mag=x.abs(); return and(x/mag, mag.r_is_true());}
		__forceinline Vect4d step(Vect4d const &x){return m_half+m_half*sgn(x);}
		__forceinline Comp4d step(Comp4d const &x){return m_half+m_half*sgn(x);}
		__forceinline Quat4d step(Quat4d const &x){return m_half+m_half*sgn(x);}
		void  r_r_step					(VectP &r, VectP const &x)					{assign(r, step(Vect4d(x)));VLEAVE;}
		void  c_c_step					(CompP &r, CompP const &x)					{assign(r, step(Comp4d(x)));VLEAVE;}
		void  q_q_step					(QuatP &r, QuatP const &x)					{assign(r, step(Quat4d(x)));VLEAVE;}

		void  r_r_rect					(VectP &r, VectP const &x)					{Vect4d rx=x; assign(r, step(rx+m_half)-step(rx-m_half));VLEAVE;}
		void  c_c_rect					(CompP &r, CompP const &x)					{Comp4d cx=x; assign(r, step(cx+m_half)-step(cx-m_half));VLEAVE;}
		void  q_q_rect					(QuatP &r, QuatP const &x)					{Quat4d qx=x; assign(r, step(qx+m_half)-step(qx-m_half));VLEAVE;}

		void  r_r_trgl					(VectP &r, VectP const &x)					{Vect4d t=Vect4d(x).abs(); assign(r, (m_one-t)&t<m_one);VLEAVE;}
		void  r_c_trgl					(VectP &r, CompP const &x)					{Vect4d t=Comp4d(x).abs(); assign(r, and(m_one-t, t<m_one));VLEAVE;}
		void  r_q_trgl					(VectP &r, QuatP const &x)					{Vect4d t=Quat4d(x).abs(); assign(r, and(m_one-t, t<m_one));VLEAVE;}

		void  r_r_sqwv					(VectP &r, VectP const &x)					{Vect4d rx=x;	assign(r, rx-rx.floor()<m_half&m_one);VLEAVE;}
		void  r_c_sqwv					(VectP &r, CompP const &x)					{Vect4d rx=x.r; assign(r, rx-rx.floor()<m_half&m_one);VLEAVE;}
		void  r_q_sqwv					(VectP &r, QuatP const &x)					{Vect4d rx=x.r; assign(r, rx-rx.floor()<m_half&m_one);VLEAVE;}
		void r_rr_sqwv					(VectP &r, VectP const &x, VectP const &y)	{Vect4d rx=x;	assign(r, rx-rx.floor()<Vect4d(y)&m_one);VLEAVE;}
		void r_rc_sqwv					(VectP &r, VectP const &x, CompP const &y)	{Vect4d rx=x;	assign(r, rx-rx.floor()<Vect4d(y.r)&m_one);VLEAVE;}
		void r_rq_sqwv					(VectP &r, VectP const &x, QuatP const &y)	{Vect4d rx=x;	assign(r, rx-rx.floor()<Vect4d(y.r)&m_one);VLEAVE;}
		void r_cr_sqwv					(VectP &r, CompP const &x, VectP const &y)	{Vect4d rx=x.r;	assign(r, rx-rx.floor()<Vect4d(y)&m_one);VLEAVE;}
		void r_cc_sqwv					(VectP &r, CompP const &x, CompP const &y)	{Vect4d rx=x.r;	assign(r, rx-rx.floor()<Vect4d(y.r)&m_one);VLEAVE;}
		void r_cq_sqwv					(VectP &r, CompP const &x, QuatP const &y)	{Vect4d rx=x.r;	assign(r, rx-rx.floor()<Vect4d(y.r)&m_one);VLEAVE;}
		void r_qr_sqwv					(VectP &r, QuatP const &x, VectP const &y)	{Vect4d rx=x.r;	assign(r, rx-rx.floor()<Vect4d(y)&m_one);VLEAVE;}
		void r_qc_sqwv					(VectP &r, QuatP const &x, CompP const &y)	{Vect4d rx=x.r;	assign(r, rx-rx.floor()<Vect4d(y.r)&m_one);VLEAVE;}
		void r_qq_sqwv					(VectP &r, QuatP const &x, QuatP const &y)	{Vect4d rx=x.r;	assign(r, rx-rx.floor()<Vect4d(y.r)&m_one);VLEAVE;}

		void  r_r_trwv					(VectP &r, VectP const &x)					{Vect4d rx=x; Vect4d t=(rx-rx.floor()-m_half).abs(); assign(r, t+t);VLEAVE;}
		void  r_c_trwv					(VectP &r, CompP const &x)					{Comp4d cx=x; Vect4d t=(cx-cx.floor()-m_half).abs(); assign(r, t+t);VLEAVE;}
		void  r_q_trwv					(VectP &r, QuatP const &x)					{Quat4d qx=x; Vect4d t=(qx-qx.floor()-m_half).abs(); assign(r, t+t);VLEAVE;}
		void r_rr_trwv					(VectP &r, VectP const &x, VectP const &y)
		{
			Vect4d rx=x;	Vect4d ry=y;
			Vect4d t=rx-rx.floor(), t2=m_one-rx;
			t2=t2-t2.floor();
			Vect4d d=ry&ry>m_zero;
			Vect4d mask=d>m_one;
			d=d&mask.complement()|m_one&mask;
			Vect4d d2=m_one-d;
			Vect4d t_d=t/d, t2_d2=t2/d2;
			assign(r, ((t_d<m_one)&t_d)+((t2_d2<m_one)&t2_d2));
			VLEAVE;
		}
		void c_cr_trwv					(CompP &r, CompP const &x, VectP const &y)
		{
			Comp4d cx=x;	Vect4d ry=y;
			auto t=cx-cx.floor(), t2=m_one-cx;
			t2=t2-t2.floor();
			auto d=and(ry, ry>m_zero);
			auto mask=d>m_one;
			d=and(d, mask.complement())|m_one&mask;
			auto d2=m_one-d;
			auto t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
			VLEAVE;
		}
		void c_cc_trwv					(CompP &r, CompP const &x, CompP const &y)
		{
			Comp4d cx=x;	Comp4d cy=y;
			Comp4d t=cx-cx.floor(), t2=m_one-cx;
			t2=t2-t2.floor();
			Comp4d d=and(cy, cy.r>m_zero);
			Vect4d mask=d.r>m_one;
			d=and(d, mask.complement())|m_one&mask;
			Comp4d d2=m_one-d;
			Comp4d t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
			VLEAVE;
		}
		void q_qq_trwv					(QuatP &r, QuatP const &x, QuatP const &y)
		{
			Quat4d qx=x;	Quat4d qy=y;
			Quat4d t=qx-qx.floor(), t2=m_one-qx;
			t2=t2-t2.floor();
			Quat4d d=and(qy, qy.r>m_zero);
			Vect4d mask=d.r>m_one;
			d=and(d, mask.complement())|m_one&mask;
			Quat4d d2=m_one-d;
			Quat4d t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
			VLEAVE;
		}

		void  r_r_saw					(VectP &r, VectP const &x)
		{
			Vect4d rx=x;
			Vect4d t=rx-rx.floor(), t2=(m_one-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t));
			VLEAVE;
		}
		void  c_c_saw					(CompP &r, CompP const &x)
		{
			Comp4d cx=x;
			auto t=cx-cx.floor(), t2=(m_one-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t));
			VLEAVE;
		}
		void  q_q_saw					(QuatP &r, QuatP const &x)
		{
			Quat4d qx=x;
			auto t=qx-qx.floor(), t2=(m_one-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t));
			VLEAVE;
		}
		void r_rr_saw					(VectP &r, VectP const &x, VectP const &y)
		{
			Vect4d rx=x;	Vect4d ry=y;
			auto t=rx-rx.floor(), t2=(ry-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
			VLEAVE;
		}
		void c_rc_saw					(CompP &r, VectP const &x, CompP const &y)
		{
			Vect4d rx=x;	Comp4d cy=y;
			auto t=rx-rx.floor();
			auto t2=(cy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
			VLEAVE;
		}
		void q_rq_saw					(QuatP &r, VectP const &x, QuatP const &y)
		{
			Vect4d rx=x;	Quat4d qy=y;
			auto t=rx-rx.floor();
			auto t2=(qy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
			VLEAVE;
		}
		void c_cr_saw					(CompP &r, CompP const &x, VectP const &y)
		{
			Comp4d cx=x;	Vect4d ry=y;
			auto t=cx-cx.floor();
			auto t2=(ry-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
			VLEAVE;
		}
		void c_cc_saw					(CompP &r, CompP const &x, CompP const &y)
		{
			Comp4d cx=x, cy=y;
			auto t=cx-cx.floor();
			auto t2=(cy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
			VLEAVE;
		}
		void q_cq_saw					(QuatP &r, CompP const &x, QuatP const &y)
		{
			Comp4d cx=x;	Quat4d qy=y;
			auto t=cx-cx.floor();
			auto t2=(qy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
			VLEAVE;
		}
		void q_qr_saw					(QuatP &r, QuatP const &x, VectP const &y)
		{
			Quat4d qx=x;	Vect4d ry=y;
			auto t=qx-qx.floor();
			auto t2=(ry-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
			VLEAVE;
		}
		void q_qc_saw					(QuatP &r, QuatP const &x, CompP const &y)
		{
			Quat4d qx=x;	Comp4d cy=y;
			auto t=qx-qx.floor();
			auto t2=(cy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
			VLEAVE;
		}
		void q_qq_saw					(QuatP &r, QuatP const &x, QuatP const &y)
		{
			Quat4d qx=x, qy=y;
			auto t=qx-qx.floor();
			auto t2=(qy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
			VLEAVE;
		}

		void r_rr_hypot					(VectP &r, VectP const &x, VectP const &y)	{Vect4d rx=x, ry=y; assign(r, sqrt(rx*rx+ry*ry));VLEAVE;}
	//	void c_cc_hypot					(CompP &r, CompP const &x, CompP const &y)	{Comp4d cx=x; assign(r, sqrt(sq(x)+sq(y));}
	//	void q_qq_hypot					(QuatP &r, QuatP const &x, QuatP const &y)	{Quat4d qx=x; assign(r, sqrt(sq(x)+sq(y));}
		
		__forceinline Vect4d mandelbrot(Vect4d const &x, Vect4d const &y)
		{
			Vect4d rez=0., sq_rez=0., m_sq_limit=16;
			Vect4d k=0., active=m_one;
			for(;;k+=active)
			{
				//condition
				active&=sq_rez<m_sq_limit;
				Vect4d stop=k>=y|active.r_is_false();
				if(stop.v0()&&stop.v1()&&stop.v2()&&stop.v3())
					break;
				//iteration
				rez=sq_rez;//calculate sq(z)

				rez+=x;//add x

				sq_rez=rez*rez;
			}
			return k;
		}
		__forceinline Vect4d mandelbrot(Comp4d const &x, Vect4d const &y)
		{
			Vect4d rez=0., imz=0., sq_rez=0., sq_imz=0., m_sq_limit=16;
			Vect4d k=0., active=m_one;
			for(;;k+=active)
			{
				//condition
				active&=sq_rez+sq_imz<m_sq_limit;
				Vect4d cont=k<y&active.r_is_true();
				if(!_mm256_movemask_ps(_mm256_castpd_ps(cont)))
					break;
				//Vect4d stop=k>=y|active.r_is_false();
				//if(stop.v0()&&stop.v1()&&stop.v2()&&stop.v3())
				//	break;

				//iteration
				imz*=rez;//calculate sq(z)
				imz+=imz;
				rez=sq_rez-sq_imz;

				rez+=x.r, imz+=x.i;//add x

				sq_rez=rez*rez, sq_imz=imz*imz;
			}
			return k;
		}
		void r_r_mandelbrot				(VectP &r, VectP const &x)					{assign(r, mandelbrot(Vect4d(x), Vect4d(200)));}
		void r_c_mandelbrot				(VectP &r, CompP const &x)					{assign(r, mandelbrot(Comp4d(x), Vect4d(200)));}
		void r_rr_mandelbrot			(VectP &r, VectP const &x, VectP const &y)	{assign(r, mandelbrot(x, y));}
		void r_cr_mandelbrot			(VectP &r, CompP const &x, VectP const &y)	{assign(r, mandelbrot(Comp4d(x), Vect4d(y)));}

		void r_rr_min					(VectP &r, VectP const &x, VectP const &y)	{Vect4d rx=x, ry=y;			assign(r, (rx+ry-(rx-ry).abs())*m_half);VLEAVE;}
		void c_cr_min					(CompP &r, CompP const &x, VectP const &y)	{Comp4d cx=x; Vect4d ry=y;	assign(r, (cx+ry-(cx-ry).abs())*m_half);VLEAVE;}
		void c_cc_min					(CompP &r, CompP const &x, CompP const &y)	{Comp4d cx=x; Comp4d cy=y;	assign(r, (cx+cy-(cx-cy).abs())*m_half);VLEAVE;}
		void q_qq_min					(QuatP &r, QuatP const &x, QuatP const &y)	{Quat4d qx=x, qy=y;			assign(r, (qx+qy-(qx-qy).abs())*m_half);VLEAVE;}

		void r_rr_max					(VectP &r, VectP const &x, VectP const &y)	{Vect4d rx=x, ry=y;			assign(r, (rx+ry+(rx-ry).abs())*m_half);VLEAVE;}
		void c_cr_max					(CompP &r, CompP const &x, VectP const &y)	{Comp4d cx=x; Vect4d ry=y;	assign(r, (cx+ry+(cx-ry).abs())*m_half);VLEAVE;}
		void c_cc_max					(CompP &r, CompP const &x, CompP const &y)	{Comp4d cx=x; Comp4d cy=y;	assign(r, (cx+cy+(cx-cy).abs())*m_half);VLEAVE;}
		void q_qq_max					(QuatP &r, QuatP const &x, QuatP const &y)	{Quat4d qx=x, qy=y;			assign(r, (qx+qy+(qx-qy).abs())*m_half);VLEAVE;}

		void r_rr_conditional_110		(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(y)&Vect4d(x)!=m_zero);VLEAVE;}
		void c_rc_conditional_110		(CompP &r, VectP const &x, CompP const &y)	{assign(r, and(Comp4d(y), Vect4d(x)!=m_zero));VLEAVE;}
		void q_rq_conditional_110		(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, and(Quat4d(y), Vect4d(x)!=m_zero));VLEAVE;}
		void r_cr_conditional_110		(VectP &r, CompP const &x, VectP const &y)	{assign(r, and(Vect4d(y), Comp4d(x)!=m_zero));VLEAVE;}
		void c_cc_conditional_110		(CompP &r, CompP const &x, CompP const &y)	{assign(r, and(Comp4d(y), Comp4d(x)!=m_zero));VLEAVE;}
		void q_cq_conditional_110		(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, and(Quat4d(y), Comp4d(x)!=m_zero));VLEAVE;}
		void r_qr_conditional_110		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, and(Vect4d(y), Quat4d(x)!=m_zero));VLEAVE;}
		void c_qc_conditional_110		(CompP &r, QuatP const &x, CompP const &y)	{assign(r, and(Comp4d(y), Quat4d(x)!=m_zero));VLEAVE;}
		void q_qq_conditional_110		(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, and(Quat4d(y), Quat4d(x)!=m_zero));VLEAVE;}
	
		void r_rr_conditional_101		(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect4d(y)&Vect4d(x)==m_zero);VLEAVE;}
		void c_rc_conditional_101		(CompP &r, VectP const &x, CompP const &y)	{assign(r, and(Comp4d(y), Vect4d(x)==m_zero));VLEAVE;}
		void q_rq_conditional_101		(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, and(Quat4d(y), Vect4d(x)==m_zero));VLEAVE;}
		void r_cr_conditional_101		(VectP &r, CompP const &x, VectP const &y)	{assign(r, and(Vect4d(y), Comp4d(x)==m_zero));VLEAVE;}
		void c_cc_conditional_101		(CompP &r, CompP const &x, CompP const &y)	{assign(r, and(Comp4d(y), Comp4d(x)==m_zero));VLEAVE;}
		void q_cq_conditional_101		(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, and(Quat4d(y), Comp4d(x)==m_zero));VLEAVE;}
		void r_qr_conditional_101		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, and(Vect4d(y), Quat4d(x)==m_zero));VLEAVE;}
		void c_qc_conditional_101		(CompP &r, QuatP const &x, CompP const &y)	{assign(r, and(Comp4d(y), Quat4d(x)==m_zero));VLEAVE;}
		void q_qq_conditional_101		(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, and(Quat4d(y), Quat4d(x)==m_zero));VLEAVE;}

		void conditional_111			(QuatP &res, QuatP const &op1, QuatP const &op2, QuatP const &op3, int idx, int op1_ms, int op_ms)
		{
			Vect4d second=
				op1_ms=='R'?	Vect4d(op1.r+idx).r_is_false()
				:op1_ms=='c'?	Comp4d(op1.r+idx, op1.i+idx).c_is_false()
				:				Quat4d(op1.r+idx, op1.i+idx, op1.j+idx, op1.k+idx).q_is_false();
			switch(op_ms)
			{
			case 0:assign(VectP(res.r+idx),											Vect4d(op2.r+idx)&second.complement()									|	Vect4d(op3.r+idx)&second);break;//r_rr
			case 1:assign(CompP(res.r+idx, res.i+idx),							and(Comp4d(op2.r+idx, op2.i+idx), second.complement())						|	Vect4d(op3.r+idx)&second);break;//c_rc
			case 2:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Quat4d(op2.r+idx, op2.i+idx, op2.j+idx, op2.k+idx), second.complement())|	Vect4d(op3.r+idx)&second);break;//q_rq
			case 3:assign(CompP(res.r+idx, res.i+idx),								Vect4d(op2.r+idx)&second.complement()									|and(Comp4d(op3.r+idx, op3.i+idx), second));break;//c_cr
			case 4:assign(CompP(res.r+idx, res.i+idx),							and(Comp4d(op2.r+idx, op2.i+idx), second.complement())						|and(Comp4d(op3.r+idx, op3.i+idx), second));break;//c_cc
			case 5:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Quat4d(op2.r+idx, op2.i+idx, op2.j+idx, op2.k+idx), second.complement())|and(Comp4d(op3.r+idx, op3.i+idx), second));break;//q_cq
			case 6:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),		Vect4d(op2.r+idx)&second.complement()									|and(Quat4d(op3.r+idx, op3.i+idx, op3.j+idx, op3.k+idx), second));break;//q_qr
			case 7:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Comp4d(op2.r+idx, op2.i+idx), second.complement())						|and(Quat4d(op3.r+idx, op3.i+idx, op3.j+idx, op3.k+idx), second));break;//q_qc
			case 8:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Quat4d(op2.r+idx, op2.i+idx, op2.j+idx, op2.k+idx), second.complement())|and(Quat4d(op3.r+idx, op3.i+idx, op3.j+idx, op3.k+idx), second));break;//q_qq
			}
			VLEAVE;
		}

		void  r_r_increment				(VectP &r, VectP const &x)					{assign(r, Vect4d(x)+m_one);VLEAVE;}
		void  c_c_increment				(CompP &r, CompP const &x)					{assign(r, Comp4d(x)+m_one);VLEAVE;}
		void  q_q_increment				(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x)+m_one);VLEAVE;}

		void  r_r_decrement				(VectP &r, VectP const &x)					{assign(r, Vect4d(x)-m_one);VLEAVE;}
		void  c_c_decrement				(CompP &r, CompP const &x)					{assign(r, Comp4d(x)-m_one);VLEAVE;}
		void  q_q_decrement				(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x)-m_one);VLEAVE;}

		void  r_r_assign				(VectP &r, VectP const &x)					{assign(r, Vect4d(x));VLEAVE;}
		void  c_c_assign				(CompP &r, CompP const &x)					{assign(r, Comp4d(x));VLEAVE;}
		void  q_q_assign				(QuatP &r, QuatP const &x)					{assign(r, Quat4d(x));VLEAVE;}
	}//namespace sse2
}//namespace G2
namespace	modes
{
	void colorFunction_bcw_avx(CompP const &_v, int *rgb)
	{
		Comp4d v=_v;

	/*	for(int k=0;k<4;++k)
		{
			double r=v.r.get(k), i=v.i.get(k);
			double const cos_pi_6=0.866025403784439, sin_pi_6=0.5, threshold=10, inv_th=1./threshold;
			double hyp=sqrt(r*r+i*i), costh=r/hyp, sinth=i/hyp,
				mag=255*exp(-hyp*G2::_ln2*inv_th);
			double
				red		=1+costh*cos_pi_6-sinth*sin_pi_6,
				green	=1+sinth,
				blue	=1+costh*-cos_pi_6-sinth*sin_pi_6;
			if(hyp<threshold)//small
				red*=mag
		}//*/

		Vect4d threshold=10;
		double inv_th=0.1;//1			//black->color->white			49.19 cycles/px
		Vect4d const cos_pi_6=0.866025403784439, sin_pi_6=0.5, _FF=255, _7F=0x7F,
			ln2_th=G2::_ln2*inv_th;
		Vect4d hyp=v.abs(), costh=v.r/hyp, sinth=v.i/hyp,
			mag=_FF*exp(-hyp*ln2_th);
		Vect4d
			red		=m_one+costh*cos_pi_6-sinth*sin_pi_6,
			green	=m_one+sinth,
			blue	=m_one+costh*-cos_pi_6-sinth*sin_pi_6;
		Vect4d small=hyp<threshold;//, small_c=small.complement();
		mag=_mm256_blendv_pd(mag, _FF-mag, small);
		Vect4d magx2=mag+mag;
		red*=mag, green*=mag, blue*=mag;
		red		=_mm256_blendv_pd(_FF-magx2+red		, red	, small);
		green	=_mm256_blendv_pd(_FF-magx2+green	, green	, small);
		blue	=_mm256_blendv_pd(_FF-magx2+blue	, blue	, small);
		//red	=_mm256_blendv_pd(_FF-mag*(m_two-red	), red	*mag, small);
		//green	=_mm256_blendv_pd(_FF-mag*(m_two-green	), green*mag, small);
		//blue	=_mm256_blendv_pd(_FF-mag*(m_two-blue	), blue	*mag, small);
		//mag=_FF-mag&small|mag&small_c;
		//red	=red	*mag&small|_FF-mag*(m_two-red	)&small_c;
		//green	=green	*mag&small|_FF-mag*(m_two-green	)&small_c;
		//blue	=blue	*mag&small|_FF-mag*(m_two-blue	)&small_c;
		Vect4d
			nan=v.r!=v.r|v.i!=v.i,// nan_c=nan.complement(),
			inf=v.r==m_inf|v.i==m_inf;//, inf_c=inf.complement();
		red		=_mm256_blendv_pd(red,		_FF, inf);
		green	=_mm256_blendv_pd(green,	_FF, inf);
		blue	=_mm256_blendv_pd(blue,		_FF, inf);
		red		=_mm256_blendv_pd(red,		_7F, nan);
		green	=_mm256_blendv_pd(green,	_7F, nan);
		blue	=_mm256_blendv_pd(blue,		_7F, nan);
		//red	=_FF&inf|red	&inf_c;
		//green	=_FF&inf|green	&inf_c;
		//blue	=_FF&inf|blue	&inf_c;
		//red	=_7F&nan|red	&nan_c;
		//green	=_7F&nan|green	&nan_c;
		//blue	=_7F&nan|blue	&nan_c;
		//red=_mm256_cvtpd_epi32(red);
		//green=_mm256_cvtpd_epi32(green);
		//blue=_mm256_cvtpd_epi32(blue);
		auto p=(unsigned char*)rgb;
		p[ 0]=unsigned char(blue.v0()), p[ 1]=unsigned char(green.v0()), p[ 2]=unsigned char(red.v0()), p[ 3]=0;
		p[ 4]=unsigned char(blue.v1()), p[ 5]=unsigned char(green.v1()), p[ 6]=unsigned char(red.v1()), p[ 7]=0;
		p[ 8]=unsigned char(blue.v2()), p[ 9]=unsigned char(green.v2()), p[10]=unsigned char(red.v2()), p[11]=0;
		p[12]=unsigned char(blue.v3()), p[13]=unsigned char(green.v3()), p[14]=unsigned char(red.v3()), p[15]=0;//*/
		VLEAVE;
	}
	void colorFunction_bc_l_avx(CompP const &_v, int *rgb)
//	void colorFunction_avx(VectP const &pr, VectP const &pi, int *c)
//	void colorFunction_avx(VectP const &pr, VectP const &pi, int &c1, int &c2, int &c3, int &c4)
	{
		Comp4d v=_v;
		Vect4d const cos_pi_6=0.866025403784439, sin_pi_6=0.5, _FF=255, _7F=0x7F,//black->color whitening loop		26.52 cycles/px
			one_ten=0.1, one_perc=0.01, ten=10, twelve=12.75, inv255=1./255;
		Vect4d hyp=v.abs(), costh=v.r/hyp, sinth=v.i/hyp,
			fh=(hyp*one_ten).floor(), c, f;
		c=fh*one_perc;
		c=_FF*(c-c.floor()), f=twelve*(hyp-fh*ten)*(m_one-c*inv255);
		Vect4d
			red		=m_one+costh*cos_pi_6-sinth*sin_pi_6,
			green	=m_one+sinth,
			blue	=m_one+costh*-cos_pi_6-sinth*sin_pi_6;
		red=c+f*red, green=c+f*green, blue=c+f*blue;
		Vect4d
			nan=v.r!=v.r|v.i!=v.i,// nan_c=nan.complement(),
			inf=v.r==m_inf|v.i==m_inf;//, inf_c=inf.complement();
		red		=_mm256_blendv_pd(red,		_FF, inf);
		green	=_mm256_blendv_pd(green,	_FF, inf);
		blue	=_mm256_blendv_pd(blue,		_FF, inf);
		red		=_mm256_blendv_pd(red,		_7F, nan);
		green	=_mm256_blendv_pd(green,	_7F, nan);
		blue	=_mm256_blendv_pd(blue,		_7F, nan);
		//red	=_FF&inf|red	&inf_c;
		//green	=_FF&inf|green	&inf_c;
		//blue	=_FF&inf|blue	&inf_c;
		//red	=_7F&nan|red	&nan_c;
		//green	=_7F&nan|green	&nan_c;
		//blue	=_7F&nan|blue	&nan_c;
		auto p=(unsigned char*)rgb;
		p[ 0]=unsigned char(blue.v0()), p[ 1]=unsigned char(green.v0()), p[ 2]=unsigned char(red.v0()), p[ 3]=0;
		p[ 4]=unsigned char(blue.v1()), p[ 5]=unsigned char(green.v1()), p[ 6]=unsigned char(red.v1()), p[ 7]=0;
		p[ 8]=unsigned char(blue.v2()), p[ 9]=unsigned char(green.v2()), p[10]=unsigned char(red.v2()), p[11]=0;
		p[12]=unsigned char(blue.v3()), p[13]=unsigned char(green.v3()), p[14]=unsigned char(red.v3()), p[15]=0;
		VLEAVE;

	//	memset(c, -1, 4*sizeof(int));
	/*	c[0]=c[1]=c[2]=c[3]=-1;
		double
			r1=r.v.m256d_f64[0], i1=i.v.m256d_f64[0],
			r2=r.v.m256d_f64[1], i2=i.v.m256d_f64[1],
			r3=r.v.m256d_f64[2], i3=i.v.m256d_f64[2],
			r4=r.v.m256d_f64[3], i4=i.v.m256d_f64[3];
		colorFunction_cases(r1, i1, c[0]);
		colorFunction_cases(r2, i2, c[1]);
		colorFunction_cases(r3, i3, c[2]);
		colorFunction_cases(r4, i4, c[3]);
		if((c[0]==-1)+(c[1]==-1)+(c[2]==-1)+(c[3]==-1)==4)
		{
			//const Vect4d m_inf=_mm256_set1_pd(_HUGE), m_minf=_mm256_set1_pd(-_HUGE);
			//const __m128i
			//	nan_color=_mm256_set1_epi32(0x007F7F7F), cinf_color=_mm256_set1_epi32(0x00FFFFFF),
			//	pinf_color=_mm256_set1_epi32(0x00ED7F11), ninf_color=_mm256_set1_epi32(0x00117FED),
			//	iinf_color=_mm256_set1_epi32(0x003FFF3F), niinf_color=_mm256_set1_epi32(0x00BF00BF);
			//Vect4d ri_nan=r!=r|i!=i;
			//Vect4d
			//	r_inf=r==m_inf, r_minf=r==m_minf,
			//	i_inf=i==m_inf, i_minf=i==m_minf, c_inf=(r_inf|r_minf)&(i_inf|i_minf), ri_zero=r==m_zero&i==m_zero;

			const Vect4d m_amp=255/G2::_pi, m_cospi6=0.866025403784439, m_sinpi6=0.5;
			Vect4d hyp=sqrt(r*r+i*i), mag=Vect4d(atan(hyp.v))*m_amp, _1_hyp=m_one/hyp,
				cos_x=r*_1_hyp, sin_x=i*_1_hyp;
			Vect4d cosx_cospi6=cos_x*m_cospi6, sinx_sinpi6=sin_x*m_sinpi6;
			Vect4d m_r=mag*(m_one+cosx_cospi6-sinx_sinpi6), m_g=mag*(m_one+sin_x), m_b=mag*(m_one-cosx_cospi6-sinx_sinpi6);
			auto p=(unsigned char*)c;
			p[0]=unsigned char(m_b.v.m256d_f64[0]), p[1]=unsigned char(m_g.v.m256d_f64[0]), p[2]=unsigned char(m_r.v.m256d_f64[0]);//argb
			p=(unsigned char*)(c+1);
			p[0]=unsigned char(m_b.v.m256d_f64[1]), p[1]=unsigned char(m_g.v.m256d_f64[1]), p[2]=unsigned char(m_r.v.m256d_f64[1]);
			p=(unsigned char*)(c+2);
			p[0]=unsigned char(m_b.v.m256d_f64[2]), p[1]=unsigned char(m_g.v.m256d_f64[2]), p[2]=unsigned char(m_r.v.m256d_f64[2]);
			p=(unsigned char*)(c+3);
			p[0]=unsigned char(m_b.v.m256d_f64[3]), p[1]=unsigned char(m_g.v.m256d_f64[3]), p[2]=unsigned char(m_r.v.m256d_f64[3]);
		}
		else
		{
			if(c[0]==-1)
				c[0]=colorFunction(r1, i1);
			if(c[1]==-1)
				c[1]=colorFunction(r2, i2);
			if(c[2]==-1)
				c[2]=colorFunction(r3, i3);
			if(c[3]==-1)
				c[3]=colorFunction(r4, i4);
		}//*/
	}
}
//#pragma warning(pop)