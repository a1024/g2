//best viewed with tab size of 4 spaces
//g2_sse2.cpp - Implementation of SSE2 version of math functions.
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

#ifdef _DEBUG
#pragma warning(push)
#pragma warning(disable:4996)
#endif
#include"g2_common.h"
#include"g2_sse2.h"

//	#define	__SSE4_2__
	#define	__SSSE3__
//	#define	__SSE2__

	#define	MAX_VECTOR_SIZE 128
#include	"vectorclass.h"
#include	"vectormath_exp.h"
#include	"vectormath_trig.h"
#include	"vectormath_hyp.h"
#include	"vectormath_common.h"
//#include	<immintrin.h>
#ifdef _DEBUG
#pragma warning(pop)
#endif
const __m128d zero=_mm_setzero_pd(), minus_one=_mm_castsi128_pd(_mm_set1_epi32(-1)), sign_mask=_mm_castsi128_pd(_mm_set_epi32(0x7FFFFFFF, -1, 0x7FFFFFFF, -1));
struct Vect2d
{
	__declspec(align(16)) __m128d v;
	Vect2d():v(zero){}
#ifdef ALIGNED_INTRINSICS
	//Vect2d(VectP const &x)
	//{
	//	if((int&)x.r&31)
	//		v=_mm_loadu_pd(x.r);
	//	else
	//		v=_mm_load_pd(x.r);
	//}
	Vect2d(VectP const &x):v(_mm_load_pd(x.r)){}
	Vect2d(double const *x):v(_mm_load_pd(x)){}
#else
	Vect2d(VectP const &x):v(_mm_loadu_pd(x.r)){}
	Vect2d(double const *x):v(_mm_loadu_pd(x)){}
#endif
	Vect2d(double v):v(_mm_set1_pd(v)){}
	Vect2d(double v0, double v1):v(_mm_set_pd(v1, v0)){}
	Vect2d(__m128d v):v(v){}
	operator __m128d()const{return v;}
	operator __m128i()const{return _mm_castpd_si128(v);}
	operator Vec2d()const{return v;}
	void set(double v0, double v1){v=_mm_set_pd(v1, v0);}
	void setzero(){v=zero;}
	double& lo(){return v.m128d_f64[0];}
	double& hi(){return v.m128d_f64[1];}
	double lo()const{return v.m128d_f64[0];}
	double hi()const{return v.m128d_f64[1];}
	double& get(int component){return v.m128d_f64[component];}
	double get(int component)const{return v.m128d_f64[component];}
	Vect2d r_is_false()const{return _mm_cmpeq_pd(v, zero);}
	Vect2d r_is_true()const{return _mm_cmpneq_pd(v, zero);}
//	Vect2d r_is_true()const{return _mm_xor_pd(_mm_cmpeq_pd(v, zero), minus_one);}
	Vect2d& operator=(__m128d v){this->v=v; return *this;}
	Vect2d& operator+=(Vect2d const &other){v=_mm_add_pd(v, other.v); return *this;}
	Vect2d& operator-=(Vect2d const &other){v=_mm_sub_pd(v, other.v); return *this;}
	Vect2d& operator*=(Vect2d const &other){v=_mm_mul_pd(v, other.v); return *this;}
	Vect2d& operator/=(Vect2d const &other){v=_mm_div_pd(v, other.v); return *this;}
	Vect2d& operator%=(Vect2d const &q)//x%q = x-floor(x/q)*q	SSE4.1
	{
		__m128d v_q=_mm_div_pd(v, q.v);
		v_q=_mm_floor_pd(v_q);
		v_q=_mm_mul_pd(v_q, q.v);
		v=_mm_sub_pd(v, v_q);
		return *this;
	}
	Vect2d& mod_sse2(Vect2d const &q)
	{
		__m128d v_q=_mm_div_pd(v, q.v);
		v_q.m128d_f64[0]=std::floor(v_q.m128d_f64[0]), v_q.m128d_f64[1]=std::floor(v_q.m128d_f64[1]);
		v_q=_mm_mul_pd(v_q, q.v);
		v=_mm_sub_pd(v, v_q);
		return *this;
	}

	Vect2d& operator&=(Vect2d const &other){v=_mm_and_pd(v, other.v); return *this;}
	Vect2d& operator|=(Vect2d const &other){v=_mm_or_pd(v, other.v); return *this;}
	Vect2d& operator^=(Vect2d const &other){v=_mm_xor_pd(v, other.v); return *this;}

	Vect2d& ceil_this(){v=_mm_ceil_pd(v);return *this;}
	Vect2d& ceil_this_sse2()
	{
		v.m128d_f64[0]=std::ceil(v.m128d_f64[0]), v.m128d_f64[1]=std::ceil(v.m128d_f64[1]);
		return *this;
	}
	Vect2d ceil()const{return _mm_ceil_pd(v);}
	Vect2d ceil_sse2()const
	{
		Vect2d x;
		x.v.m128d_f64[0]=std::ceil(v.m128d_f64[0]), x.v.m128d_f64[1]=std::ceil(v.m128d_f64[1]);
		return x;
	}
	Vect2d& floor_this(){v=_mm_floor_pd(v);return *this;}
	Vect2d& floor_this_sse2()
	{
		v.m128d_f64[0]=std::floor(v.m128d_f64[0]), v.m128d_f64[1]=std::floor(v.m128d_f64[1]);
		return *this;
	}
	Vect2d floor()const{return _mm_floor_pd(v);}
	Vect2d floor_sse2()const
	{
		Vect2d x;
		x.v.m128d_f64[0]=std::floor(v.m128d_f64[0]), x.v.m128d_f64[1]=std::floor(v.m128d_f64[1]);
		return x;
	}
	Vect2d& round_this(){v=_mm_round_pd(v, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);return *this;}
	Vect2d& round_this_sse2()
	{
		v.m128d_f64[0]=std::floor(v.m128d_f64[0]+0.5), v.m128d_f64[1]=std::floor(v.m128d_f64[1]+0.5);
		return *this;
	}
	Vect2d round()const{return _mm_round_pd(v, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);}
	Vect2d round_sse2()const
	{
		Vect2d x;
		x.v.m128d_f64[0]=std::floor(v.m128d_f64[0]+0.5), x.v.m128d_f64[1]=std::floor(v.m128d_f64[1]+0.5);
		return x;
	}
	Vect2d complement()const{return Vect2d(_mm_xor_pd(v, minus_one));}
	Vect2d& abs_this(){v=_mm_and_pd(v, sign_mask); return *this;}
	Vect2d abs()const{return Vect2d(_mm_and_pd(v, sign_mask));}
};
const Vect2d m_zero=0., m_one=1, m_two=2, m_minus_one=_mm_castsi128_pd(_mm_set1_epi32(-1)), m_pi=G2::_pi, m_pi_2=G2::_pi_2, m_qnan=G2::_qnan,
	m_sign_mask=sign_mask, m_sign_mask_complement=m_sign_mask.complement(),
	m_half=0.5, m_third=1./3, m_ln2=G2::_ln2, m_ln10=G2::_ln10, m_1_ln10=1/G2::_ln10, m_inf=_HUGE,
	m_one_percent=0.01, m_phi=G2::_phi, m_sqrt5=G2::_sqrt5;
__forceinline Vect2d operator~(Vect2d const &x){return _mm_xor_pd(x, m_minus_one);}
__forceinline Vect2d operator*(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_mul_pd(a.v, b.v));}
__forceinline Vect2d operator/(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_div_pd(a.v, b.v));}
__forceinline Vect2d operator%(Vect2d const &x, Vect2d const &q)//SSE4.1
{
	__m128d mask=_mm_xor_pd(_mm_cmpeq_pd(q.v, zero), minus_one);//zero: 0, nonzero: FFFF
	__m128d x_q=_mm_div_pd(x.v, q.v);
	x_q=_mm_floor_pd(x_q);
	x_q=_mm_mul_pd(x_q, q.v);
	x_q=_mm_sub_pd(x.v, x_q);
	return _mm_and_pd(x_q, mask);
}
__forceinline Vect2d mod_sse2(Vect2d const &x, Vect2d const &q)
{
	__m128d x_q=_mm_div_pd(x.v, q.v);
	x_q.m128d_f64[0]=std::floor(x_q.m128d_f64[0]), x_q.m128d_f64[1]=std::floor(x_q.m128d_f64[1]);
	x_q=_mm_mul_pd(x_q, q.v);
	x_q=_mm_sub_pd(x.v, x_q);
	return x_q;
}
__forceinline Vect2d operator-(Vect2d const &x){return Vect2d(_mm_xor_pd(x.v, m_sign_mask_complement.v));}
__forceinline Vect2d operator+(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_add_pd(a.v, b.v));}
__forceinline Vect2d operator-(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_sub_pd(a.v, b.v));}
__forceinline Vect2d operator>(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_cmpgt_pd(a.v, b.v));}
__forceinline Vect2d operator<(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_cmplt_pd(a.v, b.v));}
__forceinline Vect2d operator>=(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_cmpge_pd(a.v, b.v));}
__forceinline Vect2d operator<=(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_cmple_pd(a.v, b.v));}
__forceinline Vect2d operator==(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_cmpeq_pd(a.v, b.v));}
__forceinline Vect2d operator!=(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_xor_pd(_mm_cmpeq_pd(a.v, b.v), m_minus_one.v));}
__forceinline Vect2d operator&(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_and_pd(a.v, b.v));}
__forceinline Vect2d operator|(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_or_pd(a.v, b.v));}
__forceinline Vect2d operator^(Vect2d const &a, Vect2d const &b){return Vect2d(_mm_xor_pd(a.v, b.v));}
__forceinline Vect2d and(Vect2d const &a, Vect2d const &b){return a&b;}
__forceinline Vect2d sqrt(Vect2d const &x){return Vect2d(_mm_sqrt_pd(x.v));}
__forceinline Vect2d is_not_nan(Vect2d const &x){return _mm_cmpeq_pd(x.v, x.v);}
__forceinline Vect2d isinf(Vect2d const &x){return _mm_cmpeq_pd(x.abs().v, m_inf.v);}
#pragma warning(push)
#pragma warning(disable:4172)
struct Quat2d;
struct Comp2d
{
	Vect2d r, i;
	Comp2d(){}
	Comp2d(CompP const &x):r(x.r), i(x.i){}
	Comp2d(double const *r, double const *i):r(r), i(i){}
//#ifdef ALIGNED_INTRINSICS
//	Comp2d(double const *r, double const *i):r(_mm_load_pd(r)), i(_mm_load_pd(i)){}
//#else
//	Comp2d(double const *r, double const *i):r(_mm_loadu_pd(r)), i(_mm_loadu_pd(i)){}
//#endif
	Comp2d(__m128d const &r, __m128d const &i):r(r), i(i){}
	Comp2d(Vect2d const &r, Vect2d const &i):r(r), i(i){}
	//Comp2d(Vect2d *pr, Vect2d *pi):r(*pr), i(*pi){}
	//void assign(){*pr=r, *pi=i;}
	//void assign(Vect2d const &r, Vect2d const &i){*pr=r, *pi=i;}
	//void assign(Vect2d const &r){*pr=r;}
	//void assign(__m128d const &r){*pr=r;}
	//void assign(__m128d const &r, __m128d const &i){*pr=r, *pi=i;}
	//void assign(Comp2d const &x){*pr=x.r, *pi=x.i;}
	//void setself(){pr=&r, pi=&i;}
	Vect2d r_is_false()const{return r.r_is_false();}
	Vect2d c_is_false()const{return r.r_is_false()&i.r_is_false();}
	Vect2d r_is_true()const{return r!=m_zero;}
	Vect2d c_is_true()const{return ((r==m_zero)&(i==m_zero)).complement();}
	Comp2d floor()const{return Comp2d(r.floor(), i.floor());}
	Comp2d floor_sse2()const{return Comp2d(r.floor_sse2(), i.floor_sse2());}
	Comp2d ceil()const{return Comp2d(r.ceil(), i.ceil());}
	Comp2d ceil_sse2()const{return Comp2d(r.ceil_sse2(), i.ceil_sse2());}
	Comp2d round()const{return Comp2d(r.round(), i.round());}
	Comp2d round_sse2()const{return Comp2d(r.round_sse2(), i.round_sse2());}
	Vect2d abs()const{return sqrt(r*r+i*i);}
	Vect2d arg()const{return Vect2d(::atan2(i.v, r.v));}
	Comp2d& operator+=(Comp2d const &b){r+=b.r, i+=b.i; return *this;}
	Comp2d& operator-=(Comp2d const &b){r-=b.r, i-=b.i; return *this;}
	Comp2d& operator*=(Comp2d const &b)
	{
		Vect2d rr=r*b.r-i*b.i, ri=r*b.i+i*b.r;
		r=rr, i=ri;
		return *this;
	}
	Comp2d& operator*=(Vect2d const &b){r*=b, i*=b; return *this;}
	Comp2d& operator/=(Comp2d const &b)
	{
		Vect2d _1_mag_b=m_one/sqrt(b.r*b.r+b.i*b.i);
		Vect2d
			rr=(b.r*r+b.i*i)*_1_mag_b,
			ri=(b.r*i-b.i*r)*_1_mag_b;
		r=rr, i=ri;
		return *this;
	}
	Comp2d& operator/=(Vect2d const &br){r/=br, i/=br; return *this;}
	Quat2d& operator/=(Quat2d const &b);
	Comp2d& operator^=(Comp2d const &b)
	{
		Comp2d t(::log(sqrt(r*r+i*i).v), atan2(i.v, r.v));
		t*=b;
		Vect2d r0=::exp(t.r.v);
		Vec2d sin_ti, cos_ti;
		sin_ti=sincos(&cos_ti, t.i.v);
		r=r0*Vect2d(cos_ti), i=r0*Vect2d(sin_ti);
		return *this;
	}
	Comp2d& operator^=(Vect2d const &br)
	{
		Comp2d t(::log(sqrt(r*r+i*i).v), atan2(i.v, r.v));
		t*=br;
		Vect2d r0=::exp(t.r.v);
		Vec2d sin_ti, cos_ti;
		sin_ti=sincos(&cos_ti, t.i.v);
		r=r0*Vect2d(cos_ti), i=r0*Vect2d(sin_ti);
		return *this;
	}
};
__forceinline Comp2d operator*(Comp2d const &a, Comp2d const &b){return Comp2d(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r);}
__forceinline Comp2d operator*(Comp2d const &a, Vect2d const &br){return Comp2d(a.r*br, a.i*br);}
__forceinline Comp2d operator*(Vect2d const &ar, Comp2d const &b){return Comp2d(ar*b.r, ar*b.i);}
__forceinline Comp2d operator/(Comp2d const &a, Comp2d const &b)
{
	Vect2d _1_mag_b=m_one/(b.r*b.r+b.i*b.i);
	return Comp2d((b.r*a.r+b.i*a.i)*_1_mag_b, (b.r*a.i-b.i*a.r)*_1_mag_b);
}
__forceinline Comp2d operator/(Comp2d const &a, Vect2d const &br){return Comp2d(a.r/br, a.i/br);}
__forceinline Comp2d operator/(Vect2d const &ar, Comp2d const &b)
{
	Vect2d _ar_mag_b=ar/(b.r*b.r+b.i*b.i);
	return Comp2d(b.r*_ar_mag_b, -b.i*_ar_mag_b);
}
__forceinline Comp2d operator+(Comp2d const &a, Comp2d const &b){return Comp2d(a.r+b.r, a.i+b.i);}
__forceinline Comp2d operator+(Comp2d const &a, Vect2d const &b){return Comp2d(a.r+b, a.i);}
__forceinline Comp2d operator+(Vect2d const &a, Comp2d const &b){return Comp2d(a+b.r, b.i);}
__forceinline Comp2d operator-(Comp2d const &a, Comp2d const &b){return Comp2d(a.r-b.r, a.i-b.i);}
__forceinline Comp2d operator-(Comp2d const &a, Vect2d const &b){return Comp2d(a.r-b, a.i);}
__forceinline Comp2d operator-(Vect2d const &a, Comp2d const &b){return Comp2d(a-b.r, -b.i);}
__forceinline Comp2d operator-(Comp2d const &a){return Comp2d(-a.r, -a.i);}
__forceinline Vect2d operator==(Comp2d const &a, Comp2d const &b){return (a.r==b.r)&(a.i==b.i);}
__forceinline Vect2d operator==(Comp2d const &a, Vect2d const &b){return (a.r==b)&(a.i==m_zero);}
__forceinline Vect2d operator==(Vect2d const &a, Comp2d const &b){return (a==b.r)&(m_zero==b.i);}
__forceinline Vect2d operator!=(Comp2d const &a, Comp2d const &b){return (a.r!=b.r)|(a.i!=b.i);}
__forceinline Vect2d operator!=(Comp2d const &a, Vect2d const &b){return (a.r!=b)|(a.i!=m_zero);}
__forceinline Vect2d operator!=(Vect2d const &a, Comp2d const &b){return (a!=b.r)|(m_zero!=b.i);}
__forceinline Comp2d operator|(Comp2d const &a, Comp2d const &b){return Comp2d(a.r|b.r, a.i|b.i);}
__forceinline Comp2d operator|(Comp2d const &a, Vect2d const &b){return Comp2d(a.r|b, a.i);}
__forceinline Comp2d operator|(Vect2d const &a, Comp2d const &b){return Comp2d(a|b.r, b.i);}
__forceinline Comp2d and(Comp2d const &a, Vect2d const &b){return Comp2d(a.r&b, a.i&b);}
__forceinline Comp2d and(Vect2d const &a, Comp2d const &b){return Comp2d(a&b.r, a&b.i);}
__forceinline Comp2d operator%(Comp2d const &a, Comp2d const &b){return and((a-(a/b).floor()*b), a.c_is_true());}
__forceinline Comp2d operator%(Comp2d const &a, Vect2d const &b){return and((a-(a/b).floor()*b), a.c_is_true());}
__forceinline Comp2d operator%(Vect2d const &a, Comp2d const &b){return and((a-(a/b).floor()*b), a.r_is_true());}
__forceinline Comp2d mod_sse2(Comp2d const &x, Comp2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Comp2d mod_sse2(Comp2d const &x, Vect2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Comp2d mod_sse2(Vect2d const &x, Comp2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Comp2d log(Comp2d const &x)
{
	//auto LOL_3=sqrt(x.r*x.r+x.i*x.i);
	//auto LOL_1=::log(LOL_3.v), LOL_2=::atan2(x.i.v, x.r.v);
//	return Comp2d(::log(sqrt(x.r*x.r+x.i*x.i).v), ::atan((x.i/x.r).v));
	return Comp2d(::log(sqrt(x.r*x.r+x.i*x.i).v), ::atan2(x.i.v, x.r.v));
}
__forceinline Comp2d exp(Comp2d const &x)
{
	Vec2d sin_xi, cos_xi;
	sin_xi=sincos(&cos_xi, x.i.v);
	Vect2d exp_xr=exp(x.r.v);
	return Comp2d(exp_xr*Vect2d(cos_xi), exp_xr*Vect2d(sin_xi));
}
__forceinline Comp2d sqrt(Comp2d const &x)
{
	auto s=x.r+x.abs();
	s=sqrt(s+s);//3 sqrts
	auto i=sqrt(-x.r), mask=s==m_zero;
	return Comp2d(s*m_half, x.i/s&mask.complement()|i&mask);

	//auto mag=x.abs(), _1_mag=m_one/mag, cosx=x.r*_1_mag, sinx=x.i*_1_mag;//4 sqrts
	//Vect2d sign_sinx=sinx>=m_zero;
	//sign_sinx=(sign_sinx&m_one)-(sign_sinx.complement()&m_one);
	//Vect2d cosx_2=sqrt((m_one+cosx)*m_half), sinx_2=sqrt((m_one-cosx)*m_half);
	//mag=sqrt(mag);
	//return Comp2d(mag*cosx_2, mag*sign_sinx*sinx_2);

	//auto LOL_2=m_half*log(x);
	//auto LOL_1=exp(m_half*log(x));
	//return exp(m_half*log(x));//sqrt, log, atan2, exp
}
__forceinline Comp2d operator^(Comp2d const &a, Comp2d const &b)
{
	Vect2d mask=a.r==m_zero&a.i==m_zero&b.r==m_zero&b.i==m_zero, mask_c=mask.complement();
	Comp2d t(::log(sqrt(a.r*a.r+a.i*a.i).v), atan2(a.i.v, a.r.v));
	t*=b;
	Vect2d r0=::exp(t.r.v);
	Vec2d sin_ti, cos_ti;
	sin_ti=sincos(&cos_ti, t.i.v);
	return Comp2d(r0*Vect2d(cos_ti)&mask_c|m_one&mask, r0*Vect2d(sin_ti)&mask_c);
}
__forceinline Comp2d operator^(Comp2d const &a, Vect2d const &br)
{
	Vect2d mask=a.r==m_zero&a.i==m_zero&br==m_zero, mask_c=mask.complement();
	Comp2d t(::log(sqrt(a.r*a.r+a.i*a.i).v), atan2(a.i.v, a.r.v));
	t*=br;
	Vect2d r0=::exp(t.r.v);
	Vec2d sin_ti, cos_ti;
	sin_ti=sincos(&cos_ti, t.i.v);
	return Comp2d(r0*Vect2d(cos_ti)&mask_c|m_one&mask, r0*Vect2d(sin_ti)&mask_c);
}
__forceinline Comp2d operator^(Vect2d const &ar, Comp2d const &b)
{
	Vect2d mask=ar==m_zero&b.r==m_zero&b.i==m_zero, mask_c=mask.complement();
	Comp2d t(::log(ar.abs().v), atan2(zero, ar.v));
	t*=b;
	Vect2d r0=::exp(t.r.v);
	Vec2d sin_ti, cos_ti;
	sin_ti=sincos(&cos_ti, t.i.v);
	return Comp2d(r0*Vect2d(cos_ti)&mask_c|m_one&mask, r0*Vect2d(sin_ti)&mask_c);
}
struct Quat2d
{
	Vect2d r, i, j, k;
	Quat2d(){}
	Quat2d(QuatP const &x):r(x.r), i(x.i), j(x.j), k(x.k){}
	Quat2d(double const *r, double const *i, double const *j, double const *k):r(r), i(i), j(j), k(k){}
//#ifdef ALIGNED_INTRINSICS
//	Quat2d(double const *r, double const *i, double const *j, double const *k):r(_mm_load_pd(r)), i(_mm_load_pd(i)), j(_mm_load_pd(j)), k(_mm_load_pd(k)){}
//#else
//	Quat2d(double const *r, double const *i, double const *j, double const *k):r(_mm_loadu_pd(r)), i(_mm_loadu_pd(i)), j(_mm_loadu_pd(j)), k(_mm_loadu_pd(k)){}
//#endif
	Quat2d(__m128d const &r, __m128d const &i, __m128d const &j, __m128d const &k):r(r), i(i), j(j), k(k){}
	Quat2d(Vect2d const &r, Vect2d const &i, Vect2d const &j, Vect2d const &k):r(r), i(i), j(j), k(k){}
	//Quat2d(Vect2d *pr, Vect2d *pi, Vect2d *pj, Vect2d *pk):r(*pr), i(*pi), j(*pj), k(*pk), pr(pr), pi(pi), pj(pj), pk(pk){}
	Quat2d(Comp2d const &x):r(x.r), i(x.i){}
	//void assign(){*pr=r, *pi=i, *pj=j, *pk=k;}
	//void assign(Vect2d const &r){*pr=r;}
	//void assign(Vect2d const &r, Vect2d const &i){*pr=r, *pi=i;}
	//void assign(Vect2d const &r, Vect2d const &i, Vect2d const &j, Vect2d const &k){*pr=r, *pi=i, *pj=j, *pk=k;}
	//void assign(Comp2d const &x){*pr=x.r, *pi=x.i;}
	//void assign(Quat2d const &x){*pr=x.r, *pi=x.i, *pj=x.j, *pk=x.k;}
	//void assign(__m128d const &r){*pr=r;}
	//void assign(__m128d const &r, __m128d const &i){*pr=r, *pi=i;}
	//void assign(__m128d const &r, __m128d const &i, __m128d const &j, __m128d const &k){*pr=r, *pi=i, *pj=j, *pk=k;}
	void setzero(){r.setzero(), i.setzero(), j.setzero(), k.setzero();}
	Vect2d r_is_false()const{return r.r_is_false();}
	Vect2d c_is_false()const{return r.r_is_false()&i.r_is_false();}
	Vect2d q_is_false()const{return r.r_is_false()&i.r_is_false()&j.r_is_false()&k.r_is_false();}
	Vect2d r_is_true()const{return r!=m_zero;}
	Vect2d c_is_true()const{return ((r==m_zero)&(i==m_zero)).complement();}
	Vect2d q_is_true()const{return ((r==m_zero)&(i==m_zero)&(j==m_zero)&(k==m_zero)).complement();}
	Quat2d floor()const{return Quat2d(r.floor(), i.floor(), j.floor(), k.floor());}
	Quat2d floor_sse2()const{return Quat2d(r.floor_sse2(), i.floor_sse2(), j.floor_sse2(), k.floor_sse2());}
	Quat2d ceil()const{return Quat2d(r.ceil(), i.ceil(), j.ceil(), k.ceil());}
	Quat2d ceil_sse2()const{return Quat2d(r.ceil_sse2(), i.ceil_sse2(), j.ceil_sse2(), k.ceil_sse2());}
	Quat2d round()const{return Quat2d(r.round(), i.round(), j.round(), k.round());}
	Quat2d round_sse2()const{return Quat2d(r.round_sse2(), i.round_sse2(), j.round_sse2(), k.round_sse2());}
	Vect2d abs()const{return sqrt(r*r+i*i+j*j+k*k);}
	Quat2d& operator+=(Quat2d const &b){r+=b.r, i+=b.i, j+=b.j, k+=b.k;return *this;}
	Quat2d& operator+=(Comp2d const &b){r+=b.r, i+=b.i;return *this;}
	Quat2d& operator+=(Vect2d const &br){r+=br;return *this;}
	Quat2d& operator-=(Quat2d const &b){r-=b.r, i-=b.i, j-=b.j, k-=b.k;return *this;}
	Quat2d& operator-=(Comp2d const &b){r-=b.r, i-=b.i;return *this;}
	Quat2d& operator-=(Vect2d const &br){r-=br;return *this;}
	Quat2d& operator*=(Quat2d const &b)
	{
		Vect2d
			rr=r*b.r+i*b.i+j*b.j+k*b.k,
			ri=r*b.i+i*b.r+j*b.k-k*b.j,
			rj=r*b.j-i*b.k+j*b.r+k*b.i,
			rk=r*b.k+i*b.j+j*b.i+k*b.r;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat2d& operator*=(Comp2d const &b)
	{
		Vect2d
			rr=r*b.r+i*b.i,
			ri=r*b.i+i*b.r,
			rj=j*b.r+k*b.i,
			rk=j*b.i+k*b.r;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat2d& operator*=(Vect2d const &b){r*=b, i*=b, j*=b, k*=b;return *this;}
	Quat2d& operator/=(Quat2d const &b)
	{
		Vect2d _1_mag_y=m_one/sqrt(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
		Vect2d
			rr=(b.r*r+b.i*i+b.j*j+b.k*k)*_1_mag_y,
			ri=(b.r*i-b.i*r-b.j*k+b.k*j)*_1_mag_y,
			rj=(b.r*j+b.i*k-b.j*r-b.k*i)*_1_mag_y,
			rk=(b.r*k-b.i*j+b.j*i-b.k*r)*_1_mag_y;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat2d& operator/=(Comp2d const &b)
	{
		Vect2d _1_mag_y=m_one/sqrt(b.r*b.r+b.i*b.i);
		Vect2d
			rr=(b.r*r+b.i*i)*_1_mag_y,
			ri=(b.r*i-b.i*r)*_1_mag_y,
			rj=(b.r*j+b.i*k)*_1_mag_y,
			rk=(b.r*k-b.i*j)*_1_mag_y;
		r=rr, i=ri, j=rj, k=rk;
		return *this;
	}
	Quat2d& operator/=(Vect2d const &b){r/=b, i/=b, j/=b, k/=b;return *this;}
	Quat2d& operator^=(Quat2d const &b)
	{
		Vect2d mag_v=i*i+j*j+k*k;
		Vect2d ln_mag_a=log(sqrt(r*r+mag_v).v);
		mag_v=sqrt(mag_v);
		Vect2d v_mul=acos((r/ln_mag_a).v);
		v_mul/=mag_v;
		Quat2d t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
		t=Quat2d(
			b.r*t.r+b.i*t.i+b.j*t.j+b.k*t.k,
			b.r*t.i+b.i*t.r+b.j*t.k-b.k*t.j,
			b.r*t.j-b.i*t.k+b.j*t.r+b.k*t.i,
			b.r*t.k+b.i*t.j+b.j*t.i+b.k*t.r);
		Vect2d mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
		Vec2d sin_mu, cos_mu;
		sin_mu=sincos(&cos_mu, mag_u.v);
		Vect2d exp_tr=exp(t.r.v);
		v_mul=exp_tr*Vect2d(sin_mu)/mag_u;
		r=exp_tr*Vect2d(cos_mu), i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
		return *this;
	}
	Quat2d& operator^=(Comp2d const &b)
	{
		Vect2d mag_v=i*i+j*j+k*k;
		Vect2d ln_mag_a=log(sqrt(r*r+mag_v).v);
		mag_v=sqrt(mag_v);
		Vect2d v_mul=acos((r/ln_mag_a).v);
		v_mul/=mag_v;
		Quat2d t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
		t=Quat2d(
			b.r*t.r+b.i*t.i,
			b.r*t.i+b.i*t.r,
			b.r*t.j-b.i*t.k,
			b.r*t.k+b.i*t.j);
		Vect2d mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
		Vec2d sin_mu, cos_mu;
		sin_mu=sincos(&cos_mu, mag_u.v);
		Vect2d exp_tr=exp(t.r.v);
		v_mul=exp_tr*Vect2d(sin_mu)/mag_u;
		r=exp_tr*Vect2d(cos_mu), i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
		return *this;
	}
	Quat2d& operator^=(Vect2d const &br)
	{
		Vect2d mag_v=i*i+j*j+k*k;
		Vect2d ln_mag_a=log(sqrt(r*r+mag_v).v);
		mag_v=sqrt(mag_v);
		Vect2d v_mul=acos((r/ln_mag_a).v);
		v_mul/=mag_v;
		Quat2d t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
		t=Quat2d(br*t.r, br*t.i, br*t.j, br*t.k);
		Vect2d mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
		Vec2d sin_mu, cos_mu;
		sin_mu=sincos(&cos_mu, mag_u.v);
		Vect2d exp_tr=exp(t.r.v);
		v_mul=exp_tr*Vect2d(sin_mu)/mag_u;
		r=exp_tr*Vect2d(cos_mu), i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
		return *this;
	}
};
__forceinline Vect2d operator==(Quat2d const &a, Quat2d const &b){return (a.r==b.r)&(a.i==b.i)&(a.j==b.j)&(a.k==b.k);}
__forceinline Vect2d operator==(Quat2d const &a, Comp2d const &b){return (a.r==b.r)&(a.i==b.i)&(a.j==m_zero)&(a.k==m_zero);}
__forceinline Vect2d operator==(Quat2d const &a, Vect2d const &b){return (a.r==b)&(a.i==m_zero)&(a.j==m_zero)&(a.k==m_zero);}
__forceinline Vect2d operator==(Comp2d const &a, Quat2d const &b){return (a.r==b.r)&(a.i==b.i)&(m_zero==b.j)&(m_zero==b.k);}
__forceinline Vect2d operator==(Vect2d const &a, Quat2d const &b){return (a==b.r)&(m_zero==b.i)&(m_zero==b.j)&(m_zero==b.k);}

__forceinline Vect2d operator!=(Quat2d const &a, Quat2d const &b){return (a.r!=b.r)|(a.i!=b.i)|(a.j!=b.j)|(a.k!=b.k);}
__forceinline Vect2d operator!=(Quat2d const &a, Comp2d const &b){return (a.r!=b.r)|(a.i!=b.i)|(a.j!=m_zero)|(a.k!=m_zero);}
__forceinline Vect2d operator!=(Quat2d const &a, Vect2d const &b){return (a.r!=b)|(a.i!=m_zero)|(a.j!=m_zero)|(a.k!=m_zero);}
__forceinline Vect2d operator!=(Comp2d const &a, Quat2d const &b){return (a.r!=b.r)|(a.i!=b.i)|(m_zero!=b.j)|(m_zero!=b.k);}
__forceinline Vect2d operator!=(Vect2d const &a, Quat2d const &b){return (a!=b.r)|(m_zero!=b.i)|(m_zero!=b.j)|(m_zero!=b.k);}

__forceinline Quat2d operator+(Quat2d const &a , Quat2d const &b ){return Quat2d(a.r+b.r, a.i+b.i, a.j+b.j, a.k+b.k);}
__forceinline Quat2d operator+(Quat2d const &a , Comp2d const &b ){return Quat2d(a.r+b.r, a.i+b.i, a.j, a.k);}
__forceinline Quat2d operator+(Quat2d const &a , Vect2d const &br){return Quat2d(a.r+br, a.i, a.j, a.k);}
__forceinline Quat2d operator+(Comp2d const &a , Quat2d const &b ){return Quat2d(a.r+b.r, a.i+b.i, b.j, b.k);}
__forceinline Quat2d operator+(Vect2d const &ar, Quat2d const &b ){return Quat2d(ar+b.r, b.i, b.j, b.k);}

__forceinline Quat2d operator-(Quat2d const &a , Quat2d const &b ){return Quat2d(a.r-b.r, a.i-b.i, a.j-b.j, a.k-b.k);}
__forceinline Quat2d operator-(Quat2d const &a , Comp2d const &b ){return Quat2d(a.r-b.r, a.i-b.i, a.j, a.k);}
__forceinline Quat2d operator-(Quat2d const &a , Vect2d const &br){return Quat2d(a.r-br, a.i, a.j, a.k);}
__forceinline Quat2d operator-(Comp2d const &a , Quat2d const &b ){return Quat2d(a.r-b.r, a.i-b.i, -b.j, -b.k);}
__forceinline Quat2d operator-(Vect2d const &ar, Quat2d const &b ){return Quat2d(ar-b.r, -b.i, -b.j, -b.k);}
__forceinline Quat2d operator-(Quat2d const &a){return Quat2d(-a.r, -a.i, -a.j, -a.k);}
__forceinline Quat2d operator*(Quat2d const &a, Quat2d const &b)
{
	return Quat2d(
		a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,
		a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,
		a.r*b.j-a.i*b.k+a.j*b.r+a.k*b.i,
		a.r*b.k+a.i*b.j-a.j*b.i+a.k*b.r);
}
__forceinline Quat2d operator*(Quat2d const &a, Comp2d const &b)
{
	return Quat2d(
		a.r*b.r-a.i*b.i,
		a.r*b.i+a.i*b.r,
		a.j*b.r+a.k*b.i,
		-a.j*b.i+a.k*b.r);
}
__forceinline Quat2d operator*(Quat2d const &a, Vect2d const &br){return Quat2d(a.r*br, a.i*br, a.j*br, a.k*br);}
__forceinline Quat2d operator*(Comp2d const &a, Quat2d const &b)
{
	return Quat2d(
		a.r*b.r-a.i*b.i,
		a.r*b.i+a.i*b.r,
		a.r*b.j-a.i*b.k,
		a.r*b.k+a.i*b.j);
}
__forceinline Quat2d operator*(Vect2d const &ar, Quat2d const &b){return Quat2d(ar*b.r, ar*b.i, ar*b.j, ar*b.k);}
__forceinline Quat2d operator/(Quat2d const &a, Quat2d const &b)
{
	Vect2d _1_mag_y=m_one/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat2d(
		(b.r*a.r+b.i*a.i+b.j*a.j+b.k*a.k)*_1_mag_y,
		(b.r*a.i-b.i*a.r-b.j*a.k+b.k*a.j)*_1_mag_y,
		(b.r*a.j+b.i*a.k-b.j*a.r-b.k*a.i)*_1_mag_y,
		(b.r*a.k-b.i*a.j+b.j*a.i-b.k*a.r)*_1_mag_y);
}
__forceinline Quat2d& Comp2d::operator/=(Quat2d const &b)
{
	Vect2d _1_mag_y=m_one/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat2d(
		(b.r*r+b.i*i)*_1_mag_y,
		(b.r*i-b.i*r)*_1_mag_y,
		(-b.j*r-b.k*i)*_1_mag_y,
		(b.j*i-b.k*r)*_1_mag_y);
}
__forceinline Quat2d operator/(Quat2d const &a, Comp2d const &b)
{
	Vect2d _1_mag_y=m_one/(b.r*b.r+b.i*b.i);
	return Quat2d(
		(b.r*a.r+b.i*a.i)*_1_mag_y,
		(b.r*a.i-b.i*a.r)*_1_mag_y,
		(b.r*a.j+b.i*a.k)*_1_mag_y,
		(b.r*a.k-b.i*a.j)*_1_mag_y);
}
__forceinline Quat2d operator/(Quat2d const &a, Vect2d const &br)
{
	Vect2d _1_mag_y=m_one/br;
	return Quat2d(a.r*_1_mag_y, a.i*_1_mag_y, a.j*_1_mag_y, a.k*_1_mag_y);
}
__forceinline Quat2d operator/(Comp2d const &a, Quat2d const &b)
{
	Vect2d _1_mag_y=m_one/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat2d(
		(b.r*a.r+b.i*a.i)*_1_mag_y,
		(b.r*a.i-b.i*a.r)*_1_mag_y,
		(-b.j*a.r-b.k*a.i)*_1_mag_y,
		(b.j*a.i-b.k*a.r)*_1_mag_y);
}
__forceinline Quat2d operator/(Vect2d const &ar, Quat2d const &b)
{
	Vect2d _ar_mag_y=ar/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return Quat2d(b.r*_ar_mag_y, -b.i*_ar_mag_y, -b.j*_ar_mag_y, -b.k*_ar_mag_y);
}
__forceinline Quat2d log(Quat2d const &x)
{
	Vect2d mag_v=x.i*x.i+x.j*x.j+x.k*x.k;

	Vect2d real=mag_v==m_zero, real_c=real.complement(), rr=log(x.r.abs().v);

	Vect2d mag_x=sqrt(x.r*x.r+mag_v);
	mag_v=sqrt(mag_v);
	Vect2d u_mul=Vect2d(::acos((x.r/mag_x).v))/mag_v&real_c;
	Vect2d ln_mag_x=::log(mag_x.v);
	return Quat2d(rr&real|ln_mag_x&real_c, m_pi&real|x.i*u_mul, x.j*u_mul, x.k*u_mul);
}
__forceinline Quat2d exp(Quat2d const &x)
{
	Vect2d exp_r=::exp(x.r.v);
	Vect2d mag_v=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
	Vec2d sin_v, cos_v;
	sin_v=sincos(&cos_v, mag_v.v);
	Vect2d v_mul=exp_r*Vect2d(sin_v)/mag_v;
	return Quat2d(exp_r*Vect2d(cos_v), x.i*v_mul, x.j*v_mul, x.k*v_mul);
}
__forceinline Quat2d sqrt(Quat2d const &x)
{
	auto s=x.r+x.abs();
	s=sqrt(s+s);
	auto i=sqrt(-x.r), mask=s==m_zero, mask_c=mask.complement();
	auto _1_s=m_one/s&mask_c;
	return Quat2d(s*m_half, x.i*_1_s|i&mask, x.j*_1_s, x.k*_1_s);
//	return exp(m_half*log(x));
}
__forceinline Quat2d operator^(Quat2d const &a, Quat2d const &b)
{
	Vect2d mask=a.r==m_zero&a.i==m_zero&a.j==m_zero&a.k==m_zero&b.r==m_zero&b.i==m_zero&b.j==m_zero&b.k==m_zero, mask_c=mask.complement();
	Quat2d q=exp(b*log(a));
	return Quat2d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat2d operator^(Quat2d const &a, Comp2d const &b)
{
	Vect2d mask=a.r==m_zero&a.i==m_zero&a.j==m_zero&a.k==m_zero&b.r==m_zero&b.i==m_zero, mask_c=mask.complement();
	Quat2d q=exp(b*log(a));
	return Quat2d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat2d operator^(Quat2d const &a, Vect2d const &b)
{
	Vect2d mask=a.r==m_zero&a.i==m_zero&a.j==m_zero&a.k==m_zero&b==m_zero, mask_c=mask.complement();
	Quat2d q=exp(b*log(a));
	return Quat2d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat2d operator^(Comp2d const &a, Quat2d const &b)
{
	Vect2d mask=a.r==m_zero&a.i==m_zero&b.r==m_zero&b.i==m_zero&b.j==m_zero&b.k==m_zero, mask_c=mask.complement();
	Quat2d q=exp(b*log(a));
	return Quat2d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat2d operator^(Vect2d const &a, Quat2d const &b)
{
	Vect2d mask=a==m_zero&b.r==m_zero&b.i==m_zero&b.j==m_zero&b.k==m_zero, mask_c=mask.complement();
	Quat2d q=exp(b*Vect2d(::log(a.v)));
	return Quat2d(q.r&mask_c|m_one&mask, q.i&mask_c, q.j&mask_c, q.k&mask_c);
}
__forceinline Quat2d operator|(Quat2d const &a, Quat2d const &b){return Quat2d(a.r|b.r, a.i|b.i, a.j|b.j, a.k|b.k);}
__forceinline Quat2d operator|(Quat2d const &a, Comp2d const &b){return Quat2d(a.r|b.r, a.i|b.i, a.j, a.k);}
__forceinline Quat2d operator|(Quat2d const &a, Vect2d const &b){return Quat2d(a.r|b, a.i, a.j, a.k);}
__forceinline Quat2d operator|(Comp2d const &a, Quat2d const &b){return Quat2d(a.r|b.r, a.i|b.i, b.j, b.k);}
__forceinline Quat2d operator|(Vect2d const &a, Quat2d const &b){return Quat2d(a|b.r, b.i, b.j, b.k);}
__forceinline Quat2d and(Quat2d const &a, Vect2d const &b){return Quat2d(a.r&b, a.i&b, a.j&b, a.k&b);}
__forceinline Quat2d and(Vect2d const &a, Quat2d const &b){return Quat2d(a&b.r, a&b.i, a&b.j, a&b.k);}
__forceinline Quat2d operator%(Quat2d const &a, Quat2d const &b){return and((a-(a/b).floor()*b), a.q_is_true());}
__forceinline Quat2d operator%(Quat2d const &a, Comp2d const &b){return and((a-(a/b).floor()*b), a.q_is_true());}
__forceinline Quat2d operator%(Quat2d const &a, Vect2d const &b){return and((a-(a/b).floor()*b), a.q_is_true());}
__forceinline Quat2d operator%(Comp2d const &a, Quat2d const &b){return and((a-(a/b).floor()*b), a.c_is_true());}
__forceinline Quat2d operator%(Vect2d const &a, Quat2d const &b){return and((a-(a/b).floor()*b), a.r_is_true());}
__forceinline Quat2d mod_sse2(Quat2d const &x, Quat2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Quat2d mod_sse2(Quat2d const &x, Comp2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Quat2d mod_sse2(Quat2d const &x, Vect2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Quat2d mod_sse2(Comp2d const &x, Quat2d const &q){return x-(x/q).floor_sse2()*q;}
__forceinline Quat2d mod_sse2(Vect2d const &x, Quat2d const &q){return x-(x/q).floor_sse2()*q;}
#pragma warning(pop)
//struct CompRef2d
//{
//	Vect2d &r, &i;
//	CompRef2d(Vect2d &r, Vect2d &i):r(r), i(i){}
//	CompRef2d& operator=(Comp2d const &other){r=other.r, i=other.i; return *this;}
//};
//struct QuatRef2d
//{
//	Vect2d &r, &i, &j, &k;
//	QuatRef2d(Vect2d &r, Vect2d &i, Vect2d &j, Vect2d &k):r(r), i(i), j(j), k(k){}
//	QuatRef2d& operator=(Quat2d &other){r=other.r, i=other.i, j=other.j, k=other.k; return *this;}
//};
void minmax_sse2(double *a, int size, double *lo_hi)//{min, max}, even size, 16 byte aligned
{
	//if((int)a&15)
	//	return;

	__m128d min=_mm_load_pd(a), max=min;
	for(int k=2;k+1<size;k+=2)				//1.487 c/v on R7 2700 (fastest)
	{
		__m128d vk=_mm_load_pd(a+k);		//3.494 cycles/value
		min=_mm_min_pd(min, vk);
		max=_mm_max_pd(max, vk);
	}
	__m128d
		t0=_mm_shuffle_pd(min, max, 0),//{max[0], min[0]}	//{hiqword, loqword}
		t1=_mm_shuffle_pd(min, max, 3);//{max[1], min[1]}
	__m128d t_min=_mm_min_pd(t0, t1), t_max=_mm_max_pd(t0, t1),
		t_minmax=_mm_shuffle_pd(t_min, t_max, 2);
	_mm_store_pd(lo_hi, t_minmax);
	//lo=min.m128d_f64[0]<min.m128d_f64[1]?min.m128d_f64[0]:min.m128d_f64[1];
	//hi=max.m128d_f64[0]>max.m128d_f64[1]?max.m128d_f64[0]:max.m128d_f64[1];

/*	__declspec(align(16)) double lo2[2], hi2[2];//size multiple of 4, aligned by 16 bytes
		 if(a[0]<a[1])	lo2[0]=a[0], hi2[0]=a[1];
	else				lo2[0]=a[1], hi2[0]=a[0];
		 if(a[2]<a[3])	lo2[1]=a[2], hi2[1]=a[3];
	else				lo2[1]=a[3], hi2[1]=a[2];
	//__m128d min=_mm_load_pd(lo2), max=_mm_load_pd(hi2);					//17.462 cycles/value [sic]
	__m128d min=_mm_set_pd(lo2[1], lo2[0]), max=_mm_set_pd(hi2[1], hi2[0]);	// 7.220 7.277 cycles/value
	for(int k=4;k+3<size;k+=2)
	{
		if(a[k]<a[k+1])
			lo2[0]=a[k], hi2[0]=a[k+1];
		else
			lo2[0]=a[k+1], hi2[0]=a[k];
		if(a[k+2]<a[k+3])
			lo2[1]=a[k+2], hi2[1]=a[k+3];
		else
			lo2[1]=a[k+3], hi2[1]=a[k+2];
		//min=_mm_min_pd(min, _mm_load_pd(lo2));
		//max=_mm_max_pd(max, _mm_load_pd(hi2));
		min=_mm_min_pd(min, _mm_set_pd(lo2[1], lo2[0]));
		max=_mm_max_pd(max, _mm_set_pd(hi2[1], hi2[0]));
	}
	lo=min.m128d_f64[0]<min.m128d_f64[1]?min.m128d_f64[0]:min.m128d_f64[1];
	hi=max.m128d_f64[0]>max.m128d_f64[1]?max.m128d_f64[0]:max.m128d_f64[1];//*/
}
namespace	G2
{
/*	bool _2d_between(double x1, double y1, double x, double y, double x2, double y2)
	{
		if(x1==x2)
		{
			if(x==x1)
			{
				if(y1==y2)
					return x==x1&&y==y1;
				else
					return y1<y2?y1<=y&&y<=y2:y2<=y&&y<=y1;
			}
		}
		else
		{
			if(y1==y2)
			{
				if(y==y1)
					return x1<x2?x1<=x&&x<=x2:x2<=x&&x<=x1;
			}
			else if(x==x1)
				return y==y1;
			else if(x==x2)
				return y==y2;
			else if(y==y1)
				return x==x1;
			else if(y==y2)
				return x==x2;
			else if((y2-y1)*(x-x1)==(x2-x1)*(y-y1))
			{
				if(x1<x2)
				{
					if(x1<=x&&x<=x2)
						return y1<y2?y1<=y&&y<=y2:y2<=y&&y<=y1;
				}
				else
				{
					if(x2<=x&&x<=x1)
						return y1<y2?y1<=y&&y<=y2:y2<=y&&y<=y1;
				}
			}
		}
		return false;
	}
	bool _1d_int_in_range(double x0, double x1){return std::floor(x0)!=std::floor(x1)||std::ceil(x0)!=std::ceil(x1);}
	bool _1d_zero_in_range(double x0, double x1){return x0<0?x1>=0:x0==0?x1<0||x1>0:x1<0;}
	double _1d_zero_crossing(double x0, double y0, double x1, double y1){return x0+(0-y0)*(x1-x0)/(y1-y0);}
	const double ll_max=9.22337203685478e+018;
	__forceinline long long convert_d2ll		(double const x){return x!=x||x<-ll_max||x>ll_max?(long long&)x	:(long long)x;}
	__forceinline long long convert_d2ll_zero	(double const x){return x!=x||x<-ll_max||x>ll_max?0				:(long long)x;}

	std::complex<double>			floor					(std::complex<double>				const &x){return std::complex<double>(::floor(x.real()), ::floor(x.imag()));}
	boost::math::quaternion<double>	floor					(boost::math::quaternion<double>	const &x){return boost::math::quaternion<double>(::floor(x.R_component_1()), ::floor(x.R_component_2()), ::floor(x.R_component_3()), ::floor(x.R_component_4()));}
	double							arg						(boost::math::quaternion<double>	const &x){return ::acos(x.R_component_1()/boost::math::abs(x));}
	boost::math::quaternion<double>	log						(boost::math::quaternion<double>	const &x)
	{
		try
		{
			double t=boost::math::abs(boost::math::quaternion<double>(0, x.R_component_2(), x.R_component_3(), x.R_component_4()));
			if(t)
			{
				t=arg(x)/t;
				return boost::math::quaternion<double>(boost::math::log1p(abs(x)-1), t*x.R_component_2(), t*x.R_component_3(), t*x.R_component_4());
			}
			if(x.R_component_1()!=0)
				return boost::math::quaternion<double>(boost::math::log1p(abs(x)-1), 0, 0, 0);
			return boost::math::quaternion<double>(-_HUGE);
		}
		catch(std::overflow_error&)
		{
			return boost::math::quaternion<double>(-_HUGE);
		}
	}
	boost::math::quaternion<double>	pow						(boost::math::quaternion<double>	const &x, boost::math::quaternion<double>	const &y){return boost::math::exp(y*log(x));}
	boost::math::quaternion<double>	pow						(double								const &x, boost::math::quaternion<double>	const &y){return boost::math::exp(y*std::log(std::complex<double>(x)));}
	double							bitwise_xor				(double								const &x)
	{
		long long t1=convert_d2ll(x);
	//	const double ll_max=9.22337203685478e+018;
	//	long long t1=x!=x||x<-ll_max||x>ll_max?(long long&)x:(long long)x;

		//const double ll_max=9.22337203685478e+018;
		//if(x!=x||x<-ll_max||x>ll_max)
		//	return _qnan;
		//long long t1=(long long)x;//high truncation		breaks font & TabbedTextOut

		//long long t1=(long long&)x;

		t1^=t1>>32, t1^=t1>>16, t1^=t1>>8, t1^=t1>>4;
		t1&=15;
		return (0x6996>>t1)&1;
	//	return 0;
	}
	double							bitwise_xnor			(double								const &x)
	{
		long long t1=convert_d2ll(x);
		t1^=t1>>32, t1^=t1>>16, t1^=t1>>8, t1^=t1>>4;
		t1&=15;
		return ~((0x6996>>t1)&1);
	}

	double							sgn						(double								const &x){return (x>0)-(x<0);}
	std::complex<double>			sgn						(std::complex<double>				const &x)
	{
		double temp=std::abs(x);
		return temp?x/temp:std::complex<double>();
	}
	boost::math::quaternion<double>	sgn						(boost::math::quaternion<double>	const &x)
	{
		double q2=abs(x);
		if(q2!=0)
			return x/q2;
		return boost::math::quaternion<double>();
	}
	boost::math::quaternion<double>	sgnu					(boost::math::quaternion<double>	const &x){return sgn(boost::math::quaternion<double>(0, x.R_component_2(), x.R_component_3(), x.R_component_4()));}
	namespace gamma//http://en.wikipedia.org/wiki/Lanczos_approximation
	{
		const double g=7, p[]={0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	}
	double							tgamma					(double								const &x)
	{
		try
		{
			return boost::math::tgamma(x);
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
	std::complex<double>			tgamma					(std::complex<double>				const &x)
	{
		using namespace gamma;
		if(x.real()<.5)
		{
			std::complex<double> t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			std::complex<double> t2=g+.5-x;
			return _pi/(sin(_pi*x)*_sqrt_2pi*std::pow(t2, .5-x)*std::exp(-t2)*t1);
		}
		else
		{
			std::complex<double> t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			std::complex<double> t2=g+.5+x-1.;
			return _sqrt_2pi*std::pow(t2, .5+x-1.)*std::exp(-t2)*t1;
		}
		//if(x.real()<.5)
		//	return _pi/(sin(_pi*x)*tgamma(1.-x));
		//else
		//{
		//	std::complex<double> t1(p[0]);
		//	for(int k=1;k<g+2;++k)
		//		t1+=p[k]/(x-1.+double(k));
		//	std::complex<double> t2=x-1.+g+.5;
		//	return _sqrt_2pi*std::pow(t2, x-1.+.5)*std::exp(-t2)*t1;
		//}
	}
	boost::math::quaternion<double>	tgamma					(boost::math::quaternion<double>	const &x)
	{
		using namespace gamma;
		if(x.real()<.5)
		{
			boost::math::quaternion<double> t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			boost::math::quaternion<double> t2=g+.5-x;
			return _pi/(sin(_pi*x)*_sqrt_2pi*pow(t2, .5-x)*exp(-t2)*t1);
		}
		else
		{
			boost::math::quaternion<double> t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			boost::math::quaternion<double> t2=g+.5+x-1.;
			return _sqrt_2pi*pow(t2, .5+x-1.)*exp(-t2)*t1;
		}
	}
	boost::math::quaternion<double>	sqrt					(boost::math::quaternion<double>	const &x){return exp(0.5*log(x));}
	boost::math::quaternion<double>	acosh					(boost::math::quaternion<double>	const &x){return log(x+sqrt(x*x-1.));}
	boost::math::quaternion<double>	asinh					(boost::math::quaternion<double>	const &x){return log(x+sqrt(x*x+1.));}
	boost::math::quaternion<double>	atanh					(boost::math::quaternion<double>	const &x){return (log(1.+x)-log(1.-x))/2.;}
	boost::math::quaternion<double>	acos					(boost::math::quaternion<double>	const &x){return -sgnu(x)*acosh(x);}
	boost::math::quaternion<double>	asin					(boost::math::quaternion<double>	const &x)
	{
		boost::math::quaternion<double> q2=sgnu(x);
		return -q2*asinh(x*q2);
	}
	boost::math::quaternion<double>	atan					(boost::math::quaternion<double>	const &x)
	{
		boost::math::quaternion<double> q2=sgnu(x);
		return -q2*atanh(x*q2);
	}
	double							step					(double								const &x){return .5+.5*sgn(x);}
	std::complex<double>			step					(std::complex<double>				const &x){return .5+.5*sgn(x);}
	boost::math::quaternion<double>	step					(boost::math::quaternion<double>	const &x){return .5+.5*sgn(x);}//*/
	namespace sse2
	{
		//__forceinline void assign(VectP &p, Vect2d const &v)
		//{
		//	_mm_storeu_pd(p.r, v);
		//}
		////__forceinline void assign(double *p, Vect2d const &v){_mm_storeu_pd(p, v.v);}
		//__forceinline void assign(CompP &p, Comp2d const &v)
		//{
		//	_mm_storeu_pd(p.r, v.r);//VS2010 bug
		//	_mm_storeu_pd(p.i, v.i);
		//}
		//__forceinline void assign(QuatP &p, Quat2d const &v){_mm_storeu_pd(p.r, v.r), _mm_storeu_pd(p.i, v.i), _mm_storeu_pd(p.j, v.j), _mm_storeu_pd(p.k, v.k);}
		__forceinline void assign(VectP &p, Vect2d const &v)
		{
			//if((int&)p.r&31)
			//	_mm_storeu_pd(p.r, v);
			//else
			//	_mm_store_pd(p.r, v);
			_mm_store_pd(p.r, v);
		}
		//__forceinline void assign(double *p, Vect2d const &v){_mm_store_pd(p, v.v);}
		__forceinline void assign(CompP &p, Comp2d const &v)
		{
			_mm_store_pd(p.r, v.r);//VS2010 bug
			_mm_store_pd(p.i, v.i);
		}
		__forceinline void assign(QuatP &p, Quat2d const &v)
		{
			_mm_store_pd(p.r, v.r);//VS2010 bug
			_mm_store_pd(p.i, v.i);
			_mm_store_pd(p.j, v.j);
			_mm_store_pd(p.k, v.k);
		}

		void r_r_setzero				(VectP &r, VectP const&)					{assign(r, m_zero);}
		void c_c_setzero				(CompP &r, CompP const&)					{assign(r, Comp2d(_mm_setzero_pd(), _mm_setzero_pd()));}
		void q_q_setzero				(QuatP &r, QuatP const&)					{assign(r, Quat2d(_mm_setzero_pd(), _mm_setzero_pd(), _mm_setzero_pd(), _mm_setzero_pd()));}

		void  r_r_ceil					(VectP &r, VectP const &x)					{assign(r, Vect2d(x).ceil());}
		void  r_r_ceil_sse2				(VectP &r, VectP const &x)					{assign(r, Vect2d(x).ceil_sse2());}
		void  c_c_ceil					(CompP &r, CompP const &x)					{assign(r, Comp2d(x).ceil());}
		void  c_c_ceil_sse2				(CompP &r, CompP const &x)					{assign(r, Comp2d(x).ceil_sse2());}
		void  q_q_ceil					(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x).ceil());}
		void  q_q_ceil_sse2				(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x).ceil_sse2());}

		void  r_r_floor					(VectP &r, VectP const &x)					{assign(r, Vect2d(x).floor());}
		void  r_r_floor_sse2			(VectP &r, VectP const &x)					{assign(r, Vect2d(x).floor_sse2());}
		void  c_c_floor					(CompP &r, CompP const &x)					{assign(r, Comp2d(x).floor());}
		void  c_c_floor_sse2			(CompP &r, CompP const &x)					{assign(r, Comp2d(x).floor_sse2());}
		void  q_q_floor					(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x).floor());}
		void  q_q_floor_sse2			(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x).floor_sse2());}

		void  r_r_round					(VectP &r, VectP const &x)					{assign(r, Vect2d(x).round());}
		void  r_r_round_sse2			(VectP &r, VectP const &x)					{assign(r, Vect2d(x).round_sse2());}
		void  c_c_round					(CompP &r, CompP const &x)					{assign(r, Comp2d(x).round());}
		void  c_c_round_sse2			(CompP &r, CompP const &x)					{assign(r, Comp2d(x).round_sse2());}
		void  q_q_round					(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x).round());}
		void  q_q_round_sse2			(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x).round_sse2());}

		void  r_r_abs					(VectP &r, VectP const &x)					{assign(r, Vect2d(x).abs());}
		void  r_c_abs					(VectP &r, CompP const &x)					{assign(r, Comp2d(x).abs());}
		void  r_q_abs					(VectP &r, QuatP const &x)					{assign(r, Quat2d(x).abs());}

		void  r_r_arg					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, m_pi&(rx<m_zero)|m_qnan&(rx==m_zero));}
		void  r_c_arg					(VectP &r, CompP const &x)
		{
			Comp2d cx=x;
			Vect2d mask=(cx.r==m_zero)&(cx.i==m_zero);//zero: FFFF, nonzero: 0
			assign(r, Vect2d(::atan2(cx.i.v, cx.r.v))&mask.complement()|m_qnan&mask);
		}
		void  r_q_arg					(VectP &r, QuatP const &x)					{Quat2d qx=x; assign(r, Vect2d(::acos((qx.r/qx.abs()).v)));}

		void  r_c_real					(VectP &r, CompP const &x)					{assign(r, Vect2d(x.r));}

		void  r_c_imag					(VectP &r, CompP const &x)					{assign(r, Vect2d(x.i));}

		//r_conjugate: assign
		void c_c_conjugate				(CompP &r, CompP const &x)					{assign(r, Comp2d(Vect2d(x.r), -Vect2d(x.i)));}
		void q_q_conjugate				(QuatP &r, QuatP const &x)					{assign(r, Quat2d(Vect2d(x.r), -Vect2d(x.i), -Vect2d(x.j), -Vect2d(x.k)));}

		void  c_r_polar					(CompP &r, VectP const &x)					{Vect2d rx=x; assign(r, Comp2d(rx&m_sign_mask, m_pi&(rx<m_zero)|m_qnan&(rx==m_zero)));}
		void  c_c_polar					(CompP &r, CompP const &x)
		{
			Comp2d cx=x;
			Vect2d mag=cx.abs(), z_mask=mag==m_zero, arg=::atan2(cx.i.v, cx.r.v);
			assign(r, Comp2d(mag, arg&z_mask|m_qnan&z_mask.complement()));
		}
		void  c_q_polar					(CompP &r, QuatP const &x)
		{
			Quat2d qx=x;
			Vect2d mag=qx.abs();
			assign(r, Comp2d(mag, Vect2d(::acos((qx.r/mag).v))));
		}

		//r_cartesian	assign
		void  c_c_cartesian				(CompP &r, CompP const &x)
		{
			Comp2d cx=x;
			Vec2d sin_i, cos_i;
			sin_i=sincos(&cos_i, cx.i.v);
			assign(r, Comp2d(cx.r*Vect2d(cos_i), cx.r*Vect2d(sin_i)));
		}
		void  q_q_cartesian				(QuatP &r, QuatP const &x)
		{
			Quat2d qx=x;
			Vec2d sin_i, cos_i, sin_j, cos_j, sin_k, cos_k;
			sin_i=sincos(&cos_i, qx.i.v);
			sin_j=sincos(&cos_j, qx.j.v);
			sin_k=sincos(&cos_k, qx.k.v);
			cos_k*=qx.r.v;
			assign(r, Quat2d(
				Vect2d(cos_i)*Vect2d(cos_j)*Vect2d(cos_k),
				Vect2d(sin_i)*Vect2d(cos_j)*Vect2d(cos_k),
				Vect2d(sin_j)*Vect2d(cos_k),
				qx.r*Vect2d(sin_k)));
		}

		void r_rr_plus					(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x)+Vect2d(y));}
		void c_rc_plus					(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)+Comp2d(y));}
		void q_rq_plus					(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)+Quat2d(y));}
		void c_cr_plus					(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)+Vect2d(y));}
		void c_cc_plus					(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)+Comp2d(y));}
		void q_cq_plus					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)+Quat2d(y));}
		void q_qr_plus					(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)+Vect2d(y));}
		void q_qc_plus					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)+Comp2d(y));}
		void q_qq_plus					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)+Quat2d(y));}

		void  r_r_minus					(VectP &r, VectP const &x)					{assign(r, -Vect2d(x));}
		void  c_c_minus					(CompP &r, CompP const &x)					{assign(r, -Comp2d(x));}
		void  q_q_minus					(QuatP &r, QuatP const &x)					{assign(r, -Quat2d(x));}
		void r_rr_minus					(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x)-Vect2d(y));}
		void c_rc_minus					(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)-Comp2d(y));}
		void q_rq_minus					(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)-Quat2d(y));}
		void c_cr_minus					(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)-Vect2d(y));}
		void c_cc_minus					(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)-Comp2d(y));}
		void q_cq_minus					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)-Quat2d(y));}
		void q_qr_minus					(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)-Vect2d(y));}
		void q_qc_minus					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)-Comp2d(y));}
		void q_qq_minus					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)-Quat2d(y));}

		void r_rr_multiply				(VectP &r, VectP const &x, VectP const &y)
		{
			_mm_store_pd(r.r, _mm_mul_pd(_mm_load_pd(x.r), _mm_load_pd(y.r)));
		//	assign(r, Vect2d(x)*Vect2d(y));
		}
		void c_rc_multiply				(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)*Comp2d(y));}
		void q_rq_multiply				(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)*Quat2d(y));}
		void c_cr_multiply				(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)*Vect2d(y));}
		void c_cc_multiply				(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)*Comp2d(y));}//(xr+i*xi)(yr+i*yi) = xr*yr-xi*yi+i(xr*yi+xi*yr)
		void q_cq_multiply				(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)*Quat2d(y));}
		void q_qr_multiply				(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)*Vect2d(y));}
		void q_qc_multiply				(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)*Comp2d(y));}
		void q_qq_multiply				(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)*Quat2d(y));}

		void  r_r_divide				(VectP &r, VectP const &y)					{assign(r, m_one/Vect2d(y));}
		void  c_c_divide				(CompP &r, CompP const &y)					{assign(r, m_one/Comp2d(y));}
		void  q_q_divide				(QuatP &r, QuatP const &y)					{assign(r, m_one/Quat2d(y));}
		void r_rr_divide				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x)/Vect2d(y));}
		void c_rc_divide				(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)/Comp2d(y));}
		void q_rq_divide				(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)/Quat2d(y));}
		void c_cr_divide				(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)/Vect2d(y));}
		void c_cc_divide				(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)/Comp2d(y));}
		void q_cq_divide				(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)/Quat2d(y));}
		void q_qr_divide				(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)/Vect2d(y));}
		void q_qc_divide				(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)/Comp2d(y));}
		void q_qq_divide				(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)/Quat2d(y));}

		void r_rr_logic_divides			(VectP &r, VectP const &y, VectP const &x)	{auto t=Vect2d(x)/Vect2d(y); assign(r, (t==t.floor())&m_one);}//rc_divides, rq_divides: applied to each component
		void r_rc_logic_divides			(VectP &r, VectP const &y, CompP const &x)	{auto t=Comp2d(x)/Vect2d(y); assign(r, (t==t.floor())&m_one);}
		void r_rq_logic_divides			(VectP &r, VectP const &y, QuatP const &x)	{auto t=Quat2d(x)/Vect2d(y); assign(r, (t==t.floor())&m_one);}
		void r_cr_logic_divides			(VectP &r, CompP const &y, VectP const &x)	{auto t=Vect2d(x)/Comp2d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&m_one);}
		void r_cc_logic_divides			(VectP &r, CompP const &y, CompP const &x)	{auto t=Comp2d(x)/Comp2d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&m_one);}
		void r_cq_logic_divides			(VectP &r, CompP const &y, QuatP const &x)	{auto t=Quat2d(x)/Comp2d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);}
		void r_qr_logic_divides			(VectP &r, QuatP const &y, VectP const &x)	{auto t=Vect2d(x)/Quat2d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);}
		void r_qc_logic_divides			(VectP &r, QuatP const &y, CompP const &x)	{auto t=Comp2d(x)/Quat2d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);}
		void r_qq_logic_divides			(VectP &r, QuatP const &y, QuatP const &x)	{auto t=Quat2d(x)/Quat2d(y); assign(r, (t.r==t.r.floor())&(t.i==t.i.floor())&(t.j==t.j.floor())&(t.k==t.k.floor())&m_one);}
		
		void r_rr_logic_divides_sse2	(VectP &r, VectP const &y, VectP const &x)	{auto t=Vect2d(x)/Vect2d(y); assign(r, (t==t.floor_sse2())&m_one);}//rc_divides, rq_divides: applied to each component
		void r_rc_logic_divides_sse2	(VectP &r, VectP const &y, CompP const &x)	{auto t=Comp2d(x)/Vect2d(y); assign(r, (t==t.floor_sse2())&m_one);}
		void r_rq_logic_divides_sse2	(VectP &r, VectP const &y, QuatP const &x)	{auto t=Quat2d(x)/Vect2d(y); assign(r, (t==t.floor_sse2())&m_one);}
		void r_cr_logic_divides_sse2	(VectP &r, CompP const &y, VectP const &x)	{auto t=Vect2d(x)/Comp2d(y); assign(r, (t.r==t.r.floor_sse2())&(t.i==t.i.floor_sse2())&m_one);}
		void r_cc_logic_divides_sse2	(VectP &r, CompP const &y, CompP const &x)	{auto t=Comp2d(x)/Comp2d(y); assign(r, (t.r==t.r.floor_sse2())&(t.i==t.i.floor_sse2())&m_one);}
		void r_cq_logic_divides_sse2	(VectP &r, CompP const &y, QuatP const &x)	{auto t=Quat2d(x)/Comp2d(y); assign(r, (t.r==t.r.floor_sse2())&(t.i==t.i.floor_sse2())&(t.j==t.j.floor_sse2())&(t.k==t.k.floor_sse2())&m_one);}
		void r_qr_logic_divides_sse2	(VectP &r, QuatP const &y, VectP const &x)	{auto t=Vect2d(x)/Quat2d(y); assign(r, (t.r==t.r.floor_sse2())&(t.i==t.i.floor_sse2())&(t.j==t.j.floor_sse2())&(t.k==t.k.floor_sse2())&m_one);}
		void r_qc_logic_divides_sse2	(VectP &r, QuatP const &y, CompP const &x)	{auto t=Comp2d(x)/Quat2d(y); assign(r, (t.r==t.r.floor_sse2())&(t.i==t.i.floor_sse2())&(t.j==t.j.floor_sse2())&(t.k==t.k.floor_sse2())&m_one);}
		void r_qq_logic_divides_sse2	(VectP &r, QuatP const &y, QuatP const &x)	{auto t=Quat2d(x)/Quat2d(y); assign(r, (t.r==t.r.floor_sse2())&(t.i==t.i.floor_sse2())&(t.j==t.j.floor_sse2())&(t.k==t.k.floor_sse2())&m_one);}

		void c_cr_pow					(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)^Vect2d(y));}
		void c_cc_pow					(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)^Comp2d(y));}
		void q_cq_pow					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)^Quat2d(y));}
		void q_qr_pow					(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)^Vect2d(y));}
		void q_qc_pow					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)^Comp2d(y));}
		void q_qq_pow					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)^Quat2d(y));}

		void  c_c_ln					(CompP &r, CompP const &x)					{assign(r, log(Comp2d(x)));}
		void  q_q_ln					(QuatP &r, QuatP const &x)					{assign(r, log(Quat2d(x)));}
	
		const Vect2d m_inv_ln10=m_one/m_ln10;
		void  c_c_log					(CompP &r, CompP const &x)					{assign(r, log(Comp2d(x))*m_inv_ln10);}
		void  q_q_log					(QuatP &r, QuatP const &x)					{assign(r, log(Quat2d(x))*m_inv_ln10);}
		void c_cr_log					(CompP &r, CompP const &x, VectP const &y)	{assign(r, log(Comp2d(x))/log(Comp2d(Vect2d(y), m_zero)));}
		void c_cc_log					(CompP &r, CompP const &x, CompP const &y)	{assign(r, log(Comp2d(x))/log(Comp2d(y)));}
		void q_cq_log					(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, log(Comp2d(x))/log(Quat2d(y)));}
		void q_qc_log					(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, log(Quat2d(x))/log(Comp2d(y)));}
		void q_qq_log					(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, log(Quat2d(x))/log(Quat2d(y)));}
	
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
		Comp2d c_rr_tetrate					(Vect2d const &x, Vect2d const &y)
		{
			Complex1d
				lo=tetrate(x.lo(), y.lo()),
				hi=tetrate(x.hi(), y.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_rc_tetrate					(Vect2d const &xr, Comp2d const &y)
		{
			Complex1d
				lo=tetrate(xr.lo(), Complex1d(y.r.lo(), y.i.lo())),
				hi=tetrate(xr.hi(), Complex1d(y.r.hi(), y.i.hi()));
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_cr_tetrate					(Comp2d const &x, Vect2d const &y)
		{
			Complex1d
				lo=tetrate(Complex1d(x.r.lo(), x.i.lo()), y.lo()),
				hi=tetrate(Complex1d(x.r.hi(), x.i.hi()), y.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_cc_tetrate					(Comp2d const &x, Comp2d const &y)
		{
			Complex1d
				lo=tetrate(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.r.lo(), y.i.lo())),
				hi=tetrate(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.r.hi(), y.i.hi()));
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Quat2d q_qr_tetrate					(Quat2d const &x, Vect2d const &y)
		{
			Quaternion1d
				lo=tetrate(Quaternion1d(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo()), y.lo()),
				hi=tetrate(Quaternion1d(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi()), y.hi());
			return Quat2d(Vect2d(lo.R_component_1(), hi.R_component_1()), Vect2d(lo.R_component_2(), hi.R_component_2()),
				Vect2d(lo.R_component_3(), hi.R_component_3()), Vect2d(lo.R_component_4(), hi.R_component_4()));
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
		Comp2d c_rr_pentate					(Vect2d const &x, Vect2d const &y)
		{
			Complex1d
				lo=pentate(x.lo(), y.lo()),
				hi=pentate(x.hi(), y.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_cr_pentate					(Comp2d const &x, Vect2d const &y)
		{
			Complex1d
				lo=pentate(Complex1d(x.r.lo(), x.i.lo()), y.lo()),
				hi=pentate(Complex1d(x.r.hi(), x.i.hi()), y.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		bool disc_rr_pentate_i			(Value const &x0, Value const &y0, Value const &x1, Value const &y1){return false;}//
		bool disc_cr_pentate_i			(Value const &x0, Value const &y0, Value const &x1, Value const &y1){return false;}////*/

		void  r_r_bitwise_shift_left_l	(VectP &r, VectP const &x)					{assign(r, Vect2d(::exp((Vect2d(x).floor()*m_ln2).v)));}//<<x = 2^x
		void  c_c_bitwise_shift_left_l	(CompP &r, CompP const &x)					{assign(r, exp(Comp2d(x).floor()*m_ln2));}
		void  q_q_bitwise_shift_left_l	(QuatP &r, QuatP const &x)					{assign(r, exp(Quat2d(x).floor()*m_ln2));}
		void  r_r_bitwise_shift_left_r	(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, rx+rx);}//x<< = x*2
		void  c_c_bitwise_shift_left_r	(CompP &r, CompP const &x)					{Comp2d cx=x; assign(r, cx+cx);}
		void  q_q_bitwise_shift_left_r	(QuatP &r, QuatP const &x)					{Quat2d qx=x; assign(r, qx+qx);}
		void r_rr_bitwise_shift_left	(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x)*Vect2d(::exp((Vect2d(y).floor()*m_ln2).v)));}//x<<y = x*2^y
		void c_rc_bitwise_shift_left	(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)*exp(Comp2d(y).floor()*m_ln2));}
		void q_rq_bitwise_shift_left	(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)*exp(Quat2d(y).floor()*m_ln2));}
		void c_cr_bitwise_shift_left	(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)*Vect2d(::exp((Vect2d(y).floor()*m_ln2).v)));}
		void c_cc_bitwise_shift_left	(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)*exp(Comp2d(y).floor()*m_ln2));}
		void q_cq_bitwise_shift_left	(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)*exp(Quat2d(y).floor()*m_ln2));}
		void q_qr_bitwise_shift_left	(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)*Vect2d(::exp((Vect2d(y).floor()*m_ln2).v)));}
		void q_qc_bitwise_shift_left	(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)*exp(Comp2d(y).floor()*m_ln2));}
		void q_qq_bitwise_shift_left	(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)*exp(Quat2d(y).floor()*m_ln2));}

		void  r_r_bitwise_shift_right_l	(VectP &r, VectP const &x)					{assign(r, Vect2d(::exp((-Vect2d(x).floor()*m_ln2).v)));}//>>x = 2^-x = exp(-x*ln2)
		void  c_c_bitwise_shift_right_l	(CompP &r, CompP const &x)					{assign(r, exp(-Comp2d(x).floor()*m_ln2));}
		void  q_q_bitwise_shift_right_l	(QuatP &r, QuatP const &x)					{assign(r, exp(-Quat2d(x).floor()*m_ln2));}
		void  r_r_bitwise_shift_right_r	(VectP &r, VectP const &x)					{assign(r, Vect2d(x)*m_half);}//x>> = x/2
		void  c_c_bitwise_shift_right_r	(CompP &r, CompP const &x)					{assign(r, Comp2d(x)*m_half);}
		void  q_q_bitwise_shift_right_r	(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x)*m_half);}
		void r_rr_bitwise_shift_right	(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x)*Vect2d(::exp((-Vect2d(y).floor()*m_ln2).v)));}
		void c_rc_bitwise_shift_right	(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)*exp(-Comp2d(y).floor()*m_ln2));}
		void q_rq_bitwise_shift_right	(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)*exp(-Quat2d(y).floor()*m_ln2));}
		void c_cr_bitwise_shift_right	(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)*Vect2d(::exp((-Vect2d(y).floor()*m_ln2).v)));}
		void c_cc_bitwise_shift_right	(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)*exp(-Comp2d(y).floor()*m_ln2));}
		void q_cq_bitwise_shift_right	(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)*exp(-Quat2d(y).floor()*m_ln2));}
		void q_qr_bitwise_shift_right	(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)*Vect2d(::exp((-Vect2d(y).floor()*m_ln2).v)));}
		void q_qc_bitwise_shift_right	(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)*exp(-Comp2d(y).floor()*m_ln2));}
		void q_qq_bitwise_shift_right	(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)*exp(-Quat2d(y).floor()*m_ln2));}

	/*	__forceinline Vect2d bitwise_not(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			Vect2d xc((double)~convert_d2ll(x.lo()), (double)~convert_d2ll(x.hi()));
			return xc&mask;
		}
		Vect2d  r_r_bitwise_not				(Vect2d const &x)					{return bitwise_not(x);}
		Comp2d  c_c_bitwise_not				(Comp2d const &x)					{return Comp2d(bitwise_not(x.r), bitwise_not(x.i));}
		Quat2d  q_q_bitwise_not				(Quat2d const &x)					{return Quat2d(bitwise_not(x.r), bitwise_not(x.i), bitwise_not(x.j), bitwise_not(x.k));}

		__forceinline Vect2d bitwise_and(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			return Vect2d((double)!~convert_d2ll(x.lo()), (double)!~convert_d2ll(x.hi()))&mask;
		}
		__forceinline Vect2d bitwise_and(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect2d(double(convert_d2ll(x.lo())&convert_d2ll(y.lo())), double(convert_d2ll(x.hi())&convert_d2ll(y.hi())))&mask;
		}
		__forceinline Vect2d bitwise_and_ll(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect2d v;
			(long long&)v.lo()=convert_d2ll(x.lo())&convert_d2ll(y.lo());
			(long long&)v.hi()=convert_d2ll(x.hi())&convert_d2ll(y.hi());
			return v&mask;
		}
		__forceinline Vect2d convert_ll2d(Vect2d const &x){Vect2d x2=x; return Vect2d((double)(long long&)x2.lo(), (double)(long long&)x2.hi());}
		Vect2d  r_r_bitwise_and				(Vect2d const &x)					{return bitwise_and(x);}
		Comp2d  c_c_bitwise_and				(Comp2d const &x)					{return Comp2d(bitwise_and(x.r), bitwise_and(x.i));}
		Quat2d  q_q_bitwise_and				(Quat2d const &x)					{return Quat2d(bitwise_and(x.r), bitwise_and(x.i), bitwise_and(x.j), bitwise_and(x.k));}
		Vect2d r_rr_bitwise_and				(Vect2d const &x, Vect2d const &y)	{return bitwise_and(x, y);}
		Comp2d c_rc_bitwise_and				(Vect2d const &x, Comp2d const &y)	{return Comp2d(bitwise_and(x, y.r), bitwise_and(x, y.i));}
		Quat2d q_rq_bitwise_and				(Vect2d const &x, Quat2d const &y)	{return Quat2d(bitwise_and(x, y.r), bitwise_and(x, y.i), bitwise_and(x, y.j), bitwise_and(x, y.k));}
		Comp2d c_cr_bitwise_and				(Comp2d const &x, Vect2d const &y)	{return Comp2d(bitwise_and(x.r, y), bitwise_and(x.i, y));}
		Comp2d c_cc_bitwise_and				(Comp2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i);
			return Comp2d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr));
		}
		Quat2d q_cq_bitwise_and				(Comp2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i),
				xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j),
				xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k);
			return Quat2d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xr_yj-xi_yk), convert_ll2d(xr_yk+xi_yj));
		}
		Quat2d q_qr_bitwise_and				(Quat2d const &x, Vect2d const &y)	{return Quat2d(bitwise_and(x.r, y), bitwise_and(x.i, y), bitwise_and(x.j, y), bitwise_and(x.k, y));}
		Quat2d q_qc_bitwise_and				(Quat2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i);
			return Quat2d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xj_yr+xk_yi), convert_ll2d(-xj_yi+xk_yr));
		}
		Quat2d q_qq_bitwise_and				(Quat2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
				xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i),
				xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j), xj_yj=bitwise_and_ll(x.j, y.j), xk_yj=bitwise_and_ll(x.k, y.j),
				xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k), xj_yk=bitwise_and_ll(x.j, y.k), xk_yk=bitwise_and_ll(x.k, y.k);
			return Quat2d(convert_ll2d(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d(xr_yk+xi_yj-xj_yi+xk_yr));
		}

		__forceinline Vect2d bitwise_nand(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			return Vect2d((double)(~convert_d2ll(x.lo())!=0), (double)(~convert_d2ll(x.hi())!=0))&mask;
		}
		__forceinline Vect2d bitwise_nand(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect2d(double(~(convert_d2ll(x.lo())&convert_d2ll(y.lo()))), double(~(convert_d2ll(x.hi())&convert_d2ll(y.hi()))))&mask;
		}
		__forceinline Vect2d bitwise_nand_ll(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect2d v;
			(long long&)v.lo()=~(convert_d2ll(x.lo())&convert_d2ll(y.lo()));
			(long long&)v.hi()=~(convert_d2ll(x.hi())&convert_d2ll(y.hi()));
			return v&mask;
		}
		Vect2d  r_r_bitwise_nand			(Vect2d const &x)					{return bitwise_nand(x);}
		Comp2d  c_c_bitwise_nand			(Comp2d const &x)					{return Comp2d(bitwise_nand(x.r), bitwise_nand(x.i));}
		Quat2d  q_q_bitwise_nand			(Quat2d const &x)					{return Quat2d(bitwise_nand(x.r), bitwise_nand(x.i), bitwise_nand(x.j), bitwise_nand(x.k));}
		Vect2d r_rr_bitwise_nand			(Vect2d const &x, Vect2d const &y)	{return bitwise_nand(x, y);}
		Comp2d c_rc_bitwise_nand			(Vect2d const &x, Comp2d const &y)	{return Comp2d(bitwise_nand(x, y.r), bitwise_nand(x, y.i));}
		Quat2d q_rq_bitwise_nand			(Vect2d const &x, Quat2d const &y)	{return Quat2d(bitwise_nand(x, y.r), bitwise_nand(x, y.i), bitwise_nand(x, y.j), bitwise_nand(x, y.k));}
		Comp2d c_cr_bitwise_nand			(Comp2d const &x, Vect2d const &y)	{return Comp2d(bitwise_nand(x.r, y), bitwise_nand(x.i, y));}
		Comp2d c_cc_bitwise_nand			(Comp2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i);
			return Comp2d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr));
		}
		Quat2d q_cq_bitwise_nand			(Comp2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i),
				xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j),
				xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k);
			return Quat2d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xr_yj-xi_yk), convert_ll2d(xr_yk+xi_yj));
		}
		Quat2d q_qr_bitwise_nand			(Quat2d const &x, Vect2d const &y)	{return Quat2d(bitwise_nand(x.r, y), bitwise_nand(x.i, y), bitwise_nand(x.j, y), bitwise_nand(x.k, y));}
		Quat2d q_qc_bitwise_nand			(Quat2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i);
			return Quat2d(convert_ll2d(xr_yr-xi_yi), convert_ll2d(xr_yi+xi_yr), convert_ll2d(xj_yr+xk_yi), convert_ll2d(-xj_yi+xk_yr));
		}
		Quat2d q_qq_bitwise_nand			(Quat2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
				xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i),
				xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j), xj_yj=bitwise_nand_ll(x.j, y.j), xk_yj=bitwise_nand_ll(x.k, y.j),
				xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k), xj_yk=bitwise_nand_ll(x.j, y.k), xk_yk=bitwise_nand_ll(x.k, y.k);
			return Quat2d(convert_ll2d(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d(xr_yk+xi_yj-xj_yi+xk_yr));
		}
	
		__forceinline Vect2d bitwise_or(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			return Vect2d((double)(convert_d2ll(x.lo())!=0), (double)(convert_d2ll(x.hi())!=0))&mask;
		}
		__forceinline Vect2d bitwise_or(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect2d(double(convert_d2ll(x.lo())|convert_d2ll(y.lo())), double(convert_d2ll(x.hi())|convert_d2ll(y.hi())))&mask;
		}
		__forceinline Vect2d bitwise_or_ll_c(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect2d v;
			(long long&)v.lo()=~convert_d2ll(x.lo())&~convert_d2ll(y.lo());
			(long long&)v.hi()=~convert_d2ll(x.hi())&~convert_d2ll(y.hi());
			return v&mask;
		}
		__forceinline Vect2d convert_ll2d_c(Vect2d const &x){Vect2d x2=x; return Vect2d((double)~(long long)(long long&)x2.lo(), (double)~(long long)(long long&)x2.hi());}
		Vect2d  r_r_bitwise_or				(Vect2d const &x)					{return bitwise_or(x);}
		Comp2d  c_c_bitwise_or				(Comp2d const &x)					{return Comp2d(bitwise_or(x.r), bitwise_or(x.i));}
		Quat2d  q_q_bitwise_or				(Quat2d const &x)					{return Quat2d(bitwise_or(x.r), bitwise_or(x.i), bitwise_or(x.j), bitwise_or(x.k));}
		Vect2d r_rr_bitwise_or				(Vect2d const &x, Vect2d const &y)	{return bitwise_or(x, y);}
		Comp2d c_rc_bitwise_or				(Vect2d const &x, Comp2d const &y)	{return Comp2d(bitwise_or(x, y.r), bitwise_or(x, y.i));}
		Quat2d q_rq_bitwise_or				(Vect2d const &x, Quat2d const &y)	{return Quat2d(bitwise_or(x, y.r), bitwise_or(x, y.i), bitwise_or(x, y.j), bitwise_or(x, y.k));}
		Comp2d c_cr_bitwise_or				(Comp2d const &x, Vect2d const &y)	{return Comp2d(bitwise_or(x.r, y), bitwise_or(x.i, y));}
		Comp2d c_cc_bitwise_or				(Comp2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i);
			return Comp2d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr));
		}
		Quat2d q_cq_bitwise_or				(Comp2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i),
				xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j),
				xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k);
			return Quat2d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xr_yj-xi_yk), convert_ll2d_c(xr_yk+xi_yj));
		}
		Quat2d q_qr_bitwise_or				(Quat2d const &x, Vect2d const &y)	{return Quat2d(bitwise_or(x.r, y), bitwise_or(x.i, y), bitwise_or(x.j, y), bitwise_or(x.k, y));}
		Quat2d q_qc_bitwise_or				(Quat2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i);
			return Quat2d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xj_yr+xk_yi), convert_ll2d_c(-xj_yi+xk_yr));
		}
		Quat2d q_qq_bitwise_or				(Quat2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
				xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i),
				xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j), xj_yj=bitwise_or_ll_c(x.j, y.j), xk_yj=bitwise_or_ll_c(x.k, y.j),
				xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k), xj_yk=bitwise_or_ll_c(x.j, y.k), xk_yk=bitwise_or_ll_c(x.k, y.k);
			return Quat2d(convert_ll2d_c(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d_c(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d_c(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d_c(xr_yk+xi_yj-xj_yi+xk_yr));
		}
	
		__forceinline Vect2d bitwise_nor(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			return Vect2d((double)!convert_d2ll(x.lo()), (double)!convert_d2ll(x.hi()))&mask;
		}
		__forceinline Vect2d bitwise_nor(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			return Vect2d(double(~(convert_d2ll(x.lo())|convert_d2ll(y.lo()))), double(~(convert_d2ll(x.hi())|convert_d2ll(y.hi()))))&mask;
		}
		__forceinline Vect2d bitwise_nor_ll_c(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()|is_not_nan(y)&isinf(y).complement();
			Vect2d v;
			(long long&)v.lo()=~convert_d2ll(x.lo())|~convert_d2ll(y.lo());
			(long long&)v.hi()=~convert_d2ll(x.hi())|~convert_d2ll(y.hi());
			return v&mask;
		}
		Vect2d  r_r_bitwise_nor				(Vect2d const &x)					{return bitwise_nor(x);}
		Comp2d  c_c_bitwise_nor				(Comp2d const &x)					{return Comp2d(bitwise_nor(x.r), bitwise_nor(x.i));}
		Quat2d  q_q_bitwise_nor				(Quat2d const &x)					{return Quat2d(bitwise_nor(x.r), bitwise_nor(x.i), bitwise_nor(x.j), bitwise_nor(x.k));}
		Vect2d r_rr_bitwise_nor				(Vect2d const &x, Vect2d const &y)	{return bitwise_nor(x, y);}
		Comp2d c_rc_bitwise_nor				(Vect2d const &x, Comp2d const &y)	{return Comp2d(bitwise_nor(x, y.r), bitwise_nor(x, y.i));}
		Quat2d q_rq_bitwise_nor				(Vect2d const &x, Quat2d const &y)	{return Quat2d(bitwise_nor(x, y.r), bitwise_nor(x, y.i), bitwise_nor(x, y.j), bitwise_nor(x, y.k));}
		Comp2d c_cr_bitwise_nor				(Comp2d const &x, Vect2d const &y)	{return Comp2d(bitwise_nor(x.r, y), bitwise_nor(x.i, y));}
		Comp2d c_cc_bitwise_nor				(Comp2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i);
			return Comp2d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr));
		}
		Quat2d q_cq_bitwise_nor				(Comp2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i),
				xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j),
				xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k);
			return Quat2d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xr_yj-xi_yk), convert_ll2d_c(xr_yk+xi_yj));
		}
		Quat2d q_qr_bitwise_nor				(Quat2d const &x, Vect2d const &y)	{return Quat2d(bitwise_nor(x.r, y), bitwise_nor(x.i, y), bitwise_nor(x.j, y), bitwise_nor(x.k, y));}
		Quat2d q_qc_bitwise_nor				(Quat2d const &x, Comp2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i);
			return Quat2d(convert_ll2d_c(xr_yr-xi_yi), convert_ll2d_c(xr_yi+xi_yr), convert_ll2d_c(xj_yr+xk_yi), convert_ll2d_c(-xj_yi+xk_yr));
		}
		Quat2d q_qq_bitwise_nor				(Quat2d const &x, Quat2d const &y)
		{
			Vect2d
				xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
				xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i),
				xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j), xj_yj=bitwise_nor_ll_c(x.j, y.j), xk_yj=bitwise_nor_ll_c(x.k, y.j),
				xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k), xj_yk=bitwise_nor_ll_c(x.j, y.k), xk_yk=bitwise_nor_ll_c(x.k, y.k);
			return Quat2d(convert_ll2d_c(xr_yr-xi_yi-xj_yj-xk_yk), convert_ll2d_c(xr_yi+xi_yr+xj_yk-xk_yj), convert_ll2d_c(xj_yj-xi_yk+xj_yr+xk_yi), convert_ll2d_c(xr_yk+xi_yj-xj_yi+xk_yr));
		}
	
		__forceinline Vect2d bitwise_xor(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			return Vect2d((double)G2::bitwise_xor(x.lo()), (double)bitwise_xor(x.hi()))&mask;
		}
		__forceinline Vect2d bitwise_xor(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()&is_not_nan(y)&isinf(y).complement();
			return Vect2d(double(convert_d2ll(x.lo())^convert_d2ll(y.lo())), double(convert_d2ll(x.hi())^convert_d2ll(y.hi())))&mask;
		}
		Vect2d  r_r_bitwise_xor				(Vect2d const &x)					{return bitwise_xor(x);}
		Comp2d  c_c_bitwise_xor				(Comp2d const &x)					{return Comp2d(bitwise_xor(x.r), bitwise_xor(x.i));}
		Quat2d  q_q_bitwise_xor				(Quat2d const &x)					{return Quat2d(bitwise_xor(x.r), bitwise_xor(x.i), bitwise_xor(x.j), bitwise_xor(x.k));}
		Vect2d r_rr_bitwise_xor				(Vect2d const &x, Vect2d const &y)	{return bitwise_xor(x, y);}
		Comp2d c_rc_bitwise_xor				(Vect2d const &x, Comp2d const &y)	{return Comp2d(bitwise_xor(x, y.r), y.i);}
		Quat2d q_rq_bitwise_xor				(Vect2d const &x, Quat2d const &y)	{return Quat2d(bitwise_xor(x, y.r), y.i, y.j, y.k);}
		Comp2d c_cr_bitwise_xor				(Comp2d const &x, Vect2d const &y)	{return Comp2d(bitwise_xor(x.r, y), x.i);}
		Comp2d c_cc_bitwise_xor				(Comp2d const &x, Comp2d const &y)	{return Comp2d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i));}
		Quat2d q_cq_bitwise_xor				(Comp2d const &x, Quat2d const &y)	{return Quat2d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), y.j, y.k);}
		Quat2d q_qr_bitwise_xor				(Quat2d const &x, Vect2d const &y)	{return Quat2d(bitwise_xor(x.r, y), x.i, x.j, x.k);}
		Quat2d q_qc_bitwise_xor				(Quat2d const &x, Comp2d const &y)	{return Quat2d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), x.j, x.k);}
		Quat2d q_qq_bitwise_xor				(Quat2d const &x, Quat2d const &y)	{return Quat2d(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), bitwise_xor(x.j, y.j), bitwise_xor(x.k, y.k));}
	
		__forceinline Vect2d bitwise_xnor(Vect2d const &x)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement();
			return Vect2d((double)!bitwise_xor(x.lo()), (double)!bitwise_xor(x.hi()))&mask;
		}
		__forceinline Vect2d bitwise_xnor(Vect2d const &x, Vect2d const &y)
		{
			Vect2d mask=is_not_nan(x)&isinf(x).complement()&is_not_nan(y)&isinf(y).complement();
			return Vect2d((double)~(convert_d2ll(x.lo())^convert_d2ll(y.lo())), (double)~(convert_d2ll(x.hi())^convert_d2ll(y.hi())))&mask;
		}
		Vect2d  r_r_bitwise_xnor			(Vect2d const &x)					{return bitwise_xnor(x);}
		Comp2d  c_c_bitwise_xnor			(Comp2d const &x)					{return Comp2d(bitwise_xnor(x.r), bitwise_xnor(x.i));}
		Quat2d  q_q_bitwise_xnor			(Quat2d const &x)					{return Quat2d(bitwise_xnor(x.r), bitwise_xnor(x.i), bitwise_xnor(x.j), bitwise_xnor(x.k));}
		Vect2d r_rr_bitwise_xnor			(Vect2d const &x, Vect2d const &y)	{return bitwise_xnor(x, y);}
		Comp2d c_rc_bitwise_xnor			(Vect2d const &x, Comp2d const &y)	{return Comp2d(bitwise_xnor(x, y.r), y.i);}
		Quat2d q_rq_bitwise_xnor			(Vect2d const &x, Quat2d const &y)	{return Quat2d(bitwise_xnor(x, y.r), y.i, y.j, y.k);}
		Comp2d c_cr_bitwise_xnor			(Comp2d const &x, Vect2d const &y)	{return Comp2d(bitwise_xnor(x.r, y), x.i);}
		Comp2d c_cc_bitwise_xnor			(Comp2d const &x, Comp2d const &y)	{return Comp2d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i));}
		Quat2d q_cq_bitwise_xnor			(Comp2d const &x, Quat2d const &y)	{return Quat2d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), y.j, y.k);}
		Quat2d q_qr_bitwise_xnor			(Quat2d const &x, Vect2d const &y)	{return Quat2d(bitwise_xnor(x.r, y), x.i, x.j, x.k);}
		Quat2d q_qc_bitwise_xnor			(Quat2d const &x, Comp2d const &y)	{return Quat2d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), x.j, x.k);}
		Quat2d q_qq_bitwise_xnor			(Quat2d const &x, Quat2d const &y)	{return Quat2d(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), bitwise_xnor(x.j, y.j), bitwise_xnor(x.k, y.k));}//*/
	
		void  r_r_logic_equal			(VectP &r, VectP const &x)					{assign(r, (Vect2d(x)==m_zero)&m_one);}
		void  r_c_logic_equal			(VectP &r, CompP const &x)					{assign(r, Comp2d(x).c_is_true().complement()&m_one);}
		void  r_q_logic_equal			(VectP &r, QuatP const &x)					{assign(r, Quat2d(x).q_is_true().complement()&m_one);}
		void r_rr_logic_equal			(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x)==Vect2d(y))&m_one);}
		void r_rc_logic_equal			(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x)==Comp2d(y))&m_one);}
		void r_rq_logic_equal			(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x)==Quat2d(y))&m_one);}
		void r_cr_logic_equal			(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp2d(x)==Vect2d(y))&m_one);}
		void r_cc_logic_equal			(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp2d(x)==Comp2d(y))&m_one);}
		void r_cq_logic_equal			(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp2d(x)==Quat2d(y))&m_one);}
		void r_qr_logic_equal			(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat2d(x)==Vect2d(y))&m_one);}
		void r_qc_logic_equal			(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat2d(x)==Comp2d(y))&m_one);}
		void r_qq_logic_equal			(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat2d(x)==Quat2d(y))&m_one);}
	
		void  r_r_logic_not_equal		(VectP &r, VectP const &x)					{assign(r, Vect2d(x).r_is_true()&m_one);}
		void  r_c_logic_not_equal		(VectP &r, CompP const &x)					{assign(r, Comp2d(x).c_is_true()&m_one);}
		void  r_q_logic_not_equal		(VectP &r, QuatP const &x)					{assign(r, Quat2d(x).q_is_true()&m_one);}
		void r_rr_logic_not_equal		(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x)!=Vect2d(y))&m_one);}
		void r_rc_logic_not_equal		(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x)!=Comp2d(y))&m_one);}
		void r_rq_logic_not_equal		(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x)!=Quat2d(y))&m_one);}
		void r_cr_logic_not_equal		(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp2d(x)!=Vect2d(y))&m_one);}
		void r_cc_logic_not_equal		(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp2d(x)!=Comp2d(y))&m_one);}
		void r_cq_logic_not_equal		(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp2d(x)!=Quat2d(y))&m_one);}
		void r_qr_logic_not_equal		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat2d(x)!=Vect2d(y))&m_one);}
		void r_qc_logic_not_equal		(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat2d(x)!=Comp2d(y))&m_one);}
		void r_qq_logic_not_equal		(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat2d(x)!=Quat2d(y))&m_one);}
	
		void  r_r_logic_less_l			(VectP &r, VectP const &x)					{assign(r, (m_zero<Vect2d(x))&m_one);}
		void  r_c_logic_less_l			(VectP &r, CompP const &x)					{assign(r, (m_zero<Vect2d(x.r))&m_one);}
		void  r_q_logic_less_l			(VectP &r, QuatP const &x)					{assign(r, (m_zero<Vect2d(x.r))&m_one);}
		void  r_r_logic_less_r			(VectP &r, VectP const &x)					{assign(r, (Vect2d(x)<m_zero)&m_one);}
		void  r_c_logic_less_r			(VectP &r, CompP const &x)					{assign(r, (Vect2d(x.r)<m_zero)&m_one);}
		void  r_q_logic_less_r			(VectP &r, QuatP const &x)					{assign(r, (Vect2d(x.r)<m_zero)&m_one);}
		void r_rr_logic_less			(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x)<Vect2d(y))&m_one);}
		void r_rc_logic_less			(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x)<Vect2d(y.r))&m_one);}
		void r_rq_logic_less			(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x)<Vect2d(y.r))&m_one);}
		void r_cr_logic_less			(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)<Vect2d(y))&m_one);}
		void r_cc_logic_less			(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)<Vect2d(y.r))&m_one);}
		void r_cq_logic_less			(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)<Vect2d(y.r))&m_one);}
		void r_qr_logic_less			(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)<Vect2d(y))&m_one);}
		void r_qc_logic_less			(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)<Vect2d(y.r))&m_one);}
		void r_qq_logic_less			(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)<Vect2d(y.r))&m_one);}
	
		void  r_r_logic_less_equal_l	(VectP &r, VectP const &x)					{assign(r, (m_zero<=Vect2d(x))&m_one);}
		void  r_c_logic_less_equal_l	(VectP &r, CompP const &x)					{assign(r, (m_zero<=Vect2d(x.r))&m_one);}
		void  r_q_logic_less_equal_l	(VectP &r, QuatP const &x)					{assign(r, (m_zero<=Vect2d(x.r))&m_one);}
		void  r_r_logic_less_equal_r	(VectP &r, VectP const &x)					{assign(r, (Vect2d(x)<=m_zero)&m_one);}
		void  r_c_logic_less_equal_r	(VectP &r, CompP const &x)					{assign(r, (Vect2d(x.r)<=m_zero)&m_one);}
		void  r_q_logic_less_equal_r	(VectP &r, QuatP const &x)					{assign(r, (Vect2d(x.r)<=m_zero)&m_one);}
		void r_rr_logic_less_equal		(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x)<=Vect2d(y))&m_one);}
		void r_rc_logic_less_equal		(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x)<=Vect2d(y.r))&m_one);}
		void r_rq_logic_less_equal		(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x)<=Vect2d(y.r))&m_one);}
		void r_cr_logic_less_equal		(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)<=Vect2d(y))&m_one);}
		void r_cc_logic_less_equal		(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)<=Vect2d(y.r))&m_one);}
		void r_cq_logic_less_equal		(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)<=Vect2d(y.r))&m_one);}
		void r_qr_logic_less_equal		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)<=Vect2d(y))&m_one);}
		void r_qc_logic_less_equal		(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)<=Vect2d(y.r))&m_one);}
		void r_qq_logic_less_equal		(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)<=Vect2d(y.r))&m_one);}
	
		void  r_r_logic_greater_l		(VectP &r, VectP const &x)					{assign(r, (m_zero>Vect2d(x))&m_one);}
		void  r_c_logic_greater_l		(VectP &r, CompP const &x)					{assign(r, (m_zero>Vect2d(x.r))&m_one);}
		void  r_q_logic_greater_l		(VectP &r, QuatP const &x)					{assign(r, (m_zero>Vect2d(x.r))&m_one);}
		void  r_r_logic_greater_r		(VectP &r, VectP const &x)					{assign(r, (Vect2d(x)>m_zero)&m_one);}
		void  r_c_logic_greater_r		(VectP &r, CompP const &x)					{assign(r, (Vect2d(x.r)>m_zero)&m_one);}
		void  r_q_logic_greater_r		(VectP &r, QuatP const &x)					{assign(r, (Vect2d(x.r)>m_zero)&m_one);}
		void r_rr_logic_greater			(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x)>Vect2d(y))&m_one);}
		void r_rc_logic_greater			(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x)>Vect2d(y.r))&m_one);}
		void r_rq_logic_greater			(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x)>Vect2d(y.r))&m_one);}
		void r_cr_logic_greater			(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)>Vect2d(y))&m_one);}
		void r_cc_logic_greater			(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)>Vect2d(y.r))&m_one);}
		void r_cq_logic_greater			(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)>Vect2d(y.r))&m_one);}
		void r_qr_logic_greater			(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)>Vect2d(y))&m_one);}
		void r_qc_logic_greater			(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)>Vect2d(y.r))&m_one);}
		void r_qq_logic_greater			(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)>Vect2d(y.r))&m_one);}
	
		void  r_r_logic_greater_equal_l	(VectP &r, VectP const &x)					{assign(r, (m_zero>=Vect2d(x))&m_one);}
		void  r_c_logic_greater_equal_l	(VectP &r, CompP const &x)					{assign(r, (m_zero>=Vect2d(x.r))&m_one);}
		void  r_q_logic_greater_equal_l	(VectP &r, QuatP const &x)					{assign(r, (m_zero>=Vect2d(x.r))&m_one);}
		void  r_r_logic_greater_equal_r	(VectP &r, VectP const &x)					{assign(r, (Vect2d(x)>=m_zero)&m_one);}
		void  r_c_logic_greater_equal_r	(VectP &r, CompP const &x)					{assign(r, (Vect2d(x.r)>=m_zero)&m_one);}
		void  r_q_logic_greater_equal_r	(VectP &r, QuatP const &x)					{assign(r, (Vect2d(x.r)>=m_zero)&m_one);}
		void r_rr_logic_greater_equal	(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x)>=Vect2d(y))&m_one);}
		void r_rc_logic_greater_equal	(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x)>=Vect2d(y.r))&m_one);}
		void r_rq_logic_greater_equal	(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x)>=Vect2d(y.r))&m_one);}
		void r_cr_logic_greater_equal	(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)>=Vect2d(y))&m_one);}
		void r_cc_logic_greater_equal	(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)>=Vect2d(y.r))&m_one);}
		void r_cq_logic_greater_equal	(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)>=Vect2d(y.r))&m_one);}
		void r_qr_logic_greater_equal	(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Vect2d(x.r)>=Vect2d(y))&m_one);}
		void r_qc_logic_greater_equal	(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Vect2d(x.r)>=Vect2d(y.r))&m_one);}
		void r_qq_logic_greater_equal	(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Vect2d(x.r)>=Vect2d(y.r))&m_one);}
	
		void  r_r_logic_not				(VectP &r, VectP const &x)					{assign(r, (Vect2d(x)==m_zero)&m_one);}
		void  r_c_logic_not				(VectP &r, CompP const &x)					{Comp2d cx=x; assign(r, (cx.r==m_zero)&(cx.i==m_zero)&m_one);}
		void  r_q_logic_not				(VectP &r, QuatP const &x)					{Quat2d qx=x; assign(r, (qx.r==m_zero)&(qx.i==m_zero)&(qx.j==m_zero)&(qx.k==m_zero)&m_one);}
	
		void r_rr_logic_and				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x).r_is_true()&Vect2d(y).r_is_true()&m_one);}
		void r_rc_logic_and				(VectP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x).r_is_true()&Comp2d(y).c_is_true()&m_one);}
		void r_rq_logic_and				(VectP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x).r_is_true()&Quat2d(y).q_is_true()&m_one);}
		void r_cr_logic_and				(VectP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x).c_is_true()&Vect2d(y).r_is_true()&m_one);}
		void r_cc_logic_and				(VectP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x).c_is_true()&Comp2d(y).c_is_true()&m_one);}
		void r_cq_logic_and				(VectP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x).c_is_true()&Quat2d(y).q_is_true()&m_one);}
		void r_qr_logic_and				(VectP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x).q_is_true()&Vect2d(y).r_is_true()&m_one);}
		void r_qc_logic_and				(VectP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x).q_is_true()&Comp2d(y).c_is_true()&m_one);}
		void r_qq_logic_and				(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x).q_is_true()&Quat2d(y).q_is_true()&m_one);}
	
		void r_rr_logic_or				(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x).r_is_true()|Vect2d(y).r_is_true())&m_one);}
		void r_rc_logic_or				(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x).r_is_true()|Comp2d(y).c_is_true())&m_one);}
		void r_rq_logic_or				(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x).r_is_true()|Quat2d(y).q_is_true())&m_one);}
		void r_cr_logic_or				(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp2d(x).c_is_true()|Vect2d(y).r_is_true())&m_one);}
		void r_cc_logic_or				(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp2d(x).c_is_true()|Comp2d(y).c_is_true())&m_one);}
		void r_cq_logic_or				(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp2d(x).c_is_true()|Quat2d(y).q_is_true())&m_one);}
		void r_qr_logic_or				(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat2d(x).q_is_true()|Vect2d(y).r_is_true())&m_one);}
		void r_qc_logic_or				(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat2d(x).q_is_true()|Comp2d(y).c_is_true())&m_one);}
		void r_qq_logic_or				(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat2d(x).q_is_true()|Quat2d(y).q_is_true())&m_one);}
	
		void r_rr_logic_xor				(VectP &r, VectP const &x, VectP const &y)	{assign(r, (Vect2d(x).r_is_true()^Vect2d(y).r_is_true())&m_one);}
		void r_rc_logic_xor				(VectP &r, VectP const &x, CompP const &y)	{assign(r, (Vect2d(x).r_is_true()^Comp2d(y).c_is_true())&m_one);}
		void r_rq_logic_xor				(VectP &r, VectP const &x, QuatP const &y)	{assign(r, (Vect2d(x).r_is_true()^Quat2d(y).q_is_true())&m_one);}
		void r_cr_logic_xor				(VectP &r, CompP const &x, VectP const &y)	{assign(r, (Comp2d(x).c_is_true()^Vect2d(y).r_is_true())&m_one);}
		void r_cc_logic_xor				(VectP &r, CompP const &x, CompP const &y)	{assign(r, (Comp2d(x).c_is_true()^Comp2d(y).c_is_true())&m_one);}
		void r_cq_logic_xor				(VectP &r, CompP const &x, QuatP const &y)	{assign(r, (Comp2d(x).c_is_true()^Quat2d(y).q_is_true())&m_one);}
		void r_qr_logic_xor				(VectP &r, QuatP const &x, VectP const &y)	{assign(r, (Quat2d(x).q_is_true()^Vect2d(y).r_is_true())&m_one);}
		void r_qc_logic_xor				(VectP &r, QuatP const &x, CompP const &y)	{assign(r, (Quat2d(x).q_is_true()^Comp2d(y).c_is_true())&m_one);}
		void r_qq_logic_xor				(VectP &r, QuatP const &x, QuatP const &y)	{assign(r, (Quat2d(x).q_is_true()^Quat2d(y).q_is_true())&m_one);}
	
		void r_rr_condition_zero	(VectP &r, VectP const &x, VectP const &y)	{Vect2d rx=x, mask=rx.r_is_true();	assign(r, rx&mask|Vect2d(y)&mask.complement());}
		void c_rc_condition_zero	(CompP &r, VectP const &x, CompP const &y)	{Vect2d rx=x, mask=rx.r_is_true();	assign(r, rx&mask|and(Comp2d(y), mask.complement()));}
		void q_rq_condition_zero	(QuatP &r, VectP const &x, QuatP const &y)	{Vect2d rx=x, mask=rx.r_is_true();	assign(r, rx&mask|and(Quat2d(y), mask.complement()));}
		void c_cr_condition_zero	(CompP &r, CompP const &x, VectP const &y)	{Comp2d cx=x; Vect2d mask=cx.c_is_true();	assign(r, and(cx, mask)|Vect2d(y)&mask.complement());}
		void c_cc_condition_zero	(CompP &r, CompP const &x, CompP const &y)	{Comp2d cx=x; Vect2d mask=cx.c_is_true();	assign(r, and(cx, mask)|and(Comp2d(y), mask.complement()));}
		void q_cq_condition_zero	(QuatP &r, CompP const &x, QuatP const &y)	{Comp2d cx=x; Vect2d mask=cx.c_is_true();	assign(r, and(cx, mask)|and(Quat2d(y), mask.complement()));}
		void q_qr_condition_zero	(QuatP &r, QuatP const &x, VectP const &y)	{Quat2d qx=x; Vect2d mask=qx.q_is_true();	assign(r, and(qx, mask)|Vect2d(y)&mask.complement());}
		void q_qc_condition_zero	(QuatP &r, QuatP const &x, CompP const &y)	{Quat2d qx=x; Vect2d mask=qx.q_is_true();	assign(r, and(qx, mask)|and(Comp2d(y), mask.complement()));}
		void q_qq_condition_zero	(QuatP &r, QuatP const &x, QuatP const &y)	{Quat2d qx=x; Vect2d mask=qx.q_is_true();	assign(r, and(qx, mask)|and(Quat2d(y), mask.complement()));}

		void  r_r_percent				(VectP &r, VectP const &x)					{assign(r, Vect2d(x)*m_one_percent);}
		void  c_c_percent				(CompP &r, CompP const &x)					{assign(r, Comp2d(x)*m_one_percent);}
		void  q_q_percent				(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x)*m_one_percent);}
	
		void r_rr_modulo				(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(x)%Vect2d(y));}
		void c_rc_modulo				(CompP &r, VectP const &x, CompP const &y)	{assign(r, Vect2d(x)%Comp2d(y));}
		void q_rq_modulo				(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, Vect2d(x)%Quat2d(y));}
		void c_cr_modulo				(CompP &r, CompP const &x, VectP const &y)	{assign(r, Comp2d(x)%Vect2d(y));}
		void c_cc_modulo				(CompP &r, CompP const &x, CompP const &y)	{assign(r, Comp2d(x)%Comp2d(y));}
		void q_cq_modulo				(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, Comp2d(x)%Quat2d(y));}
		void q_qr_modulo				(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, Quat2d(x)%Vect2d(y));}
		void q_qc_modulo				(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, Quat2d(x)%Comp2d(y));}
		void q_qq_modulo				(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, Quat2d(x)%Quat2d(y));}
	
		void r_rr_modulo_sse2			(VectP &r, VectP const &x, VectP const &y)	{assign(r, mod_sse2(Vect2d(x), Vect2d(y)));}
		void c_rc_modulo_sse2			(CompP &r, VectP const &x, CompP const &y)	{assign(r, mod_sse2(Vect2d(x), Comp2d(y)));}
		void q_rq_modulo_sse2			(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, mod_sse2(Vect2d(x), Quat2d(y)));}
		void c_cr_modulo_sse2			(CompP &r, CompP const &x, VectP const &y)	{assign(r, mod_sse2(Comp2d(x), Vect2d(y)));}
		void c_cc_modulo_sse2			(CompP &r, CompP const &x, CompP const &y)	{assign(r, mod_sse2(Comp2d(x), Comp2d(y)));}
		void q_cq_modulo_sse2			(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, mod_sse2(Comp2d(x), Quat2d(y)));}
		void q_qr_modulo_sse2			(QuatP &r, QuatP const &x, VectP const &y)	{assign(r, mod_sse2(Quat2d(x), Vect2d(y)));}
		void q_qc_modulo_sse2			(QuatP &r, QuatP const &x, CompP const &y)	{assign(r, mod_sse2(Quat2d(x), Comp2d(y)));}
		void q_qq_modulo_sse2			(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, mod_sse2(Quat2d(x), Quat2d(y)));}

		void  r_r_sgn					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, ((rx>m_zero)&m_one)-((rx<m_zero)&m_one));}
		void  c_c_sgn					(CompP &r, CompP const &x)					{Comp2d cx=x; Vect2d mag=cx.abs(); assign(r, and((cx/mag), mag.r_is_true()));}
		void  q_q_sgn					(QuatP &r, QuatP const &x)					{Quat2d qx=x; Vect2d mag=qx.abs(); assign(r, and((qx/mag), mag.r_is_true()));}
		
		__forceinline Comp2d sq(Comp2d const &x){Vect2d ri=x.r*x.i; return Comp2d(x.r*x.r-x.i*x.i, ri+ri);}
		//__forceinline Comp2d sq(Comp2d const &x){return Comp2d(x.r*x.r-x.i*x.i, m_two*x.r*x.i);}
		__forceinline Quat2d sq(Quat2d const &x)
		{
			auto _2r=x.r+x.r;
			return Quat2d(x.r*x.r-x.i*x.i-x.j*x.j-x.k*x.k, x.i*_2r, x.j*_2r, x.k*_2r);
		}
		void  r_r_sq					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, rx*rx);}
		void  c_c_sq					(CompP &r, CompP const &x)					{assign(r, sq(Comp2d(x)));}
		void  q_q_sq					(QuatP &r, QuatP const &x)					{assign(r, sq(Quat2d(x)));}

		void  c_c_sqrt					(CompP &r, CompP const &x)					{assign(r, sqrt(Comp2d(x)));}
		void  q_q_sqrt					(QuatP &r, QuatP const &x)					{assign(r, sqrt(Quat2d(x)));}

		void  r_r_invsqrt				(VectP &r, VectP const &x)
		{
			Vect2d rx=x;
			__m128i mask=_mm_set_epi32(0x5FE6EC85, 0xE7DE30DA, 0x5FE6EC85, 0xE7DE30DA);
			const Vect2d m_one_and_half=_mm_set1_pd(1.5);
			Vect2d t0=_mm_castsi128_pd(_mm_srli_epi64(_mm_castpd_si128(rx.v), 1));
			t0.v=_mm_castsi128_pd(_mm_sub_epi64(mask, _mm_castpd_si128(t0.v)));
			assign(r, t0*(m_one_and_half-m_half*rx*t0*t0));
		}

		void  r_r_cbrt					(VectP &r, VectP const &x)					{assign(r, Vect2d(::cbrt(Vect2d(x).v)));}
		void  c_c_cbrt					(CompP &r, CompP const &x)					{assign(r, exp(m_third*log(Comp2d(x))));}//
		void  q_q_cbrt					(QuatP &r, QuatP const &x)
		{
			Quat2d qx=x;
			Vect2d mag_v=qx.i*qx.i+qx.j*qx.j+qx.k*qx.k;

			Vect2d real=mag_v==m_zero, real_c=real.complement(), rr=cbrt(qx.r.v);

			Vect2d mag_x=sqrt(qx.r*qx.r+mag_v);
			mag_v=sqrt(mag_v);
			Vect2d u_mul=Vect2d(::acos((qx.r/mag_x).v))/(mag_v*m_third);
			Vect2d ln_mag_x=::log(mag_x.v);
			Quat2d result(ln_mag_x*m_third, m_pi&real|qx.i*u_mul, qx.j*u_mul, qx.k*u_mul);
			result=exp(result);
			assign(r, Quat2d(rr&real|result.r&real_c, result.i&real_c, result.j&real_c, result.k&real_c));

		//	assign(r, exp(m_third*log(Quat2d(x))));
		}

		void  r_r_gauss					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, Vect2d(::exp((-rx*rx).v)));}
		void  c_c_gauss					(CompP &r, CompP const &x)					{assign(r, exp(-sq(Comp2d(x))));}
		void  q_q_gauss					(QuatP &r, QuatP const &x)					{assign(r, exp(-sq(Quat2d(x))));}

	/*	void  r_r_erf					(VectP &r, VectP const &x)					{assign(r, Vect2d(boost::math::erf(x.lo()), boost::math::erf(x.hi()));}

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
		Vect2d  r_r_zeta					(Vect2d const &x)					{return Vect2d(zeta(x.lo()), zeta(x.hi()));}

		Vect2d  r_r_tgamma					(Vect2d const &x)					{return Vect2d(tgamma(x.lo()), tgamma(x.hi()));}
		Comp2d  c_c_tgamma					(Comp2d const &x)
		{
			Complex1d lo(x.r.lo(), x.i.lo());	lo=tgamma(lo);
			Complex1d hi(x.r.hi(), x.i.hi());	hi=tgamma(hi);
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Quat2d  q_q_tgamma					(Quat2d const &x)
		{
			Quaternion1d lo(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo());	lo=tgamma(lo);
			Quaternion1d hi(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi());	hi=tgamma(hi);
			return Quat2d(Vect2d(lo.R_component_1(), hi.R_component_1()), Vect2d(lo.R_component_2(), hi.R_component_2()),
				Vect2d(lo.R_component_3(), hi.R_component_3()), Vect2d(lo.R_component_4(), hi.R_component_4()));
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
		Vect2d r_rr_tgamma					(Vect2d const &x, Vect2d const &y)	{return Vect2d(tgamma(x.lo(), y.lo()), tgamma(x.hi(), y.hi()));}
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
			__declspec(align(16)) const Vect2d coef[6]={_mm_set1_pd(76.18009172947146), _mm_set1_pd(-86.50532032941677), _mm_set1_pd(24.01409824083091),//https://jamesmccaffrey.wordpress.com/2013/06/19/the-log-gamma-function-with-c/
				_mm_set1_pd(-1.231739572450155), _mm_set1_pd(0.1208650973866179e-2), _mm_set1_pd(-0.5395239384953e-5)};
			Vect2d rx=x, denom=rx+m_one, y=rx+Vect2d(_mm_set1_pd(5.5)), series=Vect2d(_mm_set1_pd(1.000000000190015));
			for(int i=0;i<6;++i)
			{
				series+=coef[i]/denom;
				denom+=m_one;
			}
			assign(r, Vect2d(_mm_set1_pd(_ln_sqrt_2pi))+(rx+m_half)*Vect2d(::log(y.v))-y+Vect2d(::log((series/rx).v)));
		}

	/*	Vect2d  r_r_factorial				(Vect2d const &x)					{return Vect2d(tgamma(x.lo()+1), tgamma(x.hi()+1));}
		Comp2d  c_c_factorial				(Comp2d const &x)
		{
			Complex1d lo(x.r.lo()+1, x.i.lo());	lo=tgamma(lo);
			Complex1d hi(x.r.hi()+1, x.i.hi());	hi=tgamma(hi);
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Quat2d  q_q_factorial				(Quat2d const &x)
		{
			Quaternion1d lo(x.r.lo()+1, x.i.lo(), x.j.lo(), x.k.lo());	lo=tgamma(lo);
			Quaternion1d hi(x.r.hi()+1, x.i.hi(), x.j.hi(), x.k.hi());	hi=tgamma(hi);
			return Quat2d(Vect2d(lo.R_component_1(), hi.R_component_1()), Vect2d(lo.R_component_2(), hi.R_component_2()),
				Vect2d(lo.R_component_3(), hi.R_component_3()), Vect2d(lo.R_component_4(), hi.R_component_4()));
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
		Vect2d  r_r_permutation				(Vect2d const &x)					{return m_one;}
		Comp2d  c_c_permutation				(Comp2d const &x)					{return Comp2d(m_one, m_zero);}
		Quat2d  q_q_permutation				(Quat2d const &x)					{return Quat2d(m_one, m_zero, m_zero, m_zero);}
		Vect2d r_rr_permutation				(Vect2d const &x, Vect2d const &y)	{return Vect2d(permutation(x.lo(), y.lo()), permutation(x.hi(), y.hi()));}
		Comp2d c_cr_permutation				(Comp2d const &x, Vect2d const &y)
		{
			Complex1d lo=permutation(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.lo(), 0));
			Complex1d hi=permutation(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.hi(), 0));
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_cc_permutation				(Comp2d const &x, Comp2d const &y)
		{
			Complex1d lo=permutation(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.r.lo(), y.i.lo()));
			Complex1d hi=permutation(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.r.hi(), y.i.hi()));
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Quat2d q_qq_permutation				(Quat2d const &x, Quat2d const &y)
		{
			Quaternion1d lo=permutation(Quaternion1d(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo()), Quaternion1d(y.r.lo(), y.i.lo(), y.j.lo(), y.k.lo()));
			Quaternion1d hi=permutation(Quaternion1d(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi()), Quaternion1d(y.r.hi(), y.i.hi(), y.j.hi(), y.k.hi()));
			return Quat2d(Vect2d(lo.R_component_1(), hi.R_component_1()), Vect2d(lo.R_component_2(), hi.R_component_2()),
				Vect2d(lo.R_component_3(), hi.R_component_3()), Vect2d(lo.R_component_4(), hi.R_component_4()));
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
		Vect2d  r_r_combination				(Vect2d const &x)					{return m_one;}
		Comp2d  c_c_combination				(Comp2d const &x)					{return Comp2d(m_one, m_zero);}
		Quat2d  q_q_combination				(Quat2d const &x)					{return Quat2d(m_one, m_zero, m_zero, m_zero);}
		Vect2d r_rr_combination				(Vect2d const &x, Vect2d const &y)	{return Vect2d(combination(x.lo(), y.lo()), combination(x.hi(), y.hi()));}
		Comp2d c_cr_combination				(Comp2d const &x, Vect2d const &y)
		{
			Complex1d lo=combination(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.lo(), 0));
			Complex1d hi=combination(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.hi(), 0));
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_cc_combination				(Comp2d const &x, Comp2d const &y)
		{
			Complex1d lo=combination(Complex1d(x.r.lo(), x.i.lo()), Complex1d(y.r.lo(), y.i.lo()));
			Complex1d hi=combination(Complex1d(x.r.hi(), x.i.hi()), Complex1d(y.r.hi(), y.i.hi()));
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Quat2d q_qq_combination				(Quat2d const &x, Quat2d const &y)
		{
			Quaternion1d lo=combination(Quaternion1d(x.r.lo(), x.i.lo(), x.j.lo(), x.k.lo()), Quaternion1d(y.r.lo(), y.i.lo(), y.j.lo(), y.k.lo()));
			Quaternion1d hi=combination(Quaternion1d(x.r.hi(), x.i.hi(), x.j.hi(), x.k.hi()), Quaternion1d(y.r.hi(), y.i.hi(), y.j.hi(), y.k.hi()));
			return Quat2d(Vect2d(lo.R_component_1(), hi.R_component_1()), Vect2d(lo.R_component_2(), hi.R_component_2()),
				Vect2d(lo.R_component_3(), hi.R_component_3()), Vect2d(lo.R_component_4(), hi.R_component_4()));
		}//*/

		const Comp2d m_i(m_zero, m_one);
		__forceinline Quat2d sgnu(Quat2d const &x){return Quat2d(m_zero, x.i, x.j, x.k);}
		__forceinline Quat2d acosh(Quat2d const &x){return log(x+sqrt(sq(x)-m_one));}
		__forceinline Comp2d asinh(Comp2d const &x){return log(x+sqrt(sq(x)+m_one));}
		__forceinline Quat2d asinh(Quat2d const &x){return log(x+sqrt(sq(x)+m_one));}
		__forceinline Vect2d sinhc(Vect2d const &x){return Vect2d(::sinh(x.v))/x;}
		__forceinline Comp2d cos(Comp2d const &x)
		{
			Comp2d exp_ix=exp(m_i*x);
			return m_half*exp_ix+m_half/exp_ix;
			//Vec2d sin_xr, cos_xr;
			//sin_xr=sincos(&cos_xr, x.r.v);
			//cosh(x.i.v);
		//	return (exp(m_i*x)+exp(-m_i*x))*m_half;//sincos, exp x2
		}
		__forceinline Quat2d cos(Quat2d const &x)//boost::math
		{
			Vect2d z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
			Vec2d sin_xr, cos_xr;
			sin_xr=sincos(&cos_xr, x.r.v);
			Vect2d w=-Vect2d(sin_xr)*sinhc(z);
			return Quat2d(Vect2d(cos_xr)*Vect2d(::cosh(z.v)), w*x.i, w*x.j, w*x.k);
		}
		__forceinline Comp2d sin(Comp2d const &x)
		{
			Comp2d exp_ix=exp(m_i*x);
			return m_half*exp_ix-m_half/exp_ix;
		//	return (exp(m_i*x)-exp(-m_i*x))*(-m_half*m_i);
		}
		__forceinline Quat2d sin(Quat2d const &x)//boost::math
		{
			Vect2d z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
			Vec2d sin_xr, cos_xr;
			sin_xr=sincos(&cos_xr, x.r.v);
			Vect2d w=-Vect2d(cos_xr)*sinhc(z);
			return Quat2d(Vect2d(sin_xr)*Vect2d(::cosh(z.v)), w*x.i, w*x.j, w*x.k);
		}
		void  r_r_cos					(VectP &r, VectP const &x)					{assign(r, Vect2d(::cos(Vect2d(x).v)));}
		void  c_c_cos					(CompP &r, CompP const &x)					{assign(r, cos(Comp2d(x)));}
		void  q_q_cos					(QuatP &r, QuatP const &x)					{assign(r, cos(Quat2d(x)));}

		__forceinline Comp2d acos(Comp2d const &x){return -m_i*log(x+sqrt(sq(x)-m_one));}
		__forceinline Quat2d acos(Quat2d const &x){return -sgnu(x)*acosh(x);}
		void  c_c_acos					(CompP &r, CompP const &x)					{assign(r, acos(Comp2d(x)));}
		void  q_q_acos					(QuatP &r, QuatP const &x)					{assign(r, acos(Quat2d(x)));}

		__forceinline Comp2d cosh(Comp2d const &x)
		{
			Comp2d exp_x=exp(x);
			return m_half*exp_x+m_half/exp_x;
		//	return (exp(x)+exp(-x))*m_half;
		}
		__forceinline Quat2d cosh(Quat2d const &x)
		{
			Quat2d exp_x=exp(x);
			return m_half*exp_x+m_half/exp_x;
		//	return (exp(x)+exp(-x))*m_half;
		}
		void  r_r_cosh					(VectP &r, VectP const &x)					{assign(r, Vect2d(::cosh(Vect2d(x).v)));}
		void  c_c_cosh					(CompP &r, CompP const &x)					{assign(r, cosh(Comp2d(x)));}
		void  q_q_cosh					(QuatP &r, QuatP const &x)					{assign(r, cosh(Quat2d(x)));}

		__forceinline Comp2d acosh(Comp2d const &x){return log(x+sqrt(sq(x)-m_one));}
		void  c_c_acosh					(CompP &r, CompP const &x)					{assign(r, acosh(Comp2d(x)));}
		void  q_q_acosh					(QuatP &r, QuatP const &x)					{assign(r, acosh(Quat2d(x)));}

		void  r_r_cosc					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, Vect2d(::cos(rx.v))/rx);}
		void  c_c_cosc					(CompP &r, CompP const &x)					{Comp2d cx=x; assign(r, cos(cx)/cx);}
		void  q_q_cosc					(QuatP &r, QuatP const &x)					{Quat2d qx=x; assign(r, cos(qx)/qx);}

		void  r_r_sec					(VectP &r, VectP const &x)					{assign(r, m_one/Vect2d(::cos(Vect2d(x).v)));}
		void  c_c_sec					(CompP &r, CompP const &x)					{assign(r, m_one/cos(Comp2d(x)));}
		void  q_q_sec					(QuatP &r, QuatP const &x)					{assign(r, m_one/cos(Quat2d(x)));}

		void  c_c_asec					(CompP &r, CompP const &x)					{assign(r, acos(m_one/Comp2d(x)));}
		void  q_q_asec					(QuatP &r, QuatP const &x)					{assign(r, acos(m_one/Quat2d(x)));}

		void  r_r_sech					(VectP &r, VectP const &x)					{assign(r, m_one/Vect2d(::cosh(Vect2d(x).v)));}
		void  c_c_sech					(CompP &r, CompP const &x)					{assign(r, m_one/cosh(Comp2d(x)));}
		void  q_q_sech					(QuatP &r, QuatP const &x)					{assign(r, m_one/cosh(Quat2d(x)));}

		void  c_c_asech					(CompP &r, CompP const &x)					{assign(r, acosh(m_one/Comp2d(x)));}
		void  q_q_asech					(QuatP &r, QuatP const &x)					{assign(r, acosh(m_one/Quat2d(x)));}

		void  r_r_sin					(VectP &r, VectP const &x)					{assign(r, Vect2d(::sin(Vect2d(x).v)));}
		void  c_c_sin					(CompP &r, CompP const &x)					{assign(r, sin(Comp2d(x)));}
		void  q_q_sin					(QuatP &r, QuatP const &x)					{assign(r, sin(Quat2d(x)));}

		__forceinline Comp2d asin(Comp2d const &x){return -m_i*log(m_i*x+sqrt(m_one-sq(x)));}
		__forceinline Quat2d asin(Quat2d const &x){Quat2d t=sgnu(x); return -t*asinh(x*t);}
		void  c_c_asin					(CompP &r, CompP const &x)					{assign(r, asin(Comp2d(x)));}
		void  q_q_asin					(QuatP &r, QuatP const &x)					{assign(r, asin(Quat2d(x)));}

		__forceinline Comp2d sinh(Comp2d const &x)
		{
			Comp2d exp_x=exp(x);
			return m_half*exp_x-m_half/exp_x;
		//	return (exp(x)-exp(-x))*m_half;
		}
		__forceinline Quat2d sinh(Quat2d const &x)
		{
			Quat2d exp_x=exp(x);
			return m_half*exp_x-m_half/exp_x;
		//	return (exp(x)-exp(-x))*m_half;
		}
		void  r_r_sinh					(VectP &r, VectP const &x)					{assign(r, Vect2d(::sinh(Vect2d(x).v)));}
		void  c_c_sinh					(CompP &r, CompP const &x)					{assign(r, sinh(Comp2d(x)));}
		void  q_q_sinh					(QuatP &r, QuatP const &x)					{assign(r, sinh(Quat2d(x)));}

		void  r_r_asinh					(VectP &r, VectP const &x)					{assign(r, Vect2d(::asinh(Vect2d(x).v)));}
		void  c_c_asinh					(CompP &r, CompP const &x)					{assign(r, asinh(Comp2d(x)));}
		void  q_q_asinh					(QuatP &r, QuatP const &x)					{assign(r, asinh(Quat2d(x)));}

		void  r_r_sinc					(VectP &r, VectP const &x)					{Vect2d rx=x, mask=rx==m_zero; assign(r, Vect2d(::sin(rx.v))/rx&mask.complement()|m_one&mask);}
		void  c_c_sinc					(CompP &r, CompP const &x)					{Comp2d cx=x; Vect2d mask=cx.r==m_zero&cx.i==m_zero; assign(r, and(sin(cx)/cx, mask.complement())|m_one&mask);}
		void  q_q_sinc					(QuatP &r, QuatP const &x)					{Quat2d qx=x; Vect2d mask=qx.r==m_zero&qx.i==m_zero&qx.j==m_zero&qx.k==m_zero; assign(r, and(sin(qx)/qx, mask.complement())|m_one&mask);}

		void  r_r_sinhc					(VectP &r, VectP const &x)					{Vect2d rx=x, mask=rx==m_zero; assign(r, Vect2d(::sinh(rx.v))/rx&mask.complement()|m_one&mask);}
		void  c_c_sinhc					(CompP &r, CompP const &x)					{Comp2d cx=x; Vect2d mask=cx.r==m_zero&cx.i==m_zero; assign(r, and(sinh(cx)/cx, mask.complement())|m_one&mask);}
		void  q_q_sinhc					(QuatP &r, QuatP const &x)					{Quat2d qx=x; Vect2d mask=qx.r==m_zero&qx.i==m_zero&qx.j==m_zero&qx.k==m_zero; assign(r, and(sinh(qx)/qx, mask.complement())|m_one&mask);}

		void  r_r_csc					(VectP &r, VectP const &x)					{assign(r, m_one/Vect2d(::sin(Vect2d(x).v)));}
		void  c_c_csc					(CompP &r, CompP const &x)					{assign(r, m_one/sin(Comp2d(x)));}
		void  q_q_csc					(QuatP &r, QuatP const &x)					{assign(r, m_one/sin(Quat2d(x)));}

		void  c_c_acsc					(CompP &r, CompP const &x)					{assign(r, asin(m_one/Comp2d(x)));}
		void  q_q_acsc					(QuatP &r, QuatP const &x)					{assign(r, asin(m_one/Quat2d(x)));}

		void  r_r_csch					(VectP &r, VectP const &x)					{assign(r, m_one/Vect2d(::sinh(Vect2d(x).v)));}
		void  c_c_csch					(CompP &r, CompP const &x)					{assign(r, m_one/sinh(Comp2d(x)));}
		void  q_q_csch					(QuatP &r, QuatP const &x)					{assign(r, m_one/sinh(Quat2d(x)));}

		void  r_r_acsch					(VectP &r, VectP const &x)					{assign(r, Vect2d(::asinh((m_one/Vect2d(x)).v)));}
		void  c_c_acsch					(CompP &r, CompP const &x)					{assign(r, asinh(m_one/Comp2d(x)));}
		void  q_q_acsch					(QuatP &r, QuatP const &x)					{assign(r, asinh(m_one/Quat2d(x)));}

		__forceinline Comp2d tan(Comp2d const &x)
		{
			const Comp2d m_two_i=m_i+m_i;
			Comp2d exp_2ix=exp(m_two_i*x);
			return (exp_2ix-m_one)/((exp_2ix+m_one)*m_i);
		}
		__forceinline Quat2d tan(Quat2d const &x)
		{
			const Quat2d m_two_i=m_i+m_i;
			Quat2d exp_2ix=exp(m_two_i*x);
			return (exp_2ix-m_one)/((exp_2ix+m_one)*m_i);
		}
		void  r_r_tan					(VectP &r, VectP const &x)					{assign(r, Vect2d(::tan(Vect2d(x).v)));}
		void  c_c_tan					(CompP &r, CompP const &x)					{assign(r, tan(Comp2d(x)));}
		void  q_q_tan					(QuatP &r, QuatP const &x)					{assign(r, tan(Quat2d(x)));}

		__forceinline Comp2d atan(Comp2d const &x){return (m_i*m_half)*log((m_i+x)/(m_i-x));}
		__forceinline Quat2d atan(Quat2d const &x){return (m_i*m_half)*log((m_i+x)/(m_i-x));}
		void  r_r_atan					(VectP &r, VectP const &x)					{assign(r, Vect2d(::atan(Vect2d(x).v)));}
		void  c_c_atan					(CompP &r, CompP const &x)					{assign(r, atan(Comp2d(x)));}
		void  q_q_atan					(QuatP &r, QuatP const &x)					{assign(r, atan(Quat2d(x)));}
		void r_rr_atan					(VectP &r, VectP const &y, VectP const &x)	{assign(r, Vect2d(::atan2(Vect2d(y).v, Vect2d(x).v)));}
		void c_rc_atan					(CompP &r, VectP const &y, CompP const &x)
		{
			Comp2d cx=x;	Vect2d ry=y;
			Comp2d tr=atan(ry/cx);
			Vect2d mask_y=ry<m_zero, mask_x=cx.r<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void q_rq_atan					(QuatP &r, VectP const &y, QuatP const &x)
		{
			Quat2d qx=x;	Vect2d ry=y;
			Quat2d tr=atan(ry/qx);
			Vect2d mask_y=ry<m_zero, mask_x=qx.r<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void c_cr_atan					(CompP &r, CompP const &y, VectP const &x)
		{
			Vect2d rx=x;	Comp2d cy=y;
			Comp2d tr=atan(cy/rx);
			Vect2d mask_y=cy.r<m_zero, mask_x=rx<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void c_cc_atan					(CompP &r, CompP const &y, CompP const &x)
		{
			Comp2d cx=x;	Comp2d cy=y;
			Comp2d tr=atan(cy/cx);
			Vect2d mask_y=cy.r<m_zero, mask_x=cx.r<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void q_cq_atan					(QuatP &r, CompP const &y, QuatP const &x)
		{
			Quat2d qx=x;	Comp2d cy=y;
			auto tr=atan(cy/qx);
			Vect2d mask_y=cy.r<m_zero, mask_x=qx.r<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void q_qr_atan					(QuatP &r, QuatP const &y, VectP const &x)
		{
			Vect2d rx=x;	Quat2d qy=y;
			auto tr=atan(qy/rx);
			Vect2d mask_y=qy.r<m_zero, mask_x=rx<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void q_qc_atan					(QuatP &r, QuatP const &y, CompP const &x)
		{
			Comp2d cx=x;	Quat2d qy=y;
			auto tr=atan(qy/cx);
			Vect2d mask_y=qy.r<m_zero, mask_x=cx.r<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}
		void q_qq_atan					(QuatP &r, QuatP const &y, QuatP const &x)
		{
			Quat2d qx=x;	Quat2d qy=y;
			auto tr=atan(qy/qx);
			Vect2d mask_y=qy.r<m_zero, mask_x=qx.r<m_zero;
			Vect2d addition=m_pi&mask_y;
			Vect2d add_sign=m_sign_mask_complement&mask_x;
			addition^=add_sign;
			assign(r, tr+addition);
		}

		__forceinline Comp2d tanh(Comp2d const &x){Comp2d e2x=exp(x+x); return (e2x-m_one)/(e2x+m_one);}
		__forceinline Quat2d tanh(Quat2d const &x){Quat2d e2x=exp(x+x); return (e2x-m_one)/(e2x+m_one);}
		void  r_r_tanh					(VectP &r, VectP const &x)					{assign(r, Vect2d(::tanh(Vect2d(x).v)));}
		void  c_c_tanh					(CompP &r, CompP const &x)					{assign(r, tanh(Comp2d(x)));}
		void  q_q_tanh					(QuatP &r, QuatP const &x)					{assign(r, tanh(Quat2d(x)));}

		__forceinline Comp2d atanh(Comp2d const &x){return m_half*log((m_one+x)/(m_one-x));}
		__forceinline Quat2d atanh(Quat2d const &x){return m_half*log((m_one+x)/(m_one-x));}
		void  c_c_atanh					(CompP &r, CompP const &x)					{assign(r, atanh(Comp2d(x)));}
		void  q_q_atanh					(QuatP &r, QuatP const &x)					{assign(r, atanh(Quat2d(x)));}

		void  r_r_tanc					(VectP &r, VectP const &x)					{Vect2d rx=x, mask=rx==m_zero; assign(r, Vect2d(::tan(rx.v))/rx&mask.complement()|m_one&mask);}
		void  c_c_tanc					(CompP &r, CompP const &x)					{Comp2d cx=x; Vect2d mask=cx.r==m_zero&cx.i==m_zero; assign(r, and(tan(cx)/cx, mask.complement())|m_one&mask);}
		void  q_q_tanc					(QuatP &r, QuatP const &x)					{Quat2d qx=x; Vect2d mask=qx.r==m_zero&qx.i==m_zero&qx.j==m_zero&qx.k==m_zero; assign(r, and(tan(qx)/qx, mask.complement())|m_one&mask);}

		void  r_r_cot					(VectP &r, VectP const &x)					{assign(r, m_one/Vect2d(::tan(Vect2d(x).v)));}
		void  c_c_cot					(CompP &r, CompP const &x)					{assign(r, m_one/tan(Comp2d(x)));}
		void  q_q_cot					(QuatP &r, QuatP const &x)					{assign(r, m_one/tan(Quat2d(x)));}

		void  r_r_acot					(VectP &r, VectP const &x)
		{
			Vect2d rx=x, mask=rx==m_zero;
			assign(r, Vect2d(::atan((m_one/rx).v))&mask.complement()|m_pi_2&mask);
		}
		void  c_c_acot					(CompP &r, CompP const &x)
		{
			Comp2d cx=x;
			Vect2d mask=(cx.r==m_zero)&(cx.i==m_zero);
			assign(r, and(atan(m_one/cx), mask.complement())|m_pi_2&mask);
		}
		void  q_q_acot					(QuatP &r, QuatP const &x)
		{
			Quat2d qx=x;
			Vect2d mask=(qx.r==m_zero)&(qx.i==m_zero)&(qx.j==m_zero)&(qx.k==m_zero);
			assign(r, and(atan(m_one/qx), mask.complement())|m_pi_2&mask);
		}

		void  r_r_coth					(VectP &r, VectP const &x)					{assign(r, m_one/Vect2d(::tanh(Vect2d(x).v)));}
		void  c_c_coth					(CompP &r, CompP const &x)					{assign(r, m_one/tanh(Comp2d(x)));}
		void  q_q_coth					(QuatP &r, QuatP const &x)					{assign(r, m_one/tanh(Quat2d(x)));}

		void  c_c_acoth					(CompP &r, CompP const &x)					{assign(r, atanh(m_one/Comp2d(x)));}
		void  q_q_acoth					(QuatP &r, QuatP const &x)					{assign(r, atanh(m_one/Quat2d(x)));}

		void  r_r_exp					(VectP &r, VectP const &x)					{assign(r, Vect2d(::exp(Vect2d(x).v)));}
		void  c_c_exp					(CompP &r, CompP const &x)					{assign(r, exp(Comp2d(x)));}
		void  q_q_exp					(QuatP &r, QuatP const &x)					{assign(r, exp(Quat2d(x)));}
	
		__forceinline Vect2d pow(Vect2d const &x, Vect2d const &y){return Vect2d(::exp((y*Vect2d(::log(x.v))).v));}
		const Vect2d m_ln_phi=::log(m_phi.v), m_inv_sqrt5=m_one/m_sqrt5;
		void  r_r_fib					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, (Vect2d(::exp((rx*m_ln_phi).v))-Vect2d(::cos((m_pi*rx).v))*Vect2d(::exp((-rx*m_ln_phi).v)))*m_inv_sqrt5);}
	//	void  r_r_fib						(Vect2d const &x)					{return (pow(m_phi, x)-Vect2d(::cos((m_pi*x).v))*pow(m_phi, -x))/m_sqrt5;}
	//	void  r_r_fib						(Vect2d const &x)					{return (Vect2d(::pow(m_phi.v, x.v))-Vect2d(::cos((m_pi*x).v))*Vect2d(::pow(m_phi.v, (-x).v)))/m_sqrt5;}
		void  c_c_fib					(CompP &r, CompP const &x)					{Comp2d cx=x; assign(r, (exp(cx*m_ln_phi)-cos(m_pi*cx)*exp(-cx*m_ln_phi))*m_inv_sqrt5);}
	//	void  c_c_fib					(CompP &r, CompP const &x)					{Comp2d cx=x; assign(r, ((m_phi^cx)-cos(m_pi*cx)*(m_phi^-cx))/m_sqrt5;}
		void  q_q_fib					(QuatP &r, QuatP const &x)					{Quat2d qx=x; assign(r, (exp(qx*m_ln_phi)-cos(m_pi*qx)*exp(-qx*m_ln_phi))*m_inv_sqrt5);}
	
	/*	const double rand_norm=9.31322574615479e-010;
		Vect2d my_rand(){return Vect2d((rand()<<15|rand())*rand_norm, (rand()<<15|rand())*rand_norm);}
		Vect2d  r_r_random					(Vect2d const &x)					{return my_rand();}
		Comp2d  c_c_random					(Comp2d const &x)					{return Comp2d(my_rand(), my_rand());}
		Quat2d  q_q_random					(Quat2d const &x)					{return Quat2d(my_rand(), my_rand(), my_rand(), my_rand());}
		Vect2d r_rr_random					(Vect2d const &x, Vect2d const &y)	{return my_rand();}
		Comp2d c_cr_random					(Comp2d const &x, Vect2d const &y)	{return Comp2d(my_rand(), my_rand());}
		Comp2d c_cc_random					(Comp2d const &x, Comp2d const &y)	{return Comp2d(my_rand(), my_rand());}
		Quat2d q_qq_random					(Quat2d const &x, Quat2d const &y)	{return Quat2d(my_rand(), my_rand(), my_rand(), my_rand());}

		__forceinline double beta(double x, double y)
		{
			try
			{
				return boost::math::beta(x, y);
			}
			catch(...){return _qnan;}
		}
		Vect2d  r_r_beta					(Vect2d const &x)					{return Vect2d(beta(x.lo(), x.lo()), beta(x.hi(), x.hi()));}
		Vect2d r_rr_beta					(Vect2d const &x, Vect2d const &y)	{return Vect2d(beta(x.lo(), y.lo()), beta(x.hi(), y.hi()));}

		__forceinline double cyl_bessel_j(double x, double y)
		{
			try
			{
				return boost::math::cyl_bessel_j(x, y);
			}
			catch(...){return _qnan;}
		}
		Vect2d  r_r_cyl_bessel_j			(Vect2d const &x)					{return Vect2d(cyl_bessel_j(0, x.lo()), cyl_bessel_j(0, x.hi()));}
		Vect2d r_rr_cyl_bessel_j			(Vect2d const &x, Vect2d const &y)	{return Vect2d(cyl_bessel_j(x.lo(), y.lo()), cyl_bessel_j(x.hi(), y.hi()));}

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
		Vect2d  r_r_cyl_neumann				(Vect2d const &x)					{return Vect2d(cyl_neumann(0, x.lo()), cyl_neumann(0, x.hi()));}
		Vect2d r_rr_cyl_neumann				(Vect2d const &x, Vect2d const &y)	{return Vect2d(cyl_neumann(x.lo(), y.lo()), cyl_neumann(x.hi(), y.hi()));}

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
		Comp2d  c_r_hankel1					(Vect2d const &x)
		{
			Complex1d lo=r_hankel1(x.lo(), x.lo()), hi=r_hankel1(x.hi(), x.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d  c_c_hankel1					(Comp2d const &x)
		{
			Complex1d lo=r_hankel1(x.r.lo(), x.i.lo()), hi=r_hankel1(x.r.hi(), x.i.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}
		Comp2d c_rr_hankel1					(Vect2d const &x, Vect2d const &y)
		{
			Complex1d lo=r_hankel1(x.lo(), y.lo()), hi=r_hankel1(x.hi(), y.hi());
			return Comp2d(Vect2d(lo.real(), hi.real()), Vect2d(lo.imag(), hi.imag()));
		}//*/

		__forceinline Vect2d sgn(Vect2d const &x){return ((x>m_zero)&m_one)-((x<m_zero)&m_one);}
		__forceinline Comp2d sgn(Comp2d const &x){Vect2d mag=x.abs(); return and(x/mag, mag.r_is_true());}
		__forceinline Quat2d sgn(Quat2d const &x){Vect2d mag=x.abs(); return and(x/mag, mag.r_is_true());}
		__forceinline Vect2d step(Vect2d const &x){return m_half+m_half*sgn(x);}
		__forceinline Comp2d step(Comp2d const &x){return m_half+m_half*sgn(x);}
		__forceinline Quat2d step(Quat2d const &x){return m_half+m_half*sgn(x);}
		void  r_r_step					(VectP &r, VectP const &x)					{assign(r, step(Vect2d(x)));}
		void  c_c_step					(CompP &r, CompP const &x)					{assign(r, step(Comp2d(x)));}
		void  q_q_step					(QuatP &r, QuatP const &x)					{assign(r, step(Quat2d(x)));}

		void  r_r_rect					(VectP &r, VectP const &x)					{Vect2d rx=x; assign(r, step(rx+m_half)-step(rx-m_half));}
		void  c_c_rect					(CompP &r, CompP const &x)					{Comp2d cx=x; assign(r, step(cx+m_half)-step(cx-m_half));}
		void  q_q_rect					(QuatP &r, QuatP const &x)					{Quat2d qx=x; assign(r, step(qx+m_half)-step(qx-m_half));}

		void  r_r_trgl					(VectP &r, VectP const &x)					{Vect2d t=Vect2d(x).abs(); assign(r, (m_one-t)&t<m_one);}
		void  r_c_trgl					(VectP &r, CompP const &x)					{Vect2d t=Comp2d(x).abs(); assign(r, and(m_one-t, t<m_one));}
		void  r_q_trgl					(VectP &r, QuatP const &x)					{Vect2d t=Quat2d(x).abs(); assign(r, and(m_one-t, t<m_one));}

		void  r_r_sqwv					(VectP &r, VectP const &x)					{Vect2d rx=x;	assign(r, rx-rx.floor()<m_half&m_one);}
		void  r_c_sqwv					(VectP &r, CompP const &x)					{Vect2d rx=x.r; assign(r, rx-rx.floor()<m_half&m_one);}
		void  r_q_sqwv					(VectP &r, QuatP const &x)					{Vect2d rx=x.r; assign(r, rx-rx.floor()<m_half&m_one);}
		void r_rr_sqwv					(VectP &r, VectP const &x, VectP const &y)	{Vect2d rx=x;	assign(r, rx-rx.floor()<Vect2d(y)&m_one);}
		void r_rc_sqwv					(VectP &r, VectP const &x, CompP const &y)	{Vect2d rx=x;	assign(r, rx-rx.floor()<Vect2d(y.r)&m_one);}
		void r_rq_sqwv					(VectP &r, VectP const &x, QuatP const &y)	{Vect2d rx=x;	assign(r, rx-rx.floor()<Vect2d(y.r)&m_one);}
		void r_cr_sqwv					(VectP &r, CompP const &x, VectP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor()<Vect2d(y)&m_one);}
		void r_cc_sqwv					(VectP &r, CompP const &x, CompP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor()<Vect2d(y.r)&m_one);}
		void r_cq_sqwv					(VectP &r, CompP const &x, QuatP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor()<Vect2d(y.r)&m_one);}
		void r_qr_sqwv					(VectP &r, QuatP const &x, VectP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor()<Vect2d(y)&m_one);}
		void r_qc_sqwv					(VectP &r, QuatP const &x, CompP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor()<Vect2d(y.r)&m_one);}
		void r_qq_sqwv					(VectP &r, QuatP const &x, QuatP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor()<Vect2d(y.r)&m_one);}
		
		void  r_r_sqwv_sse2				(VectP &r, VectP const &x)					{Vect2d rx=x;	assign(r, rx-rx.floor_sse2()<m_half&m_one);}
		void  r_c_sqwv_sse2				(VectP &r, CompP const &x)					{Vect2d rx=x.r; assign(r, rx-rx.floor_sse2()<m_half&m_one);}
		void  r_q_sqwv_sse2				(VectP &r, QuatP const &x)					{Vect2d rx=x.r; assign(r, rx-rx.floor_sse2()<m_half&m_one);}
		void r_rr_sqwv_sse2				(VectP &r, VectP const &x, VectP const &y)	{Vect2d rx=x;	assign(r, rx-rx.floor_sse2()<Vect2d(y)&m_one);}
		void r_rc_sqwv_sse2				(VectP &r, VectP const &x, CompP const &y)	{Vect2d rx=x;	assign(r, rx-rx.floor_sse2()<Vect2d(y.r)&m_one);}
		void r_rq_sqwv_sse2				(VectP &r, VectP const &x, QuatP const &y)	{Vect2d rx=x;	assign(r, rx-rx.floor_sse2()<Vect2d(y.r)&m_one);}
		void r_cr_sqwv_sse2				(VectP &r, CompP const &x, VectP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor_sse2()<Vect2d(y)&m_one);}
		void r_cc_sqwv_sse2				(VectP &r, CompP const &x, CompP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor_sse2()<Vect2d(y.r)&m_one);}
		void r_cq_sqwv_sse2				(VectP &r, CompP const &x, QuatP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor_sse2()<Vect2d(y.r)&m_one);}
		void r_qr_sqwv_sse2				(VectP &r, QuatP const &x, VectP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor_sse2()<Vect2d(y)&m_one);}
		void r_qc_sqwv_sse2				(VectP &r, QuatP const &x, CompP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor_sse2()<Vect2d(y.r)&m_one);}
		void r_qq_sqwv_sse2				(VectP &r, QuatP const &x, QuatP const &y)	{Vect2d rx=x.r;	assign(r, rx-rx.floor_sse2()<Vect2d(y.r)&m_one);}

		void  r_r_trwv					(VectP &r, VectP const &x)					{Vect2d rx=x; Vect2d t=(rx-rx.floor()-m_half).abs(); assign(r, t+t);}
		void  r_c_trwv					(VectP &r, CompP const &x)					{Comp2d cx=x; Vect2d t=(cx-cx.floor()-m_half).abs(); assign(r, t+t);}
		void  r_q_trwv					(VectP &r, QuatP const &x)					{Quat2d qx=x; Vect2d t=(qx-qx.floor()-m_half).abs(); assign(r, t+t);}
		void r_rr_trwv					(VectP &r, VectP const &x, VectP const &y)
		{
			Vect2d rx=x;	Vect2d ry=y;
			Vect2d t=rx-rx.floor(), t2=m_one-rx;
			t2=t2-t2.floor();
			Vect2d d=ry&ry>m_zero;
			Vect2d mask=d>m_one;
			d=d&mask.complement()|m_one&mask;
			Vect2d d2=m_one-d;
			Vect2d t_d=t/d, t2_d2=t2/d2;
			assign(r, ((t_d<m_one)&t_d)+((t2_d2<m_one)&t2_d2));
		}
		void c_cr_trwv					(CompP &r, CompP const &x, VectP const &y)
		{
			Comp2d cx=x;	Vect2d ry=y;
			auto t=cx-cx.floor(), t2=m_one-cx;
			t2=t2-t2.floor();
			auto d=and(ry, ry>m_zero);
			auto mask=d>m_one;
			d=and(d, mask.complement())|m_one&mask;
			auto d2=m_one-d;
			auto t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
		}
		void c_cc_trwv					(CompP &r, CompP const &x, CompP const &y)
		{
			Comp2d cx=x;	Comp2d cy=y;
			Comp2d t=cx-cx.floor(), t2=m_one-cx;
			t2=t2-t2.floor();
			Comp2d d=and(cy, cy.r>m_zero);
			Vect2d mask=d.r>m_one;
			d=and(d, mask.complement())|m_one&mask;
			Comp2d d2=m_one-d;
			Comp2d t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
		}
		void q_qq_trwv					(QuatP &r, QuatP const &x, QuatP const &y)
		{
			Quat2d qx=x;	Quat2d qy=y;
			Quat2d t=qx-qx.floor(), t2=m_one-qx;
			t2=t2-t2.floor();
			Quat2d d=and(qy, qy.r>m_zero);
			Vect2d mask=d.r>m_one;
			d=and(d, mask.complement())|m_one&mask;
			Quat2d d2=m_one-d;
			Quat2d t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
		}

		void  r_r_trwv_sse2				(VectP &r, VectP const &x)					{Vect2d rx=x; Vect2d t=(rx-rx.floor_sse2()-m_half).abs(); assign(r, t+t);}
		void  r_c_trwv_sse2				(VectP &r, CompP const &x)					{Comp2d cx=x; Vect2d t=(cx-cx.floor_sse2()-m_half).abs(); assign(r, t+t);}
		void  r_q_trwv_sse2				(VectP &r, QuatP const &x)					{Quat2d qx=x; Vect2d t=(qx-qx.floor_sse2()-m_half).abs(); assign(r, t+t);}
		void r_rr_trwv_sse2				(VectP &r, VectP const &x, VectP const &y)
		{
			Vect2d rx=x;	Vect2d ry=y;
			Vect2d t=rx-rx.floor_sse2(), t2=m_one-rx;
			t2=t2-t2.floor_sse2();
			Vect2d d=ry&(ry<m_zero).complement();
			Vect2d mask=d>m_one;
			d=d&mask.complement()|m_one&mask;
			Vect2d d2=m_one-d;
			Vect2d t_d=t/d, t2_d2=t2/d2;
			assign(r, ((t_d<m_one)&t_d)+((t2_d2<m_one)&t2_d2));
		}
		void c_cr_trwv_sse2				(CompP &r, CompP const &x, VectP const &y)
		{
			Comp2d cx=x;	Vect2d ry=y;
			auto t=cx-cx.floor_sse2(), t2=m_one-cx;
			t2=t2-t2.floor_sse2();
			auto d=and(ry, (ry<m_zero).complement());
			auto mask=d>m_one;
			d=and(d, mask.complement())|m_one&mask;
			auto d2=m_one-d;
			auto t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
		}
		void c_cc_trwv_sse2				(CompP &r, CompP const &x, CompP const &y)
		{
			Comp2d cx=x;	Comp2d cy=y;
			Comp2d t=cx-cx.floor_sse2(), t2=m_one-cx;
			t2=t2-t2.floor_sse2();
			Comp2d d=and(cy, (cy.r<m_zero).complement());
			Vect2d mask=d.r>m_one;
			d=and(d, mask.complement())|m_one&mask;
			Comp2d d2=m_one-d;
			Comp2d t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
		}
		void q_qq_trwv_sse2				(QuatP &r, QuatP const &x, QuatP const &y)
		{
			Quat2d qx=x;	Quat2d qy=y;
			Quat2d t=qx-qx.floor_sse2(), t2=m_one-qx;
			t2=t2-t2.floor_sse2();
			Quat2d d=and(qy, (qy.r<m_zero).complement());
			Vect2d mask=d.r>m_one;
			d=and(d, mask.complement())|m_one&mask;
			Quat2d d2=m_one-d;
			Quat2d t_d=t/d, t2_d2=t2/d2;
			assign(r, and((t_d.r<m_one), t_d)+and((t2_d2.r<m_one), t2_d2));
		}

		void  r_r_saw					(VectP &r, VectP const &x)
		{
			Vect2d rx=x;
			Vect2d t=rx-rx.floor(), t2=(m_one-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t));
		}
		void  c_c_saw					(CompP &r, CompP const &x)
		{
			Comp2d cx=x;
			auto t=cx-cx.floor(), t2=(m_one-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t));
		}
		void  q_q_saw					(QuatP &r, QuatP const &x)
		{
			Quat2d qx=x;
			auto t=qx-qx.floor(), t2=(m_one-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t));
		}
		void r_rr_saw					(VectP &r, VectP const &x, VectP const &y)
		{
			Vect2d rx=x;	Vect2d ry=y;
			auto t=rx-rx.floor(), t2=(ry-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
		}
		void c_rc_saw					(CompP &r, VectP const &x, CompP const &y)
		{
			Vect2d rx=x;	Comp2d cy=y;
			auto t=rx-rx.floor();
			auto t2=(cy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
		}
		void q_rq_saw					(QuatP &r, VectP const &x, QuatP const &y)
		{
			Vect2d rx=x;	Quat2d qy=y;
			auto t=rx-rx.floor();
			auto t2=(qy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
		}
		void c_cr_saw					(CompP &r, CompP const &x, VectP const &y)
		{
			Comp2d cx=x;	Vect2d ry=y;
			auto t=cx-cx.floor();
			auto t2=(ry-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
		}
		void c_cc_saw					(CompP &r, CompP const &x, CompP const &y)
		{
			Comp2d cx=x, cy=y;
			auto t=cx-cx.floor();
			auto t2=(cy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
		}
		void q_cq_saw					(QuatP &r, CompP const &x, QuatP const &y)
		{
			Comp2d cx=x;	Quat2d qy=y;
			auto t=cx-cx.floor();
			auto t2=(qy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
		}
		void q_qr_saw					(QuatP &r, QuatP const &x, VectP const &y)
		{
			Quat2d qx=x;	Vect2d ry=y;
			auto t=qx-qx.floor();
			auto t2=(ry-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
		}
		void q_qc_saw					(QuatP &r, QuatP const &x, CompP const &y)
		{
			Quat2d qx=x;	Comp2d cy=y;
			auto t=qx-qx.floor();
			auto t2=(cy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
		}
		void q_qq_saw					(QuatP &r, QuatP const &x, QuatP const &y)
		{
			Quat2d qx=x, qy=y;
			auto t=qx-qx.floor();
			auto t2=(qy-t).floor();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
		}

		void  r_r_saw_sse2				(VectP &r, VectP const &x)
		{
			Vect2d rx=x;
			Vect2d t=rx-rx.floor_sse2(), t2=(m_one-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t));
		}
		void  c_c_saw_sse2				(CompP &r, CompP const &x)
		{
			Comp2d cx=x;
			auto t=cx-cx.floor_sse2(), t2=(m_one-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t));
		}
		void  q_q_saw_sse2				(QuatP &r, QuatP const &x)
		{
			Quat2d qx=x;
			auto t=qx-qx.floor_sse2(), t2=(m_one-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t));
		}
		void r_rr_saw_sse2				(VectP &r, VectP const &x, VectP const &y)
		{
			Vect2d rx=x;	Vect2d ry=y;
			auto t=rx-rx.floor_sse2(), t2=(ry-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
		}
		void c_rc_saw_sse2				(CompP &r, VectP const &x, CompP const &y)
		{
			Vect2d rx=x;	Comp2d cy=y;
			auto t=rx-rx.floor_sse2();
			auto t2=(cy-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
		}
		void q_rq_saw_sse2				(QuatP &r, VectP const &x, QuatP const &y)
		{
			Vect2d rx=x;	Quat2d qy=y;
			auto t=rx-rx.floor_sse2();
			auto t2=(qy-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
		}
		void c_cr_saw_sse2				(CompP &r, CompP const &x, VectP const &y)
		{
			Comp2d cx=x;	Vect2d ry=y;
			auto t=cx-cx.floor_sse2();
			auto t2=(ry-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
		}
		void c_cc_saw_sse2				(CompP &r, CompP const &x, CompP const &y)
		{
			Comp2d cx=x, cy=y;
			auto t=cx-cx.floor_sse2();
			auto t2=(cy-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
		}
		void q_cq_saw_sse2				(QuatP &r, CompP const &x, QuatP const &y)
		{
			Comp2d cx=x;	Quat2d qy=y;
			auto t=cx-cx.floor_sse2();
			auto t2=(qy-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
		}
		void q_qr_saw_sse2				(QuatP &r, QuatP const &x, VectP const &y)
		{
			Quat2d qx=x;	Vect2d ry=y;
			auto t=qx-qx.floor_sse2();
			auto t2=(ry-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/ry);
		}
		void q_qc_saw_sse2				(QuatP &r, QuatP const &x, CompP const &y)
		{
			Quat2d qx=x;	Comp2d cy=y;
			auto t=qx-qx.floor_sse2();
			auto t2=(cy-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/cy);
		}
		void q_qq_saw_sse2				(QuatP &r, QuatP const &x, QuatP const &y)
		{
			Quat2d qx=x, qy=y;
			auto t=qx-qx.floor_sse2();
			auto t2=(qy-t).floor_sse2();
			assign(r, (t2+m_one)*(t2*m_half+t)/qy);
		}

		void r_rr_hypot					(VectP &r, VectP const &x, VectP const &y)	{Vect2d rx=x, ry=y; assign(r, sqrt(rx*rx+ry*ry));}
	//	void c_cc_hypot					(CompP &r, CompP const &x, CompP const &y)	{Comp2d cx=x; assign(r, sqrt(sq(x)+sq(y));}
	//	void q_qq_hypot					(QuatP &r, QuatP const &x, QuatP const &y)	{Quat2d qx=x; assign(r, sqrt(sq(x)+sq(y));}
		
		__forceinline Vect2d mandelbrot(Vect2d const &x, Vect2d const &y)
		{
			Vect2d rez, sq_rez, m_sq_limit=16;
			Vect2d k, active=m_one;
			for(;;k+=active)
			{
				//condition
				active&=sq_rez<m_sq_limit;
				Vect2d cont=k<y&active.r_is_true();
				if(!_mm_movemask_pd(cont))
				//Vect2d stop=k>=y|active.r_is_false();
				//if(stop.lo()&&stop.hi())
					break;
				//iteration
				rez=sq_rez;//calculate sq(z)

				rez+=x;//add x

				sq_rez=rez*rez;
			}
			return k;
		}
		__forceinline Vect2d mandelbrot(Comp2d const &x, Vect2d const &y)
		{
			Vect2d rez, imz, sq_rez, sq_imz, m_sq_limit=16;
			Vect2d k, active=m_one;
			for(;;k+=active)
			{
				//condition
				active&=sq_rez+sq_imz<m_sq_limit;
				Vect2d cont=k<y&active.r_is_true();
				if(!_mm_movemask_pd(cont))
			//	Vect2d stop=k>=y|active.r_is_false();
			//	stop&=_mm_castsi128_pd(_mm_shuffle_epi32(_mm_castpd_si128(stop), _MM_SHUFFLE(3, 2, 3, 2)));
			//	if(stop.lo())
					break;
				//if(stop.lo()&&stop.hi())
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
		void r_r_mandelbrot				(VectP &r, VectP const &x)					{assign(r, mandelbrot(Vect2d(x), Vect2d(200)));}
		void r_c_mandelbrot				(VectP &r, CompP const &x)					{assign(r, mandelbrot(Comp2d(x), Vect2d(200)));}
		void r_rr_mandelbrot			(VectP &r, VectP const &x, VectP const &y)	{assign(r, mandelbrot(x, y));}
		void r_cr_mandelbrot			(VectP &r, CompP const &x, VectP const &y)	{assign(r, mandelbrot(Comp2d(x), Vect2d(y)));}

		void r_rr_min					(VectP &r, VectP const &x, VectP const &y)	{Vect2d rx=x, ry=y;			assign(r, (rx+ry-(rx-ry).abs())*m_half);}
		void c_cr_min					(CompP &r, CompP const &x, VectP const &y)	{Comp2d cx=x; Vect2d ry=y;	assign(r, (cx+ry-(cx-ry).abs())*m_half);}
		void c_cc_min					(CompP &r, CompP const &x, CompP const &y)	{Comp2d cx=x; Comp2d cy=y;	assign(r, (cx+cy-(cx-cy).abs())*m_half);}
		void q_qq_min					(QuatP &r, QuatP const &x, QuatP const &y)	{Quat2d qx=x, qy=y;			assign(r, (qx+qy-(qx-qy).abs())*m_half);}

		void r_rr_max					(VectP &r, VectP const &x, VectP const &y)	{Vect2d rx=x, ry=y;			assign(r, (rx+ry+(rx-ry).abs())*m_half);}
		void c_cr_max					(CompP &r, CompP const &x, VectP const &y)	{Comp2d cx=x; Vect2d ry=y;	assign(r, (cx+ry+(cx-ry).abs())*m_half);}
		void c_cc_max					(CompP &r, CompP const &x, CompP const &y)	{Comp2d cx=x; Comp2d cy=y;	assign(r, (cx+cy+(cx-cy).abs())*m_half);}
		void q_qq_max					(QuatP &r, QuatP const &x, QuatP const &y)	{Quat2d qx=x, qy=y;			assign(r, (qx+qy+(qx-qy).abs())*m_half);}

		void r_rr_conditional_110		(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(y)&Vect2d(x)!=m_zero);}
		void c_rc_conditional_110		(CompP &r, VectP const &x, CompP const &y)	{assign(r, and(Comp2d(y), Vect2d(x)!=m_zero));}
		void q_rq_conditional_110		(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, and(Quat2d(y), Vect2d(x)!=m_zero));}
		void r_cr_conditional_110		(VectP &r, CompP const &x, VectP const &y)	{assign(r, and(Vect2d(y), Comp2d(x)!=m_zero));}
		void c_cc_conditional_110		(CompP &r, CompP const &x, CompP const &y)	{assign(r, and(Comp2d(y), Comp2d(x)!=m_zero));}
		void q_cq_conditional_110		(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, and(Quat2d(y), Comp2d(x)!=m_zero));}
		void r_qr_conditional_110		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, and(Vect2d(y), Quat2d(x)!=m_zero));}
		void c_qc_conditional_110		(CompP &r, QuatP const &x, CompP const &y)	{assign(r, and(Comp2d(y), Quat2d(x)!=m_zero));}
		void q_qq_conditional_110		(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, and(Quat2d(y), Quat2d(x)!=m_zero));}
	
		void r_rr_conditional_101		(VectP &r, VectP const &x, VectP const &y)	{assign(r, Vect2d(y)&Vect2d(x)==m_zero);}
		void c_rc_conditional_101		(CompP &r, VectP const &x, CompP const &y)	{assign(r, and(Comp2d(y), Vect2d(x)==m_zero));}
		void q_rq_conditional_101		(QuatP &r, VectP const &x, QuatP const &y)	{assign(r, and(Quat2d(y), Vect2d(x)==m_zero));}
		void r_cr_conditional_101		(VectP &r, CompP const &x, VectP const &y)	{assign(r, and(Vect2d(y), Comp2d(x)==m_zero));}
		void c_cc_conditional_101		(CompP &r, CompP const &x, CompP const &y)	{assign(r, and(Comp2d(y), Comp2d(x)==m_zero));}
		void q_cq_conditional_101		(QuatP &r, CompP const &x, QuatP const &y)	{assign(r, and(Quat2d(y), Comp2d(x)==m_zero));}
		void r_qr_conditional_101		(VectP &r, QuatP const &x, VectP const &y)	{assign(r, and(Vect2d(y), Quat2d(x)==m_zero));}
		void c_qc_conditional_101		(CompP &r, QuatP const &x, CompP const &y)	{assign(r, and(Comp2d(y), Quat2d(x)==m_zero));}
		void q_qq_conditional_101		(QuatP &r, QuatP const &x, QuatP const &y)	{assign(r, and(Quat2d(y), Quat2d(x)==m_zero));}

		void conditional_111			(QuatP &res, QuatP const &op1, QuatP const &op2, QuatP const &op3, int idx, int op1_ms, int op_ms)
		{
			Vect2d second=
				op1_ms=='R'?	Vect2d(op1.r+idx).r_is_false()
				:op1_ms=='c'?	Comp2d(op1.r+idx, op1.i+idx).c_is_false()
				:				Quat2d(op1.r+idx, op1.i+idx, op1.j+idx, op1.k+idx).q_is_false();
			switch(op_ms)
			{
			case 0:assign(VectP(res.r+idx),											Vect2d(op2.r+idx)&second.complement()									|	Vect2d(op3.r+idx)&second);break;//r_rr
			case 1:assign(CompP(res.r+idx, res.i+idx),							and(Comp2d(op2.r+idx, op2.i+idx), second.complement())						|	Vect2d(op3.r+idx)&second);break;//c_rc
			case 2:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Quat2d(op2.r+idx, op2.i+idx, op2.j+idx, op2.k+idx), second.complement())|	Vect2d(op3.r+idx)&second);break;//q_rq
			case 3:assign(CompP(res.r+idx, res.i+idx),								Vect2d(op2.r+idx)&second.complement()									|and(Comp2d(op3.r+idx, op3.i+idx), second));break;//c_cr
			case 4:assign(CompP(res.r+idx, res.i+idx),							and(Comp2d(op2.r+idx, op2.i+idx), second.complement())						|and(Comp2d(op3.r+idx, op3.i+idx), second));break;//c_cc
			case 5:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Quat2d(op2.r+idx, op2.i+idx, op2.j+idx, op2.k+idx), second.complement())|and(Comp2d(op3.r+idx, op3.i+idx), second));break;//q_cq
			case 6:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),		Vect2d(op2.r+idx)&second.complement()									|and(Quat2d(op3.r+idx, op3.i+idx, op3.j+idx, op3.k+idx), second));break;//q_qr
			case 7:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Comp2d(op2.r+idx, op2.i+idx), second.complement())						|and(Quat2d(op3.r+idx, op3.i+idx, op3.j+idx, op3.k+idx), second));break;//q_qc
			case 8:assign(QuatP(res.r+idx, res.i+idx, res.j+idx, res.k+idx),	and(Quat2d(op2.r+idx, op2.i+idx, op2.j+idx, op2.k+idx), second.complement())|and(Quat2d(op3.r+idx, op3.i+idx, op3.j+idx, op3.k+idx), second));break;//q_qq
			}
		}

		void  r_r_increment				(VectP &r, VectP const &x)					{assign(r, Vect2d(x)+m_one);}
		void  c_c_increment				(CompP &r, CompP const &x)					{assign(r, Comp2d(x)+m_one);}
		void  q_q_increment				(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x)+m_one);}

		void  r_r_decrement				(VectP &r, VectP const &x)					{assign(r, Vect2d(x)-m_one);}
		void  c_c_decrement				(CompP &r, CompP const &x)					{assign(r, Comp2d(x)-m_one);}
		void  q_q_decrement				(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x)-m_one);}

		void  r_r_assign				(VectP &r, VectP const &x)					{assign(r, Vect2d(x));}
		void  c_c_assign				(CompP &r, CompP const &x)					{assign(r, Comp2d(x));}
		void  q_q_assign				(QuatP &r, QuatP const &x)					{assign(r, Quat2d(x));}
	}//namespace sse2
}//namespace G2
void		lighten_sse2(int *rgb, int n)
{
	int k=0, kEnd=n&0xFFFFFFFC;
	__m128i mask=_mm_set1_epi16(0x7F7F);
	for(;k<kEnd;k+=4)
	{
		__m128i v=_mm_loadu_si128((__m128i*)(rgb+k));
		v=_mm_xor_si128(v, m_minus_one);
		v=_mm_srli_epi16(v, 1);
		v=_mm_and_si128(v, mask);
		v=_mm_xor_si128(v, m_minus_one);
		_mm_storeu_si128((__m128i*)(rgb+k), v);
	}
	for(k-=4;k<n;++k)
	{
		auto p=(unsigned char*)(rgb+k);
		rgb[k]=~rgb[k];
		p[0]>>=1, p[1]>>=1, p[2]>>=1;
		rgb[k]=~rgb[k];
	}
}
namespace	modes
{
	void colorFunction_bcw_sse2(CompP const &_v, int *rgb)
	{
		Comp2d v=_v;
		Vect2d threshold=10;
		double inv_th=0.1;//1			//black->color->white			66.07 cycles/px
		Vect2d const cos_pi_6=0.866025403784439, sin_pi_6=0.5, _FF=255, _7F=0x7F;
		Vect2d hyp=v.abs(), costh=v.r/hyp, sinth=v.i/hyp,
			mag=_FF*exp(-hyp*Vect2d(G2::_ln2*inv_th));
		Vect2d
			red		=m_one+costh*cos_pi_6-sinth*sin_pi_6,
			green	=m_one+sinth,
			blue	=m_one+costh*-cos_pi_6-sinth*sin_pi_6;
		Vect2d small=hyp<threshold, small_c=small.complement();
		//mag=_mm_blendv_pd(mag, _FF-mag, small);
		//red	=_mm_blendv_pd(red	*mag, _FF-mag*(m_two-red	), small);//SSE4.1 CRASHES on Core 2 Duo T7500
		//green	=_mm_blendv_pd(green*mag, _FF-mag*(m_two-green	), small);
		//blue	=_mm_blendv_pd(blue	*mag, _FF-mag*(m_two-blue	), small);
		mag=_FF-mag&small|mag&small_c;
		red		=red	*mag&small|_FF-mag*(m_two-red	)&small_c;
		green	=green	*mag&small|_FF-mag*(m_two-green	)&small_c;
		blue	=blue	*mag&small|_FF-mag*(m_two-blue	)&small_c;
		Vect2d
			nan=v.r!=v.r|v.i!=v.i, nan_c=nan.complement(),
			inf=v.r==m_inf|v.i==m_inf, inf_c=inf.complement();
		//red	=_mm_blendv_pd(red,		_FF, inf);
		//green	=_mm_blendv_pd(green,	_FF, inf);
		//blue	=_mm_blendv_pd(blue,	_FF, inf);
		//red	=_mm_blendv_pd(red,		_7F, nan);
		//green	=_mm_blendv_pd(green,	_7F, nan);
		//blue	=_mm_blendv_pd(blue,	_7F, nan);
		red		=_FF&inf|red	&inf_c;
		green	=_FF&inf|green	&inf_c;
		blue	=_FF&inf|blue	&inf_c;
		red		=_7F&nan|red	&nan_c;
		green	=_7F&nan|green	&nan_c;
		blue	=_7F&nan|blue	&nan_c;
		auto p=(unsigned char*)rgb;
		p[0]=unsigned char(blue.lo()), p[1]=unsigned char(green.lo()), p[2]=unsigned char(red.lo()), p[3]=0;
		p[4]=unsigned char(blue.hi()), p[5]=unsigned char(green.hi()), p[6]=unsigned char(red.hi()), p[7]=0;
	}
	void colorFunction_bc_l_sse2(CompP const &_v, int *rgb)
	{
		Comp2d v=_v;
		Vect2d const cos_pi_6=0.866025403784439, sin_pi_6=0.5, _FF=255, _7F=0x7F;		//32.31 cycles/px
		Vect2d hyp=v.abs(), costh=v.r/hyp, sinth=v.i/hyp,
			fh=(hyp*Vect2d(0.1)).floor(), c, f;
		c=fh*Vect2d(0.01);
		c=_FF*(c-c.floor()), f=Vect2d(12.75)*(hyp-fh*Vect2d(10))*(m_one-c*Vect2d(1./255));
		Vect2d
			red=m_one+costh*cos_pi_6-sinth*sin_pi_6,
			green=m_one+sinth,
			blue=m_one+costh*-cos_pi_6-sinth*sin_pi_6;
		red=c+f*red, green=c+f*green, blue=c+f*blue;
		Vect2d
			nan=v.r!=v.r|v.i!=v.i, nan_c=nan.complement(),
			inf=v.r==m_inf|v.i==m_inf, inf_c=inf.complement();
		red		=_FF&inf|red	&inf_c;
		green	=_FF&inf|green	&inf_c;
		blue	=_FF&inf|blue	&inf_c;
		red		=_7F&nan|red	&nan_c;
		green	=_7F&nan|green	&nan_c;
		blue	=_7F&nan|blue	&nan_c;
		auto p=(unsigned char*)rgb;
		p[0]=unsigned char(blue.lo()), p[1]=unsigned char(green.lo()), p[2]=unsigned char(red.lo()), p[3]=0;
		p[4]=unsigned char(blue.hi()), p[5]=unsigned char(green.hi()), p[6]=unsigned char(red.hi()), p[7]=0;

	/*	Vect2d r=pr, i=pi;
		c1=-1, c2=-1;
		double r1=r.v.m128d_f64[0], i1=i.v.m128d_f64[0], r2=r.v.m128d_f64[1], i2=i.v.m128d_f64[1];
		colorFunction_cases(r1, i1, c1);
		colorFunction_cases(r2, i2, c2);
		if((c1==-1)+(c2==-1)==2)
		{
			//const Vect2d m_inf=_mm_set1_pd(_HUGE), m_minf=_mm_set1_pd(-_HUGE);
			//const __m128i
			//	nan_color=_mm_set1_epi32(0x007F7F7F), cinf_color=_mm_set1_epi32(0x00FFFFFF),
			//	pinf_color=_mm_set1_epi32(0x00ED7F11), ninf_color=_mm_set1_epi32(0x00117FED),
			//	iinf_color=_mm_set1_epi32(0x003FFF3F), niinf_color=_mm_set1_epi32(0x00BF00BF);
			//Vect2d ri_nan=r!=r|i!=i;
			//Vect2d
			//	r_inf=r==m_inf, r_minf=r==m_minf,
			//	i_inf=i==m_inf, i_minf=i==m_minf, c_inf=(r_inf|r_minf)&(i_inf|i_minf), ri_zero=r==m_zero&i==m_zero;

			const Vect2d m_amp=_mm_set1_pd(255/G2::_pi), m_cospi6=_mm_set1_pd(0.866025403784439), m_sinpi6=_mm_set1_pd(0.5);
			Vect2d hyp=sqrt(r*r+i*i), mag=Vect2d(atan(hyp.v))*m_amp, _1_hyp=m_one/hyp,
				cos_x=r*_1_hyp, sin_x=i*_1_hyp;
			Vect2d cosx_cospi6=cos_x*m_cospi6, sinx_sinpi6=sin_x*m_sinpi6;
			Vect2d m_r=mag*(m_one+cosx_cospi6-sinx_sinpi6), m_g=mag*(m_one+sin_x), m_b=mag*(m_one-cosx_cospi6-sinx_sinpi6);
			auto p=(unsigned char*)&c1;
			p[0]=unsigned char(m_b.v.m128d_f64[0]), p[1]=unsigned char(m_g.v.m128d_f64[0]), p[2]=unsigned char(m_r.v.m128d_f64[0]);//argb
			p=(unsigned char*)&c2;
			p[0]=unsigned char(m_b.v.m128d_f64[1]), p[1]=unsigned char(m_g.v.m128d_f64[1]), p[2]=unsigned char(m_r.v.m128d_f64[1]);
		}
		else
		{
			if(c1==-1)
				c1=colorFunction(r1, i1);
			if(c2==-1)
				c2=colorFunction(r2, i2);
		}//*/
	}
}