//g2_common.h - Include for AVX library & G2 executable
//best viewed with tab size of 4 spaces
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

#ifndef		COMMON_H
#define		COMMON_H
#ifndef		AVX_LIB
#if defined _DEBUG && !defined __GNUC__
#pragma warning(push)
#pragma warning(disable:4996)
#endif
//#include	<boost/random/mersenne_twister.hpp>
//#include	<boost/random/uniform_real.hpp>
#include	<boost/math/special_functions.hpp>
#include	<boost/math/quaternion.hpp>
//#include	<boost/math/octonion.hpp>
#include	<complex>
#if defined _DEBUG && !defined __GNUC__
#pragma warning(pop)
#endif
#endif
#include	<limits>
#include	<cmath>

	#define	CHECK_NULL_POINTERS
	#define	ALIGNED_INTRINSICS//load/store vs loadu/storeu

#ifdef __GNUC__
#define		ALIGN(mul)	__attribute__((aligned(mul)))
#else
#define		ALIGN(mul)	__declspec(align(mul))
#endif
#ifdef CHECK_NULL_POINTERS
#define				DEREF(pointer)		(pointer?*pointer:0)
#else
#define				DEREF(pointer)		(*pointer)
#endif
inline float		clamp_positive(float x){return (x+abs(x))*0.5f;}
inline int			clamp_positive(int x){return (x+abs(x))>>1;}
inline int			minimum(int a, int b){return (a+b-abs(a-b))>>1;}
inline int			maximum(int a, int b){return (a+b+abs(a-b))>>1;}
inline int			minimum(int a, int b, int c)
{
	int a2=a<<1, temp=b+c-abs(b-c);
	return (a2+temp-abs(a2-temp))>>2;
}
inline int			maximum(int a, int b, int c)
{
	int a2=a<<1, temp=b+c+abs(b-c);
	return (a2+temp+abs(a2+temp))>>2;
}
inline double		maximum(double a, double b){return (a+b+abs(a-b))*0.5;}
inline double		minimum(double a, double b, double c)
{
	double a2=a+a, temp=b+c-abs(b-c);
	return (a2+temp-abs(a2-temp))*0.5;
}
namespace	G2
{
	extern const double
		_atm,	_bbr,
		_c,		_ele,
		_e,		_g,
		_G,		_h,
		_mag,	_me,
		_mn,	_mp,
		_Me,	_Ms,
		_Na,	_phi,
		_pi,	_q,
		_R,		_qnan,

		_ln2,	_ln10,	inv_ln10,	_sqrt2,		_sqrt3,			_sqrt5,
		_pi_2,	_2pi,	_sqrt_2pi,	_ln_sqrt_2pi,	_1_2pi,		_1_pi,
		_pi_180,
		_ln_phi, inv_sqrt5;
	extern const float		s_zero, s_infinity;
//	extern const double		zero, infinity;
}
#ifndef AVX_LIB
typedef std::complex<double> Comp1d;
typedef boost::math::quaternion<double> Quat1d;
inline bool istrue(Comp1d const &x){return (x.real()!=0)|(x.imag()!=0);}
inline bool istrue(Quat1d const &x){return (x.R_component_1()!=0)|(x.R_component_2()!=0)|(x.R_component_3()!=0)|(x.R_component_4()!=0);}
#if _MSC_VER<1700 && !defined __GNUC__//before VS2012
namespace	std
{
	inline double round(double x){return floor(x+0.5);}
	inline bool signbit(double x){return (((char*)&x)[7]>>7)!=0;}
	inline bool isnan(double x){return x!=x;}
	inline bool isinf(double x){return abs(x)==_HUGE;}
	double lgamma(double x);
	const Comp1d i(0, 1), i_2(0, 0.5);
	inline double asinh(double x){return log(x+sqrt(x*x+1.));}
	inline Comp1d acos(Comp1d const &z){return -i*log(z+sqrt(z*z-1.));}
	inline Comp1d asin(Comp1d const &z){return -i*log(i*z+sqrt(1.-z*z));}
	inline Comp1d atan(Comp1d const &z){return i_2*log((i+z)/(i-z));}
	inline Comp1d acosh(Comp1d const &z){return log(z+sqrt(z*z-1.));}
	inline Comp1d asinh(Comp1d const &z){return log(z+sqrt(z*z+1.));}
	inline Comp1d atanh(Comp1d const &z){return 0.5*log((1.+z)/(1.-z));}
}
#endif
namespace	G2
{
	bool _2d_between(double x1, double y1, double x, double y, double x2, double y2);
	inline bool _1d_int_in_range(double x0, double x1){return std::floor(x0)!=std::floor(x1)||std::ceil(x0)!=std::ceil(x1);}
	inline bool _1d_zero_in_range(double x0, double x1){return x0<0?x1>=0:x0==0?x1!=0:x1<0;}
	inline double _1d_zero_crossing(double x0, double y0, double x1, double y1){return x0+(0-y0)*(x1-x0)/(y1-y0);}
	extern const double ll_max;
	inline long long convert_d2ll		(double const x){return x!=x||x<-ll_max||x>ll_max?(long long&)x	:(long long)x;}
	inline long long convert_d2ll_zero	(double const x){return x!=x||x<-ll_max||x>ll_max?0				:(long long)x;}

	//inline bool		istrue		(Comp1d const &x){return x.real()||x.imag();}
	//inline bool		istrue		(Quat1d const &x){return x.R_component_1()||x.R_component_2()||x.R_component_3()||x.R_component_4();}
	inline Comp1d	floor		(Comp1d const &x){return Comp1d(::floor(x.real()), ::floor(x.imag()));}
	inline Quat1d	floor		(Quat1d	const &x){return Quat1d(::floor(x.R_component_1()), ::floor(x.R_component_2()), ::floor(x.R_component_3()), ::floor(x.R_component_4()));}
	inline Comp1d	ceil		(Comp1d const &x){return Comp1d(::ceil(x.real()), ::ceil(x.imag()));}
	inline Quat1d	ceil		(Quat1d	const &x){return Quat1d(::ceil(x.R_component_1()), ::ceil(x.R_component_2()), ::ceil(x.R_component_3()), ::ceil(x.R_component_4()));}
	inline Comp1d	round		(Comp1d const &x){return Comp1d(std::round(x.real()), std::round(x.imag()));}
	inline Quat1d	round		(Quat1d	const &x){return Quat1d(std::round(x.R_component_1()), std::round(x.R_component_2()), std::round(x.R_component_3()), std::round(x.R_component_4()));}
	inline Comp1d	operator%	(double const &x, Comp1d const &y){return x-floor(x/y)*y;}
	inline Quat1d	operator%	(double const &x, Quat1d const &y){return x-floor(x/y)*y;}
	inline Comp1d	operator%	(Comp1d const &x, double const &y){return x-floor(x/y)*y;}
	inline Comp1d	operator%	(Comp1d const &x, Comp1d const &y){return x-floor(x/y)*y;}
	inline Quat1d	operator%	(Comp1d const &x, Quat1d const &y){return x-floor(x/y)*y;}
	inline Quat1d	operator%	(Quat1d const &x, double const &y){return x-floor(x/y)*y;}
	inline Quat1d	operator%	(Quat1d const &x, Comp1d const &y){return x-floor(x/y)*y;}
	inline Quat1d	operator%	(Quat1d const &x, Quat1d const &y){return x-floor(x/y)*y;}
//	inline double	abs_unreal	(Quat1d const &x){return ::sqrt(x.R_component_2()*x.R_component_2()+x.R_component_3()*x.R_component_3()+x.R_component_4()*x.R_component_4();}
	inline double	arg			(Quat1d	const &x){return ::acos(x.R_component_1()/boost::math::abs(x));}
		   Quat1d	log			(Quat1d	const &x);
	inline Quat1d	pow			(Quat1d	const &x, Quat1d	const &y){return boost::math::exp(y*log(x));}
	inline Quat1d	pow			(double const &x, Quat1d	const &y){return boost::math::exp(y*std::log(Comp1d(x)));}
	inline double	bitwise_xor	(double const &x)
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
	}
	inline double	bitwise_xnor(double const &x)
	{
		long long t1=convert_d2ll(x);
		t1^=t1>>32, t1^=t1>>16, t1^=t1>>8, t1^=t1>>4;
		t1&=15;
		return !((0x6996>>t1)&1);
	}

	inline double	sgn			(double const &x){return (x>0)-(x<0);}
	inline Comp1d	sgn			(Comp1d const &x)
	{
		double temp=std::abs(x);
		return temp?x/temp:Comp1d();
	}
	inline Quat1d	sgn			(Quat1d	const &x)
	{
		double temp=abs(x);
		return temp!=0?x/temp:Quat1d();
	}
	inline Quat1d	sgnu		(Quat1d	const &x){return sgn(Quat1d(0, x.R_component_2(), x.R_component_3(), x.R_component_4()));}
		   double	tgamma		(double const &x);
		   Comp1d	tgamma		(Comp1d const &x);
		   Quat1d	tgamma		(Quat1d	const &x);
	inline Comp1d	sq			(Comp1d const &x){double ri=x.real()*x.imag(); return Comp1d(x.real()*x.real()-x.imag()*x.imag(), ri+ri);}
	inline Quat1d	sq			(Quat1d const &x)
	{
		auto _2r=x.R_component_1()+x.R_component_1();
		return Quat1d(x.R_component_1()*x.R_component_1()-x.R_component_2()*x.R_component_2()-x.R_component_3()*x.R_component_3()-x.R_component_4()*x.R_component_4(),
			x.R_component_2()*_2r, x.R_component_3()*_2r, x.R_component_4()*_2r);
	}
	inline Quat1d	sqrt		(Quat1d	const &x)//2 sqrts
	{
		double s=abs(x)+x.R_component_1();
		if(s)//entire quaternion space except origin & -ve real axis
		{
			s=::sqrt(s+s);
			double inv_s=1/s;
			return Quat1d(s*0.5, x.R_component_2()*inv_s, x.R_component_3()*inv_s, x.R_component_4()*inv_s);
		}
		return Quat1d(0, ::sqrt(-x.R_component_1()), 0, 0);
		//return exp(0.5*log(x));
	}
	inline Quat1d	acosh		(Quat1d	const &x){return log(x+sqrt(sq(x)-1.));}
	inline Quat1d	asinh		(Quat1d	const &x){return log(x+sqrt(sq(x)+1.));}
	inline Quat1d	atanh		(Quat1d	const &x){return (log(1.+x)-log(1.-x))/2.;}
	inline Quat1d	acos		(Quat1d	const &x){return -sgnu(x)*acosh(x);}
	inline Quat1d	asin		(Quat1d	const &x){Quat1d q2=sgnu(x); return -q2*asinh(x*q2);}
	inline Quat1d	atan		(Quat1d	const &x){Quat1d q2=sgnu(x); return -q2*atanh(x*q2);}
	inline double	step		(double const &x){return 0.5+0.5*sgn(x);}
	inline Comp1d	step		(Comp1d const &x){return 0.5+0.5*sgn(x);}
	inline Quat1d	step		(Quat1d	const &x){return 0.5+0.5*sgn(x);}
}
#endif
struct VectP
{
	double *r;
	VectP():r(nullptr){}
	VectP(double *r):r(r){}
	void set(double *r){this->r=r;}
	bool r_is_true()const{return DEREF(r)!=0;}

	operator double()const{return DEREF(r);}

	VectP& operator=(double x){*r=x; return *this;}
	VectP& operator=(VectP const &x){*r=*x.r; return *this;}
	VectP& operator+=(VectP const &x){*r+=*x.r; return *this;}
	VectP& operator-=(VectP const &x){*r-=*x.r; return *this;}
	VectP& operator*=(VectP const &x){*r*=*x.r; return *this;}
	VectP& operator/=(VectP const &x){*r/=*x.r; return *this;}
	VectP& operator%=(VectP const &q){*r-=::floor(*r/ *q.r)**q.r;return *this;}//x%q = x-floor(x/q)*q

	VectP& operator&=(VectP const &other){*(long long*)r&=*(long long*)other.r; return *this;}//
	VectP& operator|=(VectP const &other){*(long long*)r|=*(long long*)other.r; return *this;}//
	VectP& operator^=(VectP const &other){*(long long*)r^=*(long long*)other.r; return *this;}//

	VectP& ceil_this(){*r=::ceil(*r);return *this;}
	double ceil()const{return ::ceil(*r);}
	VectP& floor_this(){*r=::floor(*r);return *this;}
	double floor()const{return ::floor(*r);}
	VectP& round_this(){*r=std::round(*r);return *this;}
	double round()const{return std::round(*r);}
	double complement()const{long long c=~*(long long*)r; return (double&)c;}
	VectP& abs_this(){*r=::abs(*r); return *this;}
	double abs()const{return ::abs(*r);}
};
//inline double operator*(VectP const &a, VectP const &b){return *a.r**b.r;}
//inline double operator/(VectP const &a, VectP const &b){return *a.r/ *b.r;}
inline double operator%(VectP const &x, VectP const &q){return *x.r-floor(*x.r/ *q.r)**q.r;}//x%q = x-floor(x/q)*q
//inline double operator+(VectP const &a, VectP const &b){return *a.r+*b.r;}
//inline double operator-(VectP const &a, VectP const &b){return *a.r-*b.r;}
//inline double operator-(VectP const &x){return -*x.r;}
//inline bool operator>(VectP const &a, VectP const &b){return *a.r>*b.r;}
//inline bool operator<(VectP const &a, VectP const &b){return *a.r<*b.r;}
//inline bool operator>=(VectP const &a, VectP const &b){return *a.r>=*b.r;}
//inline bool operator<=(VectP const &a, VectP const &b){return *a.r<=*b.r;}
//inline bool operator==(VectP const &a, VectP const &b){return *a.r==*b.r;}
//inline bool operator!=(VectP const &a, VectP const &b){return *a.r!=*b.r;}
inline double operator&(VectP const &a, VectP const &b){long long t=*(long long*)a.r&*(long long*)b.r;return (double&)t;}//
inline double operator|(VectP const &a, VectP const &b){long long t=*(long long*)a.r|*(long long*)b.r;return (double&)t;}//
inline double operator^(VectP const &a, VectP const &b){long long t=*(long long*)a.r^*(long long*)b.r;return (double&)t;}//
inline double AND(VectP const &a, VectP const &b){return a&b;}//
//inline double sqrt(VectP const &x){return sqrt(*x.r);}
inline bool	is_not_nan(VectP const &x){return *x.r==*x.r;}
inline bool	isinf(VectP const &x){return x.abs()==_HUGE;}
inline void	reduce_angle(double &th){th=th-G2::_2pi*floor(th*G2::_1_2pi);}
inline void	update_angle(double &th, double &cth, double &sth)
{
	th=th-G2::_2pi*floor(th*G2::_1_2pi);
	cth=cos(th), sth=sin(th);
}

/*struct Quat1d;
struct Comp1d
{
	double r, i;
	Comp1d():r(0), i(0){}
	Comp1d(double r, double i):r(r), i(i){}
	bool r_is_true()const{return r!=0;}
	bool c_is_true()const{return r!=0||i!=0;}
	Comp1d floor()const{return Comp1d(::floor(r), ::floor(i));}
	Comp1d ceil()const{return Comp1d(::ceil(r), ::ceil(i));}
	Comp1d round()const{return Comp1d(::round(r), ::round(i));}
	double abs()const{return sqrt(r*r+i*i);}
	double arg()const{return ::atan2(i, r);}
	Comp1d& operator+=(Comp1d const &b){r+=b.r, i+=b.i; return *this;}
	Comp1d& operator-=(Comp1d const &b){r-=b.r, i-=b.i; return *this;}
	Comp1d& operator*=(Comp1d const &b)
	{
		double rr=r*b.r-i*b.i, ri=r*b.i+i*b.r;
		r=rr, i=ri;
		return *this;
	}
	Comp1d& operator*=(double const &b){r*=b, i*=b; return *this;}
	Comp1d& operator/=(Comp1d const &b)
	{
		double _1_mag_b=1/sqrt(b.r*b.r+b.i*b.i);
		double
			rr=(b.r*r+b.i*i)*_1_mag_b,
			ri=(b.r*i-b.i*r)*_1_mag_b;
		r=rr, i=ri;
		return *this;
	}
	Comp1d& operator/=(double const &br){r/=br, i/=br; return *this;}
	Quat1d& operator/=(Quat1d const &b);
	Comp1d& operator^=(Comp1d const &b)
	{
		Comp1d t(::log(sqrt(r*r+i*i)), atan2(i, r));
		t*=b;
		double r0=::exp(t.r);
		r=r0*cos(t.i), i=r0*sin(t.i);
		return *this;
	}
	Comp1d& operator^=(double const &br)
	{
		Comp1d t(::log(sqrt(r*r+i*i)), atan2(i, r));
		t*=br;
		double r0=::exp(t.r);
		r=r0*cos(t.i), i=r0*sin(t.i);
		return *this;
	}
};//*/
struct CompP
{
	double *r, *i;
	CompP():r(nullptr), i(nullptr){}
	CompP(double *r, double *i):r(r), i(i){}
	void set(double *r, double *i){this->r=r, this->i=i;}
	bool r_is_true()const{return DEREF(r)!=0;}
	bool c_is_true()const{return DEREF(r)!=0||DEREF(i)!=0;}
#ifndef AVX_LIB
	operator Comp1d()const{return Comp1d(*r, *i);}
	CompP& operator=(Comp1d const &x){*r=x.real(), *i=x.imag(); return *this;}
	//CompP& operator=(double const &x){*r=x, *i=0; return *this;}

/*	Comp1d floor()const{return Comp1d(::floor(*r), ::floor(*i));}
	Comp1d ceil()const{return Comp1d(::ceil(*r), ::ceil(*i));}
	Comp1d round()const{return Comp1d(::round(*r), ::round(*i));}
	double abs()const{return sqrt(*r**r+*i**i);}
	double arg()const{return ::atan2(*i, *r);}
	CompP& operator+=(CompP const &b){*r+=*b.r, *i+=*b.i;}
	CompP& operator-=(CompP const &b){*r-=*b.r, *i-=*b.i; return *this;}
	CompP& operator*=(CompP const &b)
	{
		double rr=*r**b.r-*i**b.i, ri=*r**b.i+*i**b.r;
		*r=rr, *i=ri;
		return *this;
	}
	CompP& operator*=(double const &b){*r*=b, *i*=b; return *this;}
	CompP& operator/=(CompP const &b)
	{
		double _1_mag_b=1/sqrt(*b.r**b.r+*b.i**b.i);
		double
			rr=(*b.r**r+*b.i**i)*_1_mag_b,
			ri=(*b.r**i-*b.i**r)*_1_mag_b;
		*r=rr, *i=ri;
		return *this;
	}
	CompP& operator/=(double const &br){*r/=br, *i/=br; return *this;}
	QuatP& operator/=(QuatP const &b);
	CompP& operator^=(CompP const &b)
	{
		Comp1d t(::log(sqrt(*r**r+*i**i)), atan2(*i, *r));
		t*=(Comp1d)b;
		double r0=::exp(t.r);
		*r=r0*cos(t.i), *i=r0*sin(t.i);
	}
	CompP& operator^=(double const &br)
	{
		Comp1d t(::log(sqrt(*r**r+*i**i)), atan2(*i, *r));
		t*=br;
		double r0=::exp(t.r);
		*r=r0*cos(t.i), *i=r0*sin(t.i);
	}//*/
#endif
};

struct QuatP
{
	double *r, *i, *j, *k;
	QuatP():r(nullptr), i(nullptr), j(nullptr), k(nullptr){}
	QuatP(double *r, double *i, double *j, double *k):r(r), i(i), j(j), k(k){}
	void set(double *r, double *i, double *j, double *k){this->r=r, this->i=i, this->j=j, this->k=k;}
	bool r_is_true()const{return DEREF(r)!=0;}
	bool c_is_true()const{return DEREF(r)!=0||DEREF(i)!=0;}
	bool q_is_true()const{return DEREF(r)!=0||DEREF(i)!=0||DEREF(j)!=0||DEREF(k)!=0;}
#ifndef AVX_LIB
	operator Quat1d()const{return Quat1d(DEREF(r), DEREF(i), DEREF(j), DEREF(k));}
	QuatP& operator=(Quat1d const &x){*r=x.R_component_1(), *i=x.R_component_2(), *j=x.R_component_3(), *k=x.R_component_4(); return *this;}
	//QuatP& operator=(Comp1d const &x){*r=x.real(), *i=x.imag(), *j=0, *k=0; return *this;}
	//QuatP& operator=(double const &x){*r=x, *i=0, *j=0, *k=0; return *this;}
#endif
};
//namespace	modes
//{
//	int colorFunction(double &r, double &i);
//	void colorFunction_cases(double r, double i, int &c);
//}
template<typename T, int Align>struct aligned_vector
{
	typedef T *iterator;
	typedef T const *const_iterator;
	T *p;
	int n;
	aligned_vector():p(0), n(0){}
	aligned_vector(aligned_vector<T, Align> const &other)
	{
		if(&other!=this)
		{
			n=other.n;
			if(n)
			{
				p=(T*)_aligned_malloc(other.n*sizeof(T), Align);
				memcpy(p, other.p, n*sizeof(T));
			}
			else
				p=0;
		}
	}
	aligned_vector(aligned_vector<T, Align> &&other)
	{
		if(&other!=this)
		{
			n=other.n, p=other.p;
			other.n=0, other.p=0;
		}
	}
	aligned_vector(unsigned n, T const &e=T()):n(n), p((T*)_aligned_malloc(n*sizeof(T), Align))
	{
		fill(0, n, e);
	}
	~aligned_vector()
	{
		if(p)
		{
			_aligned_free(p);
			p=0;
		}
	}
	aligned_vector& operator=(aligned_vector const &other)
	{
		if(&other!=this)
		{
			n=other.n;
			realloc(n);
			for(int k=0;k<n;++k)
				p[k]=other[k];
		}
		return *this;
	}
	aligned_vector& operator=(aligned_vector &&other)
	{
		if(&other!=this)
		{
			if(p)
				_aligned_free(p);
			n=other.n, p=other.p;
			other.n=0, other.p=0;
		}
		return *this;
	}
	unsigned size()const{return n;}
	iterator begin(){return p;}
	iterator end(){return p+n;}
	const_iterator begin()const{return p;}
	const_iterator end()const{return p+n;}
	void realloc(unsigned n)
	{
		if(p)
			p=(T*)_aligned_realloc(p, n*sizeof(T), Align);
		else
			p=(T*)_aligned_malloc(n*sizeof(T), Align);
	}
	void fill(unsigned start, unsigned end, T const &e)
	{
		for(unsigned k=start;k<end;++k)
			p[k]=e;
	}
	void resize(unsigned n, T const &e=T())
	{
		if(n!=this->n)//
		{
			realloc(n);
			fill(this->n, n, e);
			this->n=n;
		}
	}
	T& operator[](unsigned k){return p[k];}
	T const& operator[](unsigned k)const{return p[k];}
	void assign(unsigned n, T const &e=T())
	{
		this->n=n;
		realloc(n);
		fill(0, n, e);
	}
	void push_back(T const &e)
	{
		++n;
		realloc(n);
		p[n-1]=e;
	}
	void pop_back()
	{
		--n;
		realloc(n);
	}
	void insert(unsigned position, T const &e, unsigned count=1)
	{
		unsigned n2=n+count;
		realloc(n2);
		n=n2;
		unsigned pos2=position+count;
		for(unsigned k=n2-1;k>=pos2;--k)
			p[k]=p[k-count];
		fill(position, pos2, e);
	}
	void erase(unsigned position, unsigned count=1)
	{
		for(unsigned k=position, kEnd=n-count;k<kEnd;++k)
			p[k]=p[k+count];
		n-=count;
		realloc(n);
	}
	void clear(){n=0; if(p)_aligned_free(p), p=0;}
};
typedef aligned_vector<double, 32> AVector_v4d;
#define				TERM_CONTENTS	bool constant; char mathSet; int varNo; AVector_v4d r, i, j, k;
struct				ArgIdx//for built-in variadic & user-defined functions
{
	short idx;
	char mathSet;
	ArgIdx():idx(0), mathSet(0){}
	ArgIdx(short idx, char mathSet):idx(idx), mathSet(mathSet){}
};
#endif//COMMON_H