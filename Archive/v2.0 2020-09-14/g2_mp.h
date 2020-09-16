#ifndef MP_H
#define MP_H
#include		<mpreal.h>

	#define		MP_PURE_QUAT

namespace		MP
{
#if 1
	extern const double inv_ln2, ln2_ln10, ln2_ln8, ln2_ln16;
	extern int bin_prec, dec_prec;
	void	recalculate_constants();
	inline void set_prec(int prec, int base=10)//supported bases: 2, 8, 10, 16
	{
		bin_prec=base==2?prec:int(prec*::log((double)base)*inv_ln2);
		dec_prec=base==10?prec:int(bin_prec*ln2_ln10);
		mpfr_set_default_prec(bin_prec);
		recalculate_constants();
	}
	typedef mpfr::mpreal Real;
	inline Real operator~(Real const &x)
	{
		int prec=x.getPrecision();
		auto mask=(Real(1)<<prec)-1;
		return mask-x;
	}
	extern Real m_pi, m_e, m_qnan, m_ln2, m_inv_ln10, m_sqrt_2pi, m_ln_phi, m_inv_sqrt5;
	struct		Quat;
	struct		Comp
	{
		Real r, i;
		Comp(Real const &r=0, Real const &i=0):r(r), i(i){}
		bool c_is_true()const{return r!=0||i!=0;}
		Comp& operator+=(Comp const &b){r+=b.r, i+=b.i; return *this;}
		Comp& operator-=(Comp const &b){r-=b.r, i-=b.i; return *this;}
		Comp& operator*=(Comp const &b)
		{
			Real rr=r*b.r-i*b.i, ri=r*b.i+i*b.r;
			r=rr, i=ri;
			return *this;
		}
		Comp& operator*=(Real const &b){r*=b, i*=b; return *this;}
		Comp& operator/=(Comp const &b)
		{
			Real inv_mag_b=1/sqrt(b.r*b.r+b.i*b.i);
			Real
				rr=(b.r*r+b.i*i)*inv_mag_b,
				ri=(b.r*i-b.i*r)*inv_mag_b;
			r=rr, i=ri;
			return *this;
		}
		Comp& operator/=(Real const &br){r/=br, i/=br; return *this;}
		Quat  operator/=(Quat const &b);
		Comp& operator^=(Comp const &b)
		{
			Comp t(log(sqrt(r*r+i*i)), atan2(i, r));
			t*=b;
			Real r0=exp(t.r), sin_ti, cos_ti;
			sin_cos(sin_ti, cos_ti, t.i);
			r=r0*cos_ti, i=r0*sin_ti;
			return *this;
		}
		Comp& operator^=(Real const &b)
		{
			Comp t(log(sqrt(r*r+i*i)), atan2(i, r));
			t*=b;
			Real r0=exp(t.r), sin_ti, cos_ti;
			sin_cos(sin_ti, cos_ti, t.i);
			r=r0*cos_ti, i=r0*sin_ti;
			return *this;
		}
	};
	__forceinline Real abs(Comp const &x){return sqrt(x.r*x.r+x.i*x.i);}
	__forceinline Comp operator*(Comp const &a, Comp const &b){return Comp(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r);}
	__forceinline Comp operator*(Comp const &a, Real const &br){return Comp(a.r*br, a.i*br);}
	__forceinline Comp operator*(Real const &a, Comp const &b){return Comp(a*b.r, a*b.i);}
	__forceinline Comp operator/(Comp const &a, Comp const &b)
	{
		Real _1_mag_b=1/(b.r*b.r+b.i*b.i);
		return Comp((b.r*a.r+b.i*a.i)*_1_mag_b, (b.r*a.i-b.i*a.r)*_1_mag_b);
	}
	__forceinline Comp operator/(Comp const &a, Real const &br){return Comp(a.r/br, a.i/br);}
	__forceinline Comp operator/(Real const &a, Comp const &b)
	{
		Real _ar_mag_b=a/(b.r*b.r+b.i*b.i);
		return Comp(b.r*_ar_mag_b, -b.i*_ar_mag_b);
	}
	__forceinline Comp operator+(Comp const &a, Comp const &b){return Comp(a.r+b.r, a.i+b.i);}
	__forceinline Comp operator+(Comp const &a, Real const &b){return Comp(a.r+b, a.i);}
	__forceinline Comp operator+(Real const &a, Comp const &b){return Comp(a+b.r, b.i);}
	__forceinline Comp operator-(Comp const &a, Comp const &b){return Comp(a.r-b.r, a.i-b.i);}
	__forceinline Comp operator-(Comp const &a, Real const &b){return Comp(a.r-b, a.i);}
	__forceinline Comp operator-(Real const &a, Comp const &b){return Comp(a-b.r, -b.i);}
	__forceinline Comp operator-(Comp const &a){return Comp(-a.r, -a.i);}
	__forceinline Real operator==(Comp const &a, Comp const &b){return (a.r==b.r)&(a.i==b.i);}
	__forceinline Real operator==(Comp const &a, Real const &b){return (a.r==b)&(a.i==0);}
	__forceinline Real operator==(Real const &a, Comp const &b){return (a==b.r)&(0==b.i);}
	__forceinline Real operator!=(Comp const &a, Comp const &b){return (a.r!=b.r)|(a.i!=b.i);}
	__forceinline Real operator!=(Comp const &a, Real const &b){return (a.r!=b)|(a.i!=0);}
	__forceinline Real operator!=(Real const &a, Comp const &b){return (a!=b.r)|(0!=b.i);}
	//__forceinline Comp operator|(Comp const &a, Comp const &b){return Comp(a.r|b.r, a.i|b.i);}
	//__forceinline Comp operator|(Comp const &a, Real const &b){return Comp(a.r|b, a.i);}
	//__forceinline Comp operator|(Real const &a, Comp const &b){return Comp(a|b.r, b.i);}
	//__forceinline Comp and(Comp const &a, Real const &b){return Comp(a.r&b, a.i&b);}
	//__forceinline Comp and(Real const &a, Comp const &b){return Comp(a&b.r, a&b.i);}
	inline Comp		floor		(Comp const &x){return Comp(mpfr::floor(x.r), mpfr::floor(x.i));}
	inline Comp		ceil		(Comp const &x){return Comp(mpfr::ceil(x.r), mpfr::ceil(x.i));}
	inline Comp		round		(Comp const &x){return Comp(mpfr::round(x.r), mpfr::round(x.i));}
	//__forceinline Comp operator%(Comp const &a, Comp const &b){return a-floor(a/b)*b;}
	//__forceinline Comp operator%(Comp const &a, Real const &b){return a-floor(a/b)*b;}
	//__forceinline Comp operator%(Real const &a, Comp const &b){return a-floor(a/b)*b;}
	__forceinline Comp log(Comp const &x){return Comp(mpfr::log(sqrt(x.r*x.r+x.i*x.i)), mpfr::atan2(x.i, x.r));}
	__forceinline Comp exp(Comp const &x)
	{
		Real sin_xi, cos_xi;
		mpfr::sin_cos(sin_xi, cos_xi, x.i);
		Real exp_xr=exp(x.r);
		return Comp(exp_xr*cos_xi, exp_xr*sin_xi);
	}
	__forceinline Comp sqrt(Comp const &x)
	{
		auto s=x.r+abs(x);
		if(s==0)
			return Comp(0, sqrt(-x.r)); 
		s=sqrt(s+s);
		return Comp(s*0.5, x.i/s);
	}
	__forceinline Comp operator^(Comp const &a, Comp const &b)//power
	{
		if(a.r==0&a.i==0&b.r==0&b.i==0)
			return Comp(1, 0);
		Comp t(log(sqrt(a.r*a.r+a.i*a.i)), atan2(a.i, a.r));
		t*=b;
		Real r0=exp(t.r);
		Real sin_ti, cos_ti;
		sin_cos(sin_ti, cos_ti, t.i);
		return Comp(r0*cos_ti, r0*sin_ti);
	}
	__forceinline Comp operator^(Comp const &a, Real const &b)
	{
		if(a.r==0&a.i==0&b==0)
			return Comp(1, 0);
		Comp t(log(sqrt(a.r*a.r+a.i*a.i)), atan2(a.i, a.r));
		t*=b;
		Real r0=exp(t.r);
		Real sin_ti, cos_ti;
		sin_cos(sin_ti, cos_ti, t.i);
		return Comp(r0*cos_ti, r0*sin_ti);
	}
	__forceinline Comp operator^(Real const &a, Comp const &b)
	{
		if(a==0&b.r==0&b.i==0)
			return Comp(1, 0);
		Comp t(log(abs(a)), atan2(0, a));
		t*=b;
		Real r0=exp(t.r);
		Real sin_ti, cos_ti;
		sin_cos(sin_ti, cos_ti, t.i);
		return Comp(r0*cos_ti, r0*sin_ti);
	}
	
	int			print_real(char *a, Real const &x, char base);//bases 2, 8, 10, 16 supported
	template<int _size>inline void print_unreal(char (&a)[_size], int &o, Real const &x, char comp, char base, bool &written)
	{
		if(x!=x)
			o+=sprintf_s(a+o, _size-o, written?"+0/0%c":"0/0%c", comp);
		else if(x.toDouble())
		{
			if(written&&x>0)
				o+=sprintf_s(a+o, _size-o, "+");
			if(x==1)
				o+=sprintf_s(a+o, _size-o, "%c", comp);
			else if(x==-1)
				o+=sprintf_s(a+o, _size-o, "-%c", comp);
			else
			{
				o+=print_real(a+o, x, base);
				o+=sprintf_s(a+o, _size-o, "%c", comp);
			}
			written=true;
		}
	}
	struct		Quat
	{
		Real r, i, j, k;
		Quat(Real const &r=0, Real const &i=0, Real const &j=0, Real const &k=0):r(r), i(i), j(j), k(k){}
		Quat(Real &&x):r(x), i(0), j(0), k(0){}
		Quat(Comp const &x):r(x.r), i(x.i), j(0), k(0){}
		Quat(Comp &&x):r(x.r), i(x.i), j(0), k(0){}
		Quat(Quat &&x):r(x.r), i(x.i), j(x.j), k(x.k){}
		void set_prec(int bin_prec, char mathSet)
		{
			r.set_prec(bin_prec);
			if(mathSet>='c')
			{
				i.set_prec(bin_prec);
				if(mathSet=='h')
				{
					j.set_prec(bin_prec);
					k.set_prec(bin_prec);
				}
			}
		}
		void set(double r=0, double i=0, double j=0, double k=0){this->r=r, this->i=i, this->j=j, this->k=k;}
		void setzero(){r.setZero(), i.setZero(), j.setZero(), k.setZero();}
		void setnan(){r.setNan(), i.setZero(), j.setZero(), k.setZero();}
		operator Real()const{return r;}
		operator Comp()const{return Comp(r, i);}
		bool r_is_true()const{return r!=0;}
		bool c_is_true()const{return r!=0||i!=0;}
		bool q_is_true()const{return r!=0||i!=0||j!=0||k!=0;}
		bool q_isTrue()const{return q_is_true();}//check g2.cpp
		Quat& operator=(Real const &x){r=x; return *this;}
		Quat& operator=(Comp const &x){r=x.r, i=x.i; return *this;}
		Quat& operator=(Real &&x){r=x; return *this;}
		Quat& operator=(Comp &&x){r=std::move(x.r), i=std::move(x.i); return *this;}
		Quat& operator+=(Quat const &b){r+=b.r, i+=b.i, j+=b.j, k+=b.k;return *this;}
		Quat& operator+=(Comp const &b){r+=b.r, i+=b.i;return *this;}
		Quat& operator+=(Real const &br){r+=br;return *this;}
		Quat& operator-=(Quat const &b){r-=b.r, i-=b.i, j-=b.j, k-=b.k;return *this;}
		Quat& operator-=(Comp const &b){r-=b.r, i-=b.i;return *this;}
		Quat& operator-=(Real const &br){r-=br;return *this;}
		Quat& operator*=(Quat const &b)
		{
			Real
				rr=r*b.r+i*b.i+j*b.j+k*b.k,
				ri=r*b.i+i*b.r+j*b.k-k*b.j,
				rj=r*b.j-i*b.k+j*b.r+k*b.i,
				rk=r*b.k+i*b.j+j*b.i+k*b.r;
			r=rr, i=ri, j=rj, k=rk;
			return *this;
		}
		Quat& operator*=(Comp const &b)
		{
			Real
				rr=r*b.r+i*b.i,
				ri=r*b.i+i*b.r,
				rj=j*b.r+k*b.i,
				rk=j*b.i+k*b.r;
			r=rr, i=ri, j=rj, k=rk;
			return *this;
		}
		Quat& operator*=(Real const &b){r*=b, i*=b, j*=b, k*=b;return *this;}
		Quat& operator/=(Quat const &b)
		{
			Real _1_mag_y=1/sqrt(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
			Real
				rr=(b.r*r+b.i*i+b.j*j+b.k*k)*_1_mag_y,
				ri=(b.r*i-b.i*r-b.j*k+b.k*j)*_1_mag_y,
				rj=(b.r*j+b.i*k-b.j*r-b.k*i)*_1_mag_y,
				rk=(b.r*k-b.i*j+b.j*i-b.k*r)*_1_mag_y;
			r=rr, i=ri, j=rj, k=rk;
			return *this;
		}
		Quat& operator/=(Comp const &b)
		{
			Real inv_mag_y=1/sqrt(b.r*b.r+b.i*b.i);
			Real
				rr=(b.r*r+b.i*i)*inv_mag_y,
				ri=(b.r*i-b.i*r)*inv_mag_y,
				rj=(b.r*j+b.i*k)*inv_mag_y,
				rk=(b.r*k-b.i*j)*inv_mag_y;
			r=rr, i=ri, j=rj, k=rk;
			return *this;
		}
		Quat& operator/=(Real const &b){r/=b, i/=b, j/=b, k/=b;return *this;}
		Quat& operator^=(Quat const &b)
		{
			Real mag_v=i*i+j*j+k*k;
			Real ln_mag_a=log(sqrt(r*r+mag_v));
			mag_v=sqrt(mag_v);
			Real v_mul=acos((r/ln_mag_a));
			v_mul/=mag_v;
			Quat t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
			t=Quat(
				b.r*t.r+b.i*t.i+b.j*t.j+b.k*t.k,
				b.r*t.i+b.i*t.r+b.j*t.k-b.k*t.j,
				b.r*t.j-b.i*t.k+b.j*t.r+b.k*t.i,
				b.r*t.k+b.i*t.j+b.j*t.i+b.k*t.r);
			Real mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k), sin_mu, cos_mu;
			sin_cos(sin_mu, cos_mu, mag_u);
			Real exp_tr=exp(t.r);
			v_mul=exp_tr*sin_mu/mag_u;
			r=exp_tr*cos_mu, i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
			return *this;
		}
		Quat& operator^=(Comp const &b)
		{
			Real mag_v=i*i+j*j+k*k;
			Real ln_mag_a=log(sqrt(r*r+mag_v));
			mag_v=sqrt(mag_v);
			Real v_mul=acos((r/ln_mag_a));
			v_mul/=mag_v;
			Quat t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
			t=Quat(
				b.r*t.r+b.i*t.i,
				b.r*t.i+b.i*t.r,
				b.r*t.j-b.i*t.k,
				b.r*t.k+b.i*t.j);
			Real mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
			Real sin_mu, cos_mu;
			sin_cos(sin_mu, cos_mu, mag_u);
			Real exp_tr=exp(t.r);
			v_mul=exp_tr*sin_mu/mag_u;
			r=exp_tr*cos_mu, i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
			return *this;
		}
		Quat& operator^=(Real const &br)
		{
			Real mag_v=i*i+j*j+k*k;
			Real ln_mag_a=log(sqrt(r*r+mag_v));
			mag_v=sqrt(mag_v);
			Real v_mul=acos(r/ln_mag_a);
			v_mul/=mag_v;
			Quat t(ln_mag_a, i*v_mul, j*v_mul, k*v_mul);
			t=Quat(br*t.r, br*t.i, br*t.j, br*t.k);
			Real mag_u=sqrt(t.i*t.i+t.j*t.j+t.k*t.k);
			Real sin_mu, cos_mu;
			sin_cos(sin_mu, cos_mu, mag_u);
			Real exp_tr=exp(t.r);
			v_mul=exp_tr*sin_mu/mag_u;
			r=exp_tr*cos_mu, i=t.i*v_mul, j=t.j*v_mul, k=t.k*v_mul;
			return *this;
		}
		template<int _size>void print(char (&a)[_size], int &o, char mathSet, int base=10)const
		{
			switch(mathSet)
			{
			case 'R':
				o+=print_real(a+o, r, base);
				break;
			case 'c':
				{
					bool written=false;
					if(r.toDouble())
						o+=print_real(a+o, r, base), written=true;
					print_unreal(a, o, i, 'i', base, written);
					if(!written)
						o+=sprintf_s(a+o, _size-o, "0");
				}
				break;
			case 'h':
				{
					bool written=false;
					if(r.toDouble())
						o+=print_real(a+o, r, base), written=true;
					print_unreal(a, o, i, 'i', base, written);
					print_unreal(a, o, j, 'j', base, written);
					print_unreal(a, o, k, 'k', base, written);
					if(!written)
						o+=sprintf_s(a+o, _size-o, "0");
				}
				break;
			}
			//switch(mathSet)
			//{
			//case 'R':
			//	{
			//		const char *sign[]={"-"}, *base_prefix=nullptr;
			//		bool sbit=signbit(r);
			//		if(base==10)
			//			o+=mpfr_sprintf(a+o, "%.*Rg", dec_prec, abs(r).mpfr_srcptr());
			//		else//2, 8, 16
			//		{
			//			bool negative=mpfr::signbit(r);
			//			long exponent=0;
			//			char *str=mpfr_get_str(0, &exponent, base, 0, abs(r).mpfr_srcptr(), MPFR_RNDNA);
			//			std::string str2=str;
			//			mpfr_free_str(str);
			//		}
			//	}
			//	break;
			//case 'c':
			//	{
			//		bool written=false;
			//		if(r)
			//			o+=mpfr_sprintf(a+o, "%.*Rg", dec_prec, r.mpfr_srcptr()), written=true;
			//	}
			//	break;
			//case 'h':
			//	printQuaternion(a, o);
			//	break;
			//}
		}
	};
	inline Real abs(Quat const &x){return sqrt(x.r*x.r+x.i*x.i+x.j*x.j+x.k*x.k);}

	__forceinline bool operator==(Quat const &a, Quat const &b){return (a.r==b.r)&(a.i==b.i)&(a.j==b.j)&(a.k==b.k);}
	__forceinline bool operator==(Quat const &a, Comp const &b){return (a.r==b.r)&(a.i==b.i)&(a.j==0)&(a.k==0);}
	__forceinline bool operator==(Quat const &a, Real const &b){return (a.r==b)&(a.i==0)&(a.j==0)&(a.k==0);}
	__forceinline bool operator==(Comp const &a, Quat const &b){return (a.r==b.r)&(a.i==b.i)&(0==b.j)&(0==b.k);}
	__forceinline bool operator==(Real const &a, Quat const &b){return (a==b.r)&(0==b.i)&(0==b.j)&(0==b.k);}

	__forceinline bool operator!=(Quat const &a, Quat const &b){return (a.r!=b.r)|(a.i!=b.i)|(a.j!=b.j)|(a.k!=b.k);}
	__forceinline bool operator!=(Quat const &a, Comp const &b){return (a.r!=b.r)|(a.i!=b.i)|(a.j!=0)|(a.k!=0);}
	__forceinline bool operator!=(Quat const &a, Real const &b){return (a.r!=b)|(a.i!=0)|(a.j!=0)|(a.k!=0);}
	__forceinline bool operator!=(Comp const &a, Quat const &b){return (a.r!=b.r)|(a.i!=b.i)|(0!=b.j)|(0!=b.k);}
	__forceinline bool operator!=(Real const &a, Quat const &b){return (a!=b.r)|(0!=b.i)|(0!=b.j)|(0!=b.k);}

	__forceinline Quat operator+(Quat const &a, Quat const &b){return Quat(a.r+b.r, a.i+b.i, a.j+b.j, a.k+b.k);}
	__forceinline Quat operator+(Quat const &a, Comp const &b){return Quat(a.r+b.r, a.i+b.i, a.j, a.k);}
	__forceinline Quat operator+(Quat const &a, Real const &b){return Quat(a.r+b, a.i, a.j, a.k);}
	__forceinline Quat operator+(Comp const &a, Quat const &b){return Quat(a.r+b.r, a.i+b.i, b.j, b.k);}
	__forceinline Quat operator+(Real const &a, Quat const &b){return Quat(a+b.r, b.i, b.j, b.k);}

	__forceinline Quat operator-(Quat const &a, Quat const &b){return Quat(a.r-b.r, a.i-b.i, a.j-b.j, a.k-b.k);}
	__forceinline Quat operator-(Quat const &a, Comp const &b){return Quat(a.r-b.r, a.i-b.i, a.j, a.k);}
	__forceinline Quat operator-(Quat const &a, Real const &b){return Quat(a.r-b, a.i, a.j, a.k);}
	__forceinline Quat operator-(Comp const &a, Quat const &b){return Quat(a.r-b.r, a.i-b.i, -b.j, -b.k);}
	__forceinline Quat operator-(Real const &a, Quat const &b){return Quat(a-b.r, -b.i, -b.j, -b.k);}
	__forceinline Quat operator-(Quat const &a){return Quat(-a.r, -a.i, -a.j, -a.k);}
	__forceinline Quat operator*(Quat const &a, Quat const &b)
	{
		return Quat(
			a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,
			a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,
			a.r*b.j-a.i*b.k+a.j*b.r+a.k*b.i,
			a.r*b.k+a.i*b.j-a.j*b.i+a.k*b.r);
	}
	__forceinline Quat operator*(Quat const &a, Comp const &b)
	{
		return Quat(
			a.r*b.r-a.i*b.i,
			a.r*b.i+a.i*b.r,
			a.j*b.r+a.k*b.i,
			-a.j*b.i+a.k*b.r);
	}
	__forceinline Quat operator*(Quat const &a, Real const &b){return Quat(a.r*b, a.i*b, a.j*b, a.k*b);}
	__forceinline Quat operator*(Comp const &a, Quat const &b)
	{
		return Quat(
			a.r*b.r-a.i*b.i,
			a.r*b.i+a.i*b.r,
			a.r*b.j-a.i*b.k,
			a.r*b.k+a.i*b.j);
	}
	__forceinline Quat operator*(Real const &a, Quat const &b){return Quat(a*b.r, a*b.i, a*b.j, a*b.k);}
	__forceinline Quat operator/(Quat const &a, Quat const &b)
	{
		Real inv_mag_y=1/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
		return Quat(
			(b.r*a.r+b.i*a.i+b.j*a.j+b.k*a.k)*inv_mag_y,
			(b.r*a.i-b.i*a.r-b.j*a.k+b.k*a.j)*inv_mag_y,
			(b.r*a.j+b.i*a.k-b.j*a.r-b.k*a.i)*inv_mag_y,
			(b.r*a.k-b.i*a.j+b.j*a.i-b.k*a.r)*inv_mag_y);
	}
	__forceinline Quat Comp::operator/=(Quat const &b)
	{
		Real inv_mag_y=1/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
		return Quat(
			( b.r*r+b.i*i)*inv_mag_y,
			( b.r*i-b.i*r)*inv_mag_y,
			(-b.j*r-b.k*i)*inv_mag_y,
			( b.j*i-b.k*r)*inv_mag_y);
	}
	__forceinline Quat operator/(Quat const &a, Comp const &b)
	{
		Real inv_mag_y=1/(b.r*b.r+b.i*b.i);
		return Quat(
			(b.r*a.r+b.i*a.i)*inv_mag_y,
			(b.r*a.i-b.i*a.r)*inv_mag_y,
			(b.r*a.j+b.i*a.k)*inv_mag_y,
			(b.r*a.k-b.i*a.j)*inv_mag_y);
	}
	__forceinline Quat operator/(Quat const &a, Real const &b)
	{
		Real inv_mag_y=1/b;
		return Quat(a.r*inv_mag_y, a.i*inv_mag_y, a.j*inv_mag_y, a.k*inv_mag_y);
	}
	__forceinline Quat operator/(Comp const &a, Quat const &b)
	{
		Real inv_mag_y=1/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
		return Quat(
			(b.r*a.r+b.i*a.i)*inv_mag_y,
			(b.r*a.i-b.i*a.r)*inv_mag_y,
			(-b.j*a.r-b.k*a.i)*inv_mag_y,
			(b.j*a.i-b.k*a.r)*inv_mag_y);
	}
	__forceinline Quat operator/(Real const &a, Quat const &b)
	{
		Real _ar_mag_y=a/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
		return Quat(b.r*_ar_mag_y, -b.i*_ar_mag_y, -b.j*_ar_mag_y, -b.k*_ar_mag_y);
	}
	__forceinline Quat log(Quat const &x)
	{
		Real mag_v=x.i*x.i+x.j*x.j+x.k*x.k;

		if(mag_v==0)
			return Quat(log(abs(x.r)), m_pi, 0, 0);

		Real mag_x=sqrt(x.r*x.r+mag_v);
		mag_v=sqrt(mag_v);
		Real u_mul=acos(x.r/mag_x)/mag_v;
		Real ln_mag_x=log(mag_x);
		return Quat(ln_mag_x, x.i*u_mul, x.j*u_mul, x.k*u_mul);
	}
	__forceinline Quat exp(Quat const &x)
	{
		Real exp_r=exp(x.r);
		Real mag_v=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
		Real sin_v, cos_v;
		sin_cos(sin_v, cos_v, mag_v);
		Real v_mul=exp_r*sin_v/mag_v;
		return Quat(exp_r*cos_v, x.i*v_mul, x.j*v_mul, x.k*v_mul);
	}
	__forceinline Quat sqrt(Quat const &x)
	{
		auto s=x.r+abs(x);
		if(s==0)
			return Quat(0, sqrt(-x.r), 0, 0);
		s=sqrt(s+s);
		auto i=sqrt(-x.r);
		auto inv_s=1/s;
		return Quat(s*0.5, x.i*inv_s, x.j*inv_s, x.k*inv_s);
	}
	__forceinline Quat operator^(Quat const &a, Quat const &b)
	{
		if(a.r==0&a.i==0&a.j==0&a.k==0&b.r==0&b.i==0&b.j==0&b.k==0)
			return Quat(1, 0, 0, 0);
		return exp(b*log(a));
	}
	__forceinline Quat operator^(Quat const &a, Comp const &b)
	{
		if(a.r==0&a.i==0&a.j==0&a.k==0&b.r==0&b.i==0)
			return Quat(1, 0, 0, 0);
		return exp(b*log(a));
	}
	__forceinline Quat operator^(Quat const &a, Real const &b)
	{
		if(a.r==0&a.i==0&a.j==0&a.k==0&b==0)
			return Quat(1, 0, 0, 0);
		return exp(b*log(a));
	}
	__forceinline Quat operator^(Comp const &a, Quat const &b)
	{
		if(a.r==0&a.i==0&b.r==0&b.i==0&b.j==0&b.k==0)
			return Quat(1, 0, 0, 0);
		return exp(b*log(a));
	}
	__forceinline Quat operator^(Real const &a, Quat const &b)
	{
		if(a==0&b.r==0&b.i==0&b.j==0&b.k==0)
			return Quat(1, 0, 0, 0);
		return exp(b*log(a));
	}
	//__forceinline Quat operator|(Quat const &a, Quat const &b){return Quat(a.r|b.r, a.i|b.i, a.j|b.j, a.k|b.k);}
	//__forceinline Quat operator|(Quat const &a, Comp const &b){return Quat(a.r|b.r, a.i|b.i, a.j, a.k);}
	//__forceinline Quat operator|(Quat const &a, Real const &b){return Quat(a.r|b, a.i, a.j, a.k);}
	//__forceinline Quat operator|(Comp const &a, Quat const &b){return Quat(a.r|b.r, a.i|b.i, b.j, b.k);}
	//__forceinline Quat operator|(Real const &a, Quat const &b){return Quat(a|b.r, b.i, b.j, b.k);}
	//__forceinline Quat and(Quat2d const &a, Vect2d const &b){return Quat2d(a.r&b, a.i&b, a.j&b, a.k&b);}
	//__forceinline Quat and(Vect2d const &a, Quat2d const &b){return Quat2d(a&b.r, a&b.i, a&b.j, a&b.k);}
	//__forceinline Quat operator%(Quat const &a, Quat const &b){return a-floor(a/b)*b;}
	//__forceinline Quat operator%(Quat const &a, Comp const &b){return a-floor(a/b)*b;}
	//__forceinline Quat operator%(Quat const &a, Real const &b){return a-floor(a/b)*b;}
	//__forceinline Quat operator%(Comp const &a, Quat const &b){return a-floor(a/b)*b;}
	//__forceinline Quat operator%(Real const &a, Quat const &b){return a-floor(a/b)*b;}

	inline Quat		floor		(Quat const &x){return Quat(mpfr::floor(x.r), mpfr::floor(x.i), mpfr::floor(x.j), mpfr::floor(x.k));}
	inline Quat		ceil		(Quat const &x){return Quat(mpfr::ceil(x.r), mpfr::ceil(x.i), mpfr::ceil(x.j), mpfr::ceil(x.k));}
	inline Quat		round		(Quat const &x){return Quat(mpfr::round(x.r), mpfr::round(x.i), mpfr::round(x.j), mpfr::round(x.k));}
	inline Real		operator%	(Real const &x, Real const &y){return x-floor(x/y)*y;}
	inline Comp		operator%	(Real const &x, Comp const &y){return x-floor(x/y)*y;}
	inline Quat		operator%	(Real const &x, Quat const &y){return x-floor(x/y)*y;}
	inline Comp		operator%	(Comp const &x, Real const &y){return x-floor(x/y)*y;}
	inline Comp		operator%	(Comp const &x, Comp const &y){return x-floor(x/y)*y;}
	inline Quat		operator%	(Comp const &x, Quat const &y){return x-floor(x/y)*y;}
	inline Quat		operator%	(Quat const &x, Real const &y){return x-floor(x/y)*y;}
	inline Quat		operator%	(Quat const &x, Comp const &y){return x-floor(x/y)*y;}
	inline Quat		operator%	(Quat const &x, Quat const &y){return x-floor(x/y)*y;}
#endif

#ifdef MP_PURE_QUAT
	void  r_r_setzero				(Quat &r, Quat const&);
	void  c_c_setzero				(Quat &r, Quat const&);
	void  q_q_setzero				(Quat &r, Quat const&);
	
	void  r_r_ceil					(Quat &r, Quat const &x);
	void  c_c_ceil					(Quat &r, Quat const &x);
	void  q_q_ceil					(Quat &r, Quat const &x);

	void  r_r_floor					(Quat &r, Quat const &x);
	void  c_c_floor					(Quat &r, Quat const &x);
	void  q_q_floor					(Quat &r, Quat const &x);

	void  r_r_round					(Quat &r, Quat const &x);
	void  c_c_round					(Quat &r, Quat const &x);
	void  q_q_round					(Quat &r, Quat const &x);

	void  r_r_int					(Quat &r, Quat const &x);
	void  c_c_int					(Quat &r, Quat const &x);
	void  q_q_int					(Quat &r, Quat const &x);

	void  r_r_frac					(Quat &r, Quat const &x);
	void  c_c_frac					(Quat &r, Quat const &x);
	void  q_q_frac					(Quat &r, Quat const &x);

	void  r_r_abs					(Quat &r, Quat const &x);
	void  r_c_abs					(Quat &r, Quat const &x);
	void  r_q_abs					(Quat &r, Quat const &x);

	void  r_r_arg					(Quat &r, Quat const &x);
	void  r_c_arg					(Quat &r, Quat const &x);
	void  r_q_arg					(Quat &r, Quat const &x);

	void  r_c_real					(Quat &r, Quat const &x);

	void  r_c_imag					(Quat &r, Quat const &x);

	//r_conjugate: assign
	void c_c_conjugate				(Quat &r, Quat const &x);
	void q_q_conjugate				(Quat &r, Quat const &x);

	void  c_r_polar					(Quat &r, Quat const &x);
	void  c_c_polar					(Quat &r, Quat const &x);
	void  c_q_polar					(Quat &r, Quat const &x);

	//r_cartesian	assign
	void  c_c_cartesian				(Quat &r, Quat const &x);
	void  q_q_cartesian				(Quat &r, Quat const &x);

	void r_rr_plus					(Quat &r, Quat const &x, Quat const &y);
	void c_rc_plus					(Quat &r, Quat const &x, Quat const &y);
	void q_rq_plus					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_plus					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_plus					(Quat &r, Quat const &x, Quat const &y);
	void q_cq_plus					(Quat &r, Quat const &x, Quat const &y);
	void q_qr_plus					(Quat &r, Quat const &x, Quat const &y);
	void q_qc_plus					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_plus					(Quat &r, Quat const &x, Quat const &y);

	void  r_r_minus					(Quat &r, Quat const &x);
	void  c_c_minus					(Quat &r, Quat const &x);
	void  q_q_minus					(Quat &r, Quat const &x);
	void r_rr_minus					(Quat &r, Quat const &x, Quat const &y);
	void c_rc_minus					(Quat &r, Quat const &x, Quat const &y);
	void q_rq_minus					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_minus					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_minus					(Quat &r, Quat const &x, Quat const &y);
	void q_cq_minus					(Quat &r, Quat const &x, Quat const &y);
	void q_qr_minus					(Quat &r, Quat const &x, Quat const &y);
	void q_qc_minus					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_minus					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_multiply				(Quat &r, Quat const &x, Quat const &y);
	void c_rc_multiply				(Quat &r, Quat const &x, Quat const &y);
	void q_rq_multiply				(Quat &r, Quat const &x, Quat const &y);
	void c_cr_multiply				(Quat &r, Quat const &x, Quat const &y);
	void c_cc_multiply				(Quat &r, Quat const &x, Quat const &y);
	void q_cq_multiply				(Quat &r, Quat const &x, Quat const &y);
	void q_qr_multiply				(Quat &r, Quat const &x, Quat const &y);
	void q_qc_multiply				(Quat &r, Quat const &x, Quat const &y);
	void q_qq_multiply				(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_divide				(Quat &r, Quat const &y);
	void  c_c_divide				(Quat &r, Quat const &y);
	void  q_q_divide				(Quat &r, Quat const &y);
	void r_rr_divide				(Quat &r, Quat const &x, Quat const &y);
	void c_rc_divide				(Quat &r, Quat const &x, Quat const &y);
	void q_rq_divide				(Quat &r, Quat const &x, Quat const &y);
	void c_cr_divide				(Quat &r, Quat const &x, Quat const &y);
	void c_cc_divide				(Quat &r, Quat const &x, Quat const &y);
	void q_cq_divide				(Quat &r, Quat const &x, Quat const &y);
	void q_qr_divide				(Quat &r, Quat const &x, Quat const &y);
	void q_qc_divide				(Quat &r, Quat const &x, Quat const &y);
	void q_qq_divide				(Quat &r, Quat const &x, Quat const &y);
	
	void r_rr_logic_divides			(Quat &r, Quat const &y, Quat const &x);//rc_divides, rq_divides: applied to each component
	void r_rc_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_rq_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_cr_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_cc_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_cq_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_qr_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_qc_logic_divides			(Quat &r, Quat const &y, Quat const &x);
	void r_qq_logic_divides			(Quat &r, Quat const &y, Quat const &x);

	void r_rr_power_real			(Quat &r, Quat const &x, Quat const &y);//trunc
	void c_cr_power_real			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_power_real			(Quat &r, Quat const &x, Quat const &y);

	void c_cr_pow					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_pow					(Quat &r, Quat const &x, Quat const &y);
	void q_cq_pow					(Quat &r, Quat const &x, Quat const &y);
	void q_qr_pow					(Quat &r, Quat const &x, Quat const &y);
	void q_qc_pow					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_pow					(Quat &r, Quat const &x, Quat const &y); 

	void  c_c_ln					(Quat &r, Quat const &x);
	void  q_q_ln					(Quat &r, Quat const &x);
	
	void  c_c_log					(Quat &r, Quat const &x);
	void  q_q_log					(Quat &r, Quat const &x);
	void c_cr_log					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_log					(Quat &r, Quat const &x, Quat const &y);
	void q_cq_log					(Quat &r, Quat const &x, Quat const &y);
	void q_qc_log					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_log					(Quat &r, Quat const &x, Quat const &y);
	
	void c_rr_tetrate				(Quat &r, Quat const &x, Quat const &y);
	void c_rc_tetrate				(Quat &r, Quat const &x, Quat const &y);
	void c_cr_tetrate				(Quat &r, Quat const &x, Quat const &y);
	void c_cc_tetrate				(Quat &r, Quat const &x, Quat const &y);
	void q_qr_tetrate				(Quat &r, Quat const &x, Quat const &y);
	
	void c_rr_pentate				(Quat &r, Quat const &x, Quat const &y);
	void c_cr_pentate				(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_shift_left_l	(Quat &r, Quat const &x);	//<<x = 2^floor(x)
	void  c_c_bitwise_shift_left_l	(Quat &r, Quat const &x);
	void  q_q_bitwise_shift_left_l	(Quat &r, Quat const &x);
	void  r_r_bitwise_shift_left_r	(Quat &r, Quat const &x);	//x<< = 2x
	void  c_c_bitwise_shift_left_r	(Quat &r, Quat const &x);
	void  q_q_bitwise_shift_left_r	(Quat &r, Quat const &x);
	void r_rr_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);//x<<y = x*2^y
	void c_rc_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_shift_right_l	(Quat &r, Quat const &x);
	void  c_c_bitwise_shift_right_l	(Quat &r, Quat const &x);
	void  q_q_bitwise_shift_right_l	(Quat &r, Quat const &x);
	void  r_r_bitwise_shift_right_r	(Quat &r, Quat const &x);
	void  c_c_bitwise_shift_right_r	(Quat &r, Quat const &x);
	void  q_q_bitwise_shift_right_r	(Quat &r, Quat const &x);
	void r_rr_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_not			(Quat &r, Quat const &x);
	void  c_c_bitwise_not			(Quat &r, Quat const &x);
	void  q_q_bitwise_not			(Quat &r, Quat const &x);

	void  r_r_bitwise_and			(Quat &r, Quat const &x);
	void  c_c_bitwise_and			(Quat &r, Quat const &x);
	void  q_q_bitwise_and			(Quat &r, Quat const &x);
	void r_rr_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_and			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_and			(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_nand			(Quat &r, Quat const &x);
	void  c_c_bitwise_nand			(Quat &r, Quat const &x);
	void  q_q_bitwise_nand			(Quat &r, Quat const &x);
	void r_rr_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_or			(Quat &r, Quat const &x);
	void  c_c_bitwise_or			(Quat &r, Quat const &x);
	void  q_q_bitwise_or			(Quat &r, Quat const &x);
	void r_rr_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_nor			(Quat &r, Quat const &x);
	void  c_c_bitwise_nor			(Quat &r, Quat const &x);
	void  q_q_bitwise_nor			(Quat &r, Quat const &x);
	void r_rr_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_xor			(Quat &r, Quat const &x);
	void  c_c_bitwise_xor			(Quat &r, Quat const &x);
	void  q_q_bitwise_xor			(Quat &r, Quat const &x);
	void r_rr_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_xnor			(Quat &r, Quat const &x);
	void  c_c_bitwise_xnor			(Quat &r, Quat const &x);
	void  q_q_bitwise_xnor			(Quat &r, Quat const &x);
	void r_rr_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void c_rc_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void q_rq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void q_cq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void q_qr_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void q_qc_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_equal			(Quat &r, Quat const &x);
	void  r_c_logic_equal			(Quat &r, Quat const &x);
	void  r_q_logic_equal			(Quat &r, Quat const &x);
	void r_rr_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_equal			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_not_equal		(Quat &r, Quat const &x);
	void  r_c_logic_not_equal		(Quat &r, Quat const &x);
	void  r_q_logic_not_equal		(Quat &r, Quat const &x);
	void r_rr_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_not_equal		(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_less_l			(Quat &r, Quat const &x);
	void  r_c_logic_less_l			(Quat &r, Quat const &x);
	void  r_q_logic_less_l			(Quat &r, Quat const &x);
	void  r_r_logic_less_r			(Quat &r, Quat const &x);
	void  r_c_logic_less_r			(Quat &r, Quat const &x);
	void  r_q_logic_less_r			(Quat &r, Quat const &x);
	void r_rr_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_less			(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_less			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_less_equal_l	(Quat &r, Quat const &x);
	void  r_c_logic_less_equal_l	(Quat &r, Quat const &x);
	void  r_q_logic_less_equal_l	(Quat &r, Quat const &x);
	void  r_r_logic_less_equal_r	(Quat &r, Quat const &x);
	void  r_c_logic_less_equal_r	(Quat &r, Quat const &x);
	void  r_q_logic_less_equal_r	(Quat &r, Quat const &x);
	void r_rr_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_less_equal		(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_greater_l		(Quat &r, Quat const &x);
	void  r_c_logic_greater_l		(Quat &r, Quat const &x);
	void  r_q_logic_greater_l		(Quat &r, Quat const &x);
	void  r_r_logic_greater_r		(Quat &r, Quat const &x);
	void  r_c_logic_greater_r		(Quat &r, Quat const &x);
	void  r_q_logic_greater_r		(Quat &r, Quat const &x);
	void r_rr_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_greater			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_greater_equal_l	(Quat &r, Quat const &x);
	void  r_c_logic_greater_equal_l	(Quat &r, Quat const &x);
	void  r_q_logic_greater_equal_l	(Quat &r, Quat const &x);
	void  r_r_logic_greater_equal_r	(Quat &r, Quat const &x);
	void  r_c_logic_greater_equal_r	(Quat &r, Quat const &x);
	void  r_q_logic_greater_equal_r	(Quat &r, Quat const &x);
	void r_rr_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_not				(Quat &r, Quat const &x);
	void  r_c_logic_not				(Quat &r, Quat const &x);
	void  r_q_logic_not				(Quat &r, Quat const &x);
	
	void r_rr_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_and				(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_and				(Quat &r, Quat const &x, Quat const &y);
	
	void r_rr_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_or				(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_or				(Quat &r, Quat const &x, Quat const &y);

	void r_rr_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_rc_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_rq_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_cr_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_cc_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_cq_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_qr_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_qc_logic_xor				(Quat &r, Quat const &x, Quat const &y);
	void r_qq_logic_xor				(Quat &r, Quat const &x, Quat const &y);

	void r_rr_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void c_rc_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void q_rq_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void c_cr_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void c_cc_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void q_cq_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void q_qr_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void q_qc_condition_zero		(Quat &r, Quat const &x, Quat const &y);
	void q_qq_condition_zero		(Quat &r, Quat const &x, Quat const &y);

	void  r_r_percent				(Quat &r, Quat const &x);
	void  c_c_percent				(Quat &r, Quat const &x);
	void  q_q_percent				(Quat &r, Quat const &x);
	
	void r_rr_modulo				(Quat &r, Quat const &x, Quat const &y);
	void c_rc_modulo				(Quat &r, Quat const &x, Quat const &y);
	void q_rq_modulo				(Quat &r, Quat const &x, Quat const &y);
	void c_cr_modulo				(Quat &r, Quat const &x, Quat const &y);
	void c_cc_modulo				(Quat &r, Quat const &x, Quat const &y);
	void q_cq_modulo				(Quat &r, Quat const &x, Quat const &y);
	void q_qr_modulo				(Quat &r, Quat const &x, Quat const &y);
	void q_qc_modulo				(Quat &r, Quat const &x, Quat const &y);
	void q_qq_modulo				(Quat &r, Quat const &x, Quat const &y);

	void  r_r_sgn					(Quat &r, Quat const &x);
	void  c_c_sgn					(Quat &r, Quat const &x);
	void  q_q_sgn					(Quat &r, Quat const &x);
	
	void  r_r_sq					(Quat &r, Quat const &x);
	void  c_c_sq					(Quat &r, Quat const &x);
	void  q_q_sq					(Quat &r, Quat const &x);
	
	void  c_c_sqrt					(Quat &r, Quat const &x);
	void  q_q_sqrt					(Quat &r, Quat const &x);

	void  r_r_invsqrt				(Quat &r, Quat const &x);

	void  r_r_cbrt					(Quat &r, Quat const &x);
	void  c_c_cbrt					(Quat &r, Quat const &x);
	void  q_q_cbrt					(Quat &r, Quat const &x);

	void  r_r_gauss					(Quat &r, Quat const &x);
	void  c_c_gauss					(Quat &r, Quat const &x);
	void  q_q_gauss					(Quat &r, Quat const &x);

	void  r_r_erf					(Quat &r, Quat const &x);

	void  r_r_zeta					(Quat &r, Quat const &x);
	
	void  r_r_tgamma				(Quat &r, Quat const &x);
	void  c_c_tgamma				(Quat &r, Quat const &x);
	void  q_q_tgamma				(Quat &r, Quat const &x);
	void r_rr_tgamma				(Quat &r, Quat const &x, Quat const &y);

	void  r_r_loggamma				(Quat &r, Quat const &x);

	void  r_r_factorial				(Quat &r, Quat const &x);
	void  c_c_factorial				(Quat &r, Quat const &x);
	void  q_q_factorial				(Quat &r, Quat const &x);

	void  r_r_permutation			(Quat &r, Quat const &x);
	void  c_c_permutation			(Quat &r, Quat const &x);
	void  q_q_permutation			(Quat &r, Quat const &x);
	void r_rr_permutation			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_permutation			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_permutation			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_permutation			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_combination			(Quat &r, Quat const &x);
	void  c_c_combination			(Quat &r, Quat const &x);
	void  q_q_combination			(Quat &r, Quat const &x);
	void r_rr_combination			(Quat &r, Quat const &x, Quat const &y);
	void c_cr_combination			(Quat &r, Quat const &x, Quat const &y);
	void c_cc_combination			(Quat &r, Quat const &x, Quat const &y);
	void q_qq_combination			(Quat &r, Quat const &x, Quat const &y);

	void  r_r_cos					(Quat &r, Quat const &x);
	void  c_c_cos					(Quat &r, Quat const &x);
	void  q_q_cos					(Quat &r, Quat const &x);

	void  c_c_acos					(Quat &r, Quat const &x);
	void  q_q_acos					(Quat &r, Quat const &x);

	void  r_r_cosh					(Quat &r, Quat const &x);
	void  c_c_cosh					(Quat &r, Quat const &x);
	void  q_q_cosh					(Quat &r, Quat const &x);
	
	void  c_c_acosh					(Quat &r, Quat const &x);
	void  q_q_acosh					(Quat &r, Quat const &x);

	void  r_r_cosc					(Quat &r, Quat const &x);
	void  c_c_cosc					(Quat &r, Quat const &x);
	void  q_q_cosc					(Quat &r, Quat const &x);

	void  r_r_sec					(Quat &r, Quat const &x);
	void  c_c_sec					(Quat &r, Quat const &x);
	void  q_q_sec					(Quat &r, Quat const &x);

	void  c_c_asec					(Quat &r, Quat const &x);
	void  q_q_asec					(Quat &r, Quat const &x);

	void  r_r_sech					(Quat &r, Quat const &x);
	void  c_c_sech					(Quat &r, Quat const &x);
	void  q_q_sech					(Quat &r, Quat const &x);

	void  c_c_asech					(Quat &r, Quat const &x);
	void  q_q_asech					(Quat &r, Quat const &x);

	void  r_r_sin					(Quat &r, Quat const &x);
	void  c_c_sin					(Quat &r, Quat const &x);
	void  q_q_sin					(Quat &r, Quat const &x);

	void  c_c_asin					(Quat &r, Quat const &x);
	void  q_q_asin					(Quat &r, Quat const &x);

	void  r_r_sinh					(Quat &r, Quat const &x);
	void  c_c_sinh					(Quat &r, Quat const &x);
	void  q_q_sinh					(Quat &r, Quat const &x);

	void  r_r_asinh					(Quat &r, Quat const &x);
	void  c_c_asinh					(Quat &r, Quat const &x);
	void  q_q_asinh					(Quat &r, Quat const &x);

	void  r_r_sinc					(Quat &r, Quat const &x);
	void  c_c_sinc					(Quat &r, Quat const &x);
	void  q_q_sinc					(Quat &r, Quat const &x);

	void  r_r_sinhc					(Quat &r, Quat const &x);
	void  c_c_sinhc					(Quat &r, Quat const &x);
	void  q_q_sinhc					(Quat &r, Quat const &x);

	void  r_r_csc					(Quat &r, Quat const &x);
	void  c_c_csc					(Quat &r, Quat const &x);
	void  q_q_csc					(Quat &r, Quat const &x);

	void  c_c_acsc					(Quat &r, Quat const &x);
	void  q_q_acsc					(Quat &r, Quat const &x);

	void  r_r_csch					(Quat &r, Quat const &x);
	void  c_c_csch					(Quat &r, Quat const &x);
	void  q_q_csch					(Quat &r, Quat const &x);

	void  r_r_acsch					(Quat &r, Quat const &x);
	void  c_c_acsch					(Quat &r, Quat const &x);
	void  q_q_acsch					(Quat &r, Quat const &x);

	void  r_r_tan					(Quat &r, Quat const &x);
	void  c_c_tan					(Quat &r, Quat const &x);
	void  q_q_tan					(Quat &r, Quat const &x);

	void  r_r_atan					(Quat &r, Quat const &x);
	void  c_c_atan					(Quat &r, Quat const &x);
	void  q_q_atan					(Quat &r, Quat const &x);
	void r_rr_atan					(Quat &r, Quat const &x, Quat const &y);
	void c_rc_atan					(Quat &r, Quat const &x, Quat const &y);
	void q_rq_atan					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_atan					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_atan					(Quat &r, Quat const &x, Quat const &y);
	void q_cq_atan					(Quat &r, Quat const &x, Quat const &y);
	void q_qr_atan					(Quat &r, Quat const &x, Quat const &y);
	void q_qc_atan					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_atan					(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_tanh					(Quat &r, Quat const &x);
	void  c_c_tanh					(Quat &r, Quat const &x);
	void  q_q_tanh					(Quat &r, Quat const &x);

	void  c_c_atanh					(Quat &r, Quat const &x);
	void  q_q_atanh					(Quat &r, Quat const &x);

	void  r_r_tanc					(Quat &r, Quat const &x);
	void  c_c_tanc					(Quat &r, Quat const &x);
	void  q_q_tanc					(Quat &r, Quat const &x);

	void  r_r_cot					(Quat &r, Quat const &x);
	void  c_c_cot					(Quat &r, Quat const &x);
	void  q_q_cot					(Quat &r, Quat const &x);

	void  r_r_acot					(Quat &r, Quat const &x);
	void  c_c_acot					(Quat &r, Quat const &x);
	void  q_q_acot					(Quat &r, Quat const &x);

	void  r_r_coth					(Quat &r, Quat const &x);
	void  c_c_coth					(Quat &r, Quat const &x);
	void  q_q_coth					(Quat &r, Quat const &x);

	void  c_c_acoth					(Quat &r, Quat const &x);
	void  q_q_acoth					(Quat &r, Quat const &x);

	void  r_r_exp					(Quat &r, Quat const &x);
	void  c_c_exp					(Quat &r, Quat const &x);
	void  q_q_exp					(Quat &r, Quat const &x);
	
	void  r_r_fib					(Quat &r, Quat const &x);
	void  c_c_fib					(Quat &r, Quat const &x);
	void  q_q_fib					(Quat &r, Quat const &x);
	
	void  r_r_random				(Quat &r, Quat const &x);
	void  c_c_random				(Quat &r, Quat const &x);
	void  q_q_random				(Quat &r, Quat const &x);
	void r_rr_random				(Quat &r, Quat const &x, Quat const &y);
	void c_cr_random				(Quat &r, Quat const &x, Quat const &y);
	void c_cc_random				(Quat &r, Quat const &x, Quat const &y);
	void q_qq_random				(Quat &r, Quat const &x, Quat const &y);

	void  r_r_beta					(Quat &r, Quat const &x);
	void r_rr_beta					(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_cyl_bessel_j			(Quat &r, Quat const &x);
	void r_rr_cyl_bessel_j			(Quat &r, Quat const &x, Quat const &y);

	void  r_r_cyl_neumann			(Quat &r, Quat const &x);
	void r_rr_cyl_neumann			(Quat &r, Quat const &x, Quat const &y);

	void  c_r_hankel1				(Quat &r, Quat const &x);
	void  c_c_hankel1				(Quat &r, Quat const &x);
	void c_rr_hankel1				(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_step					(Quat &r, Quat const &x);
	void  c_c_step					(Quat &r, Quat const &x);
	void  q_q_step					(Quat &r, Quat const &x);

	void  r_r_rect					(Quat &r, Quat const &x);
	void  c_c_rect					(Quat &r, Quat const &x);
	void  q_q_rect					(Quat &r, Quat const &x);

	void  r_r_trgl					(Quat &r, Quat const &x);
	void  r_c_trgl					(Quat &r, Quat const &x);
	void  r_q_trgl					(Quat &r, Quat const &x);

	void  r_r_sqwv					(Quat &r, Quat const &x);
	void  r_c_sqwv					(Quat &r, Quat const &x);
	void  r_q_sqwv					(Quat &r, Quat const &x);
	void r_rr_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_rc_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_rq_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_cr_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_cc_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_cq_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_qr_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_qc_sqwv					(Quat &r, Quat const &x, Quat const &y);
	void r_qq_sqwv					(Quat &r, Quat const &x, Quat const &y);

	void  r_r_trwv					(Quat &r, Quat const &x);
	void  r_c_trwv					(Quat &r, Quat const &x);
	void  r_q_trwv					(Quat &r, Quat const &x);
	void r_rr_trwv					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_trwv					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_trwv					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_trwv					(Quat &r, Quat const &x, Quat const &y);

	void  r_r_saw					(Quat &r, Quat const &x);
	void  c_c_saw					(Quat &r, Quat const &x);
	void  q_q_saw					(Quat &r, Quat const &x);
	void r_rr_saw					(Quat &r, Quat const &x, Quat const &y);
	void c_rc_saw					(Quat &r, Quat const &x, Quat const &y);
	void q_rq_saw					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_saw					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_saw					(Quat &r, Quat const &x, Quat const &y);
	void q_cq_saw					(Quat &r, Quat const &x, Quat const &y);
	void q_qr_saw					(Quat &r, Quat const &x, Quat const &y);
	void q_qc_saw					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_saw					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_hypot					(Quat &r, Quat const &x, Quat const &y);
//	void c_cc_hypot					(Quat &r, Quat const &x, Quat const &y);
//	void q_qq_hypot					(Quat &r, Quat const &x, Quat const &y);

	void r_r_mandelbrot				(Quat &r, Quat const &x);
	void r_c_mandelbrot				(Quat &r, Quat const &x);
	void r_rr_mandelbrot			(Quat &r, Quat const &x, Quat const &y);
	void r_cr_mandelbrot			(Quat &r, Quat const &x, Quat const &y);

	void r_rr_min					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_min					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_min					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_min					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_max					(Quat &r, Quat const &x, Quat const &y);
	void c_cr_max					(Quat &r, Quat const &x, Quat const &y);
	void c_cc_max					(Quat &r, Quat const &x, Quat const &y);
	void q_qq_max					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void c_rc_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void q_rq_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void r_cr_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void c_cc_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void q_cq_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void r_qr_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void c_qc_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	void q_qq_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	
	void r_rr_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void c_rc_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void q_rq_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void r_cr_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void c_cc_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void q_cq_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void r_qr_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void c_qc_conditional_101		(Quat &r, Quat const &x, Quat const &y);
	void q_qq_conditional_101		(Quat &r, Quat const &x, Quat const &y);

	void  r_r_increment				(Quat &r, Quat const &x);
	void  c_c_increment				(Quat &r, Quat const &x);
	void  q_q_increment				(Quat &r, Quat const &x);

	void  r_r_decrement				(Quat &r, Quat const &x);
	void  c_c_decrement				(Quat &r, Quat const &x);
	void  q_q_decrement				(Quat &r, Quat const &x);

	void  r_r_assign				(Quat &r, Quat const &x);
	void  c_c_assign				(Quat &r, Quat const &x);
	void  q_q_assign				(Quat &r, Quat const &x);
#else
	void  r_r_setzero				(Real &r, Real const&);
	void  c_c_setzero				(Comp &r, Comp const&);
	void  q_q_setzero				(Quat &r, Quat const&);
	
	void  r_r_ceil					(Real &r, Real const &x);
	void  c_c_ceil					(Comp &r, Comp const &x);
	void  q_q_ceil					(Quat &r, Quat const &x);

	void  r_r_floor					(Real &r, Real const &x);
	void  c_c_floor					(Comp &r, Comp const &x);
	void  q_q_floor					(Quat &r, Quat const &x);

	void  r_r_round					(Real &r, Real const &x);
	void  c_c_round					(Comp &r, Comp const &x);
	void  q_q_round					(Quat &r, Quat const &x);

	void  r_r_int					(Real &r, Real const &x);
	void  c_c_int					(Comp &r, Comp const &x);
	void  q_q_int					(Quat &r, Quat const &x);

	void  r_r_frac					(Real &r, Real const &x);
	void  c_c_frac					(Comp &r, Comp const &x);
	void  q_q_frac					(Quat &r, Quat const &x);

	void  r_r_abs					(Real &r, Real const &x);
	void  r_c_abs					(Real &r, Comp const &x);
	void  r_q_abs					(Real &r, Quat const &x);

	void  r_r_arg					(Real &r, Real const &x);
	void  r_c_arg					(Real &r, Comp const &x);
	void  r_q_arg					(Real &r, Quat const &x);

	void  r_c_real					(Real &r, Comp const &x);

	void  r_c_imag					(Real &r, Comp const &x);

	//r_conjugate: assign
	void c_c_conjugate				(Comp &r, Comp const &x);
	void q_q_conjugate				(Quat &r, Quat const &x);

	void  c_r_polar					(Comp &r, Real const &x);
	void  c_c_polar					(Comp &r, Comp const &x);
	void  c_q_polar					(Comp &r, Quat const &x);

	//r_cartesian	assign
	void  c_c_cartesian				(Comp &r, Comp const &x);
	void  q_q_cartesian				(Quat &r, Quat const &x);

	void r_rr_plus					(Real &r, Real const &x, Real const &y);
	void c_rc_plus					(Comp &r, Real const &x, Comp const &y);
	void q_rq_plus					(Quat &r, Real const &x, Quat const &y);
	void c_cr_plus					(Comp &r, Comp const &x, Real const &y);
	void c_cc_plus					(Comp &r, Comp const &x, Comp const &y);
	void q_cq_plus					(Quat &r, Comp const &x, Quat const &y);
	void q_qr_plus					(Quat &r, Quat const &x, Real const &y);
	void q_qc_plus					(Quat &r, Quat const &x, Comp const &y);
	void q_qq_plus					(Quat &r, Quat const &x, Quat const &y);

	void  r_r_minus					(Real &r, Real const &x);
	void  c_c_minus					(Comp &r, Comp const &x);
	void  q_q_minus					(Quat &r, Quat const &x);
	void r_rr_minus					(Real &r, Real const &x, Real const &y);
	void c_rc_minus					(Comp &r, Real const &x, Comp const &y);
	void q_rq_minus					(Quat &r, Real const &x, Quat const &y);
	void c_cr_minus					(Comp &r, Comp const &x, Real const &y);
	void c_cc_minus					(Comp &r, Comp const &x, Comp const &y);
	void q_cq_minus					(Quat &r, Comp const &x, Quat const &y);
	void q_qr_minus					(Quat &r, Quat const &x, Real const &y);
	void q_qc_minus					(Quat &r, Quat const &x, Comp const &y);
	void q_qq_minus					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_multiply				(Real &r, Real const &x, Real const &y);
	void c_rc_multiply				(Comp &r, Real const &x, Comp const &y);
	void q_rq_multiply				(Quat &r, Real const &x, Quat const &y);
	void c_cr_multiply				(Comp &r, Comp const &x, Real const &y);
	void c_cc_multiply				(Comp &r, Comp const &x, Comp const &y);
	void q_cq_multiply				(Quat &r, Comp const &x, Quat const &y);
	void q_qr_multiply				(Quat &r, Quat const &x, Real const &y);
	void q_qc_multiply				(Quat &r, Quat const &x, Comp const &y);
	void q_qq_multiply				(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_divide				(Real &r, Real const &y);
	void  c_c_divide				(Comp &r, Comp const &y);
	void  q_q_divide				(Quat &r, Quat const &y);
	void r_rr_divide				(Real &r, Real const &x, Real const &y);
	void c_rc_divide				(Comp &r, Real const &x, Comp const &y);
	void q_rq_divide				(Quat &r, Real const &x, Quat const &y);
	void c_cr_divide				(Comp &r, Comp const &x, Real const &y);
	void c_cc_divide				(Comp &r, Comp const &x, Comp const &y);
	void q_cq_divide				(Quat &r, Comp const &x, Quat const &y);
	void q_qr_divide				(Quat &r, Quat const &x, Real const &y);
	void q_qc_divide				(Quat &r, Quat const &x, Comp const &y);
	void q_qq_divide				(Quat &r, Quat const &x, Quat const &y);
	
	void r_rr_logic_divides			(Real &r, Real const &y, Real const &x);//rc_divides, rq_divides: applied to each component
	void r_rc_logic_divides			(Real &r, Real const &y, Comp const &x);
	void r_rq_logic_divides			(Real &r, Real const &y, Quat const &x);
	void r_cr_logic_divides			(Real &r, Comp const &y, Real const &x);
	void r_cc_logic_divides			(Real &r, Comp const &y, Comp const &x);
	void r_cq_logic_divides			(Real &r, Comp const &y, Quat const &x);
	void r_qr_logic_divides			(Real &r, Quat const &y, Real const &x);
	void r_qc_logic_divides			(Real &r, Quat const &y, Comp const &x);
	void r_qq_logic_divides			(Real &r, Quat const &y, Quat const &x);

	void r_rr_power_real			(Real &r, Real const &x, Real const &y);//trunc
	void c_cr_power_real			(Comp &r, Comp const &x, Real const &y);
	void q_qr_power_real			(Quat &r, Quat const &x, Real const &y);

	void c_cr_pow					(Comp &r, Comp const &x, Real const &y);
	void c_cc_pow					(Comp &r, Comp const &x, Comp const &y);
	void q_cq_pow					(Quat &r, Comp const &x, Quat const &y);
	void q_qr_pow					(Quat &r, Quat const &x, Real const &y);
	void q_qc_pow					(Quat &r, Quat const &x, Comp const &y);
	void q_qq_pow					(Quat &r, Quat const &x, Quat const &y); 

	void  c_c_ln					(Comp &r, Comp const &x);
	void  q_q_ln					(Quat &r, Quat const &x);
	
	void  c_c_log					(Comp &r, Comp const &x);
	void  q_q_log					(Quat &r, Quat const &x);
	void c_cr_log					(Comp &r, Comp const &x, Real const &y);
	void c_cc_log					(Comp &r, Comp const &x, Comp const &y);
	void q_cq_log					(Quat &r, Comp const &x, Quat const &y);
	void q_qc_log					(Quat &r, Quat const &x, Comp const &y);
	void q_qq_log					(Quat &r, Quat const &x, Quat const &y);
	
	void c_rr_tetrate				(Comp &r, Real const &x, Real const &y);
	void c_rc_tetrate				(Comp &r, Real const &x, Comp const &y);
	void c_cr_tetrate				(Comp &r, Comp const &x, Real const &y);
	void c_cc_tetrate				(Comp &r, Comp const &x, Comp const &y);
	void q_qr_tetrate				(Quat &r, Quat const &x, Real const &y);
	
	void c_rr_pentate				(Comp &r, Real const &x, Real const &y);
	void c_cr_pentate				(Comp &r, Comp const &x, Real const &y);

	void  r_r_bitwise_shift_left_l	(Real &r, Real const &x);	//<<x = 2^floor(x)
	void  c_c_bitwise_shift_left_l	(Comp &r, Comp const &x);
	void  q_q_bitwise_shift_left_l	(Quat &r, Quat const &x);
	void  r_r_bitwise_shift_left_r	(Real &r, Real const &x);	//x<< = 2x
	void  c_c_bitwise_shift_left_r	(Comp &r, Comp const &x);
	void  q_q_bitwise_shift_left_r	(Quat &r, Quat const &x);
	void r_rr_bitwise_shift_left	(Real &r, Real const &x, Real const &y);//x<<y = x*2^y
	void c_rc_bitwise_shift_left	(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_shift_left	(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_shift_left	(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_shift_left	(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_shift_left	(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_shift_left	(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_shift_left	(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_shift_right_l	(Real &r, Real const &x);
	void  c_c_bitwise_shift_right_l	(Comp &r, Comp const &x);
	void  q_q_bitwise_shift_right_l	(Quat &r, Quat const &x);
	void  r_r_bitwise_shift_right_r	(Real &r, Real const &x);
	void  c_c_bitwise_shift_right_r	(Comp &r, Comp const &x);
	void  q_q_bitwise_shift_right_r	(Quat &r, Quat const &x);
	void r_rr_bitwise_shift_right	(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_shift_right	(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_shift_right	(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_shift_right	(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_shift_right	(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_shift_right	(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_shift_right	(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_shift_right	(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_not			(Real &r, Real const &x);
	void  c_c_bitwise_not			(Comp &r, Comp const &x);
	void  q_q_bitwise_not			(Quat &r, Quat const &x);

	void  r_r_bitwise_and			(Real &r, Real const &x);
	void  c_c_bitwise_and			(Comp &r, Comp const &x);
	void  q_q_bitwise_and			(Quat &r, Quat const &x);
	void r_rr_bitwise_and			(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_and			(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_and			(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_and			(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_and			(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_and			(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_and			(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_and			(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_and			(Quat &r, Quat const &x, Quat const &y);

	void  r_r_bitwise_nand			(Real &r, Real const &x);
	void  c_c_bitwise_nand			(Comp &r, Comp const &x);
	void  q_q_bitwise_nand			(Quat &r, Quat const &x);
	void r_rr_bitwise_nand			(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_nand			(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_nand			(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_nand			(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_nand			(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_nand			(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_nand			(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_nand			(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_or			(Real &r, Real const &x);
	void  c_c_bitwise_or			(Comp &r, Comp const &x);
	void  q_q_bitwise_or			(Quat &r, Quat const &x);
	void r_rr_bitwise_or			(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_or			(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_or			(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_or			(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_or			(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_or			(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_or			(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_or			(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_or			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_nor			(Real &r, Real const &x);
	void  c_c_bitwise_nor			(Comp &r, Comp const &x);
	void  q_q_bitwise_nor			(Quat &r, Quat const &x);
	void r_rr_bitwise_nor			(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_nor			(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_nor			(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_nor			(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_nor			(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_nor			(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_nor			(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_nor			(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_xor			(Real &r, Real const &x);
	void  c_c_bitwise_xor			(Comp &r, Comp const &x);
	void  q_q_bitwise_xor			(Quat &r, Quat const &x);
	void r_rr_bitwise_xor			(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_xor			(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_xor			(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_xor			(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_xor			(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_xor			(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_xor			(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_xor			(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_bitwise_xnor			(Real &r, Real const &x);
	void  c_c_bitwise_xnor			(Comp &r, Comp const &x);
	void  q_q_bitwise_xnor			(Quat &r, Quat const &x);
	void r_rr_bitwise_xnor			(Real &r, Real const &x, Real const &y);
	void c_rc_bitwise_xnor			(Comp &r, Real const &x, Comp const &y);
	void q_rq_bitwise_xnor			(Quat &r, Real const &x, Quat const &y);
	void c_cr_bitwise_xnor			(Comp &r, Comp const &x, Real const &y);
	void c_cc_bitwise_xnor			(Comp &r, Comp const &x, Comp const &y);
	void q_cq_bitwise_xnor			(Quat &r, Comp const &x, Quat const &y);
	void q_qr_bitwise_xnor			(Quat &r, Quat const &x, Real const &y);
	void q_qc_bitwise_xnor			(Quat &r, Quat const &x, Comp const &y);
	void q_qq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_equal			(Real &r, Real const &x);
	void  r_c_logic_equal			(Real &r, Comp const &x);
	void  r_q_logic_equal			(Real &r, Quat const &x);
	void r_rr_logic_equal			(Real &r, Real const &x, Real const &y);
	void r_rc_logic_equal			(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_equal			(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_equal			(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_equal			(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_equal			(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_equal			(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_equal			(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_equal			(Real &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_not_equal		(Real &r, Real const &x);
	void  r_c_logic_not_equal		(Real &r, Comp const &x);
	void  r_q_logic_not_equal		(Real &r, Quat const &x);
	void r_rr_logic_not_equal		(Real &r, Real const &x, Real const &y);
	void r_rc_logic_not_equal		(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_not_equal		(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_not_equal		(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_not_equal		(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_not_equal		(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_not_equal		(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_not_equal		(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_not_equal		(Real &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_less_l			(Real &r, Real const &x);
	void  r_c_logic_less_l			(Real &r, Comp const &x);
	void  r_q_logic_less_l			(Real &r, Quat const &x);
	void  r_r_logic_less_r			(Real &r, Real const &x);
	void  r_c_logic_less_r			(Real &r, Comp const &x);
	void  r_q_logic_less_r			(Real &r, Quat const &x);
	void r_rr_logic_less			(Real &r, Real const &x, Real const &y);
	void r_rc_logic_less			(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_less			(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_less			(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_less			(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_less			(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_less			(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_less			(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_less			(Real &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_less_equal_l	(Real &r, Real const &x);
	void  r_c_logic_less_equal_l	(Real &r, Comp const &x);
	void  r_q_logic_less_equal_l	(Real &r, Quat const &x);
	void  r_r_logic_less_equal_r	(Real &r, Real const &x);
	void  r_c_logic_less_equal_r	(Real &r, Comp const &x);
	void  r_q_logic_less_equal_r	(Real &r, Quat const &x);
	void r_rr_logic_less_equal		(Real &r, Real const &x, Real const &y);
	void r_rc_logic_less_equal		(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_less_equal		(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_less_equal		(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_less_equal		(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_less_equal		(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_less_equal		(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_less_equal		(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_less_equal		(Real &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_greater_l		(Real &r, Real const &x);
	void  r_c_logic_greater_l		(Real &r, Comp const &x);
	void  r_q_logic_greater_l		(Real &r, Quat const &x);
	void  r_r_logic_greater_r		(Real &r, Real const &x);
	void  r_c_logic_greater_r		(Real &r, Comp const &x);
	void  r_q_logic_greater_r		(Real &r, Quat const &x);
	void r_rr_logic_greater			(Real &r, Real const &x, Real const &y);
	void r_rc_logic_greater			(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_greater			(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_greater			(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_greater			(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_greater			(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_greater			(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_greater			(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_greater			(Real &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_greater_equal_l	(Real &r, Real const &x);
	void  r_c_logic_greater_equal_l	(Real &r, Comp const &x);
	void  r_q_logic_greater_equal_l	(Real &r, Quat const &x);
	void  r_r_logic_greater_equal_r	(Real &r, Real const &x);
	void  r_c_logic_greater_equal_r	(Real &r, Comp const &x);
	void  r_q_logic_greater_equal_r	(Real &r, Quat const &x);
	void r_rr_logic_greater_equal	(Real &r, Real const &x, Real const &y);
	void r_rc_logic_greater_equal	(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_greater_equal	(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_greater_equal	(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_greater_equal	(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_greater_equal	(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_greater_equal	(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_greater_equal	(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_greater_equal	(Real &r, Quat const &x, Quat const &y);
	
	void  r_r_logic_not				(Real &r, Real const &x);
	void  r_c_logic_not				(Real &r, Comp const &x);
	void  r_q_logic_not				(Real &r, Quat const &x);
	
	void r_rr_logic_and				(Real &r, Real const &x, Real const &y);
	void r_rc_logic_and				(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_and				(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_and				(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_and				(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_and				(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_and				(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_and				(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_and				(Real &r, Quat const &x, Quat const &y);
	
	void r_rr_logic_or				(Real &r, Real const &x, Real const &y);
	void r_rc_logic_or				(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_or				(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_or				(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_or				(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_or				(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_or				(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_or				(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_or				(Real &r, Quat const &x, Quat const &y);

	void r_rr_logic_xor				(Real &r, Real const &x, Real const &y);
	void r_rc_logic_xor				(Real &r, Real const &x, Comp const &y);
	void r_rq_logic_xor				(Real &r, Real const &x, Quat const &y);
	void r_cr_logic_xor				(Real &r, Comp const &x, Real const &y);
	void r_cc_logic_xor				(Real &r, Comp const &x, Comp const &y);
	void r_cq_logic_xor				(Real &r, Comp const &x, Quat const &y);
	void r_qr_logic_xor				(Real &r, Quat const &x, Real const &y);
	void r_qc_logic_xor				(Real &r, Quat const &x, Comp const &y);
	void r_qq_logic_xor				(Real &r, Quat const &x, Quat const &y);

	void r_rr_condition_zero		(Real &r, Real const &x, Real const &y);
	void c_rc_condition_zero		(Comp &r, Real const &x, Comp const &y);
	void q_rq_condition_zero		(Quat &r, Real const &x, Quat const &y);
	void c_cr_condition_zero		(Comp &r, Comp const &x, Real const &y);
	void c_cc_condition_zero		(Comp &r, Comp const &x, Comp const &y);
	void q_cq_condition_zero		(Quat &r, Comp const &x, Quat const &y);
	void q_qr_condition_zero		(Quat &r, Quat const &x, Real const &y);
	void q_qc_condition_zero		(Quat &r, Quat const &x, Comp const &y);
	void q_qq_condition_zero		(Quat &r, Quat const &x, Quat const &y);

	void  r_r_percent				(Real &r, Real const &x);
	void  c_c_percent				(Comp &r, Comp const &x);
	void  q_q_percent				(Quat &r, Quat const &x);
	
	void r_rr_modulo				(Real &r, Real const &x, Real const &y);
	void c_rc_modulo				(Comp &r, Real const &x, Comp const &y);
	void q_rq_modulo				(Quat &r, Real const &x, Quat const &y);
	void c_cr_modulo				(Comp &r, Comp const &x, Real const &y);
	void c_cc_modulo				(Comp &r, Comp const &x, Comp const &y);
	void q_cq_modulo				(Quat &r, Comp const &x, Quat const &y);
	void q_qr_modulo				(Quat &r, Quat const &x, Real const &y);
	void q_qc_modulo				(Quat &r, Quat const &x, Comp const &y);
	void q_qq_modulo				(Quat &r, Quat const &x, Quat const &y);

	void  r_r_sgn					(Real &r, Real const &x);
	void  c_c_sgn					(Comp &r, Comp const &x);
	void  q_q_sgn					(Quat &r, Quat const &x);
	
	void  r_r_sq					(Real &r, Real const &x);
	void  c_c_sq					(Comp &r, Comp const &x);
	void  q_q_sq					(Quat &r, Quat const &x);
	
	void  c_c_sqrt					(Comp &r, Comp const &x);
	void  q_q_sqrt					(Quat &r, Quat const &x);

	void  r_r_invsqrt				(Real &r, Real const &x);

	void  r_r_cbrt					(Real &r, Real const &x);
	void  c_c_cbrt					(Comp &r, Comp const &x);
	void  q_q_cbrt					(Quat &r, Quat const &x);

	void  r_r_gauss					(Real &r, Real const &x);
	void  c_c_gauss					(Comp &r, Comp const &x);
	void  q_q_gauss					(Quat &r, Quat const &x);

	void  r_r_erf					(Real &r, Real const &x);

	void  r_r_zeta					(Real &r, Real const &x);
	
	void  r_r_tgamma				(Real &r, Real const &x);
	void  c_c_tgamma				(Comp &r, Comp const &x);
	void  q_q_tgamma				(Quat &r, Quat const &x);
	void r_rr_tgamma				(Real &r, Real const &x, Real const &y);

	void  r_r_loggamma				(Real &r, Real const &x);

	void  r_r_factorial				(Real &r, Real const &x);
	void  c_c_factorial				(Comp &r, Comp const &x);
	void  q_q_factorial				(Quat &r, Quat const &x);

	void  r_r_permutation			(Real &r, Real const &x);
	void  c_c_permutation			(Comp &r, Comp const &x);
	void  q_q_permutation			(Quat &r, Quat const &x);
	void r_rr_permutation			(Real &r, Real const &x, Real const &y);
	void c_cr_permutation			(Comp &r, Comp const &x, Real const &y);
	void c_cc_permutation			(Comp &r, Comp const &x, Comp const &y);
	void q_qq_permutation			(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_combination			(Real &r, Real const &x);
	void  c_c_combination			(Comp &r, Comp const &x);
	void  q_q_combination			(Quat &r, Quat const &x);
	void r_rr_combination			(Real &r, Real const &x, Real const &y);
	void c_cr_combination			(Comp &r, Comp const &x, Real const &y);
	void c_cc_combination			(Comp &r, Comp const &x, Comp const &y);
	void q_qq_combination			(Quat &r, Quat const &x, Quat const &y);

	void  r_r_cos					(Real &r, Real const &x);
	void  c_c_cos					(Comp &r, Comp const &x);
	void  q_q_cos					(Quat &r, Quat const &x);

	void  c_c_acos					(Comp &r, Comp const &x);
	void  q_q_acos					(Quat &r, Quat const &x);

	void  r_r_cosh					(Real &r, Real const &x);
	void  c_c_cosh					(Comp &r, Comp const &x);
	void  q_q_cosh					(Quat &r, Quat const &x);
	
	void  c_c_acosh					(Comp &r, Comp const &x);
	void  q_q_acosh					(Quat &r, Quat const &x);

	void  r_r_cosc					(Real &r, Real const &x);
	void  c_c_cosc					(Comp &r, Comp const &x);
	void  q_q_cosc					(Quat &r, Quat const &x);

	void  r_r_sec					(Real &r, Real const &x);
	void  c_c_sec					(Comp &r, Comp const &x);
	void  q_q_sec					(Quat &r, Quat const &x);

	void  c_c_asec					(Comp &r, Comp const &x);
	void  q_q_asec					(Quat &r, Quat const &x);

	void  r_r_sech					(Real &r, Real const &x);
	void  c_c_sech					(Comp &r, Comp const &x);
	void  q_q_sech					(Quat &r, Quat const &x);

	void  c_c_asech					(Comp &r, Comp const &x);
	void  q_q_asech					(Quat &r, Quat const &x);

	void  r_r_sin					(Real &r, Real const &x);
	void  c_c_sin					(Comp &r, Comp const &x);
	void  q_q_sin					(Quat &r, Quat const &x);

	void  c_c_asin					(Comp &r, Comp const &x);
	void  q_q_asin					(Quat &r, Quat const &x);

	void  r_r_sinh					(Real &r, Real const &x);
	void  c_c_sinh					(Comp &r, Comp const &x);
	void  q_q_sinh					(Quat &r, Quat const &x);

	void  r_r_asinh					(Real &r, Real const &x);
	void  c_c_asinh					(Comp &r, Comp const &x);
	void  q_q_asinh					(Quat &r, Quat const &x);

	void  r_r_sinc					(Real &r, Real const &x);
	void  c_c_sinc					(Comp &r, Comp const &x);
	void  q_q_sinc					(Quat &r, Quat const &x);

	void  r_r_sinhc					(Real &r, Real const &x);
	void  c_c_sinhc					(Comp &r, Comp const &x);
	void  q_q_sinhc					(Quat &r, Quat const &x);

	void  r_r_csc					(Real &r, Real const &x);
	void  c_c_csc					(Comp &r, Comp const &x);
	void  q_q_csc					(Quat &r, Quat const &x);

	void  c_c_acsc					(Comp &r, Comp const &x);
	void  q_q_acsc					(Quat &r, Quat const &x);

	void  r_r_csch					(Real &r, Real const &x);
	void  c_c_csch					(Comp &r, Comp const &x);
	void  q_q_csch					(Quat &r, Quat const &x);

	void  r_r_acsch					(Real &r, Real const &x);
	void  c_c_acsch					(Comp &r, Comp const &x);
	void  q_q_acsch					(Quat &r, Quat const &x);

	void  r_r_tan					(Real &r, Real const &x);
	void  c_c_tan					(Comp &r, Comp const &x);
	void  q_q_tan					(Quat &r, Quat const &x);

	void  r_r_atan					(Real &r, Real const &x);
	void  c_c_atan					(Comp &r, Comp const &x);
	void  q_q_atan					(Quat &r, Quat const &x);
	void r_rr_atan					(Real &r, Real const &x, Real const &y);
	void c_rc_atan					(Comp &r, Real const &x, Comp const &y);
	void q_rq_atan					(Quat &r, Real const &x, Quat const &y);
	void c_cr_atan					(Comp &r, Comp const &x, Real const &y);
	void c_cc_atan					(Comp &r, Comp const &x, Comp const &y);
	void q_cq_atan					(Quat &r, Comp const &x, Quat const &y);
	void q_qr_atan					(Quat &r, Quat const &x, Real const &y);
	void q_qc_atan					(Quat &r, Quat const &x, Comp const &y);
	void q_qq_atan					(Quat &r, Quat const &x, Quat const &y);
	
	void  r_r_tanh					(Real &r, Real const &x);
	void  c_c_tanh					(Comp &r, Comp const &x);
	void  q_q_tanh					(Quat &r, Quat const &x);

	void  c_c_atanh					(Comp &r, Comp const &x);
	void  q_q_atanh					(Quat &r, Quat const &x);

	void  r_r_tanc					(Real &r, Real const &x);
	void  c_c_tanc					(Comp &r, Comp const &x);
	void  q_q_tanc					(Quat &r, Quat const &x);

	void  r_r_cot					(Real &r, Real const &x);
	void  c_c_cot					(Comp &r, Comp const &x);
	void  q_q_cot					(Quat &r, Quat const &x);

	void  r_r_acot					(Real &r, Real const &x);
	void  c_c_acot					(Comp &r, Comp const &x);
	void  q_q_acot					(Quat &r, Quat const &x);

	void  r_r_coth					(Real &r, Real const &x);
	void  c_c_coth					(Comp &r, Comp const &x);
	void  q_q_coth					(Quat &r, Quat const &x);

	void  c_c_acoth					(Comp &r, Comp const &x);
	void  q_q_acoth					(Quat &r, Quat const &x);

	void  r_r_exp					(Real &r, Real const &x);
	void  c_c_exp					(Comp &r, Comp const &x);
	void  q_q_exp					(Quat &r, Quat const &x);
	
	void  r_r_fib					(Real &r, Real const &x);
	void  c_c_fib					(Comp &r, Comp const &x);
	void  q_q_fib					(Quat &r, Quat const &x);
	
	void  r_r_random				(Real &r, Real const &x);
	void  c_c_random				(Comp &r, Comp const &x);
	void  q_q_random				(Quat &r, Quat const &x);
	void r_rr_random				(Real &r, Real const &x, Real const &y);
	void c_cr_random				(Comp &r, Comp const &x, Real const &y);
	void c_cc_random				(Comp &r, Comp const &x, Comp const &y);
	void q_qq_random				(Quat &r, Quat const &x, Quat const &y);

	void  r_r_beta					(Real &r, Real const &x);
	void r_rr_beta					(Real &r, Real const &x, Real const &y);
	
	void  r_r_cyl_bessel_j			(Real &r, Real const &x);
	void r_rr_cyl_bessel_j			(Real &r, Real const &x, Real const &y);

	void  r_r_cyl_neumann			(Real &r, Real const &x);
	void r_rr_cyl_neumann			(Real &r, Real const &x, Real const &y);

	void  c_r_hankel1				(Comp &r, Real const &x);
	void  c_c_hankel1				(Comp &r, Comp const &x);
	void c_rr_hankel1				(Comp &r, Real const &x, Real const &y);
	
	void  r_r_step					(Real &r, Real const &x);
	void  c_c_step					(Comp &r, Comp const &x);
	void  q_q_step					(Quat &r, Quat const &x);

	void  r_r_rect					(Real &r, Real const &x);
	void  c_c_rect					(Comp &r, Comp const &x);
	void  q_q_rect					(Quat &r, Quat const &x);

	void  r_r_trgl					(Real &r, Real const &x);
	void  r_c_trgl					(Real &r, Comp const &x);
	void  r_q_trgl					(Real &r, Quat const &x);

	void  r_r_sqwv					(Real &r, Real const &x);
	void  r_c_sqwv					(Real &r, Comp const &x);
	void  r_q_sqwv					(Real &r, Quat const &x);
	void r_rr_sqwv					(Real &r, Real const &x, Real const &y);
	void r_rc_sqwv					(Real &r, Real const &x, Comp const &y);
	void r_rq_sqwv					(Real &r, Real const &x, Quat const &y);
	void r_cr_sqwv					(Real &r, Comp const &x, Real const &y);
	void r_cc_sqwv					(Real &r, Comp const &x, Comp const &y);
	void r_cq_sqwv					(Real &r, Comp const &x, Quat const &y);
	void r_qr_sqwv					(Real &r, Quat const &x, Real const &y);
	void r_qc_sqwv					(Real &r, Quat const &x, Comp const &y);
	void r_qq_sqwv					(Real &r, Quat const &x, Quat const &y);

	void  r_r_trwv					(Real &r, Real const &x);
	void  r_c_trwv					(Real &r, Comp const &x);
	void  r_q_trwv					(Real &r, Quat const &x);
	void r_rr_trwv					(Real &r, Real const &x, Real const &y);
	void c_cr_trwv					(Comp &r, Comp const &x, Real const &y);
	void c_cc_trwv					(Comp &r, Comp const &x, Comp const &y);
	void q_qq_trwv					(Quat &r, Quat const &x, Quat const &y);

	void  r_r_saw					(Real &r, Real const &x);
	void  c_c_saw					(Comp &r, Comp const &x);
	void  q_q_saw					(Quat &r, Quat const &x);
	void r_rr_saw					(Real &r, Real const &x, Real const &y);
	void c_rc_saw					(Comp &r, Real const &x, Comp const &y);
	void q_rq_saw					(Quat &r, Real const &x, Quat const &y);
	void c_cr_saw					(Comp &r, Comp const &x, Real const &y);
	void c_cc_saw					(Comp &r, Comp const &x, Comp const &y);
	void q_cq_saw					(Quat &r, Comp const &x, Quat const &y);
	void q_qr_saw					(Quat &r, Quat const &x, Real const &y);
	void q_qc_saw					(Quat &r, Quat const &x, Comp const &y);
	void q_qq_saw					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_hypot					(Real &r, Real const &x, Real const &y);
//	void c_cc_hypot					(Real &r, Comp const &x, Comp const &y);
//	void q_qq_hypot					(Real &r, Quat const &x, Quat const &y);

	void r_r_mandelbrot				(Real &r, Real const &x);
	void r_c_mandelbrot				(Real &r, Comp const &x);
	void r_rr_mandelbrot			(Real &r, Real const &x, Real const &y);
	void r_cr_mandelbrot			(Real &r, Comp const &x, Real const &y);

	void r_rr_min					(Real &r, Real const &x, Real const &y);
	void c_cr_min					(Comp &r, Comp const &x, Real const &y);
	void c_cc_min					(Comp &r, Comp const &x, Comp const &y);
	void q_qq_min					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_max					(Real &r, Real const &x, Real const &y);
	void c_cr_max					(Comp &r, Comp const &x, Real const &y);
	void c_cc_max					(Comp &r, Comp const &x, Comp const &y);
	void q_qq_max					(Quat &r, Quat const &x, Quat const &y);

	void r_rr_conditional_110		(Real &r, Real const &x, Real const &y);
	void c_rc_conditional_110		(Comp &r, Real const &x, Comp const &y);
	void q_rq_conditional_110		(Quat &r, Real const &x, Quat const &y);
	void r_cr_conditional_110		(Real &r, Comp const &x, Real const &y);
	void c_cc_conditional_110		(Comp &r, Comp const &x, Comp const &y);
	void q_cq_conditional_110		(Quat &r, Comp const &x, Quat const &y);
	void r_qr_conditional_110		(Real &r, Quat const &x, Real const &y);
	void c_qc_conditional_110		(Comp &r, Quat const &x, Comp const &y);
	void q_qq_conditional_110		(Quat &r, Quat const &x, Quat const &y);
	
	void r_rr_conditional_101		(Real &r, Real const &x, Real const &y);
	void c_rc_conditional_101		(Comp &r, Real const &x, Comp const &y);
	void q_rq_conditional_101		(Quat &r, Real const &x, Quat const &y);
	void r_cr_conditional_101		(Real &r, Comp const &x, Real const &y);
	void c_cc_conditional_101		(Comp &r, Comp const &x, Comp const &y);
	void q_cq_conditional_101		(Quat &r, Comp const &x, Quat const &y);
	void r_qr_conditional_101		(Real &r, Quat const &x, Real const &y);
	void c_qc_conditional_101		(Comp &r, Quat const &x, Comp const &y);
	void q_qq_conditional_101		(Quat &r, Quat const &x, Quat const &y);

	void  r_r_increment				(Real &r, Real const &x);
	void  c_c_increment				(Comp &r, Comp const &x);
	void  q_q_increment				(Quat &r, Quat const &x);

	void  r_r_decrement				(Real &r, Real const &x);
	void  c_c_decrement				(Comp &r, Comp const &x);
	void  q_q_decrement				(Quat &r, Quat const &x);

	void  r_r_assign				(Real &r, Real const &x);
	void  c_c_assign				(Comp &r, Comp const &x);
	void  q_q_assign				(Quat &r, Quat const &x);
#endif
}
#endif