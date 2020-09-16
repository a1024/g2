//g2_mp.cpp - Multiprecision version of math functions.
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

#include		"g2_mp.h"
#include		"g2_common.h"
#pragma			comment(lib, "mpir32.lib")
#pragma			comment(lib, "mpfr32.lib")
namespace		MP
{
	const double inv_ln2=1/::log(2.), ln2_ln10=::log(2.)/::log(10.), ln2_ln8=::log(2.)/::log(8.), ln2_ln16=0.25;
	int bin_prec=53, dec_prec=int(bin_prec*ln2_ln10);
	Real m_pi=::acos(-1.), m_e=::exp(1.),
		m_qnan=std::numeric_limits<double>::quiet_NaN(),//use _HUGE & G2::_qnan
		m_ln2=::log(2.), m_inv_ln10=1/::log(10.), m_sqrt_2pi=::sqrt(2*G2::_pi), m_ln_phi=::log(0.5*(1+::sqrt(5.))), m_inv_sqrt5=1/::sqrt(5.);

	void	recalculate_constants()
	{
		m_pi=acos(Real(-1));
		m_e=exp(Real(1));
	//	m_qnan=std::numeric_limits<double>::quiet_NaN();
		m_ln2=log(Real(2));
		m_inv_ln10=1/log(Real(10));
		m_sqrt_2pi=sqrt(2*m_pi);
		m_ln_phi=log(0.5*(1+sqrt(Real(5))));
		m_inv_sqrt5=1/sqrt(Real(5));
	}

	void	remove_trailing_zeros(std::string &str2)
	{
		int k=str2.size()-1;
		for(;k>=0&&str2[k]=='0';--k);//remove trailing zeros
		str2.resize(k+1);
	}
	void	to_upper(std::string &str)
	{
		for(int k=0, kEnd=str.size();k<kEnd;++k)
			if(str[k]>='a'&&str[k]<='f')
				str[k]&=0xDF;
	}
	int		print_real(char *a, Real const &x, char base)//bases 2, 8, 10, 16 supported
	{
		if(x!=x)
		{
			memcpy(a, "0/0", 4);
			return 3;
		}
		if(base==10)
		{
			//std::string str=x.toString();//
			//mpfr_sprintf(a, "%Rg", x.mpfr_srcptr());//
			//char *b=nullptr;
			//int length=mpfr_asprintf(&b, "%Rg", x.mpfr_srcptr());//1.1 -> 1.10156
			//memcpy(a, b, length+1);
			//mpfr_free_str(b);
			//return mpfr_sprintf(a, "%s", x.toString("%Rg").c_str());//1.1 -> 1.10156
		//	int prec=x.getPrecision();
		//	return mpfr_sprintf(a, "%s", x.toString().c_str());
			return mpfr_sprintf(a, "%.*Rg", MP::dec_prec, x.mpfr_srcptr());//1.1 -> 1.1015625
		}
		int length=0;
		Real abs_x=abs(x);
		bool negative=signbit(x);
		if(abs_x==_HUGE)
		{
			int k=0;
			if(negative)
				a[k++]='-';
			memcpy(a+k, "inf", 4);
			return k+3;
		}
		long exponent=0;
		char *str=mpfr_get_str(0, &exponent, base, 0, abs_x.mpfr_srcptr(), MPFR_RNDNA);
		std::string str2=str;
		mpfr_free_str(str);

		const int t_size=128;
		char temp[t_size]={0};
		int exp_mul=1;
		switch(base)
		{
		case 8:exp_mul=3;break;
		case 16:exp_mul=4;break;
		}
		if(exponent<-3)//print decimal point, remove trailing zeros & possibly print exponent
		{
			str2.insert(1, 1, '.');
			remove_trailing_zeros(str2);
			int len2=sprintf_s(temp, "p%03d", exp_mul*(exponent-1));
			str2.insert(str2.size(), temp, len2);
		}
		if(exponent>=-3&&exponent<=0)
		{
			str2.insert(0, "0.0000", 2-exponent);
			remove_trailing_zeros(str2);
		}
		else if(exponent>0&&exponent<(int)str2.size())//exponent==str2.size(): point omitted
		{
			str2.insert(exponent, 1, '.');
			remove_trailing_zeros(str2);
			if(*str2.rbegin()=='.')
				str2.pop_back();
		}
		else if(exponent>(int)str2.size())
		{
			str2.insert(1, 1, '.');
			remove_trailing_zeros(str2);
			int len2=sprintf_s(temp, "p+%03d", exp_mul*(exponent-1));
			str2.insert(str2.size(), temp, len2);
		}

		switch(base)//print base prefix
		{
		case 2:str2.insert(0, 1, 'B');break;
		case 8:str2.insert(0, 1, '0');break;
		case 16:
			to_upper(str2);
			str2.insert(0, "0x");
			break;
		}

		if(negative)//print sign
			str2.insert(0, 1, '-');

		memcpy(a, str2.c_str(), str2.size()+1);
		return str2.size();
	}
	
#ifdef MP_PURE_QUAT
	void  r_r_setzero				(Quat &r, Quat const&)					{r.r.setZero();}
	void  c_c_setzero				(Quat &r, Quat const&)					{r.r.setZero(), r.i.setZero();}
	void  q_q_setzero				(Quat &r, Quat const&)					{r.r.setZero(), r.i.setZero(), r.j.setZero(), r.k.setZero();}
	
	void  r_r_ceil					(Quat &r, Quat const &x)				{r.r=ceil(x.r);}
	void  c_c_ceil					(Quat &r, Quat const &x)				{r=ceil((Comp)x);}
	void  q_q_ceil					(Quat &r, Quat const &x)				{r=ceil(x);}

	void  r_r_floor					(Quat &r, Quat const &x)				{r.r=floor(x.r);}
	void  c_c_floor					(Quat &r, Quat const &x)				{r=floor((Comp)x);}
	void  q_q_floor					(Quat &r, Quat const &x)				{r=floor(x);}

	void  r_r_round					(Quat &r, Quat const &x)				{r.r=round(x.r);}
	void  c_c_round					(Quat &r, Quat const &x)				{r=round((Comp)x);}
	void  q_q_round					(Quat &r, Quat const &x)				{r=round(x);}

	void  r_r_int					(Quat &r, Quat const &x)				{r.r=trunc(x.r);}
	void  c_c_int					(Quat &r, Quat const &x)				{r.r=trunc(x.r), r.i=trunc(x.i);}
	void  q_q_int					(Quat &r, Quat const &x)				{r.r=trunc(x.r), r.i=trunc(x.i), r.j=trunc(x.j), r.k=trunc(x.k);}

	void  r_r_frac					(Quat &r, Quat const &x)				{r.r=frac(x.r);}
	void  c_c_frac					(Quat &r, Quat const &x)				{r.r=frac(x.r), r.i=frac(x.i);}
	void  q_q_frac					(Quat &r, Quat const &x)				{r.r=frac(x.r), r.i=frac(x.i), r.j=frac(x.j), r.k=frac(x.k);}

	void  r_r_abs					(Quat &r, Quat const &x)				{r=abs(x.r);}
	void  r_c_abs					(Quat &r, Quat const &x)				{r=abs((Comp)x);}
	void  r_q_abs					(Quat &r, Quat const &x)				{r=abs(x);}

	void  r_r_arg					(Quat &r, Quat const &x)				{r=x.r<0?m_pi:x.r==0?m_qnan:0;}
	void  r_c_arg					(Quat &r, Quat const &x)				{r=(x.r==0)&(x.i==0)?m_qnan:atan2(x.i, x.r);}
	void  r_q_arg					(Quat &r, Quat const &x)				{r=acos(x.r/abs(x));}

	void  r_c_real					(Quat &r, Quat const &x)				{r=x.r;}

	void  r_c_imag					(Quat &r, Quat const &x)				{r=x.i;}

	//r_conjugate: assign
	void c_c_conjugate				(Quat &r, Quat const &x)				{r=Comp(x.r, -x.i);}
	void q_q_conjugate				(Quat &r, Quat const &x)				{r=Quat(x.r, -x.i, -x.j, -x.k);}

	void  c_r_polar					(Quat &r, Quat const &x)				{r=Comp(abs(x.r), x.r<0?m_pi:x.r==0?m_qnan:0);}
	void  c_c_polar					(Quat &r, Quat const &x)
	{
		Real mag=abs((Comp)x);
		r=Comp(mag, mag==0?m_qnan:atan2(x.i, x.r));
	}
	void  c_q_polar					(Quat &r, Quat const &x)
	{
		Real mag=abs(x);
		r=Comp(mag, acos(x.r/mag));
	}

	//r_cartesian	assign
	void  c_c_cartesian				(Quat &r, Quat const &x)				{r=Comp(x.r*cos(x.i), x.r*sin(x.i));}
	void  q_q_cartesian				(Quat &r, Quat const &x)
	{
		Real cos_j=cos(x.j), xr_cos_k=x.r*cos(x.k);
		r=Quat(
			cos(x.i)*cos_j*xr_cos_k,
			sin(x.i)*cos_j*xr_cos_k,
			sin(x.j)*xr_cos_k,
			x.r*sin(x.k));
	}

	void r_rr_plus					(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x+(Real)y;}
	void c_rc_plus					(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x+(Comp)y;}
	void q_rq_plus					(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x+y;}
	void c_cr_plus					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x+(Real)y;}
	void c_cc_plus					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x+(Comp)y;}
	void q_cq_plus					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x+y;}
	void q_qr_plus					(Quat &r, Quat const &x, Quat const &y)	{r=x+(Real)y;}
	void q_qc_plus					(Quat &r, Quat const &x, Quat const &y)	{r=x+(Comp)y;}
	void q_qq_plus					(Quat &r, Quat const &x, Quat const &y)	{r=x+y;}

	void  r_r_minus					(Quat &r, Quat const &x)				{r=-(Real)x;}
	void  c_c_minus					(Quat &r, Quat const &x)				{r=-(Comp)x;}
	void  q_q_minus					(Quat &r, Quat const &x)				{r=-x;}
	void r_rr_minus					(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x-(Real)y;}
	void c_rc_minus					(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x-(Comp)y;}
	void q_rq_minus					(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x-y;}
	void c_cr_minus					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x-(Real)y;}
	void c_cc_minus					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x-(Comp)y;}
	void q_cq_minus					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x-y;}
	void q_qr_minus					(Quat &r, Quat const &x, Quat const &y)	{r=x-(Real)y;}
	void q_qc_minus					(Quat &r, Quat const &x, Quat const &y)	{r=x-(Comp)y;}
	void q_qq_minus					(Quat &r, Quat const &x, Quat const &y)	{r=x-y;}

	void r_rr_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*(Real)y;}		//MARKER
	void c_rc_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*(Comp)y;}
	void q_rq_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*y;}
	void c_cr_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*(Real)y;}
	void c_cc_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*(Comp)y;}
	void q_cq_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*y;}
	void q_qr_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=x*(Real)y;}
	void q_qc_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=x*(Comp)y;}
	void q_qq_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=x*y;}
	
	inline Comp inv(Comp const &x)
	{
		Real inv_mag=1/abs(x);
		return Comp(x.r*inv_mag, -x.i*inv_mag);
	}
	inline Quat inv(Quat const &x)
	{
		Real inv_mag=1/abs(x);
		return Quat(x.r*inv_mag, -x.i*inv_mag, -x.j*inv_mag, -x.k*inv_mag);
	}
	void  r_r_divide				(Quat &r, Quat const &y)				{r=1/(Real)y;}
	void  c_c_divide				(Quat &r, Quat const &y)				{r=inv((Comp)y);}
	void  q_q_divide				(Quat &r, Quat const &y)				{r=inv(y);}
	void r_rr_divide				(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x/(Real)y;}
	void c_rc_divide				(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x/(Comp)y;}
	void q_rq_divide				(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x/y;}
	void c_cr_divide				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x/(Real)y;}
	void c_cc_divide				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x/(Comp)y;}
	void q_cq_divide				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x/y;}
	void q_qr_divide				(Quat &r, Quat const &x, Quat const &y)	{r=x/(Real)y;}
	void q_qc_divide				(Quat &r, Quat const &x, Quat const &y)	{r=x/(Comp)y;}
	void q_qq_divide				(Quat &r, Quat const &x, Quat const &y)	{r=x/y;}
	
	void r_rr_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=(Real)x/(Real)y; r=t==floor(t);}//rc_divides, rq_divides: applied to each component
	void r_rc_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=(Comp)x/(Real)y; r=t==floor(t);}
	void r_rq_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=x/(Real)y; r=t==floor(t);}
	void r_cr_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=(Real)x/(Comp)y; r=t==floor(t);}
	void r_cc_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=(Comp)x/(Comp)y; r=t==floor(t);}
	void r_cq_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=x/(Comp)y; r=t==floor(t);}
	void r_qr_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=(Real)x/y; r=t==floor(t);}
	void r_qc_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=(Comp)x/y; r=t==floor(t);}
	void r_qq_logic_divides			(Quat &r, Quat const &y, Quat const &x)	{auto t=x/y; r=t==floor(t);}

	void r_rr_power_real			(Quat &r, Quat const &x, Quat const &y)	{r=pow((Real)x, floor((Real)y));}//trunc
	void c_cr_power_real			(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x^floor((Real)y);}
	void q_qr_power_real			(Quat &r, Quat const &x, Quat const &y)	{r=x^floor((Real)y);}

	void c_cr_pow					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x^(Real)y;}
	void c_cc_pow					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x^(Comp)y;}
	void q_cq_pow					(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x^y;}
	void q_qr_pow					(Quat &r, Quat const &x, Quat const &y)	{r=x^(Real)y;}
	void q_qc_pow					(Quat &r, Quat const &x, Quat const &y)	{r=x^(Comp)y;}
	void q_qq_pow					(Quat &r, Quat const &x, Quat const &y)	{r=x^y;}

	void  c_c_ln					(Quat &r, Quat const &x)				{r=log((Comp)x);}
	void  q_q_ln					(Quat &r, Quat const &x)				{r=log(x);}
	
	void  c_c_log					(Quat &r, Quat const &x)				{r=log((Comp)x)*m_inv_ln10;}
	void  q_q_log					(Quat &r, Quat const &x)				{r=log(x)*m_inv_ln10;}
	void c_cr_log					(Quat &r, Quat const &x, Quat const &y)	{r=log((Comp)x)/log(Comp(y, 0));}
	void c_cc_log					(Quat &r, Quat const &x, Quat const &y)	{r=log((Comp)x)/log((Comp)y);}
	void q_cq_log					(Quat &r, Quat const &x, Quat const &y)	{r=log((Comp)x)/log(y);}
	void q_qc_log					(Quat &r, Quat const &x, Quat const &y)	{r=log(x)/log((Comp)y);}
	void q_qq_log					(Quat &r, Quat const &x, Quat const &y)	{r=log(x)/log(y);}
	
	Comp tetrate(Real const &x, Real const &y)
	{
		Real log_x=log(x), ry=y;
		int h=floor(ry-.5).toLong()+1;//rounded real part
		Comp t(ry-h, 0);//real-round(real)+i*imag
		{
			auto q=sqrt(log_x);
			t=(t+1.)*(
				(
					log_x>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
						+t*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
								+t*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
										+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
						+t*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+m_ln2)/log_x)
								+t*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
										+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
									)
							)
				)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	Comp tetrate(Comp const &x, Real const &y)
	{
		Comp log_x=log(x);
		Real ry=y;
		int h=floor(ry-.5).toLong()+1;//rounded real part
		Comp t(ry-h, 0);//real-round(real)+i*imag
		{
			auto q=sqrt(log_x);
			t=(t+1.)*(
				(
					log_x.r>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
						+t*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
								+t*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
										+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
						+t*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+m_ln2)/log_x)
								+t*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
										+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
									)
							)
				)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	Quat tetrate(Quat const &x, Real const &y)
	{
		Quat qx=x;
		Real ry=y;
		if(ry<-1)
		{
			int steps=abs(ry).toLong();
			Quat t(ry-floor(ry), 0, 0, 0), lrx=log(qx);
			for(int k=0;k<steps;++k)
				t=log(t)/lrx;
			return t;
		}
		else if(ry<=0)
			return Quat(1+ry, 0, 0, 0);
		else
		{
			int h=ry.toLong()+1;
			Quat t(ry-floor(ry), 0, 0, 0);
			for(int k=0;k<h;++k)
				t=qx^t;
			return t;
		}
	}
	Comp tetrate(Real const &x, Comp const &y)
	{
		Real log_x=log((Real)x);
		Comp cy=y;
		int h=floor(cy.r-.5).toLong()+1;//rounded real part
		Comp t(cy-(Real)h);//real-round(real)+i*imag
		{
			auto q=sqrt(log_x);
			t=(t+1.)*(
				(
					log_x>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
						+t*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
								+t*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
										+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
						+t*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+m_ln2)/log_x)
								+t*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
										+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
									)
							)
				)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	Comp tetrate(Comp const &x, Comp const &y)
	{
		Comp log_x=log(x), cy=y;
	//	if(log_x.r<.03)//abs(log_x)<1.03045453395352
	//	{
	//		if(cy.r<-1.)
	//			return -30.;
	//		return 1.;
	//	}
		int h=floor(cy.r-.5).toLong()+1;
		Comp t(cy-(Real)h);//real-round(real)+i*imag
		{
		//	bool unassigned=true;
		//	if(log_x.r<.001)//abs(log_x)<1.00100050016671
		//	{
		//		if(t.r>-1)//real-round(real)>-1
		//			unassigned=false, t=1.;
		//		else if(t.r<-1)
		//			unassigned=false, t=-990.;
		//	}
		//	if(unassigned)
		//	{
				Comp q=sqrt(log_x);
				t=(t+1.)*(
					(
						log_x.r>m_ln2?
							(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
							+t*	(
									(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
									+t*	(
											(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
											+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
										)
								)
						:
							(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
							+t*	(
									(1.1-2.608785958462561*(1.-0.6663562294911147*sqrt(log_x))*sqrt(log_x)-(-0.625+m_ln2)/log_x)
									+t*	(
											(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
											+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
										)
								)
					)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		//	}
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	void c_rr_tetrate				(Quat &r, Quat const &x, Quat const &y)	{r=tetrate((Real)x, (Real)y);}
	void c_rc_tetrate				(Quat &r, Quat const &x, Quat const &y)	{r=tetrate((Real)x, (Comp)y);}
	void c_cr_tetrate				(Quat &r, Quat const &x, Quat const &y)	{r=tetrate((Comp)x, (Real)y);}
	void c_cc_tetrate				(Quat &r, Quat const &x, Quat const &y)	{r=tetrate((Comp)x, (Comp)y);}
	void q_qr_tetrate				(Quat &r, Quat const &x, Quat const &y)	{r=tetrate((Quat)x, (Real)y);}
	
	class Tetrate//http://en.citizendium.org/wiki/Fit1.cin
	{
		static Comp _fit1(Comp &x, Comp &y)
		{
			if(x.r<.001)
			{
				if(y.r>-1)
					return Comp(1, 0);
				if(y.r<-1)
					return Comp(-990, 0);
			}
			Comp q=sqrt(x);
			return (y+1.)*(
				(
					x.r>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/x-1.)
						+y*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/x-1.)
								+y*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/x-1.)
										+y*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+x*61.6566))))/(1.-8.1192*q+37.087*x)+(-131./192.+m_ln2)/x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/x)
						+y*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*sqrt(x))*sqrt(x)-(-0.625+m_ln2)/x)
								+y*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*x)*q/(1.+ 4.240467556480155*x)+(-2./3+m_ln2)/x)
										+y*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*x)*q/(1.+4.95715636660691*q + 7.70233216637738*x)-(-131./192.+m_ln2)/x)
									)
							)
				)*y+1.)+log(y+2.)/x-m_ln2/x*(1.+y);
		}
	public:
		static Comp fit1(Comp x, Comp const &y)//[sic]
		{
			x=log(x);
			if(x.r<.03)
			{
				if(y.r<-1.)
					return Comp(-30, 0);
				return Comp(1, 0);
			}
			int h=floor(y.r-.5).toLong()+1;
			Comp result=_fit1(x, y-(Real)h);
			for(;h>0;--h)
				result=exp(x*result);
			for(;h<0;++h)
				result=log(result)/x;
			return result;
		}
	};
	inline Comp pentate(Real const &x, Real const &y)
	{
		long long h=y.toLLong();
	//	long long h=y.r!=y.r||y.r<-ll_max||y.r>ll_max?0:long long(y);
	//	long long h=std::isnan(y.r)||std::isinf(y.r)?0:long long(y);
		if(h<-2)	return Comp(_HUGE, 0);//1/::sin(0);
		if(h==-2)	return Comp(-1, 0);
		if(h==-1)	return Comp(0, 0);
		if(h==0)	return Comp(1, 0);
		if(h==1)	return Comp(x, 0);
		Real rx=x;
		Comp result(rx, 0);
		for(int k=0;k<h;++k)
			result=Tetrate::fit1(Comp(rx, 0), result);
		return result;
	}
	inline Comp pentate(Comp const &x, Real const &y)
	{
		long long h=y.toLLong();
		if(h<-2)	return Comp(_HUGE, 0);//1/::sin(0);
		if(h==-2)	return Comp(-1, 0);
		if(h==-1)	return Comp(0, 0);
		if(h==0)	return Comp(1, 0);
		if(h==1)	return x;
		Comp result=(Comp)x;
		for(int k=0;k<h;++k)
			result=Tetrate::fit1(x, result);
		return result;
	}
	void c_rr_pentate				(Quat &r, Quat const &x, Quat const &y)	{r=pentate((Real)x, (Real)y);}
	void c_cr_pentate				(Quat &r, Quat const &x, Quat const &y)	{r=pentate((Comp)x, (Real)y);}

	void  r_r_bitwise_shift_left_l	(Quat &r, Quat const &x)				{r=exp(floor((Real)x)*m_ln2);}	//<<x = 2^floor(x)
	void  c_c_bitwise_shift_left_l	(Quat &r, Quat const &x)				{r=exp(floor((Comp)x)*m_ln2);}
	void  q_q_bitwise_shift_left_l	(Quat &r, Quat const &x)				{r=exp(floor(x)*m_ln2);}
	void  r_r_bitwise_shift_left_r	(Quat &r, Quat const &x)				{r=x.r+x.r;}					//x<< = 2x
	void  c_c_bitwise_shift_left_r	(Quat &r, Quat const &x)				{Comp cx=(Comp)x; r=cx+cx;}
	void  q_q_bitwise_shift_left_r	(Quat &r, Quat const &x)				{Quat qx=x; r=qx+qx;}
	void r_rr_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*exp(floor((Real)y)*m_ln2);}//x<<y = x*2^y
	void c_rc_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*exp(floor((Comp)y)*m_ln2);}
	void q_rq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*exp(floor(y)*m_ln2);}
	void c_cr_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*exp(floor((Real)y)*m_ln2);}
	void c_cc_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*exp(floor((Comp)y)*m_ln2);}
	void q_cq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*exp(floor(y)*m_ln2);}
	void q_qr_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(floor((Real)y)*m_ln2);}
	void q_qc_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(floor((Comp)y)*m_ln2);}
	void q_qq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(floor(y)*m_ln2);}

	void  r_r_bitwise_shift_right_l	(Quat &r, Quat const &x)				{r=exp(-floor((Real)x)*m_ln2);}
	void  c_c_bitwise_shift_right_l	(Quat &r, Quat const &x)				{r=exp(-floor((Comp)x)*m_ln2);}
	void  q_q_bitwise_shift_right_l	(Quat &r, Quat const &x)				{r=exp(-floor(x)*m_ln2);}
	void  r_r_bitwise_shift_right_r	(Quat &r, Quat const &x)				{r=(Real)x*0.5;}
	void  c_c_bitwise_shift_right_r	(Quat &r, Quat const &x)				{r=(Comp)x*0.5;}
	void  q_q_bitwise_shift_right_r	(Quat &r, Quat const &x)				{r=(Quat)x*0.5;}
	void r_rr_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*exp(-floor((Real)y)*m_ln2);}
	void c_rc_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*exp(-floor((Comp)y)*m_ln2);}
	void q_rq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Real)x*exp(-floor(y)*m_ln2);}
	void c_cr_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*exp(-floor((Real)y)*m_ln2);}
	void c_cc_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*exp(-floor((Comp)y)*m_ln2);}
	void q_cq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x*exp(-floor(y)*m_ln2);}
	void q_qr_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(-floor((Real)y)*m_ln2);}
	void q_qc_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(-floor((Comp)y)*m_ln2);}
	void q_qq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(-floor(y)*m_ln2);}

	inline Real bitwise_not(Real const &x){return isnan(x)|isinf(x)?0:~x.toLLong();}
	void  r_r_bitwise_not			(Quat &r, Quat const &x)				{r=bitwise_not(x.r);}
	void  c_c_bitwise_not			(Quat &r, Quat const &x)				{r=Comp(bitwise_not(x.r), bitwise_not(x.i));}
	void  q_q_bitwise_not			(Quat &r, Quat const &x)				{r=Quat(bitwise_not(x.r), bitwise_not(x.i), bitwise_not(x.j), bitwise_not(x.k));}

	inline Real bitwise_and(Real const &x){return isnan(x)|isinf(x)?0:!~x.toLLong();}
	inline Real bitwise_and(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:Real(x.toULLong()&y.toULLong());}
	inline long long bitwise_and_ll(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:x.toULLong()&y.toULLong();}
	void  r_r_bitwise_and			(Quat &r, Quat const &x)				{r=bitwise_and(x.r);}
	void  c_c_bitwise_and			(Quat &r, Quat const &x)				{r=Comp(bitwise_and(x.r), bitwise_and(x.i));}
	void  q_q_bitwise_and			(Quat &r, Quat const &x)				{r=Quat(bitwise_and(x.r), bitwise_and(x.i), bitwise_and(x.j), bitwise_and(x.k));}
	void r_rr_bitwise_and			(Quat &r, Quat const &x, Quat const &y)	{r=bitwise_and(x.r, y.r);}
	void c_rc_bitwise_and			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_and(x.r, y.r), bitwise_and(x.r, y.i));}
	void q_rq_bitwise_and			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_and(x.r, y.r), bitwise_and(x.r, y.i), bitwise_and(x, y.j), bitwise_and(x, y.k));}
	void c_cr_bitwise_and			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_and(x.r, y.r), bitwise_and(x.i, y.r));}
	void c_cc_bitwise_and			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i);
		r=Comp(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr));
	}
	void q_cq_bitwise_and			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i),
			xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j),
			xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xr_yj-xi_yk), Real(xr_yk+xi_yj));
	}
	void q_qr_bitwise_and			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_and(x.r, y.r), bitwise_and(x.i, y.r), bitwise_and(x.j, y.r), bitwise_and(x.k, y.r));}
	void q_qc_bitwise_and			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xj_yr+xk_yi), Real(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_and			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i),
			xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j), xj_yj=bitwise_and_ll(x.j, y.j), xk_yj=bitwise_and_ll(x.k, y.j),
			xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k), xj_yk=bitwise_and_ll(x.j, y.k), xk_yk=bitwise_and_ll(x.k, y.k);
		r=Quat(Real(xr_yr-xi_yi-xj_yj-xk_yk), Real(xr_yi+xi_yr+xj_yk-xk_yj), Real(xj_yj-xi_yk+xj_yr+xk_yi), Real(xr_yk+xi_yj-xj_yi+xk_yr));
	}

	inline Real bitwise_nand(Real const &x){return isnan(x)|isinf(x)?0:!x.toDouble();}
	inline Real bitwise_nand(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:~(x.toULLong()&y.toULLong());}
	inline long long bitwise_nand_ll(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:~(x.toULLong()&y.toULLong());}
	void  r_r_bitwise_nand			(Quat &r, Quat const &x)				{r=bitwise_nand(x.r);}
	void  c_c_bitwise_nand			(Quat &r, Quat const &x)				{r=Comp(bitwise_nand(x.r), bitwise_nand(x.i));}
	void  q_q_bitwise_nand			(Quat &r, Quat const &x)				{r=Quat(bitwise_nand(x.r), bitwise_nand(x.i), bitwise_nand(x.j), bitwise_nand(x.k));}
	void r_rr_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)	{r=bitwise_nand(x.r, y.r);}
	void c_rc_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_nand(x.r, y.r), bitwise_nand(x.r, y.i));}
	void q_rq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_nand(x.r, y.r), bitwise_nand(x.r, y.i), bitwise_nand(x, y.j), bitwise_nand(x, y.k));}
	void c_cr_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_nand(x.r, y.r), bitwise_nand(x.i, y.r));}
	void c_cc_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i);
		r=Comp(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr));
	}
	void q_cq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i),
			xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j),
			xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xr_yj-xi_yk), Real(xr_yk+xi_yj));
	}
	void q_qr_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_nand(x.r, y.r), bitwise_nand(x.i, y.r), bitwise_nand(x.j, y.r), bitwise_nand(x.k, y.r));}
	void q_qc_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xj_yr+xk_yi), Real(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i),
			xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j), xj_yj=bitwise_nand_ll(x.j, y.j), xk_yj=bitwise_nand_ll(x.k, y.j),
			xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k), xj_yk=bitwise_nand_ll(x.j, y.k), xk_yk=bitwise_nand_ll(x.k, y.k);
		r=Quat(Real(xr_yr-xi_yi-xj_yj-xk_yk), Real(xr_yi+xi_yr+xj_yk-xk_yj), Real(xj_yj-xi_yk+xj_yr+xk_yi), Real(xr_yk+xi_yj-xj_yi+xk_yr));
	}
	
	inline Real bitwise_or(Real const &x){return isnan(x)|isinf(x)?0:x.toULLong()!=0;}
	inline Real bitwise_or(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:Real(x.toULLong()|y.toULLong());}
	inline long long bitwise_or_ll_c(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:~x.toULLong()&~y.toULLong();}
	void  r_r_bitwise_or			(Quat &r, Quat const &x)				{r=bitwise_or(x.r);}
	void  c_c_bitwise_or			(Quat &r, Quat const &x)				{r=Comp(bitwise_or(x.r), bitwise_or(x.i));}
	void  q_q_bitwise_or			(Quat &r, Quat const &x)				{r=Quat(bitwise_or(x.r), bitwise_or(x.i), bitwise_or(x.j), bitwise_or(x.k));}
	void r_rr_bitwise_or			(Quat &r, Quat const &x, Quat const &y)	{r=bitwise_or(x.r, y.r);}
	void c_rc_bitwise_or			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_or(x.r, y.r), bitwise_or(x.r, y.i));}
	void q_rq_bitwise_or			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_or(x.r, y.r), bitwise_or(x.r, y.i), bitwise_or(x, y.j), bitwise_or(x, y.k));}
	void c_cr_bitwise_or			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_or(x.r, y.r), bitwise_or(x.i, y.r));}
	void c_cc_bitwise_or			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i);
		r=Comp(~(xr_yr-xi_yi), ~(xr_yi+xi_yr));
	}
	void q_cq_bitwise_or			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i),
			xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j),
			xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k);
		r=Quat(~(xr_yr-xi_yi), ~(xr_yi+xi_yr), ~(xr_yj-xi_yk), ~(xr_yk+xi_yj));
	}
	void q_qr_bitwise_or			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_or(x.r, y.r), bitwise_or(x.i, y.r), bitwise_or(x.j, y.r), bitwise_or(x.k, y.r));}
	void q_qc_bitwise_or			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i);
		r=Quat(~(xr_yr-xi_yi), ~(xr_yi+xi_yr), ~(xj_yr+xk_yi), ~(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_or			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i),
			xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j), xj_yj=bitwise_or_ll_c(x.j, y.j), xk_yj=bitwise_or_ll_c(x.k, y.j),
			xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k), xj_yk=bitwise_or_ll_c(x.j, y.k), xk_yk=bitwise_or_ll_c(x.k, y.k);
		r=Quat(~(xr_yr-xi_yi-xj_yj-xk_yk), ~(xr_yi+xi_yr+xj_yk-xk_yj), ~(xj_yj-xi_yk+xj_yr+xk_yi), ~(xr_yk+xi_yj-xj_yi+xk_yr));
	}
	
	inline Real bitwise_nor(Real const &x){return isnan(x)|isinf(x)?0:!x.toULLong();}
	inline Real bitwise_nor(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(x)|isinf(x)?0:~(x.toULLong()|y.toULLong());}
	inline long long bitwise_nor_ll_c(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(x)|isinf(x)?0:~x.toULLong()|~y.toULLong();}
	void  r_r_bitwise_nor			(Quat &r, Quat const &x)				{r=bitwise_nor(x.r);}
	void  c_c_bitwise_nor			(Quat &r, Quat const &x)				{r=Comp(bitwise_nor(x.r), bitwise_nor(x.i));}
	void  q_q_bitwise_nor			(Quat &r, Quat const &x)				{r=Quat(bitwise_nor(x.r), bitwise_nor(x.i), bitwise_nor(x.j), bitwise_nor(x.k));}
	void r_rr_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)	{r=bitwise_nor(x.r, y.r);}
	void c_rc_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_nor(x.r, y.r), bitwise_nor(x.r, y.i));}
	void q_rq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_nor(x.r, y.r), bitwise_nor(x.r, y.i), bitwise_nor(x.r, y.j), bitwise_nor(x.r, y.k));}
	void c_cr_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_nor(x.r, y.r), bitwise_nor(x.i, y.r));}
	void c_cc_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i);
		r=Comp(~(xr_yr-xi_yi), ~(xr_yi+xi_yr));
	}
	void q_cq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i),
			xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j),
			xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k);
		r=Quat((Real)~(xr_yr-xi_yi), (Real)~(xr_yi+xi_yr), (Real)~(xr_yj-xi_yk), (Real)~(xr_yk+xi_yj));
	}
	void q_qr_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_nor(x.r, y.r), bitwise_nor(x.i, y.r), bitwise_nor(x.j, y.r), bitwise_nor(x.k, y.r));}
	void q_qc_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i);
		r=Quat((Real)~(xr_yr-xi_yi), (Real)~(xr_yi+xi_yr), (Real)~(xj_yr+xk_yi), (Real)~(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i),
			xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j), xj_yj=bitwise_nor_ll_c(x.j, y.j), xk_yj=bitwise_nor_ll_c(x.k, y.j),
			xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k), xj_yk=bitwise_nor_ll_c(x.j, y.k), xk_yk=bitwise_nor_ll_c(x.k, y.k);
		r=Quat((Real)~(xr_yr-xi_yi-xj_yj-xk_yk), (Real)~(xr_yi+xi_yr+xj_yk-xk_yj), (Real)~(xj_yj-xi_yk+xj_yr+xk_yi), (Real)~(xr_yk+xi_yj-xj_yi+xk_yr));
	}
	
	inline Real bitwise_xor(Real const &x){return isnan(x)|isinf(x)?0:bitwise_xor(x);}
	inline Real bitwise_xor(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:Real(x.toULLong()^y.toULLong());}
	void  r_r_bitwise_xor			(Quat &r, Quat const &x)				{r=bitwise_xor(x.r);}
	void  c_c_bitwise_xor			(Quat &r, Quat const &x)				{r=Comp(bitwise_xor(x.r), bitwise_xor(x.i));}
	void  q_q_bitwise_xor			(Quat &r, Quat const &x)				{r=Quat(bitwise_xor(x.r), bitwise_xor(x.i), bitwise_xor(x.j), bitwise_xor(x.k));}
	void r_rr_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=bitwise_xor(x.r, y.r);}
	void c_rc_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_xor(x.r, y.r), y.i);}
	void q_rq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), y.i, y.j, y.k);}
	void c_cr_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_xor(x.r, y.r), x.i);}
	void c_cc_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i));}
	void q_cq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), y.j, y.k);}
	void q_qr_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), x.i, x.j, x.k);}
	void q_qc_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), x.j, x.k);}
	void q_qq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), bitwise_xor(x.j, y.j), bitwise_xor(x.k, y.k));}
	
	inline Real bitwise_xnor(Real const &x){return isnan(x)|isinf(x)?0:!bitwise_xor(x).toDouble();}
	inline Real bitwise_xnor(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:(Real)~(x.toULLong()^y.toULLong());}
	void  r_r_bitwise_xnor			(Quat &r, Quat const &x)				{r=bitwise_xnor(x.r);}
	void  c_c_bitwise_xnor			(Quat &r, Quat const &x)				{r=Comp(bitwise_xnor(x.r), bitwise_xnor(x.i));}
	void  q_q_bitwise_xnor			(Quat &r, Quat const &x)				{r=Quat(bitwise_xnor(x.r), bitwise_xnor(x.i), bitwise_xnor(x.j), bitwise_xnor(x.k));}
	void r_rr_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=bitwise_xnor(x.r, y.r);}
	void c_rc_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_xnor(x.r, y.r), y.i);}
	void q_rq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), y.i, y.j, y.k);}
	void c_cr_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_xnor(x.r, y.r), x.i);}
	void c_cc_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Comp(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i));}
	void q_cq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), y.j, y.k);}
	void q_qr_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), x.i, x.j, x.k);}
	void q_qc_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), x.j, x.k);}
	void q_qq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), bitwise_xnor(x.j, y.j), bitwise_xnor(x.k, y.k));}
	
	void  r_r_logic_equal			(Quat &r, Quat const &x)				{r=x==0;}
	void  r_c_logic_equal			(Quat &r, Quat const &x)				{r=!x.c_is_true();}
	void  r_q_logic_equal			(Quat &r, Quat const &x)				{r=!x.q_is_true();}
	void r_rr_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=x.r==y.r;}
	void r_rc_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=x.r==(Comp)y;}
	void r_rq_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=x.r==y;}
	void r_cr_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x==y.r;}
	void r_cc_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x==(Comp)y;}
	void r_cq_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x==y;}
	void r_qr_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=x==y.r;}
	void r_qc_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=x==(Comp)y;}
	void r_qq_logic_equal			(Quat &r, Quat const &x, Quat const &y)	{r=x==y;}
	
	void  r_r_logic_not_equal		(Quat &r, Quat const &x)				{r=x.r!=0;}
	void  r_c_logic_not_equal		(Quat &r, Quat const &x)				{r=x.c_is_true();}
	void  r_q_logic_not_equal		(Quat &r, Quat const &x)				{r=x.q_is_true();}
	void r_rr_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r!=y.r;}
	void r_rc_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r!=(Comp)y;}
	void r_rq_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r!=y;}
	void r_cr_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x!=y.r;}
	void r_cc_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x!=(Comp)y;}
	void r_cq_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x!=y;}
	void r_qr_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x!=y.r;}
	void r_qc_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x!=(Comp)y;}
	void r_qq_logic_not_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x!=y;}
	
	void  r_r_logic_less_l			(Quat &r, Quat const &x)				{r=0<x.r;}
	void  r_c_logic_less_l			(Quat &r, Quat const &x)				{r=0<x.r;}
	void  r_q_logic_less_l			(Quat &r, Quat const &x)				{r=0<x.r;}
	void  r_r_logic_less_r			(Quat &r, Quat const &x)				{r=x.r<0;}
	void  r_c_logic_less_r			(Quat &r, Quat const &x)				{r=x.r<0;}
	void  r_q_logic_less_r			(Quat &r, Quat const &x)				{r=x.r<0;}
	void r_rr_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_rc_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_rq_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_cr_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_cc_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_cq_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_qr_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_qc_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	void r_qq_logic_less			(Quat &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	
	void  r_r_logic_less_equal_l	(Quat &r, Quat const &x)				{r=0<=x.r;}
	void  r_c_logic_less_equal_l	(Quat &r, Quat const &x)				{r=0<=x.r;}
	void  r_q_logic_less_equal_l	(Quat &r, Quat const &x)				{r=0<=x.r;}
	void  r_r_logic_less_equal_r	(Quat &r, Quat const &x)				{r=x.r<=0;}
	void  r_c_logic_less_equal_r	(Quat &r, Quat const &x)				{r=x.r<=0;}
	void  r_q_logic_less_equal_r	(Quat &r, Quat const &x)				{r=x.r<=0;}
	void r_rr_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_rc_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_rq_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_cr_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_cc_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_cq_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_qr_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_qc_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_qq_logic_less_equal		(Quat &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	
	void  r_r_logic_greater_l		(Quat &r, Quat const &x)				{r=0>x.r;}
	void  r_c_logic_greater_l		(Quat &r, Quat const &x)				{r=0>x.r;}
	void  r_q_logic_greater_l		(Quat &r, Quat const &x)				{r=0>x.r;}
	void  r_r_logic_greater_r		(Quat &r, Quat const &x)				{r=x.r>0;}
	void  r_c_logic_greater_r		(Quat &r, Quat const &x)				{r=x.r>0;}
	void  r_q_logic_greater_r		(Quat &r, Quat const &x)				{r=x.r>0;}
	void r_rr_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_rc_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_rq_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_cr_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_cc_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_cq_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_qr_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_qc_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	void r_qq_logic_greater			(Quat &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	
	void  r_r_logic_greater_equal_l	(Quat &r, Quat const &x)				{r=0>=x.r;}
	void  r_c_logic_greater_equal_l	(Quat &r, Quat const &x)				{r=0>=x.r;}
	void  r_q_logic_greater_equal_l	(Quat &r, Quat const &x)				{r=0>=x.r;}
	void  r_r_logic_greater_equal_r	(Quat &r, Quat const &x)				{r=x.r>=0;}
	void  r_c_logic_greater_equal_r	(Quat &r, Quat const &x)				{r=x.r>=0;}
	void  r_q_logic_greater_equal_r	(Quat &r, Quat const &x)				{r=x.r>=0;}
	void r_rr_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_rc_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_rq_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_cr_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_cc_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_cq_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_qr_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_qc_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_qq_logic_greater_equal	(Quat &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	
	void  r_r_logic_not				(Quat &r, Quat const &x)				{r=x.r==0;}
	void  r_c_logic_not				(Quat &r, Quat const &x)				{r=(x.r==0)&(x.i==0);}
	void  r_q_logic_not				(Quat &r, Quat const &x)				{r=(x.r==0)&(x.i==0)&(x.j==0)&(x.k==0);}
	
	void r_rr_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0&y.r==0;}
	void r_rc_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0&y.c_is_true();}
	void r_rq_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0&y.q_is_true();}
	void r_cr_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()&y.r==0;}
	void r_cc_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()&y.c_is_true();}
	void r_cq_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()&y.q_is_true();}
	void r_qr_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()&y.r==0;}
	void r_qc_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()&y.c_is_true();}
	void r_qq_logic_and				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()&y.q_is_true();}
	
	void r_rr_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0|y.r==0;}
	void r_rc_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0|y.c_is_true();}
	void r_rq_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0|y.q_is_true();}
	void r_cr_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()|y.r==0;}
	void r_cc_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()|y.c_is_true();}
	void r_cq_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()|y.q_is_true();}
	void r_qr_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()|y.r==0;}
	void r_qc_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()|y.c_is_true();}
	void r_qq_logic_or				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()|y.q_is_true();}
	
	void r_rr_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0^y.r==0;}
	void r_rc_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0^y.c_is_true();}
	void r_rq_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0^y.q_is_true();}
	void r_cr_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()^y.r==0;}
	void r_cc_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()^y.c_is_true();}
	void r_cq_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()^y.q_is_true();}
	void r_qr_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()^y.r==0;}
	void r_qc_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()^y.c_is_true();}
	void r_qq_logic_xor				(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()^y.q_is_true();}
	
	void r_rr_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0?x.r:y.r;}
	void c_rc_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0?Comp(x.r, 0):y;}
	void q_rq_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.r==0?Quat(x.r, 0, 0, 0):y;}
	void c_cr_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?x:Comp(y.r, 0);}
	void c_cc_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?x:y;}
	void q_cq_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?Quat(x):y;}
	void q_qr_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?x:Quat(y.r, 0, 0, 0);}
	void q_qc_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?x:Quat(y.r, y.i);}
	void q_qq_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?x:y;}

	void  r_r_percent				(Quat &r, Quat const &x)				{r=x.r*0.01;}
	void  c_c_percent				(Quat &r, Quat const &x)				{r=(Comp)x*0.01;}
	void  q_q_percent				(Quat &r, Quat const &x)				{r=(Quat)x*0.01;}
	
	void r_rr_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x.r%y.r;}
	void c_rc_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x.r%(Comp)y;}
	void q_rq_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x.r%y;}
	void c_cr_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x%y.r;}
	void c_cc_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x%(Comp)y;}
	void q_cq_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=(Comp)x%y;}
	void q_qr_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x%y.r;}
	void q_qc_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x%(Comp)y;}
	void q_qq_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x%y;}

	void  r_r_sgn					(Quat &r, Quat const &x)				{r=(x.r>0)-(x.r<0);}
	void  c_c_sgn					(Quat &r, Quat const &x)				{Comp cx=(Comp)x; Real mag=abs(cx); r=mag.toDouble()?cx/mag:Comp();}
	void  q_q_sgn					(Quat &r, Quat const &x)				{Quat qx=x; Real mag=abs(qx); r=mag.toDouble()?qx/mag:Quat();}
	
	inline Comp	sq(Comp const &x){Real ri=x.r*x.i; return Comp(x.r*x.r-x.i*x.i, ri+ri);}
	inline Quat	sq(Quat const &x)
	{
		auto _2r=x.r+x.r;
		return Quat(x.r*x.r-x.i*x.i-x.j*x.j-x.k*x.k, x.i*_2r, x.j*_2r, x.k*_2r);
	}
	void  r_r_sq					(Quat &r, Quat const &x)				{r=x.r*x.r;}
	void  c_c_sq					(Quat &r, Quat const &x)				{r=sq((Comp)x);}
	void  q_q_sq					(Quat &r, Quat const &x)				{r=sq(x);}
	
	void  c_c_sqrt					(Quat &r, Quat const &x)				{r=sqrt((Comp)x);}
	void  q_q_sqrt					(Quat &r, Quat const &x)				{r=sqrt((Quat)x);}

	void  r_r_invsqrt				(Quat &r, Quat const &x)				{r=rec_sqrt(x.r);}

	void  r_r_cbrt					(Quat &r, Quat const &x)				{r=cbrt(x.r);}
	void  c_c_cbrt					(Quat &r, Quat const &x)				{r=exp(log((Comp)x)/3);}//optimize
	void  q_q_cbrt					(Quat &r, Quat const &x)
	{
		Real t=x.i*x.i+x.j*x.j+x.k*x.k;
		if(t.toDouble())
		{
			Real absx=sqrt(x.r*x.r+t);
			t=acos(x.r/absx)/(sqrt(t)*3);
			Quat logx_3((1./3)*log(absx), t*x.i, t*x.j, t*x.k);
			r=exp(logx_3);
		}
		else
			r=Quat(cbrt(x.r));
	}

	void  r_r_gauss					(Quat &r, Quat const &x)				{r=exp(-x.r*x.r);}
	void  c_c_gauss					(Quat &r, Quat const &x)				{r=exp(-sq((Comp)x));}
	void  q_q_gauss					(Quat &r, Quat const &x)				{r=exp(-sq(x));}

	void  r_r_erf					(Quat &r, Quat const &x)				{r=erf(x.r);}

	void  r_r_zeta					(Quat &r, Quat const &x)				{r=zeta(x.r);}
	
	const Comp m_i(0, 1);
	__forceinline Quat sgnu(Quat const &x){return Quat(0, x.i, x.j, x.k);}
	__forceinline Quat acosh(Quat const &x){return log(x+sqrt(sq(x)-1));}
	__forceinline Comp asinh(Comp const &x){return log(x+sqrt(sq(x)+1));}
	__forceinline Quat asinh(Quat const &x){return log(x+sqrt(sq(x)+1));}
	__forceinline Real sinhc(Real const &x){return sinh(x)/x;}
	__forceinline Comp cos(Comp const &x)
	{
		Comp exp_ix=exp(m_i*x);//sincos, exp
		return 0.5*exp_ix+0.5/exp_ix;
	}
	__forceinline Quat cos(Quat const &x)
	{
		Real z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k), sin_xr, cos_xr;//boost::math: sqrt, sincos, sh, ch
		sin_cos(sin_xr, cos_xr, x.r);
		Real w=-sin_xr*sinhc(z);
		return Quat(cos_xr*cosh(z), w*x.i, w*x.j, w*x.k);
	}
	__forceinline Comp sin(Comp const &x)
	{
		Comp exp_ix=exp(m_i*x);
		return 0.5*exp_ix-0.5/exp_ix;
	//	return (exp(m_i*x)-exp(-m_i*x))*(-m_half*m_i);
	}
	__forceinline Quat sin(Quat const &x)//boost::math
	{
		Real z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
		Real sin_xr, cos_xr;
		sin_cos(sin_xr, cos_xr, x.r);
		Real w=-cos_xr*sinhc(z);
		return Quat(sin_xr*cosh(z), w*x.i, w*x.j, w*x.k);
	}
	namespace gamma//http://en.wikipedia.org/wiki/Lanczos_approximation
	{
		const double g=7, p[]={0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	}
	Comp			tgamma		(Comp const &x)
	{
		using namespace gamma;
		if(x.r<.5)
		{
			Comp t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			Comp t2=g+.5-x;
			return m_pi/(sin(m_pi*x)*m_sqrt_2pi*(t2^0.5-x)*exp(-t2)*t1);
		}
		else
		{
			Comp t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			Comp t2=g+.5+x-1.;
			return m_sqrt_2pi*(t2^0.5+x-1)*exp(-t2)*t1;
		}
	}
	Quat			tgamma		(Quat	const &x)
	{
		using namespace gamma;
		if(x.r<.5)
		{
			Quat t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			Quat t2=g+.5-x;
			return m_pi/(sin(m_pi*x)*m_sqrt_2pi*(t2^0.5-x)*exp(-t2)*t1);
		}
		else
		{
			Quat t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			Quat t2=g+.5+x-1.;
			return m_sqrt_2pi*(t2^0.5+x-1.)*exp(-t2)*t1;
		}
	}
	void  r_r_tgamma				(Quat &r, Quat const &x)				{r=tgamma(x.r);}
	void  c_c_tgamma				(Quat &r, Quat const &x)				{r=tgamma((Comp)x);}
	void  q_q_tgamma				(Quat &r, Quat const &x)				{r=tgamma(x);}
	double tgamma(double const &x, double const &y)//upper incomplete gamma, G(x, 0) = G(x)
	{
		try
		{
			if(x>0)
			{
				if(y>=0)
					return boost::math::tgamma(x, y);
			}
			else if(x==0)
				return _HUGE;
			long long cInf=0x7FF8000000000010;//x.r>0&&y.r<0||x.r<0
			return (double&)cInf;
			//if(x.r<0)
			//{
			//	long long cInf=0x7FF8000000000010;
			//	return (Real&)cInf;
			//}
			//else if(x.r==0)
			//	return _HUGE;
			//else if(y.r<0)
			//{
			//	long long cInf=0x7FF8000000000010;
			//	return (Real&)cInf;
			//}
			//return boost::math::tgamma((Real)x, (Real)y);
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
	void r_rr_tgamma				(Quat &r, Quat const &x, Quat const &y)	{r=tgamma(x.r.toDouble(), y.r.toDouble());}

	void  r_r_loggamma				(Quat &r, Quat const &x)				{r=lgamma(x.r);}

	void  r_r_factorial				(Quat &r, Quat const &x)				{r=tgamma(x.r+1);}
	void  c_c_factorial				(Quat &r, Quat const &x)				{r=tgamma((Comp)x+1.);}
	void  q_q_factorial				(Quat &r, Quat const &x)				{r=tgamma((Quat)x+1.);}

	inline Real permutation(Real const &x, Real const &y){return tgamma(x+1)/tgamma(x-y+1);}
	Comp permutation(Comp const &x, Comp const &y){return tgamma(x+1.)/tgamma(x-y+1.);}
	Quat permutation(Quat const &x, Quat const &y){return tgamma(x+1.)/tgamma(x-y+1.);}
	void  r_r_permutation			(Quat &r, Quat const &x)				{r=1;}
	void  c_c_permutation			(Quat &r, Quat const &x)				{r=Comp(1, 0);}
	void  q_q_permutation			(Quat &r, Quat const &x)				{r=Quat(1, 0, 0, 0);}
	void r_rr_permutation			(Quat &r, Quat const &x, Quat const &y)	{r=permutation(x.r, y.r);}
	void c_cr_permutation			(Quat &r, Quat const &x, Quat const &y)	{r=permutation((Comp)x, Comp(y.r, 0));}
	void c_cc_permutation			(Quat &r, Quat const &x, Quat const &y)	{r=permutation((Comp)x, (Comp)y);}
	void q_qq_permutation			(Quat &r, Quat const &x, Quat const &y)	{r=permutation(x, y);}
	
	Real combination(Real const &x, Real const &y){return tgamma(x+1)/(tgamma(x-y+1)*tgamma(y+1));}
	Comp combination(Comp const &x, Comp const &y){return tgamma(x+1)/(tgamma(x-y+1)*tgamma(y+1));}
	Quat combination(Quat const &x, Quat const &y){return tgamma(x+1)/(tgamma(x-y+1)*tgamma(y+1));}
	void  r_r_combination			(Quat &r, Quat const &x)				{r=1;}
	void  c_c_combination			(Quat &r, Quat const &x)				{r=Comp(1, 0);}
	void  q_q_combination			(Quat &r, Quat const &x)				{r=Quat(1, 0, 0, 0);}
	void r_rr_combination			(Quat &r, Quat const &x, Quat const &y)	{r=combination(x.r, y.r);}
	void c_cr_combination			(Quat &r, Quat const &x, Quat const &y)	{r=combination((Comp)x, Comp(y.r, 0));}
	void c_cc_combination			(Quat &r, Quat const &x, Quat const &y)	{r=combination((Comp)x, (Comp)y);}
	void q_qq_combination			(Quat &r, Quat const &x, Quat const &y)	{r=combination(x, y);}

	void  r_r_cos					(Quat &r, Quat const &x)				{r=cos(x.r);}
	void  c_c_cos					(Quat &r, Quat const &x)				{r=cos((Comp)x);}
	void  q_q_cos					(Quat &r, Quat const &x)				{r=cos(x);}

	inline Comp acos(Comp const &x){return -m_i*log(x+sqrt(sq(x)-1.));}
	inline Quat acos(Quat const &x){return -m_i*log(x+sqrt(sq(x)-1.));}
	void  c_c_acos					(Quat &r, Quat const &x)				{r=acos((Comp)x);}
	void  q_q_acos					(Quat &r, Quat const &x)				{r=acos(x);}

	inline Comp cosh(Comp const &x){return (exp(x)+exp(-x))*0.5;}
	inline Quat cosh(Quat const &x){return (exp(x)+exp(-x))*0.5;}
	void  r_r_cosh					(Quat &r, Quat const &x)				{r=cosh(x.r);}
	void  c_c_cosh					(Quat &r, Quat const &x)				{r=cosh((Comp)x);}
	void  q_q_cosh					(Quat &r, Quat const &x)				{r=cosh(x);}
	
	inline Comp acosh(Comp const &x){return log(x+sqrt(sq(x)-1.));}
	void  c_c_acosh					(Quat &r, Quat const &x)				{r=acosh((Comp)x);}
	void  q_q_acosh					(Quat &r, Quat const &x)				{r=acosh(x);}

	void  r_r_cosc					(Quat &r, Quat const &x)				{r=cos(x.r)/x.r;}
	void  c_c_cosc					(Quat &r, Quat const &x)				{Comp cx=(Comp)x; r=cos(cx)/cx;}
	void  q_q_cosc					(Quat &r, Quat const &x)				{Quat qx=x; r=cos(qx)/qx;}

	void  r_r_sec					(Quat &r, Quat const &x)				{r=1/cos(x.r);}
	void  c_c_sec					(Quat &r, Quat const &x)				{r=inv(cos((Comp)x));}
	void  q_q_sec					(Quat &r, Quat const &x)				{r=inv(cos(x));}

	void  c_c_asec					(Quat &r, Quat const &x)				{r=acos(inv((Comp)x));}
	void  q_q_asec					(Quat &r, Quat const &x)				{r=acos(inv(x));}

	void  r_r_sech					(Quat &r, Quat const &x)				{r=1/cosh(x.r);}
	void  c_c_sech					(Quat &r, Quat const &x)				{r=inv(cosh((Comp)x));}
	void  q_q_sech					(Quat &r, Quat const &x)				{r=inv(cosh(x));}

	void  c_c_asech					(Quat &r, Quat const &x)				{r=acosh(inv((Comp)x));}
	void  q_q_asech					(Quat &r, Quat const &x)				{r=acosh(inv(x));}

	void  r_r_sin					(Quat &r, Quat const &x)				{r=sin(x.r);}
	void  c_c_sin					(Quat &r, Quat const &x)				{r=sin((Comp)x);}
	void  q_q_sin					(Quat &r, Quat const &x)				{r=sin(x);}

	inline Comp asin(Comp const &x){return -m_i*log(m_i*x+sqrt(1.-sq(x)));}
	inline Quat asin(Quat const &x){return -m_i*log(m_i*x+sqrt(1.-sq(x)));}
	void  c_c_asin					(Quat &r, Quat const &x)				{r=asin((Comp)x);}
	void  q_q_asin					(Quat &r, Quat const &x)				{r=asin(x);}

	inline Comp sinh(Comp const &x){return (exp(x)-exp(-x))*0.5;}
	inline Quat sinh(Quat const &x){return (exp(x)-exp(-x))*0.5;}
	void  r_r_sinh					(Quat &r, Quat const &x)				{r=sinh(x.r);}
	void  c_c_sinh					(Quat &r, Quat const &x)				{r=sinh((Comp)x);}
	void  q_q_sinh					(Quat &r, Quat const &x)				{r=sinh(x);}

	void  r_r_asinh					(Quat &r, Quat const &x)				{r=asinh(x.r);}
	void  c_c_asinh					(Quat &r, Quat const &x)				{r=asinh((Comp)x);}
	void  q_q_asinh					(Quat &r, Quat const &x)				{r=asinh(x);}

	void  r_r_sinc					(Quat &r, Quat const &x)				{r=x.r.toDouble()?sin(x.r)/x.r:1;}
	void  c_c_sinc					(Quat &r, Quat const &x)				{r=x.c_is_true()?sin(x)/x:Comp(1);}
	void  q_q_sinc					(Quat &r, Quat const &x)				{r=x.q_is_true()?sin(x)/x:Quat(1);}

	void  r_r_sinhc					(Quat &r, Quat const &x)				{r=x.r.toDouble()?sinh(x.r)/x.r:1;}
	void  c_c_sinhc					(Quat &r, Quat const &x)				{r=x.c_is_true()?sinh(x)/x:Comp(1);}
	void  q_q_sinhc					(Quat &r, Quat const &x)				{r=x.q_is_true()?sinh(x)/x:Quat(1);}

	void  r_r_csc					(Quat &r, Quat const &x)				{r=1/sin(x.r);}
	void  c_c_csc					(Quat &r, Quat const &x)				{r=inv(sin((Comp)x));}
	void  q_q_csc					(Quat &r, Quat const &x)				{r=inv(sin((Quat)x));}

	void  c_c_acsc					(Quat &r, Quat const &x)				{r=asin(inv((Comp)x));}
	void  q_q_acsc					(Quat &r, Quat const &x)				{r=asin(inv(x));}

	void  r_r_csch					(Quat &r, Quat const &x)				{r=1/sinh(x.r);}
	void  c_c_csch					(Quat &r, Quat const &x)				{r=inv(sinh((Comp)x));}
	void  q_q_csch					(Quat &r, Quat const &x)				{r=inv(sinh(x));}

	void  r_r_acsch					(Quat &r, Quat const &x)				{r=asinh(1/x.r);}
	void  c_c_acsch					(Quat &r, Quat const &x)				{r=asinh(inv((Comp)x));}
	void  q_q_acsch					(Quat &r, Quat const &x)				{r=asinh(inv(x));}

	inline Comp tan(Comp const &x)
	{
		const Comp two_i(0, 2);
		Comp exp_2ix=exp(two_i*x);
		return (exp_2ix-1.)/((exp_2ix+1.)*m_i);
	}
	inline Quat tan(Quat const &x)
	{
		const Comp two_i(0, 2);
		Quat exp_2ix=exp(two_i*x);
		return (exp_2ix-1.)/((exp_2ix+1.)*m_i);
	}
	void  r_r_tan					(Quat &r, Quat const &x)				{r=tan(x.r);}
	void  c_c_tan					(Quat &r, Quat const &x)				{r=tan((Comp)x);}
	void  q_q_tan					(Quat &r, Quat const &x)				{r=tan(x);}

	inline Comp atan(Comp const &x){return (m_i*0.5)*log((m_i+x)/(m_i-x));}
	inline Quat atan(Quat const &x){return (m_i*0.5)*log((m_i+x)/(m_i-x));}
	void  r_r_atan					(Quat &r, Quat const &x)				{r=atan(x.r);}
	void  c_c_atan					(Quat &r, Quat const &x)				{r=atan((Comp)x);}
	void  q_q_atan					(Quat &r, Quat const &x)				{r=atan(x);}
	void r_rr_atan					(Quat &r, Quat const &x, Quat const &y)	{r=atan2(x.r, y.r);}
	void c_rc_atan					(Quat &r, Quat const &x, Quat const &y)	{Comp t=atan((Comp)y/x.r);		r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_rq_atan					(Quat &r, Quat const &x, Quat const &y)	{Quat t=atan(y/x.r);			r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void c_cr_atan					(Quat &r, Quat const &x, Quat const &y)	{Comp t=atan(y.r/(Comp)x);		r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void c_cc_atan					(Quat &r, Quat const &x, Quat const &y)	{Comp t=atan((Comp)y/(Comp)x);	r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_cq_atan					(Quat &r, Quat const &x, Quat const &y)	{Quat t=atan(y/(Comp)x);		r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_qr_atan					(Quat &r, Quat const &x, Quat const &y)	{Quat t=atan(y/x);				r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_qc_atan					(Quat &r, Quat const &x, Quat const &y)	{Quat t=atan((Comp)y/x);		r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_qq_atan					(Quat &r, Quat const &x, Quat const &y)	{Quat t=atan(y/x);				r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	
	inline Comp tanh(Comp const &x){Comp e2x=exp(x+x); return (e2x-1.)/(e2x+1.);}
	inline Quat tanh(Quat const &x){Quat e2x=exp(x+x); return (e2x-1.)/(e2x+1.);}
	void  r_r_tanh					(Quat &r, Quat const &x)				{r=tanh(x.r);}
	void  c_c_tanh					(Quat &r, Quat const &x)				{r=tanh((Comp)x);}
	void  q_q_tanh					(Quat &r, Quat const &x)				{r=tanh(x);}

	inline Comp atanh(Comp const &x){return 0.5*log((1+x)/(1-x));}
	inline Quat atanh(Quat const &x){return 0.5*log((1+x)/(1-x));}
	void  c_c_atanh					(Quat &r, Quat const &x)				{r=atanh((Comp)x);}
	void  q_q_atanh					(Quat &r, Quat const &x)				{r=atanh(x);}

	void  r_r_tanc					(Quat &r, Quat const &x)				{r=x.r.toDouble()?tan(x.r)/x.r:0;}
	void  c_c_tanc					(Quat &r, Quat const &x)				{Comp cx=(Comp)x; r=x.c_is_true()?tan(cx)/cx:Comp();}
	void  q_q_tanc					(Quat &r, Quat const &x)				{r=x.q_is_true()?tan(x)/x:Quat();}

	void  r_r_cot					(Quat &r, Quat const &x)				{r=1/tan(x.r);}
	void  c_c_cot					(Quat &r, Quat const &x)				{r=inv(tan((Comp)x));}
	void  q_q_cot					(Quat &r, Quat const &x)				{r=inv(tan(x));}

	void  r_r_acot					(Quat &r, Quat const &x)				{r=x.r.toDouble()?atan(1/x.r):m_pi/2;}
	void  c_c_acot					(Quat &r, Quat const &x)				{r=atan(inv((Comp)x));}
	void  q_q_acot					(Quat &r, Quat const &x)				{r=atan(inv(x));}

	void  r_r_coth					(Quat &r, Quat const &x)				{r=1/tanh(x.r);}
	void  c_c_coth					(Quat &r, Quat const &x)				{r=inv(tanh((Comp)x));}
	void  q_q_coth					(Quat &r, Quat const &x)				{r=inv(tanh(x));}

	void  c_c_acoth					(Quat &r, Quat const &x)				{r=atanh(inv((Comp)x));}
	void  q_q_acoth					(Quat &r, Quat const &x)				{r=atanh(inv((Quat)x));}

	void  r_r_exp					(Quat &r, Quat const &x)				{r=exp(x.r);}
	void  c_c_exp					(Quat &r, Quat const &x)				{r=exp((Comp)x);}
	void  q_q_exp					(Quat &r, Quat const &x)				{r=exp(x);}
	
	void  r_r_fib					(Quat &r, Quat const &x)				{r=(exp(x.r*m_ln_phi)-cos(m_pi*x.r)*exp(-x.r*m_ln_phi))*m_inv_sqrt5;}
	void  c_c_fib					(Quat &r, Quat const &x)				{Comp cx=(Comp)x; r=(exp(cx*m_ln_phi)-cos(m_pi*cx)*exp(-cx*m_ln_phi))*m_inv_sqrt5;}
	void  q_q_fib					(Quat &r, Quat const &x)				{r=(exp(x*m_ln_phi)-cos(m_pi*x)*exp(-x*m_ln_phi))*m_inv_sqrt5;}
	
	void  r_r_random				(Quat &r, Quat const &x)				{r=mpfr::random();}
	void  c_c_random				(Quat &r, Quat const &x)				{r=Comp(mpfr::random(), mpfr::random());}
	void  q_q_random				(Quat &r, Quat const &x)				{r=Quat(mpfr::random(), mpfr::random(), mpfr::random(), mpfr::random());}
	void r_rr_random				(Quat &r, Quat const &x, Quat const &y)	{r=mpfr::random();}
	void c_cr_random				(Quat &r, Quat const &x, Quat const &y)	{r=Comp(mpfr::random(), mpfr::random());}
	void c_cc_random				(Quat &r, Quat const &x, Quat const &y)	{r=Comp(mpfr::random(), mpfr::random());}
	void q_qq_random				(Quat &r, Quat const &x, Quat const &y)	{r=Quat(mpfr::random(), mpfr::random(), mpfr::random(), mpfr::random());}

	inline Real beta(Real x, Real y){return exp(lgamma(x)+lgamma(y)-lgamma(x+y));}
	void  r_r_beta					(Quat &r, Quat const &x)				{r=beta(x.r, x.r);}
	void r_rr_beta					(Quat &r, Quat const &x, Quat const &y)	{r=beta(x.r, y.r);}
	
	void  r_r_cyl_bessel_j			(Quat &r, Quat const &x)				{r=besselj0(x.r);}
	void r_rr_cyl_bessel_j			(Quat &r, Quat const &x, Quat const &y)	{r=besseljn(x.r.toLong(), y.r);}

	void  r_r_cyl_neumann			(Quat &r, Quat const &x)				{r=bessely0(x.r);}
	void r_rr_cyl_neumann			(Quat &r, Quat const &x, Quat const &y)	{r=besselyn(x.r.toLong(), y.r);}

	inline Comp r_hankel1(Real x, Real y)
	{
		return besseljn(x.toLong(), y)+m_i*besselyn(x.toLong(), y);
		//try
		//{
		//	return boost::math::cyl_bessel_j(x, y)+Comp(0, 1)*boost::math::cyl_neumann(x, y);
		//}
		//catch(std::domain_error&){return _qnan;}
		//catch(std::overflow_error&){return -_HUGE;}
		//catch(...){return _qnan;}
	}
	void  c_r_hankel1				(Quat &r, Quat const &x)				{r=r_hankel1(x.r, x.r);}
	void  c_c_hankel1				(Quat &r, Quat const &x)				{r=r_hankel1(x.r, x.i);}
	void c_rr_hankel1				(Quat &r, Quat const &x, Quat const &y)	{r=r_hankel1(x.r, y.r);}
	
//	inline Real	sgn			(Real const &x){return (x>0)-(x<0);}
	inline Comp	sgn			(Comp const &x)
	{
		Real temp=abs(x);
		return temp.toDouble()?x/temp:Comp();
	}
	inline Quat	sgn			(Quat	const &x)
	{
		Real temp=abs(x);
		return temp!=0?x/temp:Quat();
	}
	inline Real	step		(Real const &x){return 0.5+0.5*sgn(x);}
	inline Comp	step		(Comp const &x){return 0.5+0.5*sgn(x);}
	inline Quat	step		(Quat const &x){return 0.5+0.5*sgn(x);}
	void  r_r_step					(Quat &r, Quat const &x)				{r=step(x.r);}
	void  c_c_step					(Quat &r, Quat const &x)				{r=step((Comp)x);}
	void  q_q_step					(Quat &r, Quat const &x)				{r=step(x);}

	void  r_r_rect					(Quat &r, Quat const &x)				{r=step(x.r+0.5)-step(x.r-0.5);}
	void  c_c_rect					(Quat &r, Quat const &x)				{Comp cx=(Comp)x; r=step(cx+0.5)-step(cx-0.5);}
	void  q_q_rect					(Quat &r, Quat const &x)				{Quat qx=x; r=step(qx+0.5)-step(qx-0.5);}

	void  r_r_trgl					(Quat &r, Quat const &x)				{Real t=abs(x.r);		r=t<1?1-t:0;}
	void  r_c_trgl					(Quat &r, Quat const &x)				{Real t=abs((Comp)x);	r=t<1?1-t:0;}
	void  r_q_trgl					(Quat &r, Quat const &x)				{Real t=abs(x);	r=t<1?1-t:0;}

	void  r_r_sqwv					(Quat &r, Quat const &x)				{r=x.r-floor(x.r)<0.5;}
	void  r_c_sqwv					(Quat &r, Quat const &x)				{r=x.r-floor(x.r)<0.5;}
	void  r_q_sqwv					(Quat &r, Quat const &x)				{r=x.r-floor(x.r)<0.5;}
	void r_rr_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_rc_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_rq_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_cr_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_cc_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_cq_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_qr_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_qc_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_qq_sqwv					(Quat &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}

	Real clamp_positive(Real const &x){return signbit(x)?x:0;}
	void  r_r_trwv					(Quat &r, Quat const &x)				{auto temp=abs(x.r-floor(x.r)-0.5); r=temp+temp;}
	void  r_c_trwv					(Quat &r, Quat const &x)				{Comp cx=(Comp)x; auto temp=abs(cx-floor(cx)-0.5); r=temp+temp;}
	void  r_q_trwv					(Quat &r, Quat const &x)				{Quat qx=x; r=2*abs(qx-floor(qx)-0.5);}
	void r_rr_trwv					(Quat &r, Quat const &x, Quat const &y)
	{
		Real t=x.r-floor(x.r), t2=1-x.r;
		t2-=floor(t2);
		Real dc=clamp_positive(y.r);
		dc=dc>1?1:dc;
		Real dc2=1-dc, t_d=t/dc, t2_d2=t2/dc2;
		r=(t_d<1)*t_d+(t2_d2<1)*t2_d2;
	}
	void c_cr_trwv					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cx=(Comp)x, t=cx-floor(cx), t2=1.-cx;
		t2-=floor(t2);
		Real dc=clamp_positive(y.r);
		dc=dc>1?1:dc;
		Real dc2=1-dc;
		auto t_d=t/dc, t2_d2=t2/dc2;
		r=Real(t_d.r<1)*t_d+Real(t2_d2.r<1)*t2_d2;
	}
	void c_cc_trwv					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cx=(Comp)x, cy=y, t=cx-floor(cx), t2=1.-cx;
		t2-=floor(t2);
		Comp dc=cy;
		dc.r=clamp_positive(dc.r);
		dc.r=dc.r>1?1:dc.r;
		Comp dc2=1.-dc;
		auto t_d=t/dc, t2_d2=t2/dc2;
		r=Real(t_d.r<1)*t_d+Real(t2_d2.r<1)*t2_d2;
	}
	void q_qq_trwv					(Quat &r, Quat const &x, Quat const &y)
	{
		Quat cx=x, cy=y, t=cx-floor(cx), t2=1.-cx;
		t2-=floor(t2);
		Quat dc=cy;
		dc.r=clamp_positive(dc.r);
		dc.r=dc.r>1?1:dc.r;
		Quat dc2=1.-dc;
		auto t_d=t/dc, t2_d2=t2/dc2;
		r=Real(t_d.r<1)*t_d+Real(t2_d2.r<1)*t2_d2;
	}

	void  r_r_saw					(Quat &r, Quat const &x)				{Real t=x.r-floor(x.r), t2=floor(1-t); r=(t2+1)*(t2*0.5+t);}
	void  c_c_saw					(Quat &r, Quat const &x)				{Comp cx=(Comp)x, t=cx-floor(cx), t2=floor(1-t); r=(t2+1)*(t2*0.5+t);}
	void  q_q_saw					(Quat &r, Quat const &x)				{Quat t=x-floor(x), t2=floor(1-t); r=(t2+1)*(t2*0.5+t);}
	void r_rr_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		auto t=x.r-floor(x.r), t2=floor(y.r-t);
		r=(t2+1)*(t2*0.5+t)/y;
	}
	void c_rc_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cy=(Comp)y;
		auto t=x.r-floor(x.r);
		auto t2=floor(cy-t);
		r=(t2+1.)*(t2*0.5+t)/cy;
	}
	void q_rq_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		auto t=x.r-floor(x.r);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void c_cr_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cx=(Comp)x;
		auto t=cx-floor(cx);
		auto t2=floor(y.r-t);
		r=(t2+1.)*(t2*0.5+t)/y.r;
	}
	void c_cc_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cx=(Comp)x, cy=(Comp)y;
		auto t=cx-floor(cx);
		auto t2=floor(cy-t);
		r=(t2+1.)*(t2*0.5+t)/cy;
	}
	void q_cq_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cx=(Comp)x;
		auto t=cx-floor(cx);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void q_qr_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y.r-t);
		r=(t2+1.)*(t2*0.5+t)/y.r;
	}
	void q_qc_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		Comp cy=(Comp)y;
		auto t=x-floor(x);
		auto t2=floor(cy-t);
		r=(t2+1.)*(t2*0.5+t)/cy;
	}
	void q_qq_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}

	void r_rr_hypot					(Quat &r, Quat const &x, Quat const &y)	{r=sqrt(x.r*x.r+y.r*y.r);}
	//void c_cc_hypot					(Comp const &x, Comp const &y)	{return sqrt(sq(x)+sq(y));}
	//void q_qq_hypot					(Quat const &x, Quat const &y)	{return sqrt(sq(x)+sq(y));}

	inline int mandelbrot(Comp const &x, int n_iterations)
	{
		Real rez=0, imz=0, sq_rez=0, sq_imz=0;
		int k=0;
		for(;k<n_iterations&&sq_rez+sq_imz<16;++k)
		{
			imz=rez*imz;//calculate sq(z)
			imz+=imz;
			rez=sq_rez-sq_imz;

			rez+=x.r, imz+=x.i;//add x

			sq_rez=rez*rez, sq_imz=imz*imz;
		}
		return k;
	}
	void r_r_mandelbrot				(Quat &r, Quat const &x)				{r=mandelbrot(Comp(x.r), 200);}
	void r_c_mandelbrot				(Quat &r, Quat const &x)				{r=mandelbrot((Comp)x, 200);}
	void r_rr_mandelbrot			(Quat &r, Quat const &x, Quat const &y)	{r=mandelbrot(Comp(x.r), y.r.toLong());}
	void r_cr_mandelbrot			(Quat &r, Quat const &x, Quat const &y)	{r=mandelbrot((Comp)x, y.r.toLong());}

	void r_rr_min					(Quat &r, Quat const &x, Quat const &y)	{r=(x.r+y.r-abs(x.r-y.r))*0.5;}
	void c_cr_min					(Quat &r, Quat const &x, Quat const &y)	{Comp cx=(Comp)x;				r=(cx+y.r-abs(cx-y.r))*0.5;}
	void c_cc_min					(Quat &r, Quat const &x, Quat const &y)	{Comp cx=(Comp)x, cy=(Comp)y;	r=(cx+cy-abs(cx-cy))*0.5;}
	void q_qq_min					(Quat &r, Quat const &x, Quat const &y)	{r=(x+y-abs(x-y))*0.5;}

	void r_rr_max					(Quat &r, Quat const &x, Quat const &y)	{r=(x.r+y.r+abs(x.r-y.r))*0.5;}
	void c_cr_max					(Quat &r, Quat const &x, Quat const &y)	{Comp cx=(Comp)x;				r=(cx+y.r+abs(cx-y.r))*0.5;}
	void c_cc_max					(Quat &r, Quat const &x, Quat const &y)	{Comp cx=(Comp)x, cy=(Comp)y;	r=(cx+cy+abs(cx-cy))*0.5;}
	void q_qq_max					(Quat &r, Quat const &x, Quat const &y)	{r=(x+y+abs(x-y))*0.5;}

	void r_rr_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.r.toDouble()?y.r:0;}
	void c_rc_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.r.toDouble()?(Comp)y:Comp();}
	void q_rq_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.r.toDouble()?y:Quat();}
	void r_cr_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?y.r:0;}
	void c_cc_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?(Comp)y:Comp();}
	void q_cq_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?y:Quat();}
	void r_qr_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?y.r:0;}
	void c_qc_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?(Comp)y:Comp();}
	void q_qq_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?y:Quat();}
	
	void r_rr_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.r.toDouble()?0:y.r;}
	void c_rc_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.r.toDouble()?Comp():(Comp)y;}
	void q_rq_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.r.toDouble()?Quat():y;}
	void r_cr_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?0:y.r;}
	void c_cc_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?Comp():(Comp)y;}
	void q_cq_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.c_is_true()?Quat():y;}
	void r_qr_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?0:y.r;}
	void c_qc_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?Comp():(Comp)y;}
	void q_qq_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?Quat():y;}

	void  r_r_increment				(Quat &r, Quat const &x)				{r=x.r+1;}
	void  c_c_increment				(Quat &r, Quat const &x)				{r=(Comp)x+1;}
	void  q_q_increment				(Quat &r, Quat const &x)				{r=x+1;}

	void  r_r_decrement				(Quat &r, Quat const &x)				{r=x.r-1;}
	void  c_c_decrement				(Quat &r, Quat const &x)				{r=(Comp)x-1;}
	void  q_q_decrement				(Quat &r, Quat const &x)				{r=x-1;}

	void  r_r_assign				(Quat &r, Quat const &x)				{r=x.r;}
	void  c_c_assign				(Quat &r, Quat const &x)				{r=(Comp)x;}
	void  q_q_assign				(Quat &r, Quat const &x)				{r=x;}
#else
	void  r_r_setzero				(Real &r, Real const&)					{r.setZero();}
	void  c_c_setzero				(Comp &r, Comp const&)					{r.r.setZero(), r.i.setZero();}
	void  q_q_setzero				(Quat &r, Quat const&)					{r.r.setZero(), r.i.setZero(), r.j.setZero(), r.k.setZero();}
	
	void  r_r_ceil					(Real &r, Real const &x)				{r=ceil(x);}
	void  c_c_ceil					(Comp &r, Comp const &x)				{r=ceil(x);}
	void  q_q_ceil					(Quat &r, Quat const &x)				{r=ceil(x);}

	void  r_r_floor					(Real &r, Real const &x)				{r=floor(x);}
	void  c_c_floor					(Comp &r, Comp const &x)				{r=floor(x);}
	void  q_q_floor					(Quat &r, Quat const &x)				{r=floor(x);}

	void  r_r_round					(Real &r, Real const &x)				{r=round(x);}
	void  c_c_round					(Comp &r, Comp const &x)				{r=round(x);}
	void  q_q_round					(Quat &r, Quat const &x)				{r=round(x);}

	void  r_r_int					(Real &r, Real const &x)				{r=trunc(x);}
	void  c_c_int					(Comp &r, Comp const &x)
	{
		r.r=trunc(x.r), r.i=trunc(x.i);
		//mpfr_trunc(r.r.mpfr_ptr(), x.r.mpfr_ptr());
		//mpfr_trunc(r.i.mpfr_ptr(), x.i.mpfr_ptr());
	}
	void  q_q_int					(Quat &r, Quat const &x)
	{
		r.r=trunc(x.r), r.i=trunc(x.i), r.j=trunc(x.j), r.k=trunc(x.k);
		//mpfr_trunc(r.r.mpfr_ptr(), x.r.mpfr_ptr());
		//mpfr_trunc(r.i.mpfr_ptr(), x.i.mpfr_ptr());
		//mpfr_trunc(r.j.mpfr_ptr(), x.j.mpfr_ptr());
		//mpfr_trunc(r.k.mpfr_ptr(), x.k.mpfr_ptr());
	}

	void  r_r_frac					(Real &r, Real const &x)				{r=frac(x);}
	void  c_c_frac					(Comp &r, Comp const &x)				{r.r=frac(x.r), r.i=frac(x.i);}
	void  q_q_frac					(Quat &r, Quat const &x)				{r.r=frac(x.r), r.i=frac(x.i), r.j=frac(x.j), r.k=frac(x.k);}

	void  r_r_abs					(Real &r, Real const &x)				{r=abs(x);}
	void  r_c_abs					(Real &r, Comp const &x)				{r=abs(x);}
	void  r_q_abs					(Real &r, Quat const &x)				{r=abs(x);}

	void  r_r_arg					(Real &r, Real const &x)				{r=x<0?m_pi:x==0?m_qnan:0;}
	void  r_c_arg					(Real &r, Comp const &x)				{r=(x.r==0)&(x.i==0)?m_qnan:atan2(x.i, x.r);}
	void  r_q_arg					(Real &r, Quat const &x)				{r=acos(x.r/abs(x));}

	void  r_c_real					(Real &r, Comp const &x)				{r=x.r;}

	void  r_c_imag					(Real &r, Comp const &x)				{r=x.i;}

	//r_conjugate: assign
	void c_c_conjugate				(Comp &r, Comp const &x)				{r=Comp(x.r, -x.i);}
	void q_q_conjugate				(Quat &r, Quat const &x)				{r=Quat(x.r, -x.i, -x.j, -x.k);}

	void  c_r_polar					(Comp &r, Real const &x)				{r=Comp(abs(x), x<0?m_pi:x==0?m_qnan:0);}
	void  c_c_polar					(Comp &r, Comp const &x)
	{
		Real mag=abs(x);
		r=Comp(mag, mag==0?m_qnan:atan2(x.i, x.r));
	}
	void  c_q_polar					(Comp &r, Quat const &x)
	{
		Real mag=abs(x);
		r=Comp(mag, acos(x.r/mag));
	}

	//r_cartesian	assign
	void  c_c_cartesian				(Comp &r, Comp const &x)				{r=Comp(x.r*cos(x.i), x.r*sin(x.i));}
	void  q_q_cartesian				(Quat &r, Quat const &x)
	{
		Real cos_j=cos(x.j), xr_cos_k=x.r*cos(x.k);
		r=Quat(
			cos(x.i)*cos_j*xr_cos_k,
			sin(x.i)*cos_j*xr_cos_k,
			sin(x.j)*xr_cos_k,
			x.r*sin(x.k));
	}

	void r_rr_plus					(Real &r, Real const &x, Real const &y)	{r=x+y;}
	void c_rc_plus					(Comp &r, Real const &x, Comp const &y)	{r=x+y;}
	void q_rq_plus					(Quat &r, Real const &x, Quat const &y)	{r=x+y;}
	void c_cr_plus					(Comp &r, Comp const &x, Real const &y)	{r=x+y;}
	void c_cc_plus					(Comp &r, Comp const &x, Comp const &y)	{r=x+y;}
	void q_cq_plus					(Quat &r, Comp const &x, Quat const &y)	{r=x+y;}
	void q_qr_plus					(Quat &r, Quat const &x, Real const &y)	{r=x+y;}
	void q_qc_plus					(Quat &r, Quat const &x, Comp const &y)	{r=x+y;}
	void q_qq_plus					(Quat &r, Quat const &x, Quat const &y)	{r=x+y;}

	void  r_r_minus					(Real &r, Real const &x)				{r=-x;}
	void  c_c_minus					(Comp &r, Comp const &x)				{r=-x;}
	void  q_q_minus					(Quat &r, Quat const &x)				{r=-x;}
	void r_rr_minus					(Real &r, Real const &x, Real const &y)	{r=x-y;}
	void c_rc_minus					(Comp &r, Real const &x, Comp const &y)	{r=x-y;}
	void q_rq_minus					(Quat &r, Real const &x, Quat const &y)	{r=x-y;}
	void c_cr_minus					(Comp &r, Comp const &x, Real const &y)	{r=x-y;}
	void c_cc_minus					(Comp &r, Comp const &x, Comp const &y)	{r=x-y;}
	void q_cq_minus					(Quat &r, Comp const &x, Quat const &y)	{r=x-y;}
	void q_qr_minus					(Quat &r, Quat const &x, Real const &y)	{r=x-y;}
	void q_qc_minus					(Quat &r, Quat const &x, Comp const &y)	{r=x-y;}
	void q_qq_minus					(Quat &r, Quat const &x, Quat const &y)	{r=x-y;}

	void r_rr_multiply				(Real &r, Real const &x, Real const &y)	{r=x*y;}
	void c_rc_multiply				(Comp &r, Real const &x, Comp const &y)	{r=x*y;}
	void q_rq_multiply				(Quat &r, Real const &x, Quat const &y)	{r=x*y;}
	void c_cr_multiply				(Comp &r, Comp const &x, Real const &y)	{r=x*y;}
	void c_cc_multiply				(Comp &r, Comp const &x, Comp const &y)	{r=x*y;}
	void q_cq_multiply				(Quat &r, Comp const &x, Quat const &y)	{r=x*y;}
	void q_qr_multiply				(Quat &r, Quat const &x, Real const &y)	{r=x*y;}
	void q_qc_multiply				(Quat &r, Quat const &x, Comp const &y)	{r=x*y;}
	void q_qq_multiply				(Quat &r, Quat const &x, Quat const &y)	{r=x*y;}
	
	inline Comp inv(Comp const &x)
	{
		Real inv_mag=1/abs(x);
		return Comp(x.r*inv_mag, -x.i*inv_mag);
	}
	inline Quat inv(Quat const &x)
	{
		Real inv_mag=1/abs(x);
		return Quat(x.r*inv_mag, -x.i*inv_mag, -x.j*inv_mag, -x.k*inv_mag);
	}
	void  r_r_divide				(Real &r, Real const &y)				{r=1/y;}
	void  c_c_divide				(Comp &r, Comp const &y)				{r=inv(y);}
	void  q_q_divide				(Quat &r, Quat const &y)				{r=inv(y);}
	void r_rr_divide				(Real &r, Real const &x, Real const &y)	{r=x/y;}
	void c_rc_divide				(Comp &r, Real const &x, Comp const &y)	{r=x/y;}
	void q_rq_divide				(Quat &r, Real const &x, Quat const &y)	{r=x/y;}
	void c_cr_divide				(Comp &r, Comp const &x, Real const &y)	{r=x/y;}
	void c_cc_divide				(Comp &r, Comp const &x, Comp const &y)	{r=x/y;}
	void q_cq_divide				(Quat &r, Comp const &x, Quat const &y)	{r=x/y;}
	void q_qr_divide				(Quat &r, Quat const &x, Real const &y)	{r=x/y;}
	void q_qc_divide				(Quat &r, Quat const &x, Comp const &y)	{r=x/y;}
	void q_qq_divide				(Quat &r, Quat const &x, Quat const &y)	{r=x/y;}
	
	void r_rr_logic_divides			(Real &r, Real const &y, Real const &x)	{auto t=x/y; r=t==floor(t);}//rc_divides, rq_divides: applied to each component
	void r_rc_logic_divides			(Real &r, Real const &y, Comp const &x)	{auto t=x/y; r=t==floor(t);}
	void r_rq_logic_divides			(Real &r, Real const &y, Quat const &x)	{auto t=x/y; r=t==floor(t);}
	void r_cr_logic_divides			(Real &r, Comp const &y, Real const &x)	{auto t=x/y; r=t==floor(t);}
	void r_cc_logic_divides			(Real &r, Comp const &y, Comp const &x)	{auto t=x/y; r=t==floor(t);}
	void r_cq_logic_divides			(Real &r, Comp const &y, Quat const &x)	{auto t=x/y; r=t==floor(t);}
	void r_qr_logic_divides			(Real &r, Quat const &y, Real const &x)	{auto t=x/y; r=t==floor(t);}
	void r_qc_logic_divides			(Real &r, Quat const &y, Comp const &x)	{auto t=x/y; r=t==floor(t);}
	void r_qq_logic_divides			(Real &r, Quat const &y, Quat const &x)	{auto t=x/y; r=t==floor(t);}

	void r_rr_power_real			(Real &r, Real const &x, Real const &y)	{r=pow(x, floor(y));}//trunc
	void c_cr_power_real			(Comp &r, Comp const &x, Real const &y)	{r=x^floor(y);}
	void q_qr_power_real			(Quat &r, Quat const &x, Real const &y)	{r=x^floor(y);}

	void c_cr_pow					(Comp &r, Comp const &x, Real const &y)
	{
		//if(!x.r&&!x.i&&!y)//check operator^
		//	r=1;
		//else
			r=x^y;
	}
	void c_cc_pow					(Comp &r, Comp const &x, Comp const &y)
	{
		//if(!x.r&&!x.i&&!y.r&&!y.i)
		//	r=1;
		//else
			r=x^y;
	}
	void q_cq_pow					(Quat &r, Comp const &x, Quat const &y)
	{
		//if(!x.r&&!x.i&&!y.r&&!y.i&&!y.j&&!y.k)
		//	r=Quat(1, 0, 0, 0);
		//else
			r=x^y;
	}
	void q_qr_pow					(Quat &r, Quat const &x, Real const &y)
	{
		//if(!x.r&&!x.i&&!x.j&&!x.k&&!y)
		//	r=Quat(1, 0, 0, 0);
		//else
			r=x^y;
	}
	void q_qc_pow					(Quat &r, Quat const &x, Comp const &y)
	{
		//if(!x.r&&!x.i&&!x.j&&!x.k&&!y.r&&!y.i)
		//	r=Quat(1, 0, 0, 0);
		//else
			r=x^y;
	}
	void q_qq_pow					(Quat &r, Quat const &x, Quat const &y)
	{
		//if(!x.r&&!x.i&&!x.j&&!x.k&&!y.r&&!y.i&&!y.j&&!y.k)
		//	r=Quat(1, 0, 0, 0);
		//else
			r=x^y;
	}

	void  c_c_ln					(Comp &r, Comp const &x)				{r=log(x);}
	void  q_q_ln					(Quat &r, Quat const &x)				{r=log(x);}
	
	void  c_c_log					(Comp &r, Comp const &x)				{r=log(x)*m_inv_ln10;}
	void  q_q_log					(Quat &r, Quat const &x)				{r=log(x)*m_inv_ln10;}
	void c_cr_log					(Comp &r, Comp const &x, Real const &y)	{r=log(x)/log(Comp(y, 0));}
	void c_cc_log					(Comp &r, Comp const &x, Comp const &y)	{r=log(x)/log(y);}
	void q_cq_log					(Quat &r, Comp const &x, Quat const &y)	{r=log(x)/log(y);}
	void q_qc_log					(Quat &r, Quat const &x, Comp const &y)	{r=log(x)/log(y);}
	void q_qq_log					(Quat &r, Quat const &x, Quat const &y)	{r=log(x)/log(y);}
	
	Comp tetrate(Real const &x, Real const &y)
	{
		Real log_x=log(x), ry=y;
		int h=int(floor(ry-.5)+1);//rounded real part
		Comp t(ry-h, 0);//real-round(real)+i*imag
		{
			auto q=sqrt(log_x);
			t=(t+1.)*(
				(
					log_x>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
						+t*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
								+t*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
										+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
						+t*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+m_ln2)/log_x)
								+t*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
										+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
									)
							)
				)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	Comp tetrate(Comp const &x, Real const &y)
	{
		Comp log_x=log(x);
		Real ry=y;
		int h=int(floor(ry-.5)+1);//rounded real part
		Comp t(ry-h, 0);//real-round(real)+i*imag
		{
			auto q=sqrt(log_x);
			t=(t+1.)*(
				(
					log_x.r>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
						+t*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
								+t*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
										+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
						+t*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+m_ln2)/log_x)
								+t*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
										+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
									)
							)
				)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	Quat tetrate(Quat const &x, Real const &y)
	{
		Quat qx=x;
		Real ry=y;
		if(ry<-1)
		{
			int steps=int(abs(ry));
			Quat t(ry-floor(ry), 0, 0, 0), lrx=log(qx);
			for(int k=0;k<steps;++k)
				t=log(t)/lrx;
			return t;
		}
		else if(ry<=0)
			return Quat(1+ry, 0, 0, 0);
		else
		{
			int h=int(ry)+1;
			Quat t(ry-floor(ry), 0, 0, 0);
			for(int k=0;k<h;++k)
				t=qx^t;
			return t;
		}
	}
	Comp tetrate(Real const &x, Comp const &y)
	{
		Real log_x=log((Real)x);
		Comp cy=y;
		int h=int(floor(cy.r-.5)+1);//rounded real part
		Comp t(cy-(Real)h);//real-round(real)+i*imag
		{
			auto q=sqrt(log_x);
			t=(t+1.)*(
				(
					log_x>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
						+t*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
								+t*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
										+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
						+t*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*q)*q-(-0.625+m_ln2)/log_x)
								+t*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
										+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
									)
							)
				)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	Comp tetrate(Comp const &x, Comp const &y)
	{
		Comp log_x=log(x), cy=y;
	//	if(log_x.r<.03)//abs(log_x)<1.03045453395352
	//	{
	//		if(cy.r<-1.)
	//			return -30.;
	//		return 1.;
	//	}
		int h=int(floor(cy.r-.5)+1);
		Comp t(cy-(Real)h);//real-round(real)+i*imag
		{
		//	bool unassigned=true;
		//	if(log_x.r<.001)//abs(log_x)<1.00100050016671
		//	{
		//		if(t.r>-1)//real-round(real)>-1
		//			unassigned=false, t=1.;
		//		else if(t.r<-1)
		//			unassigned=false, t=-990.;
		//	}
		//	if(unassigned)
		//	{
				Comp q=sqrt(log_x);
				t=(t+1.)*(
					(
						log_x.r>m_ln2?
							(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/log_x-1.)
							+t*	(
									(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/log_x-1.)
									+t*	(
											(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/log_x-1.)
											+t*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+log_x*61.6566))))/(1.-8.1192*q+37.087*log_x)+(-131./192.+m_ln2)/log_x-1.)
										)
								)
						:
							(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*log_x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/log_x)
							+t*	(
									(1.1-2.608785958462561*(1.-0.6663562294911147*sqrt(log_x))*sqrt(log_x)-(-0.625+m_ln2)/log_x)
									+t*	(
											(-0.96+3.0912038297987596*(1.+0.6021398048785328*log_x)*q/(1.+ 4.240467556480155*log_x)+(-2./3+m_ln2)/log_x)
											+t*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*log_x)*q/(1.+4.95715636660691*q + 7.70233216637738*log_x)-(-131./192.+m_ln2)/log_x)
										)
								)
					)*t+1.)+log(t+2.)/log_x-m_ln2/log_x*(1.+t);
		//	}
		}
		for(;h>0;--h)
			t=exp(log_x*t);
		for(;h<0;++h)
			t=log(t)/log_x;
		return t;
	}
	void c_rr_tetrate				(Comp &r, Real const &x, Real const &y)	{r=tetrate(x, y);}
	void c_rc_tetrate				(Comp &r, Real const &x, Comp const &y)	{r=tetrate(x, y);}
	void c_cr_tetrate				(Comp &r, Comp const &x, Real const &y)	{r=tetrate((Comp)x, (Real)y);}
	void c_cc_tetrate				(Comp &r, Comp const &x, Comp const &y)	{r=tetrate(x, y);}
	void q_qr_tetrate				(Quat &r, Quat const &x, Real const &y)	{r=tetrate((Quat)x, (Real)y);}
	
	class Tetrate//http://en.citizendium.org/wiki/Fit1.cin
	{
		static Comp _fit1(Comp &x, Comp &y)
		{
			if(x.r<.001)
			{
				if(y.r>-1)
					return Comp(1, 0);
				if(y.r<-1)
					return Comp(-990, 0);
			}
			Comp q=sqrt(x);
			return (y+1.)*(
				(
					x.r>m_ln2?
						(q*(0.137467+q*(4.94969+q*0.0474179))/(1.+q*(3.23171+q*0.471222))+(-.5+m_ln2)/x-1.)
						+y*	(
								(q*(-0.424278+q*(1.75166+q*(-1.46524+q*0.93347)))/(0.0312142+q*(-0.267478+q))+(-.625+m_ln2)/x-1.)
								+y*	(
										(q*(3.39255+q*(16.1046+q*(-19.5216+q*10.7458)))/(1.+q*(4.1274+q*5.25449))+(-2./3+m_ln2)/x-1.)
										+y*(0.16*q*(1.+q*(27.7934+q*(358.688+q*(-259.233+x*61.6566))))/(1.-8.1192*q+37.087*x)+(-131./192.+m_ln2)/x-1.)
									)
							)
					:
						(-1.0018+(0.15128484821526975*(1.+33.04715298851381*q-3.51771875598067*x)*q)/(1.+3.2255053261256337*q)+(-0.5+m_ln2)/x)
						+y*	(
								(1.1-2.608785958462561*(1.-0.6663562294911147*sqrt(x))*sqrt(x)-(-0.625+m_ln2)/x)
								+y*	(
										(-0.96+3.0912038297987596*(1.+0.6021398048785328*x)*q/(1.+ 4.240467556480155*x)+(-2./3+m_ln2)/x)
										+y*(1.2-10.44604984418533*(1.+0.2137568928431227*q+0.3693275254470449*x)*q/(1.+4.95715636660691*q + 7.70233216637738*x)-(-131./192.+m_ln2)/x)
									)
							)
				)*y+1.)+log(y+2.)/x-m_ln2/x*(1.+y);
		}
	public:
		static Comp fit1(Comp x, Comp const &y)//[sic]
		{
			x=log(x);
			if(x.r<.03)
			{
				if(y.r<-1.)
					return Comp(-30, 0);
				return Comp(1, 0);
			}
			int h=int(floor(y.r-.5)+1);
			Comp result=_fit1(x, y-(Real)h);
			for(;h>0;--h)
				result=exp(x*result);
			for(;h<0;++h)
				result=log(result)/x;
			return result;
		}
	};
	Comp pentate(Real const &x, Real const &y)
	{
		long long h=y.toLLong();
	//	long long h=y.r!=y.r||y.r<-ll_max||y.r>ll_max?0:long long(y);
	//	long long h=std::isnan(y.r)||std::isinf(y.r)?0:long long(y);
		if(h<-2)	return Comp(_HUGE, 0);//1/::sin(0);
		if(h==-2)	return Comp(-1, 0);
		if(h==-1)	return Comp(0, 0);
		if(h==0)	return Comp(1, 0);
		if(h==1)	return Comp(x, 0);
		Real rx=x;
		Comp result(rx, 0);
		for(int k=0;k<h;++k)
			result=Tetrate::fit1(Comp(rx, 0), result);
		return result;
	}
	Comp pentate(Comp const &x, Real const &y)
	{
		long long h=y.toLLong();
		if(h<-2)	return Comp(_HUGE, 0);//1/::sin(0);
		if(h==-2)	return Comp(-1, 0);
		if(h==-1)	return Comp(0, 0);
		if(h==0)	return Comp(1, 0);
		if(h==1)	return x;
		Comp result=x;
		for(int k=0;k<h;++k)
			result=Tetrate::fit1(x, result);
		return result;
	}
	void c_rr_pentate				(Comp &r, Real const &x, Real const &y)	{r=pentate(x, y);}
	void c_cr_pentate				(Comp &r, Comp const &x, Real const &y)	{r=pentate(x, y);}

	void  r_r_bitwise_shift_left_l	(Real &r, Real const &x)				{r=exp(floor(x)*m_ln2);}	//<<x = 2^floor(x)
	void  c_c_bitwise_shift_left_l	(Comp &r, Comp const &x)				{r=exp(floor(x)*m_ln2);}
	void  q_q_bitwise_shift_left_l	(Quat &r, Quat const &x)				{r=exp(floor(x)*m_ln2);}
	void  r_r_bitwise_shift_left_r	(Real &r, Real const &x)				{r=x+x;}					//x<< = 2x
	void  c_c_bitwise_shift_left_r	(Comp &r, Comp const &x)				{Comp cx=x; r=cx+cx;}
	void  q_q_bitwise_shift_left_r	(Quat &r, Quat const &x)				{Quat qx=x; r=qx+qx;}
	void r_rr_bitwise_shift_left	(Real &r, Real const &x, Real const &y)	{r=x*exp(floor(y)*m_ln2);}//x<<y = x*2^y
	void c_rc_bitwise_shift_left	(Comp &r, Real const &x, Comp const &y)	{r=x*exp(floor(y)*m_ln2);}
	void q_rq_bitwise_shift_left	(Quat &r, Real const &x, Quat const &y)	{r=x*exp(floor(y)*m_ln2);}
	void c_cr_bitwise_shift_left	(Comp &r, Comp const &x, Real const &y)	{r=(Comp)x*exp(floor(y)*m_ln2);}
	void c_cc_bitwise_shift_left	(Comp &r, Comp const &x, Comp const &y)	{r=(Comp)x*exp(floor(y)*m_ln2);}
	void q_cq_bitwise_shift_left	(Quat &r, Comp const &x, Quat const &y)	{r=(Comp)x*exp(floor(y)*m_ln2);}
	void q_qr_bitwise_shift_left	(Quat &r, Quat const &x, Real const &y)	{r=(Quat)x*exp(floor(y)*m_ln2);}
	void q_qc_bitwise_shift_left	(Quat &r, Quat const &x, Comp const &y)	{r=(Quat)x*exp(floor(y)*m_ln2);}
	void q_qq_bitwise_shift_left	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(floor(y)*m_ln2);}

	void  r_r_bitwise_shift_right_l	(Real &r, Real const &x)				{r=exp(-floor(x)*m_ln2);}
	void  c_c_bitwise_shift_right_l	(Comp &r, Comp const &x)				{r=exp(-floor(x)*m_ln2);}
	void  q_q_bitwise_shift_right_l	(Quat &r, Quat const &x)				{r=exp(-floor(x)*m_ln2);}
	void  r_r_bitwise_shift_right_r	(Real &r, Real const &x)				{r=x*0.5;}
	void  c_c_bitwise_shift_right_r	(Comp &r, Comp const &x)				{r=(Comp)x*0.5;}
	void  q_q_bitwise_shift_right_r	(Quat &r, Quat const &x)				{r=(Quat)x*0.5;}
	void r_rr_bitwise_shift_right	(Real &r, Real const &x, Real const &y)	{r=x*exp(-floor(y)*m_ln2);}
	void c_rc_bitwise_shift_right	(Comp &r, Real const &x, Comp const &y)	{r=x*exp(-floor(y)*m_ln2);}
	void q_rq_bitwise_shift_right	(Quat &r, Real const &x, Quat const &y)	{r=x*exp(-floor(y)*m_ln2);}
	void c_cr_bitwise_shift_right	(Comp &r, Comp const &x, Real const &y)	{r=(Comp)x*exp(-floor(y)*m_ln2);}
	void c_cc_bitwise_shift_right	(Comp &r, Comp const &x, Comp const &y)	{r=(Comp)x*exp(-floor(y)*m_ln2);}
	void q_cq_bitwise_shift_right	(Quat &r, Comp const &x, Quat const &y)	{r=(Comp)x*exp(-floor(y)*m_ln2);}
	void q_qr_bitwise_shift_right	(Quat &r, Quat const &x, Real const &y)	{r=(Quat)x*exp(-floor(y)*m_ln2);}
	void q_qc_bitwise_shift_right	(Quat &r, Quat const &x, Comp const &y)	{r=(Quat)x*exp(-floor(y)*m_ln2);}
	void q_qq_bitwise_shift_right	(Quat &r, Quat const &x, Quat const &y)	{r=(Quat)x*exp(-floor(y)*m_ln2);}

	inline Real bitwise_not(Real const &x){return isnan(x)|isinf(x)?0:~x.toLLong();}
	void  r_r_bitwise_not			(Real &r, Real const &x)				{r=bitwise_not(x);}
	void  c_c_bitwise_not			(Comp &r, Comp const &x)				{r=Comp(bitwise_not(x.r), bitwise_not(x.i));}
	void  q_q_bitwise_not			(Quat &r, Quat const &x)				{r=Quat(bitwise_not(x.r), bitwise_not(x.i), bitwise_not(x.j), bitwise_not(x.k));}

	inline Real bitwise_and(Real const &x){return isnan(x)|isinf(x)?0:!~x.toLLong();}
	inline Real bitwise_and(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:Real(x.toULLong()&y.toULLong());}
	inline long long bitwise_and_ll(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:x.toULLong()&y.toULLong();}
	void  r_r_bitwise_and			(Real &r, Real const &x)				{r=bitwise_and(x);}
	void  c_c_bitwise_and			(Comp &r, Comp const &x)				{r=Comp(bitwise_and(x.r), bitwise_and(x.i));}
	void  q_q_bitwise_and			(Quat &r, Quat const &x)				{r=Quat(bitwise_and(x.r), bitwise_and(x.i), bitwise_and(x.j), bitwise_and(x.k));}
	void r_rr_bitwise_and			(Real &r, Real const &x, Real const &y)	{r=bitwise_and(x, y);}
	void c_rc_bitwise_and			(Comp &r, Real const &x, Comp const &y)	{r=Comp(bitwise_and(x, y.r), bitwise_and(x, y.i));}
	void q_rq_bitwise_and			(Quat &r, Real const &x, Quat const &y)	{r=Quat(bitwise_and(x, y.r), bitwise_and(x, y.i), bitwise_and(x, y.j), bitwise_and(x, y.k));}
	void c_cr_bitwise_and			(Comp &r, Comp const &x, Real const &y)	{r=Comp(bitwise_and(x.r, y), bitwise_and(x.i, y));}
	void c_cc_bitwise_and			(Comp &r, Comp const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i);
		r=Comp(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr));
	}
	void q_cq_bitwise_and			(Quat &r, Comp const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i),
			xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j),
			xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xr_yj-xi_yk), Real(xr_yk+xi_yj));
	}
	void q_qr_bitwise_and			(Quat &r, Quat const &x, Real const &y)	{r=Quat(bitwise_and(x.r, y), bitwise_and(x.i, y), bitwise_and(x.j, y), bitwise_and(x.k, y));}
	void q_qc_bitwise_and			(Quat &r, Quat const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xj_yr+xk_yi), Real(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_and			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_and_ll(x.r, y.r), xi_yr=bitwise_and_ll(x.i, y.r), xj_yr=bitwise_and_ll(x.j, y.r), xk_yr=bitwise_and_ll(x.k, y.r),
			xr_yi=bitwise_and_ll(x.r, y.i), xi_yi=bitwise_and_ll(x.i, y.i), xj_yi=bitwise_and_ll(x.j, y.i), xk_yi=bitwise_and_ll(x.k, y.i),
			xr_yj=bitwise_and_ll(x.r, y.j), xi_yj=bitwise_and_ll(x.i, y.j), xj_yj=bitwise_and_ll(x.j, y.j), xk_yj=bitwise_and_ll(x.k, y.j),
			xr_yk=bitwise_and_ll(x.r, y.k), xi_yk=bitwise_and_ll(x.i, y.k), xj_yk=bitwise_and_ll(x.j, y.k), xk_yk=bitwise_and_ll(x.k, y.k);
		r=Quat(Real(xr_yr-xi_yi-xj_yj-xk_yk), Real(xr_yi+xi_yr+xj_yk-xk_yj), Real(xj_yj-xi_yk+xj_yr+xk_yi), Real(xr_yk+xi_yj-xj_yi+xk_yr));
	}

	inline Real bitwise_nand(Real const &x){return isnan(x)|isinf(x)?0:!x;}
	inline Real bitwise_nand(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:~(x.toULLong()&y.toULLong());}
	inline long long bitwise_nand_ll(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:~(x.toULLong()&y.toULLong());}
	void  r_r_bitwise_nand			(Real &r, Real const &x)				{r=bitwise_nand(x);}
	void  c_c_bitwise_nand			(Comp &r, Comp const &x)				{r=Comp(bitwise_nand(x.r), bitwise_nand(x.i));}
	void  q_q_bitwise_nand			(Quat &r, Quat const &x)				{r=Quat(bitwise_nand(x.r), bitwise_nand(x.i), bitwise_nand(x.j), bitwise_nand(x.k));}
	void r_rr_bitwise_nand			(Real &r, Real const &x, Real const &y)	{r=bitwise_nand(x, y);}
	void c_rc_bitwise_nand			(Comp &r, Real const &x, Comp const &y)	{r=Comp(bitwise_nand(x, y.r), bitwise_nand(x, y.i));}
	void q_rq_bitwise_nand			(Quat &r, Real const &x, Quat const &y)	{r=Quat(bitwise_nand(x, y.r), bitwise_nand(x, y.i), bitwise_nand(x, y.j), bitwise_nand(x, y.k));}
	void c_cr_bitwise_nand			(Comp &r, Comp const &x, Real const &y)	{r=Comp(bitwise_nand(x.r, y), bitwise_nand(x.i, y));}
	void c_cc_bitwise_nand			(Comp &r, Comp const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i);
		r=Comp(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr));
	}
	void q_cq_bitwise_nand			(Quat &r, Comp const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i),
			xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j),
			xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xr_yj-xi_yk), Real(xr_yk+xi_yj));
	}
	void q_qr_bitwise_nand			(Quat &r, Quat const &x, Real const &y)	{r=Quat(bitwise_nand(x.r, y), bitwise_nand(x.i, y), bitwise_nand(x.j, y), bitwise_nand(x.k, y));}
	void q_qc_bitwise_nand			(Quat &r, Quat const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i);
		r=Quat(Real(xr_yr-xi_yi), Real(xr_yi+xi_yr), Real(xj_yr+xk_yi), Real(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_nand			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nand_ll(x.r, y.r), xi_yr=bitwise_nand_ll(x.i, y.r), xj_yr=bitwise_nand_ll(x.j, y.r), xk_yr=bitwise_nand_ll(x.k, y.r),
			xr_yi=bitwise_nand_ll(x.r, y.i), xi_yi=bitwise_nand_ll(x.i, y.i), xj_yi=bitwise_nand_ll(x.j, y.i), xk_yi=bitwise_nand_ll(x.k, y.i),
			xr_yj=bitwise_nand_ll(x.r, y.j), xi_yj=bitwise_nand_ll(x.i, y.j), xj_yj=bitwise_nand_ll(x.j, y.j), xk_yj=bitwise_nand_ll(x.k, y.j),
			xr_yk=bitwise_nand_ll(x.r, y.k), xi_yk=bitwise_nand_ll(x.i, y.k), xj_yk=bitwise_nand_ll(x.j, y.k), xk_yk=bitwise_nand_ll(x.k, y.k);
		r=Quat(Real(xr_yr-xi_yi-xj_yj-xk_yk), Real(xr_yi+xi_yr+xj_yk-xk_yj), Real(xj_yj-xi_yk+xj_yr+xk_yi), Real(xr_yk+xi_yj-xj_yi+xk_yr));
	}
	
	inline Real bitwise_or(Real const &x){return isnan(x)|isinf(x)?0:x.toULLong()!=0;}
	inline Real bitwise_or(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:Real(x.toULLong()|y.toULLong());}
	inline long long bitwise_or_ll_c(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:~x.toULLong()&~y.toULLong();}
	void  r_r_bitwise_or			(Real &r, Real const &x)				{r=bitwise_or(x);}
	void  c_c_bitwise_or			(Comp &r, Comp const &x)				{r=Comp(bitwise_or(x.r), bitwise_or(x.i));}
	void  q_q_bitwise_or			(Quat &r, Quat const &x)				{r=Quat(bitwise_or(x.r), bitwise_or(x.i), bitwise_or(x.j), bitwise_or(x.k));}
	void r_rr_bitwise_or			(Real &r, Real const &x, Real const &y)	{r=bitwise_or(x, y);}
	void c_rc_bitwise_or			(Comp &r, Real const &x, Comp const &y)	{r=Comp(bitwise_or(x, y.r), bitwise_or(x, y.i));}
	void q_rq_bitwise_or			(Quat &r, Real const &x, Quat const &y)	{r=Quat(bitwise_or(x, y.r), bitwise_or(x, y.i), bitwise_or(x, y.j), bitwise_or(x, y.k));}
	void c_cr_bitwise_or			(Comp &r, Comp const &x, Real const &y)	{r=Comp(bitwise_or(x.r, y), bitwise_or(x.i, y));}
	void c_cc_bitwise_or			(Comp &r, Comp const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i);
		r=Comp(~(xr_yr-xi_yi), ~(xr_yi+xi_yr));
	}
	void q_cq_bitwise_or			(Quat &r, Comp const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i),
			xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j),
			xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k);
		r=Quat(~(xr_yr-xi_yi), ~(xr_yi+xi_yr), ~(xr_yj-xi_yk), ~(xr_yk+xi_yj));
	}
	void q_qr_bitwise_or			(Quat &r, Quat const &x, Real const &y)	{r=Quat(bitwise_or(x.r, y), bitwise_or(x.i, y), bitwise_or(x.j, y), bitwise_or(x.k, y));}
	void q_qc_bitwise_or			(Quat &r, Quat const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i);
		r=Quat(~(xr_yr-xi_yi), ~(xr_yi+xi_yr), ~(xj_yr+xk_yi), ~(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_or			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_or_ll_c(x.r, y.r), xi_yr=bitwise_or_ll_c(x.i, y.r), xj_yr=bitwise_or_ll_c(x.j, y.r), xk_yr=bitwise_or_ll_c(x.k, y.r),
			xr_yi=bitwise_or_ll_c(x.r, y.i), xi_yi=bitwise_or_ll_c(x.i, y.i), xj_yi=bitwise_or_ll_c(x.j, y.i), xk_yi=bitwise_or_ll_c(x.k, y.i),
			xr_yj=bitwise_or_ll_c(x.r, y.j), xi_yj=bitwise_or_ll_c(x.i, y.j), xj_yj=bitwise_or_ll_c(x.j, y.j), xk_yj=bitwise_or_ll_c(x.k, y.j),
			xr_yk=bitwise_or_ll_c(x.r, y.k), xi_yk=bitwise_or_ll_c(x.i, y.k), xj_yk=bitwise_or_ll_c(x.j, y.k), xk_yk=bitwise_or_ll_c(x.k, y.k);
		r=Quat(~(xr_yr-xi_yi-xj_yj-xk_yk), ~(xr_yi+xi_yr+xj_yk-xk_yj), ~(xj_yj-xi_yk+xj_yr+xk_yi), ~(xr_yk+xi_yj-xj_yi+xk_yr));
	}
	
	inline Real bitwise_nor(Real const &x){return isnan(x)|isinf(x)?0:!x.toULLong();}
	inline Real bitwise_nor(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(x)|isinf(x)?0:~(x.toULLong()|y.toULLong());}
	inline long long bitwise_nor_ll_c(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(x)|isinf(x)?0:~x.toULLong()|~y.toULLong();}
	void  r_r_bitwise_nor			(Real &r, Real const &x)				{r=bitwise_nor(x);}
	void  c_c_bitwise_nor			(Comp &r, Comp const &x)				{r=Comp(bitwise_nor(x.r), bitwise_nor(x.i));}
	void  q_q_bitwise_nor			(Quat &r, Quat const &x)				{r=Quat(bitwise_nor(x.r), bitwise_nor(x.i), bitwise_nor(x.j), bitwise_nor(x.k));}
	void r_rr_bitwise_nor			(Real &r, Real const &x, Real const &y)	{r=bitwise_nor(x, y);}
	void c_rc_bitwise_nor			(Comp &r, Real const &x, Comp const &y)	{r=Comp(bitwise_nor(x, y.r), bitwise_nor(x, y.i));}
	void q_rq_bitwise_nor			(Quat &r, Real const &x, Quat const &y)	{r=Quat(bitwise_nor(x, y.r), bitwise_nor(x, y.i), bitwise_nor(x, y.j), bitwise_nor(x, y.k));}
	void c_cr_bitwise_nor			(Comp &r, Comp const &x, Real const &y)	{r=Comp(bitwise_nor(x.r, y), bitwise_nor(x.i, y));}
	void c_cc_bitwise_nor			(Comp &r, Comp const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i);
		r=Comp(~(xr_yr-xi_yi), ~(xr_yi+xi_yr));
	}
	void q_cq_bitwise_nor			(Quat &r, Comp const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i),
			xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j),
			xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k);
		r=Quat((Real)~(xr_yr-xi_yi), (Real)~(xr_yi+xi_yr), (Real)~(xr_yj-xi_yk), (Real)~(xr_yk+xi_yj));
	}
	void q_qr_bitwise_nor			(Quat &r, Quat const &x, Real const &y)	{r=Quat(bitwise_nor(x.r, y), bitwise_nor(x.i, y), bitwise_nor(x.j, y), bitwise_nor(x.k, y));}
	void q_qc_bitwise_nor			(Quat &r, Quat const &x, Comp const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i);
		r=Quat((Real)~(xr_yr-xi_yi), (Real)~(xr_yi+xi_yr), (Real)~(xj_yr+xk_yi), (Real)~(-xj_yi+xk_yr));
	}
	void q_qq_bitwise_nor			(Quat &r, Quat const &x, Quat const &y)
	{
		long long
			xr_yr=bitwise_nor_ll_c(x.r, y.r), xi_yr=bitwise_nor_ll_c(x.i, y.r), xj_yr=bitwise_nor_ll_c(x.j, y.r), xk_yr=bitwise_nor_ll_c(x.k, y.r),
			xr_yi=bitwise_nor_ll_c(x.r, y.i), xi_yi=bitwise_nor_ll_c(x.i, y.i), xj_yi=bitwise_nor_ll_c(x.j, y.i), xk_yi=bitwise_nor_ll_c(x.k, y.i),
			xr_yj=bitwise_nor_ll_c(x.r, y.j), xi_yj=bitwise_nor_ll_c(x.i, y.j), xj_yj=bitwise_nor_ll_c(x.j, y.j), xk_yj=bitwise_nor_ll_c(x.k, y.j),
			xr_yk=bitwise_nor_ll_c(x.r, y.k), xi_yk=bitwise_nor_ll_c(x.i, y.k), xj_yk=bitwise_nor_ll_c(x.j, y.k), xk_yk=bitwise_nor_ll_c(x.k, y.k);
		r=Quat((Real)~(xr_yr-xi_yi-xj_yj-xk_yk), (Real)~(xr_yi+xi_yr+xj_yk-xk_yj), (Real)~(xj_yj-xi_yk+xj_yr+xk_yi), (Real)~(xr_yk+xi_yj-xj_yi+xk_yr));
	}
	
	inline Real bitwise_xor(Real const &x){return isnan(x)|isinf(x)?0:bitwise_xor(x);}
	inline Real bitwise_xor(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:Real(x.toULLong()^y.toULLong());}
	void  r_r_bitwise_xor			(Real &r, Real const &x)				{r=bitwise_xor(x);}
	void  c_c_bitwise_xor			(Comp &r, Comp const &x)				{r=Comp(bitwise_xor(x.r), bitwise_xor(x.i));}
	void  q_q_bitwise_xor			(Quat &r, Quat const &x)				{r=Quat(bitwise_xor(x.r), bitwise_xor(x.i), bitwise_xor(x.j), bitwise_xor(x.k));}
	void r_rr_bitwise_xor			(Real &r, Real const &x, Real const &y)	{r=bitwise_xor(x, y);}
	void c_rc_bitwise_xor			(Comp &r, Real const &x, Comp const &y)	{r=Comp(bitwise_xor(x, y.r), y.i);}
	void q_rq_bitwise_xor			(Quat &r, Real const &x, Quat const &y)	{r=Quat(bitwise_xor(x, y.r), y.i, y.j, y.k);}
	void c_cr_bitwise_xor			(Comp &r, Comp const &x, Real const &y)	{r=Comp(bitwise_xor(x.r, y), x.i);}
	void c_cc_bitwise_xor			(Comp &r, Comp const &x, Comp const &y)	{r=Comp(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i));}
	void q_cq_bitwise_xor			(Quat &r, Comp const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), y.j, y.k);}
	void q_qr_bitwise_xor			(Quat &r, Quat const &x, Real const &y)	{r=Quat(bitwise_xor(x.r, y), x.i, x.j, x.k);}
	void q_qc_bitwise_xor			(Quat &r, Quat const &x, Comp const &y)	{r=Quat(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), x.j, x.k);}
	void q_qq_bitwise_xor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xor(x.r, y.r), bitwise_xor(x.i, y.i), bitwise_xor(x.j, y.j), bitwise_xor(x.k, y.k));}
	
	inline Real bitwise_xnor(Real const &x){return isnan(x)|isinf(x)?0:!bitwise_xor(x);}
	inline Real bitwise_xnor(Real const &x, Real const &y){return isnan(x)|isinf(x)|isnan(y)|isinf(y)?0:(Real)~(x.toULLong()^y.toULLong());}
	void  r_r_bitwise_xnor			(Real &r, Real const &x)				{r=bitwise_xnor(x);}
	void  c_c_bitwise_xnor			(Comp &r, Comp const &x)				{r=Comp(bitwise_xnor(x.r), bitwise_xnor(x.i));}
	void  q_q_bitwise_xnor			(Quat &r, Quat const &x)				{r=Quat(bitwise_xnor(x.r), bitwise_xnor(x.i), bitwise_xnor(x.j), bitwise_xnor(x.k));}
	void r_rr_bitwise_xnor			(Real &r, Real const &x, Real const &y)	{r=bitwise_xnor(x, y);}
	void c_rc_bitwise_xnor			(Comp &r, Real const &x, Comp const &y)	{r=Comp(bitwise_xnor(x, y.r), y.i);}
	void q_rq_bitwise_xnor			(Quat &r, Real const &x, Quat const &y)	{r=Quat(bitwise_xnor(x, y.r), y.i, y.j, y.k);}
	void c_cr_bitwise_xnor			(Comp &r, Comp const &x, Real const &y)	{r=Comp(bitwise_xnor(x.r, y), x.i);}
	void c_cc_bitwise_xnor			(Comp &r, Comp const &x, Comp const &y)	{r=Comp(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i));}
	void q_cq_bitwise_xnor			(Quat &r, Comp const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), y.j, y.k);}
	void q_qr_bitwise_xnor			(Quat &r, Quat const &x, Real const &y)	{r=Quat(bitwise_xnor(x.r, y), x.i, x.j, x.k);}
	void q_qc_bitwise_xnor			(Quat &r, Quat const &x, Comp const &y)	{r=Quat(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), x.j, x.k);}
	void q_qq_bitwise_xnor			(Quat &r, Quat const &x, Quat const &y)	{r=Quat(bitwise_xnor(x.r, y.r), bitwise_xnor(x.i, y.i), bitwise_xnor(x.j, y.j), bitwise_xnor(x.k, y.k));}
	
	void  r_r_logic_equal			(Real &r, Real const &x)				{r=x==0;}
	void  r_c_logic_equal			(Real &r, Comp const &x)				{r=!x.c_is_true();}
	void  r_q_logic_equal			(Real &r, Quat const &x)				{r=!x.q_is_true();}
	void r_rr_logic_equal			(Real &r, Real const &x, Real const &y)	{r=x==y;}
	void r_rc_logic_equal			(Real &r, Real const &x, Comp const &y)	{r=x==y;}
	void r_rq_logic_equal			(Real &r, Real const &x, Quat const &y)	{r=x==y;}
	void r_cr_logic_equal			(Real &r, Comp const &x, Real const &y)	{r=x==y;}
	void r_cc_logic_equal			(Real &r, Comp const &x, Comp const &y)	{r=x==y;}
	void r_cq_logic_equal			(Real &r, Comp const &x, Quat const &y)	{r=x==y;}
	void r_qr_logic_equal			(Real &r, Quat const &x, Real const &y)	{r=x==y;}
	void r_qc_logic_equal			(Real &r, Quat const &x, Comp const &y)	{r=x==y;}
	void r_qq_logic_equal			(Real &r, Quat const &x, Quat const &y)	{r=x==y;}
	
	void  r_r_logic_not_equal		(Real &r, Real const &x)				{r=x!=0;}
	void  r_c_logic_not_equal		(Real &r, Comp const &x)				{r=x.c_is_true();}
	void  r_q_logic_not_equal		(Real &r, Quat const &x)				{r=x.q_is_true();}
	void r_rr_logic_not_equal		(Real &r, Real const &x, Real const &y)	{r=x!=y;}
	void r_rc_logic_not_equal		(Real &r, Real const &x, Comp const &y)	{r=x!=y;}
	void r_rq_logic_not_equal		(Real &r, Real const &x, Quat const &y)	{r=x!=y;}
	void r_cr_logic_not_equal		(Real &r, Comp const &x, Real const &y)	{r=x!=y;}
	void r_cc_logic_not_equal		(Real &r, Comp const &x, Comp const &y)	{r=x!=y;}
	void r_cq_logic_not_equal		(Real &r, Comp const &x, Quat const &y)	{r=x!=y;}
	void r_qr_logic_not_equal		(Real &r, Quat const &x, Real const &y)	{r=x!=y;}
	void r_qc_logic_not_equal		(Real &r, Quat const &x, Comp const &y)	{r=x!=y;}
	void r_qq_logic_not_equal		(Real &r, Quat const &x, Quat const &y)	{r=x!=y;}
	
	void  r_r_logic_less_l			(Real &r, Real const &x)				{r=0<x;}
	void  r_c_logic_less_l			(Real &r, Comp const &x)				{r=0<x.r;}
	void  r_q_logic_less_l			(Real &r, Quat const &x)				{r=0<x.r;}
	void  r_r_logic_less_r			(Real &r, Real const &x)				{r=x<0;}
	void  r_c_logic_less_r			(Real &r, Comp const &x)				{r=x.r<0;}
	void  r_q_logic_less_r			(Real &r, Quat const &x)				{r=x.r<0;}
	void r_rr_logic_less			(Real &r, Real const &x, Real const &y)	{r=x<y;}
	void r_rc_logic_less			(Real &r, Real const &x, Comp const &y)	{r=x<y.r;}
	void r_rq_logic_less			(Real &r, Real const &x, Quat const &y)	{r=x<y.r;}
	void r_cr_logic_less			(Real &r, Comp const &x, Real const &y)	{r=x.r<y;}
	void r_cc_logic_less			(Real &r, Comp const &x, Comp const &y)	{r=x.r<y.r;}
	void r_cq_logic_less			(Real &r, Comp const &x, Quat const &y)	{r=x.r<y.r;}
	void r_qr_logic_less			(Real &r, Quat const &x, Real const &y)	{r=x.r<y;}
	void r_qc_logic_less			(Real &r, Quat const &x, Comp const &y)	{r=x.r<y.r;}
	void r_qq_logic_less			(Real &r, Quat const &x, Quat const &y)	{r=x.r<y.r;}
	
	void  r_r_logic_less_equal_l	(Real &r, Real const &x)				{r=0<=x;}
	void  r_c_logic_less_equal_l	(Real &r, Comp const &x)				{r=0<=x.r;}
	void  r_q_logic_less_equal_l	(Real &r, Quat const &x)				{r=0<=x.r;}
	void  r_r_logic_less_equal_r	(Real &r, Real const &x)				{r=x<=0;}
	void  r_c_logic_less_equal_r	(Real &r, Comp const &x)				{r=x.r<=0;}
	void  r_q_logic_less_equal_r	(Real &r, Quat const &x)				{r=x.r<=0;}
	void r_rr_logic_less_equal		(Real &r, Real const &x, Real const &y)	{r=x<=y;}
	void r_rc_logic_less_equal		(Real &r, Real const &x, Comp const &y)	{r=x<=y.r;}
	void r_rq_logic_less_equal		(Real &r, Real const &x, Quat const &y)	{r=x<=y.r;}
	void r_cr_logic_less_equal		(Real &r, Comp const &x, Real const &y)	{r=x.r<=y;}
	void r_cc_logic_less_equal		(Real &r, Comp const &x, Comp const &y)	{r=x.r<=y.r;}
	void r_cq_logic_less_equal		(Real &r, Comp const &x, Quat const &y)	{r=x.r<=y.r;}
	void r_qr_logic_less_equal		(Real &r, Quat const &x, Real const &y)	{r=x.r<=y;}
	void r_qc_logic_less_equal		(Real &r, Quat const &x, Comp const &y)	{r=x.r<=y.r;}
	void r_qq_logic_less_equal		(Real &r, Quat const &x, Quat const &y)	{r=x.r<=y.r;}
	
	void  r_r_logic_greater_l		(Real &r, Real const &x)				{r=0>x;}
	void  r_c_logic_greater_l		(Real &r, Comp const &x)				{r=0>x.r;}
	void  r_q_logic_greater_l		(Real &r, Quat const &x)				{r=0>x.r;}
	void  r_r_logic_greater_r		(Real &r, Real const &x)				{r=x>0;}
	void  r_c_logic_greater_r		(Real &r, Comp const &x)				{r=x.r>0;}
	void  r_q_logic_greater_r		(Real &r, Quat const &x)				{r=x.r>0;}
	void r_rr_logic_greater			(Real &r, Real const &x, Real const &y)	{r=x>y;}
	void r_rc_logic_greater			(Real &r, Real const &x, Comp const &y)	{r=x>y.r;}
	void r_rq_logic_greater			(Real &r, Real const &x, Quat const &y)	{r=x>y.r;}
	void r_cr_logic_greater			(Real &r, Comp const &x, Real const &y)	{r=x.r>y;}
	void r_cc_logic_greater			(Real &r, Comp const &x, Comp const &y)	{r=x.r>y.r;}
	void r_cq_logic_greater			(Real &r, Comp const &x, Quat const &y)	{r=x.r>y.r;}
	void r_qr_logic_greater			(Real &r, Quat const &x, Real const &y)	{r=x.r>y;}
	void r_qc_logic_greater			(Real &r, Quat const &x, Comp const &y)	{r=x.r>y.r;}
	void r_qq_logic_greater			(Real &r, Quat const &x, Quat const &y)	{r=x.r>y.r;}
	
	void  r_r_logic_greater_equal_l	(Real &r, Real const &x)				{r=0>=x;}
	void  r_c_logic_greater_equal_l	(Real &r, Comp const &x)				{r=0>=x.r;}
	void  r_q_logic_greater_equal_l	(Real &r, Quat const &x)				{r=0>=x.r;}
	void  r_r_logic_greater_equal_r	(Real &r, Real const &x)				{r=x>=0;}
	void  r_c_logic_greater_equal_r	(Real &r, Comp const &x)				{r=x.r>=0;}
	void  r_q_logic_greater_equal_r	(Real &r, Quat const &x)				{r=x.r>=0;}
	void r_rr_logic_greater_equal	(Real &r, Real const &x, Real const &y)	{r=x>=y;}
	void r_rc_logic_greater_equal	(Real &r, Real const &x, Comp const &y)	{r=x>=y.r;}
	void r_rq_logic_greater_equal	(Real &r, Real const &x, Quat const &y)	{r=x>=y.r;}
	void r_cr_logic_greater_equal	(Real &r, Comp const &x, Real const &y)	{r=x.r>=y;}
	void r_cc_logic_greater_equal	(Real &r, Comp const &x, Comp const &y)	{r=x.r>=y.r;}
	void r_cq_logic_greater_equal	(Real &r, Comp const &x, Quat const &y)	{r=x.r>=y.r;}
	void r_qr_logic_greater_equal	(Real &r, Quat const &x, Real const &y)	{r=x.r>=y;}
	void r_qc_logic_greater_equal	(Real &r, Quat const &x, Comp const &y)	{r=x.r>=y.r;}
	void r_qq_logic_greater_equal	(Real &r, Quat const &x, Quat const &y)	{r=x.r>=y.r;}
	
	void  r_r_logic_not				(Real &r, Real const &x)				{r=x==0;}
	void  r_c_logic_not				(Real &r, Comp const &x)				{r=(x.r==0)&(x.i==0);}
	void  r_q_logic_not				(Real &r, Quat const &x)				{r=(x.r==0)&(x.i==0)&(x.j==0)&(x.k==0);}
	
	void r_rr_logic_and				(Real &r, Real const &x, Real const &y)	{r=x==0&y==0;}
	void r_rc_logic_and				(Real &r, Real const &x, Comp const &y)	{r=x==0&y.c_is_true();}
	void r_rq_logic_and				(Real &r, Real const &x, Quat const &y)	{r=x==0&y.q_is_true();}
	void r_cr_logic_and				(Real &r, Comp const &x, Real const &y)	{r=x.c_is_true()&y==0;}
	void r_cc_logic_and				(Real &r, Comp const &x, Comp const &y)	{r=x.c_is_true()&y.c_is_true();}
	void r_cq_logic_and				(Real &r, Comp const &x, Quat const &y)	{r=x.c_is_true()&y.q_is_true();}
	void r_qr_logic_and				(Real &r, Quat const &x, Real const &y)	{r=x.q_is_true()&y==0;}
	void r_qc_logic_and				(Real &r, Quat const &x, Comp const &y)	{r=x.q_is_true()&y.c_is_true();}
	void r_qq_logic_and				(Real &r, Quat const &x, Quat const &y)	{r=x.q_is_true()&y.q_is_true();}
	
	void r_rr_logic_or				(Real &r, Real const &x, Real const &y)	{r=x==0|y==0;}
	void r_rc_logic_or				(Real &r, Real const &x, Comp const &y)	{r=x==0|y.c_is_true();}
	void r_rq_logic_or				(Real &r, Real const &x, Quat const &y)	{r=x==0|y.q_is_true();}
	void r_cr_logic_or				(Real &r, Comp const &x, Real const &y)	{r=x.c_is_true()|y==0;}
	void r_cc_logic_or				(Real &r, Comp const &x, Comp const &y)	{r=x.c_is_true()|y.c_is_true();}
	void r_cq_logic_or				(Real &r, Comp const &x, Quat const &y)	{r=x.c_is_true()|y.q_is_true();}
	void r_qr_logic_or				(Real &r, Quat const &x, Real const &y)	{r=x.q_is_true()|y==0;}
	void r_qc_logic_or				(Real &r, Quat const &x, Comp const &y)	{r=x.q_is_true()|y.c_is_true();}
	void r_qq_logic_or				(Real &r, Quat const &x, Quat const &y)	{r=x.q_is_true()|y.q_is_true();}
	
	void r_rr_logic_xor				(Real &r, Real const &x, Real const &y)	{r=x==0^y==0;}
	void r_rc_logic_xor				(Real &r, Real const &x, Comp const &y)	{r=x==0^y.c_is_true();}
	void r_rq_logic_xor				(Real &r, Real const &x, Quat const &y)	{r=x==0^y.q_is_true();}
	void r_cr_logic_xor				(Real &r, Comp const &x, Real const &y)	{r=x.c_is_true()^y==0;}
	void r_cc_logic_xor				(Real &r, Comp const &x, Comp const &y)	{r=x.c_is_true()^y.c_is_true();}
	void r_cq_logic_xor				(Real &r, Comp const &x, Quat const &y)	{r=x.c_is_true()^y.q_is_true();}
	void r_qr_logic_xor				(Real &r, Quat const &x, Real const &y)	{r=x.q_is_true()^y==0;}
	void r_qc_logic_xor				(Real &r, Quat const &x, Comp const &y)	{r=x.q_is_true()^y.c_is_true();}
	void r_qq_logic_xor				(Real &r, Quat const &x, Quat const &y)	{r=x.q_is_true()^y.q_is_true();}
	
	void r_rr_condition_zero		(Real &r, Real const &x, Real const &y)	{r=x==0?x:y;}
	void c_rc_condition_zero		(Comp &r, Real const &x, Comp const &y)	{r=x==0?Comp(x, 0):y;}
	void q_rq_condition_zero		(Quat &r, Real const &x, Quat const &y)	{r=x==0?Quat(x, 0, 0, 0):y;}
	void c_cr_condition_zero		(Comp &r, Comp const &x, Real const &y)	{r=x.c_is_true()?x:Comp(y, 0);}
	void c_cc_condition_zero		(Comp &r, Comp const &x, Comp const &y)	{r=x.c_is_true()?x:y;}
	void q_cq_condition_zero		(Quat &r, Comp const &x, Quat const &y)	{r=x.c_is_true()?Quat(x):y;}
	void q_qr_condition_zero		(Quat &r, Quat const &x, Real const &y)	{r=x.q_is_true()?x:Quat(y, 0, 0, 0);}
	void q_qc_condition_zero		(Quat &r, Quat const &x, Comp const &y)	{r=x.q_is_true()?x:Quat(y);}
	void q_qq_condition_zero		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?x:y;}

	void  r_r_percent				(Real &r, Real const &x)				{r=x*0.01;}
	void  c_c_percent				(Comp &r, Comp const &x)				{r=(Comp)x*0.01;}
	void  q_q_percent				(Quat &r, Quat const &x)				{r=(Quat)x*0.01;}
	
	void r_rr_modulo				(Real &r, Real const &x, Real const &y)	{r=x%y;}
	void c_rc_modulo				(Comp &r, Real const &x, Comp const &y)	{r=x%y;}
	void q_rq_modulo				(Quat &r, Real const &x, Quat const &y)	{r=x%y;}
	void c_cr_modulo				(Comp &r, Comp const &x, Real const &y)	{r=x%y;}
	void c_cc_modulo				(Comp &r, Comp const &x, Comp const &y)	{r=x%y;}
	void q_cq_modulo				(Quat &r, Comp const &x, Quat const &y)	{r=x%y;}
	void q_qr_modulo				(Quat &r, Quat const &x, Real const &y)	{r=x%y;}
	void q_qc_modulo				(Quat &r, Quat const &x, Comp const &y)	{r=x%y;}
	void q_qq_modulo				(Quat &r, Quat const &x, Quat const &y)	{r=x%y;}

	void  r_r_sgn					(Real &r, Real const &x)				{r=(x>0)-(x<0);}
	void  c_c_sgn					(Comp &r, Comp const &x)				{Comp cx=x; Real mag=abs(cx); r=mag?cx/mag:Comp();}
	void  q_q_sgn					(Quat &r, Quat const &x)				{Quat qx=x; Real mag=abs(qx); r=mag?qx/mag:Quat();}
	
	inline Comp	sq(Comp const &x){Real ri=x.r*x.i; return Comp(x.r*x.r-x.i*x.i, ri+ri);}
	inline Quat	sq(Quat const &x)
	{
		auto _2r=x.r+x.r;
		return Quat(x.r*x.r-x.i*x.i-x.j*x.j-x.k*x.k, x.i*_2r, x.j*_2r, x.k*_2r);
	}
	void  r_r_sq					(Real &r, Real const &x)				{r=x*x;}
	void  c_c_sq					(Comp &r, Comp const &x)				{r=sq(x);}
	void  q_q_sq					(Quat &r, Quat const &x)				{r=sq(x);}
	
	void  c_c_sqrt					(Comp &r, Comp const &x)				{r=sqrt((Comp)x);}
	void  q_q_sqrt					(Quat &r, Quat const &x)				{r=sqrt((Quat)x);}

	void  r_r_invsqrt				(Real &r, Real const &x)				{r=rec_sqrt(x);}

	void  r_r_cbrt					(Real &r, Real const &x)				{r=cbrt(x);}
	void  c_c_cbrt					(Comp &r, Comp const &x)				{r=exp(1./3*log(x));}//optimize
	void  q_q_cbrt					(Quat &r, Quat const &x)
	{
		Real t=x.i*x.i+x.j*x.j+x.k*x.k;
		if(t)
		{
			Real absx=sqrt(x.r*x.r+t);
			t=acos(x.r/absx)/(sqrt(t)*3);
			Quat logx_3((1./3)*log(absx), t*x.i, t*x.j, t*x.k);
			r=exp(logx_3);
		}
		else
			r=Quat(cbrt(x.r));
	}

	void  r_r_gauss					(Real &r, Real const &x)				{Real rx=x; r=exp(-rx*rx);}
	void  c_c_gauss					(Comp &r, Comp const &x)				{r=exp(-sq(x));}
	void  q_q_gauss					(Quat &r, Quat const &x)				{r=exp(-sq(x));}

	void  r_r_erf					(Real &r, Real const &x)				{r=erf((Real)x);}

	void  r_r_zeta					(Real &r, Real const &x)				{r=zeta(x);}
	
	const Comp m_i(0, 1);
	__forceinline Quat sgnu(Quat const &x){return Quat(0, x.i, x.j, x.k);}
	__forceinline Quat acosh(Quat const &x){return log(x+sqrt(sq(x)-1));}
	__forceinline Comp asinh(Comp const &x){return log(x+sqrt(sq(x)+1));}
	__forceinline Quat asinh(Quat const &x){return log(x+sqrt(sq(x)+1));}
	__forceinline Real sinhc(Real const &x){return sinh(x)/x;}
	__forceinline Comp cos(Comp const &x)
	{
		Comp exp_ix=exp(m_i*x);//sincos, exp
		return 0.5*exp_ix+0.5/exp_ix;
	}
	__forceinline Quat cos(Quat const &x)
	{
		Real z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k), sin_xr, cos_xr;//boost::math: sqrt, sincos, sh, ch
		sin_cos(sin_xr, cos_xr, x.r);
		Real w=-sin_xr*sinhc(z);
		return Quat(cos_xr*cosh(z), w*x.i, w*x.j, w*x.k);
	}
	__forceinline Comp sin(Comp const &x)
	{
		Comp exp_ix=exp(m_i*x);
		return 0.5*exp_ix-0.5/exp_ix;
	//	return (exp(m_i*x)-exp(-m_i*x))*(-m_half*m_i);
	}
	__forceinline Quat sin(Quat const &x)//boost::math
	{
		Real z=sqrt(x.i*x.i+x.j*x.j+x.k*x.k);
		Real sin_xr, cos_xr;
		sin_cos(sin_xr, cos_xr, x.r);
		Real w=-cos_xr*sinhc(z);
		return Quat(sin_xr*cosh(z), w*x.i, w*x.j, w*x.k);
	}
	namespace gamma//http://en.wikipedia.org/wiki/Lanczos_approximation
	{
		const double g=7, p[]={0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	}
	Comp			tgamma		(Comp const &x)
	{
		using namespace gamma;
		if(x.r<.5)
		{
			Comp t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			Comp t2=g+.5-x;
			return m_pi/(sin(m_pi*x)*m_sqrt_2pi*(t2^0.5-x)*exp(-t2)*t1);
		}
		else
		{
			Comp t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			Comp t2=g+.5+x-1.;
			return m_sqrt_2pi*(t2^0.5+x-1)*exp(-t2)*t1;
		}
	}
	Quat			tgamma		(Quat	const &x)
	{
		using namespace gamma;
		if(x.r<.5)
		{
			Quat t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			Quat t2=g+.5-x;
			return m_pi/(sin(m_pi*x)*m_sqrt_2pi*(t2^0.5-x)*exp(-t2)*t1);
		}
		else
		{
			Quat t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			Quat t2=g+.5+x-1.;
			return m_sqrt_2pi*(t2^0.5+x-1.)*exp(-t2)*t1;
		}
	}
	void  r_r_tgamma				(Real &r, Real const &x)				{r=tgamma(x);}
	void  c_c_tgamma				(Comp &r, Comp const &x)				{r=tgamma(x);}
	void  q_q_tgamma				(Quat &r, Quat const &x)				{r=tgamma(x);}
	double tgamma(double const &x, double const &y)//upper incomplete gamma, G(x, 0) = G(x)
	{
		try
		{
			if(x>0)
			{
				if(y>=0)
					return boost::math::tgamma(x, y);
			}
			else if(x==0)
				return _HUGE;
			long long cInf=0x7FF8000000000010;//x.r>0&&y.r<0||x.r<0
			return (double&)cInf;
			//if(x.r<0)
			//{
			//	long long cInf=0x7FF8000000000010;
			//	return (Real&)cInf;
			//}
			//else if(x.r==0)
			//	return _HUGE;
			//else if(y.r<0)
			//{
			//	long long cInf=0x7FF8000000000010;
			//	return (Real&)cInf;
			//}
			//return boost::math::tgamma((Real)x, (Real)y);
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
	void r_rr_tgamma				(Real &r, Real const &x, Real const &y)	{r=tgamma(x.toDouble(), y.toDouble());}

	void  r_r_loggamma				(Real &r, Real const &x)				{r=lgamma(x);}

	void  r_r_factorial				(Real &r, Real const &x)				{r=tgamma(x+1);}
	void  c_c_factorial				(Comp &r, Comp const &x)				{r=tgamma((Comp)x+1.);}
	void  q_q_factorial				(Quat &r, Quat const &x)				{r=tgamma((Quat)x+1.);}

	inline Real permutation(Real const &x, Real const &y){return tgamma(x+1)/tgamma(x-y+1);}
	Comp permutation(Comp const &x, Comp const &y){return tgamma(x+1.)/tgamma(x-y+1.);}
	Quat permutation(Quat const &x, Quat const &y){return tgamma(x+1.)/tgamma(x-y+1.);}
	void  r_r_permutation			(Real &r, Real const &x)				{r=1;}
	void  c_c_permutation			(Comp &r, Comp const &x)				{r=Comp(1, 0);}
	void  q_q_permutation			(Quat &r, Quat const &x)				{r=Quat(1, 0, 0, 0);}
	void r_rr_permutation			(Real &r, Real const &x, Real const &y)	{r=permutation(x, y);}
	void c_cr_permutation			(Comp &r, Comp const &x, Real const &y)	{r=permutation(x, Comp(y, 0));}
	void c_cc_permutation			(Comp &r, Comp const &x, Comp const &y)	{r=permutation(x, y);}
	void q_qq_permutation			(Quat &r, Quat const &x, Quat const &y)	{r=permutation(x, y);}
	
	Real combination(Real const &x, Real const &y){return tgamma(x+1)/(tgamma(x-y+1)*tgamma(y+1));}
	Comp combination(Comp const &x, Comp const &y){return tgamma(x+1)/(tgamma(x-y+1)*tgamma(y+1));}
	Quat combination(Quat const &x, Quat const &y){return tgamma(x+1)/(tgamma(x-y+1)*tgamma(y+1));}
	void  r_r_combination			(Real &r, Real const &x)				{r=1;}
	void  c_c_combination			(Comp &r, Comp const &x)				{r=Comp(1, 0);}
	void  q_q_combination			(Quat &r, Quat const &x)				{r=Quat(1, 0, 0, 0);}
	void r_rr_combination			(Real &r, Real const &x, Real const &y)	{r=combination(x, y);}
	void c_cr_combination			(Comp &r, Comp const &x, Real const &y)	{r=combination(x, Comp(y, 0));}
	void c_cc_combination			(Comp &r, Comp const &x, Comp const &y)	{r=combination(x, y);}
	void q_qq_combination			(Quat &r, Quat const &x, Quat const &y)	{r=combination(x, y);}

//	const Comp _i(0, 1);
//	inline Comp asinh(Comp const &x){return log(x+sqrt(sq(x)+1.));}
//	inline Real sinhc(Real const &x){return sinh(x)/x;}
	void  r_r_cos					(Real &r, Real const &x)				{r=cos(x);}
	void  c_c_cos					(Comp &r, Comp const &x)				{r=cos(x);}
	void  q_q_cos					(Quat &r, Quat const &x)				{r=cos(x);}

	inline Comp acos(Comp const &x){return -m_i*log(x+sqrt(sq(x)-1.));}
	inline Quat acos(Quat const &x){return -m_i*log(x+sqrt(sq(x)-1.));}
	void  c_c_acos					(Comp &r, Comp const &x)				{r=acos(x);}
	void  q_q_acos					(Quat &r, Quat const &x)				{r=acos(x);}

	inline Comp cosh(Comp const &x){return (exp(x)+exp(-x))*0.5;}
	inline Quat cosh(Quat const &x){return (exp(x)+exp(-x))*0.5;}
	void  r_r_cosh					(Real &r, Real const &x)				{r=cosh(x);}
	void  c_c_cosh					(Comp &r, Comp const &x)				{r=cosh(x);}
	void  q_q_cosh					(Quat &r, Quat const &x)				{r=cosh(x);}
	
	inline Comp acosh(Comp const &x){return log(x+sqrt(sq(x)-1.));}
	void  c_c_acosh					(Comp &r, Comp const &x)				{r=acosh(x);}
	void  q_q_acosh					(Quat &r, Quat const &x)				{r=acosh(x);}

	void  r_r_cosc					(Real &r, Real const &x)				{r=cos(x)/x;}
	void  c_c_cosc					(Comp &r, Comp const &x)				{Comp cx=x; r=cos(cx)/cx;}
	void  q_q_cosc					(Quat &r, Quat const &x)				{Quat qx=x; r=cos(qx)/qx;}

	void  r_r_sec					(Real &r, Real const &x)				{r=1/cos(x);}
	void  c_c_sec					(Comp &r, Comp const &x)				{r=inv(cos(x));}
	void  q_q_sec					(Quat &r, Quat const &x)				{r=inv(cos(x));}

	void  c_c_asec					(Comp &r, Comp const &x)				{r=acos(inv(x));}
	void  q_q_asec					(Quat &r, Quat const &x)				{r=acos(inv(x));}

	void  r_r_sech					(Real &r, Real const &x)				{r=1/cosh(x);}
	void  c_c_sech					(Comp &r, Comp const &x)				{r=inv(cosh(x));}
	void  q_q_sech					(Quat &r, Quat const &x)				{r=inv(cosh(x));}

	void  c_c_asech					(Comp &r, Comp const &x)				{r=acosh(inv(x));}
	void  q_q_asech					(Quat &r, Quat const &x)				{r=acosh(inv(x));}

	void  r_r_sin					(Real &r, Real const &x)				{r=sin(x);}
	void  c_c_sin					(Comp &r, Comp const &x)				{r=sin(x);}
	void  q_q_sin					(Quat &r, Quat const &x)				{r=sin(x);}

	inline Comp asin(Comp const &x){return -m_i*log(m_i*x+sqrt(1.-sq(x)));}
	inline Quat asin(Quat const &x){return -m_i*log(m_i*x+sqrt(1.-sq(x)));}
	void  c_c_asin					(Comp &r, Comp const &x)				{r=asin(x);}
	void  q_q_asin					(Quat &r, Quat const &x)				{r=asin(x);}

	inline Comp sinh(Comp const &x){return (exp(x)-exp(-x))*0.5;}
	inline Quat sinh(Quat const &x){return (exp(x)-exp(-x))*0.5;}
	void  r_r_sinh					(Real &r, Real const &x)				{r=sinh(x);}
	void  c_c_sinh					(Comp &r, Comp const &x)				{r=sinh(x);}
	void  q_q_sinh					(Quat &r, Quat const &x)				{r=sinh(x);}

	void  r_r_asinh					(Real &r, Real const &x)				{r=asinh(x);}
	void  c_c_asinh					(Comp &r, Comp const &x)				{r=asinh(x);}
	void  q_q_asinh					(Quat &r, Quat const &x)				{r=asinh(x);}

	void  r_r_sinc					(Real &r, Real const &x)				{r=x?sin(x)/x:1;}
	void  c_c_sinc					(Comp &r, Comp const &x)				{r=x.c_is_true()?sin(x)/x:Comp(1);}
	void  q_q_sinc					(Quat &r, Quat const &x)				{r=x.q_is_true()?sin(x)/x:Quat(1);}

	void  r_r_sinhc					(Real &r, Real const &x)				{r=x?sinh(x)/x:1;}
	void  c_c_sinhc					(Comp &r, Comp const &x)				{r=x.c_is_true()?sinh(x)/x:Comp(1);}
	void  q_q_sinhc					(Quat &r, Quat const &x)				{r=x.q_is_true()?sinh(x)/x:Quat(1);}

	void  r_r_csc					(Real &r, Real const &x)				{r=1/sin(x);}
	void  c_c_csc					(Comp &r, Comp const &x)				{r=inv(sin((Comp)x));}
	void  q_q_csc					(Quat &r, Quat const &x)				{r=inv(sin((Quat)x));}

	void  c_c_acsc					(Comp &r, Comp const &x)				{r=asin(inv(x));}
	void  q_q_acsc					(Quat &r, Quat const &x)				{r=asin(inv(x));}

	void  r_r_csch					(Real &r, Real const &x)				{r=1/sinh(x);}
	void  c_c_csch					(Comp &r, Comp const &x)				{r=inv(sinh(x));}
	void  q_q_csch					(Quat &r, Quat const &x)				{r=inv(sinh(x));}

	void  r_r_acsch					(Real &r, Real const &x)				{r=asinh(1/x);}
	void  c_c_acsch					(Comp &r, Comp const &x)				{r=asinh(inv(x));}
	void  q_q_acsch					(Quat &r, Quat const &x)				{r=asinh(inv(x));}

	inline Comp tan(Comp const &x)
	{
		const Comp two_i(0, 2);
		Comp exp_2ix=exp(two_i*x);
		return (exp_2ix-1.)/((exp_2ix+1.)*m_i);
	}
	inline Quat tan(Quat const &x)
	{
		const Comp two_i(0, 2);
		Quat exp_2ix=exp(two_i*x);
		return (exp_2ix-1.)/((exp_2ix+1.)*m_i);
	}
	void  r_r_tan					(Real &r, Real const &x)				{r=tan(x);}
	void  c_c_tan					(Comp &r, Comp const &x)				{r=tan(x);}
	void  q_q_tan					(Quat &r, Quat const &x)				{r=tan(x);}

	inline Comp atan(Comp const &x){return (m_i*0.5)*log((m_i+x)/(m_i-x));}
	inline Quat atan(Quat const &x){return (m_i*0.5)*log((m_i+x)/(m_i-x));}
	void  r_r_atan					(Real &r, Real const &x)				{r=atan(x);}
	void  c_c_atan					(Comp &r, Comp const &x)				{r=atan(x);}
	void  q_q_atan					(Quat &r, Quat const &x)				{r=atan(x);}
	void r_rr_atan					(Real &r, Real const &x, Real const &y)	{r=atan2(x, y);}
	void c_rc_atan					(Comp &r, Real const &x, Comp const &y)	{Comp t=atan(y/x); r=x  <0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_rq_atan					(Quat &r, Real const &x, Quat const &y)	{Quat t=atan(y/x); r=x  <0?y.r<0?t-m_pi:t+m_pi:t;}
	void c_cr_atan					(Comp &r, Comp const &x, Real const &y)	{Comp t=atan(y/x); r=x.r<0?y  <0?t-m_pi:t+m_pi:t;}
	void c_cc_atan					(Comp &r, Comp const &x, Comp const &y)	{Comp t=atan(y/x); r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_cq_atan					(Quat &r, Comp const &x, Quat const &y)	{Quat t=atan(y/x); r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_qr_atan					(Quat &r, Quat const &x, Real const &y)	{Quat t=atan(y/x); r=x.r<0?y  <0?t-m_pi:t+m_pi:t;}
	void q_qc_atan					(Quat &r, Quat const &x, Comp const &y)	{Quat t=atan(y/x); r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	void q_qq_atan					(Quat &r, Quat const &x, Quat const &y)	{Quat t=atan(y/x); r=x.r<0?y.r<0?t-m_pi:t+m_pi:t;}
	
	inline Comp tanh(Comp const &x){Comp e2x=exp(x+x); return (e2x-1.)/(e2x+1.);}
	inline Quat tanh(Quat const &x){Quat e2x=exp(x+x); return (e2x-1.)/(e2x+1.);}
	void  r_r_tanh					(Real &r, Real const &x)				{r=tanh(x);}
	void  c_c_tanh					(Comp &r, Comp const &x)				{r=tanh(x);}
	void  q_q_tanh					(Quat &r, Quat const &x)				{r=tanh(x);}

	inline Comp atanh(Comp const &x){return 0.5*log((1+x)/(1-x));}
	inline Quat atanh(Quat const &x){return 0.5*log((1+x)/(1-x));}
	void  c_c_atanh					(Comp &r, Comp const &x)				{r=atanh(x);}
	void  q_q_atanh					(Quat &r, Quat const &x)				{r=atanh(x);}

	void  r_r_tanc					(Real &r, Real const &x)				{r=x?tan(x)/x:0;}
	void  c_c_tanc					(Comp &r, Comp const &x)				{r=x.c_is_true()?tan(x)/x:Comp();}
	void  q_q_tanc					(Quat &r, Quat const &x)				{r=x.q_is_true()?tan(x)/x:Quat();}

	void  r_r_cot					(Real &r, Real const &x)				{r=1/tan(x);}
	void  c_c_cot					(Comp &r, Comp const &x)				{r=inv(tan(x));}
	void  q_q_cot					(Quat &r, Quat const &x)				{r=inv(tan(x));}

	void  r_r_acot					(Real &r, Real const &x)				{r=x?atan(1/x):m_pi/2;}
	void  c_c_acot					(Comp &r, Comp const &x)				{r=atan(inv(x));}
	void  q_q_acot					(Quat &r, Quat const &x)				{r=atan(inv(x));}

	void  r_r_coth					(Real &r, Real const &x)				{r=1/tanh(x);}
	void  c_c_coth					(Comp &r, Comp const &x)				{r=inv(tanh(x));}
	void  q_q_coth					(Quat &r, Quat const &x)				{r=inv(tanh(x));}

	void  c_c_acoth					(Comp &r, Comp const &x)				{r=atanh(inv((Comp)x));}
	void  q_q_acoth					(Quat &r, Quat const &x)				{r=atanh(inv((Quat)x));}

	void  r_r_exp					(Real &r, Real const &x)				{r=exp(x);}
	void  c_c_exp					(Comp &r, Comp const &x)				{r=exp(x);}
	void  q_q_exp					(Quat &r, Quat const &x)				{r=exp(x);}
	
	void  r_r_fib					(Real &r, Real const &x)				{r=(exp(x*m_ln_phi)-cos(m_pi*x)*exp(-x*m_ln_phi))*m_inv_sqrt5;}
	void  c_c_fib					(Comp &r, Comp const &x)				{r=(exp(x*m_ln_phi)-cos(m_pi*x)*exp(-x*m_ln_phi))*m_inv_sqrt5;}
	void  q_q_fib					(Quat &r, Quat const &x)				{r=(exp(x*m_ln_phi)-cos(m_pi*x)*exp(-x*m_ln_phi))*m_inv_sqrt5;}
	
	void  r_r_random				(Real &r, Real const &x)				{r=mpfr::random();}
	void  c_c_random				(Comp &r, Comp const &x)				{r=Comp(mpfr::random(), mpfr::random());}
	void  q_q_random				(Quat &r, Quat const &x)				{r=Quat(mpfr::random(), mpfr::random(), mpfr::random(), mpfr::random());}
	void r_rr_random				(Real &r, Real const &x, Real const &y)	{r=mpfr::random();}
	void c_cr_random				(Comp &r, Comp const &x, Real const &y)	{r=Comp(mpfr::random(), mpfr::random());}
	void c_cc_random				(Comp &r, Comp const &x, Comp const &y)	{r=Comp(mpfr::random(), mpfr::random());}
	void q_qq_random				(Quat &r, Quat const &x, Quat const &y)	{r=Quat(mpfr::random(), mpfr::random(), mpfr::random(), mpfr::random());}

	inline Real beta(Real x, Real y){return exp(lgamma(x)+lgamma(y)-lgamma(x+y));}
	void  r_r_beta					(Real &r, Real const &x)				{r=beta(x, x);}
	void r_rr_beta					(Real &r, Real const &x, Real const &y)	{r=beta(x, y);}
	
	//inline Real cyl_bessel_j(Real x, Real y)
	//{
	//	try
	//	{
	//		return boost::math::cyl_bessel_j(x, y);
	//	}
	//	catch(...){return _qnan;}
	//}
	void  r_r_cyl_bessel_j			(Real &r, Real const &x)				{r=besselj0(x);}
	void r_rr_cyl_bessel_j			(Real &r, Real const &x, Real const &y)	{r=besseljn(x.toLong(), y);}

	//inline Real cyl_neumann(Real x, Real y)
	//{
	//	if(x>0)
	//		int LOL_1=0;
	//	if(y!=y||y<0)
	//		return _qnan;
	//	try
	//	{
	//		return boost::math::cyl_neumann(x, y);
	//	}
	//	catch(std::domain_error&){return _qnan;}
	//	catch(std::overflow_error&){return -_HUGE;}
	//	catch(...){return _qnan;}
	//}
	void  r_r_cyl_neumann			(Real &r, Real const &x)				{r=bessely0(x);}
	void r_rr_cyl_neumann			(Real &r, Real const &x, Real const &y)	{r=besselyn(x.toLong(), y);}

	inline Comp r_hankel1(Real x, Real y)
	{
		return besseljn(x.toLong(), y)+m_i*besselyn(x.toLong(), y);
		//try
		//{
		//	return boost::math::cyl_bessel_j(x, y)+Comp(0, 1)*boost::math::cyl_neumann(x, y);
		//}
		//catch(std::domain_error&){return _qnan;}
		//catch(std::overflow_error&){return -_HUGE;}
		//catch(...){return _qnan;}
	}
	void  c_r_hankel1				(Comp &r, Real const &x)				{r=r_hankel1(x, x);}
	void  c_c_hankel1				(Comp &r, Comp const &x)				{r=r_hankel1(x.r, x.i);}
	void c_rr_hankel1				(Comp &r, Real const &x, Real const &y)	{r=r_hankel1(x, y);}
	
//	inline Real	sgn			(Real const &x){return (x>0)-(x<0);}
	inline Comp	sgn			(Comp const &x)
	{
		Real temp=abs(x);
		return temp?x/temp:Comp();
	}
	inline Quat	sgn			(Quat	const &x)
	{
		Real temp=abs(x);
		return temp!=0?x/temp:Quat();
	}
	inline Real	step		(Real const &x){return 0.5+0.5*sgn(x);}
	inline Comp	step		(Comp const &x){return 0.5+0.5*sgn(x);}
	inline Quat	step		(Quat const &x){return 0.5+0.5*sgn(x);}
	void  r_r_step					(Real &r, Real const &x)				{r=step(x);}
	void  c_c_step					(Comp &r, Comp const &x)				{r=step(x);}
	void  q_q_step					(Quat &r, Quat const &x)				{r=step(x);}

	void  r_r_rect					(Real &r, Real const &x)				{r=step(x+0.5)-step(x-0.5);}
	void  c_c_rect					(Comp &r, Comp const &x)				{Comp cx=x; r=step(cx+0.5)-step(cx-0.5);}
	void  q_q_rect					(Quat &r, Quat const &x)				{Quat qx=x; r=step(qx+0.5)-step(qx-0.5);}

	void  r_r_trgl					(Real &r, Real const &x)				{Real t=abs(x);	r=t<1?1-t:0;}
	void  r_c_trgl					(Real &r, Comp const &x)				{Real t=abs(x);	r=t<1?1-t:0;}
	void  r_q_trgl					(Real &r, Quat const &x)				{Real t=abs(x);	r=t<1?1-t:0;}

	void  r_r_sqwv					(Real &r, Real const &x)				{r=x  -floor(x  )<0.5;}
	void  r_c_sqwv					(Real &r, Comp const &x)				{r=x.r-floor(x.r)<0.5;}
	void  r_q_sqwv					(Real &r, Quat const &x)				{r=x.r-floor(x.r)<0.5;}
	void r_rr_sqwv					(Real &r, Real const &x, Real const &y)	{r=x  -floor(x  )<y  ;}
	void r_rc_sqwv					(Real &r, Real const &x, Comp const &y)	{r=x  -floor(x  )<y.r;}
	void r_rq_sqwv					(Real &r, Real const &x, Quat const &y)	{r=x  -floor(x  )<y.r;}
	void r_cr_sqwv					(Real &r, Comp const &x, Real const &y)	{r=x.r-floor(x.r)<y  ;}
	void r_cc_sqwv					(Real &r, Comp const &x, Comp const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_cq_sqwv					(Real &r, Comp const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_qr_sqwv					(Real &r, Quat const &x, Real const &y)	{r=x.r-floor(x.r)<y  ;}
	void r_qc_sqwv					(Real &r, Quat const &x, Comp const &y)	{r=x.r-floor(x.r)<y.r;}
	void r_qq_sqwv					(Real &r, Quat const &x, Quat const &y)	{r=x.r-floor(x.r)<y.r;}

	Real clamp_positive(Real const &x){return signbit(x)?x:0;}
	void  r_r_trwv					(Real &r, Real const &x)				{auto temp=abs(x-floor(x)-0.5); r=temp+temp;}
	void  r_c_trwv					(Real &r, Comp const &x)				{auto temp=abs(x-floor(x)-0.5); r=temp+temp;}
	void  r_q_trwv					(Real &r, Quat const &x)				{Quat qx=x; r=2*abs(qx-floor(qx)-0.5);}
	void r_rr_trwv					(Real &r, Real const &x, Real const &y)
	{
		Real t=x-floor(x), t2=1-x;
		t2-=floor(t2);
		Real dc=clamp_positive(y);
		dc=dc>1?1:dc;
		Real dc2=1-dc, t_d=t/dc, t2_d2=t2/dc2;
		r=(t_d<1)*t_d+(t2_d2<1)*t2_d2;
	}
	void c_cr_trwv					(Comp &r, Comp const &x, Real const &y)
	{
		Comp cx=x, t=cx-floor(cx), t2=1.-cx;
		t2-=floor(t2);
		Real dc=clamp_positive(y);
		dc=dc>1?1:dc;
		Real dc2=1-dc;
		auto t_d=t/dc, t2_d2=t2/dc2;
		r=Real(t_d.r<1)*t_d+Real(t2_d2.r<1)*t2_d2;
	}
	void c_cc_trwv					(Comp &r, Comp const &x, Comp const &y)
	{
		Comp cx=x, cy=y, t=cx-floor(cx), t2=1.-cx;
		t2-=floor(t2);
		Comp dc=cy;
		dc.r=clamp_positive(dc.r);
		dc.r=dc.r>1?1:dc.r;
		Comp dc2=1.-dc;
		auto t_d=t/dc, t2_d2=t2/dc2;
		r=Real(t_d.r<1)*t_d+Real(t2_d2.r<1)*t2_d2;
	}
	void q_qq_trwv					(Quat &r, Quat const &x, Quat const &y)
	{
		Quat cx=x, cy=y, t=cx-floor(cx), t2=1.-cx;
		t2-=floor(t2);
		Quat dc=cy;
		dc.r=clamp_positive(dc.r);
		dc.r=dc.r>1?1:dc.r;
		Quat dc2=1.-dc;
		auto t_d=t/dc, t2_d2=t2/dc2;
		r=Real(t_d.r<1)*t_d+Real(t2_d2.r<1)*t2_d2;
	}

	void  r_r_saw					(Real &r, Real const &x){Real t=x-floor(x), t2=floor(1-t); r=(t2+1)*(t2*0.5+t);}
	void  c_c_saw					(Comp &r, Comp const &x){Comp t=x-floor(x), t2=floor(1-t); r=(t2+1)*(t2*0.5+t);}
	void  q_q_saw					(Quat &r, Quat const &x){Quat t=x-floor(x), t2=floor(1-t); r=(t2+1)*(t2*0.5+t);}
	void r_rr_saw					(Real &r, Real const &x, Real const &y)
	{
		auto t=x-floor(x), t2=floor(y-t);
		r=(t2+1)*(t2*0.5+t)/y;
	}
	void c_rc_saw					(Comp &r, Real const &x, Comp const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void q_rq_saw					(Quat &r, Real const &x, Quat const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void c_cr_saw					(Comp &r, Comp const &x, Real const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void c_cc_saw					(Comp &r, Comp const &x, Comp const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void q_cq_saw					(Quat &r, Comp const &x, Quat const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void q_qr_saw					(Quat &r, Quat const &x, Real const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void q_qc_saw					(Quat &r, Quat const &x, Comp const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}
	void q_qq_saw					(Quat &r, Quat const &x, Quat const &y)
	{
		auto t=x-floor(x);
		auto t2=floor(y-t);
		r=(t2+1.)*(t2*0.5+t)/y;
	}

	void r_rr_hypot					(Real &r, Real const &x, Real const &y)	{r=sqrt(x*x+y*y);}
	//void c_cc_hypot					(Comp const &x, Comp const &y)	{return sqrt(sq(x)+sq(y));}
	//void q_qq_hypot					(Quat const &x, Quat const &y)	{return sqrt(sq(x)+sq(y));}

	inline int mandelbrot(Comp const &x, int n_iterations)
	{
		Real rez=0, imz=0, sq_rez=0, sq_imz=0;
		int k=0;
		for(;k<n_iterations&&sq_rez+sq_imz<16;++k)
		{
			imz=rez*imz;//calculate sq(z)
			imz+=imz;
			rez=sq_rez-sq_imz;

			rez+=x.r, imz+=x.i;//add x

			sq_rez=rez*rez, sq_imz=imz*imz;
		}
		return k;
	}
	void r_r_mandelbrot				(Real &r, Real const &x)				{r=mandelbrot(Comp(x), 200);}
	void r_c_mandelbrot				(Real &r, Comp const &x)				{r=mandelbrot(x, 200);}
	void r_rr_mandelbrot			(Real &r, Real const &x, Real const &y)	{r=mandelbrot(Comp(x), y.toLong());}
	void r_cr_mandelbrot			(Real &r, Comp const &x, Real const &y)	{r=mandelbrot(x, y.toLong());}

	void r_rr_min					(Real &r, Real const &x, Real const &y)	{r=(x+y-abs(x-y))*0.5;}
	void c_cr_min					(Comp &r, Comp const &x, Real const &y)	{r=(x+y-abs(x-y))*0.5;}
	void c_cc_min					(Comp &r, Comp const &x, Comp const &y)	{r=(x+y-abs(x-y))*0.5;}
	void q_qq_min					(Quat &r, Quat const &x, Quat const &y)	{r=(x+y-abs(x-y))*0.5;}

	void r_rr_max					(Real &r, Real const &x, Real const &y)	{r=(x+y+abs(x-y))*0.5;}
	void c_cr_max					(Comp &r, Comp const &x, Real const &y)	{r=(x+y+abs(x-y))*0.5;}
	void c_cc_max					(Comp &r, Comp const &x, Comp const &y)	{r=(x+y+abs(x-y))*0.5;}
	void q_qq_max					(Quat &r, Quat const &x, Quat const &y)	{r=(x+y+abs(x-y))*0.5;}

	void r_rr_conditional_110		(Real &r, Real const &x, Real const &y)	{r=x?y:0;}
	void c_rc_conditional_110		(Comp &r, Real const &x, Comp const &y)	{r=x?y:Comp();}
	void q_rq_conditional_110		(Quat &r, Real const &x, Quat const &y)	{r=x?y:Quat();}
	void r_cr_conditional_110		(Real &r, Comp const &x, Real const &y)	{r=x.c_is_true()?y:0;}
	void c_cc_conditional_110		(Comp &r, Comp const &x, Comp const &y)	{r=x.c_is_true()?y:Comp();}
	void q_cq_conditional_110		(Quat &r, Comp const &x, Quat const &y)	{r=x.c_is_true()?y:Quat();}
	void r_qr_conditional_110		(Real &r, Quat const &x, Real const &y)	{r=x.q_is_true()?y:0;}
	void c_qc_conditional_110		(Comp &r, Quat const &x, Comp const &y)	{r=x.q_is_true()?y:Comp();}
	void q_qq_conditional_110		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?y:Quat();}
	
	void r_rr_conditional_101		(Real &r, Real const &x, Real const &y)	{r=x?0:y;}
	void c_rc_conditional_101		(Comp &r, Real const &x, Comp const &y)	{r=x?Comp():y;}
	void q_rq_conditional_101		(Quat &r, Real const &x, Quat const &y)	{r=x?Quat():y;}
	void r_cr_conditional_101		(Real &r, Comp const &x, Real const &y)	{r=x.c_is_true()?0:y;}
	void c_cc_conditional_101		(Comp &r, Comp const &x, Comp const &y)	{r=x.c_is_true()?Comp():y;}
	void q_cq_conditional_101		(Quat &r, Comp const &x, Quat const &y)	{r=x.c_is_true()?Quat():y;}
	void r_qr_conditional_101		(Real &r, Quat const &x, Real const &y)	{r=x.q_is_true()?0:y;}
	void c_qc_conditional_101		(Comp &r, Quat const &x, Comp const &y)	{r=x.q_is_true()?Comp():y;}
	void q_qq_conditional_101		(Quat &r, Quat const &x, Quat const &y)	{r=x.q_is_true()?Quat():y;}

	void  r_r_increment				(Real &r, Real const &x)				{r=x+1;}
	void  c_c_increment				(Comp &r, Comp const &x)				{r=x+1.;}
	void  q_q_increment				(Quat &r, Quat const &x)				{r=x+1.;}

	void  r_r_decrement				(Real &r, Real const &x)				{r=x-1;}
	void  c_c_decrement				(Comp &r, Comp const &x)				{r=x-1.;}
	void  q_q_decrement				(Quat &r, Quat const &x)				{r=x-1.;}

	void  r_r_assign				(Real &r, Real const &x)				{r=x;}
	void  c_c_assign				(Comp &r, Comp const &x)				{r=x;}
	void  q_q_assign				(Quat &r, Quat const &x)				{r=x;}
#endif
}