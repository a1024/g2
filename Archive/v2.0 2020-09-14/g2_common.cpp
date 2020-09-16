//g2_common.cpp - Implementation of common functions used by all Grapher 2 source files.
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

#include"g2_common.h"
namespace	G2
{
	const double	_atm	=101325,			_bbr	=5.670373e-8,
					_c		=299792458,			_ele	=8.854187817620e-12,
					_e		=::exp(1.),			_g		=9.80665,
					_G		=6.67384e-11,		_h		=6.62606957e-34,
					_mag	=1.2566370614e-6,	_me		=9.10938291e-31,
					_mn		=1.674927351e-27,	_mp		=1.672621777e-27,
					_Me		=5.9736e24,			_Ms		=1.9891e30,
					_Na		=6.02214129e23,		_phi	=1.6180339887498948482045868343656381177203091798057628621354486227052604628189,
					_pi		=::acos(-1.),		_q		=1.602176565e-19,
					_R		=8.3144621,			_qnan	=std::numeric_limits<double>::quiet_NaN(),

					_ln2	=::log(2.),			_ln10	=::log(10.),	inv_ln10=1/_ln10,	_sqrt2=::sqrt(2.),	_sqrt3=::sqrt(3.),	_sqrt5=::sqrt(5.),
					_pi_2=_pi*0.5,	_2pi=2*_pi,	_sqrt_2pi=::sqrt(_2pi),	_ln_sqrt_2pi=::log(_sqrt_2pi),	_1_2pi=1/_2pi, _1_pi=1/_pi, _pi_180=_pi/180;
#pragma			warning(push)
#pragma			warning(disable:4723)
	const float		s_zero=0, s_infinity=1/s_zero;
//	const double	zero=0, infinity=1/zero;
#pragma			warning(pop)
}
#if _MSC_VER<1700//__cplusplus<201103
namespace	std
{
	double lgamma(double x)//https://jamesmccaffrey.wordpress.com/2013/06/19/the-log-gamma-function-with-c/
	{
		const double coef[6]={76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
		double denom=x+1, y=x+5.5, series=1.000000000190015;
		for (int i=0;i<6;++i)
		{
			series+=coef[i]/denom;
			denom+=1;
		}
		return G2::_ln_sqrt_2pi+(x+0.5)*::log(y)-y+::log(series/x);
	}
}
#endif
namespace	G2
{
	bool _2d_between(double x1, double y1, double x, double y, double x2, double y2)
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
	const double ll_max=9.22337203685478e+018;
	Quat1d	log					(Quat1d	const &x)
	{
		double t=x.R_component_2()*x.R_component_2()+x.R_component_3()*x.R_component_3()+x.R_component_4()*x.R_component_4();
		if(t)//quat
		{
			double absx=::sqrt(x.R_component_1()*x.R_component_1()+t);
			t=::acos(x.R_component_1()/absx)/::sqrt(t);
			return Quat1d(::log(absx), t*x.R_component_2(), t*x.R_component_3(), t*x.R_component_4());
		}
		return Quat1d(::log(::abs(x.R_component_1())), _pi*(x.R_component_1()<0));

		//try
		//{
		//	double t=boost::math::abs(Quat1d(0, x.R_component_2(), x.R_component_3(), x.R_component_4()));
		//	if(t)
		//	{
		//		t=arg(x)/t;
		//		return Quat1d(boost::math::log1p(abs(x)-1), t*x.R_component_2(), t*x.R_component_3(), t*x.R_component_4());
		//	}
		//	if(x.R_component_1()!=0)
		//		return Quat1d(boost::math::log1p(abs(x)-1), 0, 0, 0);
		//	return Quat1d(-_HUGE);
		//}
		//catch(std::overflow_error&)
		//{
		//	return Quat1d(-_HUGE);
		//}
	}
	namespace gamma//http://en.wikipedia.org/wiki/Lanczos_approximation
	{
		const double g=7, p[]={0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	}
	double			tgamma		(double const &x)
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
	Comp1d			tgamma		(Comp1d const &x)
	{
		using namespace gamma;
		if(x.real()<.5)
		{
			Comp1d t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			Comp1d t2=g+.5-x;
			return _pi/(sin(_pi*x)*_sqrt_2pi*std::pow(t2, .5-x)*std::exp(-t2)*t1);
		}
		else
		{
			Comp1d t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			Comp1d t2=g+.5+x-1.;
			return _sqrt_2pi*std::pow(t2, .5+x-1.)*std::exp(-t2)*t1;
		}
		//if(x.real()<.5)
		//	return _pi/(sin(_pi*x)*tgamma(1.-x));
		//else
		//{
		//	Comp1d t1(p[0]);
		//	for(int k=1;k<g+2;++k)
		//		t1+=p[k]/(x-1.+double(k));
		//	Comp1d t2=x-1.+g+.5;
		//	return _sqrt_2pi*std::pow(t2, x-1.+.5)*std::exp(-t2)*t1;
		//}
	}
	Quat1d			tgamma		(Quat1d	const &x)
	{
		using namespace gamma;
		if(x.real()<.5)
		{
			Quat1d t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)-x);
			Quat1d t2=g+.5-x;
			return _pi/(sin(_pi*x)*_sqrt_2pi*pow(t2, .5-x)*exp(-t2)*t1);
		}
		else
		{
			Quat1d t1(p[0]);
			for(int k=1;k<g+2;++k)
				t1+=p[k]/(double(k)+x-1.);
			Quat1d t2=g+.5+x-1.;
			return _sqrt_2pi*pow(t2, .5+x-1.)*exp(-t2)*t1;
		}
	}
}
#if 0
namespace	modes
{
	inline unsigned char f(int n, double a, double H, double L)
	{
		double k=n+H/30;
		k-=floor(k/12)*12;
		return unsigned char(255*(L-a*maximum(-1, minimum(k-3, 9-k, 1))));
	}
	inline int hsl2rgb(double H, double S, double L)
	{
		double C=(1-abs(L+L-1))*S;
		double H2=1.f/60*H;
		double X=C*(1-abs(H2-floor(H2*0.5)*2-1));
		double m=L-C*0.5;
		double R, G, B;
			 if(0<=H&&H<60)		R=C, G=X, B=0;
		else if(60<=H&&H<120)	R=X, G=C, B=0;
		else if(120<=H&&H<180)	R=0, G=C, B=X;
		else if(180<=H&&H<240)	R=0, G=X, B=C;
		else if(240<=H&&H<300)	R=X, G=0, B=C;
		else if(300<=H&&H<360)	R=C, G=0, B=X;
		return unsigned char(255*(R+m))<<16|unsigned char(255*(G+m))<<8|unsigned char(255*(B+m));

		//double a=S*(0.5-abs(L-0.5));//S*minimum(L, 1-L);
		//unsigned char r=f(0, a, H, L), g=f(8, a, H, L), b=f(4, a, H, L);
		//return r<<16|g<<8|b;
	}
	int colorFunction(double &r, double &i)
	{
	/*	const double cos_pi_6=0.866025403784439, sin_pi_6=0.5;	//black->color simple loop		56.42 cyclec/px
		double hyp=sqrt(r*r+i*i), cosx=r/hyp, sinx=i/hyp,
			mag=12.75*(hyp-floor(hyp*0.1)*10);
		double red=1+cosx*cos_pi_6-sinx*sin_pi_6, green=1+sinx, blue=1+cosx*-cos_pi_6-sinx*sin_pi_6;
		red*=mag, green*=mag, blue*=mag;
		int color=unsigned char(red)<<16|unsigned char(green)<<8|unsigned char(blue);//argb
		return color;//*/

		if(r!=r||i!=i)										//81.02 cycles/px
			return 0x7F7F7F;
		if(abs(r)==_HUGE||abs(i)==_HUGE)
			return 0xFFFFFF;
		double threshold=10, inv_th=1/threshold;//1			//black->color->white	67.52 cycles/px
		const double cos_pi_6=0.866025403784439, sin_pi_6=0.5;
		double hyp=sqrt(r*r+i*i), cosx=r/hyp, sinx=i/hyp,
			mag=255*exp(-hyp*G2::_ln2*inv_th);
		double red=1+cosx*cos_pi_6-sinx*sin_pi_6, green=1+sinx, blue=1+cosx*-cos_pi_6-sinx*sin_pi_6;
		if(hyp<threshold)
			mag=255-mag, red*=mag, green*=mag, blue*=mag;
		else
			red=255-mag*(2-red), green=255-mag*(2-green), blue=255-mag*(2-blue);
		int color=unsigned char(red)<<16|unsigned char(green)<<8|unsigned char(blue);//argb
		return color;//*/

	/*	if(r!=r||i!=i)
			return 0x7F7F7F;
		if(abs(r)==_HUGE||abs(i)==_HUGE)
			return 0xFFFFFF;
		const double cos_pi_6=0.866025403784439, sin_pi_6=0.5;//black->color whitening loop		87.29 cycles/px		//89.22 cycles/px
		double hyp=sqrt(r*r+i*i), cosx=r/hyp, sinx=i/hyp,
			fh=floor(hyp*0.1), c, f;
		c=fh*0.01;
		c=255*(c-floor(c)), f=12.75*(hyp-fh*10)*(1-c*1./255);
		double red=1+cosx*cos_pi_6-sinx*sin_pi_6, green=1+sinx, blue=1+cosx*-cos_pi_6-sinx*sin_pi_6;
	//	red*=mag, green*=mag, blue*=mag;
		red=c+f*red, green=c+f*green, blue=c+f*blue;
	//	red=255-fh*(255-red), green=255-fh*(255-green), blue=255-fh*(255-blue);
		int color=unsigned char(red)<<16|unsigned char(green)<<8|unsigned char(blue);//argb
		return color;//*/
		
	/*	if(r!=r||i!=i)									//95 cyclec/px
			return 0x7F7F7F;
		if(r==_HUGE)
		{
			if(i==_HUGE||i==-_HUGE)
				return 0x00FFFFFF;
			else
				return 0x00ED7F11;//arg=0
		}
		else if(r==-_HUGE)
		{
			if(i==_HUGE||i==-_HUGE)
				return 0x00FFFFFF;
			else
				return 0x00117FED;//arg=pi
		}
		else if(i==_HUGE)
			return 0x003FFF3F;//arg=pi/2
		else if(i==-_HUGE)
			return 0x00BF00BF;//arg=-pi/2
		const double cos_pi_6=0.866025403784439, sin_pi_6=0.5;
		double hyp=sqrt(r*r+i*i), mag=255/G2::_pi*atan(hyp), cosx=r/hyp, sinx=i/hyp;
		double red=mag*(1+cosx*cos_pi_6-sinx*sin_pi_6), green=mag*(1+sinx), blue=mag*(1+cosx*-cos_pi_6-sinx*sin_pi_6);
		int color=unsigned char(mag*(1+cosx*cos_pi_6-sinx*sin_pi_6))<<16|unsigned char(mag*(1+sinx))<<8|unsigned char(mag*(1+cosx*-cos_pi_6-sinx*sin_pi_6));//argb
		return color;//*/
	//	return unsigned char(mag*(1+cosx*cos_pi_6-sinx*sin_pi_6))<<16;

	/*	const double lnbase=log(0.75);				//683.35 cycles/px
		double H=180/G2::_pi*atan2(i, r), S=1, L=1-exp(sqrt(r*r+i*i)*lnbase);
	//	double H=r, S=1, L=i;
		return hsl2rgb(H, S, L);//*/
	}
	void colorFunction_cases(double r, double i, int &c)
	{
		if(r!=r||i!=i)
			c=0x7F7F7F;
		if(r==_HUGE)
		{
			if(i==_HUGE||i==-_HUGE)
				c=0x00FFFFFF;
			else
				c=0x00ED7F11;//arg=0
		}
		else if(r==-_HUGE)
		{
			if(i==_HUGE||i==-_HUGE)
				c=0x00FFFFFF;
			else
				c=0x00117FED;//arg=pi
		}
		else if(i==_HUGE)
			c=0x003FFF3F;//arg=pi/2
		else if(i==-_HUGE)
			c=0x00BF00BF;//arg=-pi/2
	}
}
#endif