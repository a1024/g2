G2_R_R(sin){IDX; ASSIGN_R(sin(xr[idx]));}
G2_C_C(sin){IDX; float2 ret=sin_c(VEC2(x)); RET_C;}
G2_Q_Q(sin){IDX; float4 ret=sin_q(VEC4(x)); RET_Q;}

G2_C_C(asin){IDX; float2 ret=asin_c(VEC2(x)); RET_C;}
G2_Q_Q(asin){IDX; float4 ret=asin_q(VEC4(x)); RET_Q;}
DISC_C_I(asin)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.y==x0.y)
		disc[idx]=false;
	else if(x0.x==x0.x)
		disc[idx]=(x0.y<=0?x1.y>0:x1.y<=0)&&(x0.x<-1||x0.x>1);
	else if((x0.y<=0&&x1.y>0)||(x1.y<=0&&x0.y>0))
	{
		float t=_1d_zero_crossing(x0.x, x0.y, x1.x, x1.y);
		disc[idx]=t<-1||t>1;
	}
	else
		disc[idx]=false;
}
DISC_Q_I(asin){IDX; disc[idx]=false;}//TODO

G2_R_R(sinh){IDX; ASSIGN_R(sinh(xr[idx]));}
G2_C_C(sinh){IDX; float2 ret=sinh_c(VEC2(x)); RET_C;}
G2_Q_Q(sinh){IDX; float4 ret=sinh_q(VEC4(x)); RET_Q;}

G2_R_R(asinh){IDX; ASSIGN_R(asinh(xr[idx]));}
G2_C_C(asinh){IDX; float2 ret=asinh_c(VEC2(x)); RET_C;}
G2_Q_Q(asinh){IDX; float4 ret=asinh_q(VEC4(x)); RET_Q;}
DISC_C_I(asinh)
{
	IDX;
	float2 x0=(float2)(xi[idx], xr[idx]), x1=(float2)(xi[idx+offset], xr[idx+offset]);//sic
	if(x0.y==x0.y)
		disc[idx]=false;
	else if(x0.x==x0.x)
		disc[idx]=(x0.y<=0?x1.y>0:x1.y<=0)&&(x0.x<-1||x0.x>1);
	else if((x0.y<=0&&x1.y>0)||(x1.y<=0&&x0.y>0))
	{
		float t=_1d_zero_crossing(x0.x, x0.y, x1.x, x1.y);
		disc[idx]=t<-1||t>1;
	}
	else
		disc[idx]=false;
}
DISC_Q_I(asinh){IDX; disc[idx]=false;}//TODO

G2_R_R(sinc){IDX; float a=xr[idx]; ASSIGN_R(a!=0?sin(a)/a:1);}
G2_C_C(sinc){IDX; float2 a=VEC2(x), ret=istrue_c(a)?div_cc(sin_c(a), a):(float2)(1, 0); RET_C;}
G2_Q_Q(sinc){IDX; float4 a=VEC4(x), ret=istrue_q(a)?div_qq(sin_q(a), a):(float4)(1, 0, 0, 0); RET_Q;}

G2_R_R(sinhc){IDX; float a=xr[idx]; ASSIGN_R(a!=0?sinh(a)/a:1);}
G2_C_C(sinhc){IDX; float2 a=VEC2(x), ret=istrue_c(a)?div_cc(sinh_c(a), a):(float2)(1, 0); RET_C;}
G2_Q_Q(sinhc){IDX; float4 a=VEC4(x), ret=istrue_q(a)?div_qq(sinh_q(a), a):(float4)(1, 0, 0, 0); RET_Q;}

G2_R_R(csc){IDX; ASSIGN_R(1/sin(xr[idx]));}
G2_C_C(csc){IDX; float2 ret=inv_c(sin_c(VEC2(x))); RET_C;}
G2_Q_Q(csc){IDX; float4 ret=inv_q(sin_q(VEC4(x))); RET_Q;}
DISC_R_I(csc)
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset];
	if(fabs(x1-x0)>3.2f)
		disc[idx]=true;
	else
		disc[idx]=_1d_int_in_range(x0/_pi, x1/_pi);
}
DISC_C_I(csc)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.x==x1.x)
		disc[idx]=true;
	else if(x0.y==x1.y)
		disc[idx]=x0.y==0&&_1d_int_in_range(x0.x/_pi, x1.x/_pi);
	else if(signbit(x0.y)!=signbit(x1.y))
	{
		float t=_1d_zero_crossing(x0.x, x0.y, x1.x, x0.y)/_pi;
		disc[idx]=t==floor(t);
	}
	else
		disc[idx]=false;
}
DISC_Q_I(csc){IDX; disc[idx]=false;}//TODO

G2_C_C(acsc){IDX; float2 ret=asin_c(inv_c(VEC2(x))); RET_C;}
G2_Q_Q(acsc){IDX; float4 ret=asin_q(inv_q(VEC4(x))); RET_Q;}
DISC_C_I(acsc)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.y==x1.y)
		disc[idx]=x0.y==0&&(x0.x<0?x1.x>=0:x0.x>0?x1.x<=0:x1.x!=0);//x1.x<0||x1.x>0);
	else if(x0.x==x1.x)
	{
		if(x0.x<0)
			disc[idx]=x0.x>-1&&(x0.y<=0?x1.y>0:x1.y<=0);
		else if(x0.x==0)
			disc[idx]=x0.y<0?x1.y>=0:x0.y==0?x1.y<0||x1.y>0:x1.y<=0;
		else
			disc[idx]=x0.x<1&&(x0.y<0?x1.y>=0:x1.y<0);
	}
	else
		disc[idx]=false;
}
DISC_Q_I(acsc){IDX; disc[idx]=false;}//TODO

G2_R_R(csch){IDX; ASSIGN_R(1/sinh(xr[idx]));}
G2_C_C(csch){IDX; float2 ret=inv_c(sinh_c(VEC2(x))); RET_C;}
G2_Q_Q(csch){IDX; float4 ret=inv_q(sinh_q(VEC4(x))); RET_Q;}
DISC_R_I(csch){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}
DISC_C_I(csch){disc_c_csc_i(size, offset, disc, xi, xr);}//sic
DISC_Q_I(csch){IDX; disc[idx]=false;}//TODO

G2_R_R(acsch){IDX; ASSIGN_R(asinh(1/xr[idx]));}
G2_C_C(acsch){IDX; float2 ret=asinh_c(inv_c(VEC2(x))); RET_C;}
G2_Q_Q(acsch){IDX; float4 ret=asinh_q(inv_q(VEC4(x))); RET_Q;}
DISC_R_I(acsch){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}
DISC_C_I(acsch){disc_c_acsc_i(size, offset, disc, xi, xr);}//sic
DISC_Q_I(acsch){IDX; disc[idx]=false;}//TODO