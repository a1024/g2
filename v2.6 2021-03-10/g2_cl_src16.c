G2_R_R(tan){IDX; ASSIGN_R(tan(xr[idx]));}
G2_C_C(tan){IDX; float2 ret=tan_c(VEC2(x)); RET_C;}
G2_Q_Q(tan){IDX; float4 ret=tan_q(VEC4(x)); RET_Q;}
DISC_R_I(tan)
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset];
	disc[idx]=fabs(x1-x0)>3.2f||_1d_int_in_range(x0/_pi-0.5f, x1/_pi-0.5f);
}
DISC_C_I(tan)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.x==x1.x)
	{
		if(x0.y==x1.y)
			disc[idx]=false;
		else
		{
			float t=x0.x/_pi-0.5f;
			disc[idx]=t==floor(t);
		}
	}
	else if(x0.y==x1.y)
		disc[idx]=x0.y==0&&_1d_int_in_range(x0.x/_pi-0.5f, x1.x/_pi-0.5f);
	if(sign(x0.y)!=sign(x1.y))
	{
		float t=_1d_zero_crossing(x0.x, x0.y, x1.x, x0.y)/_pi-0.5f;
		disc[idx]=t==floor(t);
	}
	disc[idx]=false;
}
DISC_Q_I(tan){IDX; disc[idx]=false;}//TODO

float atan_addition(float x, float y){return x<0?y<0?-_pi:_pi:0;}
G2_R_R(atan){IDX; ASSIGN_R(atan(xr[idx]));}
G2_C_C(atan){IDX; float2 ret=atan_c(VEC2(x)); RET_C;}
G2_Q_Q(atan){IDX; float4 ret=atan_q(VEC4(x)); RET_Q;}
G2_R_RR(atan){IDX; ASSIGN_R(atan2(yr[idx], xr[idx]));}
G2_C_RC(atan)
{
	IDX;
	float a=xr[idx];
	float2 b=VEC2(y);
	float2 ret=atan_c(div_rc(a, b))+(float2)(atan_addition(a, b.x), 0);
	RET_C;
}
G2_Q_RQ(atan)
{
	IDX;
	float a=xr[idx];
	float4 b=VEC4(y);
	float4 ret=atan_q(div_rq(a, b))+(float4)(atan_addition(a, b.x), 0, 0, 0);
	RET_Q;
}
G2_C_CR(atan)
{
	IDX;
	float2 a=VEC2(x);
	float b=yr[idx];
	float2 ret=atan_c(div_cr(a, b))+(float2)(atan_addition(a.x, b), 0);
	RET_C;
}
G2_C_CC(atan)
{
	IDX;
	float2 a=VEC2(x);
	float2 b=VEC2(y);
	float2 ret=atan_c(div_cc(a, b))+(float2)(atan_addition(a.x, b.x), 0);
	RET_C;
}
G2_Q_CQ(atan)
{
	IDX;
	float2 a=VEC2(x);
	float4 b=VEC4(y);
	float4 ret=atan_q(div_cq(a, b))+(float4)(atan_addition(a.x, b.x), 0, 0, 0);
	RET_Q;
}
G2_Q_QR(atan)
{
	IDX;
	float4 a=VEC4(x);
	float b=yr[idx];
	float4 ret=atan_q(div_qr(a, b))+(float4)(atan_addition(a.x, b), 0, 0, 0);
	RET_Q;
}
G2_Q_QC(atan)
{
	IDX;
	float4 a=VEC4(x);
	float2 b=VEC2(x);
	float4 ret=atan_q(div_qc(a, b))+(float4)(atan_addition(a.x, b.x), 0, 0, 0);
	RET_Q;
}
G2_Q_QQ(atan)
{
	IDX;
	float4 a=VEC4(x);
	float4 b=VEC4(x);
	float4 ret=atan_q(div_qq(a, b))+(float4)(atan_addition(a.x, b.x), 0, 0, 0);
	RET_Q;
}
DISC_C_I(atan)
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
DISC_Q_I(atan){IDX; disc[idx]=false;}//TODO
DISC_RR_I(atan)
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset], y0=yr[idx], y1=yr[idx+offset];
	if(y0<0)
	{
		if(x0<0)
		{
				 if(y1<0)	disc[idx]=x1>=0;
			else if(y1>0)	disc[idx]=x1>0&&y0+(0-x0)*(y1-y0)/(x1-x0)<=0;
			else			disc[idx]=x1>=0;
		}
		else if(x0>0)
		{
				 if(y1<0)	disc[idx]=x1<0;
			else if(y1>0)	disc[idx]=x1<0&&y0+(0-y0)*(y1-y0)/(x1-x1)<=0;
			else			disc[idx]=x1<=0;
		}
		else
		{
				 if(y1<0)	disc[idx]=x1<0;
			else if(y1>0)	disc[idx]=x1<=0;
			else			disc[idx]=x1<=0;
		}
	}
	else if(y0>0)
	{
		if(x0<0)
		{
				 if(y1<0)	disc[idx]=x1>=0&&(x1==0||y0+(0-y0)*(y1-y0)/(x1-x1)<=0);
			else if(y1>0)	disc[idx]=false;
			else			disc[idx]=x1==0;
		}
		else if(x0>0)
		{
				 if(y1<0)	disc[idx]=x1<0&&y0+(0-y0)*(y1-y0)/(x1-x1)<=0;
			else if(y1>0)	disc[idx]=false;
			else			disc[idx]=x1==0;
		}
		else
		{
				 if(y1<0)	disc[idx]=x1==0;
			else if(y1>0)	disc[idx]=false;
			else			disc[idx]=x1==0;
		}
	}
	else
	{
		if(x0<0)
		{
				 if(y1<0)	disc[idx]=x1>=0;
			else if(y1>0)	disc[idx]=false;
			else			disc[idx]=x1>=0;
		}
		else if(x0>0)
		{
				 if(y1<0)	disc[idx]=x1<0;
			else if(y1>0)	disc[idx]=false;
			else			disc[idx]=x1<=0;
		}
		else				disc[idx]=true;
	}
}
DISC_RC_I(atan){IDX; disc[idx]=false;}//TODO
DISC_RQ_I(atan){IDX; disc[idx]=false;}
DISC_CR_I(atan){IDX; disc[idx]=false;}
DISC_CC_I(atan){IDX; disc[idx]=false;}
DISC_CQ_I(atan){IDX; disc[idx]=false;}
DISC_QR_I(atan){IDX; disc[idx]=false;}
DISC_QC_I(atan){IDX; disc[idx]=false;}
DISC_QQ_I(atan){IDX; disc[idx]=false;}

G2_R_R(tanh){IDX; ASSIGN_R(tanh(xr[idx]));}
G2_C_C(tanh){IDX; float2 ret=tanh_c(VEC2(x)); RET_C;}
G2_Q_Q(tanh){IDX; float4 ret=tanh_q(VEC4(x)); RET_Q;}

G2_C_C(atanh){IDX; float2 ret=atanh_c(VEC2(x)); RET_C;}
G2_Q_Q(atanh){IDX; float4 ret=atanh_q(VEC4(x)); RET_Q;}
DISC_C_I(atanh)
{
	IDX;
	float2 x0=(float2)(xi[idx], xr[idx]), x1=(float2)(xi[idx+offset], xr[idx+offset]);//sic
	if(x0.y==x0.y)
		disc[idx]=false;//disc c acos i
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
DISC_Q_I(atanh){IDX; disc[idx]=false;}

G2_R_R(tanc){IDX; float a=xr[idx]; ASSIGN_R(a!=0?tan(a)/a:0);}
G2_C_C(tanc){IDX; float2 a=VEC2(x), ret=istrue_c(a)?div_cc(tan_c(a), a):(float2)(0, 0); RET_C;}
G2_Q_Q(tanc){IDX; float4 a=VEC4(x), ret=istrue_q(a)?div_qq(tan_q(a), a):(float4)(0, 0, 0, 0); RET_Q;}
DISC_R_I(tanc)
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset];
	disc[idx]=fabs(x1-x0)>3.2f||_1d_int_in_range(x0/_pi-0.5f, x1/_pi-0.5f);
}
DISC_C_I(tanc)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.x==x1.x)
	{
		if(x0.y==x1.y)
			disc[idx]=false;
		else
		{
			float t=x0.x/_pi-0.5f;
			disc[idx]=t==floor(t);
		}
	}
	else if(x0.y==x1.y)
		disc[idx]=x0.y==0&&_1d_int_in_range(x0.x/_pi-0.5f, x1.x/_pi-0.5f);
	if(sign(x0.y)!=sign(x1.y))
	{
		float t=_1d_zero_crossing(x0.x, x0.y, x1.x, x0.y)/_pi-0.5f;
		disc[idx]=t==floor(t);
	}
	disc[idx]=false;
}
DISC_Q_I(tanc){IDX; disc[idx]=false;}//TODO

G2_R_R(cot){IDX; ASSIGN_R(1/tan(xr[idx]));}
G2_C_C(cot){IDX; float2 ret=inv_c(tan_c(VEC2(x))); RET_C;}
G2_Q_Q(cot){IDX; float4 ret=inv_q(tan_q(VEC4(x))); RET_Q;}
DISC_R_I(cot)
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset];
	if(fabs(x1-x0)>3.2f)
		disc[idx]=true;
	else
		disc[idx]=_1d_int_in_range(x0/_pi, x1/_pi);
}
DISC_C_I(cot)
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
DISC_Q_I(cot){IDX; disc[idx]=false;}//TODO

G2_R_R(acot){IDX; float a=xr[idx]; ASSIGN_R(a!=0?atan(1/a):_pi*0.5f);}
G2_C_C(acot){IDX; float2 a=VEC2(x), ret=istrue_c(a)?atan_c(inv_c(a)):(float2)(_pi*0.5f, 0); RET_C;}
G2_Q_Q(acot){IDX; float4 a=VEC4(x), ret=istrue_q(a)?atan_q(inv_q(a)):(float4)(_pi*0.5f, 0, 0, 0); RET_Q;}
DISC_R_I(acot){IDX; float x0=xr[idx], x1=xr[idx+offset]; disc[idx]=x0<0?x1>=0:x1<0;}
DISC_C_I(acot)
{
	IDX;
	float2 x0=(float2)(xi[idx], xr[idx]), x1=(float2)(xi[idx+offset], xr[idx+offset]);//sic
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
DISC_Q_I(acot){IDX; disc[idx]=false;}//TODO

G2_R_R(coth){IDX; ASSIGN_R(1/tanh(xr[idx]));}
G2_C_C(coth){IDX; float2 ret=inv_c(tanh_c(VEC2(x))); RET_C;}
G2_Q_Q(coth){IDX; float4 ret=inv_q(tanh_q(VEC4(x))); RET_Q;}
DISC_R_I(coth){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}//disc r csch i
DISC_C_I(coth)//disc c csch i
{
	IDX;
	float2 x0=(float2)(xi[idx], xr[idx]), x1=(float2)(xi[idx+offset], xr[idx+offset]);//sic
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
DISC_Q_I(coth){IDX; disc[idx]=false;}

G2_C_C(acoth){IDX; float2 ret=atanh_c(inv_c(VEC2(x))); RET_C;}
G2_Q_Q(acoth){IDX; float4 ret=atanh_q(inv_q(VEC4(x))); RET_Q;}
DISC_C_I(acoth)
{
	IDX;
	float2 x0=(float2)(xi[idx], xr[idx]), x1=(float2)(xi[idx+offset], xr[idx+offset]);//sic
	if(x0.y==x1.y)
		disc[idx]=x0.y==0&&(x0.x<0?x1.x>=0:x0.x>0?x1.x<=0:x1.x!=0);//x1.x<0||x1.x>0);//disc c acsc i
	else if(x0.x==x1.x)
	{
			 if(x0.x<0)		disc[idx]=x0.x>-1&&(x0.y<=0?x1.y>0:x1.y<=0);
		else if(x0.x==0)	disc[idx]=x0.y<0?x1.y>=0:x0.y==0?x1.y<0||x1.y>0:x1.y<=0;
		else				disc[idx]=x0.x<1&&(x0.y<0?x1.y>=0:x1.y<0);
	}
	else
		disc[idx]=false;
}
DISC_Q_I(acoth){IDX; disc[idx]=false;}//TODO