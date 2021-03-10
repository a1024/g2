G2_R_R(cos){IDX; ASSIGN_R(cos(xr[idx]));}
G2_C_C(cos){IDX; float2 ret=cos_c(VEC2(x)); RET_C;}
G2_Q_Q(cos){IDX; float4 ret=cos_q(VEC4(x)); RET_Q;}

G2_C_C(acos){IDX; float2 ret=acos_c(VEC2(x)); RET_C;}
G2_Q_Q(acos){IDX; float4 ret=acos_q(VEC4(x)); RET_Q;}
DISC_C_I(acos)
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
DISC_Q_I(acos){IDX; disc[idx]=false;}

G2_R_R(cosh){IDX; ASSIGN_R(cosh(xr[idx]));}
G2_C_C(cosh){IDX; float2 ret=cosh_c(VEC2(x)); RET_C;}
G2_Q_Q(cosh){IDX; float4 ret=cosh_q(VEC4(x)); RET_Q;}

G2_C_C(acosh){IDX; float2 ret=acosh_c(VEC2(x)); RET_C;}
G2_Q_Q(acosh){IDX; float4 ret=acosh_q(VEC4(x)); RET_Q;}

G2_R_R(cosc){IDX; float a=xr[idx]; ASSIGN_R(cos(a)/a);}
G2_C_C(cosc){IDX; float2 a=VEC2(x), ret=div_cc(cos_c(a), a); RET_C;}
G2_Q_Q(cosc){IDX; float4 a=VEC4(x), ret=div_qq(cos_q(a), a); RET_Q;}
DISC_R_I(cosc){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}
DISC_C_I(cosc){IDX; disc[idx]=false;}//TODO
DISC_Q_I(cosc){IDX; disc[idx]=false;}//

G2_R_R(sec){IDX; ASSIGN_R(1/cos(xr[idx]));}
G2_C_C(sec){IDX; float2 ret=inv_c(cos_c(VEC2(x))); RET_C;}
G2_Q_Q(sec){IDX; float4 ret=inv_q(cos_q(VEC4(x))); RET_Q;}
DISC_R_I(sec)
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset];
	disc[idx]=fabs(x1-x0)>3.2f||_1d_int_in_range(x0/_pi-0.5f, x1/_pi-0.5f);
}
DISC_C_I(sec)
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
DISC_Q_I(sec){IDX; disc[idx]=false;}//TODO

G2_C_C(asec){IDX; float2 ret=acos_c(inv_c(VEC2(x))); RET_C;}
G2_Q_Q(asec){IDX; float4 ret=acos_q(inv_q(VEC4(x))); RET_Q;}
DISC_C_I(asec){IDX; disc[idx]=false;}//TODO	disc c divise i
DISC_Q_I(asec){IDX; disc[idx]=false;}//		disc q divide i

G2_R_R(sech){IDX; ASSIGN_R(1/cosh(xr[idx]));}
G2_C_C(sech){IDX; float2 ret=inv_c(cosh_c(VEC2(x))); RET_C;}
G2_Q_Q(sech){IDX; float4 ret=inv_q(cosh_q(VEC4(x))); RET_Q;}
DISC_C_I(sech)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.y==x1.y)
	{
		if(x0.x==x1.x)
			disc[idx]=false;
		else
		{
			float i=x0.y/_pi+0.5f;
			disc[idx]=i==floor(i)&&(x0.x<0?x1.x>=0:x0.x>0?x1.x<=0:x1.x!=0);
		}
	}
	else if(x0.x==x1.x)
		disc[idx]=x0.x==0&&_1d_int_in_range(x0.y/_pi-0.5f, x1.y/_pi-0.5f);
	else
		disc[idx]=false;
}
DISC_Q_I(sech){IDX; disc[idx]=false;}//TODO

G2_C_C(asech){IDX; float2 ret=acosh_c(inv_c(VEC2(x))); RET_C;}
G2_Q_Q(asech){IDX; float4 ret=acosh_q(inv_q(VEC4(x))); RET_Q;}
DISC_C_I(asech)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.y==x0.y)
		disc[idx]=x0.y==0&&x0.x!=x0.x&&signbit(x0.x)!=signbit(x1.x);
	else if(x0.x==x0.x)
		disc[idx]=(x0.y<=0?x1.y>0:x1.y<=0)&&(x0.x<=0||x0.x>1);
	else if((x0.y<=0&&x1.y>0)||(x0.y<=0&&x1.y>0))
	{
		float t=_1d_zero_crossing(x0.x, x0.y, x1.x, x1.y);
		disc[idx]=t<=0||t>1;
	}
	else
		disc[idx]=false;
}
DISC_Q_I(asech){IDX; disc[idx]=false;}