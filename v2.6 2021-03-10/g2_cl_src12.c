G2_R_R(sgn){IDX; float a=xr[idx]; ASSIGN_R((a>0)-(a<0));}
G2_C_C(sgn)
{
	IDX;
	float2 a=VEC2(x);
	float absa=abs_c(a);
	if(absa)
	{
		float2 ret=div_cr(a, absa);
		RET_C;
	}
	else
		ASSIGN_C(0, 0);
}
G2_Q_Q(sgn)
{
	IDX;
	float4 a=VEC4(x);
	float absa=abs_q(a);
	if(absa)
	{
		float4 ret=div_qr(a, absa);
		RET_Q;
	}
	else
		ASSIGN_Q(0, 0, 0, 0);
}
DISC_R_I(sgn){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}
DISC_C_I(sgn){IDX; disc[idx]=false;}//TODO
DISC_Q_I(sgn){IDX; disc[idx]=false;}//

G2_R_R(sq){IDX; float a=xr[idx]; ASSIGN_R(a*a);}
G2_C_C(sq){IDX; float2 a=VEC2(x); float2 ret=sq_c(a); RET_C;}
G2_Q_Q(sq){IDX; float4 a=VEC4(x); float4 ret=sq_q(a); RET_Q;}

G2_C_C(sqrt){IDX; float2 a=VEC2(x); float2 ret=sqrt_c(a); RET_C;}
G2_Q_Q(sqrt){IDX; float4 a=VEC4(x); float4 ret=sqrt_q(a); RET_Q;}

G2_R_R(invsqrt){IDX; ASSIGN_R(rsqrt(xr[idx]));}

G2_R_R(cbrt){IDX; ASSIGN_R(cbrt(xr[idx]));}
G2_C_C(cbrt)
{
	IDX;
	float2 a=VEC2(x);
	float2 ln_a=log_c(a);
	float2 temp=mul_cr(ln_a, 1.f/3);
	float2 ret=exp_c(temp);
	RET_C;
}
G2_Q_Q(cbrt)//TODO: optimize
{
	IDX;
	float4 a=VEC4(x);
	float4 ln_a=log_q(a);
	float4 temp=mul_qr(ln_a, 1.f/3);
	float4 ret=exp_q(temp);
	RET_Q;
}

G2_R_R(gauss)
{
	IDX;
	float a=xr[idx];
	ASSIGN_R(exp(-a*a));
}
G2_C_C(gauss)
{
	IDX;
	float2 a=VEC2(x);
	a=mul_cc(a, a);
	a.x=-a.x, a.y=-a.y;
	float2 ret=exp_c(a);
	RET_C;
}
G2_Q_Q(gauss)
{
	IDX;
	float4 a=VEC4(x);
	a=mul_qq(a, a);
	a.x=-a.x, a.y=-a.y, a.z=-a.z, a.w=-a.w;
	float4 ret=exp_q(a);
	RET_Q;
}

G2_R_R(erf){IDX; ASSIGN_R(erf(xr[idx]));}

//zeta