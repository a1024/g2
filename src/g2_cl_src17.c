G2_R_R(exp){IDX; ASSIGN_R(exp(xr[idx]));}
G2_C_C(exp){IDX; float2 ret=exp_c(VEC2(x)); RET_C;}
G2_Q_Q(exp){IDX; float4 ret=exp_q(VEC4(x)); RET_Q;}

G2_R_R(fib){IDX; float a=xr[idx], phi_p_x=exp(a*_ln_phi); ASSIGN_R((phi_p_x-cos(_pi*a)/phi_p_x)*_inv_sqrt5);}
G2_C_C(fib)
{
	IDX;
	float2 a=VEC2(x);
	float2 phi_p_x=exp_c(mul_cr(a, _ln_phi));
	float2 ret=mul_cr(phi_p_x-div_cc(cos_c(mul_rc(_pi, a)), phi_p_x), _inv_sqrt5);
	RET_C;
}
G2_Q_Q(fib)
{
	IDX;
	float4 a=VEC4(x);
	float4 phi_p_x=exp_q(mul_qr(a, _ln_phi));
	float4 ret=mul_qr(phi_p_x-div_qq(cos_q(mul_rq(_pi, a)), phi_p_x), _inv_sqrt5);
	RET_Q;
}

//random		//TODO

//beta

//bessel_j (x86: cyl_bessel_j)

//bessel_y (x86: cyl_neumann)

//hankel1

//aux functions for last part
float2 sign_c(float2 a){float abs_a=abs_c(a); return abs_a!=0?div_cr(a, abs_a):(float2)(0, 0);}
float4 sign_q(float4 a){float abs_a=abs_q(a); return abs_a!=0?div_qr(a, abs_a):(float4)(0, 0, 0, 0);}
float step_r(float a){return 0.5f+0.5f*sign(a);}
float2 step_c(float2 a){return (float2)(0.5f, 0)+mul_rc(0.5f, sign_c(a));}
float4 step_q(float4 a){return (float4)(0.5f, 0, 0, 0)+mul_rq(0.5f, sign_q(a));}

//back to kernels
G2_R_R(step){IDX; ASSIGN_R(step_r(xr[idx]));}
G2_C_C(step){IDX; float2 ret=step_c(VEC2(x)); RET_C;}
G2_Q_Q(step){IDX; float4 ret=step_q(VEC4(x)); RET_Q;}
DISC_R_I(step){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}
DISC_C_I(step){IDX; disc[idx]=false;}//TODO
DISC_Q_I(step){IDX; disc[idx]=false;}//

G2_R_R(rect){IDX; float a=xr[idx]; ASSIGN_R(step_r(a+0.5f)-step_r(a-0.5f));}
G2_C_C(rect){IDX; const float2 temp=(float2)(0.5f, 0); float2 a=VEC2(x), ret=step_c(a+temp)-step_c(a-temp); RET_C;}//'half' means half precision
G2_Q_Q(rect){IDX; const float4 temp=(float4)(0.5f, 0, 0, 0); float4 a=VEC4(x), ret=step_q(a+temp)-step_q(a-temp); RET_Q;}
bool disc_rect(float x0, float x1)
{
	const float2 d=(float2)(-0.5f, 0.5f);
		 if(x0<d.x)		return x1>=d.x;
	else if(x0==d.x)	return x1<d.x||x1>d.x;
	else if(x0<d.y)		return x1>=d.y;
	else if(x0==d.y)	return x1<d.y||x1>d.y;
	else				return x1<=d.y;
}
DISC_R_I(rect){IDX; disc[idx]=disc_rect(xr[idx], xr[idx+offset]);}
DISC_C_I(rect)
{
	IDX;
	float2 x0=VEC2(x), x1=(float2)(xr[idx+offset], xi[idx+offset]);
	if(x0.y==x1.y)
		disc[idx]=x0.y==0&&disc_rect(x0.x, x1.x);
	else if(x0.x==x1.x)
		disc[idx]=(x0.x==-0.5f||x0.x==0.5f)&&(x0.y<0?x1.y>=0:x0.y==0?x1.y<0||x1.y>0:x1.y<=0);
	else
	{
		float x=_1d_zero_crossing(x0.x, x0.y, x1.x, x1.y);
		disc[idx]=x==-0.5f||x==0.5f;
	}
}
DISC_Q_I(rect){IDX; disc[idx]=false;}//TODO

G2_R_R(trgl){IDX; float abs_a=fabs(xr[idx]); ASSIGN_R((abs_a<1)*(1-abs_a));}
G2_R_C(trgl){IDX; float abs_a=abs_c(VEC2(x)); ASSIGN_R((abs_a<1)*(1-abs_a));}
G2_R_Q(trgl){IDX; float abs_a=abs_q(VEC4(x)); ASSIGN_R((abs_a<1)*(1-abs_a));}