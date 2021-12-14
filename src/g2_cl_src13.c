float2	tgamma_c(float2 a)//http://en.wikipedia.org/wiki/Lanczos_approximation
{
	const float g=7, p[]={0.99999999999980993f, 676.5203681218851f, -1259.1392167224028f, 771.32342877765313f, -176.61502916214059f, 12.507343278686905f, -0.13857109526572012f, 9.9843695780195716e-6f, 1.5056327351493116e-7f};
	if(a.x<0.5f)
	{
		float2 t1=(float2)(p[0], 0);
		for(int k=1;k<g+2;++k)
			t1+=div_rc(p[k], (float2)(k-a.x, -a.y));
		float2 t2=(float2)(g+0.5f-a.x, -a.y);

		float2 spa_s2p=mul_cr(sin_c(mul_rc(_pi, a)), _sqrt_2pi);
		float2 temp2=pow_cc(t2, (float2)(0.5f-a.x, -a.y));
		float2 temp3=mul_cc(exp_c((float2)(-t2.x, -t2.y)), t1);
		return div_rc(_pi, mul_cc(mul_cc(spa_s2p, temp2), temp3));
	//	return _pi/(sin(_pi*a)*_sqrt_2pi*pow(t2, 0.5f-a)*exp(-t2)*t1);//C++
	}
	else
	{
		float2 t1=(float2)(p[0], 0);
		for(int k=1;k<g+2;++k)
			t1+=div_rc(p[k], (float2)(k-1+a.x, a.y));
		float2 t2=(float2)(g+0.5f-1+a.x, a.y);

		float2 temp=mul_rc(_sqrt_2pi, pow_cc(t2, (float2)(0.5f-1+a.x, a.y)));
		float2 temp2=mul_cc(exp_c((float2)(-t2.x, -t2.y)), t1);
		return mul_cc(temp, temp2);
	//	return _sqrt_2pi*pow(t2, 0.5f-1.+x)*exp(-t2)*t1;//C++
	}
}
float4	tgamma_q(float4 a)
{
	const float g=7, p[]={0.99999999999980993f, 676.5203681218851f, -1259.1392167224028f, 771.32342877765313f, -176.61502916214059f, 12.507343278686905f, -0.13857109526572012f, 9.9843695780195716e-6f, 1.5056327351493116e-7f};
	if(a.x<0.5f)
	{
		float4 t1=(float4)(p[0], 0, 0, 0);
		for(int k=1;k<g+2;++k)
			t1+=div_rq(p[k], (float4)(k-a.x, -a.y, -a.z, -a.w));
		float4 t2=(float4)(g+0.5f-a.x, -a.y, -a.z, -a.w);

		float4 spa_s2p=mul_qr(sin_q(mul_rq(_pi, a)), _sqrt_2pi);
		float4 temp2=pow_qq(t2, (float4)(0.5f-a.x, -a.y, -a.z, -a.w));
		float4 temp3=mul_qq(exp_q((float4)(-t2.x, -t2.y, -t2.z, -t2.w)), t1);
		return div_rq(_pi, mul_qq(mul_qq(spa_s2p, temp2), temp3));
	//	return _pi/(sin(_pi*a)*_sqrt_2pi*pow(t2, 0.5f-a)*exp(-t2)*t1);//C++
	}
	else
	{
		float4 t1=(float4)(p[0], 0, 0, 0);
		for(int k=1;k<g+2;++k)
			t1+=div_rq(p[k], (float4)(k-1+a.x, a.y, a.z, a.w));
		float4 t2=(float4)(g+0.5f-1+a.x, a.y, a.z, a.w);

		float4 temp=mul_rq(_sqrt_2pi, pow_qq(t2, (float4)(0.5f-1+a.x, a.y, a.z, a.w)));
		float4 temp2=mul_qq(exp_q((float4)(-t2.x, -t2.y, -t2.z, -t2.w)), t1);
		return mul_qq(temp, temp2);
	//	return _sqrt_2pi*pow(t2, 0.5f-1.f+x)*exp(-t2)*t1;//C++
	}
}
G2_R_R(tgamma){IDX; ASSIGN_R(tgamma(xr[idx]));}
G2_C_C(tgamma)
{
	IDX;
	float2 a=VEC2(x);
	float2 ret=tgamma_c(a);
	RET_C;
}
G2_Q_Q(tgamma)
{
	IDX;
	float4 a=VEC4(x);
	float4 ret=tgamma_q(a);
	RET_Q;
}
//r_rr_tgamma		<-
DISC_R_I(tgamma){IDX; float x0=xr[idx], x1=xr[idx+offset]; disc[idx]=x0<=0&&x1<=0&&_1d_int_in_range(x0, x1);}
DISC_C_I(tgamma)
{
	IDX;
	float2 a=VEC2(x), b=(float2)(xr[idx+offset], xi[idx+offset]);
	if(a.x==b.x)
		disc[idx]=false;
	else if(a.y==b.y)
		disc[idx]=a.y==0&&_1d_int_in_range(a.x, b.x);
	else if(signbit(a.y)!=signbit(b.y))
	{
		float t=_1d_zero_crossing(a.x, a.y, b.x, b.y);
		disc[idx]=t<=0&&t==floor(t);
	}
	else
		disc[idx]=false;
}
DISC_Q_I(tgamma){IDX; disc[idx]=false;}//TODO

G2_R_R(loggamma){IDX; ASSIGN_R(lgamma(xr[idx]));}
DISC_R_I(loggamma){disc_r_tgamma_i(size, offset, disc, xr);}

G2_R_R(factorial){IDX; ASSIGN_R(tgamma(xr[idx]+1));}
G2_C_C(factorial)
{
	IDX;
	float2 a=(float2)(xr[idx]+1, xi[idx]);
	float2 ret=tgamma_c(a);
	RET_C;
}
G2_Q_Q(factorial)
{
	IDX;
	float4 a=(float4)(xr[idx]+1, xi[idx], xj[idx], xk[idx]);
	float4 ret=tgamma_q(a);
	RET_Q;
}
DISC_R_I(factorial)
{
	IDX;
	float x0=xr[idx]+1, x1=xr[idx+offset]+1;
	disc[idx]=x0<=0&&x1<=0&&_1d_int_in_range(x0, x1);
}
DISC_C_I(factorial)
{
	IDX;
	float2 a=VEC2(x), b=(float2)(xr[idx+offset], xi[idx+offset]);
	a.x+=1, b.x+=1;
	if(a.x==b.x)
		disc[idx]=false;
	else if(a.y==b.y)
		disc[idx]=a.y==0&&_1d_int_in_range(a.x, b.x);
	else if(signbit(a.y)!=signbit(b.y))
	{
		float t=_1d_zero_crossing(a.x, a.y, b.x, b.y);
		disc[idx]=t<=0&&t==floor(t);
	}
	else
		disc[idx]=false;
}
DISC_Q_I(factorial){IDX; disc[idx]=false;}//TODO

G2_R_R(permutation){IDX; ASSIGN_R(1);}
G2_C_C(permutation){IDX; ASSIGN_C(1, 0);}
G2_Q_Q(permutation){IDX; ASSIGN_Q(1, 0, 0, 0);}
G2_R_RR(permutation)
{
	IDX;
	float a=xr[idx]+1, b=yr[idx];
	ASSIGN_R(tgamma(a)/tgamma(a-b));
}
//G2_C_RC(permutation)
//{
//	IDX;
//	float a=xr[idx]+1;
//	float2 b=VEC2(y);
//	float2 ret=div_cc(tgamma_c(a), tgamma_c(a-b));
//	RET_C;
//}
//G2_Q_RQ(permutation)
//{
//	IDX;
//	float a=xr[idx]+1;
//	float4 b=VEC4(y);
//	float4 ret=div_qq(tgamma_q(a), tgamma_q(a-b));
//	RET_Q;
//}
G2_C_CR(permutation)
{
	IDX;
	float2 a=(float2)(xr[idx]+1, xi[idx]), b=(float2)(yr[idx], 0);
	float2 ret=div_cc(tgamma_c(a), tgamma_c(a-b));
	RET_C;
}
G2_C_CC(permutation)
{
	IDX;
	float2 a=(float2)(xr[idx]+1, xi[idx]), b=VEC2(y);
	float2 ret=div_cc(tgamma_c(a), tgamma_c(a-b));
	RET_C;
}
G2_Q_QQ(permutation)
{
	IDX;
	float4 a=(float4)(xr[idx]+1, xi[idx], xj[idx], xk[idx]), b=VEC4(y);
	float4 ret=div_qq(tgamma_q(a), tgamma_q(a-b));
	RET_Q;
}
DISC_RR_I(permutation){IDX; disc[idx]=false;}//TODO
DISC_CR_I(permutation){IDX; disc[idx]=false;}//
DISC_CC_I(permutation){IDX; disc[idx]=false;}//
DISC_QQ_I(permutation){IDX; disc[idx]=false;}//

G2_R_R(combination){IDX; ASSIGN_R(1);}
G2_C_C(combination){IDX; ASSIGN_C(1, 0);}
G2_Q_Q(combination){IDX; ASSIGN_Q(1, 0, 0, 0);}
G2_R_RR(combination)
{
	IDX;
	float a=xr[idx]+1, b=yr[idx];
	ASSIGN_R(tgamma(a)/(tgamma(a-b)*tgamma(b+1)));
}
G2_C_CR(combination)
{
	IDX;
	float2 a=VEC2(x);
	float2 b=(float2)(yr[idx], 0);
	a.x+=1;
	float2 ret=div_cc(tgamma_c(a), mul_cc(tgamma_c(a-b), tgamma_c((float2)(b.x+1, b.y))));
	RET_C;
}
G2_C_CC(combination)
{
	IDX;
	float2 a=VEC2(x), b=VEC2(y);
	a.x+=1;
	float2 ret=div_cc(tgamma_c(a), mul_cc(tgamma_c(a-b), tgamma_c((float2)(b.x+1, b.y))));
	RET_C;
}
G2_Q_QQ(combination)
{
	IDX;
	float4 a=VEC4(x), b=VEC4(y);
	a.x+=1;
	float4 ret=div_qq(tgamma_q(a), mul_qq(tgamma_q(a-b), tgamma_q((float4)(b.x+1, b.y, b.z, b.w))));
	RET_Q;
}
DISC_RR_I(combination){IDX; disc[idx]=false;}//TODO
DISC_CR_I(combination){IDX; disc[idx]=false;}
DISC_CC_I(combination){IDX; disc[idx]=false;}
DISC_QQ_I(combination){IDX; disc[idx]=false;}