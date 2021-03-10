G2_R_RR(logic_divides){IDX; float q=xr[idx]/yr[idx]; rr[idx]=q==floor(q);}
G2_R_RC(logic_divides){IDX; float2 q=div_rc(xr[idx], VEC2(y)); rr[idx]=equal_cc(q, floor_c(q));}
G2_R_RQ(logic_divides){IDX; float4 q=div_rq(xr[idx], VEC4(y)); rr[idx]=equal_qq(q, floor_q(q));}
G2_R_CR(logic_divides){IDX; float2 q=div_cr(VEC2(x), yr[idx]); rr[idx]=equal_cc(q, floor_c(q));}
G2_R_CC(logic_divides){IDX; float2 q=div_cc(VEC2(x), VEC2(y)); rr[idx]=equal_cc(q, floor_c(q));}
G2_R_CQ(logic_divides){IDX; float4 q=div_cq(VEC2(x), VEC4(y)); rr[idx]=equal_qq(q, floor_q(q));}
G2_R_QR(logic_divides){IDX; float4 q=div_qr(VEC4(x), yr[idx]); rr[idx]=equal_qq(q, floor_q(q));}
G2_R_QC(logic_divides){IDX; float4 q=div_qc(VEC4(x), VEC2(y)); rr[idx]=equal_qq(q, floor_q(q));}
G2_R_QQ(logic_divides){IDX; float4 q=div_qq(VEC4(x), VEC4(y)); rr[idx]=equal_qq(q, floor_q(q));}
DISC_R_O(logic_divides){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//for all logic_divides functions

G2_R_RR(power_real)
{
	IDX;
	float a=xr[idx];
	int b=(int)yr[idx];
	float m[]={1, 0}, r=1;//may not compile
	if(b<0)
		m[1]=1/a, b=-b;
	else
		m[1]=a;
	for(;;)
	{
		r*=m[b&1], b>>=1;
		if(!b)
			break;
		m[1]*=m[1];
	}
	rr[idx]=r;
}
G2_C_CR(power_real)
{
	IDX;
	float2 a=VEC2(x), ret=(float2)(1, 0);
	int b=(int)yr[idx], p=abs(b);
	for(;;)
	{
		if(p&1)
			ret=mul_cc(ret, a);
		p>>=1;
		if(!p)
			break;
		a=mul_cc(a, a);
	}
	if(b<0)
		ret=inv_c(ret);
	RET_C;
}
G2_Q_QR(power_real)
{
	IDX;
	float4 a=VEC4(x), ret=(float4)(1, 0, 0, 0);
	int b=(int)yr[idx], p=abs(b);
	for(;;)
	{
		if(p&1)
			ret=mul_qq(ret, a);
		p>>=1;
		if(!p)
			break;
		a=mul_qq(a, a);
	}
	if(b<0)
		ret=inv_q(ret);
	RET_Q;
}
DISC_RR_I(power_real){IDX; disc[idx]=floor(yr[idx])!=floor(yr[idx+offset]);}
DISC_CR_I(power_real){IDX; disc[idx]=floor(yr[idx])!=floor(yr[idx+offset]);}
DISC_QR_I(power_real){IDX; disc[idx]=floor(yr[idx])!=floor(yr[idx+offset]);}

G2_C_CR(pow)
{
	IDX;
	float2 a=VEC2(x);
	float b=yr[idx];
	if(!a.x&&!a.y && !b)
		ASSIGN_C(1, 0);
	else
	{
		float2 ret=pow_cr(a, b);
		RET_C;
	}
}
G2_C_CC(pow)
{
	IDX;
	float2 a=VEC2(x);
	float2 b=VEC2(y);
	if(!a.x&&!a.y && !b.x&&!b.y)
		ASSIGN_C(1, 0);
	else
	{
		float2 ret=pow_cc(a, b);
		RET_C;
	}
}
G2_Q_CQ(pow)
{
	IDX;
	float2 a=VEC2(x);
	float4 b=VEC4(y);
	if(!a.x&&!a.y && !b.x&&!b.y&&!b.z&&!b.w)
		ASSIGN_Q(1, 0, 0, 0);
	else
	{
		float4 ret=pow_cq(a, b);
		RET_Q;
	}
}
G2_Q_QR(pow)
{
	IDX;
	float4 a=VEC4(x);
	float b=yr[idx];
	if(!a.x&&!a.y&&!a.z&&!a.w && !b)
		ASSIGN_Q(1, 0, 0, 0);
	else
	{
		float4 ret=pow_qr(a, b);
		RET_Q;
	}
}
G2_Q_QC(pow)
{
	IDX;
	float4 a=VEC4(x);
	float2 b=VEC2(y);
	if(!a.x&&!a.y&&!a.z&&!a.w && !b.x&&!b.y)
		ASSIGN_Q(1, 0, 0, 0);
	else
	{
		float4 ret=pow_qc(a, b);
		RET_Q;
	}
}
G2_Q_QQ(pow)
{
	IDX;
	float4 a=VEC4(x);
	float4 b=VEC4(y);
	if(!a.x&&!a.y&&!a.z&&!a.w && !b.x&&!b.y&&!b.z&&!b.w)
		ASSIGN_Q(1, 0, 0, 0);
	else
	{
		float4 ret=pow_qq(a, b);
		RET_Q;
	}
}
DISC_CR_I(pow){IDX; disc[idx]=false;}//TODO
DISC_CC_I(pow){IDX; disc[idx]=false;}//TODO
DISC_CQ_I(pow){IDX; disc[idx]=false;}//TODO
DISC_QR_I(pow){IDX; disc[idx]=false;}//TODO
DISC_QC_I(pow){IDX; disc[idx]=false;}//TODO
DISC_QQ_I(pow){IDX; disc[idx]=false;}//TODO

G2_C_C(ln){IDX; float2 ret=log_c(VEC2(x)); RET_C;}
G2_Q_Q(ln){IDX; float4 ret=log_q(VEC4(x)); RET_Q;}
DISC_C_I(ln){IDX; disc[idx]=false;}//TODO
DISC_Q_I(ln){IDX; disc[idx]=false;}

G2_C_C(log)
{
	IDX;
	float2 log_a=log_c(VEC2(x));
	float2 ret=mul_cr(log_a, _inv_ln10);
	RET_C;
}
G2_Q_Q(log)
{
	IDX;
	float4 log_a=log_q(VEC4(x));
	float4 ret=mul_qr(log_a, _inv_ln10);
	RET_Q;
}
G2_C_CR(log)
{
	IDX;
	float2 log_a=log_c(VEC2(x));
	float2 log_b=log_c((float2)(yr[idx], 0));
	float2 ret=div_cc(log_a, log_b);
	RET_C;
}
G2_C_CC(log)
{
	IDX;
	float2 log_a=log_c(VEC2(x));
	float2 log_b=log_c(VEC2(y));
	float2 ret=div_cc(log_a, log_b);
	RET_C;
}
G2_Q_CQ(log)
{
	IDX;
	float2 log_a=log_c(VEC2(x));
	float4 log_b=log_q(VEC4(y));
	float4 ret=div_cq(log_a, log_b);
	RET_Q;
}
G2_Q_QC(log)
{
	IDX;
	float4 log_a=log_q(VEC4(x));
	float2 log_b=log_c(VEC2(y));
	float4 ret=div_qc(log_a, log_b);
	RET_Q;
}
G2_Q_QQ(log)
{
	IDX;
	float4 log_a=log_q(VEC4(x));
	float4 log_b=log_q(VEC4(y));
	float4 ret=div_qq(log_a, log_b);
	RET_Q;
}
DISC_C_I(log){IDX; disc[idx]=false;}//TODO disc c arg
DISC_Q_I(log){IDX; disc[idx]=false;}//TODO
DISC_CR_I(log){IDX; disc[idx]=false;}//TODO
DISC_CC_I(log){IDX; disc[idx]=false;}//TODO
DISC_CQ_I(log){IDX; disc[idx]=false;}//TODO
DISC_QC_I(log){IDX; disc[idx]=false;}//TODO
DISC_QQ_I(log){IDX; disc[idx]=false;}//TODO

//tetrate

//pentate