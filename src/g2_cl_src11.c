G2_R_RR(condition_zero){IDX; ASSIGN_R(xr[idx]!=0?xr[idx]:yr[idx]);}
G2_C_RC(condition_zero)
{
	IDX;
	if(xr[idx])
		ASSIGN_C(xr[idx], 0);
	else
		ASSIGN_C(yr[idx], yi[idx]);
}
G2_Q_RQ(condition_zero)
{
	IDX;
	if(xr[idx])
		ASSIGN_Q(xr[idx], 0, 0, 0);
	else
		ASSIGN_Q(yr[idx], yi[idx], yj[idx], yk[idx]);
}
G2_C_CR(condition_zero)
{
	IDX;
	if(istrue_c(VEC2(x)))
		ASSIGN_C(xr[idx], xi[idx]);
	else
		ASSIGN_C(yr[idx], 0);
}
G2_C_CC(condition_zero)
{
	IDX;
	if(istrue_c(VEC2(x)))
		ASSIGN_C(xr[idx], xi[idx]);
	else
		ASSIGN_C(yr[idx], yi[idx]);
}
G2_Q_CQ(condition_zero)
{
	IDX;
	if(istrue_c(VEC2(x)))
		ASSIGN_Q(xr[idx], xi[idx], 0, 0);
	else
		ASSIGN_Q(yr[idx], yi[idx], yj[idx], yk[idx]);
}
G2_Q_QR(condition_zero)
{
	IDX;
	if(istrue_q(VEC4(x)))
		ASSIGN_Q(xr[idx], xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(yr[idx], 0, 0, 0);
}
G2_Q_QC(condition_zero)
{
	IDX;
	if(istrue_q(VEC4(x)))
		ASSIGN_Q(xr[idx], xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(yr[idx], yi[idx], 0, 0);
}
G2_Q_QQ(condition_zero)
{
	IDX;
	if(istrue_q(VEC4(x)))
		ASSIGN_Q(xr[idx], xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(yr[idx], yi[idx], yj[idx], yk[idx]);
}
DISC_RR_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset]);}
DISC_RC_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset]);}
DISC_RQ_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset]);}
DISC_CR_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])||_1d_zero_in_range(xi[idx], xi[idx+offset]);}
DISC_CC_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])||_1d_zero_in_range(xi[idx], xi[idx+offset]);}
DISC_CQ_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])||_1d_zero_in_range(xi[idx], xi[idx+offset]);}
DISC_QR_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])||_1d_zero_in_range(xi[idx], xi[idx+offset])||_1d_zero_in_range(xj[idx], xj[idx+offset])||_1d_zero_in_range(xk[idx], xk[idx+offset]);}
DISC_QC_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])||_1d_zero_in_range(xi[idx], xi[idx+offset])||_1d_zero_in_range(xj[idx], xj[idx+offset])||_1d_zero_in_range(xk[idx], xk[idx+offset]);}
DISC_QQ_I(condition_zero){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])||_1d_zero_in_range(xi[idx], xi[idx+offset])||_1d_zero_in_range(xj[idx], xj[idx+offset])||_1d_zero_in_range(xk[idx], xk[idx+offset]);}

G2_R_R(percent){IDX; ASSIGN_R(xr[idx]*0.01f);}
G2_C_C(percent)
{
	IDX;
	float2 ret=VEC2(x)*0.01f;
	RET_C;
}
G2_Q_Q(percent)
{
	IDX;
	float4 ret=VEC4(x)*0.01f;
	RET_Q;
}

G2_R_RR(modulo)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	ASSIGN_R(a-floor(a/b)*b);
}
G2_C_RC(modulo)
{
	IDX;
	float a=xr[idx];
	float2 b=VEC2(y);
	float2 ret=mul_cc(floor_c(div_rc(a, b)), b);
	ret.x=a-ret.x;
	ret.y=-ret.y;
	RET_C;
}
G2_Q_RQ(modulo)
{
	IDX;
	float a=xr[idx];
	float4 b=VEC4(y);
	float4 ret=mul_qq(floor_q(div_rq(a, b)), b);
	ret.x=a-ret.x;
	ret.y=-ret.y;
	ret.z=-ret.z;
	ret.w=-ret.w;
	RET_Q;
}
G2_C_CR(modulo)
{
	IDX;
	float2 a=VEC2(x);
	float b=yr[idx];
	float2 ret=mul_cr(floor_c(div_cr(a, b)), b);
	ret=a-ret;
	RET_C;
}
G2_C_CC(modulo)
{
	IDX;
	float2 a=VEC2(x);
	float2 b=VEC2(y);
	float2 ret=mul_cc(floor_c(div_cc(a, b)), b);
	ret=a-ret;
	RET_C;
}
G2_Q_CQ(modulo)
{
	IDX;
	float2 a=VEC2(x);
	float4 b=VEC4(y);
	float4 ret=mul_qq(floor_q(div_cq(a, b)), b);
	ret.x=a.x-ret.x;
	ret.y=a.y-ret.y;
	ret.z=-ret.z;
	ret.w=-ret.w;
	RET_Q;
}
G2_Q_QR(modulo)
{
	IDX;
	float4 a=VEC4(x);
	float b=yr[idx];
	float4 ret=mul_qr(floor_q(div_qr(a, b)), b);
	ret=a-ret;
	RET_Q;
}
G2_Q_QC(modulo)
{
	IDX;
	float4 a=VEC4(x);
	float2 b=VEC2(y);
	float4 ret=mul_qc(floor_q(div_qc(a, b)), b);
	ret=a-ret;
	RET_Q;
}
G2_Q_QQ(modulo)
{
	IDX;
	float4 a=VEC4(x);
	float4 b=VEC4(y);
	float4 ret=mul_qq(floor_q(div_qq(a, b)), b);
	ret=a-ret;
	RET_Q;
}
DISC_RR_I(modulo){IDX; disc[idx]=floor(xr[idx]/yr[idx])!=floor(xr[idx+offset]/yr[idx+offset]);}
DISC_RC_I(modulo){IDX; disc[idx]=false;}//TODO
DISC_RQ_I(modulo){IDX; disc[idx]=false;}
DISC_CR_I(modulo){IDX; disc[idx]=false;}
DISC_CC_I(modulo){IDX; disc[idx]=false;}
DISC_CQ_I(modulo){IDX; disc[idx]=false;}
DISC_QR_I(modulo){IDX; disc[idx]=false;}
DISC_QC_I(modulo){IDX; disc[idx]=false;}
DISC_QQ_I(modulo){IDX; disc[idx]=false;}