G2_R_R(bitwise_shift_left_l){IDX; ASSIGN_R(exp(floor(xr[idx])*_ln2));}
G2_C_C(bitwise_shift_left_l)
{
	IDX;
	float2 ret=exp_c(mul_cr(floor_c(VEC2(x)), _ln2));
	RET_C;
}
G2_Q_Q(bitwise_shift_left_l)
{
	IDX;
	float4 ret=exp_q(mul_qr(floor_q(VEC4(x)), _ln2));
	RET_Q;
}
G2_R_R(bitwise_shift_left_r){IDX; ASSIGN_R(xr[idx]+xr[idx]);}
G2_C_C(bitwise_shift_left_r)
{
	IDX;
	float2 ret=VEC2(x);
	ret+=ret;
	RET_C;
}
G2_Q_Q(bitwise_shift_left_r)
{
	IDX;
	float4 ret=VEC4(x);
	ret+=ret;
	RET_Q;
}
G2_R_RR(bitwise_shift_left){IDX; ASSIGN_R(xr[idx]*exp(floor(yr[idx])*_ln2));}
G2_C_RC(bitwise_shift_left)
{
	IDX;
	float a=xr[idx];
	float2 b=VEC2(y);
	float2 ret=mul_rc(a, exp_c(mul_cr(floor_c(b), _ln2)));
	RET_C;
}
G2_Q_RQ(bitwise_shift_left)
{
	IDX;
	float a=xr[idx];
	float4 b=VEC4(y);
	float4 ret=mul_rq(a, exp_q(mul_qr(floor_q(b), _ln2)));
	RET_Q;
}
G2_C_CR(bitwise_shift_left)
{
	IDX;
	float2 a=VEC2(x);
	float b=yr[idx];
	float2 ret=mul_cr(a, exp(floor(b)*_ln2));
	RET_C;
}
G2_C_CC(bitwise_shift_left)
{
	IDX;
	float2 a=VEC2(x);
	float2 b=VEC2(y);
	float2 ret=mul_cc(a, exp_c(mul_cr(floor_c(b), _ln2)));
	RET_C;
}
G2_Q_CQ(bitwise_shift_left)
{
	IDX;
	float2 a=VEC2(x);
	float4 b=VEC4(y);
	float4 ret=mul_cq(a, exp_q(mul_qr(floor_q(b), _ln2)));
	RET_Q;
}
G2_Q_QR(bitwise_shift_left)
{
	IDX;
	float4 a=VEC4(x);
	float b=yr[idx];
	float4 ret=mul_qr(a, exp(floor(b)*_ln2));
	RET_Q;
}
G2_Q_QC(bitwise_shift_left)
{
	IDX;
	float4 a=VEC4(x);
	float2 b=VEC2(y);
	float4 ret=mul_qc(a, exp_c(mul_cr(floor_c(b), _ln2)));
	RET_Q;
}
G2_Q_QQ(bitwise_shift_left)
{
	IDX;
	float4 a=VEC4(x);
	float4 b=VEC4(y);
	float4 ret=mul_qq(a, exp_q(mul_qr(floor_q(b), _ln2)));
	RET_Q;
}
DISC_R_O(bitwise_shift_left_l){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}
DISC_C_O(bitwise_shift_left_l){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_shift_left_l){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}
DISC_RR_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset]);}
DISC_RC_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset]);}
DISC_RQ_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset]);}
DISC_CR_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset]);}
DISC_CC_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset]);}
DISC_CQ_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset]);}
DISC_QR_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset])||floor(xj[idx])!=floor(xj[idx+offset])||floor(xk[idx])!=floor(xk[idx+offset]);}
DISC_QC_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset])||floor(xj[idx])!=floor(xj[idx+offset])||floor(xk[idx])!=floor(xk[idx+offset]);}
DISC_QQ_I(bitwise_shift_left){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset])||floor(xj[idx])!=floor(xj[idx+offset])||floor(xk[idx])!=floor(xk[idx+offset]);}

G2_R_R(bitwise_shift_right_l){IDX; ASSIGN_R(exp(-floor(xr[idx])*_ln2));}
G2_C_C(bitwise_shift_right_l)
{
	IDX;
	float2 ret=exp_c(mul_cr(floor_c(VEC2(x)), -_ln2));
	RET_C;
}
G2_Q_Q(bitwise_shift_right_l)
{
	IDX;
	float4 ret=exp_q(mul_qr(floor_q(VEC4(x)), -_ln2));
	RET_Q;
}
G2_R_R(bitwise_shift_right_r){IDX; ASSIGN_R(xr[idx]*0.5f);}
G2_C_C(bitwise_shift_right_r)
{
	IDX;
	float2 ret=VEC2(x)*0.5f;
	RET_C;
}
G2_Q_Q(bitwise_shift_right_r)
{
	IDX;
	float4 ret=VEC4(x)*0.5f;
	RET_Q;
}
G2_R_RR(bitwise_shift_right){IDX; ASSIGN_R(xr[idx]*exp(-floor(yr[idx])*_ln2));}
G2_C_RC(bitwise_shift_right)
{
	IDX;
	float a=xr[idx];
	float2 b=VEC2(y);
	float2 ret=mul_rc(a, exp_c(mul_cr(floor_c(b), -_ln2)));
	RET_C;
}
G2_Q_RQ(bitwise_shift_right)
{
	IDX;
	float a=xr[idx];
	float4 b=VEC4(y);
	float4 ret=mul_rq(a, exp_q(mul_qr(floor_q(b), -_ln2)));
	RET_Q;
}
G2_C_CR(bitwise_shift_right)
{
	IDX;
	float2 a=VEC2(x);
	float b=yr[idx];
	float2 ret=mul_cr(a, exp(-floor(b)*_ln2));
	RET_C;
}
G2_C_CC(bitwise_shift_right)
{
	IDX;
	float2 a=VEC2(x);
	float2 b=VEC2(y);
	float2 ret=mul_cc(a, exp_c(mul_cr(floor_c(b), -_ln2)));
	RET_C;
}
G2_Q_CQ(bitwise_shift_right)
{
	IDX;
	float2 a=VEC2(x);
	float4 b=VEC4(y);
	float4 ret=mul_cq(a, exp_q(mul_qr(floor_q(b), -_ln2)));
	RET_Q;
}
G2_Q_QR(bitwise_shift_right)
{
	IDX;
	float4 a=VEC4(x);
	float b=yr[idx];
	float4 ret=mul_qr(a, exp(-floor(b)*_ln2));
	RET_Q;
}
G2_Q_QC(bitwise_shift_right)
{
	IDX;
	float4 a=VEC4(x);
	float2 b=VEC2(y);
	float4 ret=mul_qc(a, exp_c(mul_cr(floor_c(b), -_ln2)));
	RET_Q;
}
G2_Q_QQ(bitwise_shift_right)
{
	IDX;
	float4 a=VEC4(x);
	float4 b=VEC4(y);
	float4 ret=mul_qq(a, exp_q(mul_qr(floor_q(b), -_ln2)));
	RET_Q;
}
DISC_R_O(bitwise_shift_right_l){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}
DISC_C_O(bitwise_shift_right_l){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_shift_right_l){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}
DISC_RR_I(bitwise_shift_right){disc_rr_bitwise_shift_left_i(size, offset, disc, xr, yr);}
DISC_RC_I(bitwise_shift_right){disc_rc_bitwise_shift_left_i(size, offset, disc, xr, yr, yi);}
DISC_RQ_I(bitwise_shift_right){disc_rq_bitwise_shift_left_i(size, offset, disc, xr, yr, yi, yj, yk);}
DISC_CR_I(bitwise_shift_right){disc_cr_bitwise_shift_left_i(size, offset, disc, xr, xi, yr);}
DISC_CC_I(bitwise_shift_right){disc_cc_bitwise_shift_left_i(size, offset, disc, xr, xi, yr, yi);}
DISC_CQ_I(bitwise_shift_right){disc_cq_bitwise_shift_left_i(size, offset, disc, xr, xi, yr, yi, yj, yk);}
DISC_QR_I(bitwise_shift_right){disc_qr_bitwise_shift_left_i(size, offset, disc, xr, xi, xj, xk, yr);}
DISC_QC_I(bitwise_shift_right){disc_qc_bitwise_shift_left_i(size, offset, disc, xr, xi, xj, xk, yr, yi);}
DISC_QQ_I(bitwise_shift_right){disc_qq_bitwise_shift_left_i(size, offset, disc, xr, xi, xj, xk, yr, yi, yj, yk);}