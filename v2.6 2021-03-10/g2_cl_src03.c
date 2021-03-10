G2_R_RR(plus){IDX; ASSIGN_R(xr[idx]+yr[idx]);}
G2_C_RC(plus){IDX; ASSIGN_C(xr[idx]+yr[idx], yi[idx]);}
G2_Q_RQ(plus){IDX; ASSIGN_Q(xr[idx]+yr[idx], yi[idx], yj[idx], yk[idx]);}
G2_C_CR(plus){IDX; ASSIGN_C(xr[idx]+yr[idx], xi[idx]);}
G2_C_CC(plus){IDX; ASSIGN_C(xr[idx]+yr[idx], xi[idx]+yi[idx]);}
G2_Q_CQ(plus){IDX; ASSIGN_Q(xr[idx]+yr[idx], xi[idx]+yi[idx], yj[idx], yk[idx]);}
G2_Q_QR(plus){IDX; ASSIGN_Q(xr[idx]+yr[idx], xi[idx], xj[idx], xk[idx]);}
G2_Q_QC(plus){IDX; ASSIGN_Q(xr[idx]+yr[idx], xi[idx]+yi[idx], xj[idx], xk[idx]);}
G2_Q_QQ(plus){IDX; ASSIGN_Q(xr[idx]+yr[idx], xi[idx]+yi[idx], xj[idx]+yj[idx], xk[idx]+yk[idx]);}

G2_R_R(minus){IDX; ASSIGN_R(-xr[idx]);}
G2_C_C(minus){IDX; ASSIGN_C(-xr[idx], -xi[idx]);}
G2_Q_Q(minus){IDX; ASSIGN_Q(-xr[idx], -xi[idx], -xj[idx], -xk[idx]);}
G2_R_RR(minus){IDX; ASSIGN_R(xr[idx]-yr[idx]);}
G2_C_RC(minus){IDX; ASSIGN_C(xr[idx]-yr[idx], -yi[idx]);}
G2_Q_RQ(minus){IDX; ASSIGN_Q(xr[idx]-yr[idx], -yi[idx], -yj[idx], -yk[idx]);}
G2_C_CR(minus){IDX; ASSIGN_C(xr[idx]-yr[idx], xi[idx]);}
G2_C_CC(minus){IDX; ASSIGN_C(xr[idx]-yr[idx], xi[idx]-yi[idx]);}
G2_Q_CQ(minus){IDX; ASSIGN_Q(xr[idx]-yr[idx], xi[idx]-yi[idx], -yj[idx], -yk[idx]);}
G2_Q_QR(minus){IDX; ASSIGN_Q(xr[idx]-yr[idx], xi[idx], xj[idx], xk[idx]);}
G2_Q_QC(minus){IDX; ASSIGN_Q(xr[idx]-yr[idx], xi[idx]-yi[idx], xj[idx], xk[idx]);}
G2_Q_QQ(minus){IDX; ASSIGN_Q(xr[idx]-yr[idx], xi[idx]-yi[idx], xj[idx]-yj[idx], xk[idx]-yk[idx]);}

G2_R_RR(multiply){IDX; ASSIGN_R(xr[idx]*yr[idx]);}
G2_C_RC(multiply){IDX; float2 ret=mul_rc(xr[idx], VEC2(y)); RET_C;}
G2_Q_RQ(multiply){IDX; float4 ret=mul_rq(xr[idx], VEC4(y)); RET_Q;}
G2_C_CR(multiply){IDX; float2 ret=mul_cr(VEC2(x), yr[idx]); RET_C;}
G2_C_CC(multiply){IDX; float2 ret=mul_cc(VEC2(x), VEC2(y)); RET_C;}
G2_Q_CQ(multiply){IDX; float4 ret=mul_cq(VEC2(x), VEC4(y)); RET_Q;}
G2_Q_QR(multiply){IDX; float4 ret=mul_qr(VEC4(x), yr[idx]); RET_Q;}
G2_Q_QC(multiply){IDX; float4 ret=mul_qc(VEC4(x), VEC2(y)); RET_Q;}
G2_Q_QQ(multiply){IDX; float4 ret=mul_qq(VEC4(x), VEC4(y)); RET_Q;}

G2_R_R(divide){IDX; rr[idx]=1/xr[idx];}
G2_C_C(divide){IDX; float2 ret=inv_c(VEC2(x)); RET_C;}
G2_Q_Q(divide){IDX; float4 ret=inv_q(VEC4(x)); RET_Q;}
G2_R_RR(divide){IDX; ASSIGN_R(xr[idx]/yr[idx]);}
G2_C_RC(divide){IDX; float2 ret=div_rc(xr[idx], VEC2(y)); RET_C;}
G2_Q_RQ(divide){IDX; float4 ret=div_rq(xr[idx], VEC4(y)); RET_Q;}
G2_C_CR(divide){IDX; float2 ret=div_cr(VEC2(x), yr[idx]); RET_C;}
G2_C_CC(divide){IDX; float2 ret=div_cc(VEC2(x), VEC2(y)); RET_C;}
G2_Q_CQ(divide){IDX; float4 ret=div_cq(VEC2(x), VEC4(y)); RET_Q;}
G2_Q_QR(divide){IDX; float4 ret=div_qr(VEC4(x), yr[idx]); RET_Q;}
G2_Q_QC(divide){IDX; float4 ret=div_qc(VEC4(x), VEC2(y)); RET_Q;}
G2_Q_QQ(divide){IDX; float4 ret=div_qq(VEC4(x), VEC4(y)); RET_Q;}
DISC_R_I(divide){IDX; float x0r=xr[idx], x1r=xr[idx+offset]; disc[idx]=x0r<0?x1r>=0:x0r>0?x1r<=0:x1r!=0;}
DISC_C_I(divide){IDX; disc[idx]=false;}//TODO
DISC_Q_I(divide){IDX; disc[idx]=false;}//
DISC_RR_I(divide){IDX; float y0r=yr[idx], y1r=yr[idx+offset]; disc[idx]=y0r<0?y1r>=0:y0r>0?y1r<=0:y1r!=0;}
DISC_RC_I(divide){IDX; disc[idx]=false;}//
DISC_RQ_I(divide){IDX; disc[idx]=false;}//
DISC_CR_I(divide){IDX; float y0r=yr[idx], y1r=yr[idx+offset]; disc[idx]=y0r<0?y1r>=0:y0r>0?y1r<=0:y1r!=0;}
DISC_CC_I(divide){IDX; disc[idx]=false;}//
DISC_CQ_I(divide){IDX; disc[idx]=false;}//
DISC_QR_I(divide){IDX; float y0r=yr[idx], y1r=yr[idx+offset]; disc[idx]=y0r<0?y1r>=0:y0r>0?y1r<=0:y1r!=0;}
DISC_QC_I(divide){IDX; disc[idx]=false;}//
DISC_QQ_I(divide){IDX; disc[idx]=false;}//