float bitwise_nor1(float x){return iscastable2int(x)?!f2i(x):NAN;}
float bitwise_nor2(float a, float b){return iscastable2int(a)&&iscastable2int(b)?~(f2i(a)|f2i(b)):NAN;}
G2_R_R(bitwise_nor){IDX; ASSIGN_R(bitwise_nor1(xr[idx]));}
G2_C_C(bitwise_nor){IDX; ASSIGN_C(bitwise_nor1(xr[idx]), bitwise_nor1(xi[idx]));}
G2_Q_Q(bitwise_nor){IDX; ASSIGN_Q(bitwise_nor1(xr[idx]), bitwise_nor1(xi[idx]), bitwise_nor1(xj[idx]), bitwise_nor1(xk[idx]));}
G2_R_RR(bitwise_nor){IDX; ASSIGN_R(bitwise_nor2(xr[idx], yr[idx]));}
G2_C_RC(bitwise_nor){IDX; ASSIGN_C(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(0, yi[idx]));}
G2_Q_RQ(bitwise_nor){IDX; ASSIGN_Q(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(0, yi[idx]), bitwise_nor2(0, yj[idx]), bitwise_nor2(0, yk[idx]));}
G2_C_CR(bitwise_nor){IDX; ASSIGN_C(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(xi[idx], 0));}
G2_C_CC(bitwise_nor){IDX; ASSIGN_C(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(xi[idx], yi[idx]));}
G2_Q_CQ(bitwise_nor){IDX; ASSIGN_Q(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(xi[idx], yi[idx]), bitwise_nor2(0, yj[idx]), bitwise_nor2(0, yk[idx]));}
G2_Q_QR(bitwise_nor){IDX; ASSIGN_Q(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(xi[idx], 0), bitwise_nor2(xj[idx], 0), bitwise_nor2(xk[idx], 0));}
G2_Q_QC(bitwise_nor){IDX; ASSIGN_Q(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(xi[idx], yi[idx]), bitwise_nor2(xj[idx], 0), bitwise_nor2(xk[idx], 0));}
G2_Q_QQ(bitwise_nor){IDX; ASSIGN_Q(bitwise_nor2(xr[idx], yr[idx]), bitwise_nor2(xi[idx], yi[idx]), bitwise_nor2(xj[idx], yj[idx]), bitwise_nor2(xk[idx], yk[idx]));}
DISC_R_O(bitwise_nor){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//repeat r,c,q for R_RR...Q_QQ
DISC_C_O(bitwise_nor){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_nor){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}

float bitwise_xor1(float x)
{
	if(iscastable2int(x))
	{
		int a=f2i(x);
		a^=a>>16, a^=a>>8, a^=a>>4;
		a&=15;
		return (0x6996>>a)&1;
	}
	return NAN;
}
float bitwise_xor2(float a, float b){return iscastable2int(a)&&iscastable2int(b)?f2i(a)^f2i(b):NAN;}
G2_R_R(bitwise_xor){IDX; ASSIGN_R(bitwise_xor1(xr[idx]));}
G2_C_C(bitwise_xor){IDX; ASSIGN_C(bitwise_xor1(xr[idx]), bitwise_xor1(xi[idx]));}
G2_Q_Q(bitwise_xor){IDX; ASSIGN_Q(bitwise_xor1(xr[idx]), bitwise_xor1(xi[idx]), bitwise_xor1(xj[idx]), bitwise_xor1(xk[idx]));}
G2_R_RR(bitwise_xor){IDX; ASSIGN_R(bitwise_xor2(xr[idx], yr[idx]));}
G2_C_RC(bitwise_xor){IDX; ASSIGN_C(bitwise_xor2(xr[idx], yr[idx]), yi[idx]);}
G2_Q_RQ(bitwise_xor){IDX; ASSIGN_Q(bitwise_xor2(xr[idx], yr[idx]), yi[idx], yj[idx], yk[idx]);}
G2_C_CR(bitwise_xor){IDX; ASSIGN_C(bitwise_xor2(xr[idx], yr[idx]), xi[idx]);}
G2_C_CC(bitwise_xor){IDX; ASSIGN_C(bitwise_xor2(xr[idx], yr[idx]), bitwise_xor2(xi[idx], yi[idx]));}
G2_Q_CQ(bitwise_xor){IDX; ASSIGN_Q(bitwise_xor2(xr[idx], yr[idx]), bitwise_xor2(xi[idx], yi[idx]), yj[idx], yk[idx]);}
G2_Q_QR(bitwise_xor){IDX; ASSIGN_Q(bitwise_xor2(xr[idx], yr[idx]), xi[idx], xj[idx], xk[idx]);}
G2_Q_QC(bitwise_xor){IDX; ASSIGN_Q(bitwise_xor2(xr[idx], yr[idx]), bitwise_xor2(xi[idx], yi[idx]), xj[idx], xk[idx]);}
G2_Q_QQ(bitwise_xor){IDX; ASSIGN_Q(bitwise_xor2(xr[idx], yr[idx]), bitwise_xor2(xi[idx], yi[idx]), bitwise_xor2(xj[idx], yj[idx]), bitwise_xor2(xk[idx], yk[idx]));}
DISC_R_O(bitwise_xor){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//repeat r,c,q for R_RR...Q_QQ
DISC_C_O(bitwise_xor){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_xor){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}

float bitwise_xnor1(float x)
{
	if(iscastable2int(x))
	{
		int a=f2i(x);
		a^=a>>16, a^=a>>8, a^=a>>4;
		a&=15;
		return !((0x6996>>a)&1);
	}
	return NAN;
}
float bitwise_xnor2(float a, float b){return iscastable2int(a)&&iscastable2int(b)?~(f2i(a)^f2i(b)):NAN;}
G2_R_R(bitwise_xnor){IDX; ASSIGN_R(bitwise_xnor1(xr[idx]));}
G2_C_C(bitwise_xnor){IDX; ASSIGN_C(bitwise_xnor1(xr[idx]), bitwise_xnor1(xi[idx]));}
G2_Q_Q(bitwise_xnor){IDX; ASSIGN_Q(bitwise_xnor1(xr[idx]), bitwise_xnor1(xi[idx]), bitwise_xnor1(xj[idx]), bitwise_xnor1(xk[idx]));}
G2_R_RR(bitwise_xnor){IDX; ASSIGN_R(bitwise_xnor2(xr[idx], yr[idx]));}
G2_C_RC(bitwise_xnor){IDX; ASSIGN_C(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(0, yi[idx]));}
G2_Q_RQ(bitwise_xnor){IDX; ASSIGN_Q(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(0, yi[idx]), bitwise_xnor2(0, yj[idx]), bitwise_xnor2(0, yk[idx]));}
G2_C_CR(bitwise_xnor){IDX; ASSIGN_C(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(xi[idx], 0));}
G2_C_CC(bitwise_xnor){IDX; ASSIGN_C(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(xi[idx], yi[idx]));}
G2_Q_CQ(bitwise_xnor){IDX; ASSIGN_Q(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(xi[idx], yi[idx]), bitwise_xnor2(0, yj[idx]), bitwise_xnor2(0, yk[idx]));}
G2_Q_QR(bitwise_xnor){IDX; ASSIGN_Q(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(xi[idx], 0), bitwise_xnor2(xj[idx], 0), bitwise_xnor2(xk[idx], 0));}
G2_Q_QC(bitwise_xnor){IDX; ASSIGN_Q(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(xi[idx], yi[idx]), bitwise_xnor2(xj[idx], 0), bitwise_xnor2(xk[idx], 0));}
G2_Q_QQ(bitwise_xnor){IDX; ASSIGN_Q(bitwise_xnor2(xr[idx], yr[idx]), bitwise_xnor2(xi[idx], yi[idx]), bitwise_xnor2(xj[idx], yj[idx]), bitwise_xnor2(xk[idx], yk[idx]));}
DISC_R_O(bitwise_xnor){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//repeat r,c,q for R_RR...Q_QQ
DISC_C_O(bitwise_xnor){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_xnor){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}