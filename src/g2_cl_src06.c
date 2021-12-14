float bitwise_not(float x){return iscastable2int(x)?~f2i(x):NAN;}
G2_R_R(bitwise_not){IDX; float a=xr[idx]; ASSIGN_R(bitwise_not(a));}
G2_C_C(bitwise_not){IDX; float2 a=VEC2(x); ASSIGN_C(bitwise_not(a.x), bitwise_not(a.y));}
G2_Q_Q(bitwise_not){IDX; float4 a=VEC4(x); ASSIGN_Q(bitwise_not(a.x), bitwise_not(a.y), bitwise_not(a.z), bitwise_not(a.w));}
DISC_R_I(bitwise_not){IDX; disc[idx]=sign(xr[idx])!=sign(xr[idx+offset]);}
DISC_C_I(bitwise_not){IDX; disc[idx]=sign(xr[idx])!=sign(xr[idx+offset])||sign(xi[idx])!=sign(xi[idx+offset]);}
DISC_Q_I(bitwise_not){IDX; disc[idx]=sign(xr[idx])!=sign(xr[idx+offset])||sign(xi[idx])!=sign(xi[idx+offset])||sign(xj[idx])!=sign(xj[idx+offset])||sign(xk[idx])!=sign(xk[idx+offset]);}

float bitwise_and1(float x){return iscastable2int(x)?f2i(x)==-1:NAN;}
float bitwise_and2(float a, float b){return iscastable2int(a)&&iscastable2int(b)?f2i(a)&f2i(b):NAN;}
G2_R_R(bitwise_and){IDX; ASSIGN_R(bitwise_and1(xr[idx]));}
G2_C_C(bitwise_and){IDX; ASSIGN_C(bitwise_and1(xr[idx]), bitwise_and1(xi[idx]));}
G2_Q_Q(bitwise_and){IDX; ASSIGN_Q(bitwise_and1(xr[idx]), bitwise_and1(xi[idx]), bitwise_and1(xj[idx]), bitwise_and1(xk[idx]));}
G2_R_RR(bitwise_and){IDX; ASSIGN_R(bitwise_and2(xr[idx], yr[idx]));}
G2_C_RC(bitwise_and){IDX; ASSIGN_C(bitwise_and2(xr[idx], yr[idx]), 0);}
G2_Q_RQ(bitwise_and){IDX; ASSIGN_Q(bitwise_and2(xr[idx], yr[idx]), 0, 0, 0);}
G2_C_CR(bitwise_and){IDX; ASSIGN_C(bitwise_and2(xr[idx], yr[idx]), 0);}
G2_C_CC(bitwise_and){IDX; ASSIGN_C(bitwise_and2(xr[idx], yr[idx]), bitwise_and2(xi[idx], yi[idx]));}
G2_Q_CQ(bitwise_and){IDX; ASSIGN_Q(bitwise_and2(xr[idx], yr[idx]), bitwise_and2(xi[idx], yi[idx]), 0, 0);}
G2_Q_QR(bitwise_and){IDX; ASSIGN_Q(bitwise_and2(xr[idx], yr[idx]), 0, 0, 0);}
G2_Q_QC(bitwise_and){IDX; ASSIGN_Q(bitwise_and2(xr[idx], yr[idx]), bitwise_and2(xi[idx], yi[idx]), 0, 0);}
G2_Q_QQ(bitwise_and){IDX; ASSIGN_Q(bitwise_and2(xr[idx], yr[idx]), bitwise_and2(xi[idx], yi[idx]), bitwise_and2(xj[idx], yj[idx]), bitwise_and2(xk[idx], yk[idx]));}
DISC_R_O(bitwise_and){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//repeat r,c,q for R_RR...Q_QQ
DISC_C_O(bitwise_and){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_and){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}

float bitwise_nand1(float x){return iscastable2int(x)?f2i(x)!=-1:NAN;}
float bitwise_nand2(float a, float b){return iscastable2int(a)&&iscastable2int(b)?~(f2i(a)&f2i(b)):NAN;}
G2_R_R(bitwise_nand){IDX; ASSIGN_R(bitwise_nand1(xr[idx]));}
G2_C_C(bitwise_nand){IDX; ASSIGN_C(bitwise_nand1(xr[idx]), bitwise_nand1(xi[idx]));}
G2_Q_Q(bitwise_nand){IDX; ASSIGN_Q(bitwise_nand1(xr[idx]), bitwise_nand1(xi[idx]), bitwise_nand1(xj[idx]), bitwise_nand1(xk[idx]));}
G2_R_RR(bitwise_nand){IDX; ASSIGN_R(bitwise_nand2(xr[idx], yr[idx]));}
G2_C_RC(bitwise_nand){IDX; ASSIGN_C(bitwise_nand2(xr[idx], yr[idx]), 0);}
G2_Q_RQ(bitwise_nand){IDX; ASSIGN_Q(bitwise_nand2(xr[idx], yr[idx]), 0, 0, 0);}
G2_C_CR(bitwise_nand){IDX; ASSIGN_C(bitwise_nand2(xr[idx], yr[idx]), 0);}
G2_C_CC(bitwise_nand){IDX; ASSIGN_C(bitwise_nand2(xr[idx], yr[idx]), bitwise_nand2(xi[idx], yi[idx]));}
G2_Q_CQ(bitwise_nand){IDX; ASSIGN_Q(bitwise_nand2(xr[idx], yr[idx]), bitwise_nand2(xi[idx], yi[idx]), 0, 0);}
G2_Q_QR(bitwise_nand){IDX; ASSIGN_Q(bitwise_nand2(xr[idx], yr[idx]), 0, 0, 0);}
G2_Q_QC(bitwise_nand){IDX; ASSIGN_Q(bitwise_nand2(xr[idx], yr[idx]), bitwise_nand2(xi[idx], yi[idx]), 0, 0);}
G2_Q_QQ(bitwise_nand){IDX; ASSIGN_Q(bitwise_nand2(xr[idx], yr[idx]), bitwise_nand2(xi[idx], yi[idx]), bitwise_nand2(xj[idx], yj[idx]), bitwise_nand2(xk[idx], yk[idx]));}
DISC_R_O(bitwise_nand){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//repeat r,c,q for R_RR...Q_QQ
DISC_C_O(bitwise_nand){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_nand){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}

float bitwise_or1(float x){return iscastable2int(x)?f2i(x)!=0:NAN;}
float bitwise_or2(float a, float b){return iscastable2int(a)&&iscastable2int(b)?f2i(a)|f2i(b):NAN;}
G2_R_R(bitwise_or){IDX; ASSIGN_R(bitwise_or1(xr[idx]));}
G2_C_C(bitwise_or){IDX; ASSIGN_C(bitwise_or1(xr[idx]), bitwise_or1(xi[idx]));}
G2_Q_Q(bitwise_or){IDX; ASSIGN_Q(bitwise_or1(xr[idx]), bitwise_or1(xi[idx]), bitwise_or1(xj[idx]), bitwise_or1(xk[idx]));}
G2_R_RR(bitwise_or){IDX; ASSIGN_R(bitwise_or2(xr[idx], yr[idx]));}
G2_C_RC(bitwise_or){IDX; ASSIGN_C(bitwise_or2(xr[idx], yr[idx]), yi[idx]);}
G2_Q_RQ(bitwise_or){IDX; ASSIGN_Q(bitwise_or2(xr[idx], yr[idx]), yi[idx], yj[idx], yk[idx]);}
G2_C_CR(bitwise_or){IDX; ASSIGN_C(bitwise_or2(xr[idx], yr[idx]), xi[idx]);}
G2_C_CC(bitwise_or){IDX; ASSIGN_C(bitwise_or2(xr[idx], yr[idx]), bitwise_or2(xi[idx], yi[idx]));}
G2_Q_CQ(bitwise_or){IDX; ASSIGN_Q(bitwise_or2(xr[idx], yr[idx]), bitwise_or2(xi[idx], yi[idx]), yj[idx], yk[idx]);}
G2_Q_QR(bitwise_or){IDX; ASSIGN_Q(bitwise_or2(xr[idx], yr[idx]), xi[idx], xj[idx], xk[idx]);}
G2_Q_QC(bitwise_or){IDX; ASSIGN_Q(bitwise_or2(xr[idx], yr[idx]), bitwise_or2(xi[idx], yi[idx]), xj[idx], xk[idx]);}
G2_Q_QQ(bitwise_or){IDX; ASSIGN_Q(bitwise_or2(xr[idx], yr[idx]), bitwise_or2(xi[idx], yi[idx]), bitwise_or2(xj[idx], yj[idx]), bitwise_or2(xk[idx], yk[idx]));}
DISC_R_O(bitwise_or){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}//repeat r,c,q for R_RR...Q_QQ
DISC_C_O(bitwise_or){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(bitwise_or){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}