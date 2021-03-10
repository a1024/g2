G2_R_R(setzero){IDX; ASSIGN_R(0);}
G2_C_C(setzero){IDX; ASSIGN_C(0, 0);}
G2_Q_Q(setzero){IDX; ASSIGN_Q(0, 0, 0, 0);}

G2_R_R(ceil){IDX; ASSIGN_R(ceil(xr[idx]));}
G2_C_C(ceil){IDX; ASSIGN_C(ceil(xr[idx]), ceil(xi[idx]));}
G2_Q_Q(ceil){IDX; ASSIGN_Q(ceil(xr[idx]), ceil(xi[idx]), ceil(xj[idx]), ceil(xk[idx]));}
DISC_R_O(ceil){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}
DISC_C_O(ceil){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset];}
DISC_Q_O(ceil){IDX; disc[idx]=xr[idx]!=xr[idx+offset]||xi[idx]!=xi[idx+offset]||xj[idx]!=xj[idx+offset]||xk[idx]!=xk[idx+offset];}

G2_R_R(floor){IDX; ASSIGN_R(floor(xr[idx]));}
G2_C_C(floor){IDX; ASSIGN_C(floor(xr[idx]), floor(xi[idx]));}
G2_Q_Q(floor){IDX; ASSIGN_Q(floor(xr[idx]), floor(xi[idx]), floor(xj[idx]), floor(xk[idx]));}
DISC_R_O(floor){disc_r_ceil_o(size, offset, disc, xr);}
DISC_C_O(floor){disc_c_ceil_o(size, offset, disc, xr, xi);}
DISC_Q_O(floor){disc_q_ceil_o(size, offset, disc, xr, xi, xj, xk);}

G2_R_R(round){IDX; ASSIGN_R(round(xr[idx]));}
G2_C_C(round){IDX; ASSIGN_C(round(xr[idx]), round(xi[idx]));}
G2_Q_Q(round){IDX; ASSIGN_Q(round(xr[idx]), round(xi[idx]), round(xj[idx]), round(xk[idx]));}
DISC_R_O(round){disc_r_ceil_o(size, offset, disc, xr);}
DISC_C_O(round){disc_c_ceil_o(size, offset, disc, xr, xi);}
DISC_Q_O(round){disc_q_ceil_o(size, offset, disc, xr, xi, xj, xk);}

G2_R_R(int){IDX; ASSIGN_R((int)xr[idx]);}
G2_C_C(int){IDX; ASSIGN_C((int)xr[idx], (int)xi[idx]);}
G2_Q_Q(int){IDX; ASSIGN_Q((int)xr[idx], (int)xi[idx], (int)xj[idx], (int)xk[idx]);}
DISC_R_O(int){disc_r_ceil_o(size, offset, disc, xr);}
DISC_C_O(int){disc_c_ceil_o(size, offset, disc, xr, xi);}
DISC_Q_O(int){disc_q_ceil_o(size, offset, disc, xr, xi, xj, xk);}

G2_R_R(frac){IDX; float fxr;				ASSIGN_R(fract(xr[idx], &fxr));}
G2_C_C(frac){IDX; float fxr, fxi;			ASSIGN_C(fract(xr[idx], &fxr), fract(xi[idx], &fxi));}
G2_Q_Q(frac){IDX; float fxr, fxi, fxj, fxk;	ASSIGN_Q(fract(xr[idx], &fxr), fract(xi[idx], &fxi), fract(xj[idx], &fxj), fract(xk[idx], &fxk));}
DISC_R_I(frac){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset]);}
DISC_C_I(frac){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset]);}
DISC_Q_I(frac){IDX; disc[idx]=floor(xr[idx])!=floor(xr[idx+offset])||floor(xi[idx])!=floor(xi[idx+offset])||floor(xj[idx])!=floor(xj[idx+offset])||floor(xk[idx])!=floor(xk[idx+offset]);}