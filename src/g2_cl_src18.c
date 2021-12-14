//aux functions for last part
float2 sign_c(float2 a){float abs_a=abs_c(a); return abs_a!=0?div_cr(a, abs_a):(float2)(0, 0);}
float4 sign_q(float4 a){float abs_a=abs_q(a); return abs_a!=0?div_qr(a, abs_a):(float4)(0, 0, 0, 0);}
float step_r(float a){return 0.5f+0.5f*sign(a);}
float2 step_c(float2 a){return (float2)(0.5f, 0)+mul_rc(0.5f, sign_c(a));}
float4 step_q(float4 a){return (float4)(0.5f, 0, 0, 0)+mul_rq(0.5f, sign_q(a));}

//back to G2 kernels
G2_R_R(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<0.5f);}
G2_R_C(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<0.5f);}
G2_R_Q(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<0.5f);}
G2_R_RR(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_RC(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_RQ(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_CR(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_CC(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_CQ(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_QR(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_QC(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
G2_R_QQ(sqwv){IDX; float a=xr[idx]; ASSIGN_R(a-floor(a)<yr[idx]);}
DISC_R_O(sqwv)//for all sqwv functions
{
	IDX;
	float x0=xr[idx], x1=xr[idx+offset];
	disc[idx]=x0!=x1;
}

float clamp01(float x)
{
	float temp=x+fabs(x);//max(0, x)
	return (temp+2-fabs(temp-2))*0.25f;//min(x, 1)
}
float trwv_dc(float x, float y)
{
	float t=x-floor(x), t2=1-x;
	t2-=floor(t2);
	float dc=clamp01(y);
	float dc2=1-dc, t_d=t/dc, t2_d2=t2/dc2;
	return (t_d<1?t_d:0)+(t2_d2<1?t2_d2:0);
}
G2_R_R(trwv)
{
	IDX;
	float a=xr[idx], t=fabs(a-floor(a)-0.5f);
	ASSIGN_R(t+t);
}
G2_R_C(trwv)
{
	IDX;
	float2 a=VEC2(x);
	float t=abs_c(a-floor_c(a)-(float2)(0.5f, 0));
	ASSIGN_R(t+t);
}
G2_R_Q(trwv)
{
	IDX;
	float4 a=VEC4(x);
	float t=abs_q(a-floor_q(a)-(float4)(0.5f, 0, 0, 0));
	ASSIGN_R(t+t);
}
G2_R_RR(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_RC(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_RQ(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_CR(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_CC(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_CQ(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_QR(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_QC(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}
G2_R_QQ(trwv){IDX; ASSIGN_R(trwv_dc(xr[idx], yr[idx]));}

float sawtooth(float x)
{
	float t=x-floor(x), t2=floor(1-t);//dc=1
	return (t2+1)*(t2*0.5f+t);
}
float sawtooth_dc(float x, float y)
{
	if(!y)
		return 0;
	float t=x-floor(x), t2=floor(y-t);
	return (t2+1)*(t2*0.5f+t)/y;
}
bool sawtooth_dc_disc(float t0, float t1){return floor(t0)!=floor(t1);}
G2_R_R(saw){IDX; ASSIGN_R(sawtooth(xr[idx]));}
G2_R_C(saw){IDX; ASSIGN_R(sawtooth(xr[idx]));}
G2_R_Q(saw){IDX; ASSIGN_R(sawtooth(xr[idx]));}
G2_R_RR(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_RC(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_RQ(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_CR(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_CC(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_CQ(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_QR(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_QC(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
G2_R_QQ(saw){IDX; ASSIGN_R(sawtooth_dc(xr[idx], yr[idx]));}
DISC_R_I(saw){IDX; disc[idx]=ceil(xr[idx])!=ceil(xr[idx+offset]);}
DISC_C_I(saw){IDX; disc[idx]=ceil(xr[idx])!=ceil(xr[idx+offset]);}
DISC_Q_I(saw){IDX; disc[idx]=ceil(xr[idx])!=ceil(xr[idx+offset]);}
DISC_RR_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_RC_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_RQ_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_CR_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_CC_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_CQ_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_QR_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_QC_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}
DISC_QQ_I(saw){IDX; disc[idx]=sawtooth_dc_disc(xr[idx]-yr[idx], xr[idx+offset]-yr[idx+offset]);}