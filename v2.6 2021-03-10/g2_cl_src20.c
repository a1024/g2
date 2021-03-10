G2_R_RR(conditional_110){IDX; ASSIGN_R(xr[idx]!=0?yr[idx]:0);}
G2_C_RC(conditional_110){IDX; float2 ret=xr[idx]!=0?VEC2(y):(float2)(0, 0); RET_C;}
G2_Q_RQ(conditional_110){IDX; float4 ret=xr[idx]!=0?VEC4(y):(float4)(0, 0, 0, 0); RET_Q;}
G2_R_CR(conditional_110){IDX; ASSIGN_R(istrue_c(VEC2(x))?yr[idx]:0);}
G2_C_CC(conditional_110){IDX; float2 ret=istrue_c(VEC2(x))?VEC2(y):(float2)(0, 0); RET_C;}
G2_Q_CQ(conditional_110){IDX; float4 ret=istrue_c(VEC2(x))?VEC4(y):(float4)(0, 0, 0, 0); RET_Q;}
G2_R_QR(conditional_110){IDX; ASSIGN_R(istrue_q(VEC4(x))?yr[idx]:0);}
G2_C_QC(conditional_110){IDX; float2 ret=istrue_q(VEC4(x))?VEC2(y):(float2)(0, 0); RET_C;}
G2_Q_QQ(conditional_110){IDX; float4 ret=istrue_q(VEC4(x))?VEC4(y):(float4)(0, 0, 0, 0); RET_Q;}
DISC_RR_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset]);}
DISC_RC_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset]);}
DISC_RQ_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset]);}
DISC_CR_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])&&_1d_zero_in_range(xi[idx], xi[idx+offset]);}
DISC_CC_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])&&_1d_zero_in_range(xi[idx], xi[idx+offset]);}
DISC_CQ_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])&&_1d_zero_in_range(xi[idx], xi[idx+offset]);}
DISC_QR_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])&&_1d_zero_in_range(xi[idx], xi[idx+offset])&&_1d_zero_in_range(xj[idx], xj[idx+offset])&&_1d_zero_in_range(xk[idx], xk[idx+offset]);}
DISC_QC_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])&&_1d_zero_in_range(xi[idx], xi[idx+offset])&&_1d_zero_in_range(xj[idx], xj[idx+offset])&&_1d_zero_in_range(xk[idx], xk[idx+offset]);}
DISC_QQ_I(conditional_110){IDX; disc[idx]=_1d_zero_in_range(xr[idx], xr[idx+offset])&&_1d_zero_in_range(xi[idx], xi[idx+offset])&&_1d_zero_in_range(xj[idx], xj[idx+offset])&&_1d_zero_in_range(xk[idx], xk[idx+offset]);}

G2_R_RR(conditional_101){IDX; ASSIGN_R(xr[idx]!=0?yr[idx]:0);}
G2_C_RC(conditional_101){IDX; float2 ret=xr[idx]!=0?VEC2(y):(float2)(0, 0); RET_C;}
G2_Q_RQ(conditional_101){IDX; float4 ret=xr[idx]!=0?VEC4(y):(float4)(0, 0, 0, 0); RET_Q;}
G2_R_CR(conditional_101){IDX; ASSIGN_R(istrue_c(VEC2(x))?yr[idx]:0);}
G2_C_CC(conditional_101){IDX; float2 ret=istrue_c(VEC2(x))?VEC2(y):(float2)(0, 0); RET_C;}
G2_Q_CQ(conditional_101){IDX; float4 ret=istrue_c(VEC2(x))?VEC4(y):(float4)(0, 0, 0, 0); RET_Q;}
G2_R_QR(conditional_101){IDX; ASSIGN_R(istrue_q(VEC4(x))?yr[idx]:0);}
G2_C_QC(conditional_101){IDX; float2 ret=istrue_q(VEC4(x))?VEC2(y):(float2)(0, 0); RET_C;}
G2_Q_QQ(conditional_101){IDX; float4 ret=istrue_q(VEC4(x))?VEC4(y):(float4)(0, 0, 0, 0); RET_Q;}
DISC_RR_I(conditional_101){disc_rr_conditional_110_i(size, offset, disc, xr, yr);}
DISC_RC_I(conditional_101){disc_rc_conditional_110_i(size, offset, disc, xr, yr, yi);}
DISC_RQ_I(conditional_101){disc_rq_conditional_110_i(size, offset, disc, xr, yr, yi, yj, yk);}
DISC_CR_I(conditional_101){disc_cr_conditional_110_i(size, offset, disc, xr, xi, yr);}
DISC_CC_I(conditional_101){disc_cc_conditional_110_i(size, offset, disc, xr, xi, yr, yi);}
DISC_CQ_I(conditional_101){disc_cq_conditional_110_i(size, offset, disc, xr, xi, yr, yi, yj, yk);}
DISC_QR_I(conditional_101){disc_qr_conditional_110_i(size, offset, disc, xr, xi, xj, xk, yr);}
DISC_QC_I(conditional_101){disc_qc_conditional_110_i(size, offset, disc, xr, xi, xj, xk, yr, yi);}
DISC_QQ_I(conditional_101){disc_qq_conditional_110_i(size, offset, disc, xr, xi, xj, xk, yr, yi, yj, yk);}

//conditional_111: pass nullptr for real & complex			TODO: disc for conditional_111
#define		NPA(pointer)				(pointer?pointer[idx]:0)
#define		ASSIGN_NP(dest, source)		if(dest)dest[idx]=NPA(source)
__kernel void conditional_111(__global const int *size,
	__global float *rr, __global float *ri, __global float *rj, __global float *rk,
	__global const float *xr, __global const float *xi, __global const float *xj, __global const float *xk,
	__global const float *yr, __global const float *yi, __global const float *yj, __global const float *yk,
	__global const float *zr, __global const float *zi, __global const float *zj, __global const float *zk)
{
	IDX;
	float4 a=(float4)(NPA(xr), NPA(xi), NPA(xj), NPA(xk));
	if(istrue_q(a))
	{
		ASSIGN_NP(rr, yr);
		ASSIGN_NP(ri, yi);
		ASSIGN_NP(rj, yj);
		ASSIGN_NP(rk, yk);
	}
	else
	{
		ASSIGN_NP(rr, zr);
		ASSIGN_NP(ri, zi);
		ASSIGN_NP(rj, zj);
		ASSIGN_NP(rk, zk);
	}
}

G2_R_R(increment){IDX; ASSIGN_R(xr[idx]+1);}
G2_C_C(increment){IDX; float2 ret=VEC2(x)+(float2)(1, 0); RET_C;}
G2_Q_Q(increment){IDX; float4 ret=VEC4(x)+(float4)(1, 0, 0, 0); RET_Q;}

G2_R_R(decrement){IDX; ASSIGN_R(xr[idx]-1);}
G2_C_C(decrement){IDX; float2 ret=VEC2(x)-(float2)(1, 0); RET_C;}
G2_Q_Q(decrement){IDX; float4 ret=VEC4(x)-(float4)(1, 0, 0, 0); RET_Q;}

G2_R_R(assign){IDX; ASSIGN_R(xr[idx]);}
G2_C_C(assign){IDX; float2 ret=VEC2(x); RET_C;}
G2_Q_Q(assign){IDX; float4 ret=VEC4(x); RET_Q;}

__kernel void initialize_parameter(__global const int *size, __global float *buffer, __global const float *args)
{
	IDX;
//	uint idx=size[0]*(size[1]*get_global_id(2)+get_global_id(1))+get_global_id(0);
	//args[0]: start
	//args[1]: ratio
	//args[2]: dimension 0:x, 1:y, 2:z
	buffer[idx]=args[0]+args[1]*get_global_id((int)args[2]);
//	buffer[idx]=args[0]+args[1]*get_global_id(((int*)args)[2]);
}

float alpha_from_line(float x1, float y1, float x2, float y2)
{
	float a=x2-x1, b=y1-y2;
	return max(0.f, 1.f-fabs(a*y1+b*x1)*rsqrt(a*a+b*b));
}
float do_quadrant(float m, float R, float U, float UR)
{
	bool down=(m>0)!=(R>0), right=(R>0)!=(UR>0), up=(UR>0)!=(U>0), left=(U>0)!=(m>0);
//	return (float)(down||right||up||left);//
	float yL=m/(m-U), xD=m/(m-R), yR=R/(R-UR), xU=U/(U-UR);
	if(left)
	{
		if(down)//case 4 & 16
			return alpha_from_line(0, yL, xD, 0);
		if(right)//case 7
			return alpha_from_line(0, yL, 1, yR);
		return alpha_from_line(0, yL, xU, 1);//up	case 11
	}
	if(down)//cases 6 & 10
	{
		if(right)//	case 6
			return alpha_from_line(xD, 0, 1, yR);
		return alpha_from_line(xD, 0, xU, 1);//up	case 10
	}
	if(up)//&&right	case 13
		return alpha_from_line(xU, 1, 1, yR);
	return 0;//case 1
}
__kernel void ti2d_rgb(__global const int *size, __global const float *xr, __global const float *curvecolor,
#ifdef G2_OCL_IMAGES
	__write_only image2d_t rgb
#else
	__global int *rgb
#endif
	)
{//size{Xplaces, Yplaces}
	const uint kx=get_global_id(0), ky=get_global_id(1);
	const int w=size[0], h=size[1], idx=size[0]*ky+kx;
	const int
		xsub=1&-(kx>0), xadd=1&-(kx<w-1),//
		ysub=w&-(ky>0), yadd=w&-(ky<h-1);
	float
		Vnw=xr[idx-ysub-xsub], Vnm=xr[idx-ysub], Vne=xr[idx-ysub+xadd],
		Vmw=xr[idx     -xsub], Vmm=xr[idx     ], Vme=xr[idx     +xadd],
		Vsw=xr[idx+yadd-xsub], Vsm=xr[idx+yadd], Vse=xr[idx+yadd+xadd];
	float alpha=Vmm==0;
	if(alpha!=1)
	{
		alpha=do_quadrant(Vmm, Vme, Vnm, Vne);
		alpha=max(alpha, do_quadrant(Vmm, Vnm, Vmw, Vnw));
		alpha=max(alpha, do_quadrant(Vmm, Vmw, Vsm, Vsw));
		alpha=max(alpha, do_quadrant(Vmm, Vsm, Vme, Vse));
	}
#ifdef G2_OCL_IMAGES
	write_imagef(rgb, (int2)(kx, ky), (float4)(curvecolor[0], curvecolor[1], curvecolor[2], alpha));
//	write_imagef(rgb, (int2)(kx, ky), (float4)((1-alpha)*curvecolor[0], (1-alpha)*curvecolor[1], (1-alpha)*curvecolor[2], 1));
#else
	rgb[w*ky+kx]=(uchar)(255*alpha)<<24|(uchar)(255*curvecolor[2])<<16|(uchar)(255*curvecolor[1])<<8|(uchar)(255*curvecolor[0]);
#endif
}
#define		COS_PI_6		0.866025403784439f
#define		SIN_PI_6		0.5f
#define		THRESHOLD		10
#define		INV_THRESHOLD	0.1f
#define		COMP_MUL		0.00392156862745098f
__kernel void c2d_rgb(__global const int *size, __global const float *xr, __global const float *xi,
#ifdef G2_OCL_IMAGES
	__write_only image2d_t rgb
#else
	__global int *rgb
#endif
	)
{//size{Xplaces, Yplaces}
	const int2 coords=(int2)(get_global_id(0), get_global_id(1));
	const uint idx=size[0]*coords.y+coords.x;
	float r=xr[idx], i=xi[idx];
	float4 color;
	if(r!=r||i!=i)
		color=(float4)(1, 0.5f, 0.5f, 0.5f);
	else if(fabs(r)==INFINITY||fabs(i)==INFINITY)
		color=(float4)(1, 1, 1, 1);
	else
	{
		float hyp=sqrt(r*r+i*i), cosx=r/hyp, sinx=i/hyp,
			mag=255*exp(-hyp*_ln2*INV_THRESHOLD);
		float red=1+cosx*COS_PI_6-sinx*SIN_PI_6, green=1+sinx, blue=1+cosx*-COS_PI_6-sinx*SIN_PI_6;
		if(hyp<THRESHOLD)
			mag=255-mag, red*=mag, green*=mag, blue*=mag;
		else
			red=255-mag*(2-red), green=255-mag*(2-green), blue=255-mag*(2-blue);
		color=(float4)(red, green, blue, 255)*COMP_MUL;
	//	color.x=red*COMP_MUL, color.y=green*COMP_MUL, color.z=blue*COMP_MUL, color.w=1;
	//	color=(float4)(1, blue*COMP_MUL, green*COMP_MUL, red*COMP_MUL);
	//	rgb[idx]=0xFF<<24|(uchar)blue<<16|(uchar)green<<8|(uchar)red;
	}
//	color=(float4)(1, 1, coords.x/size[0], coords.y/size[1]);//
#ifdef G2_OCL_IMAGES
	write_imagef(rgb, coords, color);
#else
	rgb[idx]=0xFF000000|(uchar)(255*color.z)<<16|(uchar)(255*color.y)<<8|(uchar)(255*color.x);
#endif
}

__kernel void c2d_rgb2(__global const int *size, __global const float *xr, __global const float *xi, __global int *rgb)
{//size{Xplaces, Yplaces, 1, appwidth, appheight}
	const int2 coords=(int2)(get_global_id(0), get_global_id(1));
	const uint idx=size[0]*coords.y+coords.x, idx2=size[3]*coords.y+coords.x;
	float r=xr[idx], i=xi[idx];
	if(r!=r||i!=i)
		rgb[idx2]=0xFF7F7F7F;
	else if(fabs(r)==INFINITY||fabs(i)==INFINITY)
		rgb[idx2]=0xFFFFFFFF;
	else
	{
		float hyp=sqrt(r*r+i*i), cosx=r/hyp, sinx=i/hyp,
			mag=255*exp(-hyp*_ln2*INV_THRESHOLD);
		float red=1+cosx*COS_PI_6-sinx*SIN_PI_6, green=1+sinx, blue=1+cosx*-COS_PI_6-sinx*SIN_PI_6;
		if(hyp<THRESHOLD)
			mag=255-mag, red*=mag, green*=mag, blue*=mag;
		else
			red=255-mag*(2-red), green=255-mag*(2-green), blue=255-mag*(2-blue);
		rgb[idx2]=0xFF000000|(uchar)blue<<16|(uchar)green<<8|(uchar)red;
	}
}