//math constants
#define	_pi			3.14159265358979f
#define	_2pi		6.28318530717959f
#define	_sqrt_2pi	2.506628274631f
#define	_sqrt2		1.4142135623731f
#define	_ln2		0.693147180559945f
#define	_inv_ln10	0.434294481903252f
#define	_ln_phi		0.481211825059603f
#define	_inv_sqrt5	0.447213595499958f

//macros for G2 functions
#define		ARG_CI(arg)		__global const int *arg
#define		ARG_F(arg)		__global float *arg
#define		ARG_CF(arg)		__global const float *arg
#define		ARG_I(arg)		__global int *arg
//#define	ARG_C(arg)		__global char *arg//needs cl_khr_byte_addressable_store
#define		G2_R_R(func)	__kernel void  r_r_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr))
#define		G2_C_C(func)	__kernel void  c_c_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(xi))
#define		G2_Q_Q(func)	__kernel void  q_q_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_F(rj), ARG_F(rk), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk))
#define		G2_R_RR(func)	__kernel void r_rr_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(yr))
#define		G2_C_RC(func)	__kernel void c_rc_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(yr), ARG_CF(yi))
#define		G2_Q_RQ(func)	__kernel void q_rq_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_F(rj), ARG_F(rk), ARG_CF(xr), ARG_CF(yr), ARG_CF(yi), ARG_CF(yj), ARG_CF(yk))
#define		G2_C_CR(func)	__kernel void c_cr_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(xi), ARG_CF(yr))
#define		G2_C_CC(func)	__kernel void c_cc_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(xi), ARG_CF(yr), ARG_CF(yi))
#define		G2_Q_CQ(func)	__kernel void q_cq_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_F(rj), ARG_F(rk), ARG_CF(xr), ARG_CF(xi), ARG_CF(yr), ARG_CF(yi), ARG_CF(yj), ARG_CF(yk))
#define		G2_Q_QR(func)	__kernel void q_qr_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_F(rj), ARG_F(rk), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr))
#define		G2_Q_QC(func)	__kernel void q_qc_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_F(rj), ARG_F(rk), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr), ARG_CF(yi))
#define		G2_Q_QQ(func)	__kernel void q_qq_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_F(rj), ARG_F(rk), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr), ARG_CF(yi), ARG_CF(yj), ARG_CF(yk))
#define		G2_C_R(func)	__kernel void  c_r_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr))
#define		G2_C_Q(func)	__kernel void  c_q_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk))
#define		G2_R_C(func)	__kernel void  r_c_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi))
#define		G2_R_Q(func)	__kernel void  r_q_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk))
#define		G2_C_RR(func)	__kernel void c_rr_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(yr))
#define		G2_R_RC(func)	__kernel void r_rc_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(yr), ARG_CF(yi))
#define		G2_R_RQ(func)	__kernel void r_rq_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(yr), ARG_CF(yi), ARG_CF(yj), ARG_CF(yk))
#define		G2_R_CR(func)	__kernel void r_cr_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(yr))
#define		G2_R_CC(func)	__kernel void r_cc_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(yr), ARG_CF(yi))
#define		G2_R_CQ(func)	__kernel void r_cq_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(yr), ARG_CF(yi), ARG_CF(yj), ARG_CF(yk))
#define		G2_R_QR(func)	__kernel void r_qr_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr))
#define		G2_R_QC(func)	__kernel void r_qc_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr), ARG_CF(yi))
#define		G2_R_QQ(func)	__kernel void r_qq_##func(ARG_CI(size), ARG_F(rr), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr), ARG_CF(yi), ARG_CF(yj), ARG_CF(yk))
#define		G2_C_QC(func)	__kernel void c_qc_##func(ARG_CI(size), ARG_F(rr), ARG_F(ri), ARG_CF(xr), ARG_CF(xi), ARG_CF(xj), ARG_CF(xk), ARG_CF(yr), ARG_CF(yi))

#define		DISC_R_O(func)	__kernel void disc_r_##func##_o(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr))
#define		DISC_C_O(func)	__kernel void disc_c_##func##_o(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi))
#define		DISC_Q_O(func)	__kernel void disc_q_##func##_o(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(xj), ARG_F(xk))
#define		DISC_R_I(func)	__kernel void disc_r_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr))
#define		DISC_C_I(func)	__kernel void disc_c_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi))
#define		DISC_Q_I(func)	__kernel void disc_q_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(xj), ARG_F(xk))
#define		DISC_RR_I(func)	__kernel void disc_rr_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(yr))
#define		DISC_RC_I(func)	__kernel void disc_rc_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(yr), ARG_F(yi))
#define		DISC_RQ_I(func)	__kernel void disc_rq_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(yr), ARG_F(yi), ARG_F(yj), ARG_F(yk))
#define		DISC_CR_I(func)	__kernel void disc_cr_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(yr))
#define		DISC_CC_I(func)	__kernel void disc_cc_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(yr), ARG_F(yi))
#define		DISC_CQ_I(func)	__kernel void disc_cq_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(yr), ARG_F(yi), ARG_F(yj), ARG_F(yk))
#define		DISC_QR_I(func)	__kernel void disc_qr_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(xj), ARG_F(xk), ARG_F(yr))
#define		DISC_QC_I(func)	__kernel void disc_qc_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(xj), ARG_F(xk), ARG_F(yr), ARG_F(yi))
#define		DISC_QQ_I(func)	__kernel void disc_qq_##func##_i(ARG_CI(size), const int offset, ARG_I(disc), ARG_F(xr), ARG_F(xi), ARG_F(xj), ARG_F(xk), ARG_F(yr), ARG_F(yi), ARG_F(yj), ARG_F(yk))

#define		IDX						const unsigned idx=get_idx(size)
#define		ASSIGN_R(r)				rr[idx]=r
#define		ASSIGN_C(r, i)			rr[idx]=r, ri[idx]=i
#define		ASSIGN_Q(r, i, j, k)	rr[idx]=r, ri[idx]=i, rj[idx]=j, rk[idx]=k
#define		RET_C					ASSIGN_C(ret.x, ret.y)
#define		RET_Q					ASSIGN_Q(ret.x, ret.y, ret.z, ret.w)
#define		VEC2(name)				(float2)(name##r[idx], name##i[idx])
#define		VEC4(name)				(float4)(name##r[idx], name##i[idx], name##j[idx], name##k[idx])

//auxiliary functions
int		get_idx(__global const int *size){return size[0]*(size[1]*get_global_id(2)+get_global_id(1))+get_global_id(0);}
int		iscastable2int(float x){return x>=-2147483647.f&&x<=2147483647.f;}
int		f2i(float x){return (int)floor(x);}
bool	_1d_int_in_range(float a, float b){return floor(a)!=floor(b)||ceil(a)!=ceil(b);}
bool	_1d_zero_in_range(float a, float b){return a<0?b>=0:a==0?b!=0:b<0;}
float	_1d_zero_crossing(float x0, float y0, float x1, float y1){return x0+(0-y0)*(x1-x0)/(y1-y0);}
bool	istrue_c(float2 a){return a.x||a.y;}
bool	istrue_q(float4 a){return a.x||a.y||a.z||a.w;}
float	abs_c(float2 a){return sqrt(a.x*a.x+a.y*a.y);}
float	abs_q(float4 a){return sqrt(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);}
bool	equal_rc(float a, float2 b){return a==b.x&&0==b.y;}
bool	equal_rq(float a, float4 b){return a==b.x&&0==b.y&&0==b.z&&0==b.w;}
bool	equal_cr(float2 a, float b){return a.x==b&&a.y==0;}
bool	equal_cc(float2 a, float2 b){return a.x==b.x&&a.y==b.y;}
bool	equal_cq(float2 a, float4 b){return a.x==b.x&&a.y==b.y&&0==b.z&&0==b.w;}
bool	equal_qr(float4 a, float b){return a.x==b&&a.y==0&&a.z==0&&a.w==0;}
bool	equal_qc(float4 a, float2 b){return a.x==b.x&&a.y==b.y&&a.z==0&&a.w==0;}
bool	equal_qq(float4 a, float4 b){return a.x==b.x&&a.y==b.y&&a.z==b.z&&a.w==b.w;}
float2	floor_c(float2 a){return (float2)(floor(a.x), floor(a.y));}
float4	floor_q(float4 a){return (float4)(floor(a.x), floor(a.y), floor(a.z), floor(a.w));}
float2	sq_c(float2 a){float ri=a.x*a.y; return (float2)(a.x*a.x-a.y*a.y, ri+ri);}
float4	sq_q(float4 a){float _2r=a.x+a.x; return (float4)(a.x*a.x-a.y*a.y-a.z*a.z-a.w*a.w, a.y*_2r, a.z*_2r, a.w*_2r);}
float2	sqrt_c(float2 a)
{
	float s=abs_c(a)+a.x;
	if(s)//entire complex plane except origin & -ve real axis
	{
		s=sqrt(s+s);
		return (float2)(s*0.5f, a.y/s);
	}
	return (float2)(0, sqrt(-a.x));
}
float4	sqrt_q(float4 a)
{
	float s=abs_q(a)+a.x;
	if(s)
	{
		s=sqrt(s+s);
		float inv_s=1/s;
		return (float4)(s*0.5f, a.y*inv_s, a.z*inv_s, a.w*inv_s);
	}
	return (float4)(0, sqrt(-a.x), 0, 0);
}
float2	mul_rc(float a, float2 b){return (float2)(a*b.x, a*b.y);}
float4	mul_rq(float a, float4 b){return (float4)(a*b.x, a*b.y, a*b.z, a*b.w);}
float2	mul_cr(float2 a, float b){return (float2)(a.x*b, a.y*b);}
float2	mul_cc(float2 a, float2 b){return (float2)(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);}
float4	mul_cq(float2 a, float4 b)
{
	return (float4)
	(
		a.x*b.x-a.y*b.y,
		a.x*b.y+a.y*b.x,
		a.x*b.z-a.y*b.w,
		a.x*b.w+a.y*b.z
	);
}
float4	mul_qr(float4 a, float b){return (float4)(a.x*b, a.y*b, a.z*b, a.w*b);}
float4	mul_qc(float4 a, float2 b)
{
	return (float4)
	(
		 a.x*b.x-a.y*b.y,
		 a.x*b.y+a.y*b.x,
		 a.z*b.x+a.w*b.y,
		-a.z*b.y+a.w*b.x
	);
}
float4	mul_qq(float4 a, float4 b)
{
	return (float4)
	(
		a.x*b.x-a.y*b.y-a.z*b.z-a.w*b.w,
		a.x*b.y+a.y*b.x+a.z*b.w-a.w*b.z,
		a.x*b.z-a.y*b.w+a.z*b.x+a.w*b.y,
		a.x*b.w+a.y*b.z-a.z*b.y+a.w*b.x
	);
}
float2	inv_c(float2 a)
{
	float inv_abs2=1/(a.x*a.x+a.y*a.y);
	return (float2)(a.x*inv_abs2, -a.y*inv_abs2);
}
float4	inv_q(float4 a)
{
	float inv_abs2=1/(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);
	return (float4)(a.x*inv_abs2, -a.y*inv_abs2, -a.z*inv_abs2, -a.w*inv_abs2);
}
float2	div_rc(float a, float2 b)
{
	float a_absb2=a/(b.x*b.x+b.y*b.y);
	return (float2)(b.x*a_absb2, -b.y*a_absb2);
}
float4	div_rq(float a, float4 b)
{
	float a_absb2=a/(b.x*b.x+b.y*b.y+b.z*b.z+b.w*b.w);
	return (float4)(b.x*a_absb2, -b.y*a_absb2, -b.z*a_absb2, -b.w*a_absb2);
}
float2	div_cr(float2 a, float b)
{
	float inv_b=1/b;
	return (float2)(a.x*inv_b, a.y*inv_b);
}
float2	div_cc(float2 a, float2 b)
{
	float inv_absb2=1/(b.x*b.x+b.y*b.y);
	return (float2)
	(
		(a.x*b.x+a.y*b.y)*inv_absb2,
		(a.y*b.x-a.x*b.y)*inv_absb2
	);
}
float4	div_cq(float2 a, float4 b)
{
	float inv_absb2=1/(b.x*b.x+b.y*b.y+b.z*b.z+b.w*b.w);
	return (float4)
	(
		( a.x*b.x+a.y*b.y)*inv_absb2,
		( a.y*b.x-a.x*b.y)*inv_absb2,
		(-a.x*b.z-a.y*b.w)*inv_absb2,
		( a.y*b.z-a.x*b.w)*inv_absb2
	);
}
float4	div_qr(float4 a, float b)
{
	float inv_b=1/b;
	return (float4)
	(
		a.x*inv_b,
		a.y*inv_b,
		a.z*inv_b,
		a.w*inv_b
	);
}
float4	div_qc(float4 a, float2 b)
{
	float inv_absb2=1/(b.x*b.x+b.y*b.y);
	return (float4)
	(
		(a.x*b.x+a.y*b.y)*inv_absb2,
		(a.y*b.x-a.x*b.y)*inv_absb2,
		(a.z*b.x+a.w*b.y)*inv_absb2,
		(a.w*b.x-a.z*b.y)*inv_absb2
	);
}
float4	div_qq(float4 a, float4 b)
{
	float inv_absb2=1/(b.x*b.x+b.y*b.y+b.z*b.z+b.w*b.w);
	return (float4)
	(
		(a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w)*inv_absb2,
		(a.y*b.x-a.x*b.y-a.w*b.z+a.z*b.w)*inv_absb2,
		(a.z*b.x+a.w*b.y-a.x*b.z-a.y*b.w)*inv_absb2,
		(a.w*b.x-a.z*b.y+a.y*b.z-a.x*b.w)*inv_absb2
	);
}
float2	exp_c(float2 a)
{
	float exp_r=exp(a.x);
	return (float2)
	(
		exp_r*cos(a.y),
		exp_r*sin(a.y)
	);
}
float4	exp_q(float4 a)
{
	float exp_r=exp(a.x), abs_v=sqrt(a.y*a.y+a.z*a.z+a.w*a.w);
	float cos_v, sin_v=sincos(abs_v, &cos_v);
	float v_mul=exp_r*sin_v/abs_v;
	return (float4)
	(
		exp_r*cos_v,
		a.y*v_mul,
		a.z*v_mul,
		a.w*v_mul
	);
}
float2	log_c(float2 a)
{
	return (float2)
	(
		log(sqrt(a.x*a.x+a.y*a.y)),
		atan2(a.y, a.x)
	);
}
float4	log_q(float4 a)
{
	float absv2=a.y*a.y+a.z*a.z+a.w*a.w, abs_a=sqrt(a.x*a.x+absv2);
	float v_mul=acos(a.x/abs_a)*rsqrt(absv2);
	return (float4)
	(
		log(abs_a),
		a.y*v_mul,
		a.z*v_mul,
		a.w*v_mul
	);
}
//complex trig & hyp functions
float2 coshsinh(float x)
{
	float exp_x=exp(x);
	float2 e=(float2)(exp_x, 1/exp_x)*0.5f;
	return (float2)(e.x+e.y, e.x-e.y);
}
float2 cos_c(float2 a){return (float2)(cos(a.x)*cosh(a.y), -sin(a.x)*sinh(a.y));}
float4 cos_q(float4 a)
{
	float cos_r, sin_r=sincos(a.x, &cos_r);
	float abs_v=sqrt(a.y*a.y+a.z*a.z+a.w*a.w);
	float2 chsh_absv=coshsinh(abs_v);
	float v_mul=-sin_r*chsh_absv.y/abs_v;
	return (float4)(cos_r*chsh_absv.x, v_mul*a.y, v_mul*a.z, v_mul*a.w);
}
float2 acos_c(float2 a)
{
	const float2 m_i=(float2)(0, -1);
	float2 temp=sq_c(a);
	temp.x-=1;
	return mul_cc(log_c(a+sqrt_c(temp)), m_i);
}
float4 acos_q(float4 a)
{
	const float2 m_i=(float2)(0, -1);
	float4 temp=sq_q(a);
	temp.x-=1;
	return mul_qc(log_q(a+sqrt_q(temp)), m_i);
}
float2 cosh_c(float2 a){return (float2)(cos(a.x)*cosh(a.y), -sin(a.x)*sinh(a.y));}
float4 cosh_q(float4 a)
{
	float4 exp_a=exp_q(a), exp_ma=inv_q(exp_a);
	return mul_qr(exp_a+exp_ma, 0.5f);
}
float2 acosh_c(float2 a)
{
	float2 temp=sq_c(a);
	temp.x-=1;
	return log_c(a+sqrt_c(temp));
}
float4 acosh_q(float4 a)
{
	float4 temp=sq_q(a);
	temp.x-=1;
	return log_q(a+sqrt_q(temp));
}
float2 sin_c(float2 a){return (float2)(sin(a.x)*cosh(a.y), cos(a.x)*sinh(a.y));}
float4 sin_q(float4 a)
{
	float cos_r, sin_r=sincos(a.x, &cos_r);
	float abs_v=sqrt(a.y*a.y+a.z*a.z+a.w*a.w);
	float2 chsh_absv=coshsinh(abs_v);
	float v_mul=-cos_r*chsh_absv.y/abs_v;
	return (float4)(sin_r*chsh_absv.x, v_mul*a.y, v_mul*a.z, v_mul*a.w);
}
float2 asin_c(float2 a)
{
	const float2 i=(float2)(0, 1), m_i=(float2)(0, -1);
	float2 temp=sq_c(a);
	temp.x=1-temp.x, temp.y=-temp.y;
	return mul_cc(log_c(mul_cc(a, i)+sqrt_c(temp)), m_i);
}
float4 asin_q(float4 a)
{
	const float2 i=(float2)(0, 1), m_i=(float2)(0, -1);
	float4 temp=sq_q(a);
	temp=(float4)(1, 0, 0, 0)-temp;
	return mul_qc(log_q(mul_qc(a, i)+sqrt_q(temp)), m_i);
}
float2 sinh_c(float2 a){return (float2)(sinh(a.x)*cos(a.y), cosh(a.x)*sin(a.y));}
float4 sinh_q(float4 a)
{
	float4 exp_a=exp_q(a), exp_ma=inv_q(exp_a);
	return mul_qr(exp_a-exp_ma, 0.5f);
}
float2 asinh_c(float2 a)
{
	float2 temp=sq_c(a);
	temp.x+=1;
	return log_c(a+sqrt_c(temp));
}
float4 asinh_q(float4 a)
{
	float4 temp=sq_q(a);
	temp.x+=1;
	return log_q(a+sqrt_q(temp));
}
float2 tan_c(float2 a)
{
	const float2 _2i=(float2)(0, 2), one=(float2)(1, 0);
	float2 exp2ia=exp_c(mul_cc(_2i, a));
	return div_cc(exp2ia-one, exp2ia+one);
}
float4 tan_q(float4 a)
{
	const float2 _2i=(float2)(0, 2);
	const float4 one=(float4)(1, 0, 0, 0);
	float4 exp2ia=exp_q(mul_cq(_2i, a));
	return div_qq(exp2ia-one, exp2ia+one);
}
float2 atan_c(float2 a)
{
	const float2 i=(float2)(0, 1), i_2=(float2)(0, 0.5f);
	return mul_cc(i_2, log_c(div_cc(i+a, i-a)));
}
float4 atan_q(float4 a)
{
	const float4 i=(float4)(0, 1, 0, 0);
	const float2 i_2=(float2)(0, 0.5f);
	return mul_cq(i_2, log_q(div_qq(i+a, i-a)));
}
float2 tanh_c(float2 a)
{
	const float2 one=(float2)(1, 0);
	float2 exp2a=exp_c(a+a);
	return div_cc(exp2a-one, exp2a+one);
}
float4 tanh_q(float4 a)
{
	const float4 one=(float4)(1, 0, 0, 0);
	float4 exp2a=exp_q(a+a);
	return div_qq(exp2a-one, exp2a+one);
}
float2 atanh_c(float2 a)
{
	const float2 one=(float2)(1, 0);
	return mul_rc(0.5f, log_c(div_cc(one+a, one-a)));
}
float4 atanh_q(float4 a)
{
	const float4 one=(float4)(1, 0, 0, 0);
	return mul_rq(0.5f, log_q(div_qq(one+a, one-a)));
}
//end complex trig & hyp functions

float2	pow_cr(float2 a, float b)
{
	float2 lna=log_c(a);
	float2 temp=mul_cr(lna, b);
	return exp_c(temp);
}
float2	pow_cc(float2 a, float2 b)
{
	float2 lna=log_c(a);
	float2 temp=mul_cc(lna, b);
	return exp_c(temp);
}
float4	pow_cq(float2 a, float4 b)
{
	float2 lna=log_c(a);
	float4 temp=mul_cq(lna, b);
	return exp_q(temp);
}
float4	pow_qr(float4 a, float b)
{
	float4 lna=log_q(a);
	float4 temp=mul_qr(lna, b);
	return exp_q(temp);
}
float4	pow_qc(float4 a, float2 b)
{
	float4 lna=log_q(a);
	float4 temp=mul_qc(lna, b);
	return exp_q(temp);
}
float4	pow_qq(float4 a, float4 b)
{
	float4 lna=log_q(a);
	float4 temp=mul_qq(lna, b);
	return exp_q(temp);
}