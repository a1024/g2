G2_R_RR(hypot){IDX; float a=xr[idx], b=yr[idx]; ASSIGN_R(sqrt(a*a+b*b));}

int mandelbrot(float2 point, int n_iterations)
{
	float rez=0, imz=0, sq_rez=0, sq_imz=0;
	int k=0;
	for(;k<n_iterations&&sq_rez+sq_imz<16;++k)
	{
		imz=rez*imz;//calculate sq(z)
		imz+=imz;
		rez=sq_rez-sq_imz;

		rez+=point.x, imz+=point.y;//add x

		sq_rez=rez*rez, sq_imz=imz*imz;
	}
	return k;
}
G2_R_R(mandelbrot){IDX; ASSIGN_R(mandelbrot((float2)(xr[idx], 0), 200));}
G2_R_C(mandelbrot){IDX; ASSIGN_R(mandelbrot(VEC2(x), 200));}
G2_R_RR(mandelbrot){IDX; ASSIGN_R(mandelbrot((float2)(xr[idx], 0), (int)floor(yr[idx])));}
G2_R_CR(mandelbrot){IDX; ASSIGN_R(mandelbrot(VEC2(x), (int)floor(yr[idx])));}
DISC_R_O(mandelbrot){IDX; disc[idx]=xr[idx]!=xr[idx+offset];}

G2_R_RR(min){IDX; float a=xr[idx], b=yr[idx]; ASSIGN_R((a+b-fabs(a-b))*0.5f);}
G2_C_RC(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_C(a, 0);
	else
		ASSIGN_C(b, yi[idx]);
}
G2_Q_RQ(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_Q(a, 0, 0, 0);
	else
		ASSIGN_Q(b, yi[idx], yj[idx], yk[idx]);
}
G2_C_CR(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_C(a, xi[idx]);
	else
		ASSIGN_C(b, 0);
}
G2_C_CC(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_C(a, xi[idx]);
	else
		ASSIGN_C(b, yi[idx]);
}
G2_Q_CQ(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_Q(a, xi[idx], 0, 0);
	else
		ASSIGN_Q(b, yi[idx], yj[idx], yk[idx]);
}
G2_Q_QR(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_Q(a, xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(b, 0, 0, 0);
}
G2_Q_QC(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_Q(a, xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(b, yi[idx], 0, 0);
}
G2_Q_QQ(min)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a<b)
		ASSIGN_Q(a, xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(b, yi[idx], yj[idx], yk[idx]);
}

G2_R_RR(max){IDX; float a=xr[idx], b=yr[idx]; ASSIGN_R((a+b+fabs(a-b))*0.5f);}
G2_C_RC(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_C(a, 0);
	else
		ASSIGN_C(b, yi[idx]);
}
G2_Q_RQ(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_Q(a, 0, 0, 0);
	else
		ASSIGN_Q(b, yi[idx], yj[idx], yk[idx]);
}
G2_C_CR(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_C(a, xi[idx]);
	else
		ASSIGN_C(b, 0);
}
G2_C_CC(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_C(a, xi[idx]);
	else
		ASSIGN_C(b, yi[idx]);
}
G2_Q_CQ(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_Q(a, xi[idx], 0, 0);
	else
		ASSIGN_Q(b, yi[idx], yj[idx], yk[idx]);
}
G2_Q_QR(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_Q(a, xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(b, 0, 0, 0);
}
G2_Q_QC(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_Q(a, xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(b, yi[idx], 0, 0);
}
G2_Q_QQ(max)
{
	IDX;
	float a=xr[idx], b=yr[idx];
	if(a>b)
		ASSIGN_Q(a, xi[idx], xj[idx], xk[idx]);
	else
		ASSIGN_Q(b, yi[idx], yj[idx], yk[idx]);
}