G2_R_R(abs){IDX; rr[idx]=fabs(xr[idx]);}
G2_R_C(abs){IDX; float2 a=VEC2(x); rr[idx]=sqrt(a.x*a.x+a.y*a.y);}
G2_R_Q(abs){IDX; float4 a=VEC4(x); rr[idx]=sqrt(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);}

G2_R_R(arg)
{
	IDX;
	if(xr[idx]<0)
		rr[idx]=_pi;
	else if(xr[idx]==0)
		rr[idx]=NAN;
	else
		rr[idx]=0;
}
G2_R_C(arg){IDX; rr[idx]=atan2(xi[idx], xr[idx]);}
G2_R_Q(arg)
{
	IDX;
	float4 a=VEC4(x);
	float abs_a=a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w;
	if(!abs_a)
		rr[idx]=NAN;
	else
	{
		abs_a=sqrt(abs_a);
		rr[idx]=acos(a.x/abs_a);
	}
}
DISC_R_I(arg)
{
	IDX;
		 if(xr[idx]<0)	disc[idx]=xr[idx+offset]>=0;
	else if(xr[idx]>0)	disc[idx]=xr[idx+offset]<=0;
	else				disc[idx]=xr[idx+offset]!=0;
}
DISC_C_I(arg){IDX; disc[idx]=false;}//TODO
DISC_Q_I(arg){IDX; disc[idx]=false;}//TODO

G2_R_C(real){IDX; rr[idx]=xr[idx];}

G2_R_C(imag){IDX; rr[idx]=xi[idx];}

G2_C_C(conjugate){IDX; ASSIGN_C(xr[idx], -xi[idx]);}
G2_Q_Q(conjugate){IDX; ASSIGN_Q(xr[idx], -xi[idx], -xj[idx], -xk[idx]);}

G2_C_R(polar)
{
	IDX;
	float a=xr[idx];
	float2 ret=(float2)(fabs(a), a<0?_pi:a==0?NAN:0);
	RET_C;
}
G2_C_C(polar)
{
	IDX;
	float2 a=VEC2(x);
	ASSIGN_C(sqrt(a.x*a.x+a.y*a.y), atan2(a.y, a.x));
}
G2_C_Q(polar)
{
	IDX;
	float4 a=VEC4(x);
	rr[idx]=sqrt(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);
	ri[idx]=acos(a.x/rr[idx]);
}
DISC_R_I(polar)
{
	IDX;
		 if(xr[idx]<0)	disc[idx]=xr[idx+offset]>=0;
	else if(xr[idx]>0)	disc[idx]=xr[idx+offset]<=0;
	else				disc[idx]=xr[idx+offset]!=0;
}
DISC_C_I(polar){IDX; disc[idx]=false;}//TODO
DISC_Q_I(polar){IDX; disc[idx]=false;}//TODO

G2_C_C(cartesian)
{
	IDX;
	float2 a=VEC2(x);
	ASSIGN_C(a.x*cos(a.y), a.x*sin(a.y));
}
G2_Q_Q(cartesian)
{
	IDX;
	float r=xr[idx], i=xi[idx], j=xj[idx], k=xk[idx];
	float cos_j=cos(j), r_cos_k=r*cos(k);
	rr[idx]=cos(i)*cos_j*r_cos_k;
	ri[idx]=sin(i)*cos_j*r_cos_k;
	rj[idx]=sin(j)*r_cos_k;
	rk[idx]=r*sin(k);
}