enum DataPointIdx
{
	P_UNW,			P_UNE,
			P_Umm,
	P_USW,			P_USE,

			P_mNm,
	P_mmW,	P_mmm,	P_mmE,
			P_mSm,

	P_DNW,			P_DNE,
			P_Dmm,
	P_DSW,			P_DSE,
};
__constant unsigned char e2v1[54]=//edge index to first vertex (data point) index
{
	//main (new) edge indices
	P_DSW, P_DSW, P_DSW,

	P_mmW, P_mmW, P_mmW, P_mmW,
	P_mSm, P_mSm, P_mSm, P_mSm,
	P_Dmm, P_Dmm, P_Dmm, P_Dmm,

	P_mmW, P_mmW, P_mmW, P_mmW,
	P_Umm, P_mNm, P_Dmm, P_mSm,
	P_mmE, P_mmE, P_mmE, P_mmE,

	P_mmm, P_mmm, P_mmm, P_mmm, P_mmm, P_mmm,

	//extended (redundant) edge indices
	P_DSE, P_DNE, P_DSE, P_DNE, P_DNW, P_USW, P_USE, P_UNE, P_UNW,
	P_mmE, P_mmE, P_mmE, P_mmE,
	P_mNm, P_mNm, P_mNm, P_mNm,
	P_Umm, P_Umm, P_Umm, P_Umm,
};
__constant unsigned char e2v2[54]=//edge index to second vertex (data point) index
{
	//main (new) edge indices
	P_DSE, P_DNW, P_USW,

	P_DSW, P_USW, P_UNW, P_DNW,
	P_DSW, P_USW, P_USE, P_DSE,
	P_DSW, P_DSE, P_DNE, P_DNW,

	P_Dmm, P_mSm, P_Umm, P_mNm,
	P_mNm, P_Dmm, P_mSm, P_Umm,
	P_Dmm, P_mSm, P_Umm, P_mNm,

	P_Dmm, P_mSm, P_mmW, P_Umm, P_mNm, P_mmE,

	//extended (redundant) edge indices
	P_DNE, P_DNW, P_USE, P_UNE, P_UNW, P_USE, P_UNE, P_UNW, P_USW,
	P_DSE, P_USE, P_UNE, P_DNE,
	P_DNW, P_UNW, P_UNE, P_DNE,
	P_USW, P_USE, P_UNE, P_UNW,
};
#define			_pi		3.14159265358979f
#define			_sqrt3	1.73205080756888f
float			clamp01(float x){return x<0?0:x>1?1:x;}
float4			solvequadratic(float a, float b, float c)//returns (r1, i1, r2, i2)
{
	const float lim=1e6f;//float
//	const float lim=1e10f;//double
	float4 ret;
	if(a==0)
		ret=(float4)(-c/b, 0, INFINITY, 0);
	else if(fabs(b)/fabs(a)>=lim&&fabs(c)/fabs(a)>=lim)
		ret=(float4)(-c/b, 0, -b/a, 0);
	else
	{
		float first=-b, disc2=b*b-4*a*c;
		float disc=sqrt(fabs(disc2));
		float inv_den=0.5f/a;
		if(disc2<0)
		{
			first*=inv_den, disc*=inv_den;
			ret=(float4)(first, disc, first, -disc);
		}
		else
			ret=(float4)((first+disc)*inv_den, 0, (first-disc)*inv_den, 0);
	}
	return ret;
}
float2			solvecubic(float a, float b, float c, float d)
{
	float2 ret;
	float r1, r2r, r2i, r3r, r3i;
	//http://easycalculation.com/algebra/learn-cubic-equation.php
	//http://stackoverflow.com/questions/13328676/c-solving-cubic-equations
	if(a==0)
	{
		r1=INFINITY;
		float4 ret2=solvequadratic(b, c, d);
		r2r=ret2.x, r2i=ret2.y, r3r=ret2.z, r3i=ret2.w;
	}
	else if(d==0)
	{
		r1=0;
		float4 ret2=solvequadratic(a, b, c);
		r2r=ret2.x, r2i=ret2.y, r3r=ret2.z, r3i=ret2.w;
	}
	else
	{
		float inv_a=1/a, disc, q, r, dum1, s, t, term1, r13;
		b*=inv_a, c*=inv_a, d*=inv_a;
		q=(3*c-b*b)/9;
		r=-27*d+b*(9*c-2*b*b);
		r/=54;
		disc=q*q*q+r*r;
		term1=b/3;
		if(disc>0)
		{
			s=r+sqrt(disc);
			s=cbrt(s);
			t=r-sqrt(disc);
			t=cbrt(t);
			r1=-term1+s+t;//The first root is always real
			term1+=(s+t)/2;
			float term2=_sqrt3*(-t+s)/2;
			r2r=r3r=-term1, r2i=term2, r3i=-term2;
		}
		else if(disc==0)//The remaining options are all real
		{
			r13=cbrt(r);
			r1=-term1+2*r13;
			r3r=r2r=-(r13+term1);//at least two are equal
			r2i=r3i=0;
		}
		else//Only option left is that all roots are real and unequal (to get here, q < 0)
		{
			q=-q;
			dum1=q*q*q;
			dum1=acos(r/sqrt(dum1));
			r13=2*sqrt(q);
			r1=-term1+r13*cos(dum1/3);
			r2r=-term1+r13*cos((dum1+2*_pi)/3), r2i=0;
			r3r=-term1+r13*cos((dum1+4*_pi)/3), r3i=0;
		}
	}
	const float tol=1e-5f;//tolerance
		 if(r1 >=0		&&r1 <=1)		ret.x=r1, ret.y=1;
	else if(r2r>=0		&&r2r<=1)		ret.x=r2r, ret.y=1;//what
	else if(r3r>=0		&&r3r<=1)		ret.x=r3r, ret.y=1;
	else if(r1 >=-tol	&&r1 <=1+tol)	ret.x=r1, ret.y=1;
	else if(r2r>=-tol	&&r2r<=1+tol)	ret.x=r2r, ret.y=1;//what
	else if(r3r>=-tol	&&r3r<=1+tol)	ret.x=r3r, ret.y=1;
	else if(r1 >=-tol	&&r1 <=1+tol)	ret.x=r1, ret.y=1;
	else if(r2r>=-tol	&&r2r<=1+tol)	ret.x=r2r, ret.y=1;//what
	else if(r3r>=-tol	&&r3r<=1+tol)	ret.x=r3r, ret.y=1;
	else								ret.x=-1, ret.y=0;
	return ret;
}
__kernel void ti3d_zerocross(__constant int *size, __constant float *coeffs, __constant float *ndr, __constant ulong *workidx, __global float *vertices)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize, nvert_total}, coeffs: {Xstart, Xsample, Ystart, Ysample, Zstart, Zsample, isovalue}
	int id=get_global_id(0), nvert_total=size[4];
	if(id<nvert_total)
	{
	//	int kv=id*6;	//DEBUG
	//	vertices[kv  ]=42;
	//	vertices[kv+1]=42;
	//	vertices[kv+2]=42;
	//	vertices[kv+3]=42;
	//	vertices[kv+4]=42;
	//	vertices[kv+5]=42;

	//	float4 coords=(float4)(1, 2, 3, 0), normal=(float4)(4, 5, 6, 0);	//DEBUG
	//	vertices[kv  ]=coords.x;
	//	vertices[kv+1]=coords.y;
	//	vertices[kv+2]=coords.z;
	//	vertices[kv+3]=normal.x;
	//	vertices[kv+4]=normal.y;
	//	vertices[kv+5]=normal.z;

		int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2], ndrSize=size[3], XYplaces=Xplaces*Yplaces;
		ulong wi=workidx[id];
		int kx=(ushort)(wi>>16), ky=(ushort)(wi>>32), kz=(ushort)(wi>>48), ke=(ushort)wi;
		int idx=Xplaces*(Yplaces*kz+ky)+kx;
		float
			Xstart=coeffs[0], Xsample=coeffs[1],
			Ystart=coeffs[2], Ysample=coeffs[3],
			Zstart=coeffs[4], Zsample=coeffs[5], th=coeffs[6];
		float
			X_W=Xstart+kx*Xsample, X_m=X_W+Xsample*0.5f, X_E=X_W+Xsample,//space coordinates
			Y_S=Ystart+ky*Ysample, Y_m=Y_S+Ysample*0.5f, Y_N=Y_S+Ysample,
			Z_D=Zstart+kz*Zsample, Z_m=Z_D+Zsample*0.5f, Z_U=Z_D+Zsample;
		float
			V_UNW=ndr[idx+XYplaces+Xplaces],	V_UNE=ndr[idx+XYplaces+Xplaces+1],
			V_USW=ndr[idx+XYplaces],			V_USE=ndr[idx+XYplaces+1],

			V_DNW=ndr[idx+Xplaces],				V_DNE=ndr[idx+Xplaces+1],
			V_DSW=ndr[idx],						V_DSE=ndr[idx+1];
		float
			V_Dmm=ndr[ndrSize+( idx				<<2)  ],
			V_mSm=ndr[ndrSize+( idx				<<2)+1],
			V_mmW=ndr[ndrSize+( idx				<<2)+2],
			V_mmm=ndr[ndrSize+( idx				<<2)+3],
			V_mmE=ndr[ndrSize+((idx+1)			<<2)+2],
			V_mNm=ndr[ndrSize+((idx+Xplaces)	<<2)+1],
			V_Umm=ndr[ndrSize+((idx+XYplaces)	<<2)  ];
		float4 dp[15]=//data points
		{
			{X_W, Y_N, Z_U, V_UNW},							{X_E, Y_N, Z_U, V_UNE},	//	0 P_UNW			1 P_UNE
									{X_m, Y_m, Z_U, V_Umm},							//			2 P_Umm
			{X_W, Y_S, Z_U, V_USW},							{X_E, Y_S, Z_U, V_USE},	//	3 P_USW			4 P_USE

									{X_m, Y_N, Z_m, V_mNm},							//			5 P_mNm
			{X_W, Y_m, Z_m, V_mmW}, {X_m, Y_m, Z_m, V_mmm}, {X_E, Y_m, Z_m, V_mmE},	//	6 P_mmW	7 P_mmm	8 P_mmE
									{X_m, Y_S, Z_m, V_mSm},							//			9 P_mSm

			{X_W, Y_N, Z_D, V_DNW},							{X_E, Y_N, Z_D, V_DNE},	//	10 P_DNW		11 P_DNE
									{X_m, Y_m, Z_D, V_Dmm},							//			12 P_Dmm
			{X_W, Y_S, Z_D, V_DSW},							{X_E, Y_S, Z_D, V_DSE},	//	13 P_DSW		14 P_DSE
		};
		float4 A=dp[e2v1[ke]], B=dp[e2v2[ke]];
		float
			Xd=B.x-A.x, X_aE=X_E-A.x, X_Wa=A.x-X_W,//see trilinear interpolation
			Yd=B.y-A.y, Y_aN=Y_N-A.y, Y_Sa=A.y-Y_S,
			Zd=B.z-A.z, Z_aU=Z_U-A.z, Z_Da=A.z-Z_D,
		
			A_DS=(V_DSE-V_DSW)*Xd, B_DS=X_aE*V_DSW+X_Wa*V_DSE,
			A_DN=(V_DNE-V_DNW)*Xd, B_DN=X_aE*V_DNW+X_Wa*V_DNE,
			A_US=(V_USE-V_USW)*Xd, B_US=X_aE*V_USW+X_Wa*V_USE,
			A_UN=(V_UNE-V_UNW)*Xd, B_UN=X_aE*V_UNW+X_Wa*V_UNE,
		
			C_D=(A_DN-A_DS)*Yd, D_D=(B_DN-B_DS)*Yd+Y_aN*A_DS+Y_Sa*A_DN, E_D=Y_aN*B_DS+Y_Sa*B_DN,
			C_U=(A_UN-A_US)*Yd, D_U=(B_UN-B_US)*Yd+Y_aN*A_US+Y_Sa*A_UN, E_U=Y_aN*B_US+Y_Sa*B_UN,
		
			a=(C_U-C_D)*Zd, b=(D_U-D_D)*Zd+Z_aU*C_D+Z_Da*C_U, c=(E_U-E_D)*Zd+Z_aU*D_D+Z_Da*D_U, d=Z_aU*E_D+Z_Da*E_U-(X_E-X_W)*(Y_N-Y_S)*(Z_U-Z_D)*th;
		float2 ret;
		if(a==0&&b==0&&c==0)//degenerate equation d=0
			ret=(float2)((th-A.w)/(B.w-A.w), 1);
		else
			ret=solvecubic(a, b, c, d);
	//	float3 pa={A.x, A.y, A.z}, pb={B.x, B.y, B.z};
		if(!ret.y)//solvecubic failed, use linear interpolation
			ret=(float2)((th-A.w)/(B.w-A.w), 1);
		float4 coords=mix(A, B, ret.x);
	//	float4 coords=mix(pa, pb, ret.x);
	//	vstore3(coords, id*6, vertices);
		int kv=id*6;
		vertices[kv  ]=coords.x;
		vertices[kv+1]=coords.y;
		vertices[kv+2]=coords.z;

		float4 normal=(float4)(0, 0, 0, 0);//calculate normal
		{
			float4 a=(float4)((coords.x-X_W)/Xsample, (coords.y-Y_S)/Ysample, (coords.z-Z_D)/Zsample, 0);
			float
				V_DSW_DSE=V_DSE-V_DSW,
				V_DNW_DNE=V_DNE-V_DNW,
				V_USW_USE=V_USE-V_USW,
				V_UNW_UNE=V_UNE-V_UNW,
				V_DS=V_DSW+V_DSW_DSE*a.x,
				V_DN=V_DNW+V_DNW_DNE*a.x,
				V_US=V_USW+V_USW_USE*a.x,
				V_UN=V_UNW+V_UNW_UNE*a.x,
				TY=V_DSW_DSE+(V_DNW_DNE-V_DSW_DSE)*a.y,
				V_D_SN=V_DN-V_DS, V_U_SN=V_UN-V_US;
			normal=(float4)
			(
				TY+((V_USW_USE+(V_UNW_UNE-V_USW_USE)*a.y)-TY)*a.z,
				V_D_SN+(V_U_SN-V_D_SN)*a.z,
				(V_US+V_U_SN*a.y)-(V_DS+V_D_SN*a.y),
				0
			);
			float inv_abs=sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
			if(inv_abs)
			{
				inv_abs=1/inv_abs;
				normal.x*=inv_abs, normal.y*=inv_abs, normal.z*=inv_abs;
			}
		}
	//	vstore3(normal, id*6+3, vertices);
		vertices[kv+3]=normal.x;
		vertices[kv+4]=normal.y;
		vertices[kv+5]=normal.z;	//*/
	}
}