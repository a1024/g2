__kernel void ti3d_av(__constant int *size, __constant float *result, __global float *ndr)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize}
	int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2];
	int kx=get_global_id(0), ky=get_global_id(1), kz=get_global_id(2);
	if(kx<Xplaces-1&&ky<Yplaces-1&&kz<Zplaces-1)
	{
		int ndrSize=size[3], XYplaces=Xplaces*Yplaces, idx=Xplaces*(Yplaces*kz+ky)+kx;
		float
			V_UNW=result[idx+XYplaces+Xplaces	],	V_UNE=result[idx+XYplaces+Xplaces	+1],
			V_USW=result[idx+XYplaces			],	V_USE=result[idx+XYplaces			+1],

			V_DNW=result[idx+Xplaces			],	V_DNE=result[idx+Xplaces			+1],
			V_DSW=result[idx					],	V_DSE=result[idx					+1];

		ndr[idx+XYplaces+Xplaces	]=V_UNW;	ndr[idx+XYplaces+Xplaces	+1]=V_UNE;
		ndr[idx+XYplaces			]=V_USW;	ndr[idx+XYplaces			+1]=V_USE;

		ndr[idx+Xplaces	]=V_DNW;				ndr[idx+Xplaces	+1]=V_DNE;
		ndr[idx			]=V_DSW;				ndr[idx			+1]=V_DSE;
		ndr[ndrSize+(idx<<2)  ]=(V_DNW+V_DSW+V_DNE+V_DSE)*0.25f;//V_Dmm
		ndr[ndrSize+(idx<<2)+1]=(V_DSW+V_USW+V_USE+V_DSE)*0.25f;//V_mSm
		ndr[ndrSize+(idx<<2)+2]=(V_UNW+V_DNW+V_DSW+V_USW)*0.25f;//V_mmW
		ndr[ndrSize+(idx<<2)+3]=(V_UNW+V_DNW+V_DSW+V_USW+V_UNE+V_DNE+V_USE+V_DSE)*0.125f;//V_mmm
	}
}
__kernel void ti3d_av_east(__constant int *size, __constant float *result, __global float *ndr)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize}
	int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2];
	int kx=Xplaces-1, ky=get_global_id(1), kz=get_global_id(2);
	if(ky<Yplaces-1&&kz<Zplaces-1)
	{
		int ndrSize=size[3], XYplaces=Xplaces*Yplaces, idx=Xplaces*(Yplaces*kz+ky)+kx;
		float
			V_UNW=result[idx+XYplaces+Xplaces],
			V_USW=result[idx+XYplaces],

			V_DNW=result[idx+Xplaces],
			V_DSW=result[idx];

		ndr[idx+XYplaces+Xplaces]=V_UNW;
		ndr[idx+XYplaces		]=V_USW;

		ndr[idx+Xplaces			]=V_DNW;
		ndr[idx					]=V_DSW;
		ndr[ndrSize+(idx<<2)+2]=(V_UNW+V_DNW+V_DSW+V_USW)*0.25f;//V_mmW
	}
}
__kernel void ti3d_av_north(__constant int *size, __constant float *result, __global float *ndr)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize}
	int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2];
	int kx=get_global_id(0), ky=Yplaces-1, kz=get_global_id(2);
	if(kx<Xplaces-1&&kz<Zplaces-1)
	{
		int ndrSize=size[3], XYplaces=Xplaces*Yplaces, idx=Xplaces*(Yplaces*kz+ky)+kx;
		float
			V_USW=result[idx+XYplaces],			V_USE=result[idx+XYplaces+1],

			V_DSW=result[idx],					V_DSE=result[idx+1];

		ndr[idx+XYplaces]=V_USW;			ndr[idx+XYplaces+1]=V_USE;

		ndr[idx			]=V_DSW;			ndr[idx			+1]=V_DSE;
		ndr[ndrSize+(idx<<2)+1]=(V_DSW+V_USW+V_USE+V_DSE)*0.25f;//V_mSm
	}
}
__kernel void ti3d_av_up(__constant int *size, __constant float *result, __global float *ndr)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize}
	int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2];
	int kx=get_global_id(0), ky=get_global_id(1), kz=Zplaces-1;
	if(kx<Xplaces-1&&ky<Yplaces-1)
	{
		int ndrSize=size[3], idx=Xplaces*(Yplaces*kz+ky)+kx;
		float
			V_DNW=result[idx+Xplaces],				V_DNE=result[idx+Xplaces+1],
			V_DSW=result[idx		],				V_DSE=result[idx		+1];

		ndr[idx+Xplaces	]=V_DNW;				ndr[idx+Xplaces	+1]=V_DNE;
		ndr[idx			]=V_DSW;				ndr[idx			+1]=V_DSE;
		ndr[ndrSize+(idx<<2)  ]=(V_DNW+V_DSW+V_DNE+V_DSE)*0.25f;//V_Dmm
	}
}
__constant unsigned char nt[16]=//number of triangles produced from tetrahedron given 4 bit vertex condition
{//	00	01	10	11	lo/hi
	0,	1,	1,	2,//00
	1,	2,	2,	1,//01
	1,	2,	2,	1,//10
	2,	1,	1,	0,//11
};
int hammingweight(ulong x)//https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
{
	x-=x>>1&0x5555555555555555;
	x=(x&0x3333333333333333)+(x>>2&0x3333333333333333);
	return ((x+(x>>4))&0x0F0F0F0F0F0F0F0F)*0x0101010101010101>>56;//0xc1200000, 0xc0f00000

//	x=(x&0x5555555555555555ULL)+(x>> 1&0x5555555555555555ULL);
//	x=(x&0x3333333333333333ULL)+(x>> 2&0x3333333333333333ULL);
//	x=(x&0x0F0F0F0F0F0F0F0FULL)+(x>> 4&0x0F0F0F0F0F0F0F0FULL);
//	x=(x&0x00FF00FF00FF00FFULL)+(x>> 8&0x00FF00FF00FF00FFULL);
//	x=(x&0x0000FFFF0000FFFFULL)+(x>>16&0x0000FFFF0000FFFFULL);//0xc1200000, 0xc0f00000
//	x=(uint)x+(uint)(x>>32);
//	return (int)x;

//	int sum=0;
//	for(int k=0;k<64;++k)
//		sum+=x>>k&1;
//	return sum;
}
enum ExEdgeBitIdx//54edges
{
	//main (new) edge indices
	DSW_DSE, DSW_DNW, DSW_USW,

	DSW_mmW, USW_mmW, UNW_mmW, DNW_mmW,
	DSW_mSm, USW_mSm, USE_mSm, DSE_mSm,
	DSW_Dmm, DSE_Dmm, DNE_Dmm, DNW_Dmm,

	mmW_Dmm, mmW_mSm, mmW_Umm, mmW_mNm,
	Umm_mNm, mNm_Dmm, Dmm_mSm, mSm_Umm,
	mmE_Dmm, mmE_mSm, mmE_Umm, mmE_mNm,

	mmm_Dmm, mmm_mSm, mmm_mmW, mmm_Umm, mmm_mNm, mmm_mmE,

	//extended (redundant) edge indices
	DSE_DNE, DNE_DNW, DSE_USE, DNE_UNE, DNW_UNW, USW_USE, USE_UNE, UNE_UNW, UNW_USW,
	mmE_DSE, mmE_USE, mmE_UNE, mmE_DNE,
	mNm_DNW, mNm_UNW, mNm_UNE, mNm_DNE,
	Umm_USW, Umm_USE, Umm_UNE, Umm_UNW,
};
enum ExEdgeBitIdxRev//54edges
{
	//main (new) edge indices
	DSE_DSW, DNW_DSW, USW_DSW,

	mmW_DSW, mmW_USW, mmW_UNW, mmW_DNW,
	mSm_DSW, mSm_USW, mSm_USE, mSm_DSE,
	Dmm_DSW, Dmm_DSE, Dmm_DNE, Dmm_DNW,

	Dmm_mmW, mSm_mmW, Umm_mmW, mNm_mmW,
	mNm_Umm, Dmm_mNm, mSm_Dmm, Umm_mSm,
	Dmm_mmE, mSm_mmE, Umm_mmE, mNm_mmE,

	Dmm_mmm, mSm_mmm, mmW_mmm, Umm_mmm, mNm_mmm, mmE_mmm,

	//extended (redundant) edge indices
	DNE_DSE, DNW_DNE, USE_DSE, UNE_DNE, UNW_DNW, USE_USW, UNE_USE, UNW_UNE, USW_UNW,
	DSE_mmE, USE_mmE, UNE_mmE, DNE_mmE,
	DNW_mNm, UNW_mNm, UNE_mNm, DNE_mNm,
	USW_Umm, USE_Umm, UNE_Umm, UNW_Umm,
};
__constant ulong obsmask[8]=//original bit selection mask
{
	0x00000001FFFFFFFF,//all coords inside
	(1ULL<<DSE_USE)|(1ULL<<DSE_DNE)|(1ULL<<mmE_DSE)|(1ULL<<mmE_USE)|(1ULL<<mmE_UNE)|(1ULL<<mmE_DNE)|0x00000001FFFFFFFF,//x at end: DSE's originals
	(1ULL<<DNW_UNW)|(1ULL<<DNW_DNE)|(1ULL<<mNm_DNW)|(1ULL<<mNm_UNW)|(1ULL<<mNm_UNE)|(1ULL<<mNm_DNE)|0x00000001FFFFFFFF,//y at end: DNW's originals
	(1ULL<<DNE_UNE) | (1ULL<<DSE_USE)|(1ULL<<DSE_DNE)|(1ULL<<mmE_DSE)|(1ULL<<mmE_USE)|(1ULL<<mmE_UNE)|(1ULL<<mmE_DNE) | (1ULL<<DNW_UNW)|(1ULL<<DNW_DNE)|(1ULL<<mNm_DNW)|(1ULL<<mNm_UNW)|(1ULL<<mNm_UNE)|(1ULL<<mNm_DNE)|0x00000001FFFFFFFF,//x&y at end: (DSE, DNW & DNE)'s originals

	(1ULL<<USW_USE)|(1ULL<<USW_UNW)|(1ULL<<Umm_USW)|(1ULL<<Umm_USE)|(1ULL<<Umm_UNE)|(1ULL<<Umm_UNW)|0x00000001FFFFFFFF,//z at end: USW's originals
	(1ULL<<USE_UNE) | (1ULL<<DSE_USE)|(1ULL<<DSE_DNE)|(1ULL<<mmE_DSE)|(1ULL<<mmE_USE)|(1ULL<<mmE_UNE)|(1ULL<<mmE_DNE) | (1ULL<<USW_USE)|(1ULL<<USW_UNW)|(1ULL<<Umm_USW)|(1ULL<<Umm_USE)|(1ULL<<Umm_UNE)|(1ULL<<Umm_UNW)|0x00000001FFFFFFFF,//x&z at end: (DSE, USW & USE)'s originals
	(1ULL<<UNW_UNE) | (1ULL<<DNW_UNW)|(1ULL<<DNW_DNE)|(1ULL<<mNm_DNW)|(1ULL<<mNm_UNW)|(1ULL<<mNm_UNE)|(1ULL<<mNm_DNE) | (1ULL<<USW_USE)|(1ULL<<USW_UNW)|(1ULL<<Umm_USW)|(1ULL<<Umm_USE)|(1ULL<<Umm_UNE)|(1ULL<<Umm_UNW)|0x00000001FFFFFFFF,//y&z at end: (DNW, USW & UNW)'s originals
	0x003FFFFFFFFFFFFF,//x,y&z at end: entire cube
};
__kernel void ti3d_classifyedges(__constant int *size, __constant float *coeffs, __constant float *ndr, __global ulong *edgeinfo, __global int *nvert, __global int *ntrgl)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize}, coeffs: {Xstart, Xsample, Ystart, Ysample, Zstart, Zsample, isovalue} (just for isovalue)
	int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2];
	int kx=get_global_id(0), ky=get_global_id(1), kz=get_global_id(2);
	if(kx<Xplaces-1&&ky<Yplaces-1&&kz<Zplaces-1)
	{
		float th=coeffs[6];//isovalue
		int ndrSize=size[3], XYplaces=Xplaces*Yplaces;
		int idx=Xplaces*(Yplaces*kz+ky)+kx;
		//all vertices relevant to current cube
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
		char
			C_UNW=V_UNW>th, C_UNE=V_UNE>th,
			C_USW=V_USW>th, C_USE=V_USE>th,
			C_DNW=V_DNW>th, C_DNE=V_DNE>th,
			C_DSW=V_DSW>th, C_DSE=V_DSE>th,

			C_Umm=V_Umm>th, C_mNm=V_mNm>th,
			C_mmW=V_mmW>th, C_mmE=V_mmE>th,
			C_mSm=V_mSm>th, C_Dmm=V_Dmm>th,
			C_mmm=V_mmm>th;
		char ec[54]=//edge conditions
		{
			C_DSW!=C_DSE, C_DSW!=C_DNW, C_DSW!=C_USW,

			C_mmW!=C_DSW, C_mmW!=C_USW, C_mmW!=C_UNW, C_mmW!=C_DNW,
			C_mSm!=C_DSW, C_mSm!=C_USW, C_mSm!=C_USE, C_mSm!=C_DSE,
			C_Dmm!=C_DSW, C_Dmm!=C_DSE, C_Dmm!=C_DNE, C_Dmm!=C_DNW,

			C_mmW!=C_Dmm, C_mmW!=C_mSm, C_mmW!=C_Umm, C_mmW!=C_mNm,
			C_Umm!=C_mNm, C_mNm!=C_Dmm, C_Dmm!=C_mSm, C_mSm!=C_Umm,
			C_mmE!=C_Dmm, C_mmE!=C_mSm, C_mmE!=C_Umm, C_mmE!=C_mNm,

			C_mmm!=C_Dmm, C_mmm!=C_mSm, C_mmm!=C_mmW, C_mmm!=C_Umm, C_mmm!=C_mNm, C_mmm!=C_mmE,

			C_DSE!=C_DNE, C_DNE!=C_DNW, C_DSE!=C_USE, C_DNE!=C_UNE, C_DNW!=C_UNW, C_USW!=C_USE, C_USE!=C_UNE, C_UNE!=C_UNW, C_UNW!=C_USW,
			C_mmE!=C_DSE, C_mmE!=C_USE, C_mmE!=C_UNE, C_mmE!=C_DNE, 
			C_mNm!=C_DNW, C_mNm!=C_UNW, C_mNm!=C_UNE, C_mNm!=C_DNE, 
			C_Umm!=C_USW, C_Umm!=C_USE, C_Umm!=C_UNE, C_Umm!=C_UNW,
		};
#define		EC(number)	((ulong)ec[number]<<number)
		edgeinfo[idx]=
			EC( 0)|EC( 1)|EC( 2)|EC( 3)|EC( 4)|EC( 5)|EC( 6)|EC( 7)|EC( 8)|
			EC( 9)|EC(10)|EC(11)|EC(12)|EC(13)|EC(14)|EC(15)|EC(16)|EC(17)|
			EC(18)|EC(19)|EC(20)|EC(21)|EC(22)|EC(23)|EC(24)|EC(25)|EC(26)|
			EC(27)|EC(28)|EC(29)|EC(30)|EC(31)|EC(32)|EC(33)|EC(34)|EC(35)|
			EC(36)|EC(37)|EC(38)|EC(39)|EC(40)|EC(41)|EC(42)|EC(43)|EC(44)|
			EC(45)|EC(46)|EC(47)|EC(48)|EC(49)|EC(50)|EC(51)|EC(52)|EC(53);
#undef		EC
		//if(edgeinfo[idx]&&kx==3&&ky==0&&kz==0)
		//{
		//	printf("th=%f\n\n", th);
		//	printf("V_UNW=%f, V_UNE=%f,\nV_USW=%f, V_USE=%f,\n\nV_DNW=%f, V_DNE=%f,\nV_DSW=%f, V_DSE=%f\n\n", V_UNW, V_UNE, V_USW, V_USE, V_DNW, V_DNE, V_DSW, V_DSE);
		//	printf("V_Dmm=%f,\nV_mSm=%f,\nV_mmW=%f,\nV_mmm=%f,\nV_mmE=%f,\nV_mNm=%f,\nV_Umm=%f\n\n", V_Dmm, V_mSm, V_mmW, V_mmm, V_mmE, V_mNm, V_Umm);
		//	printf("C_UNW=%d, C_UNE=%d,\nC_USW=%d, C_USE=%d,\nC_DNW=%d, C_DNE=%d,\nC_DSW=%d, C_DSE=%d,\n\nC_Umm=%d, C_mNm=%d,\nC_mmW=%d, C_mmE=%d,\nC_mSm=%d, C_Dmm=%d,\nC_mmm=%d\n\n",
		//		(int)C_UNW, (int)C_UNE,
		//		(int)C_USW, (int)C_USE,
		//		(int)C_DNW, (int)C_DNE,
		//		(int)C_DSW, (int)C_DSE,
		//
		//		(int)C_Umm, (int)C_mNm,
		//		(int)C_mmW, (int)C_mmE,
		//		(int)C_mSm, (int)C_Dmm,
		//		(int)C_mmm);
		//	printf("edgeinfo[%d] = 0x%016LLX\n\n", idx, edgeinfo[idx]);
		//}
		ulong mask=obsmask[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)];
	//	nvert[idx]=idx;
	//	nvert[idx]=hammingweight(0xc0f00000c0f00000&mask);//
		nvert[idx]=hammingweight(edgeinfo[idx]&mask);
#define		NT(V0, V1, V2, V3)		nt[V3<<3|V2<<2|V1<<1|V0]
		ntrgl[idx]=
			NT(C_USW, C_DSW, C_mmW, C_mSm)+
			NT(C_UNW, C_USW, C_mmW, C_Umm)+
			NT(C_DNW, C_UNW, C_mmW, C_mNm)+
			NT(C_DSW, C_DNW, C_mmW, C_Dmm)+

			NT(C_USE, C_DSE, C_mmE, C_mSm)+
			NT(C_UNE, C_USE, C_mmE, C_Umm)+
			NT(C_DNE, C_UNE, C_mmE, C_mNm)+
			NT(C_DSE, C_DNE, C_mmE, C_Dmm)+

			NT(C_DNE, C_DNW, C_mNm, C_Dmm)+
			NT(C_DSE, C_DSW, C_Dmm, C_mSm)+
			NT(C_USE, C_USW, C_mSm, C_Umm)+
			NT(C_UNE, C_UNW, C_Umm, C_mNm)+

			NT(C_UNW, C_mmW, C_Umm, C_mNm)+
			NT(C_USW, C_Umm, C_mmW, C_mSm)+
			NT(C_USE, C_mSm, C_mmE, C_Umm)+
			NT(C_UNE, C_mmE, C_Umm, C_mNm)+

			NT(C_DNW, C_mmW, C_Dmm, C_mNm)+
			NT(C_DSW, C_Dmm, C_mmW, C_mSm)+
			NT(C_DSE, C_mSm, C_mmE, C_Dmm)+
			NT(C_DNE, C_mmE, C_Dmm, C_mNm)+

			NT(C_mmm, C_mmW, C_Umm, C_mNm)+
			NT(C_mmm, C_Umm, C_mmW, C_mSm)+
			NT(C_mmm, C_mSm, C_mmE, C_Umm)+
			NT(C_mmm, C_mmE, C_Umm, C_mNm)+

			NT(C_mmm, C_mmW, C_Dmm, C_mNm)+
			NT(C_mmm, C_Dmm, C_mmW, C_mSm)+
			NT(C_mmm, C_mSm, C_mmE, C_Dmm)+
			NT(C_mmm, C_mmE, C_Dmm, C_mNm);
#undef		NT
	}
}