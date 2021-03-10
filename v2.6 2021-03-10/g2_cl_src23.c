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
#define		TEHE(V0, V1, V2, V3)		V0##_##V1, V0##_##V2, V0##_##V3, V1##_##V2, V1##_##V3, V2##_##V3
__constant unsigned char ti[28*6]=//tetrahedron edge indices
{
	TEHE(USW, DSW, mmW, mSm),
	TEHE(UNW, USW, mmW, Umm),
	TEHE(DNW, UNW, mmW, mNm),
	TEHE(DSW, DNW, mmW, Dmm),

	TEHE(USE, DSE, mmE, mSm),
	TEHE(UNE, USE, mmE, Umm),
	TEHE(DNE, UNE, mmE, mNm),
	TEHE(DSE, DNE, mmE, Dmm),

	TEHE(DNE, DNW, mNm, Dmm),
	TEHE(DSE, DSW, Dmm, mSm),
	TEHE(USE, USW, mSm, Umm),
	TEHE(UNE, UNW, Umm, mNm),

	TEHE(UNW, mmW, Umm, mNm),
	TEHE(USW, Umm, mmW, mSm),
	TEHE(USE, mSm, mmE, Umm),
	TEHE(UNE, mmE, Umm, mNm),

	TEHE(DNW, mmW, Dmm, mNm),
	TEHE(DSW, Dmm, mmW, mSm),
	TEHE(DSE, mSm, mmE, Dmm),
	TEHE(DNE, mmE, Dmm, mNm),

	TEHE(mmm, mmW, Umm, mNm),
	TEHE(mmm, Umm, mmW, mSm),
	TEHE(mmm, mSm, mmE, Umm),
	TEHE(mmm, mmE, Umm, mNm),

	TEHE(mmm, mmW, Dmm, mNm),
	TEHE(mmm, Dmm, mmW, mSm),
	TEHE(mmm, mSm, mmE, Dmm),
	TEHE(mmm, mmE, Dmm, mNm),
};
#undef			TEHE
struct WE_Offset//work-edge offset, used by CPU-side part 2
{
	short dx, dy, dz;
	char ke2[8];//{normal jump, (xyz, yz, xz, z, xy, y, x)@end}
};
__constant struct WE_Offset we_offsets[54-33]=//work-edge offsets, for redundant bits, access: [ke-33]={Xoffset, Yomask, Zomask, ke2[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)]}
{//	(dx, dy, dz)	(xyz,		yz,			xz,			z,			xy,			y,			x)@end		normal jump,			
	{1, 0, 0,		{DSE_DNE,	DSW_DNW,	DSE_DNE,	DSW_DNW,	DSE_DNE,	DSW_DNW,	DSE_DNE,	DSW_DNW}},//DSE_DNE
	{0, 1, 0,		{DNE_DNW,	DNE_DNW,	DSW_DSE,	DSW_DSE,	DNE_DNW,	DNE_DNW,	DSW_DSE,	DSW_DSE}},//DNE_DNW
	{1, 0, 0,		{DSE_USE,	DSW_USW,	DSE_USE,	DSW_USW,	DSE_USE,	DSW_USW,	DSE_USE,	DSW_USW}},//DSE_USE
	{1, 1, 0,		{DNE_UNE,	DNW_UNW,	DSE_USE,	DSW_USW,	DNE_UNE,	DNW_UNW,	DSE_USE,	DSW_USW}},//DNE_UNE//2D jump
	{0, 1, 0,		{DNW_UNW,	DNW_UNW,	DSW_USW,	DSW_USW,	DNW_UNW,	DNW_UNW,	DSW_USW,	DSW_USW}},//DNW_UNW
	{0, 0, 1,		{USW_USE,	USW_USE,	USW_USE,	USW_USE,	DSW_DSE,	DSW_DSE,	DSW_DSE,	DSW_DSE}},//USW_USE
	{1, 0, 1,		{USE_UNE,	USW_UNW,	USE_UNE,	USW_UNW,	DSE_DNE,	DSW_DNW,	DSE_DNE,	DSW_DNW}},//USE_UNE//2D jump
	{0, 1, 1,		{UNE_UNW,	UNE_UNW,	USW_USE,	USW_USE,	DNE_DNW,	DNE_DNW,	DSW_DSE,	DSW_DSE}},//UNE_UNW//2D jump
	{0, 0, 1,		{UNW_USW,	UNW_USW,	UNW_USW,	UNW_USW,	DSW_DNW,	DSW_DNW,	DSW_DNW,	DSW_DNW}},//UNW_USW

	{1, 0, 0,		{mmE_DSE,	mmW_DSW,	mmE_DSE,	mmW_DSW,	mmE_DSE,	mmW_DSW,	mmE_DSE,	mmW_DSW}},//mmE_DSE
	{1, 0, 0,		{mmE_USE,	mmW_USW,	mmE_USE,	mmW_USW,	mmE_USE,	mmW_USW,	mmE_USE,	mmW_USW}},//mmE_USE
	{1, 0, 0,		{mmE_UNE,	mmW_UNW,	mmE_UNE,	mmW_UNW,	mmE_UNE,	mmW_UNW,	mmE_UNE,	mmW_UNW}},//mmE_UNE
	{1, 0, 0,		{mmE_DNE,	mmW_DNW,	mmE_DNE,	mmW_DNW,	mmE_DNE,	mmW_DNW,	mmE_DNE,	mmW_DNW}},//mmE_DNE

	{0, 1, 0,		{mNm_DNW,	mNm_DNW,	mSm_DSW,	mSm_DSW,	mNm_DNW,	mNm_DNW,	mSm_DSW,	mSm_DSW}},//mNm_DNW
	{0, 1, 0,		{mNm_UNW,	mNm_UNW,	mSm_USW,	mSm_USW,	mNm_UNW,	mNm_UNW,	mSm_USW,	mSm_USW}},//mNm_UNW
	{0, 1, 0,		{mNm_UNE,	mNm_UNE,	mSm_USE,	mSm_USE,	mNm_UNE,	mNm_UNE,	mSm_USE,	mSm_USE}},//mNm_UNE
	{0, 1, 0,		{mNm_DNE,	mNm_DNE,	mSm_DSE,	mSm_DSE,	mNm_DNE,	mNm_DNE,	mSm_DSE,	mSm_DSE}},//mNm_DNE

	{0, 0, 1,		{Umm_USW,	Umm_USW,	Umm_USW,	Umm_USW,	Dmm_DSW,	Dmm_DSW,	Dmm_DSW,	Dmm_DSW}},//Umm_USW
	{0, 0, 1,		{Umm_USE,	Umm_USE,	Umm_USE,	Umm_USE,	Dmm_DSE,	Dmm_DSE,	Dmm_DSE,	Dmm_DSE}},//Umm_USE
	{0, 0, 1,		{Umm_UNE,	Umm_UNE,	Umm_UNE,	Umm_UNE,	Dmm_DNE,	Dmm_DNE,	Dmm_DNE,	Dmm_DNE}},//Umm_UNE
	{0, 0, 1,		{Umm_UNW,	Umm_UNW,	Umm_UNW,	Umm_UNW,	Dmm_DNW,	Dmm_DNW,	Dmm_DNW,	Dmm_DNW}},//Umm_UNW
};
char getbit(ulong x, char bit){return x>>bit&1;}
__kernel void ti3d_trgl_indices(__constant int *size, __constant ulong *ac_workidx, __constant ulong *edgeinfo, __constant ulong *workidx, __global int *indices)
{//size: {Xplaces, Yplaces, Zplaces, ndrSize, nvert_total, active_ndrsize, ntrgl_total}, workidx: {short ke, kx, ky, kz; (kz<<48|ky<<32|kx<<16|ke)}
	int id=get_global_id(0), active_ndrsize=size[5];
	if(id<active_ndrsize)
	{
		const int Xplaces=size[0], Yplaces=size[1], Zplaces=size[2], nvert_total=size[4], nindices=size[6]*3, XYplaces=Xplaces*Yplaces;

	//	indices[id]=42;//DEBUG
	//	indices[id]=id%nvert_total;//DEBUG

		ulong twi=ac_workidx[id];//active cube work index
		unsigned kw=(uint)twi, ki=(uint)(twi>>32)*3;
		int kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
		ulong work=edgeinfo[kw];
		ulong mask=obsmask[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)];

	//	indices[id]=ki<<9|kz<<6|ky<<3|kx;//DEBUG
	//	indices[id]=kz<<6|ky<<3|kx;//DEBUG

		for(int kth=0;kth<28;++kth)//for each tetrahedron kth of the 28 tetrahedra in the data cube
		{
			int kt6=kth*6;
			int e[]=//edge states
			{
				getbit(work, ti[kt6  ]),
				getbit(work, ti[kt6+1]),
				getbit(work, ti[kt6+2]),
				getbit(work, ti[kt6+3]),
				getbit(work, ti[kt6+4]),
				getbit(work, ti[kt6+5]),
			},
				esum=0;
			int tcvi[6]={};//tetrahedron cut vertex indices
			for(char kb=0;kb<6;++kb)//for each tetrahedron edge
			{
				if(e[kb])
				{
					int ke=ti[kth*6+kb];//relative edge index
					ulong idx=0;
					int original=getbit(mask, ke), not_original_mask=-!original;
					int kv=0;

					__constant struct WE_Offset *weo=we_offsets+((ke-33)&not_original_mask);
				//	__constant struct WE_Offset *weo=we_offsets+select(ke-33, 0, -original);//select(a, b, c)==c?b:a
					int
						xjm=kx<Xplaces-2, yjm=ky<Yplaces-2, zjm=kz<Zplaces-2,//jump masks: control jump in each dimension by boundary conditions: *places-2 is last cube (pre-last datapoint)

						kxF=kx+(weo->dx&-xjm&not_original_mask),//
						kyF=ky+(weo->dy&-yjm&not_original_mask),
						kzF=kz+(weo->dz&-zjm&not_original_mask),

					//	kxF=select(kx+(weo->dx&-xjm), kx, original),//TOSH: succeeds, MSI: succeeds
					//	kyF=select(ky+(weo->dy&-yjm), ky, original),
					//	kzF=select(kz+(weo->dz&-zjm), kz, original),

						keF=select((int)weo->ke2[zjm<<2|yjm<<1|xjm], ke, -original);
					idx=(ulong)kzF<<48|(ulong)kyF<<32|(ulong)kxF<<16|keF;

					ulong v=0;//bin search
					int lk=0, rk=nvert_total-1;
					//if(kx==4&&ky==0&&kz==0&&ke==8)//
					//	printf("lk,\trk,\tkv,\tv,\tworkidx[kv], (idx=0x%016LLX)\n", idx);//
					for(int ik=0;ik<31;++ik)
					{
						kv=(lk+rk)>>1;
						v=select(workidx[0], workidx[kv], (ulong)(kv>=0&&kv<nvert_total));
						//if(kx==4&&ky==0&&kz==0&&ke==8)//
						//	printf("%d, %d, %d, 0x%016LLX, 0x%016LLX\n", lk, rk, kv, v, workidx[kv]);//
						//	printf("lk=%d, rk=%d, kv=%d, v=0x%016LLX, workidx[kv]=0x%016LLX\n", lk, rk, kv, v, workidx[kv]);//
					//	v=workidx[kv];
						lk=select(lk, kv+1, v<idx);
						rk=select(rk, kv-1, v>idx);
					}
				//	tcvi[esum]=((kzF&7)<<12|(kyF&7)<<9|(kxF&7)<<6|(0x3F&keF))<<16|(kz&7)<<12|(ky&7)<<9|(kx&7)<<6|(0x3F&ke);//DEBUG
				//	if(kx==4&&ky==0&&kz==0&&ke==8)
				//	{
				//		printf("nvert_total=%d, active_ndrsize=%d, ntrgl_total=%d, nindices=%d\n", nvert_total, active_ndrsize, nindices/3, nindices);
				//		printf("(%d, %d, %d, %d), kth=%d, kb=%d\n", kx, ky, kz, ke, kth, kb);
				//		printf("original = %d != 0x%08X\n", original, not_original_mask);
				//		printf("weo-we_offsets=%d, select: %d, and: %d, MSB select=%d\n", (int)(weo-we_offsets), select(ke-33, 0, original), (ke-33)&not_original_mask, select(ke-33, 0, -original));
				//		printf("-xjm=0x%08X, -yjm=0x%08X, -zjm=0x%08X\n", -xjm, -yjm, -zjm);
				//		printf("F(%d, %d, %d, %d)\n", kxF, kyF, kzF, keF);
				//		printf("lk=%d, rk=%d, kv=%d\n\n", lk, rk, kv);
				//	}
				//	tcvi[esum]=(kzF&7)<<12|(kyF&7)<<9|(kxF&7)<<6|(0x3F&keF);//DEBUG
				//	tcvi[esum]=select(-42, kv, lk==rk);
				//	if(tcvi[esum]==-1)//RESULTS 0 or -43
				//		tcvi[esum]=-43;
				//	tcvi[esum]=select(-42, tcvi[esum], tcvi[esum]!=-1);//DEBUG
					tcvi[esum]=kv;
					++esum;
				}
			}
			if(esum)//tetrahedron has active edges (produced triangles)
			{
				int //one_triangle=esum==3&&ki+2<nindices,
					two_triangles=esum==4&&ki+5<nindices;
			//	indices[ki]=tcvi[0];//

			//	indices[select(ki  , ki+3, two_triangles)]=-id;//DEBUG
			//	indices[select(ki+1, ki+4, two_triangles)]=-id;
			//	indices[select(ki+2, ki+5, two_triangles)]=-id;
			//	indices[ki  ]=id;
			//	indices[ki+1]=id;
			//	indices[ki+2]=id;

			//	indices[select(ki  , ki+3, two_triangles)]=kx;//DEBUG
			//	indices[select(ki+1, ki+4, two_triangles)]=ky;
			//	indices[select(ki+2, ki+5, two_triangles)]=kz;
			//	indices[ki  ]=-ki;
			//	indices[ki+1]=-kw;
			//	indices[ki+2]=-kth;

			//	indices[select(ki  , ki+3, two_triangles)]=43;//DEBUG
			//	indices[select(ki+1, ki+4, two_triangles)]=id;
			//	indices[select(ki+2, ki+5, two_triangles)]=ki;
			//	indices[ki  ]=size[4+id%3];
			//	indices[ki+1]=size[4+(id+1)%3];
			//	indices[ki+2]=size[4+(id+2)%3];

			//	indices[select(ki  , ki+3, two_triangles)]=43;//DEBUG
			//	indices[select(ki+1, ki+4, two_triangles)]=43;
			//	indices[select(ki+2, ki+5, two_triangles)]=43;
			//	indices[ki  ]=42;
			//	indices[ki+1]=42;
			//	indices[ki+2]=42;

			//	indices[ki  ]=nvert_total-2;//DEBUG
			//	indices[ki+1]=nvert_total-2;
			//	indices[ki+2]=nvert_total-2;

				indices[select(ki  , ki+3, two_triangles)]=tcvi[1];//2nd triangle (optional)
				indices[select(ki+1, ki+4, two_triangles)]=tcvi[2];
				indices[select(ki+2, ki+5, two_triangles)]=tcvi[3];
				indices[ki  ]=tcvi[0];//1st triangle
				indices[ki+1]=tcvi[1];
				indices[ki+2]=tcvi[2];
				ki+=3+3*two_triangles;
			}
		}//end tetrahedron loop		//*/
	}
}