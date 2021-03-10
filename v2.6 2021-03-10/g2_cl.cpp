//best viewed with tab size of 4 spaces
//g2_cl.cpp - Implementation of G2 OpenCL component.
//Copyright (C) 2012-2020  Ayman Wagih Mohsen, unless source link provided.
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.
#include		"g2_cl.h"
#include		"g2_resources.h"
#include		"g2_error.h"
#include		"g2_file.h"
#include		"g2_expr.h"
#include		"g2_graphics.h"
#define			CL_TARGET_OPENCL_VERSION 120
#include		<CL/opencl.h>
int				OCL_state=CL_NOTHING, OCL_version=0;
bool			cl_gl_interop=false;

//	#define		COPY_INITIALIZATION
	#define		BLOCKING_INITIALIZATION
//	#define		DEBUG3
	#define		V6_GPU
//	#define		V6_CPU
//	#define		V5_CPU
//	#define		DEBUG2//performance impact

	const bool	loadbinary=true;
//	const bool	loadbinary=false;//DEBUG

int				*rgb=nullptr;
std::vector<DebugInfo> debug_info;
std::vector<float> debug_vertices;
std::vector<int> debug_indices;
unsigned		g_nvert=0, g_ntrgl=0, g_max_trgl_cube=0;
#if 1
#define			UNW		0, 2, 2
#define			UNE		2, 2, 2
#define			USW		0, 0, 2
#define			USE		2, 0, 2
#define			DNW		0, 2, 0
#define			DNE		2, 2, 0
#define			DSW		0, 0, 0
#define			DSE		2, 0, 0
#define			Dmm		1, 1, 0
#define			mSm		1, 0, 1
#define			mmW		0, 1, 1
#define			mmm		1, 1, 1
#define			mmE		2, 1, 1
#define			mNm		1, 2, 1
#define			Umm		1, 1, 2
const char		debug_edges[54*7]=//x1, y1, z1, x2, y2, z2, color_idx		0: no addition, 1: half sample, 2: full sample
{
	//main (new) edge indices
	DSW, DSE, 0,	DSW, DNW, 0,	DSW, USW, 0,

	DSW, mmW, 1,	USW, mmW, 1,	UNW, mmW, 1,	DNW, mmW, 1,
	DSW, mSm, 1,	USW, mSm, 1,	USE, mSm, 1,	DSE, mSm, 1,
	DSW, Dmm, 1,	DSE, Dmm, 1,	DNE, Dmm, 1,	DNW, Dmm, 1,

	mmW, Dmm, 2,	mmW, mSm, 2,	mmW, Umm, 2,	mmW, mNm, 2,
	Umm, mNm, 2,	mNm, Dmm, 2,	Dmm, mSm, 2,	mSm, Umm, 2,
	mmE, Dmm, 2,	mmE, mSm, 2,	mmE, Umm, 2,	mmE, mNm, 2,

	mmm, Dmm, 3,	mmm, mSm, 3,	mmm, mmW, 3,	mmm, Umm, 3,	mmm, mNm, 3,	mmm, mmE, 3,

	//extended (redundant) edge indices
	DSE, DNE, 0,	DNE, DNW, 0,	DSE, USE, 0,	DNE, UNE, 0,	DNW, UNW, 0,	USW, USE, 0,	USE, UNE, 0,	UNE, UNW, 0,	UNW, USW, 0,
	mmE, DSE, 1,	mmE, USE, 1,	mmE, UNE, 1,	mmE, DNE, 1,
	mNm, DNW, 1,	mNm, UNW, 1,	mNm, UNE, 1,	mNm, DNE, 1,
	Umm, USW, 1,	Umm, USE, 1,	Umm, UNE, 1,	Umm, UNW, 1,
};
#undef			UNW
#undef			UNE
#undef			USW
#undef			USE
#undef			DNW
#undef			DNE
#undef			DSW
#undef			DSE
#undef			Dmm
#undef			mSm
#undef			mmW
#undef			mmm
#undef			mmE
#undef			mNm
#undef			Umm
#endif

//OpenCL API
#define 		DECL_CL_FUNC(clFunc)	decltype(clFunc) *p_##clFunc=nullptr
const int 		cl_api_decl_start=__LINE__;
DECL_CL_FUNC(clGetPlatformIDs);
DECL_CL_FUNC(clGetPlatformInfo);
DECL_CL_FUNC(clGetDeviceIDs);
DECL_CL_FUNC(clGetDeviceInfo);
DECL_CL_FUNC(clCreateContext);
DECL_CL_FUNC(clReleaseContext);
DECL_CL_FUNC(clRetainContext);
DECL_CL_FUNC(clGetContextInfo);
DECL_CL_FUNC(clCreateCommandQueue);
DECL_CL_FUNC(clCreateProgramWithSource);
DECL_CL_FUNC(clBuildProgram);
DECL_CL_FUNC(clGetProgramBuildInfo);
DECL_CL_FUNC(clGetProgramInfo);
DECL_CL_FUNC(clCreateProgramWithBinary);
DECL_CL_FUNC(clCreateBuffer);
DECL_CL_FUNC(clCreateKernel);
DECL_CL_FUNC(clSetKernelArg);
//DECL_CL_FUNC(clEnqueueFillBuffer);//OpenCL 1.2+
DECL_CL_FUNC(clEnqueueWriteBuffer);
DECL_CL_FUNC(clEnqueueNDRangeKernel);
DECL_CL_FUNC(clEnqueueReadBuffer);
DECL_CL_FUNC(clFlush);
DECL_CL_FUNC(clFinish);
DECL_CL_FUNC(clCreateFromGLBuffer);
DECL_CL_FUNC(clCreateFromGLTexture);//OpenCL 1.2+?
DECL_CL_FUNC(clReleaseMemObject);
const int 		cl_api_decl_end=__LINE__;
#undef			DECL_CL_FUNC
HMODULE			hOpenCL=nullptr;
void 			load_OpenCL_API()
{
	if(OCL_state<CL_API_LOADED)
	{
		static const char *opencl_so_paths[]=//https://stackoverflow.com/questions/31611790/using-opencl-in-the-new-android-studio
		{
			"OpenCL.dll",
		//	//Android
		//	"OpenCL.so",
		//	"/system/vendor/lib64/libOpenCL.so",
		//	"/system/lib/libOpenCL.so",
		//	"/system/vendor/lib/libOpenCL.so",
		//	"/system/vendor/lib/egl/libGLES_mali.so",
		//	"/system/vendor/lib/libPVROCL.so",
		//	"/data/data/org.pocl.libs/files/lib/libpocl.so",
		//	//Linux
		//	"/usr/lib/libOpenCL.so",
		//	"/usr/local/lib/libOpenCL.so",
		//	"/usr/local/lib/libpocl.so",
		//	"/usr/lib64/libOpenCL.so",
		//	"/usr/lib32/libOpenCL.so",
		};
		const int npaths=sizeof(opencl_so_paths)/sizeof(char*);
		for(int k=0;k<npaths&&!hOpenCL;++k)
			hOpenCL=LoadLibraryA(opencl_so_paths[k]);
		if(!hOpenCL)
		{
			OCL_state=CL_NOTHING;
			LOGERROR("Cannot find OpenCL library");
		}
		else
		{
			OCL_state=CL_LOADING_API;
#define		GET_CL_FUNC(handle, clFunc)				p_##clFunc=(decltype(p_##clFunc))GetProcAddress(handle, #clFunc), p_check((void(*)())p_##clFunc, __FILE__, __LINE__, #clFunc)
#define		GET_CL_FUNC_UNCHECKED(handle, clFunc)	p_##clFunc=(decltype(p_##clFunc))GetProcAddress(handle, #clFunc)
			const int cl_api_init_start=__LINE__;
			GET_CL_FUNC(hOpenCL, clGetPlatformIDs);
			GET_CL_FUNC(hOpenCL, clGetPlatformInfo);
			GET_CL_FUNC(hOpenCL, clGetDeviceIDs);
			GET_CL_FUNC(hOpenCL, clGetDeviceInfo);
			GET_CL_FUNC(hOpenCL, clCreateContext);
			GET_CL_FUNC(hOpenCL, clReleaseContext);
			GET_CL_FUNC(hOpenCL, clRetainContext);
			GET_CL_FUNC(hOpenCL, clGetContextInfo);
			GET_CL_FUNC(hOpenCL, clCreateCommandQueue);
			GET_CL_FUNC(hOpenCL, clCreateProgramWithSource);
			GET_CL_FUNC(hOpenCL, clBuildProgram);
			GET_CL_FUNC(hOpenCL, clGetProgramBuildInfo);
			GET_CL_FUNC(hOpenCL, clGetProgramInfo);
			GET_CL_FUNC(hOpenCL, clCreateProgramWithBinary);
			GET_CL_FUNC(hOpenCL, clCreateBuffer);
			GET_CL_FUNC(hOpenCL, clCreateKernel);
			GET_CL_FUNC(hOpenCL, clSetKernelArg);
		//	GET_CL_FUNC(hOpenCL, clEnqueueFillBuffer);
			GET_CL_FUNC(hOpenCL, clEnqueueWriteBuffer);
			GET_CL_FUNC(hOpenCL, clEnqueueNDRangeKernel);
			GET_CL_FUNC(hOpenCL, clEnqueueReadBuffer);
			GET_CL_FUNC(hOpenCL, clFlush);
			GET_CL_FUNC(hOpenCL, clFinish);
			GET_CL_FUNC(hOpenCL, clCreateFromGLBuffer);
			GET_CL_FUNC_UNCHECKED(hOpenCL, clCreateFromGLTexture);
			GET_CL_FUNC(hOpenCL, clReleaseMemObject);
			const int cl_api_init_end=__LINE__;
#undef		GET_CL_FUNC
#undef		GET_CL_FUNC_UNCHECKED
#if CL_API_DECL_END-(CL_API_DECL_START+1)!=CL_API_INIT_END-(CL_API_INIT_START+1)
#error "Number of declared functions not equal to number of initialized functions."
#endif
			const int n_functions=cl_api_decl_end-(cl_api_decl_start+1), n_initialized=cl_api_init_end-(cl_api_init_start+1);
			MY_ASSERT(n_functions==n_initialized, "Number of declared functions not equal to number of initialized functions.");
			if(OCL_state!=CL_NOTHING)
				OCL_state=CL_API_LOADED;
		}
	}
}
void 			unload_OpenCL_API()
{
	if(hOpenCL)//when finishing opencl session
	{
		FreeLibrary(hOpenCL), hOpenCL=nullptr;
		OCL_state=CL_NOTHING;
	}
}

enum 			CLDiscType
{
	DISC_C,//continuous
	DISC_I,//depends on input
	DISC_O,//depends on output
};
struct 			CLKernel
{
	int idx;//enum CLKernelIdx
	short
		signature,//enum CLKernelSignature
		disc_type;//enum CLDiscType
	const char
		*name,//kernel name string, if 0 then software
		*disc_name;//if 0 then always continuous
};
struct			KernelDB
{
	CLKernel *kernels;
	int nkernels;
};
#if 0
//#define		STRING_LITERAL(...)		#__VA_ARGS__//one line: all whitespace & newlines replaced by single space
namespace		CLSource
{
	static const char program_common[]=R"CLSRC(
)CLSRC";
	static const char program00[]=R"CLSRC(
)CLSRC";
	static const char program01[]=R"CLSRC(
)CLSRC";
	static const char program02[]=R"CLSRC(
)CLSRC";
	static const char program03[]=R"CLSRC(
)CLSRC";
	static const char program04[]=R"CLSRC(
)CLSRC";
	static const char program05[]=R"CLSRC(
)CLSRC";
	static const char program06[]=R"CLSRC(
)CLSRC";
	static const char program07[]=R"CLSRC(
)CLSRC";
	static const char program08[]=R"CLSRC(
)CLSRC";
	static const char program09[]=R"CLSRC(
)CLSRC";
	static const char program10[]=R"CLSRC(
)CLSRC";
	static const char program11[]=R"CLSRC(
)CLSRC";
	static const char program12[]=R"CLSRC(
)CLSRC";
	static const char program13[]=R"CLSRC(
//back to G2 kernels
)CLSRC";
	static const char program14[]=R"CLSRC(
//back to G2 kernels
)CLSRC";
	static const char program15[]=R"CLSRC(
//back to G2 kernels
)CLSRC";
	static const char program16[]=R"CLSRC(
)CLSRC";
	static const char program17[]=R"CLSRC(
)CLSRC";
	static const char program18[]=R"CLSRC(
)CLSRC";
	static const char program19[]=R"CLSRC(
)CLSRC";
	static const char program20[]=R"CLSRC(
)CLSRC";
	static const char program21[]=R"CLSRC(
)CLSRC";
	static const char program22[]=R"CLSRC(
)CLSRC";

	const char *programs[]=
	{
		CLSource::program00,
		CLSource::program01,
		CLSource::program02,
		CLSource::program03,
		CLSource::program04,
		CLSource::program05,
		CLSource::program06,
		CLSource::program07,
		CLSource::program08,
		CLSource::program09,
		CLSource::program10,
		CLSource::program11,
		CLSource::program12,
		CLSource::program13,
		CLSource::program14,
		CLSource::program15,
		CLSource::program16,
		CLSource::program17,
		CLSource::program18,
		CLSource::program19,

		CLSource::program20,//TI3D
		CLSource::program21,
		CLSource::program22,
	};
}//end CLSource
//#undef		STRING_LITERAL
const int		nprograms=sizeof(CLSource::programs)/sizeof(const char*), n_g2programs=20;
#endif
const int		nprograms=G2_CL_NSOURCES, n_g2programs=G2_CL_NG2SRC;
struct 			ProgramBinary
{
	unsigned char *bin;
	size_t size;
	void free()
	{
		if(bin)
			::free(bin);
		bin=nullptr, size=0;
	}
};
cl_platform_id	platform=nullptr;//0x00000000de763ed3
cl_device_id	device=nullptr;//changes when app is launched
cl_context		context=nullptr;//changes when resuming
cl_command_queue commandqueue=nullptr;//changes when resuming
cl_kernel		kernels[N_KERNELS]={nullptr};//all kernels
size_t 			g_maxlocalsize=0, g_maxlocalX=0, g_maxlocalY=0, g_maxlocalZ=0;
namespace 		G2_CL
{
	cl_program programs[nprograms]={nullptr};

//declare continuous kernel
#define	DECL_C(SIG, ret, arg, NAME, name)	{SIG##_##NAME, SIG_##SIG, DISC_C, #ret "_" #arg "_" #name, nullptr}
//declare kernel with discontinuities depending on input
#define	DECL_I(SIG, ret, arg, NAME, name)	{SIG##_##NAME, SIG_##SIG, DISC_I, #ret "_" #arg "_" #name, "disc_" #arg "_" #name"_i"}
//declare kernel with discontinuities depending on output
#define	DECL_O(SIG, ret, arg, NAME, name)	{SIG##_##NAME, SIG_##SIG, DISC_O, #ret "_" #arg "_" #name, "disc_" #ret "_" #name"_o"}

//kernel signatures, DISCTYPE: C continuous, I depends on input, or O depends on output
#define	DECL_R_R(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_R, r, r, NAME, name)
#define	DECL_C_C(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_C, c, c, NAME, name)
#define	DECL_Q_Q(NAME, name, DISCTYPE)		DECL_##DISCTYPE(Q_Q, q, q, NAME, name)
#define	DECL_R_RR(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_RR, r, rr, NAME, name)
#define	DECL_C_RC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_RC, c, rc, NAME, name)
#define	DECL_Q_RQ(NAME, name, DISCTYPE)		DECL_##DISCTYPE(Q_RQ, q, rq, NAME, name)

#define	DECL_C_CR(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_CR, c, cr, NAME, name)
#define	DECL_C_CC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_CC, c, cc, NAME, name)
#define	DECL_Q_CQ(NAME, name, DISCTYPE)		DECL_##DISCTYPE(Q_CQ, q, cq, NAME, name)

#define	DECL_Q_QR(NAME, name, DISCTYPE)		DECL_##DISCTYPE(Q_QR, q, qr, NAME, name)
#define	DECL_Q_QC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(Q_QC, q, qc, NAME, name)
#define	DECL_Q_QQ(NAME, name, DISCTYPE)		DECL_##DISCTYPE(Q_QQ, q, qq, NAME, name)

#define	DECL_C_R(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_R, c, r, NAME, name)
#define	DECL_C_Q(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_Q, c, q, NAME, name)
#define	DECL_R_C(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_C, r, c, NAME, name)
#define	DECL_R_Q(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_Q, r, q, NAME, name)
#define	DECL_C_RR(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_RR, c, rr, NAME, name)

#define	DECL_R_RC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_RC, r, rc, NAME, name)
#define	DECL_R_RQ(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_RQ, r, rq, NAME, name)
#define	DECL_R_CR(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_CR, r, cr, NAME, name)
#define	DECL_R_CC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_CC, r, cc, NAME, name)
#define	DECL_R_CQ(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_CQ, r, cq, NAME, name)
#define	DECL_R_QR(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_QR, r, qr, NAME, name)
#define	DECL_R_QC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_QC, r, qc, NAME, name)
#define	DECL_R_QQ(NAME, name, DISCTYPE)		DECL_##DISCTYPE(R_QQ, r, qq, NAME, name)

#define	DECL_C_QC(NAME, name, DISCTYPE)		DECL_##DISCTYPE(C_QC, c, qc, NAME, name)

//declare function implemented in software
#define	DECL_SW(SIG, NAME, DISCTYPE)		{SIG##_##NAME, SIG_##SIG, DISC_##DISCTYPE, nullptr, nullptr}

//declare kernel - outdated
//#define	DECL_K(signature, name, disctype, namestr, discstr)		{signature##_##name, CL_##signature, DISC_##disctype, namestr, discstr}

	CLKernel kernels_p00[]=
	{
		//G2 functions
		DECL_R_R(SETZERO, setzero, C),
		DECL_C_C(SETZERO, setzero, C),
		DECL_Q_Q(SETZERO, setzero, C),

		DECL_R_R(CEIL, ceil, O),
		DECL_C_C(CEIL, ceil, O),
		DECL_Q_Q(CEIL, ceil, O),

		DECL_R_R(FLOOR, floor, O),
		DECL_C_C(FLOOR, floor, O),
		DECL_Q_Q(FLOOR, floor, O),

		DECL_R_R(ROUND, round, O),
		DECL_C_C(ROUND, round, O),
		DECL_Q_Q(ROUND, round, O),

		DECL_R_R(INT, int, O),
		DECL_C_C(INT, int, O),
		DECL_Q_Q(INT, int, O),

		DECL_R_R(FRAC, frac, I),
		DECL_C_C(FRAC, frac, I),
		DECL_Q_Q(FRAC, frac, I),
	};
	CLKernel kernels_p01[]=
	{
		DECL_R_R(ABS, abs, C),
		DECL_R_C(ABS, abs, C),
		DECL_R_Q(ABS, abs, C),

		DECL_R_R(ARG, arg, I),
		DECL_R_C(ARG, arg, I),
		DECL_R_Q(ARG, arg, I),

		DECL_R_C(REAL, real, C),

		DECL_R_C(IMAG, imag, C),

		DECL_C_C(CONJUGATE, conjugate, C),
		DECL_Q_Q(CONJUGATE, conjugate, C),

		DECL_C_R(POLAR, polar, I),
		DECL_C_C(POLAR, polar, I),
		DECL_C_Q(POLAR, polar, I),

		DECL_C_C(CARTESIAN, cartesian, C),
		DECL_Q_Q(CARTESIAN, cartesian, C),
	};
	CLKernel kernels_p02[]=
	{
		DECL_R_RR(PLUS, plus, C),
		DECL_C_RC(PLUS, plus, C),
		DECL_Q_RQ(PLUS, plus, C),
		DECL_C_CR(PLUS, plus, C),
		DECL_C_CC(PLUS, plus, C),
		DECL_Q_CQ(PLUS, plus, C),
		DECL_Q_QR(PLUS, plus, C),
		DECL_Q_QC(PLUS, plus, C),
		DECL_Q_QQ(PLUS, plus, C),

		DECL_R_R(MINUS, minus, C),
		DECL_C_C(MINUS, minus, C),
		DECL_Q_Q(MINUS, minus, C),
		DECL_R_RR(MINUS, minus, C),
		DECL_C_RC(MINUS, minus, C),
		DECL_Q_RQ(MINUS, minus, C),
		DECL_C_CR(MINUS, minus, C),
		DECL_C_CC(MINUS, minus, C),
		DECL_Q_CQ(MINUS, minus, C),
		DECL_Q_QR(MINUS, minus, C),
		DECL_Q_QC(MINUS, minus, C),
		DECL_Q_QQ(MINUS, minus, C),
		
		DECL_R_RR(MULTIPLY, multiply, C),
		DECL_C_RC(MULTIPLY, multiply, C),
		DECL_Q_RQ(MULTIPLY, multiply, C),
		DECL_C_CR(MULTIPLY, multiply, C),
		DECL_C_CC(MULTIPLY, multiply, C),
		DECL_Q_CQ(MULTIPLY, multiply, C),
		DECL_Q_QR(MULTIPLY, multiply, C),
		DECL_Q_QC(MULTIPLY, multiply, C),
		DECL_Q_QQ(MULTIPLY, multiply, C),
		
		DECL_R_R(DIVIDE, divide, I),
		DECL_C_C(DIVIDE, divide, I),
		DECL_Q_Q(DIVIDE, divide, I),
		DECL_R_RR(DIVIDE, divide, I),
		DECL_C_RC(DIVIDE, divide, I),
		DECL_Q_RQ(DIVIDE, divide, I),
		DECL_C_CR(DIVIDE, divide, I),
		DECL_C_CC(DIVIDE, divide, I),
		DECL_Q_CQ(DIVIDE, divide, I),
		DECL_Q_QR(DIVIDE, divide, I),
		DECL_Q_QC(DIVIDE, divide, I),
		DECL_Q_QQ(DIVIDE, divide, I),
	};
	CLKernel kernels_p03[]=
	{
		DECL_R_RR(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_RC(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_RQ(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_CR(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_CC(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_CQ(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_QR(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_QC(LOGIC_DIVIDES, logic_divides, O),
		DECL_R_QQ(LOGIC_DIVIDES, logic_divides, O),

		DECL_R_RR(POWER_REAL, power_real, I),
		DECL_C_CR(POWER_REAL, power_real, I),
		DECL_Q_QR(POWER_REAL, power_real, I),

		DECL_C_CR(POW, pow, I),
		DECL_C_CC(POW, pow, I),
		DECL_Q_CQ(POW, pow, I),
		DECL_Q_QR(POW, pow, I),
		DECL_Q_QC(POW, pow, I),
		DECL_Q_QQ(POW, pow, I),

		DECL_C_C(LN, ln, I),
		DECL_Q_Q(LN, ln, I),

		DECL_C_C(LOG, log, I),
		DECL_Q_Q(LOG, log, I),
		DECL_C_CR(LOG, log, I),
		DECL_C_CC(LOG, log, I),
		DECL_Q_CQ(LOG, log, I),
		DECL_Q_QC(LOG, log, I),
		DECL_Q_QQ(LOG, log, I),
		
		DECL_SW(C_RR, TETRATE, I),
		DECL_SW(C_RC, TETRATE, I),
		DECL_SW(C_CR, TETRATE, I),
		DECL_SW(C_CC, TETRATE, I),
		DECL_SW(Q_QR, TETRATE, I),
		//{C_RR_TETRATE, CL_, 0, nullptr, nullptr},
		//{C_RC_TETRATE, CL_, 0, nullptr, nullptr},
		//{C_CR_TETRATE, CL_, 0, nullptr, nullptr},
		//{C_CC_TETRATE, CL_, 0, nullptr, nullptr},
		//{Q_QR_TETRATE, CL_, 0, nullptr, nullptr},
		
		DECL_SW(C_RR, PENTATE, I),
		DECL_SW(C_CR, PENTATE, I),
	};
	CLKernel kernels_p04[]=
	{
		DECL_R_R(BITWISE_SHIFT_LEFT_L, bitwise_shift_left_l, O),
		DECL_C_C(BITWISE_SHIFT_LEFT_L, bitwise_shift_left_l, O),
		DECL_Q_Q(BITWISE_SHIFT_LEFT_L, bitwise_shift_left_l, O),
		DECL_R_R(BITWISE_SHIFT_LEFT_R, bitwise_shift_left_r, C),
		DECL_C_C(BITWISE_SHIFT_LEFT_R, bitwise_shift_left_r, C),
		DECL_Q_Q(BITWISE_SHIFT_LEFT_R, bitwise_shift_left_r, C),
		DECL_R_RR(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_C_RC(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_Q_RQ(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_C_CR(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_C_CC(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_Q_CQ(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_Q_QR(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_Q_QC(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),
		DECL_Q_QQ(BITWISE_SHIFT_LEFT, bitwise_shift_left, I),

		DECL_R_R(BITWISE_SHIFT_RIGHT_L, bitwise_shift_right_l, O),
		DECL_C_C(BITWISE_SHIFT_RIGHT_L, bitwise_shift_right_l, O),
		DECL_Q_Q(BITWISE_SHIFT_RIGHT_L, bitwise_shift_right_l, O),
		DECL_R_R(BITWISE_SHIFT_RIGHT_R, bitwise_shift_right_r, C),
		DECL_C_C(BITWISE_SHIFT_RIGHT_R, bitwise_shift_right_r, C),
		DECL_Q_Q(BITWISE_SHIFT_RIGHT_R, bitwise_shift_right_r, C),
		DECL_R_RR(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_C_RC(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_Q_RQ(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_C_CR(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_C_CC(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_Q_CQ(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_Q_QR(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_Q_QC(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
		DECL_Q_QQ(BITWISE_SHIFT_RIGHT, bitwise_shift_right, I),
	};
	CLKernel kernels_p05[]=
	{
		DECL_R_R(BITWISE_NOT, bitwise_not, I),
		DECL_C_C(BITWISE_NOT, bitwise_not, I),
		DECL_Q_Q(BITWISE_NOT, bitwise_not, I),

		DECL_R_R(BITWISE_AND, bitwise_and, O),
		DECL_C_C(BITWISE_AND, bitwise_and, O),
		DECL_Q_Q(BITWISE_AND, bitwise_and, O),
		DECL_R_RR(BITWISE_AND, bitwise_and, O),
		DECL_C_RC(BITWISE_AND, bitwise_and, O),
		DECL_Q_RQ(BITWISE_AND, bitwise_and, O),
		DECL_C_CR(BITWISE_AND, bitwise_and, O),
		DECL_C_CC(BITWISE_AND, bitwise_and, O),
		DECL_Q_CQ(BITWISE_AND, bitwise_and, O),
		DECL_Q_QR(BITWISE_AND, bitwise_and, O),
		DECL_Q_QC(BITWISE_AND, bitwise_and, O),
		DECL_Q_QQ(BITWISE_AND, bitwise_and, O),

		DECL_R_R(BITWISE_NAND, bitwise_nand, O),
		DECL_C_C(BITWISE_NAND, bitwise_nand, O),
		DECL_Q_Q(BITWISE_NAND, bitwise_nand, O),
		DECL_R_RR(BITWISE_NAND, bitwise_nand, O),
		DECL_C_RC(BITWISE_NAND, bitwise_nand, O),
		DECL_Q_RQ(BITWISE_NAND, bitwise_nand, O),
		DECL_C_CR(BITWISE_NAND, bitwise_nand, O),
		DECL_C_CC(BITWISE_NAND, bitwise_nand, O),
		DECL_Q_CQ(BITWISE_NAND, bitwise_nand, O),
		DECL_Q_QR(BITWISE_NAND, bitwise_nand, O),
		DECL_Q_QC(BITWISE_NAND, bitwise_nand, O),
		DECL_Q_QQ(BITWISE_NAND, bitwise_nand, O),

		DECL_R_R(BITWISE_OR, bitwise_or, O),
		DECL_C_C(BITWISE_OR, bitwise_or, O),
		DECL_Q_Q(BITWISE_OR, bitwise_or, O),
		DECL_R_RR(BITWISE_OR, bitwise_or, O),
		DECL_C_RC(BITWISE_OR, bitwise_or, O),
		DECL_Q_RQ(BITWISE_OR, bitwise_or, O),
		DECL_C_CR(BITWISE_OR, bitwise_or, O),
		DECL_C_CC(BITWISE_OR, bitwise_or, O),
		DECL_Q_CQ(BITWISE_OR, bitwise_or, O),
		DECL_Q_QR(BITWISE_OR, bitwise_or, O),
		DECL_Q_QC(BITWISE_OR, bitwise_or, O),
		DECL_Q_QQ(BITWISE_OR, bitwise_or, O),
	};
	CLKernel kernels_p06[]=
	{
		DECL_R_R(BITWISE_NOR, bitwise_nor, O),
		DECL_C_C(BITWISE_NOR, bitwise_nor, O),
		DECL_Q_Q(BITWISE_NOR, bitwise_nor, O),
		DECL_R_RR(BITWISE_NOR, bitwise_nor, O),
		DECL_C_RC(BITWISE_NOR, bitwise_nor, O),
		DECL_Q_RQ(BITWISE_NOR, bitwise_nor, O),
		DECL_C_CR(BITWISE_NOR, bitwise_nor, O),
		DECL_C_CC(BITWISE_NOR, bitwise_nor, O),
		DECL_Q_CQ(BITWISE_NOR, bitwise_nor, O),
		DECL_Q_QR(BITWISE_NOR, bitwise_nor, O),
		DECL_Q_QC(BITWISE_NOR, bitwise_nor, O),
		DECL_Q_QQ(BITWISE_NOR, bitwise_nor, O),

		DECL_R_R(BITWISE_XOR, bitwise_xor, O),
		DECL_C_C(BITWISE_XOR, bitwise_xor, O),
		DECL_Q_Q(BITWISE_XOR, bitwise_xor, O),
		DECL_R_RR(BITWISE_XOR, bitwise_xor, O),
		DECL_C_RC(BITWISE_XOR, bitwise_xor, O),
		DECL_Q_RQ(BITWISE_XOR, bitwise_xor, O),
		DECL_C_CR(BITWISE_XOR, bitwise_xor, O),
		DECL_C_CC(BITWISE_XOR, bitwise_xor, O),
		DECL_Q_CQ(BITWISE_XOR, bitwise_xor, O),
		DECL_Q_QR(BITWISE_XOR, bitwise_xor, O),
		DECL_Q_QC(BITWISE_XOR, bitwise_xor, O),
		DECL_Q_QQ(BITWISE_XOR, bitwise_xor, O),

		DECL_R_R(BITWISE_XNOR, bitwise_xnor, O),
		DECL_C_C(BITWISE_XNOR, bitwise_xnor, O),
		DECL_Q_Q(BITWISE_XNOR, bitwise_xnor, O),
		DECL_R_RR(BITWISE_XNOR, bitwise_xnor, O),
		DECL_C_RC(BITWISE_XNOR, bitwise_xnor, O),
		DECL_Q_RQ(BITWISE_XNOR, bitwise_xnor, O),
		DECL_C_CR(BITWISE_XNOR, bitwise_xnor, O),
		DECL_C_CC(BITWISE_XNOR, bitwise_xnor, O),
		DECL_Q_CQ(BITWISE_XNOR, bitwise_xnor, O),
		DECL_Q_QR(BITWISE_XNOR, bitwise_xnor, O),
		DECL_Q_QC(BITWISE_XNOR, bitwise_xnor, O),
		DECL_Q_QQ(BITWISE_XNOR, bitwise_xnor, O),
	};
	CLKernel kernels_p07[]=
	{
		DECL_R_R(LOGIC_EQUAL, logic_equal, O),
		DECL_R_C(LOGIC_EQUAL, logic_equal, O),
		DECL_R_Q(LOGIC_EQUAL, logic_equal, O),
		DECL_R_RR(LOGIC_EQUAL, logic_equal, O),
		DECL_R_RC(LOGIC_EQUAL, logic_equal, O),
		DECL_R_RQ(LOGIC_EQUAL, logic_equal, O),
		DECL_R_CR(LOGIC_EQUAL, logic_equal, O),
		DECL_R_CC(LOGIC_EQUAL, logic_equal, O),
		DECL_R_CQ(LOGIC_EQUAL, logic_equal, O),
		DECL_R_QR(LOGIC_EQUAL, logic_equal, O),
		DECL_R_QC(LOGIC_EQUAL, logic_equal, O),
		DECL_R_QQ(LOGIC_EQUAL, logic_equal, O),

		DECL_R_R(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_C(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_Q(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_RR(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_RC(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_RQ(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_CR(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_CC(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_CQ(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_QR(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_QC(LOGIC_NOT_EQUAL, logic_not_equal, O),
		DECL_R_QQ(LOGIC_NOT_EQUAL, logic_not_equal, O),

		DECL_R_R(LOGIC_LESS_L, logic_less_l, O),
		DECL_R_C(LOGIC_LESS_L, logic_less_l, O),
		DECL_R_Q(LOGIC_LESS_L, logic_less_l, O),
		DECL_R_R(LOGIC_LESS_R, logic_less_r, O),
		DECL_R_C(LOGIC_LESS_R, logic_less_r, O),
		DECL_R_Q(LOGIC_LESS_R, logic_less_r, O),
		DECL_R_RR(LOGIC_LESS, logic_less, O),
		DECL_R_RC(LOGIC_LESS, logic_less, O),
		DECL_R_RQ(LOGIC_LESS, logic_less, O),
		DECL_R_CR(LOGIC_LESS, logic_less, O),
		DECL_R_CC(LOGIC_LESS, logic_less, O),
		DECL_R_CQ(LOGIC_LESS, logic_less, O),
		DECL_R_QR(LOGIC_LESS, logic_less, O),
		DECL_R_QC(LOGIC_LESS, logic_less, O),
		DECL_R_QQ(LOGIC_LESS, logic_less, O),
	};
	CLKernel kernels_p08[]=
	{
		DECL_R_R(LOGIC_LESS_EQUAL_L, logic_less_equal_l, O),
		DECL_R_C(LOGIC_LESS_EQUAL_L, logic_less_equal_l, O),
		DECL_R_Q(LOGIC_LESS_EQUAL_L, logic_less_equal_l, O),
		DECL_R_R(LOGIC_LESS_EQUAL_R, logic_less_equal_r, O),
		DECL_R_C(LOGIC_LESS_EQUAL_R, logic_less_equal_r, O),
		DECL_R_Q(LOGIC_LESS_EQUAL_R, logic_less_equal_r, O),
		DECL_R_RR(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_RC(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_RQ(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_CR(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_CC(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_CQ(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_QR(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_QC(LOGIC_LESS_EQUAL, logic_less_equal, O),
		DECL_R_QQ(LOGIC_LESS_EQUAL, logic_less_equal, O),

		DECL_R_R(LOGIC_GREATER_L, logic_greater_l, O),
		DECL_R_C(LOGIC_GREATER_L, logic_greater_l, O),
		DECL_R_Q(LOGIC_GREATER_L, logic_greater_l, O),
		DECL_R_R(LOGIC_GREATER_R, logic_greater_r, O),
		DECL_R_C(LOGIC_GREATER_R, logic_greater_r, O),
		DECL_R_Q(LOGIC_GREATER_R, logic_greater_r, O),
		DECL_R_RR(LOGIC_GREATER, logic_greater, O),
		DECL_R_RC(LOGIC_GREATER, logic_greater, O),
		DECL_R_RQ(LOGIC_GREATER, logic_greater, O),
		DECL_R_CR(LOGIC_GREATER, logic_greater, O),
		DECL_R_CC(LOGIC_GREATER, logic_greater, O),
		DECL_R_CQ(LOGIC_GREATER, logic_greater, O),
		DECL_R_QR(LOGIC_GREATER, logic_greater, O),
		DECL_R_QC(LOGIC_GREATER, logic_greater, O),
		DECL_R_QQ(LOGIC_GREATER, logic_greater, O),

		DECL_R_R(LOGIC_GREATER_EQUAL_L, logic_greater_equal_l, O),
		DECL_R_C(LOGIC_GREATER_EQUAL_L, logic_greater_equal_l, O),
		DECL_R_Q(LOGIC_GREATER_EQUAL_L, logic_greater_equal_l, O),
		DECL_R_R(LOGIC_GREATER_EQUAL_R, logic_greater_equal_r, O),
		DECL_R_C(LOGIC_GREATER_EQUAL_R, logic_greater_equal_r, O),
		DECL_R_Q(LOGIC_GREATER_EQUAL_R, logic_greater_equal_r, O),
		DECL_R_RR(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_RC(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_RQ(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_CR(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_CC(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_CQ(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_QR(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_QC(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
		DECL_R_QQ(LOGIC_GREATER_EQUAL, logic_greater_equal, O),
	};
	CLKernel kernels_p09[]=
	{
		DECL_R_R(LOGIC_NOT, logic_not, O),
		DECL_R_C(LOGIC_NOT, logic_not, O),
		DECL_R_Q(LOGIC_NOT, logic_not, O),

		DECL_R_RR(LOGIC_AND, logic_and, O),
		DECL_R_RC(LOGIC_AND, logic_and, O),
		DECL_R_RQ(LOGIC_AND, logic_and, O),
		DECL_R_CR(LOGIC_AND, logic_and, O),
		DECL_R_CC(LOGIC_AND, logic_and, O),
		DECL_R_CQ(LOGIC_AND, logic_and, O),
		DECL_R_QR(LOGIC_AND, logic_and, O),
		DECL_R_QC(LOGIC_AND, logic_and, O),
		DECL_R_QQ(LOGIC_AND, logic_and, O),

		DECL_R_RR(LOGIC_OR, logic_or, O),
		DECL_R_RC(LOGIC_OR, logic_or, O),
		DECL_R_RQ(LOGIC_OR, logic_or, O),
		DECL_R_CR(LOGIC_OR, logic_or, O),
		DECL_R_CC(LOGIC_OR, logic_or, O),
		DECL_R_CQ(LOGIC_OR, logic_or, O),
		DECL_R_QR(LOGIC_OR, logic_or, O),
		DECL_R_QC(LOGIC_OR, logic_or, O),
		DECL_R_QQ(LOGIC_OR, logic_or, O),

		DECL_R_RR(LOGIC_XOR, logic_xor, O),
		DECL_R_RC(LOGIC_XOR, logic_xor, O),
		DECL_R_RQ(LOGIC_XOR, logic_xor, O),
		DECL_R_CR(LOGIC_XOR, logic_xor, O),
		DECL_R_CC(LOGIC_XOR, logic_xor, O),
		DECL_R_CQ(LOGIC_XOR, logic_xor, O),
		DECL_R_QR(LOGIC_XOR, logic_xor, O),
		DECL_R_QC(LOGIC_XOR, logic_xor, O),
		DECL_R_QQ(LOGIC_XOR, logic_xor, O),
	};
	CLKernel kernels_p10[]=
	{
		DECL_R_RR(CONDITION_ZERO, condition_zero, I),
		DECL_C_RC(CONDITION_ZERO, condition_zero, I),
		DECL_Q_RQ(CONDITION_ZERO, condition_zero, I),
		DECL_C_CR(CONDITION_ZERO, condition_zero, I),
		DECL_C_CC(CONDITION_ZERO, condition_zero, I),
		DECL_Q_CQ(CONDITION_ZERO, condition_zero, I),
		DECL_Q_QR(CONDITION_ZERO, condition_zero, I),
		DECL_Q_QC(CONDITION_ZERO, condition_zero, I),
		DECL_Q_QQ(CONDITION_ZERO, condition_zero, I),

		DECL_R_R(PERCENT, percent, C),
		DECL_C_C(PERCENT, percent, C),
		DECL_Q_Q(PERCENT, percent, C),

		DECL_R_RR(MODULO, modulo, I),
		DECL_C_RC(MODULO, modulo, I),
		DECL_Q_RQ(MODULO, modulo, I),
		DECL_C_CR(MODULO, modulo, I),
		DECL_C_CC(MODULO, modulo, I),
		DECL_Q_CQ(MODULO, modulo, I),
		DECL_Q_QR(MODULO, modulo, I),
		DECL_Q_QC(MODULO, modulo, I),
		DECL_Q_QQ(MODULO, modulo, I),
	};
	CLKernel kernels_p11[]=
	{
		DECL_R_R(SGN, sgn, I),
		DECL_C_C(SGN, sgn, I),
		DECL_Q_Q(SGN, sgn, I),

		DECL_R_R(SQ, sq, C),
		DECL_C_C(SQ, sq, C),
		DECL_Q_Q(SQ, sq, C),

		DECL_C_C(SQRT, sqrt, C),
		DECL_Q_Q(SQRT, sqrt, C),

		DECL_R_R(INVSQRT, invsqrt, C),

		DECL_R_R(CBRT, cbrt, C),
		DECL_C_C(CBRT, cbrt, C),
		DECL_Q_Q(CBRT, cbrt, C),

		DECL_R_R(GAUSS, gauss, C),
		DECL_C_C(GAUSS, gauss, C),
		DECL_Q_Q(GAUSS, gauss, C),

		DECL_R_R(ERF, erf, C),

		DECL_SW(R_R, ZETA, I),
	};
	CLKernel kernels_p12[]=
	{
		DECL_R_R(TGAMMA, tgamma, I),
		DECL_C_C(TGAMMA, tgamma, I),
		DECL_Q_Q(TGAMMA, tgamma, I),
		DECL_SW(R_RR, TGAMMA, I),

		DECL_R_R(LOGGAMMA, loggamma, I),

		DECL_R_R(FACTORIAL, factorial, I),
		DECL_C_C(FACTORIAL, factorial, I),
		DECL_Q_Q(FACTORIAL, factorial, I),

		DECL_R_R(PERMUTATION, permutation, I),
		DECL_C_C(PERMUTATION, permutation, I),
		DECL_Q_Q(PERMUTATION, permutation, I),
		DECL_R_RR(PERMUTATION, permutation, I),
		DECL_C_CR(PERMUTATION, permutation, I),
		DECL_C_CC(PERMUTATION, permutation, I),
		DECL_Q_QQ(PERMUTATION, permutation, I),

		DECL_R_R(COMBINATION, combination, I),
		DECL_C_C(COMBINATION, combination, I),
		DECL_Q_Q(COMBINATION, combination, I),
		DECL_R_RR(COMBINATION, combination, I),
		DECL_C_CR(COMBINATION, combination, I),
		DECL_C_CC(COMBINATION, combination, I),
		DECL_Q_QQ(COMBINATION, combination, I),
	};
	CLKernel kernels_p13[]=
	{
		DECL_R_R(COS, cos, C),
		DECL_C_C(COS, cos, C),
		DECL_Q_Q(COS, cos, C),

		DECL_C_C(ACOS, acos, I),
		DECL_Q_Q(ACOS, acos, I),

		DECL_R_R(COSH, cosh, C),
		DECL_C_C(COSH, cosh, C),
		DECL_Q_Q(COSH, cosh, C),

		DECL_C_C(ACOSH, acosh, C),
		DECL_Q_Q(ACOSH, acosh, C),

		DECL_R_R(COSC, cosc, I),
		DECL_C_C(COSC, cosc, I),
		DECL_Q_Q(COSC, cosc, I),

		DECL_R_R(SEC, sec, I),
		DECL_C_C(SEC, sec, I),
		DECL_Q_Q(SEC, sec, I),

		DECL_C_C(ASEC, asec, I),
		DECL_Q_Q(ASEC, asec, I),

		DECL_R_R(SECH, sech, C),
		DECL_C_C(SECH, sech, I),
		DECL_Q_Q(SECH, sech, I),

		DECL_C_C(ASECH, asech, I),
		DECL_Q_Q(ASECH, asech, I),
	};
	CLKernel kernels_p14[]=
	{
		DECL_R_R(SIN, sin, C),
		DECL_C_C(SIN, sin, C),
		DECL_Q_Q(SIN, sin, C),

		DECL_C_C(ASIN, asin, I),
		DECL_Q_Q(ASIN, asin, I),

		DECL_R_R(SINH, sinh, C),
		DECL_C_C(SINH, sinh, C),
		DECL_Q_Q(SINH, sinh, C),

		DECL_R_R(ASINH, asinh, C),
		DECL_C_C(ASINH, asinh, I),
		DECL_Q_Q(ASINH, asinh, I),

		DECL_R_R(SINC, sinc, C),
		DECL_C_C(SINC, sinc, C),
		DECL_Q_Q(SINC, sinc, C),

		DECL_R_R(SINHC, sinhc, C),
		DECL_C_C(SINHC, sinhc, C),
		DECL_Q_Q(SINHC, sinhc, C),

		DECL_R_R(CSC, csc, I),
		DECL_C_C(CSC, csc, I),
		DECL_Q_Q(CSC, csc, I),

		DECL_C_C(ACSC, acsc, I),
		DECL_Q_Q(ACSC, acsc, I),

		DECL_R_R(CSCH, csch, I),
		DECL_C_C(CSCH, csch, I),
		DECL_Q_Q(CSCH, csch, I),

		DECL_R_R(ACSCH, acsch, I),
		DECL_C_C(ACSCH, acsch, I),
		DECL_Q_Q(ACSCH, acsch, I),
	};
	CLKernel kernels_p15[]=
	{
		DECL_R_R(TAN, tan, I),
		DECL_C_C(TAN, tan, I),
		DECL_Q_Q(TAN, tan, I),

		DECL_R_R(ATAN, atan, C),
		DECL_C_C(ATAN, atan, I),
		DECL_Q_Q(ATAN, atan, I),
		DECL_R_RR(ATAN, atan, I),
		DECL_C_RC(ATAN, atan, I),
		DECL_Q_RQ(ATAN, atan, I),
		DECL_C_CR(ATAN, atan, I),
		DECL_C_CC(ATAN, atan, I),
		DECL_Q_CQ(ATAN, atan, I),
		DECL_Q_QR(ATAN, atan, I),
		DECL_Q_QC(ATAN, atan, I),
		DECL_Q_QQ(ATAN, atan, I),

		DECL_R_R(TANH, tanh, C),
		DECL_C_C(TANH, tanh, C),
		DECL_Q_Q(TANH, tanh, C),

		DECL_C_C(ATANH, atanh, I),
		DECL_Q_Q(ATANH, atanh, I),

		DECL_R_R(TANC, tanc, I),
		DECL_C_C(TANC, tanc, I),
		DECL_Q_Q(TANC, tanc, I),

		DECL_R_R(COT, cot, I),
		DECL_C_C(COT, cot, I),
		DECL_Q_Q(COT, cot, I),

		DECL_R_R(ACOT, acot, I),
		DECL_C_C(ACOT, acot, I),
		DECL_Q_Q(ACOT, acot, I),

		DECL_R_R(COTH, coth, I),
		DECL_C_C(COTH, coth, I),
		DECL_Q_Q(COTH, coth, I),

		DECL_C_C(ACOTH, acoth, I),
		DECL_Q_Q(ACOTH, acoth, I),
	};
	CLKernel kernels_p16[]=
	{
		DECL_R_R(EXP, exp, C),
		DECL_C_C(EXP, exp, C),
		DECL_Q_Q(EXP, exp, C),

		DECL_R_R(FIB, fib, C),
		DECL_C_C(FIB, fib, C),
		DECL_Q_Q(FIB, fib, C),

		DECL_SW(R_R, RANDOM, O),
		DECL_SW(C_C, RANDOM, O),
		DECL_SW(Q_Q, RANDOM, O),
		DECL_SW(R_RR, RANDOM, O),
		DECL_SW(C_CR, RANDOM, O),
		DECL_SW(C_CC, RANDOM, O),
		DECL_SW(Q_QQ, RANDOM, O),

		DECL_SW(R_R, BETA, I),
		DECL_SW(R_RR, BETA, I),

		DECL_SW(R_R, BESSEL_J, I),
		DECL_SW(R_RR, BESSEL_J, I),

		DECL_SW(R_R, BESSEL_Y, I),
		DECL_SW(R_RR, BESSEL_Y, I),

		DECL_SW(C_R, HANKEL1, I),
		DECL_SW(C_C, HANKEL1, I),
		DECL_SW(C_RR, HANKEL1, I),

		DECL_R_R(STEP, step, I),
		DECL_C_C(STEP, step, I),
		DECL_Q_Q(STEP, step, I),

		DECL_R_R(RECT, rect, I),
		DECL_C_C(RECT, rect, I),
		DECL_Q_Q(RECT, rect, I),

		DECL_R_R(TRGL, trgl, C),
		DECL_R_C(TRGL, trgl, C),
		DECL_R_Q(TRGL, trgl, C),
	};
	CLKernel kernels_p17[]=
	{
		DECL_R_R(SQWV, sqwv, O),
		DECL_R_C(SQWV, sqwv, O),
		DECL_R_Q(SQWV, sqwv, O),
		DECL_R_RR(SQWV, sqwv, O),
		DECL_R_RC(SQWV, sqwv, O),
		DECL_R_RQ(SQWV, sqwv, O),
		DECL_R_CR(SQWV, sqwv, O),
		DECL_R_CC(SQWV, sqwv, O),
		DECL_R_CQ(SQWV, sqwv, O),
		DECL_R_QR(SQWV, sqwv, O),
		DECL_R_QC(SQWV, sqwv, O),
		DECL_R_QQ(SQWV, sqwv, O),

		DECL_R_R(TRWV, trwv, C),
		DECL_R_C(TRWV, trwv, C),
		DECL_R_Q(TRWV, trwv, C),
		DECL_R_RR(TRWV, trwv, C),
		DECL_R_RC(TRWV, trwv, C),
		DECL_R_RQ(TRWV, trwv, C),
		DECL_R_CR(TRWV, trwv, C),
		DECL_R_CC(TRWV, trwv, C),
		DECL_R_CQ(TRWV, trwv, C),
		DECL_R_QR(TRWV, trwv, C),
		DECL_R_QC(TRWV, trwv, C),
		DECL_R_QQ(TRWV, trwv, C),

		DECL_R_R(SAW, saw, I),
		DECL_R_C(SAW, saw, I),
		DECL_R_Q(SAW, saw, I),
		DECL_R_RR(SAW, saw, I),
		DECL_R_RC(SAW, saw, I),
		DECL_R_RQ(SAW, saw, I),
		DECL_R_CR(SAW, saw, I),
		DECL_R_CC(SAW, saw, I),
		DECL_R_CQ(SAW, saw, I),
		DECL_R_QR(SAW, saw, I),
		DECL_R_QC(SAW, saw, I),
		DECL_R_QQ(SAW, saw, I),
	};
	CLKernel kernels_p18[]=
	{
		DECL_R_RR(HYPOT, hypot, C),

		DECL_R_R(MANDELBROT, mandelbrot, O),
		DECL_R_C(MANDELBROT, mandelbrot, O),
		DECL_R_RR(MANDELBROT, mandelbrot, O),
		DECL_R_CR(MANDELBROT, mandelbrot, O),

		DECL_R_RR(MIN, min, C),
		DECL_C_CR(MIN, min, C),
		DECL_C_CC(MIN, min, C),
		DECL_Q_QQ(MIN, min, C),

		DECL_R_RR(MAX, max, C),
		DECL_C_CR(MAX, max, C),
		DECL_C_CC(MAX, max, C),
		DECL_Q_QQ(MAX, max, C),
	};
	CLKernel kernels_p19[]=
	{
		DECL_R_RR(CONDITIONAL_110, conditional_110, I),
		DECL_C_RC(CONDITIONAL_110, conditional_110, I),
		DECL_Q_RQ(CONDITIONAL_110, conditional_110, I),
		DECL_R_CR(CONDITIONAL_110, conditional_110, I),
		DECL_C_CC(CONDITIONAL_110, conditional_110, I),
		DECL_Q_CQ(CONDITIONAL_110, conditional_110, I),
		DECL_R_QR(CONDITIONAL_110, conditional_110, I),
		DECL_C_QC(CONDITIONAL_110, conditional_110, I),
		DECL_Q_QQ(CONDITIONAL_110, conditional_110, I),

		DECL_R_RR(CONDITIONAL_101, conditional_101, I),
		DECL_C_RC(CONDITIONAL_101, conditional_101, I),
		DECL_Q_RQ(CONDITIONAL_101, conditional_101, I),
		DECL_R_CR(CONDITIONAL_101, conditional_101, I),
		DECL_C_CC(CONDITIONAL_101, conditional_101, I),
		DECL_Q_CQ(CONDITIONAL_101, conditional_101, I),
		DECL_R_QR(CONDITIONAL_101, conditional_101, I),
		DECL_C_QC(CONDITIONAL_101, conditional_101, I),
		DECL_Q_QQ(CONDITIONAL_101, conditional_101, I),

		{CONDITIONAL_111, 0, 0, "conditional_111", nullptr},//TODO: disc for conditional_111

		DECL_R_R(INCREMENT, increment, C),
		DECL_C_C(INCREMENT, increment, C),
		DECL_Q_Q(INCREMENT, increment, C),

		DECL_R_R(DECREMENT, decrement, C),
		DECL_C_C(DECREMENT, decrement, C),
		DECL_Q_Q(DECREMENT, decrement, C),

		DECL_R_R(ASSIGN, assign, C),
		DECL_C_C(ASSIGN, assign, C),
		DECL_Q_Q(ASSIGN, assign, C),

		//special kernels
		{V_INITIALIZE_PARAMETER, 0, 0, "initialize_parameter", nullptr},
		{V_TI2D_RGB, 0, 0, "ti2d_rgb", nullptr},
		{V_C2D_RGB, 0, 0, "c2d_rgb", nullptr},
		{V_C2D_RGB2, 0, 0, "c2d_rgb2", nullptr},
		//{V_INITIALIZE_CONSTANTS, 0, 0, "initialize_constants", nullptr},
	};
	CLKernel kernels_p20[]=
	{
		//special kernels
		{V_TI3D_AV, 0, 0, "ti3d_av", nullptr},
		{V_TI3D_AV_EAST, 0, 0, "ti3d_av_east", nullptr},
		{V_TI3D_AV_NORTH, 0, 0, "ti3d_av_north", nullptr},
		{V_TI3D_AV_UP, 0, 0, "ti3d_av_up", nullptr},
		{TI3D_CLASSIFYEDGES, 0, 0, "ti3d_classifyedges", nullptr},
	};
	CLKernel kernels_p21[]=
	{
		//special kernels
		{TI3D_ZEROCROSS, 0, 0, "ti3d_zerocross", nullptr},
	};
	CLKernel kernels_p22[]=
	{
		//special kernels
		{TI3D_TRGL_INDICES, 0, 0, "ti3d_trgl_indices", nullptr},
	};
	KernelDB kernel_db[]=
	{
		{kernels_p00, sizeof(kernels_p00)/sizeof(CLKernel)},
		{kernels_p01, sizeof(kernels_p01)/sizeof(CLKernel)},
		{kernels_p02, sizeof(kernels_p02)/sizeof(CLKernel)},
		{kernels_p03, sizeof(kernels_p03)/sizeof(CLKernel)},
		{kernels_p04, sizeof(kernels_p04)/sizeof(CLKernel)},
		{kernels_p05, sizeof(kernels_p05)/sizeof(CLKernel)},
		{kernels_p06, sizeof(kernels_p06)/sizeof(CLKernel)},
		{kernels_p07, sizeof(kernels_p07)/sizeof(CLKernel)},
		{kernels_p08, sizeof(kernels_p08)/sizeof(CLKernel)},
		{kernels_p09, sizeof(kernels_p09)/sizeof(CLKernel)},
		{kernels_p10, sizeof(kernels_p10)/sizeof(CLKernel)},
		{kernels_p11, sizeof(kernels_p11)/sizeof(CLKernel)},
		{kernels_p12, sizeof(kernels_p12)/sizeof(CLKernel)},
		{kernels_p13, sizeof(kernels_p13)/sizeof(CLKernel)},
		{kernels_p14, sizeof(kernels_p14)/sizeof(CLKernel)},
		{kernels_p15, sizeof(kernels_p15)/sizeof(CLKernel)},
		{kernels_p16, sizeof(kernels_p16)/sizeof(CLKernel)},
		{kernels_p17, sizeof(kernels_p17)/sizeof(CLKernel)},
		{kernels_p18, sizeof(kernels_p18)/sizeof(CLKernel)},
		{kernels_p19, sizeof(kernels_p19)/sizeof(CLKernel)},

		{kernels_p20, sizeof(kernels_p20)/sizeof(CLKernel)},//TI3D
		{kernels_p21, sizeof(kernels_p21)/sizeof(CLKernel)},
		{kernels_p22, sizeof(kernels_p22)/sizeof(CLKernel)},
	};
//	const int CLKernel_size=sizeof(CLKernel), nkernels=sizeof(kernels)/CLKernel_size;
}
double			cl_progress(std::string &ret)
{
	switch(OCL_state)
	{
	case CL_NOTHING:			ret="OpenCL API not loaded.";			break;
	case CL_LOADING_API:		ret="Loading OpenCL API...";			break;
	case CL_API_LOADED:			ret="OpenCL API loaded.";				break;
	case CL_CREATING_CONTEXT:	ret="Creating OpenCL context...";		break;
	case CL_CONTEXT_CREATED:	ret="OpenCL context created.";			break;
	case CL_COMPILING_ENTRY:	ret="Starting OpenCL compilation...";	break;

	case CL_COMPILING_PROGRAM00:
	case CL_COMPILING_PROGRAM01:
	case CL_COMPILING_PROGRAM02:
	case CL_COMPILING_PROGRAM03:
	case CL_COMPILING_PROGRAM04:
	case CL_COMPILING_PROGRAM05:
	case CL_COMPILING_PROGRAM06:
	case CL_COMPILING_PROGRAM07:
	case CL_COMPILING_PROGRAM08:
	case CL_COMPILING_PROGRAM09:
	case CL_COMPILING_PROGRAM10:
	case CL_COMPILING_PROGRAM11:
	case CL_COMPILING_PROGRAM12:
	case CL_COMPILING_PROGRAM13:
	case CL_COMPILING_PROGRAM14:
	case CL_COMPILING_PROGRAM15:
	case CL_COMPILING_PROGRAM16:
	case CL_COMPILING_PROGRAM17:
	case CL_COMPILING_PROGRAM18:
	case CL_COMPILING_PROGRAM19:
		ret="Compiling OpenCL program "+std::to_string((long long)(OCL_state+1-CL_COMPILING_PROGRAM00))
				+" of "+std::to_string((long long)nprograms)+"...";
		break;
	case CL_PROGRAMS_COMPILED:
		ret="OpenCL programs compiled.";
		break;
	case CL_RETRIEVING_BINARIES:
		ret="Retrieving OpenCL binaries.";
		break;
	case CL_LOADING_PROGRAM00:
	case CL_LOADING_PROGRAM01:
	case CL_LOADING_PROGRAM02:
	case CL_LOADING_PROGRAM03:
	case CL_LOADING_PROGRAM04:
	case CL_LOADING_PROGRAM05:
	case CL_LOADING_PROGRAM06:
	case CL_LOADING_PROGRAM07:
	case CL_LOADING_PROGRAM08:
	case CL_LOADING_PROGRAM09:
	case CL_LOADING_PROGRAM10:
	case CL_LOADING_PROGRAM11:
	case CL_LOADING_PROGRAM12:
	case CL_LOADING_PROGRAM13:
	case CL_LOADING_PROGRAM14:
	case CL_LOADING_PROGRAM15:
	case CL_LOADING_PROGRAM16:
	case CL_LOADING_PROGRAM17:
	case CL_LOADING_PROGRAM18:
	case CL_LOADING_PROGRAM19:
		ret="Loading OpenCL program "+std::to_string((long long)(OCL_state+1-CL_LOADING_PROGRAM00))
				+" of "+std::to_string((long long)nprograms)+"...";
		break;
	case CL_PROGRAMS_LOADED:
		ret="OpenCL programs loaded.";
		break;
	case CL_READY_UNTESTED:
		ret="OpenCL ready, untested.";
		break;
	case CL_TESTING:
		ret="Testing OpenCL";
		break;
	case CL_READY:
		ret="OpenCL is ready.";
		break;
	default:
		ret="Unknown state";
		return 0.5;
	}
	return (double)OCL_state/CL_READY;
}
void			cl_test()
{
	if(OCL_state<CL_READY_UNTESTED)
		return;
	//	LOGERROR("cl_test: OpenCL is not ready.");
	else
	{
		OCL_state=CL_TESTING;
#if 1
		cl_int error=0;
		//static cl_platform_id platform0=nullptr;
		//static cl_device_id device0=nullptr;
		//static cl_context context0=nullptr;
		//static cl_command_queue commandqueue0=nullptr;
		//platform0=platform;
		//device0=device;
		//context0=context;
		//commandqueue0=commandqueue;
		size_t host_sizes[3]={1, 1, 1}, host_sizes_local[3]={1, 1, 1};//DEBUG		fails after lock
		float host_args[3]={42, 1, 0};
		float testvalue=0;

		cl_mem size_buf=p_clCreateBuffer(context, CL_MEM_READ_WRITE, 2*sizeof(int), nullptr, &error);	CL_CHECK(error);
		error=p_clEnqueueWriteBuffer(commandqueue, size_buf, CL_FALSE, 0, 2*sizeof(int), host_sizes, 0, nullptr, nullptr);	CL_CHECK(error);
		//error=p_clFinish(commandqueue), CL_CHECK(error);//succeeds second time
		cl_mem args_buf=p_clCreateBuffer(context, CL_MEM_READ_WRITE, 2*sizeof(float)+sizeof(int), nullptr, &error);	CL_CHECK(error);
		error=p_clEnqueueWriteBuffer(commandqueue, args_buf, CL_FALSE, 0, 2*sizeof(float)+sizeof(int), host_args, 0, nullptr, nullptr);	CL_CHECK(error);
		//error=p_clFinish(commandqueue), CL_CHECK(error);//

		cl_mem buffer=p_clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float), nullptr, &error);	CL_CHECK(error);

		cl_kernel k_fill=kernels[V_INITIALIZE_PARAMETER];
		error=p_clSetKernelArg(k_fill, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
		error=p_clSetKernelArg(k_fill, 1, sizeof(cl_mem), &buffer);		CL_CHECK(error);
		error=p_clSetKernelArg(k_fill, 2, sizeof(cl_mem), &args_buf);	CL_CHECK(error);
		error=p_clEnqueueNDRangeKernel(commandqueue, k_fill, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
		error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_FALSE, 0, sizeof(float), &testvalue, 0, nullptr, nullptr);	CL_CHECK(error);
		error=p_clFinish(commandqueue);		CL_CHECK(error);
	//	LOGI("G2_CL: TESTVALUE = %f", testvalue);//
		if(testvalue!=42)
			LOGERROR("testvalue = %f", testvalue);
		float testvalue2=testvalue;//END DEBUG
		//testvalue+=100;
		//error=p_clEnqueueFillBuffer(commandqueue, buffer, &testvalue, sizeof(float), 0, sizeof(float), 0, nullptr, nullptr); CL_CHECK(error);//CRASH SIGSEGV
		//testvalue=0;
		//error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_FALSE, 0, sizeof(float), &testvalue, 0, nullptr, nullptr);	CL_CHECK(error);
		//error=p_clFinish(commandqueue);		CL_CHECK(error);
		//LOGI("G2_CL: TESTVALUE = %f", testvalue);//
#endif
		if(!error)
			OCL_state=CL_READY;
		//	OCL_state=CL_KERNELS_LOADED;
	}
}
void			programname2g_buf(int kp)
{
	init_directories();
	sprintf_s(g_buf, g_buf_size, "%scl_program%02d.bin", statedir, kp);
//	sprintf_s(g_buf, g_buf_size, "%s/cl_program%02d.bin", statedir, kp);
}
void			checkforbuildfailure(int error, cl_program *programs, int prog_idx, cl_device_id device, std::string &err_msg)
{
	if(error)
	{
		size_t length=0;
		error=p_clGetProgramBuildInfo(programs[prog_idx], device, CL_PROGRAM_BUILD_LOG, g_buf_size, g_buf, &length);	CL_CHECK(error);
		err_msg+="\n\n\tPROGRAM "+std::to_string((long long)prog_idx)+"\n\n"+g_buf;
		copy_to_clipboard(g_buf, length);
		messageboxa(ghWnd, "Information", "OpenCL program compilation failed. Output copied to clipboard.");
	//	LOGE("PROGRAM %d: %s", prog_idx, g_buf);
	}
}
void			cl_compile()
{
	if(OCL_state<CL_CONTEXT_CREATED)
		LOGERROR("Can't compile without OpenCL context.");
	else
	{
		OCL_state=CL_COMPILING_ENTRY;
		using namespace G2_CL;
		cl_int error=0;
		size_t retlen=0;
		std::string err_msg;

	//	const auto language=__cplusplus;//MSVC: 199711, android studio: 201402
		//std::string statedir=appdatapath;
		//statedir+='/';
		//statedir+=statefoldername;
		//mkdir_firsttime(statedir.c_str());
		char *text=nullptr;
		size_t textlen=0;
		statefolder_readstate(text, textlen);
		//{
		//	std::string str;
		//	for(unsigned k=0;k<textlen;++k)
		//	{
		//		char c=text[k];
		//		if(c>=32&&c<127)
		//			str+=c;
		//		else
		//			str+='<'+std::to_string((int)c)+'>';
		//	}
		//	//LOGI("%s: %s", statefilename, str.c_str());
		//}
		unsigned version=hexstr2uint(text);
		if(version!=g2_version.id)
		{
			messageboxa(ghWnd, "Warning", "State version mismatch.\nAbout to delete everything in the folder:\n%s", statedir);
		//	LOGI("G2_CL: Wrong state version, deleting state.");
			statefolder_deletecontents();
			statefolder_writestate();
			//const char *vstr=version2str();
			//saveFile((statedir+'/'+statefilename).c_str(), vstr, strlen(vstr));
		}
#if 0
		bool statedirexists=false;
	//	snprintf(g_buf, g_buf_size, "%s/", appdatapath);
	//	snprintf(g_buf, g_buf_size, "%s", appdatapath);
	//	snprintf(g_buf, g_buf_size, ".");
		DIR *d=opendir(appdatapath);	SYS_CHECK();//https://stackoverflow.com/questions/4204666/how-to-list-files-in-a-directory-in-a-c-program/17683417
		if(d)
		{
			struct dirent *dir;
			while(dir=readdir(d))
			{
				SYS_CHECK();
				LOGI("G2_CL: %s", dir->d_name);
				if(!strcmp(dir->d_name, statedir))
				{
					statedirexists=true;
					break;
				}
			}
			SYS_CHECK();
			closedir(d);	SYS_CHECK();
		}
		if(!statedirexists)
		{
			snprintf(g_buf, g_buf_size, "%s/%s", appdatapath, statedir);
			mode_t x;
			auto result=mkdir(g_buf, );
		}
#endif
#if 0
		if(binaries_valid)
		{
			for(int kp=0;kp<nprograms;++kp)
			{
				OCL_state=CL_COMPILING_PROGRAM00+kp;
				auto &bin=binaries[kp];
				int status=0;
				programs[kp]=p_clCreateProgramWithBinary(context, 1, &device, &bin.size, (const unsigned char**)&bin.bin, &status, &error);	CL_CHECK(error);
				error=p_clBuildProgram(programs[kp], 0, nullptr, "", nullptr, nullptr);		CL_CHECK(error);
				checkforbuildfailure(error, programs, kp, device, err_msg);
			}
			OCL_state=CL_PROGRAMS_COMPILED;
		}
		else//compile from source first
#endif
		{
			double t_start=time_sec();
			int success=0;
			int common_size=0;
			const char *rc_common=loadresource(IDR_CL_COMMON, TEXTFILE, common_size);
			std::string common(rc_common, common_size);
			for(int kp=0;kp<nprograms;++kp)
		//	for(int kp=19;kp<nprograms;++kp)//
			{
				OCL_state=CL_COMPILING_PROGRAM00+kp;
				ProgramBinary bin={};
				programname2g_buf(kp);
				std::string programpath=g_buf;
			//	snprintf(g_buf, g_buf_size, "%s/cl_program%02d.bin", statedir, kp);
				double t2=time_sec();
				if(loadbinary&&!loadFile(programpath.c_str(), (char*&)bin.bin, bin.size))//binary file is there
				{
					set_window_title("Loading program %d/%d... (%.2lf sec)", kp+1, nprograms, t2-t_start);
					int status=0;
					programs[kp]=p_clCreateProgramWithBinary(context, 1, &device, &bin.size, (const unsigned char**)&bin.bin, &status, &error);	CL_CHECK(error);
					error=p_clBuildProgram(programs[kp], 0, nullptr, "", nullptr, nullptr);					CL_CHECK(error);
					checkforbuildfailure(error, programs, kp, device, err_msg);
					bin.free();
				}
				else//compile from source
				{
					set_window_title("Compiling program %d/%d... (%.2lf sec)", kp+1, nprograms, t2-t_start);
					const char *sources[2]={};
					size_t srclen[2]={};
					int nsources=0;

					int srck_size=0;
					auto rc_srck=loadresource(IDR_CL_SRC01+kp, TEXTFILE, srck_size);
					std::string srck(rc_srck, srck_size);
					if(kp<n_g2programs)//add header for G2 main programs
						sources[0]=common.c_str(),	srclen[0]=common.size(),		++nsources;
					sources[nsources]=srck.c_str(),	srclen[nsources]=srck.size(),	++nsources;
				//	if(kp<n_g2programs)//add header for G2 main programs
				//		sources[0]=CLSource::program_common,	srclen[0]=strlen(CLSource::program_common), ++nsources;
				//	sources[nsources]=CLSource::programs[kp],	srclen[nsources]=strlen(CLSource::programs[kp]), ++nsources;
					//if(kp<n_g2programs)
					//{
					//	sources[0]=CLSource::program_common,	srclen[0]=strlen(CLSource::program_common);
					//	sources[1]=CLSource::programs[kp],		srclen[1]=strlen(CLSource::programs[kp]);
					//	nsources=2;
					//}
					//else
					//{
					//	sources[1]=CLSource::programs[kp],		srclen[1]=strlen(CLSource::programs[kp]);
					//	nsources=1;
					//}
				//	if(kp==19)
				//	//if(kp==22)
				//		int LOL_1=0;
					const char *c_options="";//compilation options
					if(kp==19&&cl_gl_interop&&OCL_version>120&&p_clCreateFromGLTexture)
						c_options=" -D G2_OCL_IMAGES ";
					G2_CL::programs[kp]=p_clCreateProgramWithSource(context, nsources, sources, srclen, &error);	CL_CHECK(error);
					//G2_CL::programs[kp]=p_clCreateProgramWithSource(context, 1, sources, srclen, &error);			CL_CHECK(error);
					error=p_clBuildProgram(programs[kp], 0, nullptr, c_options, nullptr, nullptr);					CL_CHECK(error);//calls abort(1)
				//	error=p_clBuildProgram(programs[kp], 0, nullptr, " -cl-opt-disable -Werror -cl-std=CL1.1 ", nullptr, nullptr);			CL_CHECK(error);//amdocl.dll: U.E. A.V. reading 0x14D3A98C, 0x147AF43C
				//	error=p_clBuildProgram(programs[kp], 0, nullptr, "-cl-opt-disable", nullptr, nullptr);			CL_CHECK(error);//amdocl.dll: unhandled exception: access violation, unrelated to try/catch
				//	error=p_clBuildProgram(programs[kp], 0, nullptr, "", nullptr, nullptr);							CL_CHECK(error);//calls abort(1)
					checkforbuildfailure(error, programs, kp, device, err_msg);
					if(!error)//save binary
					{
						error=p_clGetProgramInfo(programs[kp], CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin.size, &retlen);	CL_CHECK(error);
						bin.bin=(unsigned char*)malloc(bin.size);
						error=p_clGetProgramInfo(programs[kp], CL_PROGRAM_BINARIES, sizeof(unsigned char*), &bin.bin, &retlen);	CL_CHECK(error);

						programname2g_buf(kp);
					//	snprintf(g_buf, g_buf_size, "%s/cl_program%02d.bin", statedir, kp);
						error=saveFile(g_buf, (char*&)bin.bin, bin.size);		CL_CHECK(error);
						bin.free();
					}
				}
			}
			OCL_state=CL_PROGRAMS_COMPILED;
			set_window_title("");
			if(!err_msg.empty())
				copy_to_clipboard(err_msg.c_str(), err_msg.size());
#if 0
			else//build successful, retrieve binaries
			{
				OCL_state=CL_RETRIEVING_BINARIES;
				for(int kp=0;kp<nprograms;++kp)
				{
					auto &bin=binaries[kp];
					error=p_clGetProgramInfo(programs[kp], CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin.size, &retlen);	CL_CHECK(error);
					bin.bin=(unsigned char*)malloc(bin.size);
					error=p_clGetProgramInfo(programs[kp], CL_PROGRAM_BINARIES, sizeof(unsigned char*), &bin.bin, &retlen);	CL_CHECK(error);

					snprintf(g_buf, g_buf_size, "cl_program%02d.bin", kp);
					saveFile(g_buf, bin.bin, bin.size);
					//std::ofstream file(g_buf, std::ios::out|std::ios::binary);
					//if(file.is_open())
					//{
					//	file.write((char*)bin.bin, bin.size);
					//	file.close();
					//}
				}
				if(!error)
					binaries_valid=true;
			}
#endif
		}
		if(!error)//build successful: extract kernel handles
		{
			for(int kp=0;kp<nprograms;++kp)
			{
				OCL_state=CL_LOADING_PROGRAM00+kp;
				auto &kernel_batch=kernel_db[kp];
				for(int kk=0;kk<kernel_batch.nkernels;++kk)
				{
					auto &kernel=kernel_batch.kernels[kk];
					if(!kernels[kernel.idx]&&kernel.name)
					{
						kernels[kernel.idx]=p_clCreateKernel(programs[kp], kernel.name, &error);		CL_CHECK(error);
						if(error)
							LOGERROR("Failed to retrieve kernel %d - %s from program %d.", kk, kernel.name, kp);
						//	LOGE("Error line %d:\n\n\tFailed to retrieve kernel %d from program %d:\t%s\n\n", __LINE__, kk, kp, kernel.name);
					}
					//TODO: disc_idx
				}
			}
			OCL_state=CL_PROGRAMS_LOADED;
		}
		if(!error)
			OCL_state=CL_READY_UNTESTED;
	}
}
void 			cl_reset(){}//TODO: delete, recompile & save cl programs when user types reset
void 			cl_initiate()
{
	load_OpenCL_API();
	using namespace G2_CL;
	static_assert(sizeof(kernel_db)/sizeof(KernelDB)==nprograms, "kernel_db size is wrong");

	if(OCL_state<CL_API_LOADED)
		LOGERROR("Can't create context without OpenCL API");
	else
	{
		OCL_state=CL_CREATING_CONTEXT;
		cl_int error=0;
		size_t retlen=0;
		//static bool firsttime=true;
		//if(firsttime)
		//{
		//	firsttime=false;
		unsigned n_platforms=0, n_devices=0;
		// Fetch the Platform and Device IDs; we only want one.		//https://donkey.vernier.se/~yann/hellocl.c
		if(!platform)
			{error=p_clGetPlatformIDs(1, &platform, &n_platforms);							CL_CHECK(error);}
		if(!device)
			{error=p_clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, &n_devices);	CL_CHECK(error);}//was CL_DEVICE_TYPE_ALL

		error=p_clGetPlatformInfo(platform, CL_PLATFORM_VERSION, g_buf_size, g_buf, &retlen);	CL_CHECK(error);
		if(!OCL_version)
		{
			int k=0;
			for(;g_buf[k]>='A'&&g_buf[k]<='Z'||g_buf[k]>='a'&&g_buf[k]<='z';++k);//skip 'OpenCL'
			for(;g_buf[k]==' '||g_buf[k]=='\t'||g_buf[k]=='\r'||g_buf[k]=='\n';++k);//skip whitespace
			if(g_buf[k]>='0'&&g_buf[k]<='9')
			{
				OCL_version=(g_buf[k]-'0')*100;
				if(g_buf[k+1]=='.'&&k+2<(int)retlen)
					OCL_version+=(g_buf[k+2]-'0')*10;
			}
		}
	//	LOGI("%s", g_buf);
		error=p_clGetDeviceInfo(device, CL_DEVICE_NAME, g_buf_size, g_buf, &retlen);	CL_CHECK(error);	//query device info
	//	LOGI("\n\n\t%*s\n\n", (int)retlen, g_buf);
		error=p_clGetDeviceInfo(device, CL_DEVICE_VENDOR, g_buf_size, g_buf, &retlen);	CL_CHECK(error);
	//	LOGI("\n\n\t%*s\n\n", (int)retlen, g_buf);
		error=p_clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &g_maxlocalsize, &retlen);	CL_CHECK(error);
	//	LOGI("\n\n\tApp resolution: w=%d, h=%d\n\n", w, h);
	//	LOGI("\n\n\tCL_DEVICE_MAX_WORK_GROUP_SIZE = %d\n\n", (int)g_maxlocalsize);
		error=p_clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, g_buf_size, g_buf, &retlen);	CL_CHECK(error);
	//	LOGI("\n\n\t%*s\n\n", (int)retlen, g_buf);
		std::string extensions=g_buf;
		auto extpos=extensions.find("cl_khr_gl_sharing");
		cl_gl_interop=extpos!=std::string::npos;
	//	if(!cl_gl_interop)
	//		LOGI("cl_khr_gl_sharing not supported");
		//}

		if(!context)
		{
			cl_context_properties properties[8]={};
			if(cl_gl_interop)
			{
				auto gl_context=wglGetCurrentContext();
				auto hDC=wglGetCurrentDC();
				properties[0]=CL_GL_CONTEXT_KHR,	properties[1]=(cl_context_properties)gl_context,	//Win32		https://stackoverflow.com/questions/26802905/getting-opengl-buffers-using-opencl
				properties[2]=CL_WGL_HDC_KHR,		properties[3]=(cl_context_properties)hDC,
				properties[4]=CL_CONTEXT_PLATFORM,	properties[5]=(cl_context_properties)platform,//OpenCL platform object
				properties[6]=0,					properties[7]=0;
				
			//	auto gl_context=eglGetCurrentContext();//changes when resuming
			//	auto egl_display=eglGetCurrentDisplay();
			//	properties[0]=CL_GL_CONTEXT_KHR,	properties[1]=(cl_context_properties)gl_context;	//Android
			//	properties[2]=CL_EGL_DISPLAY_KHR,	properties[3]=(cl_context_properties)egl_display;
			//	properties[4]=CL_CONTEXT_PLATFORM,	properties[5]=(cl_context_properties)platform;
			//	properties[6]=0, properties[7]=0;
			}
			else
			{
				properties[0]=CL_CONTEXT_PLATFORM, properties[1]=(cl_context_properties)platform;
				properties[2]=0, properties[3]=0;
			}
			context=p_clCreateContext(properties, 1, &device, nullptr, nullptr, &error);	CL_CHECK(error);
		}
		else
		{
			error=p_clRetainContext(context);												CL_CHECK(error);
		}
		commandqueue=p_clCreateCommandQueue(context, device, 0, &error);				CL_CHECK(error);
		OCL_state=CL_CONTEXT_CREATED;
		cl_compile();
	//	std::thread thr_cl_init(cl_compile);
	//	thr_cl_init.detach();
#if 0
		{
			error=p_clGetPlatformIDs(1, &platform, &n_platforms);							CL_CHECK(error);//
			error=p_clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, &n_devices);	CL_CHECK(error);//
			cl_context_properties properties[]=//https://stackoverflow.com/questions/26802905/getting-opengl-buffers-using-opencl
			{
				CL_GL_CONTEXT_KHR,   (cl_context_properties)eglGetCurrentContext(),
				CL_EGL_DISPLAY_KHR,  (cl_context_properties)eglGetCurrentDisplay(),
				CL_CONTEXT_PLATFORM, (cl_context_properties)platform, // OpenCL platform object
				0, 0,
			};
			context=p_clCreateContext(properties, 1, &device, nullptr, nullptr, &error);	CL_CHECK(error);//
			commandqueue=p_clCreateCommandQueue(context, device, 0, &error);				CL_CHECK(error);//
		}
#endif
		//if(!error)
		//	OCL_state=CL_CONTEXT_CREATED;
#if 0
		float testvalue=42;//DEBUG		always succeeds
		cl_mem mem_in1=p_clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float), nullptr, &error);	CL_CHECK(error);
		error=p_clEnqueueWriteBuffer(commandqueue, mem_in1, CL_FALSE, 0, sizeof(float), &testvalue, 0, nullptr, nullptr);	CL_CHECK(error);
		testvalue=0;
		error=p_clEnqueueReadBuffer(commandqueue, mem_in1, CL_FALSE, 0, sizeof(float), &testvalue, 0, nullptr, nullptr);	CL_CHECK(error);
		error=p_clFinish(commandqueue);		CL_CHECK(error);//END DEBUG
		//if(NEWPOINTER(platform)||NEWPOINTER(device)||NEWPOINTER(context)||NEWPOINTER(commandqueue))
		//	int LOL_1=0;
#endif

#if 0
		//hellocl.c		//https://donkey.vernier.se/~yann/hellocl.c
		cl_int error=0;
		cl_platform_id platform;
		cl_device_id device;
		unsigned n_platforms, n_devices;
		// Fetch the Platform and Device IDs; we only want one.
		error=p_clGetPlatformIDs(1, &platform, &n_platforms);						CL_CHECK(error);//get available OpenCL platforms & devices
		error=p_clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, &n_devices);	CL_CHECK(error);
		cl_context_properties properties[]=
		{
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)platform,
			0
		};
		cl_context context=p_clCreateContext(properties, 1, &device, nullptr, nullptr, &error);	CL_CHECK(error);//create OpenCL context
		cl_command_queue cq=p_clCreateCommandQueue(context, device, 0, &error);					CL_CHECK(error);

		//compile program & extract kernels
		const char src[]=R"CLSRC(
__kernel void myfunc(__global const float *in1, __global const float *in2, __global float *out)
{
	const unsigned idx=get_global_id(0);
	out[idx]=in1[idx]+in2[idx];
//	out[idx]=get_global_id(0);
}
)CLSRC";
		const char *sources[]=
		{
			src,
		};
		size_t srclen[]=
		{
			strlen(src),
		};
		cl_program program=p_clCreateProgramWithSource(context, 1, sources, srclen, &error);CL_CHECK(error);
		error=p_clBuildProgram(program, 0, nullptr, "", nullptr, nullptr);					CL_CHECK(error);
		if(error)
		{
			size_t length=0;
			error=p_clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, g_buf_size, g_buf, &length);	CL_CHECK(error);
		}
		cl_kernel k_add=p_clCreateKernel(program, "myfunc", &error);						CL_CHECK(error);

		//use kernel
	//	glFlush();
		for(int k=0;k<worksize;++k)//initialize	TODO: ititialize in kernel
		{
			in1[k]=k;
			in2[k]=worksize-1-k;
			out[k]=rand();//
		}
		cl_mem mem_in1=p_clCreateBuffer(context, CL_MEM_READ_ONLY, worksize*sizeof(float), nullptr, &error);	CL_CHECK(error);//create OpenCL buffers (or use OpenGL buffers)
		cl_mem mem_in2=p_clCreateBuffer(context, CL_MEM_READ_ONLY, worksize*sizeof(float), nullptr, &error);	CL_CHECK(error);
		cl_mem mem_out=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, worksize*sizeof(float), nullptr, &error);	CL_CHECK(error);
		error=p_clSetKernelArg(k_add, 0, sizeof(cl_mem), &mem_in1); CL_CHECK(error);//set arguments
		error=p_clSetKernelArg(k_add, 1, sizeof(cl_mem), &mem_in2); CL_CHECK(error);
		error=p_clSetKernelArg(k_add, 2, sizeof(cl_mem), &mem_out); CL_CHECK(error);
		error=p_clEnqueueWriteBuffer(commandqueue, mem_in1, CL_FALSE, 0, worksize*sizeof(float), in1, 0, nullptr, nullptr);	CL_CHECK(error);//send arguments
		error=p_clEnqueueWriteBuffer(commandqueue, mem_in2, CL_FALSE, 0, worksize*sizeof(float), in2, 0, nullptr, nullptr);	CL_CHECK(error);
		error=p_clEnqueueNDRangeKernel(commandqueue, k_add, 1, nullptr, &worksize, &worksize, 0, nullptr, nullptr);			CL_CHECK(error);//execute
		error=p_clEnqueueReadBuffer(commandqueue, mem_out, CL_FALSE, 0, worksize*sizeof(float), out, 0, nullptr, nullptr);	CL_CHECK(error);
		error=p_clFinish(commandqueue);	CL_CHECK(error);
		unload_OpenCL_API();//
#endif
	}
}
void			cl_terminate()
{
	if(OCL_state>=CL_API_LOADED)
	{
		int error=p_clReleaseContext(context);	CL_CHECK(error);
	}
}

struct			CLTerm
{
	char mathSet;
	cl_mem r, i, j, k;
	CLTerm():mathSet(0), r(nullptr), i(nullptr), j(nullptr), k(nullptr){}
};
struct 			DebugBuffer
{
	cl_mem buffer;
	const char *name;
};
void			debug_printbuffer(cl_mem buffer, int Xplaces, int Yplaces, int Zplaces, int stride, const char *bufname=nullptr)
{
	int nfloats=Xplaces*Yplaces*Zplaces;
	auto result=(float*)malloc(nfloats*sizeof(float));
	int error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_TRUE, 0, nfloats*sizeof(float), result, 0, nullptr, nullptr);	CL_CHECK(error);
	std::stringstream LOL_1;
	if(bufname)
		LOL_1<<bufname<<"\r\n";
	for(int kz=0;kz<Zplaces;kz+=stride)
	{
		LOL_1<<"z="<<kz<<'\t';
		for(int kx=0;kx<Xplaces;kx+=stride)
			LOL_1<<"\t["<<kx<<']';
		LOL_1<<"\r\n";
		for(int ky=0;ky<Yplaces;ky+=stride)
		{
			LOL_1<<"["<<ky<<"] ";
			for(int kx=0;kx<Xplaces;kx+=stride)
			{
				int idx=Xplaces*(Yplaces*kz+ky)+kx;
				LOL_1<<result[idx]<<"="<<std::hex<<(int&)result[idx]<<", ";
			}
			LOL_1<<"\r\n";
		}
		LOL_1<<"\r\n";
	}
	//for(int k=0;k<nfloats;k+=stride)
	//	str+=bufname?bufname:"buffer", str+='[', str+=std::to_string(k), str+="] = ", str+=std::to_string((double)result[k]), str+="\r\n";
	//	//LOGI("G2_CL: %s[%d] = %g", bufname?bufname:"buffer", k, (double)result[k]);
	copy_to_clipboard(LOL_1.str());
	free(result);
}
void			debug_printbuffers(DebugBuffer *buffers, int nbuffers, int Xplaces, int Yplaces, int stride)
//void			debug_printbuffers(DebugBuffer *buffers, int nbuffers, int nfloats, int stride)
{
	int nfloats=Xplaces*Yplaces;
	auto result=(float*)malloc(nbuffers*nfloats*sizeof(float));
	for(int k=0;k<nbuffers;++k)
	{
		int error=p_clEnqueueReadBuffer(commandqueue, buffers[k].buffer, CL_TRUE, 0, nfloats*sizeof(float), result+k*nfloats, 0, nullptr, nullptr);
		CL_CHECK(error);
	}
	std::stringstream LOL_1;
	for(int kx=0;kx<Xplaces;kx+=stride)
		LOL_1<<"\t["<<kx<<']';
	LOL_1<<"\r\n";
	for(int ky=0;ky<Yplaces;ky+=stride)
	{
		LOL_1<<"["<<ky<<"] ";
		for(int kx=0;kx<Xplaces;kx+=stride)
		{
			int idx=Xplaces*ky+kx;
		//	LOL_1<<"["<<kx<<", "<<ky<<"] ";
		//	LOL_1<<"G2_CL: ["<<k<<"] ";
			for(int kb=0;kb<nbuffers;++kb)
			{
				const char *name=buffers[kb].name;
				LOL_1<<(name?name:"buffer")<<": "<<result[nfloats*kb+idx]<<"="<<std::hex<<(int&)result[nfloats*kb+idx]<<", ";
			}
		//	LOL_1<<"\r\n";
		//	LOGI("%s", LOL_1.str().c_str());
		}
		LOL_1<<"\r\n";
	}
	copy_to_clipboard(LOL_1.str());
	free(result);
}
static size_t	host_sizes[3]={},//{globalwidth, globalheight, globaldepth}
				host_sizes_local[3]={};//{localwidth, localheight, localdepth}
//static float 	Xstart=0, Xsample=0,
//				Yend=0, mYsample=0;//negative Ysample for C2D
void 			create_resize_buffer(cl_mem &buffer, unsigned nfloats)
{
	int error=0;
	buffer=p_clCreateBuffer(context, CL_MEM_READ_WRITE, nfloats*sizeof(float), nullptr, &error);
	CL_CHECK(error);
}
#ifdef DEBUG3
int				g_Xplaces=0, g_Yplaces=0, g_Zplaces=0;//DEBUG
//const char	*g_bufname=nullptr;
#endif
void			initialize_const_comp(cl_mem &buffer, unsigned nfloats, float value, cl_mem size_buf, cl_mem args_buf)
{
	int error=0;
	if(!buffer)
		create_resize_buffer(buffer, nfloats);

	float host_args[3]={value, 0, 0};
	error=p_clEnqueueWriteBuffer(commandqueue, args_buf, CL_FALSE, 0, 3*sizeof(float), host_args, 0, nullptr, nullptr);	CL_CHECK(error);
	cl_kernel k_fill=kernels[V_INITIALIZE_PARAMETER];
	error=p_clSetKernelArg(k_fill, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
	error=p_clSetKernelArg(k_fill, 1, sizeof(cl_mem), &buffer);		CL_CHECK(error);
	error=p_clSetKernelArg(k_fill, 2, sizeof(cl_mem), &args_buf);	CL_CHECK(error);
	error=p_clEnqueueNDRangeKernel(commandqueue, k_fill, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
#ifdef BLOCKING_INITIALIZATION
//	error=p_clFlush(commandqueue);		CL_CHECK(error);//wrong results
	error=p_clFinish(commandqueue);		CL_CHECK(error);//correct results
#endif

//	float pattern[4]={value};
//	error=p_clEnqueueFillBuffer(commandqueue, buffer, &value, sizeof(float), 0, nfloats*sizeof(float), 0, nullptr, nullptr); CL_CHECK(error);//CRASH SIGSEGV on Galaxy S5
	//error=p_clEnqueueFillBuffer(commandqueue, buffer, pattern, sizeof(float), 0, sizeof(float), 0, nullptr, nullptr); CL_CHECK(error);
#if 0
	std::vector<float> args(3);
	error=p_clEnqueueReadBuffer(commandqueue, args_buf, CL_TRUE, 0, 3*sizeof(float), args.data(), 0, nullptr, nullptr);	CL_CHECK(error);
	args[0]=args[0];
	std::vector<float> debug_buf(nfloats);
	error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_TRUE, 0, nfloats*sizeof(float), debug_buf.data(), 0, nullptr, nullptr);	CL_CHECK(error);
	debug_buf[0]=debug_buf[0];
	std::stringstream LOL_1;
	LOL_1<<value<<"\r\n";
	LOL_1<<args[0]<<' '<<args[1]<<' '<<args[2]<<"\r\n";
	LOL_1<<debug_buf[0]<<"\r\n";
	copy_to_clipboard(LOL_1.str());
#endif
#ifdef DEBUG2
	float testval=0;
	error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_TRUE, 0, 1*sizeof(float), &testval, 0, nullptr, nullptr);	CL_CHECK(error);
	LOGI("G2_CL: const = %f", testval);
#endif
}
void 			initialize_component(ModeParameters const &mp, cl_mem &buffer, unsigned nfloats, char varType, float value, float time, cl_mem size_buf, cl_mem args_buf)
{
	int dimension=-1;//-1: fill with constant 'value'
	float data=0;
	switch(varType)
	{
	case 'x':dimension=0;break;
	case 'y':dimension=1;break;
	case 'z':dimension=2;break;
	case 'c':data=value;break;
	case 't':data=time;break;
	default:dimension=3;break;
	}
	if(dimension==-1)
		initialize_const_comp(buffer, nfloats, data, size_buf, args_buf);
	else if(dimension<3)
//	else if(dimension<2)//C2D: exactly 2 input dimensions
	{
		if(!buffer)
			create_resize_buffer(buffer, nfloats);
		int error=0;
		float host_args[3];
		if(dimension==0)
			host_args[0]=(float)mp.cx, host_args[1]=(float)mp.mx;
		else if(dimension==1)
			host_args[0]=(float)mp.cy, host_args[1]=(float)mp.my;
		else
			host_args[0]=(float)mp.cz, host_args[1]=(float)mp.mz;
		//if(dimension==0)
		//	host_args[0]=Xstart, host_args[1]=Xsample;
		//else
		//	host_args[0]=Yend, host_args[1]=mYsample;
		host_args[2]=(float)dimension;
		error=p_clEnqueueWriteBuffer(commandqueue, args_buf, CL_FALSE, 0, 3*sizeof(float), host_args, 0, nullptr, nullptr);	CL_CHECK(error);
		cl_kernel k_fill=kernels[V_INITIALIZE_PARAMETER];
		error=p_clSetKernelArg(k_fill, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
		error=p_clSetKernelArg(k_fill, 1, sizeof(cl_mem), &buffer);		CL_CHECK(error);
		error=p_clSetKernelArg(k_fill, 2, sizeof(cl_mem), &args_buf);	CL_CHECK(error);
		int error0=error;
		error=p_clEnqueueNDRangeKernel(commandqueue, k_fill, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
#ifdef BLOCKING_INITIALIZATION
	//	error=p_clFlush(commandqueue);		CL_CHECK(error);//wrong results
		error=p_clFinish(commandqueue);		CL_CHECK(error);//correct results
#endif
#ifdef COPY_INITIALIZATION
		const char bufname[]={varType, '\0'};
		debug_printbuffer(buffer, g_Xplaces, g_Yplaces, 1, 1, bufname);
		messageboxa(ghWnd, "Information", "%s sent to clipboard", bufname?bufname:"Data");//
#endif
#ifdef DEBUG2
		LOGI("G2_CL: error: %d -> %d, lsizes={%d, %d, %d}", error0, error, (int)host_sizes_local[0], (int)host_sizes_local[1], (int)host_sizes_local[2]);
#endif
		//if(blocking)
		//	{error=p_clFinish(commandqueue);		CL_CHECK(error);}
#ifdef DEBUG2
		float testval=0;
		error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_TRUE, 0, 1*sizeof(float), &testval, 0, nullptr, nullptr);	CL_CHECK(error);
		LOGI("G2_CL: buffer[0] = %f", testval);
#endif
	}
}
void			cl_free_buffer(cl_mem &buffer)
{
	if(buffer)
	{
		int error=p_clReleaseMemObject(buffer);
		CL_CHECK(error);
		if(!error)
			buffer=nullptr;
	}
}
#if 0
inline float	max(float a, float b){return (a+b+abs(a-b))*0.5f;}
inline float 	rsqrt(float x){return 1.f/sqrt(x);}
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
void ti2d_rgb(int kx, int ky, const int *size, const float *xr, const float *curvecolor, int *rgb)
{//size{Xplaces, Yplaces, 1, appwidth, appheight}
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
//	//write_imagef(rgb, (int2)(kx, ky), (float4)(curvecolor[0], curvecolor[1], curvecolor[2], alpha));
//	write_imagef(rgb, (int2)(kx, ky), (float4)(alpha*curvecolor[0], alpha*curvecolor[1], alpha*curvecolor[2], 1));//plain white
//	//write_imagef(rgb, (int2)(kx, ky), (float4)((1-alpha)*curvecolor[0], (1-alpha)*curvecolor[1], (1-alpha)*curvecolor[2], 1));//points exactly on curve
	unsigned char c=(unsigned char)(255*(1-alpha));
//	unsigned char c=(unsigned char)(255*alpha);
	rgb[size[3]*ky+kx]=c<<16|c<<8|c;
}
#endif
#if 1
void			colorFunction_bcw(Comp1d const &v, int *color)
{
	const double r=v.real(), i=v.imag();
	if(r!=r||i!=i)
		*color=0x7F7F7F;
	else if(abs(r)==_HUGE||abs(i)==_HUGE)
		*color=0xFFFFFF;
	else
	{
		double threshold=10, inv_th=1/threshold;//1			//black->color->white		81.98 cycles/px
		const double cos_pi_6=0.866025403784439, sin_pi_6=0.5;
		double hyp=sqrt(r*r+i*i), cosx=r/hyp, sinx=i/hyp,
			mag=255*exp(-hyp*G2::_ln2*inv_th);
		double red=1+cosx*cos_pi_6-sinx*sin_pi_6, green=1+sinx, blue=1+cosx*-cos_pi_6-sinx*sin_pi_6;
		if(hyp<threshold)
			mag=255-mag, red*=mag, green*=mag, blue*=mag;
		else
			red=255-mag*(2-red), green=255-mag*(2-green), blue=255-mag*(2-blue);
		auto p=(unsigned char*)color;
		p[0]=(unsigned char)red, p[1]=(unsigned char)green, p[2]=(unsigned char)blue, p[3]=0xFF;//0xAABBGGRR	//argb
	}
}
#endif
inline unsigned	dimension_work_size(unsigned places, size_t maxlocaldim)
{
	places-=places%maxlocaldim;//when maxlocaldim isn't a power of 2
	return places>0?places:maxlocaldim;

//	return places-places%maxlocaldim;
//	return places&~(g_maxlocaldim-1);
}
inline unsigned	dimension_work_size(unsigned places, size_t maxlocaldim, char reschange, unsigned maxdim)
{
	places+=maxlocaldim*reschange;
	if(places>maxdim)
		places=maxdim;
	places-=places%maxlocaldim;//when maxlocaldim isn't a power of 2
	return places>0?places:maxlocaldim;
}
inline void		read_buffer(cl_mem buffer, float *&data, unsigned nfloats)
{
	if(buffer)
	{
		data=(float*)malloc(nfloats*sizeof(float));
		int error=p_clEnqueueReadBuffer(commandqueue, buffer, CL_TRUE, 0, nfloats*sizeof(float), data, 0, nullptr, nullptr);	CL_CHECK(error);
	}
}
namespace		SW
{
	double
		dpxr=0, dpxi=0, dpxj=0, dpxk=0,//op1
		dpyr=0, dpyi=0, dpyj=0, dpyk=0,//op2
		dprr=0, dpri=0, dprj=0, dprk=0;//res
	inline VectP in1(float *xr)
	{
		dpxr=*xr;
		return &dpxr;
	}
	inline CompP in1(float *xr, float *xi)
	{
		dpxr=*xr, dpxi=*xi;
		return CompP(&dpxr, &dpxi);
	}
	inline QuatP in1(float *xr, float *xi, float *xj, float *xk)
	{
		dpxr=*xr, dpxi=*xi, dpxj=*xj, dpxk=*xk;
		return QuatP(&dpxr, &dpxi, &dpxj, &dpxk);
	}
	inline VectP in2(float *yr)
	{
		dpyr=*yr;
		return &dpyr;
	}
	inline CompP in2(float *yr, float *yi)
	{
		dpyr=*yr, dpyi=*yi;
		return CompP(&dpyr, &dpyi);
	}
	inline QuatP in2(float *yr, float *yi, float *yj, float *yk)
	{
		dpyr=*yr, dpyi=*yi, dpyj=*yj, dpyk=*yk;
		return QuatP(&dpyr, &dpyi, &dpyj, &dpyk);
	}
	inline void out(float *rr){*rr=(float)dprr;}
	inline void out(float *rr, float *ri){*rr=(float)dprr, *ri=(float)dpri;}
	inline void out(float *rr, float *ri, float *rj, float *rk){*rr=(float)dprr, *ri=(float)dpri, *rj=(float)dprj, *rk=(float)dprk;}
}
void			cl_setsizes(int mode_idx, int *Xplaces, int *Yplaces, int *Zplaces, char reschange)
{
	int nd=0, maxX=0, maxY=0, maxZ=0;
	switch(mode_idx)
	{
	case MODE_I1D:
	case MODE_T1D:
	case MODE_T1D_C:
	case MODE_T1D_H:
		{
			int log_mld=floor_log2(g_maxlocalsize);
			g_maxlocalX=1<<(log_mld>>1), g_maxlocalY=1, g_maxlocalZ=1;
			maxX=(1<<30)-1, maxY=1, maxZ=1;
			nd=1;
		}
		break;
	case MODE_I2D:
	case MODE_T2D:
	case MODE_C2D:
	case MODE_T2D_H:
		{
			int log_mld=floor_log2(g_maxlocalsize);
			g_maxlocalY=g_maxlocalX=1<<(log_mld>>1), g_maxlocalZ=1;
			maxX=65536, maxY=65536, maxZ=1;
			nd=2;
		}
		break;
	case MODE_I3D:
	case MODE_C3D:
		{
			int log_mld=floor_log2(g_maxlocalsize);
			int log_xy=log_mld/3+(log_mld%3!=0), log_z=log_mld-(log_xy<<1);//log_xy >= log_z
			g_maxlocalY=g_maxlocalX=1<<log_xy, g_maxlocalZ=1<<log_z;//g_maxlocalX & g_maxlocalY are multiples of g_maxlocalZ
			maxX=512, maxY=512, maxZ=512;
			nd=3;
		}
		break;
	}
	if(nd)
	{
		*Xplaces=dimension_work_size(*Xplaces, g_maxlocalX, reschange, maxX);
		if(nd>=2)
		{
			*Yplaces=dimension_work_size(*Yplaces, g_maxlocalY, reschange, maxY);
			if(nd>=3)
				*Zplaces=dimension_work_size(*Zplaces, g_maxlocalZ, reschange, maxZ);
		}
	}
}
typedef unsigned long long ulong;
union			WorkIdx
{
	struct{unsigned short ke, kx, ky, kz;};
	struct{unsigned kw, kt;};
	unsigned long long idx;
	WorkIdx():idx(0){}
	WorkIdx(short kx, short ky, short kz, short ke):kx(kx), ky(ky), kz(kz), ke(ke){}
	WorkIdx(unsigned kw, unsigned kt):kw(kw), kt(kt){}
	void set(short kx, short ky, short kz, short ke){this->kx=kx, this->ky=ky, this->kz=kz, this->ke=ke;}
	void set(unsigned kw, unsigned kt){this->kw=kw, this->kt=kt;}
};
enum			ExEdgeBitIdx//54edges
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
enum			ExEdgeBitIdxRev//54edges
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
const ulong		obsmask[8]=//original bit selection mask
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
const char*		edgenumber2str(int ke)
{
	const char *a="???";
	switch(ke)
	{
		//case-string
#define	CS(name)	case name:a=#name;break;
		//main (new) edge indices
		CS(DSW_DSE)	CS(DSW_DNW)	CS(DSW_USW)
		
		CS(DSW_mmW)	CS(USW_mmW)	CS(UNW_mmW)	CS(DNW_mmW)
		CS(DSW_mSm)	CS(USW_mSm)	CS(USE_mSm)	CS(DSE_mSm)
		CS(DSW_Dmm)	CS(DSE_Dmm)	CS(DNE_Dmm)	CS(DNW_Dmm)
		
		CS(mmW_Dmm)	CS(mmW_mSm)	CS(mmW_Umm)	CS(mmW_mNm)
		CS(Umm_mNm)	CS(mNm_Dmm)	CS(Dmm_mSm)	CS(mSm_Umm)
		CS(mmE_Dmm)	CS(mmE_mSm)	CS(mmE_Umm)	CS(mmE_mNm)
		
		CS(mmm_Dmm)	CS(mmm_mSm)	CS(mmm_mmW)	CS(mmm_Umm)	CS(mmm_mNm)	CS(mmm_mmE)

		//extended (redundant) edge indices
		CS(DSE_DNE)	CS(DNE_DNW)	CS(DSE_USE)	CS(DNE_UNE)	CS(DNW_UNW)	CS(USW_USE)	CS(USE_UNE)	CS(UNE_UNW)	CS(UNW_USW)

		CS(mmE_DSE)	CS(mmE_USE)	CS(mmE_UNE)	CS(mmE_DNE)
		CS(mNm_DNW)	CS(mNm_UNW)	CS(mNm_UNE)	CS(mNm_DNE)
		CS(Umm_USW)	CS(Umm_USE)	CS(Umm_UNE)	CS(Umm_UNW)
#undef	CS
	}
	return a;
}
//int hammingweight(ulong x)//https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
//{
//	x-=x>>1&0x5555555555555555;
//	x=(x&0x3333333333333333)+(x>>2&0x3333333333333333);
//	return ((x+(x>>4))&0x0F0F0F0F0F0F0F0F)*0x0101010101010101>>56;
//}
#ifdef V6_CPU
typedef unsigned long long ulong;
static const unsigned char nt[16]=//number of triangles produced from tetrahedron given 4 bit vertex condition
{//	00	01	10	11	lo/hi
	0,	1,	1,	2,//00
	1,	2,	2,	1,//01
	1,	2,	2,	1,//10
	2,	1,	1,	0,//11
};
char hammingweight(unsigned long long x)//https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
{
	x-=x>>1&0x5555555555555555;
	x=(x&0x3333333333333333)+(x>>2&0x3333333333333333);
	return ((x+(x>>4))&0x0F0F0F0F0F0F0F0F)*0x0101010101010101>>56;
}
//union			WorkIdx
//{
//	struct{unsigned short ke, kx, ky, kz;};
//	struct{unsigned kw, kt;};
//	unsigned long long idx;
//	WorkIdx():idx(0){}
//	WorkIdx(short kx, short ky, short kz, short ke):kx(kx), ky(ky), kz(kz), ke(ke){}
//	WorkIdx(unsigned kw, unsigned kt):kw(kw), kt(kt){}
//	void set(short kx, short ky, short kz, short ke){this->kx=kx, this->ky=ky, this->kz=kz, this->ke=ke;}
//	void set(unsigned kw, unsigned kt){this->kw=kw, this->kt=kt;}
//};
//bool operator<(WorkIdx const &a, WorkIdx const &b){return a.idx<b.idx;}
//enum ExEdgeBitIdx//54edges
//{
//	//main (new) edge indices
//	DSW_DSE, DSW_DNW, DSW_USW,
//
//	DSW_mmW, USW_mmW, UNW_mmW, DNW_mmW,
//	DSW_mSm, USW_mSm, USE_mSm, DSE_mSm,
//	DSW_Dmm, DSE_Dmm, DNE_Dmm, DNW_Dmm,
//
//	mmW_Dmm, mmW_mSm, mmW_Umm, mmW_mNm,
//	Umm_mNm, mNm_Dmm, Dmm_mSm, mSm_Umm,
//	mmE_Dmm, mmE_mSm, mmE_Umm, mmE_mNm,
//
//	mmm_Dmm, mmm_mSm, mmm_mmW, mmm_Umm, mmm_mNm, mmm_mmE,
//
//	//extended (redundant) edge indices
//	DSE_DNE, DNE_DNW, DSE_USE, DNE_UNE, DNW_UNW, USW_USE, USE_UNE, UNE_UNW, UNW_USW,
//	mmE_DSE, mmE_USE, mmE_UNE, mmE_DNE,
//	mNm_DNW, mNm_UNW, mNm_UNE, mNm_DNE,
//	Umm_USW, Umm_USE, Umm_UNE, Umm_UNW,
//};
//enum ExEdgeBitIdxRev//54edges
//{
//	//main (new) edge indices
//	DSE_DSW, DNW_DSW, USW_DSW,
//
//	mmW_DSW, mmW_USW, mmW_UNW, mmW_DNW,
//	mSm_DSW, mSm_USW, mSm_USE, mSm_DSE,
//	Dmm_DSW, Dmm_DSE, Dmm_DNE, Dmm_DNW,
//
//	Dmm_mmW, mSm_mmW, Umm_mmW, mNm_mmW,
//	mNm_Umm, Dmm_mNm, mSm_Dmm, Umm_mSm,
//	Dmm_mmE, mSm_mmE, Umm_mmE, mNm_mmE,
//
//	Dmm_mmm, mSm_mmm, mmW_mmm, Umm_mmm, mNm_mmm, mmE_mmm,
//
//	//extended (redundant) edge indices
//	DNE_DSE, DNW_DNE, USE_DSE, UNE_DNE, UNW_DNW, USE_USW, UNE_USE, UNW_UNE, USW_UNW,
//	DSE_mmE, USE_mmE, UNE_mmE, DNE_mmE,
//	DNW_mNm, UNW_mNm, UNE_mNm, DNE_mNm,
//	USW_Umm, USE_Umm, UNE_Umm, UNW_Umm,
//};
//const ulong		obsmask[8]=//original bit selection mask
//{
//	0x00000001FFFFFFFF,//all coords inside
//	(1ULL<<DSE_USE)|(1ULL<<DSE_DNE)|(1ULL<<mmE_DSE)|(1ULL<<mmE_USE)|(1ULL<<mmE_UNE)|(1ULL<<mmE_DNE)|0x00000001FFFFFFFF,//x at end: DSE's originals
//	(1ULL<<DNW_UNW)|(1ULL<<DNW_DNE)|(1ULL<<mNm_DNW)|(1ULL<<mNm_UNW)|(1ULL<<mNm_UNE)|(1ULL<<mNm_DNE)|0x00000001FFFFFFFF,//y at end: DNW's originals
//	(1ULL<<DNE_UNE) | (1ULL<<DSE_USE)|(1ULL<<DSE_DNE)|(1ULL<<mmE_DSE)|(1ULL<<mmE_USE)|(1ULL<<mmE_UNE)|(1ULL<<mmE_DNE) | (1ULL<<DNW_UNW)|(1ULL<<DNW_DNE)|(1ULL<<mNm_DNW)|(1ULL<<mNm_UNW)|(1ULL<<mNm_UNE)|(1ULL<<mNm_DNE)|0x00000001FFFFFFFF,//x&y at end: (DSE, DNW & DNE)'s originals
//
//	(1ULL<<USW_USE)|(1ULL<<USW_UNW)|(1ULL<<Umm_USW)|(1ULL<<Umm_USE)|(1ULL<<Umm_UNE)|(1ULL<<Umm_UNW)|0x00000001FFFFFFFF,//z at end: USW's originals
//	(1ULL<<USE_UNE) | (1ULL<<DSE_USE)|(1ULL<<DSE_DNE)|(1ULL<<mmE_DSE)|(1ULL<<mmE_USE)|(1ULL<<mmE_UNE)|(1ULL<<mmE_DNE) | (1ULL<<USW_USE)|(1ULL<<USW_UNW)|(1ULL<<Umm_USW)|(1ULL<<Umm_USE)|(1ULL<<Umm_UNE)|(1ULL<<Umm_UNW)|0x00000001FFFFFFFF,//x&z at end: (DSE, USW & USE)'s originals
//	(1ULL<<UNW_UNE) | (1ULL<<DNW_UNW)|(1ULL<<DNW_DNE)|(1ULL<<mNm_DNW)|(1ULL<<mNm_UNW)|(1ULL<<mNm_UNE)|(1ULL<<mNm_DNE) | (1ULL<<USW_USE)|(1ULL<<USW_UNW)|(1ULL<<Umm_USW)|(1ULL<<Umm_USE)|(1ULL<<Umm_UNE)|(1ULL<<Umm_UNW)|0x00000001FFFFFFFF,//y&z at end: (DNW, USW & UNW)'s originals
//	0x003FFFFFFFFFFFFF,//x,y&z at end: entire cube
//};
struct WE_Offset//work-edge offset, used by CPU-side part 2
{
	short dx, dy, dz;
	char ke2[8];//{normal jump, (xyz, yz, xz, z, xy, y, x)@end}
//	char ke2[2];//{normal jump, @end (use original)}
};
static const WE_Offset we_offsets[54-33]=//work-edge offsets, for redundant bits, access: [ke-33]={Xoffset, Yomask, Zomask, ke2[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)]}
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
//static const WE_Offset we_offsets[54-33]=//work-edge offsets, for redundant bits, access: [ke-33]={Xoffset, Yomask, Zomask, ke2[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)]}
//{//(dx, ym, zm)	normal jump, (x,		y,			xy,			z,			xz,			yz,			xyz)@end
//	{1,  0,  0,		{}},//{DSW_DNW,	DSE_DNE,	DSW_DNW,	DSE_DNE,	DSW_DNW,	DSE_DNE,	DSW_DNW,	DSE_DNE}},//DSE_DNE
//	{0, -1,	 0,		{}},//{DSW_DSE,	DSW_DSE,	DNE_DNW,	DNE_DNW,	DSW_DSE,	DSW_DSE,	DNE_DNW,	DNE_DNW}},//DNE_DNW
//	{1,  0,  0,		{}},//{DSW_USW,	DSE_USE,	DSW_USW,	DSE_USE,	DSW_USW,	DSE_USE,	DSW_USW,	DSE_USE}},//DSE_USE
//	{1, -1,	 0,		{}},//{DSW_USW,	DSE_USE,	DNW_UNW,	DNE_UNE,	DSW_USW,	DSE_USE,	DNW_UNW,	DNE_UNE}},//DNE_UNE//2D jump
//	{0, -1,	 0,		{}},//{DSW_USW,	DSW_USW,	DNW_UNW,	DNW_UNW,	DSW_USW,	DSW_USW,	DNW_UNW,	DNW_UNW}},//DNW_UNW
//	{0,  0, -1,		{}},//{DSW_DSE,	DSW_DSE,	DSW_DSE,	DSW_DSE,	USW_USE,	USW_USE,	USW_USE,	USW_USE}},//USW_USE
//	{1,  0, -1,		{}},//{DSW_DNW,	DSE_DNE,	DSW_DNW,	DSE_DNE,	USW_UNW,	USE_UNE,	USW_UNW,	USE_UNE}},//USE_UNE//2D jump
//	{0, -1,	-1,		{}},//{DSW_DSE,	DSW_DSE,	DNE_DNW,	DNE_DNW,	USW_USE,	USW_USE,	UNE_UNW,	UNE_UNW}},//UNE_UNW//2D jump
//	{0,  0, -1,		{}},//{DSW_DNW,	DSW_DNW,	DSW_DNW,	DSW_DNW,	UNW_USW,	UNW_USW,	UNW_USW,	UNW_USW}},//UNW_USW
//
//	{1,  0,  0,		{}},//{mmW_DSW,	mmE_DSE,	mmW_DSW,	mmE_DSE,	mmW_DSW,	mmE_DSE,	mmW_DSW,	mmE_DSE}},//mmE_DSE
//	{1,  0,  0,		{}},//{mmW_USW,	mmE_USE,	mmW_USW,	mmE_USE,	mmW_USW,	mmE_USE,	mmW_USW,	mmE_USE}},//mmE_USE
//	{1,  0,  0,		{}},//{mmW_UNW,	mmE_UNE,	mmW_UNW,	mmE_UNE,	mmW_UNW,	mmE_UNE,	mmW_UNW,	mmE_UNE}},//mmE_UNE
//	{1,  0,  0,		{}},//{mmW_DNW,	mmE_DNE,	mmW_DNW,	mmE_DNE,	mmW_DNW,	mmE_DNE,	mmW_DNW,	mmE_DNE}},//mmE_DNE
//
//	{0, -1,  0,		{}},//{mSm_DSW,	mSm_DSW,	mNm_DNW,	mNm_DNW,	mSm_DSW,	mSm_DSW,	mNm_DNW,	mNm_DNW}},//mNm_DNW
//	{0, -1,  0,		{}},//{mSm_USW,	mSm_USW,	mNm_UNW,	mNm_UNW,	mSm_USW,	mSm_USW,	mNm_UNW,	mNm_UNW}},//mNm_UNW
//	{0, -1,  0,		{}},//{mSm_USE,	mSm_USE,	mNm_UNE,	mNm_UNE,	mSm_USE,	mSm_USE,	mNm_UNE,	mNm_UNE}},//mNm_UNE
//	{0, -1,  0,		{}},//{mSm_DSE,	mSm_DSE,	mNm_DNE,	mNm_DNE,	mSm_DSE,	mSm_DSE,	mNm_DNE,	mNm_DNE}},//mNm_DNE
//
//	{0,  0, -1,		{}},//{Dmm_DSW,	Dmm_DSW,	Dmm_DSW,	Dmm_DSW,	Umm_USW,	Umm_USW,	Umm_USW,	Umm_USW}},//Umm_USW
//	{0,  0, -1,		{}},//{Dmm_DSE,	Dmm_DSE,	Dmm_DSE,	Dmm_DSE,	Umm_USE,	Umm_USE,	Umm_USE,	Umm_USE}},//Umm_USE
//	{0,  0, -1,		{}},//{Dmm_DNE,	Dmm_DNE,	Dmm_DNE,	Dmm_DNE,	Umm_UNE,	Umm_UNE,	Umm_UNE,	Umm_UNE}},//Umm_UNE
//	{0,  0, -1,		{}},//{Dmm_DNW,	Dmm_DNW,	Dmm_DNW,	Dmm_DNW,	Umm_UNW,	Umm_UNW,	Umm_UNW,	Umm_UNW}},//Umm_UNW
//};
//static const WE_Offset we_offsets[54-33]=//work-edge offsets, for redundant bits, access: [ke-33]={Xoffset, Yomask, Zomask, ke2[?]}
//{//	(dx, ym, zm)	normal jump, @end (original)
//	{1,  0,  0,		{DSW_DNW,	DSE_DNE}},//DSE_DNE
//	{0, -1,	 0,		{DSW_DSE,	DNE_DNW}},//DNE_DNW
//	{1,  0,  0,		{DSW_USW,	DSE_USE}},//DSE_USE
//	{1, -1,	 0,		{DSW_USW,	DNE_UNE}},//DNE_UNE
//	{0, -1,	 0,		{DSW_USW,	DNW_UNW}},//DNW_UNW
//	{0,  0, -1,		{DSW_DSE,	USW_USE}},//USW_USE
//	{1,  0, -1,		{DSW_DNW,	USE_UNE}},//USE_UNE//<-
//	{0, -1,	-1,		{DSW_DSE,	UNE_UNW}},//UNE_UNW
//	{0,  0, -1,		{DSW_DNW,	UNW_USW}},//UNW_USW
//
//	{1,  0,  0,		{mmW_DSW,	mmE_DSE}},//mmE_DSE
//	{1,  0,  0,		{mmW_USW,	mmE_USE}},//mmE_USE
//	{1,  0,  0,		{mmW_UNW,	mmE_UNE}},//mmE_UNE
//	{1,  0,  0,		{mmW_DNW,	mmE_DNE}},//mmE_DNE
//
//	{0, -1,  0,		{mSm_DSW,	mNm_DNW}},//mNm_DNW
//	{0, -1,  0,		{mSm_USW,	mNm_UNW}},//mNm_UNW
//	{0, -1,  0,		{mSm_USE,	mNm_UNE}},//mNm_UNE
//	{0, -1,  0,		{mSm_DSE,	mNm_DNE}},//mNm_DNE
//
//	{0,  0, -1,		{Dmm_DSW,	Umm_USW}},//Umm_USW
//	{0,  0, -1,		{Dmm_DSE,	Umm_USE}},//Umm_USE
//	{0,  0, -1,		{Dmm_DNE,	Umm_UNE}},//Umm_UNE
//	{0,  0, -1,		{Dmm_DNW,	Umm_UNW}},//Umm_UNW
//};
#define		TEHE(V0, V1, V2, V3)		V0##_##V1, V0##_##V2, V0##_##V3, V1##_##V2, V1##_##V3, V2##_##V3
static const unsigned char ti[28*6]=//tetrahedron edge indices
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
inline char		getbit(long long x, char bit){return x>>bit&1;}
struct			float2
{
	float x, y;
	void set(float x, float y){this->x=x, this->y=y;}
};
struct			float3
{
	float x, y, z;
	void set(float x, float y, float z){this->x=x, this->y=y, this->z=z;}
};
struct			float4
{
	float x, y, z, w;
	void set(float x, float y, float z, float w){this->x=x, this->y=y, this->z=z, this->w=w;}
};
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
static const unsigned char e2v1[54]=//edge index to first vertex (data point) index
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
static const unsigned char e2v2[54]=//edge index to second vertex (data point) index
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
		ret.set(-c/b, 0, INFINITY, 0);
	//	ret=(float4)(-c/b, 0, INFINITY, 0);
	else if(abs(b)/abs(a)>=lim&&abs(c)/abs(a)>=lim)
		ret.set(-c/b, 0, -b/a, 0);
	//	ret=(float4)(-c/b, 0, -b/a, 0);
	else
	{
		float first=-b, disc2=b*b-4*a*c;
		float disc=sqrt(abs(disc2));
		float inv_den=0.5f/a;
		if(disc2<0)
		{
			first*=inv_den, disc*=inv_den;
			ret.set(first, disc, first, -disc);
			//ret=(float4)(first, disc, first, -disc);
		}
		else
			ret.set((first+disc)*inv_den, 0, (first-disc)*inv_den, 0);
			//ret=(float4)(first+disc, 0, first-disc, 0);

		//b/=a, c/=a;
		//float first=-0.5f*b, disc2=b*b-4*c;
		//float disc=0.5f*sqrt(abs(disc2));
		//if(disc2<0)
		//	ret.set(first, disc, first, -disc);
		//	//ret=(float4)(first, disc, first, -disc);
		//else
		//	ret.set(first+disc, 0, first-disc, 0);
		//	//ret=(float4)(first+disc, 0, first-disc, 0);
	}
	return ret;
}
float2			solvecubic(float a, float b, float c, float d)
{
	float2 ret;
	//if(a==0&&b==0&&c==0)
	//{
	//	ret.x=(Vx-A.V)/(B.V-A.V);
	//	ret.y=1;
	//	return;
	//}
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
		//	s=s<0?-pow(-s, 1.f/3):pow(s, 1.f/3);
			t=r-sqrt(disc);
			t=cbrt(t);
		//	t=t<0?-pow(-t, 1.f/3):pow(t, 1.f/3);
			r1=-term1+s+t;//The first root is always real
			term1+=(s+t)/2;
			float term2=_sqrt3*(-t+s)/2;
			r2r=r3r=-term1, r2i=term2, r3i=-term2;
		}
		else if(disc==0)//The remaining options are all real
		{
			r13=cbrt(r);
		//	float sign=(r>0)-(r<0);
		//	r13=sign*pow(sign*r, 1.f/3);
		//	r13=r<0?-pow(-r, 1.f/3):pow(r, 1.f/3);
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
	else
	{
		float results[3]={r1, r2r, r3r};
		int nresults=1+(r2i==0)+(r3i==0);
		int smallest=0;
		for(int k=1;k<nresults;++k)
			if(abs(results[smallest])>abs(results[k]))
				smallest=k;
		ret.set(clamp01(results[smallest]), 0);
	}
//	else								ret.x=-1, ret.y=0;
	return ret;
}
float3			mix(float3 const &a, float3 const &b, float alpha)
{
	float3 ret;
	ret.x=a.x+(b.x-a.x)*alpha;
	ret.y=a.y+(b.y-a.y)*alpha;
	ret.z=a.z+(b.z-a.z)*alpha;
	return ret;
}
void			vstore3(float3 const &result, int idx, float *a){a[idx]=result.x, a[idx+1]=result.y, a[idx+2]=result.z;}
#endif
void			exndr2clipboard(const float *ndr5, int Xplaces, int Yplaces, int Zplaces)
{
	const char *endl="\r\n";
	std::stringstream LOL_1;
	LOL_1<<
		"ExNDR\r\n"
		"kx = \t";
	for(int kx=0;kx<Xplaces;++kx)
		LOL_1<<'\t'<<kx;
	LOL_1<<endl<<endl;
	for(int kz=0;kz<Zplaces;++kz)
	{
		LOL_1<<"kz = "<<kz<<endl;
		for(int ky=0;ky<Yplaces;++ky)
		{
			LOL_1<<"(kx, "<<ky<<", "<<kz<<')';
			for(int kx=0;kx<Xplaces;++kx)
				LOL_1<<'\t'<<ndr5[Xplaces*(Yplaces*kz+ky)+kx];
			LOL_1<<endl;
		}
		LOL_1<<endl;
	}
	int XYplaces=Xplaces*Yplaces, ndrSize=Xplaces*Yplaces*Zplaces, exndrSize=ndrSize*5;
	LOL_1<<"\r\nAverages: down face, south face, west face, center\r\n";
	for(int k=ndrSize;k+3<exndrSize;k+=4)
	{
		int kw=(k-ndrSize)>>2, kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
		LOL_1<<kw<<'('<<kx<<", "<<ky<<", "<<kz<<"):\t"<<ndr5[kw]<<'\t'<<ndr5[kw+1]<<'\t'<<ndr5[kw+2]<<'\t'<<ndr5[kw+3]<<endl;
	}
	copy_to_clipboard(LOL_1.str());
	messageboxa(ghWnd, "Information", "Extended NDR copied to clipboard.");
}
void			classifyedges2clipboard(const ulong *edgeinfo, const int *nvert, const int *ntrgl, int Xplaces, int Yplaces, int Zplaces)
{
	const char *endl="\r\n";
	std::stringstream LOL_1;
	LOL_1<<
		"classifyedges output\r\n"
		"kw\t(kx, ky, kz)\tedgeinfo[kw]\t\tnvert[kw]\tntrgl[kw]\r\n";
	int XYplaces=Xplaces*Yplaces, ndrSize=XYplaces*Zplaces;
	for(int kw=0;kw<ndrSize;++kw)
	{
		int kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
		int length=0;
		if(kx>=Xplaces-1||ky>=Yplaces-1||kz>=Zplaces-1)
			length+=sprintf_s(g_buf+length, g_buf_size-length, "\t");
	//	int length=sprintf_s(g_buf, g_buf_size, "%d\t0x%016LLX\t", kw, edgeinfo[kw]);
		length+=sprintf_s(g_buf+length, g_buf_size-length, "%d\t(%d, %d, %d)\t0x%016llX\t", kw, kx, ky, kz, edgeinfo[kw]);
		int nvk=nvert[kw];
		length+=sprintf_s(g_buf+length, g_buf_size-length, nvk<0||nvk>54?"0x%08X\t":"%d\t", nvk);
		int ntk=ntrgl[kw];
		length+=sprintf_s(g_buf+length, g_buf_size-length, ntk<0||ntk>100?"0x%08X":"%d", ntk);
	//	length+=sprintf_s(g_buf, g_buf_size, "%d\t0x%016LLX\t%d\t%d\r\n", kw, edgeinfo[kw], nvert[kw], ntrgl[kw]);
		LOL_1<<g_buf;
		if(kx<Xplaces-1||ky<Yplaces-1||kz<Zplaces-1)
		{
			LOL_1<<'\t';
			auto eik=edgeinfo[kw];
			for(int ke=0;ke<64;++ke)
				if(eik>>ke&1)
					LOL_1<<' '<<edgenumber2str(ke)<<',';
		}
		LOL_1<<endl;
	//	LOL_1<<g_buf<<endl;
		if(!((kw+1)%Xplaces))
		{
			LOL_1<<endl;
			if(!((kw+1)%XYplaces))
				LOL_1<<endl;
		}
	}
	copy_to_clipboard(LOL_1.str());
	messageboxa(ghWnd, "Information", "ti3d_classifyedges() output copied to clipboard.");
}
void 			cl_solve(Expression const &ex, ModeParameters const &mp, double time, ...)
{//expression -> OpenGL texture
	if(OCL_state<CL_READY_UNTESTED)
		return;
	//	LOGERROR("Solve: OCL_state = %s", cl_state2str(OCL_state));
	else
	{
		if(OCL_state==CL_READY_UNTESTED)
			cl_test();
		//{//size_t								//GA70	GS5	Win32
		//	auto sp=sizeof(void*);				//8		4	4
		//	auto ssize=sizeof(size_t);			//8		4	4
		//	auto schar=sizeof(char);			//1		1	1
		//	auto sshort=sizeof(short);			//2		2	2
		//	auto swchar_t=sizeof(wchar_t);		//4		4	2
		//	auto sint=sizeof(int);				//4		4	4
		//	auto suint=sizeof(unsigned);		//4		4	4
		//	auto slong=sizeof(long);			//8		4	4
		//	auto sllong=sizeof(long long);		//8		8	8
		//	auto sfloat=sizeof(float);			//4		4	4
		//	auto sdouble=sizeof(double);		//8		8	8
		//	auto sldouble=sizeof(long double);	//16	8	8
		//	int LOL_1=0;
		//}//
#ifdef DEBUG3
		g_Xplaces=mp.Xplaces, g_Yplaces=mp.Yplaces, g_Zplaces=mp.Zplaces;
#endif
		prof_add("cl_solve entry");
		glFlush();
		prof_add("glFlush");
		unsigned ndrSize=mp.Xplaces*mp.Yplaces*mp.Zplaces;
		//unsigned ndrSize=Xplaces*Yplaces;
		//Xstart=float(VX-DX*0.5), Xsample	=float(DX/Xplaces);
		//Yend	=float(VY+DY*0.5), mYsample	=float(-DY/Yplaces);
		auto ftime=(float)time;
		int error=0;
		cl_mem size_buf=p_clCreateBuffer(context, CL_MEM_READ_ONLY, 5*sizeof(int), nullptr, &error);	CL_CHECK(error);
		cl_mem args_buf=p_clCreateBuffer(context, CL_MEM_READ_ONLY, 3*sizeof(float), nullptr, &error);	CL_CHECK(error);
		unsigned host_sizes32[6]={};
		switch(mp.mode_idx)//determine local work dimensions
		{
		case MODE_I1D:
		case MODE_T1D:
		case MODE_T1D_C:
		case MODE_T1D_H:
			{
				int log_mld=floor_log2(g_maxlocalsize);
				g_maxlocalX=1<<(log_mld>>1), g_maxlocalY=1, g_maxlocalZ=1;
			}
		//	//g_maxlocaldim=g_maxlocaldim;
		//	g_maxlocaldim=(size_t)sqrt(g_maxlocalsize);
		//	g_maxlocaldim=1<<floor_log2(g_maxlocaldim);
			host_sizes32[0]=dimension_work_size(mp.Xplaces, g_maxlocalX), host_sizes32[1]=1, host_sizes32[2]=1;
			host_sizes32[3]=mp.Xplaces, host_sizes32[4]=1, host_sizes32[5]=1;
			host_sizes[0]=host_sizes32[0], host_sizes[1]=host_sizes32[1], host_sizes[2]=host_sizes32[2];
			host_sizes_local[0]=g_maxlocalX, host_sizes_local[1]=g_maxlocalY, host_sizes_local[2]=g_maxlocalZ;
			break;
		case MODE_I2D:
		case MODE_T2D:
		case MODE_C2D:
		case MODE_T2D_H:
			{
				int log_mld=floor_log2(g_maxlocalsize);
				g_maxlocalY=g_maxlocalX=1<<(log_mld>>1), g_maxlocalZ=1;
			}
		//	g_maxlocaldim=(size_t)sqrt(g_maxlocalsize);
		//	g_maxlocaldim=1<<floor_log2(g_maxlocaldim);
			host_sizes32[0]=dimension_work_size(mp.Xplaces, g_maxlocalX), host_sizes32[1]=dimension_work_size(mp.Yplaces, g_maxlocalY), host_sizes32[2]=1;
			host_sizes32[3]=mp.Xplaces, host_sizes32[4]=mp.Yplaces, host_sizes32[5]=1;
			host_sizes[0]=host_sizes32[0], host_sizes[1]=host_sizes32[1], host_sizes[2]=host_sizes32[2];
			host_sizes_local[0]=g_maxlocalX, host_sizes_local[1]=g_maxlocalY, host_sizes_local[2]=g_maxlocalZ;
			break;
		case MODE_I3D:
		case MODE_C3D:
			{
				int log_mld=floor_log2(g_maxlocalsize);
				int log_xy=log_mld/3+(log_mld%3!=0), log_z=log_mld-(log_xy<<1);
				g_maxlocalY=g_maxlocalX=1<<log_xy, g_maxlocalZ=1<<log_z;
			}
		//	g_maxlocaldim=(size_t)cbrt(g_maxlocalsize);
		//	//g_maxlocaldim=1<<floor_log2(g_maxlocaldim);
			host_sizes32[0]=dimension_work_size(mp.Xplaces, g_maxlocalX), host_sizes32[1]=dimension_work_size(mp.Yplaces, g_maxlocalY), host_sizes32[2]=dimension_work_size(mp.Zplaces, g_maxlocalZ);
			host_sizes32[3]=mp.Xplaces, host_sizes32[4]=mp.Yplaces, host_sizes32[5]=mp.Zplaces;
			host_sizes[0]=host_sizes32[0], host_sizes[1]=host_sizes32[1], host_sizes[2]=host_sizes32[2];
			host_sizes_local[0]=g_maxlocalX, host_sizes_local[1]=g_maxlocalY, host_sizes_local[2]=g_maxlocalZ;
			break;
		}
		host_sizes32[0]=host_sizes32[3], host_sizes32[1]=host_sizes32[4], host_sizes32[2]=host_sizes32[5];
		error=p_clEnqueueWriteBuffer(commandqueue, size_buf, CL_FALSE, 0, 5*sizeof(int), host_sizes32, 0, nullptr, nullptr);	CL_CHECK(error);//send size buffer
#ifdef DEBUG2
		LOGI("G2_CL: hsizes={%d, %d, %d, %d, %d}", host_sizes32[0], host_sizes32[1], host_sizes32[2], host_sizes32[3], host_sizes32[4]);
		int sizes[5]={};
		error=p_clEnqueueReadBuffer(commandqueue, size_buf, CL_TRUE, 0, 5*sizeof(int), sizes, 0, nullptr, nullptr);	CL_CHECK(error);
		LOGI("G2_CL: sizes={%d, %d, %d, %d, %d}", sizes[0], sizes[1], sizes[2], sizes[3], sizes[4]);
#endif
		int nterms=ex.n.size();
		std::vector<CLTerm> terms(nterms);
		prof_add("prep");
		for(int kn=0;kn<nterms;++kn)//initialize terms
		{
			auto &term=terms[kn];
			auto &n=ex.n[kn];
			auto &val=ex.data[kn];
			if(n.constant)
			{
				term.mathSet=n.mathSet;
				initialize_const_comp(term.r, ndrSize, val.r.toFloat(), size_buf, args_buf);
				if(n.mathSet>='c')
				{
					initialize_const_comp(term.i, ndrSize, val.i.toFloat(), size_buf, args_buf);
					if(n.mathSet>='h')
					{
						initialize_const_comp(term.j, ndrSize, val.j.toFloat(), size_buf, args_buf);
						initialize_const_comp(term.k, ndrSize, val.k.toFloat(), size_buf, args_buf);
					}
				}
			}
			else//variable term
			{
				auto &var=ex.variables[n.varNo];
				term.mathSet=var.mathSet;
				initialize_component(mp, term.r, ndrSize, var.varTypeR, val.r.toFloat(), ftime, size_buf, args_buf);
				if(n.mathSet>='c')
				{
					initialize_component(mp, term.i, ndrSize, var.varTypeI, val.i.toFloat(), ftime, size_buf, args_buf);
					if(n.mathSet>='h')
					{
						initialize_component(mp, term.j, ndrSize, var.varTypeJ, val.j.toFloat(), ftime, size_buf, args_buf);
						initialize_component(mp, term.k, ndrSize, var.varTypeK, val.k.toFloat(), ftime, size_buf, args_buf);
					}
				}
			}
		}
#if 0
		if(terms.size()==3)
		{
			DebugBuffer buffers[]=
			{
				{terms[0].r, "x"},
				{terms[1].r, "y"},
				{terms[2].r, "z"},

			//	{terms[0].r, "x"},
			//	{terms[1].i, "i"},
			//	{terms[2].r, "y"},
			};
			debug_printbuffers(buffers, sizeof(buffers)/sizeof(DebugBuffer), mp.Xplaces, mp.Yplaces, 1);
		//	debug_printbuffers(buffers, sizeof(buffers)/sizeof(DebugBuffer), mp.Xplaces, mp.Yplaces, mp.Xplaces/10);
		}
#endif
		prof_add("initialization");
		for(int i=0, nInstr=ex.i.size();i<nInstr;++i)//solve
		{
			auto &in=ex.i[i];
			auto kernel=kernels[in.cl_idx];
			if(mp.mode_idx==MODE_I2D&&i==nInstr-1&&ex.resultLogicType>=2)//logic result: subtract NDRs at the end
			{
				switch(in.type)
				{
				case SIG_NOOP:								kernel=kernels[CL_NOOP];	break;
				case SIG_R_R:case SIG_C_R:					kernel=kernels[R_R_MINUS];	break;
				case SIG_C_C:case SIG_R_C:					kernel=kernels[C_C_MINUS];	break;
				case SIG_Q_Q:case SIG_C_Q:case SIG_R_Q:		kernel=kernels[Q_Q_MINUS];	break;
				case SIG_R_RR:case SIG_C_RR:				kernel=kernels[R_RR_MINUS];	break;
				case SIG_C_RC:case SIG_R_RC:				kernel=kernels[C_RC_MINUS];	break;
				case SIG_Q_RQ:case SIG_R_RQ:				kernel=kernels[Q_RQ_MINUS];	break;
				case SIG_C_CR:case SIG_R_CR:				kernel=kernels[C_CR_MINUS];	break;
				case SIG_C_CC:case SIG_R_CC:				kernel=kernels[C_CC_MINUS];	break;
				case SIG_Q_CQ:case SIG_R_CQ:				kernel=kernels[Q_CQ_MINUS];	break;
				case SIG_Q_QR:case SIG_R_QR:				kernel=kernels[Q_QR_MINUS];	break;
				case SIG_Q_QC:case SIG_R_QC:case SIG_C_QC:	kernel=kernels[Q_QC_MINUS];	break;
				case SIG_Q_QQ:case SIG_R_QQ:				kernel=kernels[Q_QQ_MINUS];	break;

				default://unreachable
				//case SIG_INLINE_IF:
				//case SIG_CALL:
				//case SIG_BIF:
				//case SIG_BIN:
				//case SIG_JUMP:
				//case SIG_RETURN:
					kernel=kernels[CL_NOOP];
					break;
				}
			}
			if(kernel)//TODO: user function -> OpenCL kernel		TODO: inline if
			{
				if(in.type<SIG_INLINE_IF)
				{
					int arg_number=0;
					error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &size_buf);	CL_CHECK(error);	++arg_number;
					auto &res=terms[in.result], &op1=terms[in.op1];
					error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &res.r);	CL_CHECK(error);	++arg_number;
					if(in.r_ms>='c')
					{
						error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &res.i);	CL_CHECK(error);	++arg_number;
						if(in.r_ms>='h')
						{
							error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &res.j);	CL_CHECK(error);	++arg_number;
							error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &res.j);	CL_CHECK(error);	++arg_number;
						}
					}
					error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op1.r);	CL_CHECK(error);	++arg_number;
					if(in.op1_ms>='c')
					{
						error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op1.i);	CL_CHECK(error);	++arg_number;
						if(in.op1_ms>='h')
						{
							error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op1.j);	CL_CHECK(error);	++arg_number;
							error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op1.j);	CL_CHECK(error);	++arg_number;
						}
					}
					if(in.op2_ms>='R')
					{
						auto &op2=terms[in.op2];
						error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op2.r);	CL_CHECK(error);	++arg_number;
						if(in.op2_ms>='c')
						{
							error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op2.i);	CL_CHECK(error);	++arg_number;
							if(in.op2_ms>='h')
							{
								error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op2.j);	CL_CHECK(error);	++arg_number;
								error=p_clSetKernelArg(kernel, arg_number, sizeof(cl_mem), &op2.j);	CL_CHECK(error);//	++arg_number;
							}
						}
					}
				}
				else if(in.type==SIG_INLINE_IF)
				{
					auto &res=terms[in.result], &op1=terms[in.op1], &op2=terms[in.op2], &op3=terms[in.op3];
					error=p_clSetKernelArg(kernel,  0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);

					error=p_clSetKernelArg(kernel,  1, sizeof(cl_mem), &res.r);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel,  2, sizeof(cl_mem), &res.i);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel,  3, sizeof(cl_mem), &res.j);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel,  4, sizeof(cl_mem), &res.k);		CL_CHECK(error);

					error=p_clSetKernelArg(kernel,  5, sizeof(cl_mem), &op1.r);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel,  6, sizeof(cl_mem), &op1.i);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel,  7, sizeof(cl_mem), &op1.j);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel,  8, sizeof(cl_mem), &op1.k);		CL_CHECK(error);

					error=p_clSetKernelArg(kernel,  9, sizeof(cl_mem), &op2.r);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 10, sizeof(cl_mem), &op2.i);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 11, sizeof(cl_mem), &op2.j);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 12, sizeof(cl_mem), &op2.k);		CL_CHECK(error);

					error=p_clSetKernelArg(kernel, 13, sizeof(cl_mem), &op3.r);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 14, sizeof(cl_mem), &op3.i);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 15, sizeof(cl_mem), &op3.j);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 16, sizeof(cl_mem), &op3.k);		CL_CHECK(error);
				}
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				//if(blocking)
				//	{error=p_clFinish(commandqueue);		CL_CHECK(error);}
			}
			else//kernel missing, solve in software
#if 1
			{
				float
					*rr=nullptr, *ri=nullptr, *rj=nullptr, *rk=nullptr,
					*xr=nullptr, *xi=nullptr, *xj=nullptr, *xk=nullptr,
					*yr=nullptr, *yi=nullptr, *yj=nullptr, *yk=nullptr;
				auto &res=terms[in.result], &op1=terms[in.op1];
				read_buffer(res.r, rr, ndrSize);
				read_buffer(res.i, ri, ndrSize);
				read_buffer(res.j, rj, ndrSize);
				read_buffer(res.k, rk, ndrSize);

				read_buffer(op1.r, xr, ndrSize);
				read_buffer(op1.i, xi, ndrSize);
				read_buffer(op1.j, xj, ndrSize);
				read_buffer(op1.k, xk, ndrSize);
				if(in.is_binary())
				{
					auto &op2=terms[in.op2];
					read_buffer(op2.r, yr, ndrSize);
					read_buffer(op2.i, yi, ndrSize);
					read_buffer(op2.j, yj, ndrSize);
					read_buffer(op2.k, yk, ndrSize);
				}
				int Xplaces=(int)mp.Xplaces, Yplaces=(int)mp.Yplaces, Zplaces=(int)mp.Zplaces;
				int workSize=(int)ndrSize, offset=0;
				int x1=0, x2=Xplaces, y1=0, y2=Yplaces, z1=0, z2=Zplaces;//
				int dx=x2-x1, dy=y2-y1;
				using namespace SW;
				auto res_r=VectP(&dprr);
				auto res_c=CompP(&dprr, &dpri);
				auto res_q=QuatP(&dprr, &dpri, &dprj, &dprk);
#define			IDX		int x=x1+k%dx, y=y1+(k/dx)%dy, z=z1+k/(dy*dx), idx=offset+Xplaces*(Yplaces*z+y)+x
#define			IN1R	in1(xr+idx)
#define			IN1C	in1(xr+idx, xi+idx)
#define			IN1Q	in1(xr+idx, xi+idx, xj+idx, xk+idx)
#define			IN2R	in1(yr+idx)
#define			IN2C	in1(yr+idx, yi+idx)
#define			IN2Q	in1(yr+idx, yi+idx, yj+idx, yk+idx)
#define			OUTR	out(rr+idx)
#define			OUTC	out(rr+idx, ri+idx)
#define			OUTQ	out(rr+idx, ri+idx, rj+idx, rk+idx)
				switch(in.type)
				{
				//case 'c':
				//	{
				//		Solve_UserFunction uf(ex, in, false);
				//		for(int k=0;k<workSize;++k)
				//		{
				//			IDX;
				//			uf(idx);
				//		}
				//	}
				//	break;
				case SIG_R_R://r_r
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_r(res_r, IN1R);
						OUTR;
					//	rr[idx]=in.r_r(xr[idx]);
					}
					break;
				case SIG_C_C://c_c
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_c(res_c, IN1C);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_c(Comp1d(xr[idx], xi[idx]));
					}
					break;
				case SIG_Q_Q://q_q
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.q_q(res_q, IN1Q);
						OUTQ;
					//	QuatRef32(rr[idx], ri[idx], rj[idx], rj[idx])=
					//		in.q_q(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]));
					}
					break;
				case SIG_R_RR://r_rr
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_rr(res_r, IN1R, IN2R);
						OUTR;
					//	rr[idx]=in.r_rr(xr[idx], yr[idx]);
					}
					break;
				case SIG_C_RC://c_rc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_rc(res_c, IN1R, IN2C);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_rc(xr[idx], Comp1d(yr[idx], yi[idx]));
					}
					break;
				case SIG_Q_RQ://q_rq
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.q_rq(res_q, IN1R, IN2Q);
						OUTQ;
					//	QuatRef32(rr[idx], ri[idx], rj[idx], rk[idx])=
					//		in.q_rq(xr[idx], Quat1d(yr[idx], yi[idx], yj[idx], yk[idx]));
					}
					break;
				case SIG_C_CR://c_cr
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_cr(res_c, IN1C, IN2R);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_cr(Comp1d(xr[idx], xi[idx]), yr[idx]);
					}
					break;
				case SIG_C_CC://c_cc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_cc(res_c, IN1C, IN2C);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_cc(Comp1d(xr[idx], xi[idx]), Comp1d(yr[idx], yi[idx]));
					}
					break;
				case SIG_Q_CQ://q_cq
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.q_cq(res_q, IN1C, IN2Q);
						OUTQ;
					//	QuatRef32(rr[idx], ri[idx], rj[idx], rk[idx])=
					//		in.q_cq(Comp1d(xr[idx], xi[idx]), Quat1d(yr[idx], yi[idx], yj[idx], yk[idx]));
					}
					break;
				case SIG_Q_QR://q_qr
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.q_qr(res_q, IN1Q, IN2R);
						OUTQ;
					//	QuatRef32(rr[idx], ri[idx], rj[idx], rk[idx])=
					//		in.q_qr(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), yr[idx]);
					}
					break;
				case SIG_Q_QC://q_qc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.q_qc(res_q, IN1Q, IN2C);
						OUTQ;
					//	QuatRef32(rr[idx], ri[idx], rj[idx], rk[idx])=
					//		in.q_qc(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), Comp1d(yr[idx], yi[idx]));
					}
					break;
				case SIG_Q_QQ://q_qq
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.q_qq(res_q, IN1Q, IN2Q);
						OUTQ;
					//	QuatRef32(rr[idx], ri[idx], rj[idx], rk[idx])=
					//		in.q_qq(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), Quat1d(yr[idx], yi[idx], yj[idx], yk[idx]));
					}
					break;

				case SIG_C_R://c_r
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_r(res_c, IN1R);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_r(xr[idx]);
					}
					break;
				case SIG_C_Q://c_q
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_q(res_c, IN1Q);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_q(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]));
					}
					break;

				case SIG_R_C://r_c
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_c(res_r, IN1C);
						OUTR;
					//	rr[idx]=in.r_c(Comp1d(xr[idx], xi[idx]));
					}
					break;
				case SIG_R_Q://r_q
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_q(res_r, IN1Q);
						OUTR;
					//	rr[idx]=in.r_q(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]));
					}
					break;

				case SIG_C_RR://c_rr
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_rr(res_c, IN1R, IN2R);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_rr(xr[idx], yr[idx]);
					}
					break;

				case SIG_R_RC://r_rc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_rc(res_r, IN1R, IN2C);
						OUTR;
					//	rr[idx]=in.r_rc(xr[idx], Comp1d(yr[idx], yi[idx]));
					}
					break;
				case SIG_R_RQ://r_rq
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_rq(res_r, IN1R, IN2Q);
						OUTR;
					//	rr[idx]=in.r_rq(xr[idx], Quat1d(yr[idx], yi[idx], yj[idx], yk[idx]));
					}
					break;
				case SIG_R_CR://r_cr
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_cr(res_r, IN1C, IN2R);
						OUTR;
					//	rr[idx]=in.r_cr(Comp1d(xr[idx], xi[idx]), yr[idx]);
					}
					break;
				case SIG_R_CC://r_cc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_cc(res_r, IN1C, IN2C);
						OUTR;
					//	rr[idx]=in.r_cc(Comp1d(xr[idx], xi[idx]), Comp1d(yr[idx], yi[idx]));
					}
					break;
				case SIG_R_CQ://r_cq
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_cq(res_r, IN1C, IN2Q);
						OUTR;
					//	rr[idx]=in.r_cq(Comp1d(xr[idx], xi[idx]), Quat1d(yr[idx], yi[idx], yj[idx], yk[idx]));
					}
					break;
				case SIG_R_QR://r_qr
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_qr(res_r, IN1Q, IN2R);
						OUTR;
					//	rr[idx]=in.r_qr(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), yr[idx]);
					}
					break;
				case SIG_R_QC://r_qc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_qc(res_r, IN1Q, IN2C);
						OUTR;
					//	rr[idx]=in.r_qc(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), Comp1d(yr[idx], yi[idx]));
					}
					break;
				case SIG_R_QQ://r_qq
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.r_qq(res_r, IN1Q, IN2Q);
						OUTR;
					//	rr[idx]=in.r_qq(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), Quat1d(yr[idx], yi[idx], yj[idx], yk[idx]));
					}
					break;

				case SIG_C_QC://c_qc
					for(int k=0;k<workSize;++k)
					{
						IDX;
						in.ia32.c_qc(res_c, IN1Q, IN2C);
						OUTC;
					//	CompRef32(rr[idx], ri[idx])=in.c_qc(Quat1d(xr[idx], xi[idx], xj[idx], xk[idx]), Comp1d(yr[idx], yi[idx]));
					}
					break;
				}
			}
#endif
		}//end solve
#if 0
		if(terms.size()==3)
		{
			DebugBuffer buffers[]=
			{
				{terms[0].r, "x"},
				{terms[1].r, "y"},
				{terms[2].r, "z"},

			//	{terms[0].r, "x"},
			//	{terms[1].i, "i"},
			//	{terms[2].r, "y"},
			};
			debug_printbuffers(buffers, sizeof(buffers)/sizeof(DebugBuffer), mp.Xplaces, mp.Yplaces, 1);
		//	debug_printbuffers(buffers, sizeof(buffers)/sizeof(DebugBuffer), mp.Xplaces, mp.Yplaces, mp.Xplaces/10);
		}
#endif
		prof_add("solve");
		unsigned gl_texture=0;
		GPUBuffer *gl_buf=nullptr;
	//	unsigned gl_vertex_buf=0, gl_idx_buf=0;
		{
			auto p=(unsigned*)(&time+1);
			switch(mp.mode_idx)
			{
			case MODE_I2D:
			case MODE_C2D:
				gl_texture=p[0];
				break;
			case MODE_C3D:
				gl_buf=(GPUBuffer*)*p;
			//	gl_vertex_buf=p[1], gl_idx_buf=p[2];
				break;
			}
		}
		static cl_context context0=nullptr;
		static cl_mem image=nullptr;
		auto &result=terms[ex.resultTerm];
		switch(mp.mode_idx)
		{
		case MODE_I2D:
			//if(!p_clCreateFromGLTexture)
			//	LOGERROR("OpenCL-TI2D is not supported yet w/o clCreateFromGLTexture.");
			//else
			{
#if 0//DEBUG
				auto xr=(float*)malloc(ndrSize*sizeof(float));
				error=p_clEnqueueReadBuffer(commandqueue, result.r, CL_TRUE, 0, ndrSize*sizeof(float), xr, 0, nullptr, nullptr);	CL_CHECK(error);
				rgb=(int*)realloc(rgb, ndrSize*sizeof(int));
				int size[]={(int)dimension_work_size(mp.Xplaces), (int)dimension_work_size(mp.Yplaces), 1, (int)mp.Xplaces, (int)mp.Yplaces, 1};
				float curvecolor[]={1, 1, 1};
				for(int ky=0;ky<mp.Yplaces;++ky)
					for(int kx=0;kx<mp.Xplaces;++kx)
						ti2d_rgb(kx, ky, size, xr, curvecolor, rgb);
				free(xr);
#endif
#if 1//RELEASE
				cl_kernel kernel=nullptr;
				kernel=kernels[V_TI2D_RGB];
				bool use_texture=cl_gl_interop&&OCL_version>120&&p_clCreateFromGLTexture;
				if(use_texture)
				{
				//	if(context0!=context)
				//	{
						context0=context;
						cl_free_buffer(image);
						image=p_clCreateFromGLTexture(context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, gl_texture, &error);	CL_CHECK(error);
				//	}
				}
				else//use int buffer
				{
					image=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndrSize*sizeof(int), nullptr, &error);	CL_CHECK(error);
				}
				cl_mem color_buf=p_clCreateBuffer(context, CL_MEM_READ_ONLY, 3*sizeof(float), nullptr, &error);	CL_CHECK(error);
				float color[3];
				if(mp.nExpr>1)
				{
					int c=ex.getColor();
					color[0]=(float)(unsigned char)c*inv255, color[1]=(float)(unsigned char)(c>>8)*inv255, color[2]=(float)(unsigned char)(c>>16)*inv255;
				}
				else
					color[0]=color[1]=color[2]=0;
			//	float color[]={1, 1, 1};//white for debug purposes
				error=p_clEnqueueWriteBuffer(commandqueue, color_buf, CL_FALSE, 0, 3*sizeof(float), color, 0, nullptr, nullptr);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &result.r);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &color_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 3, sizeof(cl_mem), &image);		CL_CHECK(error);
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 2, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				if(!use_texture)
				{
					rgb=(int*)realloc(rgb, ndrSize*sizeof(int));
					error=p_clEnqueueReadBuffer(commandqueue, image, CL_TRUE, 0, ndrSize*sizeof(int), rgb, 0, nullptr, nullptr);	CL_CHECK(error);
					cl_free_buffer(image);
				}
#endif
				prof_add("raster");
			}
			break;
		case MODE_C2D:					//result -> rgb
			{
				if(result.mathSet=='c')//always true for C2D
				{
					cl_kernel kernel=nullptr;
					bool use_texture=cl_gl_interop&&OCL_version>120&&p_clCreateFromGLTexture;
					if(use_texture)
					{
						kernel=kernels[V_C2D_RGB];
						//if(context0!=context)
						//{
							context0=context;
							cl_free_buffer(image);
						//	image=p_clCreateFromGLTexture(context, CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, gl_texture, &error);	CL_CHECK(error);//
							image=p_clCreateFromGLTexture(context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, gl_texture, &error);	CL_CHECK(error);
						//}
					}
					else
					{
						kernel=kernels[V_C2D_RGB2];
						image=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndrSize*sizeof(int), nullptr, &error);	CL_CHECK(error);
					}
					error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &result.r);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &result.i);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 3, sizeof(cl_mem), &image);		CL_CHECK(error);
					error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 2, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
					if(!use_texture)
					{
						rgb=(int*)realloc(rgb, ndrSize*sizeof(int));
						error=p_clEnqueueReadBuffer(commandqueue, image, CL_FALSE, 0, ndrSize*sizeof(int), rgb, 0, nullptr, nullptr);	CL_CHECK(error);
						cl_free_buffer(image);
#ifdef DEBUG2
						LOGI("G2_CL: rgb[0] = 0x%08X", rgb[0]);
#endif
					}
				//	if(blocking)
				//		{error=p_clFinish(commandqueue);		CL_CHECK(error);}
#if 0
					error=p_clFinish(commandqueue);	CL_CHECK(error);	//DEBUG
					prof_add("clFinish");
					auto
						ndr_r=(float*)malloc(ndrSize*sizeof(float)),
						ndr_i=(float*)malloc(ndrSize*sizeof(float));
					error=p_clEnqueueReadBuffer(commandqueue, result.r, CL_TRUE, 0, ndrSize*sizeof(float), ndr_r, 0, nullptr, nullptr);	CL_CHECK(error);
					error=p_clEnqueueReadBuffer(commandqueue, result.i, CL_TRUE, 0, ndrSize*sizeof(float), ndr_i, 0, nullptr, nullptr);	CL_CHECK(error);
					prof_add("DEBUG read");
					rgb=(int*)realloc(rgb, ndrSize*sizeof(int));
					for(unsigned k=0;k<ndrSize;++k)
					{
						colorFunction_bcw(Comp1d(ndr_r[k], ndr_i[k]), rgb+k);
						//if(!(k&15))
						//	LOGE("\nXY(%g, %g)\t0x%08X", ndr_r[k], ndr_i[k], rgb[k]);
					}
					prof_add("DEBUG rgb");
					LOGERROR("G2_CL: TOPRIGHT CORNER = %f + %f i, rgb[0] = 0x%08X", ndr_r[0], ndr_i[0], rgb[0]);//DEBUG2
				//	display_texture(0, Xplaces, 0, Yplaces, rgb, Xplaces, Yplaces);
				//	glFlush();
					free(ndr_r);
					free(ndr_i);
					//free(rgb);
#endif
				}
				prof_add("rgb");
			}
			break;
		case MODE_C3D:
			{
				float th=1;//isovalue
#ifdef V6_GPU
				int Xplaces=mp.Xplaces, Yplaces=mp.Yplaces, Zplaces=mp.Zplaces;
				int XYplaces=Xplaces*Yplaces;
				int sizes_host[7]=
				{
					Xplaces, Yplaces, Zplaces, (int)ndrSize
				};
				cl_free_buffer(size_buf);
				size_buf=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, 7*sizeof(int), nullptr, &error);	CL_CHECK(error);
				error=p_clEnqueueWriteBuffer(commandqueue, size_buf, CL_FALSE, 0, 7*sizeof(int), sizes_host, 0, nullptr, nullptr);	CL_CHECK(error);
				cl_mem ndr_buf=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndrSize*5*sizeof(float), nullptr, &error);	CL_CHECK(error);
				auto kernel=kernels[V_TI3D_AV];
				error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &result.r);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &ndr_buf);	CL_CHECK(error);
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				prof_add("K1.1 av");

				kernel=kernels[V_TI3D_AV_EAST];
				error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &result.r);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &ndr_buf);	CL_CHECK(error);
				size_t temp_g=host_sizes[0], temp_l=host_sizes_local[0];
				host_sizes[0]=1, host_sizes_local[0]=1;
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				host_sizes[0]=temp_g, host_sizes_local[0]=temp_l;
				prof_add("K1.2 av yz");
				
				kernel=kernels[V_TI3D_AV_NORTH];
				error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &result.r);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &ndr_buf);	CL_CHECK(error);
				temp_g=host_sizes[1], temp_l=host_sizes_local[1];
				host_sizes[1]=1, host_sizes_local[1]=1;
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				host_sizes[1]=temp_g, host_sizes_local[1]=temp_l;
				prof_add("K1.3 av xz");
				
				kernel=kernels[V_TI3D_AV_UP];
				error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &result.r);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &ndr_buf);	CL_CHECK(error);
				temp_g=host_sizes[2], temp_l=host_sizes_local[2];
				host_sizes[2]=1, host_sizes_local[2]=1;
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				host_sizes[2]=temp_g, host_sizes_local[2]=temp_l;
				prof_add("K1.3 av xy");
				for(int kn=0;kn<nterms;++kn)//free memory
				{
					auto &term=terms[kn];
					cl_free_buffer(term.r);
					cl_free_buffer(term.i);
					cl_free_buffer(term.j);
					cl_free_buffer(term.k);
				}
				prof_add("free terms");
				error=p_clFinish(commandqueue);	CL_CHECK(error);
				prof_add("clFinish");
				
				cl_mem coeffs_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, 7*sizeof(float), nullptr, &error);		CL_CHECK(error);
				cl_mem edgeinfo_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndrSize*sizeof(ulong), nullptr, &error);	CL_CHECK(error);
				cl_mem nvert_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndrSize*sizeof(int), nullptr, &error);	CL_CHECK(error);
				cl_mem ntrgl_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndrSize*sizeof(int), nullptr, &error);	CL_CHECK(error);
				float coeffs_host[7]=
				{
					(float)mp.cx, (float)mp.mx,
					(float)mp.cy, (float)mp.my,
					(float)mp.cz, (float)mp.mz, th//threshold (isovalue)
				};
				error=p_clEnqueueWriteBuffer(commandqueue, coeffs_buf,		CL_FALSE, 0, 7*sizeof(float), coeffs_host, 0, nullptr, nullptr);	CL_CHECK(error);
				prof_add("alloc info");
				
				kernel=kernels[TI3D_CLASSIFYEDGES];
				//__kernel void ti3d_classifyedges(__constant int *size, __constant float *coeffs, __constant float *ndr, __global ulong *edgeinfo, __global int *nvert, __global int *ntrgl)
				error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);		CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &coeffs_buf);		CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &ndr_buf);		CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 3, sizeof(cl_mem), &edgeinfo_buf);	CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 4, sizeof(cl_mem), &nvert_buf);		CL_CHECK(error);
				error=p_clSetKernelArg(kernel, 5, sizeof(cl_mem), &ntrgl_buf);		CL_CHECK(error);
				error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 3, nullptr, host_sizes, host_sizes_local, 0, nullptr, nullptr);	CL_CHECK(error);
				prof_add("K2 classify edges");
				error=p_clFinish(commandqueue);	CL_CHECK(error);//just a debug check
				prof_add("clFinish");//

				std::vector<ulong> edgeinfo(ndrSize);
				std::vector<int> nvert(ndrSize), ntrgl(ndrSize);
				error=p_clEnqueueReadBuffer(commandqueue, edgeinfo_buf, CL_TRUE, 0, ndrSize*sizeof(ulong), edgeinfo.data(), 0, nullptr, nullptr);	CL_CHECK(error);
				error=p_clEnqueueReadBuffer(commandqueue, nvert_buf,	CL_TRUE, 0, ndrSize*sizeof(int), nvert.data(),	0, nullptr, nullptr);		CL_CHECK(error);//TODO: 1 buffer nvert_ntrgl_buf[kw]=ntrgl[kw]<<8|nvert[kw];
				error=p_clEnqueueReadBuffer(commandqueue, ntrgl_buf,	CL_TRUE, 0, ndrSize*sizeof(int), ntrgl.data(),	0, nullptr, nullptr);		CL_CHECK(error);
				prof_add("read");

				//CPU side code: determine work size and space indices for kernel 3 & reverse workidx 'bitindices' for CPU part 2
				int nvert_total=0, ntrgl_total=0;
				std::vector<WorkIdx> workidx, tworkidx;
			//	std::vector<ulong> bitindices;
			//	std::map<WorkIdx, int> bitindices;
				g_max_trgl_cube=0;
				for(unsigned kw=0;kw<ndrSize;++kw)
				{
					//fill 'workidx' & 'tworkidx'
					short kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
					if(kx<Xplaces-1&&ky<Yplaces-1&&kz<Zplaces-1)
					{
						int nvert_k=nvert[kw], ntrgl_k=ntrgl[kw];
						if(nvert_k)
						{
							auto mask=obsmask[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)];//boundary mask: selects original bits at boundary
							auto uwork=edgeinfo[kw]&mask;//unique work
							for(unsigned short ke=0;ke<54;++ke)//for each original bit
								if(uwork>>ke&1)
									workidx.push_back(WorkIdx(kx, ky, kz, ke));
						}
						if(ntrgl_k)
						{
							tworkidx.push_back(WorkIdx(kw, ntrgl_total));
							if((int)g_max_trgl_cube<ntrgl_k)
								g_max_trgl_cube=ntrgl_k;
						}
						nvert_total+=nvert_k, ntrgl_total+=ntrgl_k;
					}
				}
				int active_ndrSize=tworkidx.size();
#if 0//DEBUG
				//{//get extended NDR
				//	std::vector<float> ndr5(ndrSize*5);
				//	error=p_clEnqueueReadBuffer(commandqueue, ndr_buf, CL_TRUE, 0, ndrSize*5*sizeof(float), ndr5.data(), 0, nullptr, nullptr);	CL_CHECK(error);
				//	exndr2clipboard(ndr5.data(), Xplaces, Yplaces, Zplaces);//identical to V6 CPU
				//}
				//classifyedges2clipboard(edgeinfo.data(), nvert.data(), ntrgl.data(), Xplaces, Yplaces, Zplaces);//get edgeinfo, nvert, ntrgl
				{
					std::stringstream LOL_1;//get tworkidx
					LOL_1<<"V6 GPU tworkidx (active cubes)\r\n"
						"idx\tkw(kx, ky, kz)\tktr\tki=ktr*3\r\n";
					for(int k=0;k<active_ndrSize;++k)
					{
						auto &idx=tworkidx[k];
						int kx=idx.kw%Xplaces, ky=idx.kw/Xplaces%Yplaces, kz=idx.kw/XYplaces;
						LOL_1<<k<<"\t"<<idx.kw<<'('<<kx<<", "<<ky<<", "<<kz<<")\t"<<idx.kt<<'\t'<<idx.kt*3<<"\r\n";
					}
					copy_to_clipboard(LOL_1.str());
					messageboxa(ghWnd, "Information", "Active cubes copied to clipboard.");
				}
#endif
#if 0//DEBUG: get unique values from nvert
				std::stringstream LOL_1;
				std::vector<int> unique_nvert;
				for(int k=0;k<ndrSize;++k)
				{
					int nvk=nvert[k];
					bool found=false;
					for(int k2=0;k2<(int)unique_nvert.size();++k2)
					{
						if(nvk==unique_nvert[k2])
						{
							found=true;
							break;
						}
					}
					if(!found)
					{
						unique_nvert.push_back(nvk);
						sprintf_s(g_buf, g_buf_size, "%d\t0x%08X\r\n", k, nvk);
						LOL_1<<g_buf;
					}
				}
				copy_to_clipboard(LOL_1.str());
				messageboxa(ghWnd, "Information", "Unique nvert copied to clipboard.");
				//int hw=hammingweight(0xc0f00000c0f00000);//
#endif
				prof_add("CPU workidx");
				
				if(nvert_total&&ntrgl_total&&active_ndrSize)
				{
					cl_mem workidx_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, nvert_total*sizeof(ulong), nullptr, &error);		CL_CHECK(error);
					cl_mem tworkidx_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, active_ndrSize*sizeof(ulong), nullptr, &error);	CL_CHECK(error);
					cl_mem vertices_buf=nullptr, indices_buf=nullptr;
					bool use_gl_buffers=cl_gl_interop&&OCL_version>120;
					if(use_gl_buffers)
					{
						gl_buf->create_VN_CL(nvert_total*6, ntrgl_total*3);
						vertices_buf	=p_clCreateFromGLBuffer(context, CL_MEM_WRITE_ONLY, gl_buf->VBO, &error);						CL_CHECK(error);
						indices_buf		=p_clCreateFromGLBuffer(context, CL_MEM_WRITE_ONLY, gl_buf->EBO, &error);						CL_CHECK(error);
					}
					else
					{
						vertices_buf	=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, nvert_total*6*sizeof(float), nullptr, &error);	CL_CHECK(error);
						indices_buf		=p_clCreateBuffer(context, CL_MEM_WRITE_ONLY, ntrgl_total*3*sizeof(int), nullptr, &error);		CL_CHECK(error);
					}
					prof_add("alloc VN");
					//size: {Xplaces, Yplaces, Zplaces, ndrSize, nvert_total, active_ndrsize, ntrgl_total}
					sizes_host[4]=nvert_total, sizes_host[5]=active_ndrSize, sizes_host[6]=ntrgl_total;
					error=p_clEnqueueWriteBuffer(commandqueue, size_buf,		CL_FALSE, 0, 7*sizeof(int), sizes_host, 0, nullptr, nullptr);	CL_CHECK(error);
					error=p_clEnqueueWriteBuffer(commandqueue, workidx_buf,		CL_FALSE, 0, nvert_total*sizeof(ulong), workidx.data(), 0, nullptr, nullptr);		CL_CHECK(error);
					error=p_clEnqueueWriteBuffer(commandqueue, tworkidx_buf,	CL_FALSE, 0, active_ndrSize*sizeof(ulong), tworkidx.data(), 0, nullptr, nullptr);	CL_CHECK(error);
					prof_add("send workidx");
				
				//	log_start(LL_PROGRESS);//
					kernel=kernels[TI3D_ZEROCROSS];//3) ti3d_zerocross
					//__kernel void ti3d_zerocross(__constant int *size, __constant float *coeffs, __constant float *ndr, __constant ulong *workidx, __global float *vertices)
					error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &coeffs_buf);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &ndr_buf);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 3, sizeof(cl_mem), &workidx_buf);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 4, sizeof(cl_mem), &vertices_buf);	CL_CHECK(error);
					size_t l_worksize=g_maxlocalsize,
				//	size_t l_worksize=1,
						g_worksize=nvert_total-nvert_total%g_maxlocalsize+g_maxlocalsize;
					//	g_worksize=1;
					error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 1, nullptr, &g_worksize, &l_worksize, 0, nullptr, nullptr);	CL_CHECK(error);
					prof_add("K3 zerocross");

					kernel=kernels[TI3D_TRGL_INDICES];//4) ti3d_trgl_indices
					//__kernel void ti3d_trgl_indices(__constant int *size, __constant ulong *ac_workidx, __constant ulong *edgeinfo, __constant int *workidx, __global int *indices)
					error=p_clSetKernelArg(kernel, 0, sizeof(cl_mem), &size_buf);		CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 1, sizeof(cl_mem), &tworkidx_buf);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 2, sizeof(cl_mem), &edgeinfo_buf);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 3, sizeof(cl_mem), &workidx_buf);	CL_CHECK(error);
					error=p_clSetKernelArg(kernel, 4, sizeof(cl_mem), &indices_buf);	CL_CHECK(error);
					g_worksize=active_ndrSize-active_ndrSize%g_maxlocalsize+g_maxlocalsize;
				//	g_worksize=1;
					error=p_clEnqueueNDRangeKernel(commandqueue, kernel, 1, nullptr, &g_worksize, &l_worksize, 0, nullptr, nullptr);	CL_CHECK(error);
					prof_add("K4 indices");
#if 0//DEBUG
				std::vector<float> vertices(nvert_total*6);
				std::vector<int> indices(ntrgl_total*3);
				error=p_clEnqueueReadBuffer(commandqueue, vertices_buf, CL_TRUE, 0, vertices.size()*sizeof(float), vertices.data(), 0, nullptr, nullptr);	CL_CHECK(error);
				error=p_clEnqueueReadBuffer(commandqueue, indices_buf, CL_TRUE, 0, indices.size()*sizeof(int), indices.data(), 0, nullptr, nullptr);		CL_CHECK(error);
			//	cl_test();//
#if 0
				{
					error=p_clEnqueueReadBuffer(commandqueue, workidx_buf, CL_TRUE, 0, nvert_total*sizeof(ulong), workidx.data(), 0, nullptr, nullptr);	CL_CHECK(error);
					std::stringstream LOL_1;//get workidx calculated by CPU
					LOL_1<<"V6 GPU workidx\r\n";
					LOL_1<<"e/v\tkx\tky\tkz\tke\r\n";
					for(int k=0;k<(int)workidx.size();++k)
					{
						auto &idx=workidx[k];
						LOL_1<<k<<"\t("<<idx.kx<<", "<<idx.ky<<", "<<idx.kz<<", "<<idx.ke<<' '<<edgenumber2str(idx.ke)<<")\r\n";
					}
					copy_to_clipboard(LOL_1.str());
					messageboxa(ghWnd, "Information", "workidx (vertex locations) copied to clipboard.");
				}
			/*	{//
					std::stringstream LOL_1;//DEBUG test 9: get just active cubes from ti3d_trgl_indices
					LOL_1<<"V6 GPU Active cubes reported by ti3d_trgl_indices\r\n";
					LOL_1<<"id\t(kx, ky, kz)\tstarting ki\r\n";
					for(int k=0;k<(int)indices.size();++k)
					{
						int idx=indices[k];
						LOL_1<<k<<"\t("<<(idx&7)<<", "<<(idx>>3&7)<<", "<<(idx>>6&7)<<")\t"<<(idx>>9)<<"\r\n";
					}
					copy_to_clipboard(LOL_1.str());
					messageboxa(ghWnd, "Information", "Active cubes copied to clipboard.");
				}//*/
				std::vector<int> indices2(indices.size());
				for(int k=0;k<(int)indices.size();++k)
				{
					int idx2=indices[k];
					WorkIdx idx(idx2>>6&7, idx2>>9&7, idx2>>12&7, idx2&0x3F);
					int kv=0;
					int lk=0, rk=nvert_total-1;
					ulong v=0;
					for(int ik=0;ik<31;++ik)
					{
						kv=(lk+rk)>>1;
						if(kv>0&&kv<nvert_total)
							v=workidx[kv].idx;
						else
							v=workidx[0].idx;
					//	v=workidx[kv];
						lk=v<idx.idx?kv+1:lk;
						rk=v>idx.idx?kv-1:rk;
					//	lk=select(lk, kv+1, v<idx);
					//	rk=select(rk, kv-1, v>idx);
					}
					indices2[k]=kv;
				}
				{
					std::stringstream LOL_1;//DEBUG test 7 special
					LOL_1<<
						"indices\r\n"
						"idx-idx\t(kx, ky, kz, ke)\tjump(kx, ky, kz, ke)\tCPU idx\r\n";
					//LOL_1<<
					//	"indices\r\n"
					//	"idx-idx\tkx\tky\tkz\tke\tCPU idx\r\n";
					for(int k=0;k<(int)indices.size();++k)
					{
						int idx=(short)indices[k], ke=idx&0x3F,
							idx2=short(indices[k]>>16), ke2=idx2&0x3F;

						LOL_1<<k<<"\t("<<(idx >>6&7)<<", "<<(idx >>9&7)<<", "<<(idx >>12&7)<<", "<<ke <<' '<<edgenumber2str(ke )<<")\t";
						if(idx==idx2)
							LOL_1<<"same\t\t";
						else
							LOL_1<<"("<<(idx2>>6&7)<<", "<<(idx2>>9&7)<<", "<<(idx2>>12&7)<<", "<<ke2<<' '<<edgenumber2str(ke2)<<")\t";
						LOL_1<<indices2[k]<<"\r\n";

					//	const char *comparison=idx==idx2?"same":"";
					//	LOL_1<<k<<"\t("
					//				<<(idx >>6&7)<<", "<<(idx >>9&7)<<", "<<(idx >>12&7)<<", "<<ke <<' '<<edgenumber2str(ke )<<")\t"<<comparison
					//		<<"("	<<(idx2>>6&7)<<", "<<(idx2>>9&7)<<", "<<(idx2>>12&7)<<", "<<ke2<<' '<<edgenumber2str(ke2)<<")\t"
					//		<<indices2[k]<<"\r\n";

					//	LOL_1<<k<<"\t("<<(idx>>6&7)<<", "<<(idx>>9&7)<<", "<<(idx>>12&7)<<", "<<ke<<' '<<edgenumber2str(ke)<<")\t"<<indices2[k]<<"\r\n";
					//	LOL_1<<k<<'\t'<<(idx>>6&7)<<'\t'<<(idx>>9&7)<<'\t'<<(idx>>12&7)<<'\t'<<(idx&0x3F)<<":\t"<<indices2[k]<<"\r\n";
					//	LOL_1<<k<<'\t'<<(idx>>6&7)<<'\t'<<(idx>>9&7)<<'\t'<<(idx>>12&7)<<'\t'<<(idx&0x3F)<<"\r\n";
					}
					copy_to_clipboard(LOL_1.str());
					messageboxa(ghWnd, "Information", "indices + info copied to clipboard.");
				}//*/
#endif
#endif
					if(!use_gl_buffers)
					{
						std::vector<float> vertices(nvert_total*6);
						std::vector<int> indices(ntrgl_total*3);
						error=p_clEnqueueReadBuffer(commandqueue, vertices_buf,	CL_TRUE, 0, nvert_total*6*sizeof(float), vertices.data(), 0, nullptr, nullptr);	CL_CHECK(error);
						error=p_clEnqueueReadBuffer(commandqueue, indices_buf,	CL_TRUE, 0, ntrgl_total*3*sizeof(int), indices.data(), 0, nullptr, nullptr);	CL_CHECK(error);
						gl_buf->create_VN_I(vertices.data(), nvert_total*6, indices.data(), ntrgl_total*3);
					}
					cl_free_buffer(workidx_buf), cl_free_buffer(tworkidx_buf);
					cl_free_buffer(vertices_buf);
					cl_free_buffer(indices_buf);
				}
				g_nvert=nvert_total, g_ntrgl=ntrgl_total;

				cl_free_buffer(ndr_buf);
				cl_free_buffer(coeffs_buf), cl_free_buffer(edgeinfo_buf), cl_free_buffer(nvert_buf), cl_free_buffer(ntrgl_buf);
				prof_add("free");

				error=p_clFinish(commandqueue);	CL_CHECK(error);
				prof_add("clFinish");
#endif
#if defined V5_CPU||defined V6_CPU
				int Xplaces=mp.Xplaces, Yplaces=mp.Yplaces, Zplaces=mp.Zplaces;
				int XYplaces=Xplaces*Yplaces;
				auto ndr=(float*)malloc(ndrSize*5*sizeof(float));
			//	memset(ndr, 0, ndrSize*5*sizeof(float));//
				memset(ndr, 0xCD, ndrSize*5*sizeof(float));//
#ifdef DEBUG3
				debug_printbuffer(result.r, Xplaces, Yplaces, Zplaces, 1, "result");
				messageboxa(ghWnd, "Information", "Result sent to clipboard");
#endif
				error=p_clEnqueueReadBuffer(commandqueue, result.r, CL_TRUE, 0, ndrSize*sizeof(float), ndr, 0, nullptr, nullptr);	CL_CHECK(error);
				prof_add("read");
				for(int kz=0, zend=Zplaces-1;kz<zend;++kz)//kernel 1.1: compute average (interpolated) values
				{
					for(int ky=0, yend=Yplaces-1;ky<yend;++ky)
					{
						for(int kx=0, xend=Xplaces-1;kx<xend;++kx)//for each data cube
						{
							int idx=Xplaces*(Yplaces*kz+ky)+kx;
							float
								V_UNW=ndr[idx+XYplaces+Xplaces],	V_UNE=ndr[idx+XYplaces+Xplaces+1],
								V_USW=ndr[idx+XYplaces],			V_USE=ndr[idx+XYplaces+1],

								V_DNW=ndr[idx+Xplaces],				V_DNE=ndr[idx+Xplaces+1],
								V_DSW=ndr[idx],						V_DSE=ndr[idx+1];
							float
								V_mmW=(V_UNW+V_DNW+V_DSW+V_USW)*0.25f,
								V_mSm=(V_DSW+V_USW+V_USE+V_DSE)*0.25f,
								V_Dmm=(V_DNW+V_DSW+V_DNE+V_DSE)*0.25f,
								V_mmm=(V_UNW+V_DNW+V_DSW+V_USW+V_UNE+V_DNE+V_USE+V_DSE)*0.125f;
							ndr[ndrSize+(idx<<2)  ]=V_Dmm;
							ndr[ndrSize+(idx<<2)+1]=V_mSm;
							ndr[ndrSize+(idx<<2)+2]=V_mmW;
							ndr[ndrSize+(idx<<2)+3]=V_mmm;
						}
					}
				}
				prof_add("k1 av1");
				for(int kz=0, zend=Zplaces-1;kz<zend;++kz)//kernel 1.2: compute average (interpolated) values on east NDR face
				{
					for(int ky=0, yend=Yplaces-1;ky<yend;++ky)
					{
						int kx=Xplaces-1;
					//	for(int kx=0, xend=Xplaces-1;kx<xend;++kx)
						{
							int idx=Xplaces*(Yplaces*kz+ky)+kx;
							float
								V_UNW=ndr[idx+XYplaces+Xplaces],
								V_USW=ndr[idx+XYplaces],

								V_DNW=ndr[idx+Xplaces],
								V_DSW=ndr[idx];
							float
								V_mmW=(V_UNW+V_DNW+V_DSW+V_USW)*0.25f;
							ndr[ndrSize+(idx<<2)+2]=V_mmW;
						}
					}
				}
				prof_add("k1 av2");
				for(int kz=0, zend=Zplaces-1;kz<zend;++kz)//kernel 1.3: compute average (interpolated) values on north NDR face
				{
					int ky=Yplaces-1;
				//	for(int ky=0, yend=Yplaces-1;ky<yend;++ky)
					{
						for(int kx=0, xend=Xplaces-1;kx<xend;++kx)
						{
							int idx=Xplaces*(Yplaces*kz+ky)+kx;
							float
								V_USW=ndr[idx+XYplaces],			V_USE=ndr[idx+XYplaces+1],

								V_DSW=ndr[idx],						V_DSE=ndr[idx+1];
							float
								V_mSm=(V_DSW+V_USW+V_USE+V_DSE)*0.25f;
							ndr[ndrSize+(idx<<2)+1]=V_mSm;
						}
					}
				}
				prof_add("k1 av3");
			//	for(int kz=0, zend=Zplaces-1;kz<zend;++kz)//kernel 1.4: compute average (interpolated) values on up NDR face
				{
					int kz=Zplaces-1;
					for(int ky=0, yend=Yplaces-1;ky<yend;++ky)
					{
						for(int kx=0, xend=Xplaces-1;kx<xend;++kx)
						{
							int idx=Xplaces*(Yplaces*kz+ky)+kx;
							float
								V_DNW=ndr[idx+Xplaces],				V_DNE=ndr[idx+Xplaces+1],
								V_DSW=ndr[idx],						V_DSE=ndr[idx+1];
							float
								V_Dmm=(V_DNW+V_DSW+V_DNE+V_DSE)*0.25f;
							ndr[ndrSize+(idx<<2)  ]=V_Dmm;
						}
					}
				}
				prof_add("k1 av4");
				auto edgeinfo=(ulong*)malloc(ndrSize*sizeof(ulong));
				memset(edgeinfo, 0, ndrSize*sizeof(ulong));
				auto nvert=(int*)malloc(ndrSize*sizeof(int));//number of edges producing a vertex
				memset(nvert, 0, ndrSize*sizeof(int));
				auto ntrgl=(int*)malloc(ndrSize*sizeof(int));
				memset(ntrgl, 0, ndrSize*sizeof(int));
				//kernel 2: classify edges: fills edgeinfo, nvert, ntrgl
				//arguments: __global const float *ndr, __global ulong *edgeinfo, __global char *nvert, __global char *ntrgl
				for(int kz=0, zend=Zplaces-1;kz<zend;++kz)
				{
					for(int ky=0, yend=Yplaces-1;ky<yend;++ky)
					{
						for(int kx=0, xend=Xplaces-1;kx<xend;++kx)//for each data cube
						{
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
#define						EC(number)	((ulong)ec[number]<<number)
							edgeinfo[idx]=
								EC( 0)|EC( 1)|EC( 2)|EC( 3)|EC( 4)|EC( 5)|EC( 6)|EC( 7)|EC( 8)|
								EC( 9)|EC(10)|EC(11)|EC(12)|EC(13)|EC(14)|EC(15)|EC(16)|EC(17)|
								EC(18)|EC(19)|EC(20)|EC(21)|EC(22)|EC(23)|EC(24)|EC(25)|EC(26)|
								EC(27)|EC(28)|EC(29)|EC(30)|EC(31)|EC(32)|EC(33)|EC(34)|EC(35)|
								EC(36)|EC(37)|EC(38)|EC(39)|EC(40)|EC(41)|EC(42)|EC(43)|EC(44)|
								EC(45)|EC(46)|EC(47)|EC(48)|EC(49)|EC(50)|EC(51)|EC(52)|EC(53);
#undef						EC
							if(kx==3&&ky==0&&kz==0)
								int LOL_1=0;
							auto mask=obsmask[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)];
							nvert[idx]=hammingweight(edgeinfo[idx]&mask);
#define						NT(V0, V1, V2, V3)		nt[V3<<3|V2<<2|V1<<1|V0]
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
#undef						NT
						}
					}
				}
				prof_add("k2 condition");
				//blocking read: edgeinfo, nvert, ntrgl
#if 0
				std::stringstream LOL_1;
				LOL_1<<"edgeinfo,\t\tkw(kx, ky, kz),\tnvert,\tntrgl\r\n";
				for(int kw=0;kw<ndrSize;++kw)
				{
					if(edgeinfo[kw])
					{
						auto p=(unsigned short*)(edgeinfo+kw);
						short kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
						sprintf_s(g_buf, g_buf_size, "0x%04X %04X %04X %04X, %d(%d, %d, %d),\t%d,\t%d\r\n", (int)p[3], (int)p[2], (int)p[1], (int)p[0], kw, kx, ky, kz, (int)nvert[kw], (int)ntrgl[kw]);
						LOL_1<<g_buf;
					}
				}
				copy_to_clipboard(LOL_1.str());
#endif
			//	std::stringstream LOL_1;
				//CPU side code 1: determine work size and space indices for kernel 3 & reverse workidx 'bitindices' for CPU part 2
				int nvert_total=0, ntrgl_total=0;
				std::vector<WorkIdx> workidx, tworkidx;
#ifdef V6_CPU
				std::vector<ulong> bitindices;//same as workidx
#else
				std::map<WorkIdx, int> bitindices;
#endif
			//	std::map<long long, int> bitindices;
			//	std::vector<int> bitindices(ndrSize*33);
				for(unsigned kw=0;kw<ndrSize;++kw)
				{
					//fill 'workidx' & 'tworkidx'
					if(nvert[kw])
					{
						short kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
						auto mask=obsmask[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)];//boundary mask: selects original bits at boundary
						auto uwork=edgeinfo[kw]&mask;//unique work
						for(unsigned short ke=0;ke<54;++ke)//for each original bit
						{
							if(uwork>>ke&1)
							{
								WorkIdx idx(kx, ky, kz, ke);
#ifdef V6_CPU
								bitindices.push_back((ulong&)idx);//same as workidx
#else
								bitindices[idx]=workidx.size();
#endif
								workidx.push_back(idx);
							}
						}
					}
					if(ntrgl[kw])
						tworkidx.push_back(WorkIdx(kw, ntrgl_total));
					nvert_total+=nvert[kw], ntrgl_total+=ntrgl[kw];
				}
				prof_add("CPU1 workidx");
			//	exndr2clipboard(ndr, Xplaces, Yplaces, Zplaces);//
				classifyedges2clipboard(edgeinfo, nvert, ntrgl, Xplaces, Yplaces, Zplaces);//

				debug_info.clear();//
				auto vertices=(float*)malloc(nvert_total*6*sizeof(float));
				memset(vertices, 0, nvert_total*6*sizeof(float));
				float
					Xstart=(float)mp.cx, Xsample=(float)mp.mx,
					Ystart=(float)mp.cy, Ysample=(float)mp.my,
					Zstart=(float)mp.cz, Zsample=(float)mp.mz;
				//kernel 3: zero cross (per-edge trilinear interpolation)
				//arguments: __constant int *size, __constant float *coeffs, __global const float *ndr, __constant ulong *workidx, __global float *vertices
				//int breakpoint=7;//nvert_total>>4;
				for(int id=0;id<nvert_total;++id)
				{
					//if(id==breakpoint)
					//	int LOL_1=0;
					auto &wi=workidx[id];
					int kx=wi.kx, ky=wi.ky, kz=wi.kz, ke=wi.ke;
				//	int id=get_global_id(0);
				//	ulong winfo=workidx[id];
				//	int kx=(unsigned short)winfo, ky=(unsigned short)(winfo>>16), kz=(unsigned short)(winfo>>32), ke=(unsigned short)(winfo>>48);
					int idx=Xplaces*(Yplaces*kz+ky)+kx;
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
						V_mmE=ndr[ndrSize+((idx+1)			<<2)+2],//x+y+z, id=7 (6, 4, 0, 26), x+ face 0xCDCDCDCD
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
					//float LOL_1=-431602080.;//0xCDCDCDCD
					//int LOL_2=(int&)LOL_1;
					//a=-0, b=2.39067e-010, c=-0.000199733, d=0.000196106;
					if(a==0&&b==0&&c==0)//degenerate equation d=0
						ret.set((th-A.w)/(B.w-A.w), 1);
					else
						ret=solvecubic(a, b, c, d);
					float3 coords;
					if(ret.y)
					{
						float3 pa={A.x, A.y, A.z}, pb={B.x, B.y, B.z};
						coords=mix(pa, pb, ret.x);
						vstore3(coords, id*6, vertices);
					}
					else//solvecubic failed
					{
						ret.set((th-A.w)/(B.w-A.w), 1);
						float3 pa={A.x, A.y, A.z}, pb={B.x, B.y, B.z};
						coords=mix(pa, pb, ret.x);
						vstore3(coords, id*6, vertices);

						debug_info.push_back(DebugInfo(kx, ky, kz, ke, id, ret.x, coords.x, coords.y, coords.z));
					//	LOL_1<<'('<<kx<<", "<<ky<<", "<<kz<<", "<<ke<<") "<<a<<"xxx+"<<b<<"xx+"<<c<<"x+"<<d<<"\r\n";
						if(LOGERROR("solvecubic(%f, %f, %f, %f) failed at (%d, %d, %d, %d).", a, b, c, d, kx, ky, kz, ke))//
						{
							std::stringstream LOL_1;
							LOL_1<<a<<"xxx+"<<b<<"xx+"<<c<<"x+"<<d;
							copy_to_clipboard(LOL_1.str());
						}
					//	LOGERROR("solvecubic(%f, %f, %f, %f) failed.", a, b, c, d);//
					}
					//V_DSW_DSE = V_DSE-V_DSW
					//V_DNW_DNE = V_DNE-V_DNW
					//V_USW_USE = V_USE-V_USW
					//V_UNW_UNE = V_UNE-V_UNW
					//V_DS = V_DSW+V_DSW_DSE*ax
					//V_DN = V_DNW+V_DNW_DNE*ax
					//V_US = V_USW+V_USW_USE*ax
					//V_UN = V_UNW+V_UNW_UNE*ax
					//TY = [V_DSW_DSE+(V_DNW_DNE-V_DSW_DSE)*ay]
					//<
					//	TY+([V_USW_USE+(V_UNW_UNE-V_USW_USE)*ay]-TY)*az,
					//	[V_DN-V_DS]+([V_UN-V_US]-[V_DN-V_DS])*az,
					//	[V_US+(V_UN-V_US)*ay]-[V_DS+(V_DN-V_DS)*ay]
					//>
					float3 normal;
					{
						float3 a;
						a.set((coords.x-X_W)/Xsample, (coords.y-Y_S)/Ysample, (coords.z-Z_D)/Zsample);
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
						normal.set
						(
							TY+((V_USW_USE+(V_UNW_UNE-V_USW_USE)*a.y)-TY)*a.z,
							V_D_SN+(V_U_SN-V_D_SN)*a.z,
							(V_US+V_U_SN*a.y)-(V_DS+V_D_SN*a.y)
						);
						float inv_abs=sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
						if(inv_abs)
						{
							inv_abs=1/inv_abs;
							normal.x*=inv_abs, normal.y*=inv_abs, normal.z*=inv_abs;
						}
					}
					vstore3(normal, id*6+3, vertices);
				}
				//auto &str=LOL_1.str();
				//if(str.size())
				//{
				//	copy_to_clipboard(str);
				//	messageboxa(ghWnd, "Information", "Coordinates & equations copied to clipboard.");
				//}
				prof_add("k3 vertices");

				//kernel 4: fill 'indices'
				std::vector<int> indices(ntrgl_total*3, 0);
				for(unsigned ka=0;ka<(int)tworkidx.size();++ka)//for each active data cube
				{
					auto &twi=tworkidx[ka];//triangle work index
					unsigned kw=twi.kw, ki=twi.kt*3;
					unsigned short kx=kw%Xplaces, ky=kw/Xplaces%Yplaces, kz=kw/XYplaces;
					auto work=edgeinfo[kw];
					unsigned ktr=0;//number of triangles in this tetrahedron
					auto mask=obsmask[(kz>=Zplaces-2)<<2|(ky>=Yplaces-2)<<1|(kx>=Xplaces-2)];
					for(int kth=0;kth<28;++kth)//for each tetrahedron kth of the 28 tetrahedra in the data cube
					{
						//if(kx==6&&ky==4&&kz==0&&kth==1)//
						//	int LOL_1=0;//
						int kt6=kth*6;
						char e[]=//edge states
						{
							getbit(work, ti[kt6]),
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
								char ke=ti[kth*6+kb];//relative edge index
								WorkIdx idx;
								char original=getbit(mask, ke);
								if(original)//original bit, in this data point		note: all bits are original at boundary
									idx.set(kx, ky, kz, ke);
								else//redundant bit, look up in another data point (jump)
								{
									auto &weo=we_offsets[ke-33];
									short
										xbc=kx<Xplaces-2, ybc=ky<Yplaces-2, zbc=kz<Zplaces-2,//boundary conditions: *places-2 is last cube (pre-last datapoint)
										kx2=kx+(weo.dx&-xbc),
										ky2=ky+(weo.dy&-ybc),
										kz2=kz+(weo.dz&-zbc),
										ke2=weo.ke2[zbc<<2|ybc<<1|xbc];
									idx.set(kx2, ky2, kz2, ke2);
								}
#ifdef V6_CPU
								int kv=0;
								char found=0;

								for(int lk=0, rk=nvert_total-1, ik=0;ik<31;++ik)
								{
									kv=(lk+rk)>>1;
									ulong v=workidx[kv].idx;
									lk=v<idx.idx?kv+1:lk;
									rk=v>idx.idx?kv-1:rk;
									//lk=select(lk, kv+1, v<idx);
									//rk=select(rk, kv-1, v>idx);
								}
								found=kv>=0&&kv<nvert_total;

							//	for(int lk=0, rk=bitindices.size()-1;lk<=rk;)
							//	{
							//		kv=(lk+rk)>>1;
							//		ulong v=bitindices[kv];
							//		if(v<idx.idx)
							//			lk=kv+1;
							//		else if(v>idx.idx)
							//			rk=kv-1;
							//		else
							//		{
							//			found=1;
							//			break;
							//		}
							//	}
								if(found)
									tcvi[esum]=kv;
#else
								auto it=bitindices.find(idx);
								if(it!=bitindices.end())
									tcvi[esum]=it->second;
#endif
								else
									LOGERROR("(%d, %d, %d, %d) not found.", idx.kx, idx.ky, idx.kz, idx.ke);
								++esum;
							}
						}
						if(esum==4)//there are edges producing vertices
						{
							if(ki+5>=indices.size())
							{
								LOGERROR("indices size=%d, ki=%d, want to add 6 indices.", indices.size(), ki);
								break;
							}
							//if(ktr==8||ktr==9)//
							//	int LOL_1=0;//
							++ktr;//
							indices[ki]=tcvi[0], ++ki;
							indices[ki]=tcvi[1], ++ki;
							indices[ki]=tcvi[2], ++ki;

							//if(ktr==8||ktr==9)//
							//	int LOL_1=0;//
							++ktr;//
							indices[ki]=tcvi[1], ++ki;
							indices[ki]=tcvi[2], ++ki;
							indices[ki]=tcvi[3], ++ki;
						}
						else if(esum==3)
						{
							if(ki+2>=indices.size())
							{
								LOGERROR("indices size=%d, ki=%d, want to add 3 indices.", indices.size(), ki);
							//	LOGERROR("indices[%d+2] is OOB", ki);
								break;
							}
							//if(ktr==8||ktr==9)//
							//	int LOL_1=0;//
							++ktr;//
							indices[ki]=tcvi[0], ++ki;
							indices[ki]=tcvi[1], ++ki;
							indices[ki]=tcvi[2], ++ki;
						}
						else if(esum)
							LOGERROR("Tetrahedron produced %d vertices", esum);
						//	messageboxa(ghWnd, "Error", "Tetrahedron produced %d vertices", esum);
						//char emask=getbit(work, ti[kt6+5])<<5|getbit(work, ti[kt6+4])<<4|getbit(work, ti[kt6+3])<<3|getbit(work, ti[kt6+2])<<2|getbit(work, ti[kt6+1])<<1|getbit(work, ti[kt6]);//6 edge codition bits from the tetrahedron
						//if(emask)//there are edges producing vertices
						//	for(int kb=0;kb<6;++kb)			//X no triangle fans, just GL_TRIANGLES
						//		if(emask>>kb&1)
						//			indices[ki]=offset+bi[ti[kth*6+kb]], ++ki;//original bit indices
					}//end tetrahedron loop
					if(g_max_trgl_cube<ktr)
						g_max_trgl_cube=ktr;
				}//end work loop
				prof_add("k4 indices");
#if 0//DEBUG - results to clipboard
				std::stringstream LOL_1;
				LOL_1<<
					"vertices:\r\n"
					"vertex#\tfloat#:\tx\ty\tz\tnx\tny\tnz\r\n";
					//"kv\tidx:\tx\ty\tz\tnx\tny\tnz\r\n";
				for(int kv=0, nf=nvert_total*6;kv+5<nf;kv+=6)
					LOL_1<<kv/6<<'\t'<<kv<<":\t\t"<<vertices[kv]<<",\t"<<vertices[kv+1]<<",\t"<<vertices[kv+2]<<",\t"<<vertices[kv+3]<<",\t"<<vertices[kv+4]<<",\t"<<vertices[kv+5]<<"\r\n";
				LOL_1<<
					"\r\n"
					"indices (triangles):\r\n"
					"kt\tidx:\tv0\tv1\tv2\t\r\n";
				for(int ki=0;ki+2<(int)indices.size();ki+=3)
					LOL_1<<ki/3<<'\t'<<ki<<":\t"<<indices[ki]<<'\t'<<indices[ki+1]<<'\t'<<indices[ki+2]<<"\r\n";
				copy_to_clipboard(LOL_1.str());
				messageboxa(ghWnd, "Information", "Vertices & indices copied to clipboard.");
#endif
#if 0//DEBUG
				//indices.resize(nvert_total-nvert_total%3);//just show vertices
				//for(unsigned k=0;k<indices.size();++k)
				//	indices[k]=k;

				static int triangle_counter=0;
				++triangle_counter;

				//int idx_idx=triangle_counter*3;				//one at a time
				//if(idx_idx+3>(int)indices.size())
				//	idx_idx=indices.size()-3;
				//std::vector<int> indices2(indices.begin()+idx_idx, indices.begin()+idx_idx+3);
				//indices=std::move(indices2);

				int ntriangles=triangle_counter,//2 3 100		//first n triangles
					nvertices=3*ntriangles;
				if((int)indices.size()>nvertices)
					indices.resize(nvertices);


				//float vert2[]=//DEBUG
				//{
				//	0, 0, 0,	0, 0, 1,
				//	1, 0, 0,	0, 0, 1,
				//	0, 1, 0,	0, 0, 1,
				//};
				//nvert_total=3;
				//vertices=(float*)realloc(vertices, nvert_total*6*sizeof(float));
				//memcpy(vertices, vert2, nvert_total*6*sizeof(float));
				//indices.clear();
				//indices.push_back(0);
				//indices.push_back(1);
				//indices.push_back(2);
#endif
#if 1//DEBUG - expose vertices & indices
				debug_vertices.assign(vertices, vertices+nvert_total*6);
				debug_indices=indices;
#endif
				g_nvert=nvert_total, g_ntrgl=indices.size()/3;
				gl_buf->create_VN_I(vertices, nvert_total*6, indices.data(), indices.size());
				prof_add("send");

				free(edgeinfo), free(nvert), free(ntrgl), free(vertices);
				free(ndr);
				prof_add("free");
				prof_sum("CPU V6", 11);
#endif
			}
			break;
		}
#if 0
		{
			DebugBuffer buffers[]=
			{
				{image, "image"},
			};
			debug_printbuffers(buffers, sizeof(buffers)/sizeof(DebugBuffer), ndrSize, 512);
		}
#endif
		for(int kn=0;kn<nterms;++kn)//free memory
		{
			auto &term=terms[kn];
			cl_free_buffer(term.r);
			cl_free_buffer(term.i);
			cl_free_buffer(term.j);
			cl_free_buffer(term.k);
		}
		cl_free_buffer(size_buf);
		cl_free_buffer(args_buf);
		prof_add("free");
#if 0//DEBUG
		error=p_clFinish(commandqueue);	CL_CHECK(error);//DEBUG
		prof_add("clFinish");
		rgb=(int*)realloc(rgb, ndrSize*sizeof(int));
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, gl_texture);
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb);
#endif
	}
}
void			cl_finish()
{
	if(OCL_state<CL_READY)
		return;
	//	LOGERROR("cl_finish: OCL_state = %d", OCL_state);
//	if(OCL_API_not_loaded)
//		LOGERROR("OpenCL API not loaded.");
	else
	{
		int error=p_clFinish(commandqueue);
		CL_CHECK(error);
		prof_add("clFinish");
	}
}