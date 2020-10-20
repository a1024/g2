//best viewed with tab size of 4 spaces
//g2_expr.h - Include for G2 Expression class and its dependencies.
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
#ifndef				G2_EXPR_H
#define				G2_EXPR_H
#include			"g2_common.h"
#include			"g2_error.h"
#include			"g2_mp.h"
#include			"g2_graphics.h"//usingOpenGL

	#define	PROFILER
//	#define PROFILER_CYCLES

	#define TIMING_USE_QueryPerformanceCounter
//	#define	TIMING_USE_rdtsc
//	#define TIMING_USE_GetProcessTimes	//~16ms resolution
//	#define TIMING_USE_GetTickCount		//~16ms resolution
//	#define TIMING_USE_timeGetTime		//~16ms resolution

	#define	MULTIPRECISION

extern bool			showBenchmark;
union				G2Version
{
	struct
	{
		short cl_comp, g2;
	};
	int id;
	G2Version(short g2, short cl_component):g2(g2), cl_comp(cl_component){}
};
extern const G2Version g2_version;//hi: g2 version, lo: cl component version

//time-measuring functions
inline double		time_sec()
{
#ifdef TIMING_USE_QueryPerformanceCounter
	static long long t=0;
	static LARGE_INTEGER li={};
	QueryPerformanceFrequency(&li);
	t=li.QuadPart;
	QueryPerformanceCounter(&li);
	return (double)li.QuadPart/t;
#elif defined TIMING_USE_rdtsc
	static LARGE_INTEGER li={};
	QueryPerformanceFrequency(&li);
	return (double)__rdtsc()*0.001/li.QuadPart;//pre-multiplied by 1000
#elif defined TIMING_USE_GetProcessTimes
	FILETIME create, exit, kernel, user;
	int success=GetProcessTimes(GetCurrentProcess(), &create, &exit, &kernel, &user);
	if(success)
//#ifdef PROFILER_CYCLES
	{
		const auto hns2sec=100e-9;
		return hns2ms*(unsigned long long&)user;
	//	return hns2ms*(unsigned long long&)kernel;
	}
//#else
//	{
//		SYSTEMTIME t;
//		success=FileTimeToSystemTime(&user, &t);
//		if(success)
//			return t.wHour*3600000.+t.wMinute*60000.+t.wSecond*1000.+t.wMilliseconds;
//		//	return t.wHour*3600.+t.wMinute*60.+t.wSecond+t.wMilliseconds*0.001;
//	}
//#endif
	SYS_CHECK();
	return -1;
#elif defined TIMING_USE_GetTickCount
	return (double)GetTickCount()*0.001;//the number of milliseconds that have elapsed since the system was started
#elif defined TIMING_USE_timeGetTime
	return (double)timeGetTime()*0.001;//system time, in milliseconds
#endif
}
inline double		time_ms()
{
#ifdef TIMING_USE_QueryPerformanceCounter
	static long long t=0;
	static LARGE_INTEGER li={};
	QueryPerformanceFrequency(&li);
	t=li.QuadPart;
	QueryPerformanceCounter(&li);
	return (double)li.QuadPart*1000./t;
#elif defined TIMING_USE_rdtsc
	static LARGE_INTEGER li={};
	QueryPerformanceFrequency(&li);
	return (double)__rdtsc()/li.QuadPart;//pre-multiplied by 1000
#elif defined TIMING_USE_GetProcessTimes
	FILETIME create, exit, kernel, user;
	int success=GetProcessTimes(GetCurrentProcess(), &create, &exit, &kernel, &user);
	if(success)
//#ifdef PROFILER_CYCLES
	{
		const auto hns2ms=100e-9*1000.;
		return hns2ms*(unsigned long long&)user;
	//	return hns2ms*(unsigned long long&)kernel;
	}
//#else
//	{
//		SYSTEMTIME t;
//		success=FileTimeToSystemTime(&user, &t);
//		if(success)
//			return t.wHour*3600000.+t.wMinute*60000.+t.wSecond*1000.+t.wMilliseconds;
//		//	return t.wHour*3600.+t.wMinute*60.+t.wSecond+t.wMilliseconds*0.001;
//	}
//#endif
	SYS_CHECK();
	return -1;
#elif defined TIMING_USE_GetTickCount
	return (double)GetTickCount();//the number of milliseconds that have elapsed since the system was started
#elif defined TIMING_USE_timeGetTime
	return (double)timeGetTime();//system time, in milliseconds
#endif
}
inline double		elapsed_ms()//since last call
{
	static double t1=0;
	double t2=time_ms(), diff=t2-t1;
	t1=t2;
	return diff;
}
inline double		elapsed_cycles()//since last call
{
	static long long t1=0;
	long long t2=__rdtsc();
	double diff=double(t2-t1);
	t1=t2;
	return diff;
}
#ifdef TIMING_USE_QueryPerformanceCounter
	#undef TIMING_USE_QueryPerformanceCounter
#endif
#ifdef TIMING_USE_GetTickCount
	#undef TIMING_USE_GetTickCount
#endif
#ifdef TIMING_USE_timeGetTime
	#undef TIMING_USE_timeGetTime
#endif

//profiler
typedef std::pair<std::string, double> ProfInfo;
void				prof_start();
void				prof_add(const char *label, int divisor=1);
void				prof_sum(const char *label, int count);
void				prof_loop_start(const char **labels, int n);
void				prof_add_loop(int idx);
void				prof_print();

namespace			G2
{
	enum			Map
	{
		M_IGNORED,
		
		M_N,
		
		M_LPR, M_RPR,																//(		)
		M_COMMA,																	//,
		M_QUESTION_MARK, M_COLON,													//?		:

			M_PROCEDURAL_START,
		M_IF, M_ELSE, M_FOR, M_DO, M_WHILE,											//if else for do while
		M_CONTINUE, M_BREAK, M_RETURN,												//continue break return
		M_LBRACE, M_RBRACE,															//{		}
		M_SEMICOLON,																//;

			M_PROCEDURAL_ASSIGN_START,
		M_ASSIGN, M_ASSIGN_MULTIPLY, M_ASSIGN_DIVIDE, M_ASSIGN_MOD,					//= *= /= %=
		M_ASSIGN_PLUS, M_ASSIGN_MINUS, M_ASSIGN_LEFT, M_ASSIGN_RIGHT,				//+= -= <<= >>=
		M_ASSIGN_AND, M_ASSIGN_XOR, M_ASSIGN_OR,									//&= #= |=
			M_PROCEDURAL_ASSIGN_END,

			M_PROCEDURAL_END,

		M_INCREMENT, M_DECREMENT,													//++ --
		M_FACTORIAL_LOGIC_NOT,														//!
		M_MODULO_PERCENT,															//%
		M_BITWISE_NOT,																//~
		M_PENTATE,																	//^^^
		M_TETRATE,																	//^^
		M_POWER, M_POWER_REAL,														//^		**
		M_MULTIPLY, M_DIVIDE, M_LOGIC_DIVIDES,										//*		/		@
		M_PLUS, M_MINUS,															//+		-
		M_BITWISE_SHIFT_LEFT, M_BITWISE_SHIFT_RIGHT,								//<<	>>
		M_LOGIC_LESS, M_LOGIC_LESS_EQUAL, M_LOGIC_GREATER, M_LOGIC_GREATER_EQUAL,	//<		<=		>		>=
		M_LOGIC_EQUAL, M_LOGIC_NOT_EQUAL,											//==	!=
		M_BITWISE_AND, M_BITWISE_NAND,												//&		~& |~
		M_BITWISE_XOR, M_BITWISE_XNOR,												//#		~# #~
		M_VERTICAL_BAR, M_BITWISE_NOR,												//|		~| &~
		M_LOGIC_AND,																//&&
		M_LOGIC_XOR,																//##
		M_LOGIC_OR,																	//||
		M_CONDITION_ZERO,															//??

		M_S_EQUAL_ASSIGN, M_S_NOT_EQUAL,											//= _!=				not supported
		M_S_LESS, M_S_LESS_EQUAL, M_S_GREATER, M_S_GREATER_EQUAL,					//= _< =< _> =>		not supported

			M_FSTART,

		M_COS, M_ACOS, M_COSH, M_ACOSH, M_COSC,
		M_SEC, M_ASEC, M_SECH, M_ASECH,
		M_SIN, M_ASIN, M_SINH, M_ASINH, M_SINC, M_SINHC,
		M_CSC, M_ACSC, M_CSCH, M_ACSCH,
		M_TAN,         M_TANH, M_ATANH, M_TANC,
		M_COT, M_ACOT, M_COTH, M_ACOTH,
		M_EXP, M_LN, M_SQRT, M_CBRT, M_INVSQRT, M_SQ,
		M_GAUSS, M_ERF, M_FIB, M_ZETA, M_LNGAMMA,
		M_STEP, M_SGN, M_RECT, M_TENT,
		M_CEIL, M_FLOOR, M_ROUND, M_INT, M_FRAC,
		M_ABS, M_ARG, M_REAL, M_IMAG, M_CONJUGATE, M_POLAR, M_CARTESIAN,
		
			M_BFSTART,

		M_RAND,
		M_ATAN,
		M_LOG,
		M_BETA, M_GAMMA, M_PERMUTATION, M_COMBINATION,
		M_BESSEL_J, M_BESSEL_Y, M_HANKEL1,
		M_SQWV, M_TRWV, M_SAW, M_MIN, M_MAX, M_HYPOT, M_MANDELBROT,

		M_USER_FUNCTION
	};
}
enum 				ResultMode
{
	MODE_NO_EXPR,
	MODE_N0D,
	MODE_T1D,
	MODE_T1D_C,
	MODE_T1D_H,
	MODE_T2D,
	MODE_C2D,
	MODE_L2D,
	MODE_T2D_H,
	MODE_C3D,
	MODE_I1D,
	MODE_I2D,

	MODE_I3D,
//	MODE_N1D,
	
	//MODE_NO_EXPR,
	//MODE_N0D,
	//MODE_I1D, MODE_N1D,
	//MODE_T1D, MODE_T1D_C, MODE_T1D_H,
	//MODE_I2D,
	//MODE_T2D,
	//MODE_C2D,
	//MODE_L2D,
	//MODE_T2D_H,
	//MODE_I3D,
	//MODE_C3D,
};
struct 				ModeParameters
{
	int mode_idx,
		nExpr;//number of expressions of this type (for color condition)

	unsigned Xplaces, Yplaces, Zplaces;//product = ndrSize

	double cx, mx, cy, my, cz, mz;//example: ndr[kx]=cx+mx*kx
//	double VX, DX, VY, DY, VZ, DZ;

	//output: mode-specific
	//MODE_I2D, MODE_C2D: {gl_texture}
	//MODE_C3D: {gl_vertexbuffer, gl_indexbuffer}
	//void *output;

	int *shift_args;//Xoffset, Yoffset, Zoffset
	int *range_args;//x1, x2, y1, y2, z1, z2
};
enum				InstructionSignature
{
	SIG_NOOP,

	SIG_R_R,	SIG_C_C,	SIG_Q_Q,

	SIG_R_RR,	SIG_C_RC,	SIG_Q_RQ,
	SIG_C_CR,	SIG_C_CC,	SIG_Q_CQ,
	SIG_Q_QR,	SIG_Q_QC,	SIG_Q_QQ,

	SIG_C_R,	SIG_C_Q,
	SIG_R_C,	SIG_R_Q,

	SIG_C_RR,

				SIG_R_RC,	SIG_R_RQ,
	SIG_R_CR,	SIG_R_CC,	SIG_R_CQ,
	SIG_R_QR,	SIG_R_QC,	SIG_R_QQ,

	SIG_C_QC,

	SIG_INLINE_IF,

	SIG_CALL='c',
	SIG_BIF='b',
	SIG_BIN='B',
	SIG_JUMP='j',
	SIG_RETURN='r',
};
char 				returnMathSet_from_signature(int signature, char op1_ms, char op2_ms=0, char op3_ms=0);
struct				Map
{
	G2::Map _0;
	int _1;
	int pos, len;
	Map(int pos=0, int len=0, G2::Map _0=G2::M_IGNORED, int _1=0):_0(_0), _1(_1), pos(pos), len(len){}
};
template<int buffer_size>inline bool printValue_real		(				char (&buffer)[buffer_size], int &offset, double const &value)
{
	if(value)
	{
		offset+=sprintf_s(buffer+offset, buffer_size-offset,

			value!=value?//NAN
				(long long&)value==0x7FF8000000000010?	"+-inf"
				:										"0/0"
			:value==_HUGE?	"inf"
			:value==-_HUGE?	"-inf"
			:				"%.15g"

			, value);//
		return true;
	}
	return false;
}
template<int buffer_size>inline bool printValue_real		(				char (&buffer)[buffer_size], int &offset, double const &value, int const base)
{
	if(value)
	{
		if(value!=value)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, (long long&)value==0x7FF8000000000010?"+-inf":"0/0");
		else if(value==_HUGE)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "inf");
		else if(value==-_HUGE)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "-inf");
		else if(base==10)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "%.15g", value);
		else
		{
			const int l2int[]={0xFF, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4};
			int l2base=l2int[base];
			if(value<0)
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "-");
			const char *prefix[]={"", "B", "", "0", "0x"};
			offset+=sprintf_s(buffer+offset, buffer_size-offset, prefix[l2base]);
			int exponent=(((int*)&value)[1]>>20&0x000007FF)-1023;
			auto mantissa=(long long&)value&0x000FFFFFFFFFFFFF|0x0010000000000000;//1.<52bits>*2^<11bits>		1.m<<e
			long long digit;
			if(exponent>=-13&&exponent<=50)
			{
				if(exponent<0)
				{
					int nzeros=-(exponent+1)/l2base;
					if(nzeros)
						offset+=sprintf_s(buffer+offset, buffer_size-offset, "0.%0*d", nzeros, 0);
					else
						offset+=sprintf_s(buffer+offset, buffer_size-offset, "0.");
					for(int k=52-(l2base+exponent%l2base)%l2base;k>=0&&mantissa;k-=l2base)
					{
						digit=mantissa>>k, mantissa-=digit<<k;
						offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
					}
				}
				else for(int k=52-exponent%l2base;k>=0;k-=l2base)//>=0 <=50
				{
					digit=mantissa>>k, mantissa-=digit<<k;
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
					if(k==52-exponent)
					{
						if(mantissa)
						{
							offset+=sprintf_s(buffer+offset, buffer_size-offset, ".");
							for(k-=l2base;k>=0&&mantissa;k-=l2base)
							{
								digit=mantissa>>k, mantissa-=digit<<k;
								offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
							}
							if(mantissa)
							{
								digit=mantissa<<((l2base-(52-exponent)%l2base)%l2base), mantissa-=digit<<k;
								offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
							}
						}
						break;
					}
				}
			}
			else
			{
				int k=52;//p=2
			//	int k=52-exponent%l2base;//p=16
			//	int k=52-(l2base+exponent%l2base)%l2base;
				digit=mantissa>>k, mantissa-=digit<<k;
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X.", digit);
				for(k-=l2base;k>=0&&mantissa;k-=l2base)
				{
					digit=mantissa>>k, mantissa-=digit<<k;
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
				}
				if(exponent<0)
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "p-");
				else
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "p+");
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "%03d", exponent);
			//	exponent/=l2base;
				//for(int k=10-10%l2base, digit;k>=0&&exponent;k-=l2base)
				//{
				//	if(digit=exponent>>k)
				//	{
				//		exponent-=digit<<k;
				//		offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
				//		for(k-=l2base;k>=0&&exponent;k-=l2base)
				//		{
				//			digit=exponent>>k, exponent-=digit<<k;
				//			offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
				//		}
				//		break;
				//	}
				//}
			}
		}
		return true;
	}
	return false;
}
template<int buffer_size>inline void printValue_unreal		(bool &written, char (&buffer)[buffer_size], int &offset, double const &value, char const *component)
{
	if(value)
	{
		offset+=sprintf_s(buffer+offset, buffer_size-offset,
			
			value!=value?//NAN
				(long long&)value==0x7FF8000000000010?	"+-inf"
				:										written?"+0/0":"0/0"
			:value==_HUGE?			written?"+inf":"inf"
			:value==-_HUGE?			"-inf"
			:value==-1?				"-"
			:value==1?				written?"+":""
			:						written&&value>0?"+%.15g":"%.15g"

			, value);//
		offset+=sprintf_s(buffer+offset, buffer_size-offset, component);
		written=true;
	}
}
template<int buffer_size>inline void printValue_unreal		(bool &written, char (&buffer)[buffer_size], int &offset, double const &value, char const *component, int const base)
{
	if(value)
	{
		if(value!=value)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, (long long&)value==0x7FF8000000000010?"+-inf":written?"+0/0":"0/0");
		else if(value==_HUGE)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, written?"+inf":"inf");
		else if(value==-_HUGE)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "-inf");
		else if(value==1)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, written?"+":"");
		else if(value==-1)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "-");
		else if(base==10)
			offset+=sprintf_s(buffer+offset, buffer_size-offset, written&&value>0?"+%.15g":"%.15g", value);
		else
		{
			const int l2int[]={0x80000000, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4};
			int l2base=l2int[base];
			if(written&&value>0)
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "+");
			if(value<0)
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "-");
			const char *prefix[]={"", "B", "", "0", "0x"};
			offset+=sprintf_s(buffer+offset, buffer_size-offset, prefix[l2base]);
			int exponent=(((int*)&value)[1]>>20&0x000007FF)-1023;
			auto mantissa=(long long&)value&0x000FFFFFFFFFFFFF|0x0010000000000000;//1.<52bits>*2^<11bits>		1.m<<e
			long long digit;
			if(exponent>=-13&&exponent<=50)
			{
				if(exponent<0)
				{
					int nzeros=-(exponent+1)/l2base;
					if(nzeros)
						offset+=sprintf_s(buffer+offset, buffer_size-offset, "0.%0*d", nzeros, 0);
					else
						offset+=sprintf_s(buffer+offset, buffer_size-offset, "0.");
					for(int k=52-exponent%l2base;k>=0&&mantissa;k-=l2base)
					{
						digit=mantissa>>k, mantissa-=digit<<k;
						offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
					}
				}
				else for(int k=52-exponent%l2base;k>=0;k-=l2base)//>=0 <=50
				{
					digit=mantissa>>k, mantissa-=digit<<k;
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
					if(k==52-exponent)
					{
						if(mantissa)
						{
							offset+=sprintf_s(buffer+offset, buffer_size-offset, ".");
							for(k-=l2base;k>=0&&mantissa;k-=l2base)
							{
								digit=mantissa>>k, mantissa-=digit<<k;
								offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
							}
							if(mantissa)
							{
								digit=mantissa<<((l2base-(52-exponent)%l2base)%l2base), mantissa-=digit<<k;
								offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
							}
						}
						break;
					}
				}
			}
			else
			{
				int k=52;
			//	int k=52-(l2base+exponent%l2base)%l2base;
				digit=mantissa>>k, mantissa-=digit<<k;
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X.", digit);
				for(k-=l2base;k>=0&&mantissa;k-=l2base)
				{
					digit=mantissa>>k, mantissa-=digit<<k;
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
				}
				if(exponent<0)
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "p-");
				else
					offset+=sprintf_s(buffer+offset, buffer_size-offset, "p+");
				offset+=sprintf_s(buffer+offset, buffer_size-offset, "%03d", exponent);
			//	exponent/=l2base;
				//for(int k=10-10%l2base, digit;k>=0&&exponent;k-=l2base)
				//{
				//	if(digit=exponent>>k)
				//	{
				//		exponent-=digit<<k;
				//		offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
				//		for(k-=l2base;k>=0&&exponent;k-=l2base)
				//		{
				//			digit=exponent>>k, exponent-=digit<<k;
				//			offset+=sprintf_s(buffer+offset, buffer_size-offset, "%X", digit);
				//		}
				//		break;
				//	}
				//}
			}
		}
		offset+=sprintf_s(buffer+offset, buffer_size-offset, component);
		written=true;
	}
}
template<typename T, int Align>struct aligned_vector
{
	typedef T *iterator;
	typedef T const *const_iterator;
	T *p;
	int n;
	aligned_vector():p(0), n(0){}
	aligned_vector(aligned_vector<T, Align> const &other)
	{
		if(&other!=this)
		{
			n=other.n;
			if(n)
			{
				p=(T*)_aligned_malloc(other.n*sizeof(T), Align);
				memcpy(p, other.p, n*sizeof(T));
			}
			else
				p=0;
		}
	}
	aligned_vector(aligned_vector<T, Align> &&other)
	{
		if(&other!=this)
		{
			n=other.n, p=other.p;
			other.n=0, other.p=0;
		}
	}
	aligned_vector(unsigned n, T const &e=T()):n(n), p((T*)_aligned_malloc(n*sizeof(T), Align))
	{
		fill(0, n, e);
	}
	~aligned_vector()
	{
		if(p)
		{
			_aligned_free(p);
			p=0;
		}
	}
	aligned_vector& operator=(aligned_vector const &other)
	{
		if(&other!=this)
		{
			n=other.n;
			realloc(n);
			for(int k=0;k<n;++k)
				p[k]=other[k];
		}
		return *this;
	}
	aligned_vector& operator=(aligned_vector &&other)
	{
		if(&other!=this)
		{
			if(p)
				_aligned_free(p);
			n=other.n, p=other.p;
			other.n=0, other.p=0;
		}
		return *this;
	}
	unsigned size()const{return n;}
	iterator begin(){return p;}
	iterator end(){return p+n;}
	const_iterator begin()const{return p;}
	const_iterator end()const{return p+n;}
	void realloc(unsigned n)
	{
		if(p)
			p=(T*)_aligned_realloc(p, n*sizeof(T), Align);
		else
			p=(T*)_aligned_malloc(n*sizeof(T), Align);
	}
	void fill(unsigned start, unsigned end, T const &e)
	{
		for(unsigned k=start;k<end;++k)
			p[k]=e;
	}
	void resize(unsigned n, T const &e=T())
	{
		if(n!=this->n)//
		{
			realloc(n);
			fill(this->n, n, e);
			this->n=n;
		}
	}
	T& operator[](unsigned k){return p[k];}
	T const& operator[](unsigned k)const{return p[k];}
	void assign(unsigned n, T const &e=T())
	{
		this->n=n;
		realloc(n);
		fill(0, n, e);
	}
	void push_back(T const &e)
	{
		++n;
		realloc(n);
		p[n-1]=e;
	}
	void pop_back()
	{
		--n;
		realloc(n);
	}
	void insert(unsigned position, T const &e, unsigned count=1)
	{
		unsigned n2=n+count;
		realloc(n2);
		n=n2;
		unsigned pos2=position+count;
		for(unsigned k=n2-1;k>=pos2;--k)
			p[k]=p[k-count];
		fill(position, pos2, e);
	}
	void erase(unsigned position, unsigned count=1)
	{
		for(unsigned k=position, kEnd=n-count;k<kEnd;++k)
			p[k]=p[k+count];
		n-=count;
		realloc(n);
	}
	void clear(){n=0; if(p)_aligned_free(p), p=0;}
};
typedef aligned_vector<double, 32> AVector_v4d;
struct				CompRef
{
	double &r, &i;
	CompRef(double &r, double &i):r(r), i(i){}
	CompRef& operator=(std::complex<double> &other){r=other.real(), i=other.imag(); return *this;}
	CompRef& operator+=(std::complex<double> &other){r+=other.real(), i+=other.imag(); return *this;}
	CompRef& operator-=(std::complex<double> &other){r-=other.real(), i-=other.imag(); return *this;}
	CompRef& operator*=(double other){r*=other, i*=other; return *this;}
	CompRef& operator/=(double other){r/=other, i/=other; return *this;}
};
inline std::complex<double> operator+(CompRef const &a, CompRef const &b){return std::complex<double>(a.r+b.r, a.i+b.i);}
inline std::complex<double> operator-(CompRef const &a, CompRef const &b){return std::complex<double>(a.r-b.r, a.i-b.i);}
inline std::complex<double> operator*(CompRef const &a, CompRef const &b){return std::complex<double>(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r);}
inline std::complex<double> operator/(CompRef const &a, CompRef const &b)
{
	double _1_mag_b=1/sqrt(b.r*b.r+b.i*b.i);
	return std::complex<double>((b.r*a.r+b.i*a.i)*_1_mag_b, (b.r*a.i-b.i*a.r)*_1_mag_b);
}
struct				QuatRef
{
	double &r, &i, &j, &k;
	QuatRef(double &r, double &i, double &j, double &k):r(r), i(i), j(j), k(k){}
	QuatRef& operator=(boost::math::quaternion<double> &other){r=other.R_component_1(), i=other.R_component_2(), j=other.R_component_3(), k=other.R_component_4(); return *this;}
	QuatRef& operator+=(boost::math::quaternion<double> &other){r+=other.R_component_1(), i+=other.R_component_2(), j+=other.R_component_3(), k+=other.R_component_4(); return *this;}
	QuatRef& operator-=(boost::math::quaternion<double> &other){r-=other.R_component_1(), i-=other.R_component_2(), j-=other.R_component_3(), k-=other.R_component_4(); return *this;}
};
inline boost::math::quaternion<double> operator+(QuatRef const &a, QuatRef const &b){return boost::math::quaternion<double>(a.r+b.r, a.i+b.i, a.j+b.j, a.k+b.k);}
inline boost::math::quaternion<double> operator-(QuatRef const &a, QuatRef const &b){return boost::math::quaternion<double>(a.r-b.r, a.i-b.i, a.j-b.j, a.k-b.k);}
inline boost::math::quaternion<double> operator*(QuatRef const &a, QuatRef const &b)
{
	return boost::math::quaternion<double>(
		a.r*b.r+a.i*b.i+a.j*b.j+a.k*b.k,
		a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,
		a.r*b.j-a.i*b.k+a.j*b.r+a.k*b.i,
		a.r*b.k+a.i*b.j+a.j*b.i+a.k*b.r);
}
inline boost::math::quaternion<double> operator/(QuatRef const &a, QuatRef const &b)
{
	double _1_mag_y=1/(b.r*b.r+b.i*b.i+b.j*b.j+b.k*b.k);
	return boost::math::quaternion<double>(
		(b.r*a.r+b.i*a.i+b.j*a.j+b.k*a.k)*_1_mag_y,
		(b.r*a.i-b.i*a.r-b.j*a.k+b.k*a.j)*_1_mag_y,
		(b.r*a.j+b.i*a.k-b.j*a.r-b.k*a.i)*_1_mag_y,
		(b.r*a.k-b.i*a.j+b.j*a.i-b.k*a.r)*_1_mag_y);
}
struct				Value
{
	double r, i, j, k;//in order
	Value(double r=0, double i=0, double j=0, double k=0):r(r), i(i), j(j), k(k){}
	Value(std::complex<double> const &x):r(x.real()), i(x.imag()), j(0), k(0){}
	Value(boost::math::quaternion<double> const &x):r(x.R_component_1()), i(x.R_component_2()), j(x.R_component_3()), k(x.R_component_4()){}
	Value(MP::Quat const &x):r(x.r.toDouble()), i(x.i.toDouble()), j(x.j.toDouble()), k(x.k.toDouble()){}
	void set(double r, double i=0, double j=0, double k=0){this->r=r, this->i=i, this->j=j, this->k=k;}
	void setzero(){r=i=j=k=0;}
	void setnan(){r=G2::_qnan, i=j=k=0;}
	operator double&						(){return r;}
	operator double							()const{return r;}
	operator std::complex<double>			()const{return std::complex<double>(r, i);}
	operator boost::math::quaternion<double>()const{return boost::math::quaternion<double>(r, i, j, k);}
	operator MP::Quat						()const{return MP::Quat(r, i, j, k);}
	Value& operator=	(double							const &x){r=x, i=j=k=0;																				return *this;}
	Value& operator=	(std::complex<double>			const &x){r=x.real(), i=x.imag(), j=k=0;															return *this;}
	Value& operator=	(boost::math::quaternion<double>const &x){r=x.R_component_1(), i=x.R_component_2(), j=x.R_component_3(), k=x.R_component_4();		return *this;}
	Value& operator=	(MP::Quat const &x						){r=x.r.toDouble(), i=x.i.toDouble(), j=x.j.toDouble(), k=x.k.toDouble();					return *this;}
	Value& operator+=	(double							const &x){r+=x;																						return *this;}
	Value& operator+=	(std::complex<double>			const &x){r+=x.real(), i+=x.imag();																	return *this;}
	Value& operator+=	(boost::math::quaternion<double>const &x){r+=x.R_component_1(), i+=x.R_component_2(), j+=x.R_component_3(), k+=x.R_component_4();	return *this;}
	Value& operator-=	(double							const &x){r-=x;																						return *this;}
	Value& operator-=	(std::complex<double>			const &x){r-=x.real(), i-=x.imag();																	return *this;}
	Value& operator-=	(boost::math::quaternion<double>const &x){r-=x.R_component_1(), i-=x.R_component_2(), j-=x.R_component_3(), k-=x.R_component_4();	return *this;}
	bool r_isTrue()const{return r!=0;}
	bool c_isTrue()const{return r!=0||i!=0;}
	bool q_isTrue()const{return r!=0||i!=0||j!=0||k!=0;}
	template<int buffer_size>void printReal(char (&buffer)[buffer_size], int &offset)const
	{
		if(!printValue_real(buffer, offset, r))
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "0");
	}
	template<int buffer_size>void printReal(char (&buffer)[buffer_size], int &offset, int base)const
	{
		if(!printValue_real(buffer, offset, r, base))
			offset+=sprintf_s(buffer+offset, buffer_size-offset, "0");
	}
	template<int _size>void printComplex(char(&)[_size], int&)const;
	template<int _size>void printComplex(char(&)[_size], int&, int)const;
	template<int _size>void printQuaternion(char(&)[_size], int&)const;
	template<int _size>void printQuaternion(char(&)[_size], int&, int)const;
	template<int _size>void print(char (&a)[_size], int &o, char mathSet)const
	{
		switch(mathSet)
		{
		case 'R':
			printReal(a, o);
			break;
		case 'c':
			printComplex(a, o);
			break;
		case 'h':
			printQuaternion(a, o);
			break;
		}
	}
	template<int _size>void print(char (&a)[_size], int &o, char mathSet, int base)const
	{
		switch(mathSet)
		{
		case 'R':
			printReal(a, o, base);
			break;
		case 'c':
			printComplex(a, o, base);
			break;
		case 'h':
			printQuaternion(a, o, base);
			break;
		}
	}
};
template<int buffer_size>inline void Value::printComplex	(				char (&buffer)[buffer_size], int &offset			)const
{
	bool written=printValue_real(buffer, offset, r);
	printValue_unreal(written, buffer, offset, i, "i");
	if(!written)
		offset+=sprintf_s(buffer+offset, buffer_size-offset, "0");
}
template<int buffer_size>inline void Value::printComplex	(				char (&buffer)[buffer_size], int &offset, int base	)const
{
	bool written=printValue_real(buffer, offset, r, base);
	printValue_unreal(written, buffer, offset, i, "i", base);
	if(!written)
		offset+=sprintf_s(buffer+offset, buffer_size-offset, "0");
}
template<int buffer_size>inline void Value::printQuaternion(				char (&buffer)[buffer_size], int &offset			)const
{
	bool written=printValue_real(buffer, offset, r);
	printValue_unreal(written, buffer, offset, i, "i");
	printValue_unreal(written, buffer, offset, j, "j");
	printValue_unreal(written, buffer, offset, k, "k");
	if(!written)
		offset+=sprintf_s(buffer+offset, buffer_size-offset, "0");
}
template<int buffer_size>inline void Value::printQuaternion(				char (&buffer)[buffer_size], int &offset, int base	)const
{
	bool written=printValue_real(buffer, offset, r, base);
	printValue_unreal(written, buffer, offset, i, "i", base);
	printValue_unreal(written, buffer, offset, j, "j", base);
	printValue_unreal(written, buffer, offset, k, "k", base);
	if(!written)
		offset+=sprintf_s(buffer+offset, buffer_size-offset, "0");
}
struct				DiscontinuityFunction
{
	bool disc_in, disc_out;//evaluated before/after the function
	union
	{
		bool
			(*d_o)(Value const&, Value const&),
			(*ud_i)(Value const&, Value const&),
			(*bd_i)(Value const&, Value const&, Value const&, Value const&),
			(*td_i)(Value const&, Value const&, Value const&, Value const&, Value const&, Value const&);
	};

	//continuous function
	DiscontinuityFunction():
		disc_in(false), disc_out(false){}

	//discontinuity takes unary function argument or any function's output
	DiscontinuityFunction(bool (*d)(Value const&, Value const&), bool disc_in):
		disc_in(disc_in), disc_out(!disc_in), ud_i(d){}

	//discontinuity takes binary function arguments
	DiscontinuityFunction(bool (*bd_i)(Value const&, Value const&, Value const&, Value const&)):
		disc_in(true), disc_out(false), bd_i(bd_i){}

	//discontinuity takes triary function arguments
	DiscontinuityFunction(bool (*td_i)(Value const&, Value const&, Value const&, Value const&, Value const&, Value const&)):
		disc_in(true), disc_out(false), td_i(td_i){}

	//continuous function
	void operator()()
	{disc_out=disc_in=false;}

	//discontinuity takes unary function argument or any function's output
	void operator()(bool (*d)(Value const&, Value const&), bool disc_in=true)
	{disc_out=!(this->disc_in=disc_in), ud_i=d;}

	//discontinuity takes binary function arguments
	void operator()(bool (*bd_i)(Value const&, Value const&, Value const&, Value const&))
	{disc_out=!(disc_in=true), this->bd_i=bd_i;}

	//discontinuity takes triary function arguments
	void operator()(bool (*td_i)(Value const&, Value const&, Value const&, Value const&, Value const&, Value const&))
	{disc_out=!(disc_in=true), this->td_i=td_i;}
};
namespace			MP
{
	union			FP//multiprecision function pointer
	{
		//enum	Type
		//{
		//	NOOP_FUNCTION,
		//	UNARY__R_R_FUNC, UNARY__C_C_FUNC, UNARY__Q_Q_FUNC,
		//	BINARY_R_RR_FUNC, BINARY_C_RC_FUNC, BINARY_Q_RQ_FUNC,
		//	BINARY_C_CR_FUNC, BINARY_C_CC_FUNC, BINARY_Q_CQ_FUNC,
		//	BINARY_Q_QR_FUNC, BINARY_Q_QC_FUNC, BINARY_Q_QQ_FUNC,
		//	UNARY__C_R_FUNC, UNARY__C_Q_FUNC,
		//	UNARY__R_C_FUNC, UNARY__R_Q_FUNC,
		//	BINARY_C_RR_FUNC, BINARY_R_RC_FUNC, BINARY_R_RQ_FUNC,
		//	BINARY_R_CR_FUNC, BINARY_R_CC_FUNC, BINARY_R_CQ_FUNC,
		//	BINARY_R_QR_FUNC, BINARY_R_QC_FUNC, BINARY_R_QQ_FUNC,
		//	BINARY_C_QC_FUNC,
		//	INLINE_IF
		//};
		typedef void (*UF)(Quat &r, Quat const &x);//1, 2, 3, 13, 14, 15, 16
		typedef void (*BF)(Quat &r, Quat const &x, Quat const &y);//4,5,6,7,8,9,10,11,12, 17, 18,19,20,21,22,23,24,25, 26
		UF uf;
		BF bf;
		FP():uf(nullptr){}
		FP(void *p):uf((void(*)(Quat&, const Quat&))p){}
		FP& operator=(void *p){uf=(void(*)(Quat&, const Quat&))p; return *this;}//*/
	/*	enum	Type
		{
			NOOP_FUNCTION,
			_R_R_FUNCTION, _C_C_FUNCTION, _Q_Q_FUNCTION,
			R_RR_FUNCTION, C_RC_FUNCTION, Q_RQ_FUNCTION,
			C_CR_FUNCTION, C_CC_FUNCTION, Q_CQ_FUNCTION,
			Q_QR_FUNCTION, Q_QC_FUNCTION, Q_QQ_FUNCTION,
			_C_R_FUNCTION, _C_Q_FUNCTION,
			_R_C_FUNCTION, _R_Q_FUNCTION,
			C_RR_FUNCTION, R_RC_FUNCTION, R_RQ_FUNCTION,
			R_CR_FUNCTION, R_CC_FUNCTION, R_CQ_FUNCTION,
			R_QR_FUNCTION, R_QC_FUNCTION, R_QQ_FUNCTION,
			C_QC_FUNCTION,
			INLINE_IF
		};
		typedef void (*R_R)(Real &r, Real const &x);//1
		typedef void (*C_C)(Comp &r, Comp const &x);//2
		typedef void (*Q_Q)(Quat &r, Quat const &x);//3

		typedef void (*R_RR)(Real &r, Real const &x, Real const &y);//4
		typedef void (*C_RC)(Comp &r, Real const &x, Comp const &y);//5
		typedef void (*Q_RQ)(Quat &r, Real const &x, Quat const &y);//6
		typedef void (*C_CR)(Comp &r, Comp const &x, Real const &y);//7
		typedef void (*C_CC)(Comp &r, Comp const &x, Comp const &y);//8
		typedef void (*Q_CQ)(Quat &r, Comp const &x, Quat const &y);//9
		typedef void (*Q_QR)(Quat &r, Quat const &x, Real const &y);//10
		typedef void (*Q_QC)(Quat &r, Quat const &x, Comp const &y);//11
		typedef void (*Q_QQ)(Quat &r, Quat const &x, Quat const &y);//12

		typedef void (*C_R)(Comp &r, Real const &x);//13
		typedef void (*C_Q)(Comp &r, Quat const &x);//14

		typedef void (*R_C)(Real &r, Comp const &x);//15
		typedef void (*R_Q)(Real &r, Quat const &x);//16

		typedef void (*C_RR)(Comp &r, Real const &x, Real const &y);//17	//pow, tetrate
													 
		typedef void (*R_RC)(Real &r, Real const &x, Comp const &y);//18	//&& || ##
		typedef void (*R_RQ)(Real &r, Real const &x, Quat const &y);//19
		typedef void (*R_CR)(Real &r, Comp const &x, Real const &y);//20
		typedef void (*R_CC)(Real &r, Comp const &x, Comp const &y);//21
		typedef void (*R_CQ)(Real &r, Comp const &x, Quat const &y);//22
		typedef void (*R_QR)(Real &r, Quat const &x, Real const &y);//23
		typedef void (*R_QC)(Real &r, Quat const &x, Comp const &y);//24
		typedef void (*R_QQ)(Real &r, Quat const &x, Quat const &y);//25
													 
		typedef void (*C_QC)(Comp &r, Quat const &x, Comp const &y);//26	//conditional 110 & 101
		R_R r_r;//1
		C_C c_c;//2
		Q_Q q_q;//3

		R_RR r_rr;//4
		C_RC c_rc;//5
		Q_RQ q_rq;//6
		C_CR c_cr;//7
		C_CC c_cc;//8
		Q_CQ q_cq;//9
		Q_QR q_qr;//10
		Q_QC q_qc;//11
		Q_QQ q_qq;//12

		C_R c_r;//13
		C_Q c_q;//14

		R_C r_c;//15
		R_Q r_q;//16

		C_RR c_rr;//17	//pow, tetrate

		R_RC r_rc;//18	//&& || ##
		R_RQ r_rq;//19
		R_CR r_cr;//20
		R_CC r_cc;//21
		R_CQ r_cq;//22
		R_QR r_qr;//23
		R_QC r_qc;//24
		R_QQ r_qq;//25

		C_QC c_qc;//26	//conditional 110 & 101

		FP():r_r(nullptr){}
		FP(void *p):r_r((void(*)(Real&, const Real&))p){}
		FP& operator=(void *p){r_r=(void(*)(Real&, const Real&))p; return *this;}//*/
	};
}
union				FP//function pointer
{
	//enum			Type//TODO: remove this, use enum InstructionSignature instead
	//{
	//	NOOP_FUNCTION,
	//	_R_R_FUNCTION, _C_C_FUNCTION, _Q_Q_FUNCTION,
	//	R_RR_FUNCTION, C_RC_FUNCTION, Q_RQ_FUNCTION,
	//	C_CR_FUNCTION, C_CC_FUNCTION, Q_CQ_FUNCTION,
	//	Q_QR_FUNCTION, Q_QC_FUNCTION, Q_QQ_FUNCTION,
	//	_C_R_FUNCTION, _C_Q_FUNCTION,
	//	_R_C_FUNCTION, _R_Q_FUNCTION,
	//	C_RR_FUNCTION, R_RC_FUNCTION, R_RQ_FUNCTION,
	//	R_CR_FUNCTION, R_CC_FUNCTION, R_CQ_FUNCTION,
	//	R_QR_FUNCTION, R_QC_FUNCTION, R_QQ_FUNCTION,
	//	C_QC_FUNCTION,
	//	INLINE_IF
	//};
	typedef void (*R_R)(VectP &r, VectP const &x);//1
	typedef void (*C_C)(CompP &r, CompP const &x);//2
	typedef void (*Q_Q)(QuatP &r, QuatP const &x);//3

	typedef void (*R_RR)(VectP &r, VectP const &x, VectP const &y);//4
	typedef void (*C_RC)(CompP &r, VectP const &x, CompP const &y);//5
	typedef void (*Q_RQ)(QuatP &r, VectP const &x, QuatP const &y);//6
	typedef void (*C_CR)(CompP &r, CompP const &x, VectP const &y);//7
	typedef void (*C_CC)(CompP &r, CompP const &x, CompP const &y);//8
	typedef void (*Q_CQ)(QuatP &r, CompP const &x, QuatP const &y);//9
	typedef void (*Q_QR)(QuatP &r, QuatP const &x, VectP const &y);//10
	typedef void (*Q_QC)(QuatP &r, QuatP const &x, CompP const &y);//11
	typedef void (*Q_QQ)(QuatP &r, QuatP const &x, QuatP const &y);//12

	typedef void (*C_R)(CompP &r, VectP const &x);//13
	typedef void (*C_Q)(CompP &r, QuatP const &x);//14

	typedef void (*R_C)(VectP &r, CompP const &x);//15
	typedef void (*R_Q)(VectP &r, QuatP const &x);//16

	typedef void (*C_RR)(CompP &r, VectP const &x,	VectP const &y);//17	//pow, tetrate

	typedef void (*R_RC)(VectP &r, VectP const &x,	CompP const &y);//18	//&& || ##
	typedef void (*R_RQ)(VectP &r, VectP const &x,	QuatP const &y);//19
	typedef void (*R_CR)(VectP &r, CompP const &x,	VectP const &y);//20
	typedef void (*R_CC)(VectP &r, CompP const &x,	CompP const &y);//21
	typedef void (*R_CQ)(VectP &r, CompP const &x,	QuatP const &y);//22
	typedef void (*R_QR)(VectP &r, QuatP const &x,	VectP const &y);//23
	typedef void (*R_QC)(VectP &r, QuatP const &x,	CompP const &y);//24
	typedef void (*R_QQ)(VectP &r, QuatP const &x,	QuatP const &y);//25

	typedef void (*C_QC)(CompP &r, QuatP const &x,	CompP const &y);//26	//conditional 110 & 101
	R_R r_r;//1
	C_C c_c;//2
	Q_Q q_q;//3

	R_RR r_rr;//4
	C_RC c_rc;//5
	Q_RQ q_rq;//6
	C_CR c_cr;//7
	C_CC c_cc;//8
	Q_CQ q_cq;//9
	Q_QR q_qr;//10
	Q_QC q_qc;//11
	Q_QQ q_qq;//12

	C_R c_r;//13
	C_Q c_q;//14

	R_C r_c;//15
	R_Q r_q;//16

	C_RR c_rr;//17	//pow, tetrate

	R_RC r_rc;//18	//&& || ##
	R_RQ r_rq;//19
	R_CR r_cr;//20
	R_CC r_cc;//21
	R_CQ r_cq;//22
	R_QR r_qr;//23
	R_QC r_qc;//24
	R_QQ r_qq;//25

	C_QC c_qc;//26	//conditional 110 & 101

	//Vect2d (*r_r)(Vect2d const&);//1
	//Comp2d (*c_c)(Comp2d const&);//2
	//Quat2d (*q_q)(Quat2d const&);//3
	//
	//Vect2d (*r_rr)(Vect2d const&, Vect2d const&);//4
	//Comp2d (*c_rc)(Vect2d const&, Comp2d const&);//5
	//Quat2d (*q_rq)(Vect2d const&, Quat2d const&);//6
	//Comp2d (*c_cr)(Comp2d const&, Vect2d const&);//7
	//Comp2d (*c_cc)(Comp2d const&, Comp2d const&);//8
	//Quat2d (*q_cq)(Comp2d const&, Quat2d const&);//9
	//Quat2d (*q_qr)(Quat2d const&, Vect2d const&);//10
	//Quat2d (*q_qc)(Quat2d const&, Comp2d const&);//11
	//Quat2d (*q_qq)(Quat2d const&, Quat2d const&);//12
	//
	//Comp2d (*c_r)(Vect2d const &x);//13
	//Comp2d (*c_q)(Quat2d const &x);//14
	//
	//Vect2d (*r_c)(Comp2d const &x);//15
	//Vect2d (*r_q)(Quat2d const &x);//16
	//
	//Comp2d (*c_rr)(Vect2d const &x,	Vect2d const &y);//17	//pow, tetrate
	//
	//Vect2d (*r_rc)(Vect2d const &x,	Comp2d const &y);//18	//&& || ##
	//Vect2d (*r_rq)(Vect2d const &x,	Quat2d const &y);//19
	//Vect2d (*r_cr)(Comp2d const &x,	Vect2d const &y);//20
	//Vect2d (*r_cc)(Comp2d const &x,	Comp2d const &y);//21
	//Vect2d (*r_cq)(Comp2d const &x,	Quat2d const &y);//22
	//Vect2d (*r_qr)(Quat2d const &x,	Vect2d const &y);//23
	//Vect2d (*r_qc)(Quat2d const &x,	Comp2d const &y);//24
	//Vect2d (*r_qq)(Quat2d const &x,	Quat2d const &y);//25
	//
	//Comp2d (*c_qc)(Quat2d const &x,	Comp2d const &y);//26	//conditional 110 & 101

	FP():r_r(nullptr){}
	FP(void *p):r_r((void(*)(VectP &,const VectP &))p){}
	FP& operator=(void *p){r_r=(void(*)(VectP &,const VectP &))p; return *this;}
};
struct				FPSetter
{
	char type;
	bool simd;
	FP ia32, sse2, avx;
	MP::FP mp;
	FPSetter():ia32(0), sse2(0), avx(0), mp(0), type(0){}
	void set(){ia32=sse2=avx=0, mp=0, type=SIG_NOOP;}

	//non-simd setters
	void set(FP::R_R r_r,	MP::FP::UF mp_r_r)	{ia32.r_r=sse2.r_r=avx.r_r=r_r,		mp=mp_r_r,	type=SIG_R_R, simd=false;}
	void set(FP::C_C c_c,	MP::FP::UF mp_c_c)	{ia32.c_c=sse2.c_c=avx.c_c=c_c,		mp=mp_c_c,	type=SIG_C_C, simd=false;}
	void set(FP::Q_Q q_q,	MP::FP::UF mp_q_q)	{ia32.q_q=sse2.q_q=avx.q_q=q_q,		mp=mp_q_q,	type=SIG_Q_Q, simd=false;}
	void set(FP::R_RR r_rr,	MP::FP::BF mp_r_rr)	{ia32.r_rr=sse2.r_rr=avx.r_rr=r_rr, mp=mp_r_rr, type=SIG_R_RR, simd=false;}
	void set(FP::C_RC c_rc,	MP::FP::BF mp_c_rc)	{ia32.c_rc=sse2.c_rc=avx.c_rc=c_rc, mp=mp_c_rc, type=SIG_C_RC, simd=false;}
	void set(FP::Q_RQ q_rq,	MP::FP::BF mp_q_rq)	{ia32.q_rq=sse2.q_rq=avx.q_rq=q_rq, mp=mp_q_rq, type=SIG_Q_RQ, simd=false;}
	void set(FP::C_CR c_cr,	MP::FP::BF mp_c_cr)	{ia32.c_cr=sse2.c_cr=avx.c_cr=c_cr, mp=mp_c_cr, type=SIG_C_CR, simd=false;}
	void set(FP::C_CC c_cc,	MP::FP::BF mp_c_cc)	{ia32.c_cc=sse2.c_cc=avx.c_cc=c_cc, mp=mp_c_cc, type=SIG_C_CC, simd=false;}
	void set(FP::Q_CQ q_cq,	MP::FP::BF mp_q_cq)	{ia32.q_cq=sse2.q_cq=avx.q_cq=q_cq, mp=mp_q_cq, type=SIG_Q_CQ, simd=false;}
	void set(FP::Q_QR q_qr,	MP::FP::BF mp_q_qr)	{ia32.q_qr=sse2.q_qr=avx.q_qr=q_qr, mp=mp_q_qr, type=SIG_Q_QR, simd=false;}
	void set(FP::Q_QC q_qc,	MP::FP::BF mp_q_qc)	{ia32.q_qc=sse2.q_qc=avx.q_qc=q_qc, mp=mp_q_qc, type=SIG_Q_QC, simd=false;}
	void set(FP::Q_QQ q_qq,	MP::FP::BF mp_q_qq)	{ia32.q_qq=sse2.q_qq=avx.q_qq=q_qq, mp=mp_q_qq, type=SIG_Q_QQ, simd=false;}

	void set(FP::C_R c_r,	MP::FP::UF mp_c_r)	{ia32.c_r=sse2.c_r=avx.c_r=c_r,		mp=mp_c_r,	type=SIG_C_R, simd=false;}
	void set(FP::C_Q c_q,	MP::FP::UF mp_c_q)	{ia32.c_q=sse2.c_q=avx.c_q=c_q,		mp=mp_c_q,	type=SIG_C_Q, simd=false;}

	void set(FP::R_C r_c,	MP::FP::UF mp_r_c)	{ia32.r_c=sse2.r_c=avx.r_c=r_c,		mp=mp_r_c,	type=SIG_R_C, simd=false;}
	void set(FP::R_Q r_q,	MP::FP::UF mp_r_q)	{ia32.r_q=sse2.r_q=avx.r_q=r_q,		mp=mp_r_q,	type=SIG_R_Q, simd=false;}

	void set(FP::C_RR c_rr,	MP::FP::BF mp_c_rr)	{ia32.c_rr=sse2.c_rr=avx.c_rr=c_rr, mp=mp_c_rr,	type=SIG_C_RR, simd=false;}
	
	void set(FP::R_RC r_rc,	MP::FP::BF mp_r_rc)	{ia32.r_rc=sse2.r_rc=avx.r_rc=r_rc, mp=mp_r_rc,	type=SIG_C_RC, simd=false;}
	void set(FP::R_RQ r_rq,	MP::FP::BF mp_r_rq)	{ia32.r_rq=sse2.r_rq=avx.r_rq=r_rq, mp=mp_r_rq,	type=SIG_R_RQ, simd=false;}
	void set(FP::R_CR r_cr,	MP::FP::BF mp_r_cr)	{ia32.r_cr=sse2.r_cr=avx.r_cr=r_cr, mp=mp_r_cr,	type=SIG_R_CR, simd=false;}
	void set(FP::R_CC r_cc,	MP::FP::BF mp_r_cc)	{ia32.r_cc=sse2.r_cc=avx.r_cc=r_cc, mp=mp_r_cc,	type=SIG_R_CC, simd=false;}
	void set(FP::R_CQ r_cq,	MP::FP::BF mp_r_cq)	{ia32.r_cq=sse2.r_cq=avx.r_cq=r_cq, mp=mp_r_cq,	type=SIG_R_CQ, simd=false;}
	void set(FP::R_QR r_qr,	MP::FP::BF mp_r_qr)	{ia32.r_qr=sse2.r_qr=avx.r_qr=r_qr, mp=mp_r_qr,	type=SIG_R_QR, simd=false;}
	void set(FP::R_QC r_qc,	MP::FP::BF mp_r_qc)	{ia32.r_qc=sse2.r_qc=avx.r_qc=r_qc, mp=mp_r_qc,	type=SIG_R_QC, simd=false;}
	void set(FP::R_QQ r_qq,	MP::FP::BF mp_r_qq)	{ia32.r_qq=sse2.r_qq=avx.r_qq=r_qq, mp=mp_r_qq,	type=SIG_R_QQ, simd=false;}

	void set(FP::C_QC c_qc,	MP::FP::BF mp_c_qc)	{ia32.c_qc=sse2.c_qc=avx.c_qc=c_qc, mp=mp_c_qc,	type=SIG_C_QC, simd=false;}

	//simd setters
	void set(FP::R_R ia32_r_r,		FP::R_R sse2_r_r,	FP::R_R avx_r_r,	MP::FP::UF mp_r_r)	{ia32.r_r=ia32_r_r,		sse2.r_r=sse2_r_r,		avx.r_r=avx_r_r,	mp=mp_r_r,	type=SIG_R_R, simd=true;}
	void set(FP::C_C ia32_c_c,		FP::C_C sse2_c_c,	FP::C_C avx_c_c,	MP::FP::UF mp_c_c)	{ia32.c_c=ia32_c_c,		sse2.c_c=sse2_c_c,		avx.c_c=avx_c_c,	mp=mp_c_c,	type=SIG_C_C, simd=true;}
	void set(FP::Q_Q ia32_q_q,		FP::Q_Q sse2_q_q,	FP::Q_Q avx_q_q,	MP::FP::UF mp_q_q)	{ia32.q_q=ia32_q_q,		sse2.q_q=sse2_q_q,		avx.q_q=avx_q_q,	mp=mp_q_q,	type=SIG_Q_Q, simd=true;}
	void set(FP::R_RR ia32_r_rr,	FP::R_RR sse2_r_rr, FP::R_RR avx_r_rr,	MP::FP::BF mp_r_rr)	{ia32.r_rr=ia32_r_rr,	sse2.r_rr=sse2_r_rr,	avx.r_rr=avx_r_rr,	mp=mp_r_rr,	type=SIG_R_RR, simd=true;}
	void set(FP::C_RC ia32_c_rc,	FP::C_RC sse2_c_rc, FP::C_RC avx_c_rc,	MP::FP::BF mp_c_rc)	{ia32.c_rc=ia32_c_rc,	sse2.c_rc=sse2_c_rc,	avx.c_rc=avx_c_rc,	mp=mp_c_rc,	type=SIG_C_RC, simd=true;}
	void set(FP::Q_RQ ia32_q_rq,	FP::Q_RQ sse2_q_rq, FP::Q_RQ avx_q_rq,	MP::FP::BF mp_q_rq)	{ia32.q_rq=ia32_q_rq,	sse2.q_rq=sse2_q_rq,	avx.q_rq=avx_q_rq,	mp=mp_q_rq,	type=SIG_Q_RQ, simd=true;}
	void set(FP::C_CR ia32_c_cr,	FP::C_CR sse2_c_cr, FP::C_CR avx_c_cr,	MP::FP::BF mp_c_cr)	{ia32.c_cr=ia32_c_cr,	sse2.c_cr=sse2_c_cr,	avx.c_cr=avx_c_cr,	mp=mp_c_cr,	type=SIG_C_CR, simd=true;}
	void set(FP::C_CC ia32_c_cc,	FP::C_CC sse2_c_cc, FP::C_CC avx_c_cc,	MP::FP::BF mp_c_cc)	{ia32.c_cc=ia32_c_cc,	sse2.c_cc=sse2_c_cc,	avx.c_cc=avx_c_cc,	mp=mp_c_cc,	type=SIG_C_CC, simd=true;}
	void set(FP::Q_CQ ia32_q_cq,	FP::Q_CQ sse2_q_cq, FP::Q_CQ avx_q_cq,	MP::FP::BF mp_q_cq)	{ia32.q_cq=ia32_q_cq,	sse2.q_cq=sse2_q_cq,	avx.q_cq=avx_q_cq,	mp=mp_q_cq,	type=SIG_Q_CQ, simd=true;}
	void set(FP::Q_QR ia32_q_qr,	FP::Q_QR sse2_q_qr, FP::Q_QR avx_q_qr,	MP::FP::BF mp_q_qr)	{ia32.q_qr=ia32_q_qr,	sse2.q_qr=sse2_q_qr,	avx.q_qr=avx_q_qr,	mp=mp_q_qr,	type=SIG_Q_QR, simd=true;}
	void set(FP::Q_QC ia32_q_qc,	FP::Q_QC sse2_q_qc, FP::Q_QC avx_q_qc,	MP::FP::BF mp_q_qc)	{ia32.q_qc=ia32_q_qc,	sse2.q_qc=sse2_q_qc,	avx.q_qc=avx_q_qc,	mp=mp_q_qc,	type=SIG_Q_QC, simd=true;}
	void set(FP::Q_QQ ia32_q_qq,	FP::Q_QQ sse2_q_qq, FP::Q_QQ avx_q_qq,	MP::FP::BF mp_q_qq)	{ia32.q_qq=ia32_q_qq,	sse2.q_qq=sse2_q_qq,	avx.q_qq=avx_q_qq,	mp=mp_q_qq,	type=SIG_Q_QQ, simd=true;}

	void set(FP::C_R ia32_c_r,		FP::C_R sse2_c_r,	FP::C_R avx_c_r,	MP::FP::UF mp_c_r)	{ia32.c_r=ia32_c_r,		sse2.c_r=sse2_c_r,		avx.c_r=avx_c_r,	mp=mp_c_r,	type=SIG_C_R, simd=true;}
	void set(FP::C_Q ia32_c_q,		FP::C_Q sse2_c_q,	FP::C_Q avx_c_q,	MP::FP::UF mp_c_q)	{ia32.c_q=ia32_c_q,		sse2.c_q=sse2_c_q,		avx.c_q=avx_c_q,	mp=mp_c_q,	type=SIG_C_Q, simd=true;}

	void set(FP::R_C ia32_r_c,		FP::R_C sse2_r_c,	FP::R_C avx_r_c,	MP::FP::UF mp_r_c)	{ia32.r_c=ia32_r_c,		sse2.r_c=sse2_r_c,		avx.r_c=avx_r_c,	mp=mp_r_c,	type=SIG_R_C, simd=true;}
	void set(FP::R_Q ia32_r_q,		FP::R_Q sse2_r_q,	FP::R_Q avx_r_q,	MP::FP::UF mp_r_q)	{ia32.r_q=ia32_r_q,		sse2.r_q=sse2_r_q,		avx.r_q=avx_r_q,	mp=mp_r_q,	type=SIG_R_Q, simd=true;}

	void set(FP::C_RR ia32_c_rr,	FP::C_RR sse2_c_rr,	FP::C_RR avx_c_rr,	MP::FP::BF mp_c_rr)	{ia32.c_rr=ia32_c_rr,	sse2.c_rr=sse2_c_rr,	avx.c_rr=avx_c_rr,	mp=mp_c_rr,	type=SIG_C_RR, simd=true;}
	
	void set(FP::R_RC ia32_r_rc,	FP::R_RC sse2_r_rc,	FP::R_RC avx_r_rc,	MP::FP::BF mp_r_rc)	{ia32.r_rc=ia32_r_rc,	sse2.r_rc=sse2_r_rc,	avx.r_rc=avx_r_rc,	mp=mp_r_rc,	type=SIG_R_RC, simd=true;}
	void set(FP::R_RQ ia32_r_rq,	FP::R_RQ sse2_r_rq,	FP::R_RQ avx_r_rq,	MP::FP::BF mp_r_rq)	{ia32.r_rq=ia32_r_rq,	sse2.r_rq=sse2_r_rq,	avx.r_rq=avx_r_rq,	mp=mp_r_rq,	type=SIG_R_RQ, simd=true;}
	void set(FP::R_CR ia32_r_cr,	FP::R_CR sse2_r_cr,	FP::R_CR avx_r_cr,	MP::FP::BF mp_r_cr)	{ia32.r_cr=ia32_r_cr,	sse2.r_cr=sse2_r_cr,	avx.r_cr=avx_r_cr,	mp=mp_r_cr,	type=SIG_R_CR, simd=true;}
	void set(FP::R_CC ia32_r_cc,	FP::R_CC sse2_r_cc,	FP::R_CC avx_r_cc,	MP::FP::BF mp_r_cc)	{ia32.r_cc=ia32_r_cc,	sse2.r_cc=sse2_r_cc,	avx.r_cc=avx_r_cc,	mp=mp_r_cc,	type=SIG_R_CC, simd=true;}
	void set(FP::R_CQ ia32_r_cq,	FP::R_CQ sse2_r_cq,	FP::R_CQ avx_r_cq,	MP::FP::BF mp_r_cq)	{ia32.r_cq=ia32_r_cq,	sse2.r_cq=sse2_r_cq,	avx.r_cq=avx_r_cq,	mp=mp_r_cq,	type=SIG_R_CQ, simd=true;}
	void set(FP::R_QR ia32_r_qr,	FP::R_QR sse2_r_qr,	FP::R_QR avx_r_qr,	MP::FP::BF mp_r_qr)	{ia32.r_qr=ia32_r_qr,	sse2.r_qr=sse2_r_qr,	avx.r_qr=avx_r_qr,	mp=mp_r_qr,	type=SIG_R_QR, simd=true;}
	void set(FP::R_QC ia32_r_qc,	FP::R_QC sse2_r_qc,	FP::R_QC avx_r_qc,	MP::FP::BF mp_r_qc)	{ia32.r_qc=ia32_r_qc,	sse2.r_qc=sse2_r_qc,	avx.r_qc=avx_r_qc,	mp=mp_r_qc,	type=SIG_R_QC, simd=true;}
	void set(FP::R_QQ ia32_r_qq,	FP::R_QQ sse2_r_qq,	FP::R_QQ avx_r_qq,	MP::FP::BF mp_r_qq)	{ia32.r_qq=ia32_r_qq,	sse2.r_qq=sse2_r_qq,	avx.r_qq=avx_r_qq,	mp=mp_r_qq,	type=SIG_R_QQ, simd=true;}

	void set(FP::C_QC ia32_c_qc,	FP::C_QC sse2_c_qc,	FP::C_QC avx_c_qc,	MP::FP::BF mp_c_qc)	{ia32.c_qc=ia32_c_qc,	sse2.c_qc=sse2_c_qc,	avx.c_qc=avx_c_qc,	mp=mp_c_qc,	type=SIG_C_QC, simd=true;}
};
namespace			MP
{
	struct			Instruction//multiprecision instruction
	{
		//type
		// 1	 r_r(Real &r, Real const &x)
		// 2	 c_c(Comp &r, Comp const &x)
		// 3	 q_q(Quat &r, Quat const &x)
		//
		// 4	r_rr(Real &r, Real const &x, Real const &y)
		// 5 *	c_rc(Comp &r, Real const &x, Comp const &y)
		// 6 *	q_rq(Quat &r, Real const &x, Quat const &y)
		// 7	c_cr(Comp &r, Comp const &x, Real const &y)
		// 8	c_cc(Comp &r, Comp const &x, Comp const &y)
		// 9 *	q_cq(Quat &r, Comp const &x, Quat const &y)
		//10	q_qr(Quat &r, Quat const &x, Real const &y)
		//11	q_qc(Quat &r, Quat const &x, Comp const &y)
		//12	q_qq(Quat &r, Quat const &x, Quat const &y)
		//
		//13	 c_r(Comp &r, Real const &x);				//sqrt, ln, polar
		//14 +	 c_q(Comp &r, Quat const &x);				//polar
		//
		//15 +	 r_c(Real &r, Comp const &x);				//re, im, arg
		//16 +	 r_q(Real &r, Quat const &x);				//abs
		//
		//17 +	c_rr(Comp &r, Real const &x, Real const &y);//pow, tetrate
		//
		//18	r_rc(Real &r, Real const &x, Comp const &y);//&& || ##
		//19	r_rq(Real &r, Real const &x, Quat const &y);
		//20	r_cr(Real &r, Comp const &x, Real const &y);
		//21	r_cc(Real &r, Comp const &x, Comp const &y);
		//22	r_cq(Real &r, Comp const &x, Quat const &y);
		//23	r_qr(Real &r, Quat const &x, Real const &y);
		//24	r_qc(Real &r, Quat const &x, Comp const &y);
		//25	r_qq(Real &r, Quat const &x, Quat const &y);
		//26	c_qc(Comp &r, Quat const &x, Comp const &y);
		//
		//27	a ? b : c
		//'c' call						n[result]=ufd[op1]
		//'b' branch if					if(data[op1])i=result
		//'B' branch if not				if(!data[op1])i=result
		//'j' jump						i=result
		//'r' return					data[result]
		char type;
		int result, op1, op2, op3;//data idx
		char r_ms, op1_ms, op2_ms, op3_ms;//math set: R: real, c: complex, h: quaternion
		FP fp;

		std::vector<int> args;

		bool is_binary()const{return type>=4&&type<=12||type>=17&&type<=26;}
		bool has_nonreal_args()const{return op1_ms!='R'||is_binary()&&op2_ms!='R';}

		Instruction(int function, std::vector<int> const &args, int n_result){type='c', op1=function, this->args=args, result=n_result;}
	
		Instruction(char type, int n_condition){this->type=type, op1=n_condition;}//branch:		'b' branch if		'B' branch if not		op1 condition, result dest
		Instruction(){type='j';}//jump		result dest
		Instruction(int n_result){type='r', result=n_result;}//return		result result
	
		//unary function
		Instruction(FPSetter &fps, int op1, char op1_ms, int result, char r_ms):
			fp(fps.mp), type(fps.type), op1(op1), op1_ms(op1_ms), op2(-1), op2_ms(-1), op3(-1), op3_ms(-1), result(result), r_ms(r_ms){}

		//binary function
		Instruction(FPSetter &fps, int op1, char op1_ms, int op2, char op2_ms, int result, char r_ms):
			fp(fps.mp), type(fps.type), op1(op1), op1_ms(op1_ms), op2(op2), op2_ms(op2_ms), op3(-1), op3_ms(-1), result(result), r_ms(r_ms){}

		//inline if
		Instruction(int op1, char op1_ms, int op2, char op2_ms, int op3, char op3_ms, int result, char r_ms):
			fp(0), type(SIG_INLINE_IF), op1(op1), op1_ms(op1_ms), op2(op2), op2_ms(op2_ms), op3(op3), op3_ms(op3_ms), result(result), r_ms(r_ms){}
	};
}
struct				Instruction
{
	// 1	 r_r(VectP &r, VectP const &x)
	// 2	 c_c(CompP &r, CompP const &x)
	// 3	 q_q(QuatP &r, QuatP const &x)
	// 4	r_rr(VectP &r, VectP const &x, VectP const &y)
	// 5 *	c_rc(CompP &r, VectP const &x, CompP const &y)
	// 6 *	q_rq(QuatP &r, VectP const &x, QuatP const &y)
	// 7	c_cr(CompP &r, CompP const &x, VectP const &y)
	// 8	c_cc(CompP &r, CompP const &x, CompP const &y)
	// 9 *	q_cq(QuatP &r, CompP const &x, QuatP const &y)
	//10	q_qr(QuatP &r, QuatP const &x, VectP const &y)
	//11	q_qc(QuatP &r, QuatP const &x, CompP const &y)
	//12	q_qq(QuatP &r, QuatP const &x, QuatP const &y)
	//
	//13	 c_r(CompP &r, VectP const &x);				//sqrt, ln, polar
	//14 +	 c_q(CompP &r, QuatP const &x);				//polar
	//
	//15 +	 r_c(VectP &r, CompP const &x);				//re, im, arg
	//16 +	 r_q(VectP &r, QuatP const &x);				//abs
	//
	//17 +	c_rr(CompP &r, VectP const &x, VectP const &y);//pow, tetrate
	//
	//18	r_rc(VectP &r, VectP const &x, CompP const &y);//&& || ##
	//19	r_rq(VectP &r, VectP const &x, QuatP const &y);
	//20	r_cr(VectP &r, CompP const &x, VectP const &y);
	//21	r_cc(VectP &r, CompP const &x, CompP const &y);
	//22	r_cq(VectP &r, CompP const &x, QuatP const &y);
	//23	r_qr(VectP &r, QuatP const &x, VectP const &y);
	//24	r_qc(VectP &r, QuatP const &x, CompP const &y);
	//25	r_qq(VectP &r, QuatP const &x, QuatP const &y);
	//
	//26	c_qc(CompP &r, QuatP const &x, CompP const &y);
	//
	//27	a ? b : c
	//'c' call						n[result]=ufd[op1]
	//'b' branch if					if(data[op1])i=result
	//'B' branch if not				if(!data[op1])i=result
	//'j' jump						i=result
	//'r' return					data[result]
	char type;
	bool simd;
	int cl_idx, cl_disc_idx;

	int result, op1, op2, op3;
	char r_ms, op1_ms, op2_ms, op3_ms;
	FP ia32, sse2, avx;//in order
	DiscontinuityFunction d;
	
	std::vector<int> args;//arg positions

	bool is_binary()const{return type>=4&&type<=12||type>=17&&type<=26;}
	bool has_nonreal_args()const{return op1_ms!='R'||is_binary()&&op2_ms!='R';}

	//call
	Instruction(int function, std::vector<int> const &args, int n_result):type('c'), op1(function), args(args), result(n_result), simd(false),
		cl_idx(0), cl_disc_idx(0), op2(-1), op3(-1), r_ms(0), op1_ms(0), op2_ms(0), op3_ms(0){}
	
	Instruction(char type, int n_condition):type(type), op1(n_condition), simd(false),//branch:		'b' branch if		'B' branch if not		op1 condition, result dest
		cl_idx(0), cl_disc_idx(0), result(0), op2(-1), op3(-1), r_ms(0), op1_ms(0), op2_ms(0), op3_ms(0){}
	Instruction():type('j'), simd(false),//jump		result dest
		cl_idx(0), cl_disc_idx(0), result(0), op1(-1), op2(-1), op3(-1), r_ms(0), op1_ms(0), op2_ms(0), op3_ms(0){}
	Instruction(int n_result):type('r'), result(n_result), simd(false),//return		result result
		cl_idx(0), cl_disc_idx(0), op1(-1), op2(-1), op3(-1), r_ms(0), op1_ms(0), op2_ms(0), op3_ms(0){}
	
	//unary function
	Instruction(FPSetter &fps, int op1, char op1_ms, int result, char r_ms, DiscontinuityFunction &d, int cl_idx, int cl_disc_idx):
		ia32(fps.ia32), sse2(fps.sse2), avx(fps.avx), type(fps.type), simd(fps.simd), d(d), op1(op1), op1_ms(op1_ms), op2(-1), op2_ms(0), op3(-1), op3_ms(0), result(result), r_ms(r_ms), cl_idx(cl_idx), cl_disc_idx(cl_disc_idx){}

	//binary function
	Instruction(FPSetter &fps, int op1, char op1_ms, int op2, char op2_ms, int result, char r_ms, DiscontinuityFunction &d, int cl_idx, int cl_disc_idx):
		ia32(fps.ia32), sse2(fps.sse2), avx(fps.avx), type(fps.type), simd(fps.simd), d(d), op1(op1), op1_ms(op1_ms), op2(op2), op2_ms(op2_ms), op3(-1), op3_ms(0), result(result), r_ms(r_ms), cl_idx(cl_idx), cl_disc_idx(cl_disc_idx){}

	//inline if
	Instruction(int op1, char op1_ms, int op2, char op2_ms, int op3, char op3_ms, int result, char r_ms, DiscontinuityFunction &d, int cl_idx, int cl_disc_idx):
		d(d), type(SIG_INLINE_IF), simd(true), op1(op1), op1_ms(op1_ms), op2(op2), op2_ms(op2_ms), op3(op3), op3_ms(op3_ms), result(result), r_ms(r_ms), cl_idx(cl_idx), cl_disc_idx(cl_disc_idx){}
};
struct				Variable
{
	std::string name;
	char mathSet,//R real, c complex, h quaternion
		varTypeR, varTypeI, varTypeJ, varTypeK;//x y z space, t time, c constant
	Value val;

	//real variable
	Variable(char const *a, int len, int varTypeR):
		name(a, len), mathSet('R'), varTypeR(varTypeR){}
	//	name(a, len), mathSet('R'), varTypeR(varTypeR){}

	//complex variable
	Variable(char const *a, int len, int varTypeR, int varTypeI):
		name(a, len), mathSet('c'), varTypeR(varTypeR), varTypeI(varTypeI){}
	//	name(a, len), mathSet('C'), varTypeR(varTypeR), varTypeI(varTypeI){}

	//quaternion variable
	Variable(char const *a, int len, int varTypeR, int varTypeI, int varTypeJ, int varTypeK):
		name(a, len), mathSet('h'), varTypeR(varTypeR), varTypeI(varTypeI), varTypeJ(varTypeJ), varTypeK(varTypeK){}

	//user function variable
	Variable(std::string const &name, char mathSet):name(name), mathSet(mathSet), varTypeR(-1), varTypeI(-1), varTypeJ(-1), varTypeK(-1){}
};
struct				UFVariableName
{
	std::string name;
	int scopeLevel, data_idx;
	UFVariableName(char const *a, int len, int scopeLevel, int data_idx):name(a, len), scopeLevel(scopeLevel), data_idx(data_idx){}
};
struct				Term
{
	bool constant;
	char mathSet;//'R' real, 'c' complex, 'h' quaternion	larger value = superset	//'C'==67, 'H'==72, ['R'==82, 'c'==99, 'h'==104], 'r'==114
	int varNo;

	AVector_v4d r, i, j, k;
	void assign(int position, Value const &v, char mathSet)
	{
		r[position]=v.r;
		if(mathSet>='c')
		{
			i[position]=v.i;
			if(mathSet=='h')
				j[position]=v.j, k[position]=v.k;
		}
		//switch(mathSet)
		//{
		//case 'R':r[position]=v.r;break;
		//case 'c':r[position]=v.r, i[position]=v.i;break;
		//case 'h':r[position]=v.r, i[position]=v.i, j[position]=v.j, k[position]=v.k;break;
		//}
	}
	void assign(int position, Value const &v)
	{
		r[position]=v.r;
		if(mathSet>='c')
		{
			i[position]=v.i;
			if(mathSet=='h')
				j[position]=v.j, k[position]=v.k;
		}
		//switch(mathSet)
		//{
		//case 'R':r[position]=v.r;break;
		//case 'c':r[position]=v.r, i[position]=v.i;break;
		//case 'h':r[position]=v.r, i[position]=v.i, j[position]=v.j, k[position]=v.k;break;
		//}
	}
	void assign(int position, double r, double i=0, double j=0, double k=0)
	{
		this->r[position]=r;
		if(mathSet>='c')
		{
			this->i[position]=i;
			if(mathSet=='h')
				this->j[position]=j, this->k[position]=k;
		}
		//switch(mathSet)
		//{
		//case 'R':this->r[position]=r;break;
		//case 'c':this->r[position]=r, this->i[position]=i;break;
		//case 'h':this->r[position]=r, this->i[position]=i, this->j[position]=j, this->k[position]=k;break;
		//}
	}

	//constant
	Term(char mathSet='R'):constant(true), mathSet(mathSet){}

	//function variable constant==false//
	Term(char mathSet, bool constant, int varNo):mathSet(mathSet), constant(constant), varNo(varNo){}

	//expr variable
	Term(char mathSet, int varNo):mathSet(mathSet), constant(false), varNo(varNo){}
};
//struct			DecimalInfo
//{
//	int data_idx, text_start, text_end;
//};
struct				Expression
{
	std::vector<std::pair<int, int>> syntaxErrors;//highlight text[first, second[
	void insertSyntaxError(int first, int second)
	{
		int k=0, kEnd=syntaxErrors.size();
		for(;k<kEnd&&syntaxErrors[k].first<first;++k);
		if(k==kEnd||syntaxErrors[k].first!=first)
			syntaxErrors.insert(syntaxErrors.begin()+k, std::pair<int, int>(first, second));
	}

	std::vector<bool> discontinuities;

	int lineNo, endLineNo, boundNo,
		winColor;//0bgr, weird alpha behavior
private:
	int color,//0rgb, weird alpha behavior
		glColor;//abgr
public:
	int nx, nZ, nQ,
		nISD;//3 space dimentions max
	bool nITD;//1 time dimention max
	char resultMathSet;//'R' real, 'c' complex, 'h' quaternion	larger value = superset	//'C'==67, 'H'==72, ['R'==82, 'c'==99, 'h'==104], 'r'==114
	int resultTerm;

	std::vector<Map> m;
	std::vector<Variable> variables;
	std::vector<Term> n;
#ifdef MULTIPRECISION
	std::vector<MP::Quat> data;//source data
	MP::Quat nresult;//numeric data
	std::vector<int> pi_idx, e_idx;
	typedef std::pair<int, std::string> IdxStr;
	std::vector<IdxStr> decimals;//data idx, source string
#else
	std::vector<Value> data;
#endif
	std::vector<Instruction> i;//graph instructions
	std::vector<MP::Instruction> ni;//numeric instructions
	int resultLogicType;//0: not bool, 1: && ## || < <= > >= logic/inequality, 2: = @ equation, 3: != anti-equation
	int lastInstruction;
	std::vector<int> rmode;

	//user function
	bool valid;
	int name_id;//unique name_id > M_USER_FUNCTION_START
	int nArgs;//, lineNo;
	bool functionStuck;//true: user function returns nan without evaluation when markFunctionsStuck=true

	Expression():color(0), winColor(0), nx(0), nZ(0), nQ(0), nISD(0), nITD(false), rmode(1, 0),
		valid(true), nArgs(0), resultTerm(0){}
	
	void free()
	{
		syntaxErrors.clear();
		discontinuities.clear();
		lineNo=0, endLineNo=0, boundNo=0, color=0, winColor=0, glColor=0, nx=0, nZ=0, nQ=0, nISD=0, nITD=0;
		resultMathSet=0, resultTerm=0;
		m.clear();
		variables.clear();
		//if(n.size()&&n[0].r.size())//
		//{
		//	AVector_v4d test=n[0].r;
		//	_aligned_free(n[0].r.p), n[0].r.p=0;
		//	n[0].r=test;
		//}
		//if(n.size()&&n[0].i.size())//
		//{
		//	AVector_v4d test=n[0].i;
		//	_aligned_free(n[0].i.p), n[0].i.p=0;
		//	n[0].i=test;
		//}
		n.clear();
		data.clear();
	//	ndata.clear();
		pi_idx.clear();
		e_idx.clear();
		decimals.clear();
		i.clear();
		ni.clear();
		resultLogicType=0, lastInstruction=0;
		rmode.assign(1, 0);
		valid=false;
		name_id=0;
		nArgs=0;
		functionStuck=0;
	}
	void insertMap(int pos, int len, G2::Map _0, int _1=0){m.push_back(Map(pos, len, _0, _1));}
//	void insertMap(G2::Map _0, int _1=0){m.push_back(Map(_0, _1));}
	void insertMapData(int pos, int len, MP::Real const &x)
	{
		m.push_back(Map(pos, len, G2::M_N, data.size()));
		n.push_back(Term('R'));
		data.push_back(MP::Quat(x));
	}
	void insertMapPi(int pos, int len)
	{
		m.push_back(Map(pos, len, G2::M_N, data.size()));
		n.push_back(Term('R'));
		pi_idx.push_back(data.size());
		data.push_back(MP::Quat(MP::m_pi));
	}
	void insertMapEuler(int pos, int len)
	{
		m.push_back(Map(pos, len, G2::M_N, data.size()));
		n.push_back(Term('R'));
		e_idx.push_back(data.size());
		data.push_back(MP::Quat(MP::m_e));
	}
	void insertMapData(int pos, int len, char mathSet, double r=0, double i=0, double j=0, double k=0)
	{
		m.push_back(Map(pos, len, G2::M_N, data.size()));
		n.push_back(Term(mathSet));
#ifdef MULTIPRECISION
		data.push_back(MP::Quat(r, i, j, k));
#else
		data.push_back(Value(r, i, j, k));
#endif
	}//*/
/*	void insertMapData(char mathSet, double r=0, double i=0, double j=0, double k=0)
	{
		m.push_back(Map(G2::M_N, data.size()));
		n.push_back(Term(mathSet)), data.push_back(Value(r, i, j, k));
	}//*/
//	void insertData(char mathSet, Value &x){n.push_back(Term(mathSet)), data.push_back(x);}//ambiguous
#ifdef MULTIPRECISION
	void insertData(char mathSet, MP::Quat const &x=MP::Quat(0, 0, 0, 0))
#else
	void insertData(char mathSet, Value x=Value())
#endif
	{
		n.push_back(Term(mathSet));
		data.push_back(x);
	}

//	void insertData(char mathSet, Vect2d const &x, int component){n.push_back(Term(mathSet)), data.push_back(Value(x.get(component)));}
//	void insertData(char mathSet, Comp2d const &c, int component){n.push_back(Term(mathSet)), data.push_back(Value(c.r.get(component), c.i.get(component)));}

	void insertRVar(int, int, char const*, int);
	void insertCVar(int, int, char const*);
	void insertHVar(int, int, char const*);
//	void insertRVar(char const*, int, int);
//	void insertCVar(char const*, int);
//	void insertHVar(char const*, int);
	void setColor_random()
	{
		auto c=(unsigned char*)&color, wc=(unsigned char*)&winColor, gc=(unsigned char*)&glColor;
		gc[0]=c[2]=wc[0]=rand();//red
		gc[1]=c[1]=wc[1]=rand();//green
		gc[2]=c[0]=wc[2]=rand();//blue
		c[3]=wc[3]=0, gc[3]=0xFF;
	}
	void setColor_black(){color=winColor=0, glColor=0xFF000000;}
	int& getColor(){return usingOpenGL?glColor:color;}
	int getColor()const{return usingOpenGL?glColor:color;}
	
	//user function
	void insertFVar(char mathSet, std::string const &name)
	{
		n.push_back(Term(mathSet, false, variables.size()));
#ifdef MULTIPRECISION
		data.push_back(MP::Quat());
#else
		data.push_back(Value());
#endif
		variables.push_back(Variable(name, mathSet));
	}
//	void insertFVarRef(int idx, int varNo){n[idx].varNo=varNo;}
};
#endif