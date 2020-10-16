//best viewed with tab size of 4 spaces
//g2_expr.cpp - Implementation of G2 Expression class and its dependencies.
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
#include			"g2_expr.h"

//cl component version 2:
//state.txt:	version in hexadecimal \n saved user-typed expressions
//cl_program[00~19].bin:	OpenCL programs compiled for GPU
const G2Version		g2_version(250, 1);//hi: g2 version, lo: cl component version

//ProfInfo			longest;
std::vector<ProfInfo> prof;
int					prof_array_start_idx=0;
long long			prof_t1=0;
void				prof_start()
{
#ifdef PROFILER
	if(showBenchmark)
	{
		LARGE_INTEGER li;
		QueryPerformanceCounter(&li);
		prof_t1=li.QuadPart;

	//	prof_t1=__rdtsc();
	}
#endif
}
void				prof_add(const char *label, int divisor)
{
#ifdef PROFILER
	if(showBenchmark)
	{
		LARGE_INTEGER li;
		QueryPerformanceCounter(&li);
		long long t2=li.QuadPart;
		QueryPerformanceFrequency(&li);
		prof.push_back(ProfInfo(std::string(label), 1000.*double(t2-prof_t1)/(li.QuadPart*divisor)));
		QueryPerformanceCounter(&li);
		prof_t1=li.QuadPart;

		//long long t2=__rdtsc();
		//prof.push_back(ProfInfo(std::string(label), double(t2-prof_t1)/divisor));
		//prof_t1=__rdtsc();
	}
#endif
}
void				prof_loop_start(const char **labels, int n)
{
#ifdef PROFILER
	if(showBenchmark)
	{
		prof_array_start_idx=prof.size();
		for(int k=0;k<n;++k)
			prof.push_back(ProfInfo(labels[k], 0));
	}
#endif
}
void				prof_add_loop(int idx)
{
#ifdef PROFILER
	if(showBenchmark&&prof_array_start_idx+idx<(int)prof.size())
	{
		LARGE_INTEGER li;
		QueryPerformanceCounter(&li);
		long long t2=li.QuadPart;
		QueryPerformanceFrequency(&li);
		prof[prof_array_start_idx+idx].second+=1000.*double(t2-prof_t1)/(li.QuadPart);
		QueryPerformanceCounter(&li);
		prof_t1=li.QuadPart;
	}
#endif
}
void				prof_print()
{
#ifdef PROFILER
	if(showBenchmark)
	{
		int xpos=w-400, xpos2=w-200;
		int k=0;
		for(int kEnd=prof.size();k<kEnd;++k)
		{
			auto &p=prof[k];
			//if(longest.second<p.second)
			//	longest=p;
			int ypos=k<<4;
		//	int ypos=k*18;
			GUIPrint(xpos, ypos, p.first.c_str());
			GUIPrint(xpos2, ypos, "%lf", p.second);
		//	GUIPrint(xpos2, ypos, "%g", p.second);
		}
		//GUIPrint(xpos, k*18, longest.first.c_str());
		//GUIPrint(xpos2, k*18, "%lf", longest.second);
		//copy to clipboard?
		prof.clear();
		prof_start();
	}
#endif
}

char 				returnMathSet_from_signature(int signature, char op1_ms, char op2_ms, char op3_ms)
{
	switch(signature)
	{
	case SIG_R_R:
	case SIG_R_RR:
	case SIG_R_C:
	case SIG_R_Q:
	case SIG_R_RC:
	case SIG_R_RQ:
	case SIG_R_CR:
	case SIG_R_CC:
	case SIG_R_CQ:
	case SIG_R_QR:
	case SIG_R_QC:
	case SIG_R_QQ:
		return 'R';
	case SIG_C_C:
	case SIG_C_RC:
	case SIG_C_CR:
	case SIG_C_CC:
	case SIG_C_R:
	case SIG_C_Q:
	case SIG_C_RR:
	case SIG_C_QC:
		return 'c';
	case SIG_Q_Q:
	case SIG_Q_RQ:
	case SIG_Q_CQ:
	case SIG_Q_QR:
	case SIG_Q_QC:
	case SIG_Q_QQ:
		return 'q';
	case SIG_INLINE_IF:
		return maximum(op1_ms, op2_ms, op3_ms);
	}
	return 0;
}

void				Expression::insertRVar(int pos, int len, char const *a, int varType)
{
	std::string str(a, len);
	int k=std::find_if(variables.begin(), variables.end(), [&](Variable &v){return v.name==str;})-variables.begin();
	m.push_back(Map(pos, len, G2::M_N, n.size()));
//	n.push_back(Term('r', k));
	n.push_back(Term('R', k));
#ifdef MULTIPRECISION
	data.push_back(MP::Quat());
#else
	data.push_back(Value());
#endif
	if(k==variables.size())
	{
		++nx;
		switch(varType)
		{
		case 's':
			switch(nISD)
			{
			case 0:
			case 1:
			case 2:
				varType='x'+nISD, ++nISD;
				break;
			case 3:
				if(nITD)
					varType='c';
				else
					varType='t', nITD=true;
				break;
			}
			break;
		case 't':
			if(nITD)
				switch(nISD)
				{
				case 0:
				case 1:
				case 2:
					varType='x'+nISD, ++nISD;
					break;
				case 3:
					varType='c';
					break;
				}
			else
				nITD=true;
			break;
		}
		variables.push_back(Variable(a, len, varType));
	}
}
void				Expression::insertCVar(int pos, int len, char const *a)
{
	std::string str(a, len);
	int k=std::find_if(variables.begin(), variables.end(), [&](Variable &v){return v.name==str;})-variables.begin();
	m.push_back(Map(pos, len, G2::M_N, n.size()));
	n.push_back(Term('c', k));
//	n.push_back(Term('C', k));
#ifdef MULTIPRECISION
	data.push_back(MP::Quat());
#else
	data.push_back(Value());
#endif
	if(k==variables.size())
	{
		++nZ;
		switch(nISD)
		{
		case 0:
		case 1:
			variables.push_back(Variable(a, len,
				'x'+nISD,
				'x'+nISD+1)), nISD+=2;
			break;
		case 2:
			variables.push_back(Variable(a, len,
				'z',
				nITD?'c':'t')), ++nISD, nITD=true;
			break;
		case 3:
			variables.push_back(Variable(a, len,
				nITD?'c':'t',
				'c')), nITD=true;
			break;
		}
	}
}
void				Expression::insertHVar(int pos, int len, char const *a)
{
	std::string str(a, len);
	int k=std::find_if(variables.begin(), variables.end(), [&](Variable &v){return v.name==str;})-variables.begin();
	m.push_back(Map(pos, len, G2::M_N, n.size()));
	n.push_back(Term('h', k));
//	n.push_back(Term('C', k));
#ifdef MULTIPRECISION
	data.push_back(MP::Quat());
#else
	data.push_back(Value());
#endif
	if(k==variables.size())
	{
		++nQ;
		switch(nISD)
		{
		case 0:
			variables.push_back(Variable(a, len,
				nITD?'c':'t',
				'x'+nISD,
				'x'+nISD+1,
				'x'+nISD+2)), nISD+=3, nITD=true;
			break;
		case 1:
			variables.push_back(Variable(a, len,
				'x'+nISD,
				'x'+nISD+1,
				nITD?'c':'t', 'c')), nISD+=2, nITD=true;
			break;
		case 2:
			variables.push_back(Variable(a, len,
				'x'+nISD,
				nITD?'c':'t',
				'c',
				'c')), ++nISD, nITD=true;
			break;
		case 3:
			variables.push_back(Variable(a, len,
				nITD?'c':'t',
				'c',
				'c',
				'c')), nITD=true;
			break;
		}
	}
}