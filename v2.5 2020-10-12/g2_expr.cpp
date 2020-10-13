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
const unsigned		g2_version=250<<16|1;//hi: g2 version, lo: cl component version