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
extern const unsigned g2_version;
enum 				G2ModeIdx
{
	MODE_N0D,
	MODE_I1D, MODE_N1D,
	MODE_T1D, MODE_T1D_C, MODE_T1D_H,
	MODE_I2D,
	MODE_T2D,
	MODE_C2D,
	MODE_L2D,
	MODE_T2D_H,
	MODE_I3D,
	MODE_C3D,
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
#endif