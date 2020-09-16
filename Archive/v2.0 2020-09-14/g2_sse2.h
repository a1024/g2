//g2_sse2.h - Include for SSE2 version of math functions.
//Copyright (C) 2012-2020  Ayman Wagih Mohsen
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

#ifndef		SSE2_H
#define		SSE2_H
#include	"g2_common.h"
void		minmax_sse2(double *a, int size, double *lo_hi);//{min, max}, even size, 16 byte aligned
namespace	G2
{
	namespace sse2
	{
		void r_r_setzero				(VectP &r, VectP const&);
		void c_c_setzero				(CompP &r, CompP const&);
		void q_q_setzero				(QuatP &r, QuatP const&);

		void  r_r_ceil					(VectP &r, VectP const &x);
		void  r_r_ceil_sse2				(VectP &r, VectP const &x);
		void  c_c_ceil					(CompP &r, CompP const &x);
		void  c_c_ceil_sse2				(CompP &r, CompP const &x);
		void  q_q_ceil					(QuatP &r, QuatP const &x);
		void  q_q_ceil_sse2				(QuatP &r, QuatP const &x);

		void  r_r_floor					(VectP &r, VectP const &x);
		void  r_r_floor_sse2			(VectP &r, VectP const &x);
		void  c_c_floor					(CompP &r, CompP const &x);
		void  c_c_floor_sse2			(CompP &r, CompP const &x);
		void  q_q_floor					(QuatP &r, QuatP const &x);
		void  q_q_floor_sse2			(QuatP &r, QuatP const &x);

		void  r_r_round					(VectP &r, VectP const &x);
		void  r_r_round_sse2			(VectP &r, VectP const &x);
		void  c_c_round					(CompP &r, CompP const &x);
		void  c_c_round_sse2			(CompP &r, CompP const &x);
		void  q_q_round					(QuatP &r, QuatP const &x);
		void  q_q_round_sse2			(QuatP &r, QuatP const &x);

		void  r_r_abs					(VectP &r, VectP const &x);
		void  r_c_abs					(VectP &r, CompP const &x);
		void  r_q_abs					(VectP &r, QuatP const &x);

		void  r_r_arg					(VectP &r, VectP const &x);
		void  r_c_arg					(VectP &r, CompP const &x);
		void  r_q_arg					(VectP &r, QuatP const &x);

		void  r_c_real					(VectP &r, CompP const &x);

		void  r_c_imag					(VectP &r, CompP const &x);

		//r_conjugate: assign
		void c_c_conjugate				(CompP &r, CompP const &x);
		void q_q_conjugate				(QuatP &r, QuatP const &x);

		void  c_r_polar					(CompP &r, VectP const &x);
		void  c_c_polar					(CompP &r, CompP const &x);
		void  c_q_polar					(CompP &r, QuatP const &x);

		//r_cartesian	assign
		void  c_c_cartesian				(CompP &r, CompP const &x);
		void  q_q_cartesian				(QuatP &r, QuatP const &x);

		void r_rr_plus					(VectP &r, VectP const &x, VectP const &y);
		void c_rc_plus					(CompP &r, VectP const &x, CompP const &y);
		void q_rq_plus					(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_plus					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_plus					(CompP &r, CompP const &x, CompP const &y);
		void q_cq_plus					(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_plus					(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_plus					(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_plus					(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_minus					(VectP &r, VectP const &x);
		void  c_c_minus					(CompP &r, CompP const &x);
		void  q_q_minus					(QuatP &r, QuatP const &x);
		void r_rr_minus					(VectP &r, VectP const &x, VectP const &y);
		void c_rc_minus					(CompP &r, VectP const &x, CompP const &y);
		void q_rq_minus					(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_minus					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_minus					(CompP &r, CompP const &x, CompP const &y);
		void q_cq_minus					(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_minus					(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_minus					(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_minus					(QuatP &r, QuatP const &x, QuatP const &y);

		void r_rr_multiply				(VectP &r, VectP const &x, VectP const &y);
		void c_rc_multiply				(CompP &r, VectP const &x, CompP const &y);
		void q_rq_multiply				(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_multiply				(CompP &r, CompP const &x, VectP const &y);
		void c_cc_multiply				(CompP &r, CompP const &x, CompP const &y);//(xr+i*xi)(yr+i*yi) = xr*yr-xi*yi+i(xr*yi+xi*yr)
		void q_cq_multiply				(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_multiply				(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_multiply				(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_multiply				(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_divide				(VectP &r, VectP const &y);
		void  c_c_divide				(CompP &r, CompP const &y);
		void  q_q_divide				(QuatP &r, QuatP const &y);
		void r_rr_divide				(VectP &r, VectP const &x, VectP const &y);
		void c_rc_divide				(CompP &r, VectP const &x, CompP const &y);
		void q_rq_divide				(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_divide				(CompP &r, CompP const &x, VectP const &y);
		void c_cc_divide				(CompP &r, CompP const &x, CompP const &y);
		void q_cq_divide				(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_divide				(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_divide				(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_divide				(QuatP &r, QuatP const &x, QuatP const &y);

		void r_rr_logic_divides			(VectP &r, VectP const &y, VectP const &x);
		void r_rc_logic_divides			(VectP &r, VectP const &y, CompP const &x);
		void r_rq_logic_divides			(VectP &r, VectP const &y, QuatP const &x);
		void r_cr_logic_divides			(VectP &r, CompP const &y, VectP const &x);
		void r_cc_logic_divides			(VectP &r, CompP const &y, CompP const &x);
		void r_cq_logic_divides			(VectP &r, CompP const &y, QuatP const &x);
		void r_qr_logic_divides			(VectP &r, QuatP const &y, VectP const &x);
		void r_qc_logic_divides			(VectP &r, QuatP const &y, CompP const &x);
		void r_qq_logic_divides			(VectP &r, QuatP const &y, QuatP const &x);
		
		void r_rr_logic_divides_sse2	(VectP &r, VectP const &y, VectP const &x);
		void r_rc_logic_divides_sse2	(VectP &r, VectP const &y, CompP const &x);
		void r_rq_logic_divides_sse2	(VectP &r, VectP const &y, QuatP const &x);
		void r_cr_logic_divides_sse2	(VectP &r, CompP const &y, VectP const &x);
		void r_cc_logic_divides_sse2	(VectP &r, CompP const &y, CompP const &x);
		void r_cq_logic_divides_sse2	(VectP &r, CompP const &y, QuatP const &x);
		void r_qr_logic_divides_sse2	(VectP &r, QuatP const &y, VectP const &x);
		void r_qc_logic_divides_sse2	(VectP &r, QuatP const &y, CompP const &x);
		void r_qq_logic_divides_sse2	(VectP &r, QuatP const &y, QuatP const &x);

		void c_cr_pow					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_pow					(CompP &r, CompP const &x, CompP const &y);
		void q_cq_pow					(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_pow					(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_pow					(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_pow					(QuatP &r, QuatP const &x, QuatP const &y);

		void  c_c_ln					(CompP &r, CompP const &x);
		void  q_q_ln					(QuatP &r, QuatP const &x);
	
		void  c_c_log					(CompP &r, CompP const &x);
		void  q_q_log					(QuatP &r, QuatP const &x);
		void c_cr_log					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_log					(CompP &r, CompP const &x, CompP const &y);
		void q_cq_log					(QuatP &r, CompP const &x, QuatP const &y);
		void q_qc_log					(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_log					(QuatP &r, QuatP const &x, QuatP const &y);
	
	/*	Comp2d c_rr_tetrate					(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_tetrate					(Vect2d const &xr, Comp2d const &y);
		Comp2d c_cr_tetrate					(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_tetrate					(Comp2d const &x, Comp2d const &y);
		Quat2d q_qr_tetrate					(Quat2d const &x, Vect2d const &y);
	
		Comp2d c_rr_pentate					(Vect2d const &x, Vect2d const &y);
		Comp2d c_cr_pentate					(Comp2d const &x, Vect2d const &y);//*/

		void  r_r_bitwise_shift_left_l	(VectP &r, VectP const &x);
		void  c_c_bitwise_shift_left_l	(CompP &r, CompP const &x);
		void  q_q_bitwise_shift_left_l	(QuatP &r, QuatP const &x);
		void  r_r_bitwise_shift_left_r	(VectP &r, VectP const &x);
		void  c_c_bitwise_shift_left_r	(CompP &r, CompP const &x);
		void  q_q_bitwise_shift_left_r	(QuatP &r, QuatP const &x);
		void r_rr_bitwise_shift_left	(VectP &r, VectP const &x, VectP const &y);
		void c_rc_bitwise_shift_left	(CompP &r, VectP const &x, CompP const &y);
		void q_rq_bitwise_shift_left	(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_bitwise_shift_left	(CompP &r, CompP const &x, VectP const &y);
		void c_cc_bitwise_shift_left	(CompP &r, CompP const &x, CompP const &y);
		void q_cq_bitwise_shift_left	(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_bitwise_shift_left	(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_bitwise_shift_left	(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_bitwise_shift_left	(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_bitwise_shift_right_l	(VectP &r, VectP const &x);
		void  c_c_bitwise_shift_right_l	(CompP &r, CompP const &x);
		void  q_q_bitwise_shift_right_l	(QuatP &r, QuatP const &x);
		void  r_r_bitwise_shift_right_r	(VectP &r, VectP const &x);
		void  c_c_bitwise_shift_right_r	(CompP &r, CompP const &x);
		void  q_q_bitwise_shift_right_r	(QuatP &r, QuatP const &x);
		void r_rr_bitwise_shift_right	(VectP &r, VectP const &x, VectP const &y);
		void c_rc_bitwise_shift_right	(CompP &r, VectP const &x, CompP const &y);
		void q_rq_bitwise_shift_right	(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_bitwise_shift_right	(CompP &r, CompP const &x, VectP const &y);
		void c_cc_bitwise_shift_right	(CompP &r, CompP const &x, CompP const &y);
		void q_cq_bitwise_shift_right	(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_bitwise_shift_right	(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_bitwise_shift_right	(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_bitwise_shift_right	(QuatP &r, QuatP const &x, QuatP const &y);

	/*	Vect2d  r_r_bitwise_not				(Vect2d const &x);
		Comp2d  c_c_bitwise_not				(Comp2d const &x);
		Quat2d  q_q_bitwise_not				(Quat2d const &x);

		Vect2d  r_r_bitwise_and				(Vect2d const &x);
		Comp2d  c_c_bitwise_and				(Comp2d const &x);
		Quat2d  q_q_bitwise_and				(Quat2d const &x);
		Vect2d r_rr_bitwise_and				(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_bitwise_and				(Vect2d const &x, Comp2d const &y);
		Quat2d q_rq_bitwise_and				(Vect2d const &x, Quat2d const &y);
		Comp2d c_cr_bitwise_and				(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_bitwise_and				(Comp2d const &x, Comp2d const &y);
		Quat2d q_cq_bitwise_and				(Comp2d const &x, Quat2d const &y);
		Quat2d q_qr_bitwise_and				(Quat2d const &x, Vect2d const &y);
		Quat2d q_qc_bitwise_and				(Quat2d const &x, Comp2d const &y);
		Quat2d q_qq_bitwise_and				(Quat2d const &x, Quat2d const &y);

		Vect2d  r_r_bitwise_nand			(Vect2d const &x);
		Comp2d  c_c_bitwise_nand			(Comp2d const &x);
		Quat2d  q_q_bitwise_nand			(Quat2d const &x);
		Vect2d r_rr_bitwise_nand			(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_bitwise_nand			(Vect2d const &x, Comp2d const &y);
		Quat2d q_rq_bitwise_nand			(Vect2d const &x, Quat2d const &y);
		Comp2d c_cr_bitwise_nand			(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_bitwise_nand			(Comp2d const &x, Comp2d const &y);
		Quat2d q_cq_bitwise_nand			(Comp2d const &x, Quat2d const &y);
		Quat2d q_qr_bitwise_nand			(Quat2d const &x, Vect2d const &y);
		Quat2d q_qc_bitwise_nand			(Quat2d const &x, Comp2d const &y);
		Quat2d q_qq_bitwise_nand			(Quat2d const &x, Quat2d const &y);
		
		Vect2d  r_r_bitwise_or				(Vect2d const &x);
		Comp2d  c_c_bitwise_or				(Comp2d const &x);
		Quat2d  q_q_bitwise_or				(Quat2d const &x);
		Vect2d r_rr_bitwise_or				(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_bitwise_or				(Vect2d const &x, Comp2d const &y);
		Quat2d q_rq_bitwise_or				(Vect2d const &x, Quat2d const &y);
		Comp2d c_cr_bitwise_or				(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_bitwise_or				(Comp2d const &x, Comp2d const &y);
		Quat2d q_cq_bitwise_or				(Comp2d const &x, Quat2d const &y);
		Quat2d q_qr_bitwise_or				(Quat2d const &x, Vect2d const &y);
		Quat2d q_qc_bitwise_or				(Quat2d const &x, Comp2d const &y);
		Quat2d q_qq_bitwise_or				(Quat2d const &x, Quat2d const &y);
	
		Vect2d  r_r_bitwise_nor				(Vect2d const &x);
		Comp2d  c_c_bitwise_nor				(Comp2d const &x);
		Quat2d  q_q_bitwise_nor				(Quat2d const &x);
		Vect2d r_rr_bitwise_nor				(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_bitwise_nor				(Vect2d const &x, Comp2d const &y);
		Quat2d q_rq_bitwise_nor				(Vect2d const &x, Quat2d const &y);
		Comp2d c_cr_bitwise_nor				(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_bitwise_nor				(Comp2d const &x, Comp2d const &y);
		Quat2d q_cq_bitwise_nor				(Comp2d const &x, Quat2d const &y);
		Quat2d q_qr_bitwise_nor				(Quat2d const &x, Vect2d const &y);
		Quat2d q_qc_bitwise_nor				(Quat2d const &x, Comp2d const &y);
		Quat2d q_qq_bitwise_nor				(Quat2d const &x, Quat2d const &y);
	
		Vect2d  r_r_bitwise_xor				(Vect2d const &x);
		Comp2d  c_c_bitwise_xor				(Comp2d const &x);
		Quat2d  q_q_bitwise_xor				(Quat2d const &x);
		Vect2d r_rr_bitwise_xor				(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_bitwise_xor				(Vect2d const &x, Comp2d const &y);
		Quat2d q_rq_bitwise_xor				(Vect2d const &x, Quat2d const &y);
		Comp2d c_cr_bitwise_xor				(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_bitwise_xor				(Comp2d const &x, Comp2d const &y);
		Quat2d q_cq_bitwise_xor				(Comp2d const &x, Quat2d const &y);
		Quat2d q_qr_bitwise_xor				(Quat2d const &x, Vect2d const &y);
		Quat2d q_qc_bitwise_xor				(Quat2d const &x, Comp2d const &y);
		Quat2d q_qq_bitwise_xor				(Quat2d const &x, Quat2d const &y);
	
		Vect2d  r_r_bitwise_xnor			(Vect2d const &x);
		Comp2d  c_c_bitwise_xnor			(Comp2d const &x);
		Quat2d  q_q_bitwise_xnor			(Quat2d const &x);
		Vect2d r_rr_bitwise_xnor			(Vect2d const &x, Vect2d const &y);
		Comp2d c_rc_bitwise_xnor			(Vect2d const &x, Comp2d const &y);
		Quat2d q_rq_bitwise_xnor			(Vect2d const &x, Quat2d const &y);
		Comp2d c_cr_bitwise_xnor			(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_bitwise_xnor			(Comp2d const &x, Comp2d const &y);
		Quat2d q_cq_bitwise_xnor			(Comp2d const &x, Quat2d const &y);
		Quat2d q_qr_bitwise_xnor			(Quat2d const &x, Vect2d const &y);
		Quat2d q_qc_bitwise_xnor			(Quat2d const &x, Comp2d const &y);
		Quat2d q_qq_bitwise_xnor			(Quat2d const &x, Quat2d const &y);//*/
	
		void  r_r_logic_equal			(VectP &r, VectP const &x);
		void  r_c_logic_equal			(VectP &r, CompP const &x);
		void  r_q_logic_equal			(VectP &r, QuatP const &x);
		void r_rr_logic_equal			(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_equal			(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_equal			(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_equal			(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_equal			(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_equal			(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_equal			(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_equal			(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_equal			(VectP &r, QuatP const &x, QuatP const &y);
	
		void  r_r_logic_not_equal		(VectP &r, VectP const &x);
		void  r_c_logic_not_equal		(VectP &r, CompP const &x);
		void  r_q_logic_not_equal		(VectP &r, QuatP const &x);
		void r_rr_logic_not_equal		(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_not_equal		(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_not_equal		(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_not_equal		(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_not_equal		(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_not_equal		(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_not_equal		(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_not_equal		(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_not_equal		(VectP &r, QuatP const &x, QuatP const &y);
	
		void  r_r_logic_less_l			(VectP &r, VectP const &x);
		void  r_c_logic_less_l			(VectP &r, CompP const &x);
		void  r_q_logic_less_l			(VectP &r, QuatP const &x);
		void  r_r_logic_less_r			(VectP &r, VectP const &x);
		void  r_c_logic_less_r			(VectP &r, CompP const &x);
		void  r_q_logic_less_r			(VectP &r, QuatP const &x);
		void r_rr_logic_less			(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_less			(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_less			(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_less			(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_less			(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_less			(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_less			(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_less			(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_less			(VectP &r, QuatP const &x, QuatP const &y);
	
		void  r_r_logic_less_equal_l	(VectP &r, VectP const &x);
		void  r_c_logic_less_equal_l	(VectP &r, CompP const &x);
		void  r_q_logic_less_equal_l	(VectP &r, QuatP const &x);
		void  r_r_logic_less_equal_r	(VectP &r, VectP const &x);
		void  r_c_logic_less_equal_r	(VectP &r, CompP const &x);
		void  r_q_logic_less_equal_r	(VectP &r, QuatP const &x);
		void r_rr_logic_less_equal		(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_less_equal		(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_less_equal		(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_less_equal		(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_less_equal		(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_less_equal		(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_less_equal		(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_less_equal		(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_less_equal		(VectP &r, QuatP const &x, QuatP const &y);
	
		void  r_r_logic_greater_l		(VectP &r, VectP const &x);
		void  r_c_logic_greater_l		(VectP &r, CompP const &x);
		void  r_q_logic_greater_l		(VectP &r, QuatP const &x);
		void  r_r_logic_greater_r		(VectP &r, VectP const &x);
		void  r_c_logic_greater_r		(VectP &r, CompP const &x);
		void  r_q_logic_greater_r		(VectP &r, QuatP const &x);
		void r_rr_logic_greater			(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_greater			(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_greater			(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_greater			(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_greater			(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_greater			(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_greater			(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_greater			(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_greater			(VectP &r, QuatP const &x, QuatP const &y);
	
		void  r_r_logic_greater_equal_l	(VectP &r, VectP const &x);
		void  r_c_logic_greater_equal_l	(VectP &r, CompP const &x);
		void  r_q_logic_greater_equal_l	(VectP &r, QuatP const &x);
		void  r_r_logic_greater_equal_r	(VectP &r, VectP const &x);
		void  r_c_logic_greater_equal_r	(VectP &r, CompP const &x);
		void  r_q_logic_greater_equal_r	(VectP &r, QuatP const &x);
		void r_rr_logic_greater_equal	(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_greater_equal	(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_greater_equal	(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_greater_equal	(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_greater_equal	(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_greater_equal	(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_greater_equal	(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_greater_equal	(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_greater_equal	(VectP &r, QuatP const &x, QuatP const &y);
	
		void  r_r_logic_not				(VectP &r, VectP const &x);
		void  r_c_logic_not				(VectP &r, CompP const &x);
		void  r_q_logic_not				(VectP &r, QuatP const &x);
	
		void r_rr_logic_and				(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_and				(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_and				(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_and				(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_and				(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_and				(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_and				(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_and				(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_and				(VectP &r, QuatP const &x, QuatP const &y);
	
		void r_rr_logic_or				(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_or				(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_or				(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_or				(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_or				(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_or				(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_or				(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_or				(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_or				(VectP &r, QuatP const &x, QuatP const &y);
	
		void r_rr_logic_xor				(VectP &r, VectP const &x, VectP const &y);
		void r_rc_logic_xor				(VectP &r, VectP const &x, CompP const &y);
		void r_rq_logic_xor				(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_logic_xor				(VectP &r, CompP const &x, VectP const &y);
		void r_cc_logic_xor				(VectP &r, CompP const &x, CompP const &y);
		void r_cq_logic_xor				(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_logic_xor				(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_logic_xor				(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_logic_xor				(VectP &r, QuatP const &x, QuatP const &y);
	
		void r_rr_condition_zero	(VectP &r, VectP const &x, VectP const &y);
		void c_rc_condition_zero	(CompP &r, VectP const &x, CompP const &y);
		void q_rq_condition_zero	(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_condition_zero	(CompP &r, CompP const &x, VectP const &y);
		void c_cc_condition_zero	(CompP &r, CompP const &x, CompP const &y);
		void q_cq_condition_zero	(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_condition_zero	(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_condition_zero	(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_condition_zero	(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_percent				(VectP &r, VectP const &x);
		void  c_c_percent				(CompP &r, CompP const &x);
		void  q_q_percent				(QuatP &r, QuatP const &x);
	
		void r_rr_modulo				(VectP &r, VectP const &x, VectP const &y);
		void c_rc_modulo				(CompP &r, VectP const &x, CompP const &y);
		void q_rq_modulo				(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_modulo				(CompP &r, CompP const &x, VectP const &y);
		void c_cc_modulo				(CompP &r, CompP const &x, CompP const &y);
		void q_cq_modulo				(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_modulo				(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_modulo				(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_modulo				(QuatP &r, QuatP const &x, QuatP const &y);
	
		void r_rr_modulo_sse2			(VectP &r, VectP const &x, VectP const &y);
		void c_rc_modulo_sse2			(CompP &r, VectP const &x, CompP const &y);
		void q_rq_modulo_sse2			(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_modulo_sse2			(CompP &r, CompP const &x, VectP const &y);
		void c_cc_modulo_sse2			(CompP &r, CompP const &x, CompP const &y);
		void q_cq_modulo_sse2			(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_modulo_sse2			(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_modulo_sse2			(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_modulo_sse2			(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_sgn					(VectP &r, VectP const &x);
		void  c_c_sgn					(CompP &r, CompP const &x);
		void  q_q_sgn					(QuatP &r, QuatP const &x);
		
		void  r_r_sq					(VectP &r, VectP const &x);
		void  c_c_sq					(CompP &r, CompP const &x);
		void  q_q_sq					(QuatP &r, QuatP const &x);

		void  c_c_sqrt					(CompP &r, CompP const &x);
		void  q_q_sqrt					(QuatP &r, QuatP const &x);

		void  r_r_invsqrt				(VectP &r, VectP const &x);

		void  r_r_cbrt					(VectP &r, VectP const &x);
		void  c_c_cbrt					(CompP &r, CompP const &x);
		void  q_q_cbrt					(QuatP &r, QuatP const &x);

		void  r_r_gauss					(VectP &r, VectP const &x);
		void  c_c_gauss					(CompP &r, CompP const &x);
		void  q_q_gauss					(QuatP &r, QuatP const &x);

	/*	void  r_r_erf					(VectP &r, VectP const &x);

		Vect2d  r_r_zeta					(Vect2d const &x);

		Vect2d  r_r_tgamma					(Vect2d const &x);
		Comp2d  c_c_tgamma					(Comp2d const &x);
		Quat2d  q_q_tgamma					(Quat2d const &x);
		Vect2d r_rr_tgamma					(Vect2d const &x, Vect2d const &y);//*/

		void  r_r_loggamma				(VectP &r, VectP const &x);

	/*	Vect2d  r_r_factorial				(Vect2d const &x);
		Comp2d  c_c_factorial				(Comp2d const &x);
		Quat2d  q_q_factorial				(Quat2d const &x);

		Vect2d  r_r_permutation				(Vect2d const &x);
		Comp2d  c_c_permutation				(Comp2d const &x);
		Quat2d  q_q_permutation				(Quat2d const &x);
		Vect2d r_rr_permutation				(Vect2d const &x, Vect2d const &y);
		Comp2d c_cr_permutation				(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_permutation				(Comp2d const &x, Comp2d const &y);
		Quat2d q_qq_permutation				(Quat2d const &x, Quat2d const &y);
	
		Vect2d  r_r_combination				(Vect2d const &x);
		Comp2d  c_c_combination				(Comp2d const &x);
		Quat2d  q_q_combination				(Quat2d const &x);
		Vect2d r_rr_combination				(Vect2d const &x, Vect2d const &y);
		Comp2d c_cr_combination				(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_combination				(Comp2d const &x, Comp2d const &y);
		Quat2d q_qq_combination				(Quat2d const &x, Quat2d const &y);//*/

		void  r_r_cos					(VectP &r, VectP const &x);
		void  c_c_cos					(CompP &r, CompP const &x);
		void  q_q_cos					(QuatP &r, QuatP const &x);

		void  c_c_acos					(CompP &r, CompP const &x);
		void  q_q_acos					(QuatP &r, QuatP const &x);

		void  r_r_cosh					(VectP &r, VectP const &x);
		void  c_c_cosh					(CompP &r, CompP const &x);
		void  q_q_cosh					(QuatP &r, QuatP const &x);

		void  c_c_acosh					(CompP &r, CompP const &x);
		void  q_q_acosh					(QuatP &r, QuatP const &x);

		void  r_r_cosc					(VectP &r, VectP const &x);
		void  c_c_cosc					(CompP &r, CompP const &x);
		void  q_q_cosc					(QuatP &r, QuatP const &x);

		void  r_r_sec					(VectP &r, VectP const &x);
		void  c_c_sec					(CompP &r, CompP const &x);
		void  q_q_sec					(QuatP &r, QuatP const &x);

		void  c_c_asec					(CompP &r, CompP const &x);
		void  q_q_asec					(QuatP &r, QuatP const &x);

		void  r_r_sech					(VectP &r, VectP const &x);
		void  c_c_sech					(CompP &r, CompP const &x);
		void  q_q_sech					(QuatP &r, QuatP const &x);

		void  c_c_asech					(CompP &r, CompP const &x);
		void  q_q_asech					(QuatP &r, QuatP const &x);

		void  r_r_sin					(VectP &r, VectP const &x);
		void  c_c_sin					(CompP &r, CompP const &x);
		void  q_q_sin					(QuatP &r, QuatP const &x);

		void  c_c_asin					(CompP &r, CompP const &x);
		void  q_q_asin					(QuatP &r, QuatP const &x);

		void  r_r_sinh					(VectP &r, VectP const &x);
		void  c_c_sinh					(CompP &r, CompP const &x);
		void  q_q_sinh					(QuatP &r, QuatP const &x);

		void  r_r_asinh					(VectP &r, VectP const &x);
		void  c_c_asinh					(CompP &r, CompP const &x);
		void  q_q_asinh					(QuatP &r, QuatP const &x);

		void  r_r_sinc					(VectP &r, VectP const &x);
		void  c_c_sinc					(CompP &r, CompP const &x);
		void  q_q_sinc					(QuatP &r, QuatP const &x);

		void  r_r_sinhc					(VectP &r, VectP const &x);
		void  c_c_sinhc					(CompP &r, CompP const &x);
		void  q_q_sinhc					(QuatP &r, QuatP const &x);

		void  r_r_csc					(VectP &r, VectP const &x);
		void  c_c_csc					(CompP &r, CompP const &x);
		void  q_q_csc					(QuatP &r, QuatP const &x);

		void  c_c_acsc					(CompP &r, CompP const &x);
		void  q_q_acsc					(QuatP &r, QuatP const &x);

		void  r_r_csch					(VectP &r, VectP const &x);
		void  c_c_csch					(CompP &r, CompP const &x);
		void  q_q_csch					(QuatP &r, QuatP const &x);

		void  r_r_acsch					(VectP &r, VectP const &x);
		void  c_c_acsch					(CompP &r, CompP const &x);
		void  q_q_acsch					(QuatP &r, QuatP const &x);

		void  r_r_tan					(VectP &r, VectP const &x);
		void  c_c_tan					(CompP &r, CompP const &x);
		void  q_q_tan					(QuatP &r, QuatP const &x);

		void  r_r_atan					(VectP &r, VectP const &x);
		void  c_c_atan					(CompP &r, CompP const &x);
		void  q_q_atan					(QuatP &r, QuatP const &x);
		void r_rr_atan					(VectP &r, VectP const &y, VectP const &x);
		void c_rc_atan					(CompP &r, VectP const &y, CompP const &x);
		void q_rq_atan					(QuatP &r, VectP const &y, QuatP const &x);
		void c_cr_atan					(CompP &r, CompP const &y, VectP const &x);
		void c_cc_atan					(CompP &r, CompP const &y, CompP const &x);
		void q_cq_atan					(QuatP &r, CompP const &y, QuatP const &x);
		void q_qr_atan					(QuatP &r, QuatP const &y, VectP const &x);
		void q_qc_atan					(QuatP &r, QuatP const &y, CompP const &x);
		void q_qq_atan					(QuatP &r, QuatP const &y, QuatP const &x);

		void  r_r_tanh					(VectP &r, VectP const &x);
		void  c_c_tanh					(CompP &r, CompP const &x);
		void  q_q_tanh					(QuatP &r, QuatP const &x);

		void  c_c_atanh					(CompP &r, CompP const &x);
		void  q_q_atanh					(QuatP &r, QuatP const &x);

		void  r_r_tanc					(VectP &r, VectP const &x);
		void  c_c_tanc					(CompP &r, CompP const &x);
		void  q_q_tanc					(QuatP &r, QuatP const &x);

		void  r_r_cot					(VectP &r, VectP const &x);
		void  c_c_cot					(CompP &r, CompP const &x);
		void  q_q_cot					(QuatP &r, QuatP const &x);

		void  r_r_acot					(VectP &r, VectP const &x);
		void  c_c_acot					(CompP &r, CompP const &x);
		void  q_q_acot					(QuatP &r, QuatP const &x);

		void  r_r_coth					(VectP &r, VectP const &x);
		void  c_c_coth					(CompP &r, CompP const &x);
		void  q_q_coth					(QuatP &r, QuatP const &x);

		void  c_c_acoth					(CompP &r, CompP const &x);
		void  q_q_acoth					(QuatP &r, QuatP const &x);

		void  r_r_exp					(VectP &r, VectP const &x);
		void  c_c_exp					(CompP &r, CompP const &x);
		void  q_q_exp					(QuatP &r, QuatP const &x);
	
		void  r_r_fib					(VectP &r, VectP const &x);
		void  c_c_fib					(CompP &r, CompP const &x);
		void  q_q_fib					(QuatP &r, QuatP const &x);
	
	/*	Vect2d  r_r_random					(Vect2d const &x);
		Comp2d  c_c_random					(Comp2d const &x);
		Quat2d  q_q_random					(Quat2d const &x);
		Vect2d r_rr_random					(Vect2d const &x, Vect2d const &y);
		Comp2d c_cr_random					(Comp2d const &x, Vect2d const &y);
		Comp2d c_cc_random					(Comp2d const &x, Comp2d const &y);
		Quat2d q_qq_random					(Quat2d const &x, Quat2d const &y);

		Vect2d  r_r_beta					(Vect2d const &x);
		Vect2d r_rr_beta					(Vect2d const &x, Vect2d const &y);
		
		Vect2d  r_r_cyl_bessel_j			(Vect2d const &x);
		Vect2d r_rr_cyl_bessel_j			(Vect2d const &x, Vect2d const &y);

		Vect2d  r_r_cyl_neumann				(Vect2d const &x);
		Vect2d r_rr_cyl_neumann				(Vect2d const &x, Vect2d const &y);

		Comp2d  c_r_hankel1					(Vect2d const &x);
		Comp2d  c_c_hankel1					(Comp2d const &x);
		Comp2d c_rr_hankel1					(Vect2d const &x, Vect2d const &y);//*/

		void  r_r_step					(VectP &r, VectP const &x);
		void  c_c_step					(CompP &r, CompP const &x);
		void  q_q_step					(QuatP &r, QuatP const &x);

		void  r_r_rect					(VectP &r, VectP const &x);
		void  c_c_rect					(CompP &r, CompP const &x);
		void  q_q_rect					(QuatP &r, QuatP const &x);

		void  r_r_trgl					(VectP &r, VectP const &x);
		void  r_c_trgl					(VectP &r, CompP const &x);
		void  r_q_trgl					(VectP &r, QuatP const &x);

		void  r_r_sqwv					(VectP &r, VectP const &x);
		void  r_c_sqwv					(VectP &r, CompP const &x);
		void  r_q_sqwv					(VectP &r, QuatP const &x);
		void r_rr_sqwv					(VectP &r, VectP const &x, VectP const &y);
		void r_rc_sqwv					(VectP &r, VectP const &x, CompP const &y);
		void r_rq_sqwv					(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_sqwv					(VectP &r, CompP const &x, VectP const &y);
		void r_cc_sqwv					(VectP &r, CompP const &x, CompP const &y);
		void r_cq_sqwv					(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_sqwv					(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_sqwv					(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_sqwv					(VectP &r, QuatP const &x, QuatP const &y);
		
		void  r_r_sqwv_sse2				(VectP &r, VectP const &x);
		void  r_c_sqwv_sse2				(VectP &r, CompP const &x);
		void  r_q_sqwv_sse2				(VectP &r, QuatP const &x);
		void r_rr_sqwv_sse2				(VectP &r, VectP const &x, VectP const &y);
		void r_rc_sqwv_sse2				(VectP &r, VectP const &x, CompP const &y);
		void r_rq_sqwv_sse2				(VectP &r, VectP const &x, QuatP const &y);
		void r_cr_sqwv_sse2				(VectP &r, CompP const &x, VectP const &y);
		void r_cc_sqwv_sse2				(VectP &r, CompP const &x, CompP const &y);
		void r_cq_sqwv_sse2				(VectP &r, CompP const &x, QuatP const &y);
		void r_qr_sqwv_sse2				(VectP &r, QuatP const &x, VectP const &y);
		void r_qc_sqwv_sse2				(VectP &r, QuatP const &x, CompP const &y);
		void r_qq_sqwv_sse2				(VectP &r, QuatP const &x, QuatP const &y);

		void  r_r_trwv					(VectP &r, VectP const &x);
		void  r_c_trwv					(VectP &r, CompP const &x);
		void  r_q_trwv					(VectP &r, QuatP const &x);
		void r_rr_trwv					(VectP &r, VectP const &x, VectP const &y);
		void c_cr_trwv					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_trwv					(CompP &r, CompP const &x, CompP const &y);
		void q_qq_trwv					(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_trwv_sse2				(VectP &r, VectP const &x);
		void  r_c_trwv_sse2				(VectP &r, CompP const &x);
		void  r_q_trwv_sse2				(VectP &r, QuatP const &x);
		void r_rr_trwv_sse2				(VectP &r, VectP const &x, VectP const &y);
		void c_cr_trwv_sse2				(CompP &r, CompP const &x, VectP const &y);
		void c_cc_trwv_sse2				(CompP &r, CompP const &x, CompP const &y);
		void q_qq_trwv_sse2				(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_saw					(VectP &r, VectP const &x);
		void  c_c_saw					(CompP &r, CompP const &x);
		void  q_q_saw					(QuatP &r, QuatP const &x);
		void r_rr_saw					(VectP &r, VectP const &x, VectP const &y);
		void c_rc_saw					(CompP &r, VectP const &x, CompP const &y);
		void q_rq_saw					(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_saw					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_saw					(CompP &r, CompP const &x, CompP const &y);
		void q_cq_saw					(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_saw					(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_saw					(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_saw					(QuatP &r, QuatP const &x, QuatP const &y);

		void  r_r_saw_sse2				(VectP &r, VectP const &x);
		void  c_c_saw_sse2				(CompP &r, CompP const &x);
		void  q_q_saw_sse2				(QuatP &r, QuatP const &x);
		void r_rr_saw_sse2				(VectP &r, VectP const &x, VectP const &y);
		void c_rc_saw_sse2				(CompP &r, VectP const &x, CompP const &y);
		void q_rq_saw_sse2				(QuatP &r, VectP const &x, QuatP const &y);
		void c_cr_saw_sse2				(CompP &r, CompP const &x, VectP const &y);
		void c_cc_saw_sse2				(CompP &r, CompP const &x, CompP const &y);
		void q_cq_saw_sse2				(QuatP &r, CompP const &x, QuatP const &y);
		void q_qr_saw_sse2				(QuatP &r, QuatP const &x, VectP const &y);
		void q_qc_saw_sse2				(QuatP &r, QuatP const &x, CompP const &y);
		void q_qq_saw_sse2				(QuatP &r, QuatP const &x, QuatP const &y);

		void r_rr_hypot					(VectP &r, VectP const &x, VectP const &y);
	//	void c_cc_hypot					(CompP &r, CompP const &x, CompP const &y);
	//	void q_qq_hypot					(QuatP &r, QuatP const &x, QuatP const &y);
		
		void r_r_mandelbrot				(VectP &r, VectP const &x);
		void r_c_mandelbrot				(VectP &r, CompP const &x);
		void r_rr_mandelbrot			(VectP &r, VectP const &x, VectP const &y);
		void r_cr_mandelbrot			(VectP &r, CompP const &x, VectP const &y);

		void r_rr_min					(VectP &r, VectP const &x, VectP const &y);
		void c_cr_min					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_min					(CompP &r, CompP const &x, CompP const &y);
		void q_qq_min					(QuatP &r, QuatP const &x, QuatP const &y);

		void r_rr_max					(VectP &r, VectP const &x, VectP const &y);
		void c_cr_max					(CompP &r, CompP const &x, VectP const &y);
		void c_cc_max					(CompP &r, CompP const &x, CompP const &y);
		void q_qq_max					(QuatP &r, QuatP const &x, QuatP const &y);

		void r_rr_conditional_110		(VectP &r, VectP const &x, VectP const &y);
		void c_rc_conditional_110		(CompP &r, VectP const &x, CompP const &y);
		void q_rq_conditional_110		(QuatP &r, VectP const &x, QuatP const &y);
		void r_cr_conditional_110		(VectP &r, CompP const &x, VectP const &y);
		void c_cc_conditional_110		(CompP &r, CompP const &x, CompP const &y);
		void q_cq_conditional_110		(QuatP &r, CompP const &x, QuatP const &y);
		void r_qr_conditional_110		(VectP &r, QuatP const &x, VectP const &y);
		void c_qc_conditional_110		(CompP &r, QuatP const &x, CompP const &y);
		void q_qq_conditional_110		(QuatP &r, QuatP const &x, QuatP const &y);
	
		void r_rr_conditional_101		(VectP &r, VectP const &x, VectP const &y);
		void c_rc_conditional_101		(CompP &r, VectP const &x, CompP const &y);
		void q_rq_conditional_101		(QuatP &r, VectP const &x, QuatP const &y);
		void r_cr_conditional_101		(VectP &r, CompP const &x, VectP const &y);
		void c_cc_conditional_101		(CompP &r, CompP const &x, CompP const &y);
		void q_cq_conditional_101		(QuatP &r, CompP const &x, QuatP const &y);
		void r_qr_conditional_101		(VectP &r, QuatP const &x, VectP const &y);
		void c_qc_conditional_101		(CompP &r, QuatP const &x, CompP const &y);
		void q_qq_conditional_101		(QuatP &r, QuatP const &x, QuatP const &y);
		
		void conditional_111			(QuatP &res, QuatP const &op1, QuatP const &op2, QuatP const &op3, int idx, int op1_ms, int op_ms);

		void  r_r_increment				(VectP &r, VectP const &x);
		void  c_c_increment				(CompP &r, CompP const &x);
		void  q_q_increment				(QuatP &r, QuatP const &x);

		void  r_r_decrement				(VectP &r, VectP const &x);
		void  c_c_decrement				(CompP &r, CompP const &x);
		void  q_q_decrement				(QuatP &r, QuatP const &x);

		void  r_r_assign				(VectP &r, VectP const &x);
		void  c_c_assign				(CompP &r, CompP const &x);
		void  q_q_assign				(QuatP &r, QuatP const &x);
	}//namespace sse2
}//namespace G2
void		lighten_sse2(int *rgb, int n);
namespace	modes
{
	void colorFunction_bcw_sse2(CompP const &v, int *c);
	void colorFunction_bc_l_sse2(CompP const &v, int *c);
//	void colorFunction_sse2(VectP const &r, VectP const &i, int &c1, int &c2);
}
#endif//SSE2_H