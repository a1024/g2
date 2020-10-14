//best viewed with tab size of 4 spaces
//g2_error.h - Include for G2 error handling.
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
#ifndef				G2_ERROR_H
#define				G2_ERROR_H
#include			<Windows.h>
extern HWND__		*ghWnd;
extern const int	g_buf_size;//1MB
extern char			g_buf[0x100000];

void				copy_to_clipboard(const char *a, int size);

static const int	e_msg_size=2048;
extern char			first_error_msg[e_msg_size], latest_error_msg[e_msg_size];
void				messageboxa(HWND hWnd, const char *title, const char *format, ...);
void 				log_error(const char *file, int line, const char *format, ...);
#define 			LOGERROR(...)				log_error(__FILE__, __LINE__, __VA_ARGS__)
#define 			LOGERROR_LINE(LINE, ...)	log_error(__FILE__, LINE, __VA_ARGS__)
void				my_assert(int condition, const char *file, int line, const char *msg);
#define				MY_ASSERT(CONDITION, MESSAGE)		my_assert(CONDITION, __FILE__, __LINE__, MESSAGE)

void				sys_check(int line);
#define				SYS_CHECK()		sys_check(__LINE__)

void 				gl_check(int line);
#define				GL_CHECK()		gl_check(__LINE__)
void				gl_error(int line);
#define				GL_ERROR()		gl_error(__LINE__)

void				cl_check(int err, int line);
#define 			CL_CHECK(error)	cl_check(error, __LINE__)
const char*			cl_state2str(int state);
void 				p_check(void *p, int line, const char *func_name);
//#define 			P_CHECK(pointer)	p_check(pointer, __LINE__, #pointer)
#endif