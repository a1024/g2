//best viewed with tab size of 4 spaces
//g2_error.cpp - Implementation of G2 error handling.
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
#include			"g2_error.h"
#include			"g2_cl.h"
#define				CL_TARGET_OPENCL_VERSION 120
#include			<CL/opencl.h>
#include			<GL/gl.h>
#include			<GL/glu.h>
#include			<stdio.h>
#include			<stdarg.h>
#include			<string.h>
#include			<string>

#include			<io.h>//for console
#include			<fcntl.h>
#include			<conio.h>
const int			g_buf_size=0x100000;//1MB
char				g_buf[g_buf_size]={};

void				copy_to_clipboard(const char *a, int size)
{
	char *clipboard=(char*)LocalAlloc(LMEM_FIXED, (size+1)*sizeof(char));
	memcpy(clipboard, a, (size+1)*sizeof(char));
	clipboard[size]='\0';
	OpenClipboard(ghWnd), EmptyClipboard(), SetClipboardData(CF_OEMTEXT, (void*)clipboard), CloseClipboard();
}
void				set_window_title(const char *format, ...)
{
	vsprintf_s(g_buf, g_buf_size, format, (char*)(&format+1));
	int success=SetWindowTextA(ghWnd, g_buf);
	if(!success)
		SYS_CHECK();
}

int					loglevel=LL_PROGRESS;
bool				consoleactive=false;
void				RedirectIOToConsole()//https://stackoverflow.com/questions/191842/how-do-i-get-console-output-in-c-with-a-windows-program
{
	if(!consoleactive)
	{
		consoleactive=true;
		int hConHandle;
		long lStdHandle;
		CONSOLE_SCREEN_BUFFER_INFO coninfo;
		FILE *fp;

		// allocate a console for this app
		AllocConsole();

		// set the screen buffer to be big enough to let us scroll text
		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &coninfo);
		coninfo.dwSize.Y=1000;
		SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), coninfo.dwSize);

		// redirect unbuffered STDOUT to the console
		lStdHandle=(long)GetStdHandle(STD_OUTPUT_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "w");
		*stdout=*fp;
		setvbuf(stdout, nullptr, _IONBF, 0);

		// redirect unbuffered STDIN to the console
		lStdHandle=(long)GetStdHandle(STD_INPUT_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "r");
		*stdin=*fp;
		setvbuf(stdin, nullptr, _IONBF, 0);

		// redirect unbuffered STDERR to the console
		lStdHandle=(long)GetStdHandle(STD_ERROR_HANDLE);
		hConHandle=_open_osfhandle(lStdHandle, _O_TEXT);
		fp=_fdopen(hConHandle, "w");
		*stderr=*fp;
		setvbuf(stderr, nullptr, _IONBF, 0);

		// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog
		// point to console as well
		std::ios::sync_with_stdio();
	}
}
void				freeconsole()
{
	if(consoleactive)
	{
		FreeConsole();
		consoleactive=false;
	}
}
void				print_closewarning(){printf("\n\tWARNING: CLOSING THIS WINDOW WILL CLOSE GRAPHER 2\n\n");}
void				log_start(int priority)
{
	if(loglevel>=priority)
	{
		RedirectIOToConsole();
		print_closewarning();
	}
}
void				log(int priority, const char *format, ...)
{
	if(loglevel>=priority)
		vprintf(format, (char*)(&format+1));
}
void				log_pause(int priority)
{
	if(loglevel>=priority)
		system("pause");
}
#define				log_end freeconsole

char				first_error_msg[e_msg_size]={}, latest_error_msg[e_msg_size]={};
void				messageboxa(HWND hWnd, const char *title, const char *format, ...)
{
	vsprintf_s(g_buf, g_buf_size, format, (char*)(&format+1));
	MessageBoxA(hWnd, g_buf, title, MB_OK);
}
bool 				log_error(const char *file, int line, const char *format, ...)
{
	bool firsttime=first_error_msg[0]=='\0';
	char *buf=first_error_msg[0]?latest_error_msg:first_error_msg;
	va_list args;
	va_start(args, format);
	vsprintf_s(g_buf, e_msg_size, format, args);
	va_end(args);
	int size=strlen(file), start=size-1;
	for(;start>=0&&file[start]!='/'&&file[start]!='\\';--start);
	start+=start==-1||file[start]=='/'||file[start]=='\\';
//	int length=snprintf(buf, e_msg_size, "%s (%d)%s", g_buf, line, file+start);
//	int length=snprintf(buf, e_msg_size, "%s\n%s(%d)", g_buf, file, line);
//	int length=snprintf(buf, e_msg_size, "%s(%d):\n\t%s", file, line, g_buf);
	int length=sprintf_s(buf, e_msg_size, "%s(%d): %s", file+start, line, g_buf);
	if(firsttime)
//	if(buf!=latest_error_msg)
	{
		memcpy(latest_error_msg, first_error_msg, length);
		messageboxa(ghWnd, "Error", latest_error_msg);//redundant, since report_error/emergencyPrint prints both
	}
//	LOGE("%s", latest_error_msg);
	return firsttime;
}

void				my_assert(int condition, const char *file, int line, const char *msg)
{
	if(!condition)
		log_error(file, line, msg);
}

void				sys_check(const char *file, int line)
{
	int error=GetLastError();
	if(error)
	{
		char *messageBuffer=nullptr;
		size_t size=FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER|FORMAT_MESSAGE_FROM_SYSTEM|FORMAT_MESSAGE_IGNORE_INSERTS, NULL, error, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);
	//	std::string message(messageBuffer, size);
		log_error(file, line, "System %d: %s", error, messageBuffer);
		LocalFree(messageBuffer);
	}
}

const char*			glerr2str(int error)
{
#define 			EC(x)	case x:a=(const char*)#x;break
	const char *a=nullptr;
	switch(error)
	{
	case 0:a="SUCCESS";break;
	EC(GL_INVALID_ENUM);
	EC(GL_INVALID_VALUE);
	EC(GL_INVALID_OPERATION);
	case 0x0503:a="GL_STACK_OVERFLOW";break;
	case 0x0504:a="GL_STACK_UNDERFLOW";break;
	EC(GL_OUT_OF_MEMORY);
	case 0x0506:a="GL_INVALID_FRAMEBUFFER_OPERATION";break;
	case 0x0507:a="GL_CONTEXT_LOST";break;
	case 0x8031:a="GL_TABLE_TOO_LARGE";break;
	default:a="???";break;
	}
	return a;
#undef				EC
}
void 				gl_check(const char *file, int line)
{
	int err=glGetError();
	if(err)
		log_error(file, line, "GL %d: %s", err, glerr2str(err));
	//	LOGERROR_LINE(line, "GL Error %d: %s", err, glerr2str(err));
}
void				gl_error(const char *file, int line)
{
	int err=glGetError();
	log_error(file, line, "GL %d: %s", err, glerr2str(err));
}

#define 			EC(x)		case x:a=(const char*)#x;break;
#define 			EC2(n, x)	case n:a=(const char*)#x;break;
const char*			clerr2str(int error)
{
	const char *a=nullptr;
	switch(error)
	{
	EC(CL_SUCCESS)
	EC(CL_DEVICE_NOT_FOUND)
	EC(CL_DEVICE_NOT_AVAILABLE)
	EC(CL_COMPILER_NOT_AVAILABLE)
	EC(CL_MEM_OBJECT_ALLOCATION_FAILURE)
	EC(CL_OUT_OF_RESOURCES)
	EC(CL_OUT_OF_HOST_MEMORY)
	EC(CL_PROFILING_INFO_NOT_AVAILABLE)
	EC(CL_MEM_COPY_OVERLAP)
	EC(CL_IMAGE_FORMAT_MISMATCH)
	EC(CL_IMAGE_FORMAT_NOT_SUPPORTED)
	EC(CL_BUILD_PROGRAM_FAILURE)
	EC(CL_MAP_FAILURE)
//#ifdef CL_VERSION_1_1
	EC(CL_MISALIGNED_SUB_BUFFER_OFFSET)
	EC(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
//#endif
//#ifdef CL_VERSION_1_2
	EC(CL_COMPILE_PROGRAM_FAILURE)
	EC(CL_LINKER_NOT_AVAILABLE)
	EC(CL_LINK_PROGRAM_FAILURE)
	EC(CL_DEVICE_PARTITION_FAILED)
	EC(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
//#endif
	EC(CL_INVALID_VALUE)
	EC(CL_INVALID_DEVICE_TYPE)
	EC(CL_INVALID_PLATFORM)
	EC(CL_INVALID_DEVICE)
	EC(CL_INVALID_CONTEXT)
	EC(CL_INVALID_QUEUE_PROPERTIES)
	EC(CL_INVALID_COMMAND_QUEUE)
	EC(CL_INVALID_HOST_PTR)
	EC(CL_INVALID_MEM_OBJECT)
	EC(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
	EC(CL_INVALID_IMAGE_SIZE)
	EC(CL_INVALID_SAMPLER)
	EC(CL_INVALID_BINARY)
	EC(CL_INVALID_BUILD_OPTIONS)
	EC(CL_INVALID_PROGRAM)
	EC(CL_INVALID_PROGRAM_EXECUTABLE)
	EC(CL_INVALID_KERNEL_NAME)
	EC(CL_INVALID_KERNEL_DEFINITION)
	EC(CL_INVALID_KERNEL)
	EC(CL_INVALID_ARG_INDEX)
	EC(CL_INVALID_ARG_VALUE)
	EC(CL_INVALID_ARG_SIZE)
	EC(CL_INVALID_KERNEL_ARGS)
	EC(CL_INVALID_WORK_DIMENSION)
	EC(CL_INVALID_WORK_GROUP_SIZE)
	EC(CL_INVALID_WORK_ITEM_SIZE)
	EC(CL_INVALID_GLOBAL_OFFSET)
	EC(CL_INVALID_EVENT_WAIT_LIST)
	EC(CL_INVALID_EVENT)
	EC(CL_INVALID_OPERATION)
	EC(CL_INVALID_GL_OBJECT)
	EC(CL_INVALID_BUFFER_SIZE)
	EC(CL_INVALID_MIP_LEVEL)
	EC(CL_INVALID_GLOBAL_WORK_SIZE)
//#ifdef CL_VERSION_1_1
	EC(CL_INVALID_PROPERTY)
//#endif
//#ifdef CL_VERSION_1_2
	EC(CL_INVALID_IMAGE_DESCRIPTOR)
	EC(CL_INVALID_COMPILER_OPTIONS)
	EC(CL_INVALID_LINKER_OPTIONS)
	EC(CL_INVALID_DEVICE_PARTITION_COUNT)
//#endif
//#ifdef CL_VERSION_2_0
	EC2(-69, CL_INVALID_PIPE_SIZE)
	EC2(-70, CL_INVALID_DEVICE_QUEUE)
//#endif
//#ifdef CL_VERSION_2_2
	EC2(-71, CL_INVALID_SPEC_ID)
	EC2(-72, CL_MAX_SIZE_RESTRICTION_EXCEEDED)
//#endif
	EC(CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR)
	EC(CL_PLATFORM_NOT_FOUND_KHR)
	EC2(-1002, CL_INVALID_D3D10_DEVICE_KHR)
	EC2(-1003, CL_INVALID_D3D10_RESOURCE_KHR)
	EC2(-1004, CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR)
	EC2(-1005, CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1006, CL_INVALID_D3D11_DEVICE_KHR)
	EC2(-1007, CL_INVALID_D3D11_RESOURCE_KHR)
	EC2(-1008, CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR)
	EC2(-1009, CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1010, CL_INVALID_D3D9_DEVICE_NV_or_CL_INVALID_DX9_DEVICE_INTEL)
	EC2(-1011, CL_INVALID_D3D9_RESOURCE_NV_or_CL_INVALID_DX9_RESOURCE_INTEL)
	EC2(-1012, CL_D3D9_RESOURCE_ALREADY_ACQUIRED_NV_or_CL_DX9_RESOURCE_ALREADY_ACQUIRED_INTEL)
	EC2(-1013, CL_D3D9_RESOURCE_NOT_ACQUIRED_NV_or_CL_DX9_RESOURCE_NOT_ACQUIRED_INTEL)
	EC2(-1092, CL_EGL_RESOURCE_NOT_ACQUIRED_KHR)
	EC2(-1093, CL_INVALID_EGL_OBJECT_KHR)
	EC2(-1094, CL_INVALID_ACCELERATOR_INTEL)
	EC2(-1095, CL_INVALID_ACCELERATOR_TYPE_INTEL)
	EC2(-1096, CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL)
	EC2(-1097, CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL)
	EC2(-1098, CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL)
	EC2(-1099, CL_INVALID_VA_API_MEDIA_SURFACE_INTEL)
	EC2(-1101, CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL)
	case 1:a="File failure";break;//
	default:
		a="???";
		break;
	}
	return a;
}
const char*			cl_state2str(int state)
{
	const char *a=nullptr;
	switch(state)
	{
	EC(CL_NOTHING)
	EC(CL_LOADING_API)
	EC(CL_API_LOADED)
	EC(CL_CREATING_CONTEXT)
	EC(CL_CONTEXT_CREATED)
	EC(CL_COMPILING_ENTRY)
	EC(CL_COMPILING_PROGRAM00)
	EC(CL_COMPILING_PROGRAM01)
	EC(CL_COMPILING_PROGRAM02)
	EC(CL_COMPILING_PROGRAM03)
	EC(CL_COMPILING_PROGRAM04)
	EC(CL_COMPILING_PROGRAM05)
	EC(CL_COMPILING_PROGRAM06)
	EC(CL_COMPILING_PROGRAM07)
	EC(CL_COMPILING_PROGRAM08)
	EC(CL_COMPILING_PROGRAM09)
	EC(CL_COMPILING_PROGRAM10)
	EC(CL_COMPILING_PROGRAM11)
	EC(CL_COMPILING_PROGRAM12)
	EC(CL_COMPILING_PROGRAM13)
	EC(CL_COMPILING_PROGRAM14)
	EC(CL_COMPILING_PROGRAM15)
	EC(CL_COMPILING_PROGRAM16)
	EC(CL_COMPILING_PROGRAM17)
	EC(CL_COMPILING_PROGRAM18)
	EC(CL_COMPILING_PROGRAM19)
	EC(CL_PROGRAMS_COMPILED)
	EC(CL_RETRIEVING_BINARIES)
	EC(CL_LOADING_PROGRAM00)
	EC(CL_LOADING_PROGRAM01)
	EC(CL_LOADING_PROGRAM02)
	EC(CL_LOADING_PROGRAM03)
	EC(CL_LOADING_PROGRAM04)
	EC(CL_LOADING_PROGRAM05)
	EC(CL_LOADING_PROGRAM06)
	EC(CL_LOADING_PROGRAM07)
	EC(CL_LOADING_PROGRAM08)
	EC(CL_LOADING_PROGRAM09)
	EC(CL_LOADING_PROGRAM10)
	EC(CL_LOADING_PROGRAM11)
	EC(CL_LOADING_PROGRAM12)
	EC(CL_LOADING_PROGRAM13)
	EC(CL_LOADING_PROGRAM14)
	EC(CL_LOADING_PROGRAM15)
	EC(CL_LOADING_PROGRAM16)
	EC(CL_LOADING_PROGRAM17)
	EC(CL_LOADING_PROGRAM18)
	EC(CL_LOADING_PROGRAM19)
	EC(CL_PROGRAMS_LOADED)
	EC(CL_READY_UNTESTED)
	EC(CL_TESTING)
	EC(CL_READY)
	}
	return a;
}
#undef				EC
#undef				EC2
void				cl_check(const char *file, int line, int err)
{
	if(err)
		log_error(file, line, "CL %d: %s", err, clerr2str(err));
	//	LOGERROR_LINE(line, "CL Error %d: %s", err, clerr2str(err));
}
void 				p_check(void *p, const char *file, int line, const char *func_name)
{
	if(!p)
	{
		log_error(file, line, "Error: %s is 0x%08llx", func_name, (long long)p);
	//	LOGE("Line %d error: %s is %lld.", line, func_name, (long long)p);
		OCL_state=CL_NOTHING;
	}
}