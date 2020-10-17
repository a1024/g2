//best viewed with tab size of 4 spaces
//g2_graphics.cpp - Implementation of G2 graphics API.
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

#include		"g2_graphics.h"
#include		"g2_common.h"
#include		"g2_error.h"
#include		"g2_expr.h"//profiler
#include		<Windows.h>
#include		<stdio.h>
#pragma			comment(lib, "OpenGL32.lib")
#pragma			comment(lib, "GLU32.lib")
int				w=0, h=0, X0=0, Y0=0;
char			usingOpenGL=false;
HDC__			*ghMemDC;

HGLRC__			*hRC;
const unsigned char *GLversion;
int				glMajorVer, glMinorVer;
//long long		broken=false;

int				fontH, *fontW;
//const float	_pi=acos(-1.f), _2pi=2*_pi, pi_2=_pi*0.5f, inv_2pi=1/_2pi, sqrt2=sqrt(2.f), torad=_pi/180, infinity=(float)_HUGE, inv255=1.f/255, inv256=1.f/256, inv128=1.f/128;

void		Bitmap::set(int w, int h)
{
	this->w=w, this->h=h;
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0);
}
void		Bitmap::resize(int w, int h)
{
	this->w=w, this->h=h;
	DeleteObject(hBitmap);
	BITMAPINFO bmpInfo={{sizeof(BITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0);
}
void		Bitmap::finish(){DeleteObject(hBitmap);}
void		Bitmap::use	(){hBitmap=(HBITMAP)SelectObject(ghMemDC, hBitmap);}
void		Bitmap::drop(){hBitmap=(HBITMAP)SelectObject(ghMemDC, hBitmap);}

//region
//Region::Region(int x1, int y1, int x2, int y2){CreateRectRgn(x1, y1, x2, y2);}
//Region::~Region(){DeleteObject(hRgn);}
const Region *Region::current=nullptr;
void		Region::create(int x1, int y1, int x2, int y2)
{
	hRgn=CreateRectRgn(x1, y1, x2, y2);
	bx1=x1, bx2=x2, by1=y1, by2=y2, bw=x2-x1, bh=y2-y1, X0=bw>>1, Y0=bh>>1;
}
void		Region::destroy(){DeleteObject(hRgn);}
void		Region::use(){SelectClipRgn(ghMemDC, hRgn);}
void		Region::drop(){SelectClipRgn(ghMemDC, 0);}

void		capture_window(int *rgb){std::copy(gBitmap.rgb, gBitmap.rgb+w*h, rgb);}
void		display_texture_fullwindow(int *rgb){std::copy(rgb, rgb+w*h, gBitmap.rgb);}

//text API:
class		Font
{
	static const int	NFONTS;
	static HFONT		_hFont_[];
	static int			_charW_[][128], _fontH_[];
	char				font;
	void				assign(){hFont=_hFont_[font], fontW=W=_charW_[font], fontH=H=_fontH_[font];}
public:
	HFONT__				*hFont;
	int					*W, H;
	void				createFonts();
	double				change(int);
	double				larger();
	double				smaller();
	char				size(){return font;}
	double				assign(int newFont)
	{
		hFont=_hFont_[newFont], fontW=W=_charW_[newFont], fontH=H=_fontH_[newFont];
		double A=double(H)/_fontH_[font];
		font=newFont;
		return A;
	}
	int					getTextW(const char *a, int i, int f){return short(GetTabbedTextExtentA(ghMemDC, &a[i], f-i, 0, 0));}
	~					Font();
} font;
HFONT Font::_hFont_[]=
{
	CreateFontA(-  7, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "Small Fonts"),
	CreateFontA(-  8, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "Small Fonts"),
	CreateFontA(- 10, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "Small Fonts"),
	CreateFontA(   1, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System"),
	CreateFontA(  29, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System"),
	CreateFontA(  47, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System"),
	CreateFontA(  64, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System"),
	CreateFontA(  94, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System"),
	CreateFontA( 112, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System"),
	CreateFontA( 128, 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH|FF_ROMAN, "System")
};
const int Font::NFONTS=sizeof(_hFont_)/sizeof(HFONT);
int Font::_charW_[NFONTS][128], Font::_fontH_[NFONTS];
void			Font::createFonts()
{
	for(int k=0;k<NFONTS;++k)
	{
		_hFont_[k]=(HFONT__*)SelectObject(ghMemDC, _hFont_[k]);
		for(int k2=0;k2<128;++k2)
		{
			auto LOL=(char*)&k2;
			_charW_[k][k2]=(short)GetTabbedTextExtentA(ghMemDC, (char*)&k2, 1, 0, 0);//
		}
		_fontH_[k]=GetTabbedTextExtentA(ghMemDC, "text", 4, 0, 0)>>16;
		_hFont_[k]=(HFONT__*)SelectObject(ghMemDC, _hFont_[k]);
	}
	font=3;
	assign();
}
double			Font::larger	(){return font<NFONTS-1	?assign(font+1):1;}
double			Font::smaller	(){return font>0		?assign(font-1):1;}
double			Font::change(int wParam)
{
	if(((short*)&wParam)[1]<0)//wheel back
	{
		if(font>0)
			return assign(font-1);
	}
	else//wheel forward
	{
		if(font<NFONTS-1)
			return assign(font+1);
	}
	return 1;
}
				Font::~Font()
{
	for(int k=0;k<NFONTS;++k)
		DeleteObject(_hFont_[k]);
}
//double		changeFont(unsigned wParam){return font.change(wParam);}
//double		largerFont(){return font.larger();}
//double		smallerFont(){return font.smaller();}
//double		setFont(int newFont){return font.assign(newFont);}
//void		selectFont	(){font.hFont=(HFONT__*)SelectObject(ghMemDC, font.hFont);}
//void		deselectFont(){font.hFont=(HFONT__*)SelectObject(ghMemDC, font.hFont);}
//int			getTextWidth(const char *a, int length)
//{
//	return font.getTextW(a, 0, length);
//}
//int			getTextWidth(const char *a, int i, int f)
//{
//	return font.getTextW(a, i, f);
//}

//int			getBkMode(){return GetBkMode(ghMemDC);}
//int			setBkMode(int mode){return SetBkMode(ghMemDC, mode);}
//int			getBkColor(){return GetBkColor(ghMemDC);}
//int			setBkColor(int color){return SetBkColor(ghMemDC, color);}
//int			getTextColor(){return GetTextColor(ghMemDC);}
//int			setTextColor(int color){return SetTextColor(ghMemDC, color);}
//const int	g_buf_size=65536;
////const int	g_buf_size=1024;
//char		g_buf[g_buf_size]={0};
void		emergencyPrint(int x, int y, const char* format, ...)
{
	int linelen=vsprintf_s(g_buf, format, (char*)(&format+1));
	if(linelen>0)
		TextOutA(ghDC, x, y, g_buf, linelen);
}
//void		GUIPrint(int x, int y, const char* format, ...)
//{
//	int linelen=vsprintf_s(g_buf, format, (char*)(&format+1));
//	if(linelen>0)
//		TextOutA(ghMemDC, x, y, g_buf, linelen);
//}
//void		GUIPrint(int x, int y, int value)
//{
//	int linelen=sprintf_s(g_buf, 1024, "%d", value);
//	if(linelen>0)
//		TextOutA(ghMemDC, x, y, g_buf, linelen);
//}
//int			print(int x, int y, int tab_origin, const char *a, int length)
//{
//	return TabbedTextOutA(ghMemDC, x, y, a, length, 0, 0, tab_origin);
//}
//int			print(int x, int y, const char *a, int length)
//{
//	return TextOutA(ghMemDC, x, y, a, length);
//}

////2D API:
//PenBrush::PenBrush(int color)
//{
//	this->color=color;
//	hPen=CreatePen(PS_SOLID, 1, color);
//	hBrush=CreateSolidBrush(color);
//}
//PenBrush::~PenBrush()
//{
//	DeleteObject(hPen);
//	DeleteObject(hBrush);
//}
//void		PenBrush::use()
//{
//	hPen=(HPEN__*)SelectObject(ghMemDC, hPen);
//	hBrush=(HBRUSH__*)SelectObject(ghMemDC, hBrush);
//}
//void		PenBrush::drop()
//{
//	hPen=(HPEN__*)SelectObject(ghMemDC, hPen);
//	hBrush=(HBRUSH__*)SelectObject(ghMemDC, hBrush);
//}

//Pen::Pen(int color){hPen=CreatePen(PS_SOLID, 1, color);}
//Pen::~Pen()
//{
//	if(hPen)
//		DeleteObject(hPen);
//}
//void		Pen::set(int color){hPen=CreatePen(PS_SOLID, 1, color);}
//void		Pen::destroy(){DeleteObject(hPen), hPen=0;}
//void		Pen::use(){hPen=(HPEN__*)SelectObject(ghMemDC, hPen);}
//void		Pen::drop(){hPen=(HPEN__*)SelectObject(ghMemDC, hPen);}

//void		line(int x1, int y1, int x2, int y2){MoveToEx(ghMemDC, x1, y1, 0), LineTo(ghMemDC, x2, y2);}
//void		moveTo(int x, int y){MoveToEx(ghMemDC, x, y, 0);}
//void		lineTo(int x, int y){LineTo(ghMemDC, x, y);}
//void		rectangle(int x1, int y1, int x2, int y2){Rectangle(ghMemDC, x1, y1, x2, y2);}

//3D API:
double			*wbuffer=0;
int				*libuffer, *lfbuffer;
Camera const	*cam;
//void			clear_depth_buffer()
//{
//	memset(wbuffer, 0, w*h*sizeof(double));
////	memset(wbuffer, 0, bw*bh*sizeof(double));
//}
void			resize_3D(int w, int h)
{
	if(!usingOpenGL)
	{
		X0=w>>1, Y0=h>>1;
		wbuffer=(double*)realloc(wbuffer, w*h*sizeof(double));
		libuffer=(int*)realloc(libuffer, h*sizeof(int)), lfbuffer=(int*)realloc(lfbuffer, h*sizeof(int));
	}
}
inline void		line_A_coeff_x(double Xcp1, double Ycp1, double Zcp1, double Xcp2, double Ycp2, double Zcp2, double &a, double &b)
{
	double t=(Xcp2-Xcp1)*Zcp1-Xcp1*(Zcp2-Zcp1);
	a=(Zcp1-Zcp2)*cam->tanfov/(X0*t), b=((Zcp2-Zcp1)*cam->tanfov+Xcp2-Xcp1)/t;
}
inline void		line_A_coeff_y(double Xcp1, double Ycp1, double Zcp1, double Xcp2, double Ycp2, double Zcp2, double &a, double &b)
{
	double t=(Ycp2-Ycp1)*Zcp1-Ycp1*(Zcp2-Zcp1);
	a=(Zcp1-Zcp2)*cam->tanfov/(X0*t), b=((Zcp2-Zcp1)*Y0*cam->tanfov/X0+Ycp2-Ycp1)/t;
}
void			line_3D(Camera const &cam, dvec2 const &s1, dvec3 const &cp1, dvec2 const &s2, dvec3 const &cp2, int lineColor)
//void			line_3D(Camera const &cam, Region const &r, dvec2 const &s1, dvec3 const &cp1, dvec2 const &s2, dvec3 const &cp2, int lineColor)
{
	::cam=&cam;
//	Region::current=&r;
	const double &x1=s1.x, &y1=s1.y, &Xcp1=cp1.x, &Ycp1=cp1.y, &Zcp1=cp1.z, &x2=s2.x, &y2=s2.y, &Xcp2=cp2.x, &Ycp2=cp2.y, &Zcp2=cp2.z;
//	const int &bx1=r.bx1, &bx2=r.bx2, &by1=r.by1, &by2=r.by2, &bw=r.bw, &bh=r.bh, &X0=r.X0, &Y0=r.Y0;
	const int bx1=0, bx2=w, by1=0, by2=h, bw=w, bh=h;
	double xa, ya, xb, yb, dx=x2-x1, dy=y2-y1;
	if(abs(dx)>abs(dy))//horizontal
	{
		double dy_dx=dy/dx;
		if(x1<x2)
		{
			if(x1<bx1)
				xa=bx1, ya=y1+dy_dx*(bx1-x1);
			else
				xa=x1, ya=y1;
			if(x2>bx2-1)
				xb=bx2-1, yb=y1+dy_dx*(bx2-1-x1);
			else
				xb=x2, yb=y2;
		}
		else
		{
			if(x2<bx1)
				xa=bx1, ya=y1+dy_dx*(bx1-x1);
			else
				xa=x2, ya=y2;
			if(x1>bx2-1)
				xb=bx2-1, yb=y1+dy_dx*(bx2-1-x1);
			else
				xb=x1, yb=y1;
		}
		double a, b;
		line_A_coeff_x(Xcp1, Ycp1, Zcp1, Xcp2, Ycp2, Zcp2, a, b);
		//{
		//	double t=(Xcp2-Xcp1)*Zcp1-Xcp1*(Zcp2-Zcp1);
		//	a=(Zcp1-Zcp2)*tanfov/(X0*t), b=((Zcp2-Zcp1)*tanfov+Xcp2-Xcp1)/t;
		//}
		if(SSE4_1)
		{
			int xEnd=int(xb-xa)-(int(xb-xa)&7);
			for(int x=0;x<xEnd;x+=4)
			{
				__m128 xf=_mm_set_ps(float(x+3), float(x+2), float(x+1), (float)x);
				xf=_mm_add_ps(xf, _mm_set1_ps((float)xa));
				__m128 yf=_mm_sub_ps(xf, _mm_set1_ps((float)x1));
				yf=_mm_mul_ps(yf, _mm_set1_ps((float)dy_dx));
				yf=_mm_add_ps(yf, _mm_set1_ps((float)y1));
				__m128i m_y=_mm_cvtps_epi32(yf);
				__m128i c1=_mm_cmpgt_epi32(m_y, _mm_set1_epi32(by1));//
				__m128i c2=_mm_cmplt_epi32(m_y, _mm_set1_epi32(by2));
				c1=_mm_and_si128(c1, c2);

				c2=_mm_srli_si128(c1, 8);
				c2=_mm_or_si128(c2, c1);
				__m128i c3=_mm_srli_si128(c2, 4);
				c2=_mm_or_si128(c2, c3);
				if(c2.m128i_i32[0])
				{
					__m128i m_x=_mm_cvtps_epi32(xf);
					__m128i pos=_mm_mullo_epi32(m_y, _mm_set1_epi32(w));//SSE4.1
					pos=_mm_add_epi32(pos, m_x);
					__m128 A=_mm_mul_ps(xf, _mm_set1_ps((float)a));
					A=_mm_add_ps(A, _mm_set1_ps((float)b));
				//	__m128 wbk=_mm_set_ps(wbuffer+pos.m128i_i32[0], wbuffer+pos.m128i_i32[0]
					if(c1.m128i_i32[0]&&A.m128_f32[0]>wbuffer[pos.m128i_i32[0]])
						gBitmap.rgb[pos.m128i_i32[0]]=lineColor, wbuffer[pos.m128i_i32[0]]=A.m128_f32[0];
					if(c1.m128i_i32[1]&&A.m128_f32[1]>wbuffer[pos.m128i_i32[1]])
						gBitmap.rgb[pos.m128i_i32[1]]=lineColor, wbuffer[pos.m128i_i32[1]]=A.m128_f32[1];
					if(c1.m128i_i32[2]&&A.m128_f32[2]>wbuffer[pos.m128i_i32[2]])
						gBitmap.rgb[pos.m128i_i32[2]]=lineColor, wbuffer[pos.m128i_i32[2]]=A.m128_f32[2];
					if(c1.m128i_i32[3]&&A.m128_f32[3]>wbuffer[pos.m128i_i32[3]])
						gBitmap.rgb[pos.m128i_i32[3]]=lineColor, wbuffer[pos.m128i_i32[3]]=A.m128_f32[3];
				}
			}
			for(int x=xEnd<0?0:xEnd, xEnd2=int(xb-xa);x<=xEnd2;++x)//horizontal
			{
				int xx=x+(int)xa;
				int y=int(std::floor(y1+dy_dx*(xx-x1)));//-0.5 truncated as 0
				if(y>=by1&&y<by2)
				{
					int pos=w*y+xx;
					double A=a*xx+b;
					if(A>wbuffer[pos])
						gBitmap.rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&gBitmap.rgb[pos], c=(unsigned char*)&lineColor;//little endian
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
			}
		}
		else
		{
			for(int x=int(xa), xEnd=int(xb);x<=xEnd;++x)//horizontal
			{
				//int y;
				//{
				//	double Y=y1+dy_dx*(x-x1);
				//	y=int(Y)-(Y<0);
				//}
				int y=int(std::floor(y1+dy_dx*(x-x1)));//-0.5 truncated as 0
				if(y>=by1&&y<by2)
				{
					int pos=w*y+x;
					double A=a*x+b;
					if(A>wbuffer[pos])
						gBitmap.rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&gBitmap.rgb[pos], c=(unsigned char*)&lineColor;//little endian
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
					//{
					//	((unsigned char*)&rgb[pos])[0]=((unsigned char*)&rgb[pos])[0]+((unsigned char*)&lineColor)[0]>>1;
					//	((unsigned char*)&rgb[pos])[1]=((unsigned char*)&rgb[pos])[1]+((unsigned char*)&lineColor)[1]>>1;
					//	((unsigned char*)&rgb[pos])[2]=((unsigned char*)&rgb[pos])[2]+((unsigned char*)&lineColor)[2]>>1;
					//}
					//	rgb[pos]=0xFFC0CB;
				}
			}
		}
	}
	else//vertical
	{
		double dx_dy=dx/dy;
		if(y1<y2)
		{
			if(y1<by1)
				xa=x1+dx_dy*(by1-y1), ya=by1;
			else
				xa=x1, ya=y1;
			if(y2>by2-1)
				xb=x1+dx_dy*(by2-1-y1), yb=by2-1;
			else
				xb=x2, yb=y2;
		}
		else
		{
			if(y2<by1)
				xa=x1+dx_dy*(by1-y1), ya=by1;
			else
				xa=x2, ya=y2;
			if(y1>by2-1)
				xb=x1+dx_dy*(by2-1-y1), yb=by2-1;
			else
				xb=x1, yb=y1;
		}
		double a, b;
		line_A_coeff_y(Xcp1, Ycp1, Zcp1, Xcp2, Ycp2, Zcp2, a, b);
		//{
		//	double t=(Ycp2-Ycp1)*Zcp1-Ycp1*(Zcp2-Zcp1);
		//	a=(Zcp1-Zcp2)*tanfov/(X0*t), b=((Zcp2-Zcp1)*Y0*tanfov/X0+Ycp2-Ycp1)/t;
		//}
		if(SSE4_1)
		{
			int yEnd=int(yb-ya)-(int(yb-ya)&7);
			for(int y=0;y<yEnd;y+=4)
			{
				__m128 yf=_mm_set_ps(float(y+3), float(y+2), float(y+1), (float)y);
				yf=_mm_add_ps(yf, _mm_set1_ps((float)ya));
				__m128 xf=_mm_sub_ps(yf, _mm_set1_ps((float)y1));
				xf=_mm_mul_ps(xf, _mm_set1_ps((float)dx_dy));
				xf=_mm_add_ps(xf, _mm_set1_ps((float)x1));
				__m128i m_x=_mm_cvtps_epi32(xf);
				__m128i c1=_mm_cmpgt_epi32(m_x, _mm_set1_epi32(bx1));//
				__m128i c2=_mm_cmplt_epi32(m_x, _mm_set1_epi32(bx2));
				c1=_mm_and_si128(c1, c2);

				c2=_mm_srli_si128(c1, 8);
				c2=_mm_or_si128(c2, c1);
				__m128i c3=_mm_srli_si128(c2, 4);
				c2=_mm_or_si128(c2, c3);
				if(c2.m128i_i32[0])
				{
					__m128i m_y=_mm_cvtps_epi32(yf);
					__m128i pos=_mm_mullo_epi32(m_y, _mm_set1_epi32(w));
					pos=_mm_add_epi32(pos, m_x);
					__m128 A=_mm_mul_ps(yf, _mm_set1_ps((float)a));
					A=_mm_add_ps(A, _mm_set1_ps((float)b));
					if(c1.m128i_i32[0]&&A.m128_f32[0]>wbuffer[pos.m128i_i32[0]])
						gBitmap.rgb[pos.m128i_i32[0]]=lineColor, wbuffer[pos.m128i_i32[0]]=A.m128_f32[0];
					if(c1.m128i_i32[1]&&A.m128_f32[1]>wbuffer[pos.m128i_i32[1]])
						gBitmap.rgb[pos.m128i_i32[1]]=lineColor, wbuffer[pos.m128i_i32[1]]=A.m128_f32[1];
					if(c1.m128i_i32[2]&&A.m128_f32[2]>wbuffer[pos.m128i_i32[2]])
						gBitmap.rgb[pos.m128i_i32[2]]=lineColor, wbuffer[pos.m128i_i32[2]]=A.m128_f32[2];
					if(c1.m128i_i32[3]&&A.m128_f32[3]>wbuffer[pos.m128i_i32[3]])
						gBitmap.rgb[pos.m128i_i32[3]]=lineColor, wbuffer[pos.m128i_i32[3]]=A.m128_f32[3];
				}
			}
			for(int y=yEnd<0?0:yEnd, yEnd2=int(yb-ya);y<=yEnd2;++y)//vertical
			{
				int yy=y+(int)ya;
				int x=int(std::floor(x1+dx_dy*(yy-y1)));//-0.5 truncated as 0
				if(x>=bx1&&x<bx2)
				{
					int pos=w*yy+x;
					double A=a*yy+b;
					if(A>wbuffer[pos])
						gBitmap.rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&gBitmap.rgb[pos], c=(unsigned char*)&lineColor;
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
			}
		}
		else
		{
			for(int y=int(ya), yEnd=int(yb);y<yEnd;++y)//vertical
			{
				int x=int(std::floor(x1+dx_dy*(y-y1)));//-0.5 truncated as 0
				if(x>=bx1&&x<bx2)
				{
					int pos=w*y+x;
					double A=a*y+b;
					if(A>wbuffer[pos])
						gBitmap.rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&gBitmap.rgb[pos], c=(unsigned char*)&lineColor;
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
			}
		}
	}
}
//void			point_3D(Camera const &cam, dvec3 &p_world, int lineColor)
//{
//	dvec3 cp;
//	cam.world2cam(p_world, cp);
//	const int bx1=0, bx2=w, by1=0, by2=h, bw=w, bh=h;
//	if(cp.z>0)
//	{
//		double Acp=1/cp.z;
//		dvec2 s;
//		cam.cam2screen(cp, s);
//		//double cpt=X0/(cp.z*cam.tanfov);
//		//dvec2 s(X0+cp.x*cpt, Y0+cp.y*cpt);
//		if((s.x>=bx1)&(s.x<bx2)&(s.y>=by1)&(s.y<by2))
//		{
//			int pos=w*int(s.y)+int(s.x);
//			if(Acp>wbuffer[pos])
//				gBitmap.rgb[pos]=lineColor, wbuffer[pos]=Acp;
//		}
//	}
//}
void			point_3D_cam(Camera const &cam, dvec3 &cp, int lineColor)
{
	const int bx1=0, bx2=w, by1=0, by2=h, bw=w, bh=h;
	if(cp.z>0)
	{
		double Acp=1/cp.z;
		dvec2 s;
		cam.cam2screen(cp, s);
		//double cpt=X0/(cp.z*cam.tanfov);
		//dvec2 s(X0+cp.x*cpt, Y0+cp.y*cpt);
		if((s.x>=bx1)&(s.x<bx2)&(s.y>=by1)&(s.y<by2))
		{
			int pos=w*int(s.y)+int(s.x);
			if(Acp>wbuffer[pos])
				gBitmap.rgb[pos]=lineColor, wbuffer[pos]=Acp;
		}
	}
}
void			_2dSet3dPoint(int x, int y, double a, int c)
{
	const int bx1=0, bx2=w, by1=0, by2=h, bw=w, bh=h;
	if((x>=bx1)&(x<bx2)&(y>=by1)&(y<by2))
//	if(x>=bx1&&x<bx2&&y>=by1&&y<by2)
	{
		int pos=int(bw*y+x);
		if(a>wbuffer[pos])
			gBitmap.rgb[pos]=c, wbuffer[pos]=a;
	}
}
//void			point_3D_2x2(Camera const &cam, dvec3 &p_world, int Rcolor)
//{
//	dvec3 cp;
//	cam.world2cam(p_world, cp);
////	double dx=x-camx, dy=y-camy, dz=z-camz, cpt=dx*cax+dy*sax, Xcp1=dx*sax-dy*cax, Ycp1=cpt*say-dz*cay;
////	double Zcp1=cpt*cay+dz*say;
//	if(cp.z>0)
//	{
//		double Acp=1/cp.z;
//		dvec2 s;
//		cam.cam2screen(cp, s);
//		//double cpt=X0/(cp.z*cam.tanfov);
//		//dvec2 s(X0+cp.x*cpt, Y0+cp.y*cpt);
//	//	GUIPrint(ghMemDC, Xs1, Ys1, Rcolor);
//		int xs1=int(s.x)-(s.x<0), ys1=int(s.y)-(s.y<0);
//		_2dSet3dPoint(xs1, ys1	, Acp, Rcolor), _2dSet3dPoint(xs1+1, ys1	, Acp, Rcolor);
//		_2dSet3dPoint(xs1, ys1+1, Acp, Rcolor), _2dSet3dPoint(xs1+1, ys1+1	, Acp, Rcolor);
//	}
//}
//void			point_3D_2x2(Camera const &cam, dvec3 &p_world, int Rcolor, int Icolor, int Jcolor, int Kcolor)
//{
//	dvec3 cp;
//	cam.world2cam(p_world, cp);
//	//double dx=x-camx, dy=y-camy, dz=z-camz, cpt=dx*cax+dy*sax, Xcp1=dx*sax-dy*cax, Ycp1=cpt*say-dz*cay;
//	//double Zcp1=cpt*cay+dz*say;
//	if(cp.z>0)
//	{
//		double Acp=1/cp.z;
//		dvec2 s;
//		cam.cam2screen(cp, s);
//		//double cpt=X0/(cp.z*cam.tanfov);
//		//double Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
//		int xs1=int(s.x)-(s.x<0), ys1=int(s.y)-(s.y<0);
//		_2dSet3dPoint(xs1		, ys1		, Acp, Rcolor);
//		_2dSet3dPoint(xs1+1		, ys1		, Acp, Rcolor);
//		_2dSet3dPoint(xs1+1		, ys1+1		, Acp, Rcolor);
//		_2dSet3dPoint(xs1		, ys1+1		, Acp, Rcolor);
//
//		_2dSet3dPoint(xs1+3		, ys1		, Acp, Icolor);
//		_2dSet3dPoint(xs1+3+1	, ys1		, Acp, Icolor);
//		_2dSet3dPoint(xs1+3+1	, ys1+1		, Acp, Icolor);
//		_2dSet3dPoint(xs1+3		, ys1+1		, Acp, Icolor);
//
//		_2dSet3dPoint(xs1		, ys1+3		, Acp, Jcolor);
//		_2dSet3dPoint(xs1+1		, ys1+3		, Acp, Jcolor);
//		_2dSet3dPoint(xs1+1		, ys1+3+1	, Acp, Jcolor);
//		_2dSet3dPoint(xs1		, ys1+3+1	, Acp, Jcolor);
//
//		_2dSet3dPoint(xs1+3		, ys1+3		, Acp, Kcolor);
//		_2dSet3dPoint(xs1+3+1	, ys1+3		, Acp, Kcolor);
//		_2dSet3dPoint(xs1+3+1	, ys1+3+1	, Acp, Kcolor);
//		_2dSet3dPoint(xs1+3		, ys1+3+1	, Acp, Kcolor);
//	}
//}

bool			lsw_transparency_multiply=false;
void			draft_start		(vec2 &s1, vec2 &s2)
{
	for(int k2=0;k2<h;++k2)libuffer[k2]=w; memset(lfbuffer, 0, h*sizeof(int));
	int k2;
	float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)clamp_positive(s1.y)+1-s1.y);
		//k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);
		if(s1.x<s2.x)
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
					k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	 k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
		else
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	}
	else
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
		}
		else
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	clamp_positive(int(s2.x))
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
		}
	}
	//	 if(s1.y<s2.y){	k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);			 if(s1.x<s2.x)	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	//else				{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
}
void			draft			(vec2 &s1, vec2 &s2)
{
	int k2;float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);
		if(s1.x<s2.x)
		{
			int kEnd=minimum(int(ceil(s2.y)), h);
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
		else
		{
			int kEnd=minimum(int(ceil(s2.y)), h);
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	clamp_positive((int)s2.x)
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
	}
	else
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
		else
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
	}
	//	 if(s1.y<s2.y)	{k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);			 if(s1.x<s2.x)	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	//else				{k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);			 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
}
void			draft_crit_start(vec2 &s1, vec2 &s2)
{
	for(int k2=0;k2<h;++k2)libuffer[k2]=w; memset(lfbuffer, 0, h*sizeof(int));
	int k2;float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
		else
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=lfbuffer[k]=k2, k3+=A;
				//k2=int(k3), k2=k2<0
				//	?	0
				//	:	k2>w
				//		?	minimum(int(s2.x), w)
				//		:	minimum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
	}
	else
	{
		k3=s1.x-A*((long long)s1.y+1);
		if(s1.x<s2.x)
			for(long long k=					0					, kEnd=minimum(int(ceil(s2.y)), h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
		else
			for(long long k=					0					, kEnd=minimum(int(ceil(s2.y)), h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=lfbuffer[k]=k2, k3+=A;
	}
	//	 if(s1.y<s2.y)	{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	//else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)	for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}

	//	 if(s1.y<s2.y)	{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	//else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)	for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
}
void			draft_crit		(vec2 &s1, vec2 &s2)
{
	int k2;float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		else
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
	}
	else
	{
		k3=s1.x-A*((long long)s1.y+1);
		if(s1.x<s2.x)
			for(long long k=					0					, kEnd=minimum(int(s2.y)+1, h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		else
			for(long long k=					0					, kEnd=minimum(int(s2.y)+1, h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
	}
	//	 if(s1.y<s2.y)	{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	//else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)	for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	//if(s1.y<s2.y)
	//{
	//	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
	//	if(s1.x<s2.x)
	//		for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)
	//			k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//	else
	//		for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)
	//			k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//}
	//else
	//{
	//	k3=s1.x-A*((long long)s1.y+1);
	//	if(s1.x<s2.x)for(long long k=			0					+1;k<h&&k<s2.y	;++k)
	//		k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//	else
	//		for(long long k=					0					+1;k<h&&k<s2.y	;++k)
	//			k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//}
//		 if(s1.y<s2.y){	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
//																				else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
//	else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
//																				else				for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
}
void			draft_1behind	(vec2 &s1, vec2 &s2, vec2 &s3)//1 & 2 forward		3 behind
{
	draft_start(s1, s2);
	draft_crit(s3, s1);
	draft_crit(s3, s2);
	if(s3.y>s1.y&&s3.y<s2.y||s3.y>s2.y&&s3.y<s1.y){		 if(s3.x<(s2.x-s1.x)*(s3.y-s1.y)/(s2.y-s1.y)+s1.x)	for(int k=0;k<h;++k)lfbuffer[k]=w;
													else													memset(libuffer, 0, h*sizeof(int));}
}
void			draft_2behind	(vec2 &s1, vec2 &s2, vec2 &s3)//1 forward		2 & 3 behind
{
	draft_crit_start(s2, s1);
	draft_crit(s3, s1);
	if(s1.y>s2.y&&s1.y<s3.y||s1.y>s3.y&&s1.y<s2.y){		 if(s1.x<(s3.x-s2.x)*(s1.y-s2.y)/(s3.y-s2.y)+s2.x)	memset(libuffer, 0, h*sizeof(int));
													else													for(int k=0;k<h;++k)lfbuffer[k]=w;}
}
#if 0
void			render_solid_transparent(Camera const &cam, dvec3 const &p1, dvec3 const &p2, dvec3 const &p3, int color)
//void			render_solid_transparent(Camera const &cam, vec3 const &p1, vec3 const &p2, vec3 const &p3, int color)
{
	vec2 s1, s2, s3;
	vec3 c1, c2, c3;
	float X0_tfov=float(X0/cam.tanfov);
	cam.world2camscreen(p1, c1, s1), cam.world2camscreen(p2, c2, s2), cam.world2camscreen(p3, c3, s3);
//	world_to_screen(p1, s1, c1), world_to_screen(p2, s2, c2), world_to_screen(p3, s3, c3);
			if(c1.z<0){		 if(c2.z<0){		 if(c3.z<0)	return;
											else			draft_2behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_2behind(s2, s3, s1);
											else			draft_1behind(s2, s3, s1);}}
	else{					 if(c2.z<0){		 if(c3.z<0)	draft_2behind(s1, s2, s3);
											else			draft_1behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_1behind(s1, s2, s3);
											else			draft_start(s1, s2), draft(s2, s3), draft(s3, s1);}}
	vec3 upy=c2-c1,		//u12	=<12>
		upx=c3-c1,			//ux3	=<13>
		a=upy.cross(upx);	//abc	=<n>	=<12>x<13>
	float t=a.dot(c1);
	if(!t)return;
	float B7=a.x, B8=a.y, B9=a.z*X0_tfov-a.x*X0-a.y*Y0;	float cpt=t*X0_tfov;
	float A1=B7/cpt, A2=B8/cpt, A3=B9/cpt;

	double *wb_y=wbuffer;
	int lib_y, lfb_y
		, *rgb_y=gBitmap.rgb
		;
	__m128i const sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
	__m128i const sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
	const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
	const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
	__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
	__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
	__m128i m_color=_mm_set1_epi32(color);
//#ifdef _DEBUG
	for(int ys=0;ys<h;++ys)
//#else
//			Concurrency::parallel_for(0, h, [&](int ys)
//#endif
	{
	//	double *wb_y=wbuffer+w*ys;
	//	int *rgb_y=rgb+w*ys;
		lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
		float admittance=A1*lib_y+A2*ys+A3;

		int xs=lib_y;
		__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
		__m128 adm_increment=_mm_set1_ps(A1*4);
		for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)
		{
		//	__m128 wb=_mm_loadu_ps(wb_y+xs);
			__m128 wb=_mm_set_ps((float)wb_y[xs+3], (float)wb_y[xs+2], (float)wb_y[xs+1], (float)wb_y[xs]);
			__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
			__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
			__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
			c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
			c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
			if(c_or.m128i_i32[0])
			{
				c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));

				__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
				__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
				__m128i m_c2=_mm_and_si128(m_color, c_result);
				m_rgb=_mm_and_si128(m_rgb, c_result);
				
				__m128i zero=_mm_setzero_si128();
				__m128i m_c_lo=_mm_unpacklo_epi8(m_c2, zero);//{0, tx[7], ..., 0, tx[0]}
				__m128i m_c_hi=_mm_unpackhi_epi8(m_c2, zero);//{0, tx[15], ..., 0, tx[8]}
				__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
				__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
				if(lsw_transparency_multiply)
				{
					m_c_lo=_mm_mullo_epi16(m_c_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
					m_c_hi=_mm_mullo_epi16(m_c_hi, m_rgb_hi);
					m_c_lo=_mm_shuffle_epi8(m_c_lo, sh_mul_lo);
					m_c_hi=_mm_shuffle_epi8(m_c_hi, sh_mul_hi);
				}
				else
				{
					m_c_lo=_mm_add_epi16(m_c_lo, m_rgb_lo);
					m_c_hi=_mm_add_epi16(m_c_hi, m_rgb_hi);
					m_c_lo=_mm_srli_epi16(m_c_lo, 1);
					m_c_hi=_mm_srli_epi16(m_c_hi, 1);
					m_c_lo=_mm_shuffle_epi8(m_c_lo, sh_add_lo);
					m_c_hi=_mm_shuffle_epi8(m_c_hi, sh_add_hi);
				}
				m_c_lo=_mm_and_si128(m_c_lo, mask_lo);
				m_c_hi=_mm_and_si128(m_c_hi, mask_hi);
				m_c2=_mm_or_si128(m_c_lo, m_c_hi);

				m_c2=_mm_or_si128(m_c2, m_rgb2);
				_mm_storeu_si128((__m128i*)(rgb_y+xs), m_c2);
			}
			adm=_mm_add_ps(adm, adm_increment);
		}
		admittance=adm.m128_f32[0];
		for(;xs<lfb_y;++xs)
		{
			if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
				admittance>=wb_y[xs])
			{
				auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&color;
				if(lsw_transparency_multiply)
				{
					p[0]=p[0]*c[0]>>8;//b
					p[1]=p[1]*c[1]>>8;//g
					p[2]=p[2]*c[2]>>8;//r
				}
				else
				{
					p[0]=(p[0]+c[0])>>1;//b
					p[1]=(p[1]+c[1])>>1;//g
					p[2]=(p[2]+c[2])>>1;//r
				}
			}
			admittance+=A1;
		}
		wb_y=wb_y+w, rgb_y=rgb_y+w;
	}
//#ifndef _DEBUG
//			);
//#endif
}
#endif
void			render_textured_transparent(Camera const &cam, dvec3 const &p1, dvec3 const &p2, dvec3 const &p3, int *texture, int txw, int txh, vec2 const &tx1, mat2 const &txm)
//void			render_textured_transparent(Camera const &cam, vec3 const &p1, vec3 const &p2, vec3 const &p3, int *texture, int txw, int txh, vec2 const &tx1, mat2 const &txm)
{
	vec2 s1, s2, s3;
	vec3 c1, c2, c3;
	float X0_tfov=float(X0/cam.tanfov);//
	cam.world2camscreen(p1, c1, s1), cam.world2camscreen(p2, c2, s2), cam.world2camscreen(p3, c3, s3);
//	world_to_screen(p1, s1, c1), world_to_screen(p2, s2, c2), world_to_screen(p3, s3, c3);
			if(c1.z<0){		 if(c2.z<0){		 if(c3.z<0)	return;
											else			draft_2behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_2behind(s2, s3, s1);
											else			draft_1behind(s2, s3, s1);}}
	else{					 if(c2.z<0){		 if(c3.z<0)	draft_2behind(s1, s2, s3);
											else			draft_1behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_1behind(s1, s2, s3);
											else			draft_start(s1, s2), draft(s2, s3), draft(s3, s1);}}
	vec3 upy=c2-c1,		//u12	=<12>
		upx=c3-c1,			//ux3	=<13>
		a=upy.cross(upx);	//abc	=<n>	=<12>x<13>
	float t=a.dot(c1);
	if(!t)return;
	float B7=a.x, B8=a.y, B9=a.z*X0_tfov-a.x*X0-a.y*Y0;	float cpt=t*X0_tfov;
	float A1=B7/cpt, A2=B8/cpt, A3=B9/cpt;

	//textured transparent triangle
	float	ux3_1=upy.magnitude();	upy/=ux3_1;			//u12	=<u12>	=<12>/|12|
			ux3_1=upx.dot(upy),		upx-=ux3_1*upy;		//ux3	=<x3>	=<13>-<13>.<u12><u12>
			ux3_1=upx.magnitude(),	upx/=ux3_1;			//ux3	=<ux3>	=<x3>/|x3|
			ux3_1=upx.dot(c1);							//ux3_1			=<ux3>.<1>
	float	u12_1=upy.dot(c1);							//u12_1			=<u12>.<1>
	vec2 C1=txm*vec2(upx.x, upy.x), C2=txm*vec2(upx.y, upy.y), C3=txm*vec2(upx.z, upy.z);
	float	D1=C1.x*t, D2=C2.x*t, D3=(C3.x*X0_tfov-C1.x*X0-C2.x*Y0)*t;
	float	D4=C1.y*t, D5=C2.y*t, D6=(C3.y*X0_tfov-C1.y*X0-C2.y*Y0)*t;
	vec2 E1=tx1-txm*vec2(ux3_1, u12_1);
	float	B1=D1+E1.x*B7, B2=D2+E1.x*B8, B3=D3+E1.x*B9;
	float	B4=D4+E1.y*B7, B5=D5+E1.y*B8, B6=D6+E1.y*B9;
	double *wb_y=wbuffer;
	int lib_y, lfb_y, *rgb_y=gBitmap.rgb, *txk=texture, txhk=txh, txwk=txw;
	if(SSE4_1)
	{
		__m128i const sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
		__m128i const sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
		const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
		const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
		__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
		__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
		for(int ys=0;ys<h;++ys)
		{
			lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
			float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

			int xs=lib_y;
			__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
			__m128 adm_increment=_mm_set1_ps(A1*4);
			__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
			__m128 m_a_increment=_mm_set1_ps(B1*4);
			__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
			__m128 m_b_increment=_mm_set1_ps(B4*4);
			__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
			__m128 m_c_increment=_mm_set1_ps(B7*4);

			__m128 s_txwk=_mm_set1_ps((float)txwk);
			__m128 s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk);
			__m128i i_txhk=_mm_set1_epi32(txhk);
			for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
			{
			//	__m128 wb=_mm_loadu_ps(wb_y+xs);
				__m128 wb=_mm_set_ps((float)wb_y[xs+3], (float)wb_y[xs+2], (float)wb_y[xs+1], (float)wb_y[xs]);
				__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
				__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
				c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				c_or=_mm_or_si128(c_or, c_result2);			//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(c_or.m128i_i32[0])
				{
					c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
					__m128 a_c=_mm_div_ps(m_a, m_c);
					__m128 b_c=_mm_div_ps(m_b, m_c);

					//mod operation		x%q = x-q*floor(x/q)
					__m128 ac_w=_mm_div_ps(a_c, s_txwk);
					__m128 bc_h=_mm_div_ps(b_c, s_txhk);

					__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
					__m128 f_bch=_mm_floor_ps(bc_h);

					f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
					f_bch=_mm_mul_ps(f_bch, s_txhk);

					a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
					b_c=_mm_sub_ps(b_c, f_bch);
				
					a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
					b_c=_mm_floor_ps(b_c);

					__m128i xtx=_mm_cvtps_epi32(a_c);
					__m128i ytx=_mm_cvtps_epi32(b_c);

					__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
					tx_idx=_mm_add_epi32(tx_idx, xtx);

					__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
					__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
					__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
					m_tx=_mm_and_si128(m_tx, c_result);
					m_rgb=_mm_and_si128(m_rgb, c_result);

					const __m128i zero=_mm_setzero_si128();
					__m128i m_tx_lo=_mm_unpacklo_epi8(m_tx, zero);//{0, tx[7], ..., 0, tx[0]}
					__m128i m_tx_hi=_mm_unpackhi_epi8(m_tx, zero);//{0, tx[15], ..., 0, tx[8]}
					__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
					__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
					if(lsw_transparency_multiply)
					{
						m_tx_lo=_mm_mullo_epi16(m_tx_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
						m_tx_hi=_mm_mullo_epi16(m_tx_hi, m_rgb_hi);
						m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_mul_lo);
						m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_mul_hi);
					}
					else
					{
						m_tx_lo=_mm_add_epi16(m_tx_lo, m_rgb_lo);
						m_tx_hi=_mm_add_epi16(m_tx_hi, m_rgb_hi);
						m_tx_lo=_mm_srli_epi16(m_tx_lo, 1);
						m_tx_hi=_mm_srli_epi16(m_tx_hi, 1);
						m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_add_lo);
						m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_add_hi);
					}
					m_tx_lo=_mm_and_si128(m_tx_lo, mask_lo);
					m_tx_hi=_mm_and_si128(m_tx_hi, mask_hi);
					m_tx=_mm_or_si128(m_tx_lo, m_tx_hi);

					m_tx=_mm_or_si128(m_tx, m_rgb2);
					_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
				}
				adm=_mm_add_ps(adm, adm_increment);
				m_a=_mm_add_ps(m_a, m_a_increment);
				m_b=_mm_add_ps(m_b, m_b_increment);
				m_c=_mm_add_ps(m_c, m_c_increment);
			}
			admittance=adm.m128_f32[0];
			a=m_a.m128_f32[0];
			b=m_b.m128_f32[0];
			c=m_c.m128_f32[0];
			for(;xs<lfb_y;++xs)
			{
				if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
					admittance>wb_y[xs])
				{
					int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
					int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
					auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&txk[Xtx+Ytx*txwk];
					if(lsw_transparency_multiply)
					{
						p[0]=p[0]*c[0]>>8;//b
						p[1]=p[1]*c[1]>>8;//g
						p[2]=p[2]*c[2]>>8;//r
					}
					else
					{
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
				admittance+=A1, a+=B1, b+=B4, c+=B7;
			}
			wb_y=wb_y+w, rgb_y=rgb_y+w;
		}
	}
	else//textured transparent triangle - no SSE4.1
	{
		const __m128 half=_mm_set1_ps(0.5f);
		const __m128i mask2_lo=_mm_set_epi32(0, -1, 0, -1), mask2_hi=_mm_set_epi32(-1, 0, -1, 0);
		const __m128i sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
		const __m128i sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
		const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
		const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
		const __m128i mask_lo=_mm_set_epi32(0, 0, -1, -1);
		const __m128i mask_hi=_mm_set_epi32(-1, -1, 0, 0);

		__m128 adm_increment=_mm_set1_ps(A1*4);
		__m128 m_a_increment=_mm_set1_ps(B1*4);
		__m128 m_b_increment=_mm_set1_ps(B4*4);
		__m128 m_c_increment=_mm_set1_ps(B7*4);

		__m128 s_txwk=_mm_set1_ps((float)txwk);
		__m128 s_txhk=_mm_set1_ps((float)txhk);
		__m128i i_txwk=_mm_set1_epi32(txwk), i_txwk_1=_mm_set1_epi32(txwk-1);
		__m128i i_txhk=_mm_set1_epi32(txhk), i_txhk_1=_mm_set1_epi32(txhk-1);
		for(int ys=0;ys<h;++ys)
		{
			lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
			float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

			int xs=lib_y;
			__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
			__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
			__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
			__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
			for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
			{
			//	__m128 wb=_mm_loadu_ps(wb_y+xs);
				__m128 wb=_mm_set_ps((float)wb_y[xs+3], (float)wb_y[xs+2], (float)wb_y[xs+1], (float)wb_y[xs]);
				__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
				__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
				c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				c_or=_mm_or_si128(c_or, c_result2);				//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(c_or.m128i_i32[0])
				{
					c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
					__m128 a_c=_mm_div_ps(m_a, m_c);
					__m128 b_c=_mm_div_ps(m_b, m_c);

					//mod operation		x%q = x-q*floor(x/q)
					__m128 ac_w=_mm_div_ps(a_c, s_txwk);
					__m128 bc_h=_mm_div_ps(b_c, s_txhk);

					//__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
					//__m128 f_bch=_mm_floor_ps(bc_h);
					ac_w=_mm_sub_ps(ac_w, half);//floor(x/q)
					bc_h=_mm_sub_ps(bc_h, half);
					__m128i i_acw=_mm_cvtps_epi32(ac_w);
					__m128i i_bch=_mm_cvtps_epi32(bc_h);
					__m128 f_acw=_mm_cvtepi32_ps(i_acw);
					__m128 f_bch=_mm_cvtepi32_ps(i_bch);

					f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
					f_bch=_mm_mul_ps(f_bch, s_txhk);

					a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
					b_c=_mm_sub_ps(b_c, f_bch);
				
					//a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
					//b_c=_mm_floor_ps(b_c);
					//__m128i xtx=_mm_cvtps_epi32(a_c);
					//__m128i ytx=_mm_cvtps_epi32(b_c);
					a_c=_mm_sub_ps(a_c, half);
					b_c=_mm_sub_ps(b_c, half);
					__m128i xtx=_mm_cvtps_epi32(a_c);
					__m128i ytx=_mm_cvtps_epi32(b_c);

					__m128i cmp_x=_mm_cmplt_epi32(xtx, _mm_setzero_si128());//index guard
					__m128i cmp_y=_mm_cmplt_epi32(ytx, _mm_setzero_si128());
					cmp_x=_mm_and_si128(cmp_x, i_txwk);
					cmp_y=_mm_and_si128(cmp_y, i_txhk);
					xtx=_mm_add_epi32(xtx, cmp_x);
					ytx=_mm_add_epi32(ytx, cmp_y);
					cmp_x=_mm_cmpgt_epi32(xtx, i_txwk_1);
					cmp_y=_mm_cmpgt_epi32(ytx, i_txhk_1);
					cmp_x=_mm_and_si128(cmp_x, i_txwk);
					cmp_y=_mm_and_si128(cmp_y, i_txhk);
					xtx=_mm_sub_epi32(xtx, cmp_x);
					ytx=_mm_sub_epi32(ytx, cmp_y);

					//__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
					__m128i tx_idx_lo=_mm_mul_epu32(ytx, i_txwk);
					ytx=_mm_shuffle_epi32(ytx, _MM_SHUFFLE(2, 3, 0, 1));//y2 y3 y0 y1
					__m128i tx_idx_hi=_mm_mul_epu32(ytx, i_txwk);
					tx_idx_hi=_mm_shuffle_epi32(tx_idx_hi, _MM_SHUFFLE(2, 3, 0, 1));
					tx_idx_lo=_mm_and_si128(tx_idx_lo, mask2_lo);
					tx_idx_hi=_mm_and_si128(tx_idx_hi, mask2_hi);
					__m128i tx_idx=_mm_or_si128(tx_idx_lo, tx_idx_hi);

					tx_idx=_mm_add_epi32(tx_idx, xtx);

					__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
					__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
					__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
					m_tx=_mm_and_si128(m_tx, c_result);
					m_rgb=_mm_and_si128(m_rgb, c_result);
				
					__m128i const zero=_mm_setzero_si128();
					__m128i m_tx_lo=_mm_unpacklo_epi8(m_tx, zero);//{0, tx[7], ..., 0, tx[0]}
					__m128i m_tx_hi=_mm_unpackhi_epi8(m_tx, zero);//{0, tx[15], ..., 0, tx[8]}
					__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
					__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
					if(lsw_transparency_multiply)
					{
						m_tx_lo=_mm_mullo_epi16(m_tx_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
						m_tx_hi=_mm_mullo_epi16(m_tx_hi, m_rgb_hi);
						m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_mul_lo);
						m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_mul_hi);
						m_tx_lo=_mm_and_si128(m_tx_lo, mask_lo);
						m_tx_hi=_mm_and_si128(m_tx_hi, mask_hi);
					}
					else
					{
						m_tx_lo=_mm_add_epi16(m_tx_lo, m_rgb_lo);
						m_tx_hi=_mm_add_epi16(m_tx_hi, m_rgb_hi);
						m_tx_lo=_mm_srli_epi16(m_tx_lo, 1);
						m_tx_hi=_mm_srli_epi16(m_tx_hi, 1);
						m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_add_lo);
						m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_add_hi);
					}
					m_tx=_mm_or_si128(m_tx_lo, m_tx_hi);

					m_tx=_mm_or_si128(m_tx, m_rgb2);
					_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
				}
				adm=_mm_add_ps(adm, adm_increment);
				m_a=_mm_add_ps(m_a, m_a_increment);
				m_b=_mm_add_ps(m_b, m_b_increment);
				m_c=_mm_add_ps(m_c, m_c_increment);
			}
			admittance=adm.m128_f32[0];
			a=m_a.m128_f32[0];
			b=m_b.m128_f32[0];
			c=m_c.m128_f32[0];
			for(;xs<lfb_y;++xs)
			{
				if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
					admittance>wb_y[xs])
				{
					int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
					int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
					auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&txk[Xtx+Ytx*txwk];
					if(lsw_transparency_multiply)
					{
						p[0]=p[0]*c[0]>>8;//b
						p[1]=p[1]*c[1]>>8;//g
						p[2]=p[2]*c[2]>>8;//r
					}
					else
					{
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
				admittance+=A1, a+=B1, b+=B4, c+=B7;
			}
			wb_y=wb_y+w, rgb_y=rgb_y+w;
		}
	/*	for(int ys=0;ys<h;++ys)
		{
			lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
			float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;
			for(int xs=lib_y;xs<lfb_y;++xs)
			{
				if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
					admittance>wb_y[xs])
				{
					int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
					int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
					auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&txk[Xtx+Ytx*txwk];
					p[0]=p[0]*c[0]>>8;//b
					p[1]=p[1]*c[1]>>8;//g
					p[2]=p[2]*c[2]>>8;//r
				}
				admittance+=A1, a+=B1, b+=B4, c+=B7;
			}
			wb_y=wb_y+w, rgb_y=rgb_y+w;
		}//*/
	}
}

bool		fonts_not_created=true;
void		initiate()
{
	ghMemDC=CreateCompatibleDC(ghDC);
	if(fonts_not_created)
		font.createFonts(), fonts_not_created=false;
	resize_3D(w, h);
}
void		finish(){DeleteDC(ghMemDC);}
//void		show(){BitBlt(ghDC, 0, 0, w, h, ghMemDC, 0, 0, SRCCOPY);}

//
//OpenGL API
//
inline int		clamp(int lo, int x, int hi)
{
	hi<<=1;
	int temp=x+lo+abs(x-lo);
	return (temp+hi-abs(temp-hi))>>2;
}
namespace		resources
{
	const char sf10_height=16, sf8_height=12, sf7_height=10, sf6_height=8;
	const char sf10_widths[]=
	{
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,4,4,6,8,8,11,9,4,
		4,4,6,8,4,4,4,4,8,8,
		8,8,8,8,8,8,8,8,4,4,
		8,8,8,8,14,8,10,9,10,9,
		8,10,10,4,7,9,8,12,10,10,
		9,10,10,9,8,10,8,14,9,10,
		9,4,4,4,5,8,5,8,8,7,
		8,8,4,8,8,4,4,7,4,12,
		8,8,8,8,5,8,4,8,8,10,
		8,8,8,5,4,5,5,0
	};
	const char sf8_widths[]=
	{
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,2,3,4,6,6,9,7,2,
		3,3,4,6,3,5,3,3,6,4,
		6,6,6,6,6,6,6,6,3,3,
		6,6,6,6,10,7,7,7,7,6,
		6,8,7,2,5,7,6,8,7,8,
		7,8,7,7,6,7,7,9,7,7,
		6,3,3,3,3,6,3,6,6,6,
		6,6,3,6,6,2,2,5,2,8,
		6,6,6,6,3,5,3,6,5,7,
		5,5,5,3,3,3,3,0,
	};
	const char sf7_widths[]=
	{
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,2,2,4,6,4,5,5,2,
		3,3,3,4,2,3,2,3,5,3,
		5,5,5,4,5,4,5,5,2,2,
		4,4,4,5,8,6,6,6,6,5,
		5,6,5,2,4,5,4,8,6,6,
		5,6,6,5,5,6,6,8,6,6,
		5,3,3,3,4,5,3,4,5,4,
		5,4,3,5,4,2,2,4,2,6,
		4,5,5,5,3,4,3,4,4,6,
		4,4,4,3,2,3,4,0,
	};
	char sf6_widths[]=
	{
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,
		0,0,1,2,3,5,4,5,5,2,
		3,3,3,4,2,3,2,3,4,3,
		4,4,4,4,4,4,4,4,2,2,
		3,3,3,3,7,5,5,5,5,4,
		4,5,5,2,4,5,4,6,5,5,
		5,5,5,5,4,5,4,6,4,4,
		4,3,3,3,4,4,3,4,4,4,
		4,4,3,4,4,2,2,4,2,6,
		4,4,4,4,3,3,3,4,4,5,
		4,4,4,3,2,3,4,0,
	};
	const unsigned char sf10[]=//system font size 10+
	{
		25,209,255,159,199,24,205,243,60,99,80,217,216,252,217,216,108,108,254,108,108,50,182,140,55,207,131,7,195,179,199,96,86,29,188,177,205,108,195,13,176,195,54,179,141,61,184,204,41,199,
		102,179,113,28,219,249,204,159,17,250,207,76,123,219,182,109,155,49,211,179,109,219,182,61,140,21,243,51,158,148,178,99,24,230,103,24,6,59,99,15,72,241,89,145,199,16,51,51,155,153,
		205,204,100,76,189,121,158,231,121,158,103,47,67,202,207,204,204,204,156,49,245,230,25,198,24,99,24,254,140,169,55,195,48,7,195,240,236,101,76,97,156,231,109,123,254,97,56,99,250,15,
		195,240,205,48,60,123,25,83,111,30,134,111,158,231,217,203,152,254,48,134,49,12,99,24,70,198,212,155,231,217,155,231,121,246,50,166,222,60,207,179,15,195,179,151,18,124,128,71,77,218,
		0,96,47,100,9,99,140,49,24,12,6,199,140,126,0,248,67,150,131,193,96,48,198,24,67,198,212,155,103,24,99,24,128,97,34,220,240,192,57,6,102,110,179,61,219,155,189,205,219,108,
		123,6,192,113,240,129,65,197,192,224,225,33,49,51,243,27,30,206,160,254,13,15,15,255,13,15,15,15,255,101,78,121,230,225,96,48,24,12,205,188,12,234,207,216,240,240,240,240,240,240,
		216,79,230,244,63,24,12,126,131,193,96,240,207,156,254,7,131,193,111,48,24,12,6,25,84,62,227,193,129,129,249,225,97,99,222,12,234,240,240,240,240,255,240,240,240,240,112,70,244,255,
		255,49,166,48,12,195,48,12,207,179,151,65,29,155,217,120,56,120,216,152,25,27,206,156,14,6,131,193,96,48,24,12,254,25,213,129,7,62,252,240,231,159,223,123,239,153,103,206,160,14,
		31,63,63,111,111,207,207,143,15,103,80,121,204,134,135,135,135,135,135,205,120,50,168,127,195,195,195,195,127,3,3,3,3,25,84,30,179,225,225,225,225,225,121,51,254,76,234,159,97,195,
		134,13,251,51,108,216,176,193,153,83,223,120,60,112,224,192,227,177,15,131,250,143,129,129,129,129,129,129,129,129,145,65,29,30,30,30,30,30,30,30,54,227,193,160,14,15,155,153,153,145,
		240,240,96,96,96,92,135,225,97,120,24,54,207,204,51,150,134,231,193,48,48,12,12,3,147,58,120,176,49,54,56,112,176,49,54,120,48,70,117,224,129,13,99,6,15,24,96,128,1,6,
		24,48,169,255,1,3,3,3,3,3,3,3,3,254,103,166,191,109,219,182,237,24,226,153,25,51,51,102,102,204,244,109,219,182,109,31,83,67,220,198,99,252,71,104,206,56,102,231,205,188,
		61,207,126,198,116,24,134,111,158,231,121,254,98,118,222,60,12,195,236,101,76,97,24,246,231,121,158,103,63,102,231,205,255,48,204,30,134,148,205,222,204,204,44,102,234,207,243,60,207,62,
		60,123,25,211,97,24,190,121,158,231,121,206,136,62,255,63,102,218,134,109,219,182,151,49,29,134,97,222,158,227,217,230,140,232,255,255,199,244,254,155,121,230,153,103,158,121,230,152,221,55,
		207,243,60,207,49,59,111,158,231,121,246,98,166,223,60,207,243,252,13,195,16,51,245,231,121,158,103,31,134,225,24,221,223,204,204,196,236,188,121,240,224,217,3,81,154,189,153,153,113,204,
		110,158,231,121,158,125,24,222,240,176,153,25,15,6,6,76,111,224,153,109,179,205,63,102,152,129,225,13,155,241,96,240,152,13,195,80,135,135,205,204,120,120,48,48,24,12,98,118,63,140,
		49,198,240,103,168,217,204,204,198,204,204,56,35,253,255,255,63,134,122,204,204,140,205,204,108,48,37,183,3,
	};
	const unsigned char sf8[]=//small fonts size 8
	{
		145,184,111,100,164,141,212,41,245,213,87,42,82,72,92,53,172,58,146,185,27,142,32,8,130,56,28,49,52,146,196,72,163,200,141,132,28,145,180,170,50,34,153,170,90,68,228,67,85,132,
		124,66,146,34,22,163,248,100,196,136,160,106,85,164,208,197,24,99,140,46,50,200,147,36,25,41,116,17,34,34,194,143,20,186,8,25,132,209,69,10,17,83,41,125,132,34,133,63,132,7,
		97,116,145,66,23,195,139,49,186,72,225,135,8,17,66,72,164,208,197,232,98,140,46,82,232,98,244,16,70,151,138,138,169,176,56,67,135,36,33,132,41,154,135,207,208,133,16,146,68,164,
		208,69,136,8,1,36,131,200,35,148,89,89,153,38,192,19,49,100,24,73,146,254,48,140,24,126,97,248,133,97,248,69,12,189,48,8,130,32,244,34,134,95,24,134,97,24,126,145,194,31,
		194,11,33,252,72,225,15,225,133,16,66,228,208,23,12,4,226,193,160,47,98,24,134,225,31,134,97,24,9,252,143,16,34,34,98,166,69,12,195,40,57,142,36,10,35,133,33,132,16,66,
		248,145,195,241,184,90,77,38,131,193,136,225,56,150,101,154,198,113,228,208,23,12,6,131,193,160,47,98,248,133,97,248,5,65,16,185,244,5,131,193,96,52,236,3,70,12,191,48,12,191,
		48,12,35,134,94,24,120,32,24,122,145,194,79,8,33,132,144,136,97,24,134,97,24,134,94,196,48,12,67,73,146,24,38,130,24,24,24,40,164,165,69,66,18,49,12,67,137,97,164,48,
		132,28,6,131,34,10,2,129,64,68,10,63,68,68,132,240,35,146,87,85,103,228,170,213,136,228,170,234,144,17,21,106,241,71,68,50,85,229,160,143,62,82,24,66,120,49,198,151,170,114,
		49,68,23,41,132,16,250,24,163,79,85,185,248,131,139,8,90,87,165,234,124,140,209,195,23,41,12,33,180,25,99,204,196,245,51,145,253,35,132,17,145,53,149,145,192,255,212,213,47,153,
		76,38,83,85,47,198,24,83,85,46,198,232,82,117,47,198,248,66,72,213,249,24,163,135,48,37,117,85,138,202,195,240,66,98,93,153,170,138,49,102,155,162,202,76,179,84,85,172,85,169,
		20,85,166,101,166,232,50,83,73,34,69,245,36,241,17,73,213,84,35,145,255,35,146,85,86,65,70,60,
	};
	const unsigned char sf7[]=//small fonts size 7
	{
		145,176,55,50,210,70,202,212,87,234,171,200,156,62,249,138,144,205,146,52,71,200,164,164,75,35,33,71,4,173,42,35,130,169,106,17,145,15,205,232,202,133,156,146,56,23,49,34,166,86,
		17,50,203,204,180,136,152,171,70,200,44,209,226,35,100,150,132,105,17,50,50,171,167,200,216,51,114,17,50,139,203,180,200,216,41,73,17,50,75,203,180,8,153,165,99,90,40,40,67,81,
		57,51,69,69,12,205,28,207,76,69,84,34,100,150,36,32,145,59,142,40,89,205,9,56,145,50,66,148,139,49,82,246,226,139,241,69,202,92,12,33,186,72,217,139,49,198,23,33,251,184,
		136,143,144,125,92,68,68,202,92,12,57,250,8,89,230,103,102,36,236,143,140,145,100,23,33,203,154,169,140,140,37,73,30,57,27,143,171,213,100,50,82,54,103,173,57,71,202,92,140,49,
		186,8,217,101,94,68,164,204,197,88,147,141,148,189,24,95,140,17,50,75,161,52,72,217,39,132,16,18,41,139,49,198,232,34,101,49,42,69,72,228,44,24,84,169,40,148,72,89,84,132,
		168,24,41,139,138,16,66,34,100,143,36,241,17,185,171,58,34,86,169,17,185,85,61,50,162,2,43,126,66,36,67,67,99,29,33,139,184,204,11,13,57,113,132,12,209,51,61,52,164,198,
		17,49,187,10,145,121,166,163,69,198,146,181,141,132,245,145,192,254,200,88,210,181,145,176,63,84,244,106,173,161,161,181,13,17,89,166,133,200,46,243,34,66,100,158,233,136,33,161,171,208,
		144,195,101,164,186,12,13,181,117,104,168,149,66,69,181,42,21,26,170,212,208,88,91,41,161,161,51,143,200,169,169,145,184,127,68,174,178,130,144,104,1,
	};
	const unsigned char sf6[]=//small fonts size 6
	{
		145,168,27,17,249,8,149,189,189,69,198,244,249,202,16,37,73,70,168,164,180,52,18,114,68,204,42,35,98,169,22,17,249,204,140,174,88,200,33,137,99,17,35,82,122,69,166,212,86,17,
		41,87,35,83,163,242,200,212,40,46,50,69,247,140,76,61,231,34,83,202,170,200,212,41,41,50,165,170,138,76,169,167,66,49,13,5,117,72,198,76,73,28,146,201,34,82,25,69,204,24,
		169,237,5,38,66,101,150,159,17,170,203,203,139,80,89,70,90,132,234,50,243,34,83,207,242,200,212,179,36,66,229,209,233,17,170,204,207,140,68,253,200,20,201,46,66,149,53,149,145,169,
		36,121,164,42,238,26,99,132,42,183,51,35,84,150,153,22,161,186,188,136,8,149,101,150,70,168,46,47,51,66,229,97,120,145,169,75,82,132,42,51,211,34,83,109,165,72,85,172,85,169,
		200,84,171,54,50,213,74,138,76,157,202,35,98,87,29,145,170,53,34,182,234,145,17,21,71,241,17,145,12,205,184,142,76,37,235,66,51,142,35,83,164,235,208,140,186,136,148,93,161,33,
		215,140,76,37,107,35,81,61,18,214,71,166,146,174,145,168,31,170,121,181,134,102,214,134,102,222,135,134,214,37,52,228,154,33,153,11,52,99,45,35,212,25,154,105,29,154,105,21,162,201,
		252,208,76,213,208,80,165,132,102,174,71,196,180,52,18,245,35,98,165,21,25,241,0,
	};
	void			uncompress_bitmap_v4(const unsigned char *data, int data_size, char fontH, const char *widths, int *&rgb2, int &w2, int &h2)
	{
		const int red=0xFF0000FF, green=0xFF00FF00, blue=0xFFFF0000;//swapped in WinAPI bitmaps
		w2=128, h2=256;
		const int size=w2*h2;
		rgb2=(int*)malloc(size<<2);
		for(int k=0;k<size;++k)
			rgb2[k]=blue;
		int nbits=data_size<<3;
		std::vector<bool> cdata(nbits);
		for(int k=0;k<nbits;k+=8)
		{
			for(int bit=0;bit<8;++bit)
				cdata[k+bit]=data[k>>3]>>bit&1;
		}
		{
			int c=' ';
			int xstart=(c&7)<<4, ystart=(c>>3)<<4;
			int bkwidth=widths[c];
			for(int ky=0;ky<fontH;++ky)
			{
				for(int kx=0;kx<bkwidth;++kx)
				{
					rgb2[(ystart+ky)<<7|(xstart+kx)]=green;
				}
			}
		}
		for(int c='!', pos=0;c<127;++c)//95 printable characters, space is all-bk
		{
			int xstart=(c&7)<<4, ystart=(c>>3)<<4;
			int bkwidth=widths[c];
			for(int ky=0;ky<fontH;++ky)
			{
				for(int kx=0;kx<bkwidth;++kx)
				{
					rgb2[(ystart+ky)<<7|(xstart+kx)]=green;
				}
			}
			int	xoffset=(int)cdata[pos+2]<<2|(int)cdata[pos+1]<<1|(int)cdata[pos],
				yoffset=(int)cdata[pos+6]<<3|(int)cdata[pos+5]<<2|(int)cdata[pos+4]<<1|(int)cdata[pos+3],
				width=(int)cdata[pos+10]<<3|(int)cdata[pos+9]<<2|(int)cdata[pos+8]<<1|(int)cdata[pos+7],
				height=(int)cdata[pos+14]<<3|(int)cdata[pos+13]<<2|(int)cdata[pos+12]<<1|(int)cdata[pos+11];
			xstart+=xoffset, ystart+=yoffset;
			pos+=15;
			for(int ky=0;ky<height;++ky)
			{
				for(int kx=0;kx<width;++kx)
				{
					if(cdata[pos])
						rgb2[(ystart+ky)<<7|(xstart+kx)]=red;
					++pos;
				}
			}
		}
	}
}
const float	_pi=acos(-1.f), _2pi=2*_pi, pi_2=_pi*0.5f, inv_2pi=1/_2pi, sqrt2=sqrt(2.f), torad=_pi/180, infinity=(float)_HUGE, inv255=1.f/255, inv256=1.f/256, inv128=1.f/128;
bool			API_not_loaded=true, one_shot=true;
//void			(__stdcall *glGenVertexArrays)(int n, unsigned *arrays)=nullptr;//OpenGL 3.0
//void			(__stdcall *glDeleteVertexArrays)(int n, unsigned *arrays)=nullptr;//OpenGL 3.0
void			(__stdcall *glBindVertexArray)(unsigned array)=nullptr;
void			(__stdcall *glGenBuffers)(int n, unsigned *buffers)=nullptr;
void			(__stdcall *glBindBuffer)(unsigned target, unsigned buffer)=nullptr;
void			(__stdcall *glBufferData)(unsigned target, int size, const void *data, unsigned usage)=nullptr;
void			(__stdcall *glBufferSubData)(unsigned target, int offset, int size, const void *data)=nullptr;
//void			(__stdcall *glGetBufferSubData)(unsigned target, int offset, int size, void *data)=nullptr;
void			(__stdcall *glEnableVertexAttribArray)(unsigned index)=nullptr;
void			(__stdcall *glVertexAttribPointer)(unsigned index, int size, unsigned type, unsigned char normalized, int stride, const void *pointer)=nullptr;
void			(__stdcall *glDisableVertexAttribArray)(unsigned index)=nullptr;
unsigned		(__stdcall *glCreateShader)(unsigned shaderType)=nullptr;
void			(__stdcall *glShaderSource)(unsigned shader, int count, const char **string, const int *length)=nullptr;
void			(__stdcall *glCompileShader)(unsigned shader)=nullptr;
void			(__stdcall *glGetShaderiv)(unsigned shader, unsigned pname, int *params)=nullptr;
void			(__stdcall *glGetShaderInfoLog)(unsigned shader, int maxLength, int *length, char *infoLog)=nullptr;
unsigned		(__stdcall *glCreateProgram)()=nullptr;
void			(__stdcall *glAttachShader)(unsigned program, unsigned shader)=nullptr;
void			(__stdcall *glLinkProgram)(unsigned program)=nullptr;
void			(__stdcall *glGetProgramiv)(unsigned program, unsigned pname, int *params)=nullptr;
void			(__stdcall *glGetProgramInfoLog)(unsigned program, int maxLength, int *length, char *infoLog)=nullptr;
void			(__stdcall *glDetachShader)(unsigned program, unsigned shader)=nullptr;
void			(__stdcall *glDeleteShader)(unsigned shader)=nullptr;
void			(__stdcall *glUseProgram)(unsigned program)=nullptr;
int				(__stdcall *glGetAttribLocation)(unsigned program, const char *name)=nullptr;
void			(__stdcall *glDeleteProgram)(unsigned program)=nullptr;
void			(__stdcall *glDeleteBuffers)(int n, const unsigned *buffers)=nullptr;
int				(__stdcall *glGetUniformLocation)(unsigned program, const char *name)=nullptr;
void			(__stdcall *glUniform1f)(int location, float v0)=nullptr;
void			(__stdcall *glUniformMatrix3fv)(int location, int count, unsigned char transpose, const float *value)=nullptr;
void			(__stdcall *glUniformMatrix4fv)(int location, int count, unsigned char transpose, const float *value)=nullptr;
void			(__stdcall *glGetBufferParameteriv)(unsigned target, unsigned value, int *data)=nullptr;
void			(__stdcall *glActiveTexture)(unsigned texture)=nullptr;
void			(__stdcall *glUniform1i)(int location, int v0)=nullptr;
void			(__stdcall *glUniform2f)(int location, float v0, float v1)=nullptr;
void			(__stdcall *glUniform3f)(int location, float v0, float v1, float v2)=nullptr;
void			(__stdcall *glUniform3fv)(int location, int count, const float *value)=nullptr;
void			(__stdcall *glUniform4f)(int location, float v0, float v1, float v2, float v3)=nullptr;
struct			ShaderVar
{
	int *pvar;
	const char *name;
	int lineNo;//__LINE__
};
unsigned		CompileShader(const char *src, unsigned type, int line)
{
	unsigned shaderID=glCreateShader(type);
	glShaderSource(shaderID, 1, &src, 0);
	glCompileShader(shaderID);
	int success=0;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		int infoLogLength;
		glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::vector<char> errorMessage(infoLogLength+1);
		glGetShaderInfoLog(shaderID, infoLogLength, 0, &errorMessage[0]);
		copy_to_clipboard(&errorMessage[0], infoLogLength);
		log_error(__FILE__, line, "Shader compilation failed. Output copied to clipboard.");
	//	GL_ERROR();
		return 0;
	}
	return shaderID;
}
unsigned		LoadShaders(const char *vertSrc, const char *fragSrc, ShaderVar *attributes, int n_attrib, ShaderVar *uniforms, int n_unif, int line)
{
	unsigned
		vertShaderID=CompileShader(vertSrc, GL_VERTEX_SHADER, line),
		fragShaderID=CompileShader(fragSrc, GL_FRAGMENT_SHADER, line);
//	prof_add("compile sh");
	if(!vertShaderID||!fragShaderID)
	{
		GL_ERROR();
		return 0;
	}
	unsigned ProgramID=glCreateProgram();
	glAttachShader(ProgramID, vertShaderID);
	glAttachShader(ProgramID, fragShaderID);
	glLinkProgram(ProgramID);
//	prof_add("link");
	int success=0;
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &success);
	if(!success)
	{
		int infoLogLength;
		glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::vector<char> errorMessage(infoLogLength+1);
		glGetProgramInfoLog(ProgramID, infoLogLength, 0, &errorMessage[0]);
		copy_to_clipboard(&errorMessage[0], infoLogLength);
		GL_ERROR();
		return 0;
	}
	glDetachShader(ProgramID, vertShaderID);
	glDetachShader(ProgramID, fragShaderID);
	glDeleteShader(vertShaderID);
	glDeleteShader(fragShaderID);
//	prof_add("delete");
	GL_CHECK();
	for(int ka=0;ka<n_attrib;++ka)
		if((*attributes[ka].pvar=glGetAttribLocation(ProgramID, attributes[ka].name))==-1)
			gl_error(__FILE__, attributes[ka].lineNo);
	for(int ku=0;ku<n_unif;++ku)
		if((*uniforms[ku].pvar=glGetUniformLocation(ProgramID, uniforms[ku].name))==-1)
			gl_error(__FILE__, uniforms[ku].lineNo);
	//if(broken)//
	//	return 0;//
	return ProgramID;
}
unsigned		make_gpu_buffer(unsigned target, const void *pointer, int size_bytes)
{
	unsigned buffer_id=0;
	glGenBuffers(1, &buffer_id);
	glBindBuffer(target, buffer_id);
	glBufferData(target, size_bytes, pointer, GL_STATIC_DRAW);
	return buffer_id;
}
void			LoadFontTexture(const unsigned char *data, int data_size, char fontH, const char *widths, unsigned &tx_id)
{
	int w2, h2;
	int *rgb2=nullptr;
	resources::uncompress_bitmap_v4(data, data_size, fontH, widths, rgb2, w2, h2);

	glGenTextures(1, &tx_id);
	glBindTexture(GL_TEXTURE_2D, tx_id);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w2, h2, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb2);
	GL_CHECK();

	free(rgb2);
}
void			calculate_text_txcoords(const char *font_widths, int gl_fontH, float *txcoord)
{
	for(int c=0, idx=0;c<128;++c, idx+=4)
	{
		int width=font_widths[c];
		int xpos=(c&0x7)<<4, ypos=(c>>3&0xF)<<4;
		txcoord[idx  ]=xpos*inv128, txcoord[idx+1]=(xpos+width)*inv128;
		txcoord[idx+2]=ypos*inv256, txcoord[idx+3]=(ypos+gl_fontH)*inv256;
	}
}
unsigned		current_program=0;
inline void		gl_setProgram(unsigned program)
{
	if(current_program!=program)
	{
	//	unsigned current=current_program;
	//	long long broken0=broken;
		glUseProgram(current_program=program);
	//	GL_CHECK();
	}
}
void			send_color(unsigned location, int color)
{
	auto p=(unsigned char*)&color;
	static const __m128 m_255=_mm_set1_ps(inv255);
	__m128 c=_mm_castsi128_ps(_mm_set_epi32(p[3], p[2], p[1], p[0]));
	c=_mm_cvtepi32_ps(_mm_castps_si128(c));
	c=_mm_mul_ps(c, m_255);
	glUniform4f(location, c.m128_f32[0], c.m128_f32[1], c.m128_f32[2], c.m128_f32[3]);
	//glUniform4f(location, p[0]*inv255, p[1]*inv255, p[2]*inv255, p[3]*inv255);
}
void			send_color_rgb(unsigned location, int color)
{
	auto p=(unsigned char*)&color;
	static const __m128 m_255=_mm_set1_ps(inv255);
	__m128 c=_mm_castsi128_ps(_mm_set_epi32(p[3], p[2], p[1], p[0]));
	c=_mm_cvtepi32_ps(_mm_castps_si128(c));
	c=_mm_mul_ps(c, m_255);
	glUniform3f(location, c.m128_f32[0], c.m128_f32[1], c.m128_f32[2]);
}
void			select_texture(unsigned tx_id, int u_location)
{
	glActiveTexture(GL_TEXTURE0);		GL_CHECK();
	glBindTexture(GL_TEXTURE_2D, tx_id);GL_CHECK();//select texture
	glUniform1i(u_location, 0);			GL_CHECK();
}
float			g_fbuf[16]={0};
int				pen_color=0xFF000000, brush_color=0xFFFFFFFF;
namespace		GL2_2D
{
	unsigned	program=0;
	int			u_color=-1, a_vertices=-1;
	unsigned	vertex_buffer=0;
	ivec4		region, current_region;//x1, y1, dx, dy
	bool		continuous=true;
	//void		set_pen_color(int color)
	//{
	//	gl_setProgram(program);
	//	send_color(u_color, color);
	//}
	struct		DrawInfo
	{
		int count, color;
		DrawInfo(int count, int color):count(count), color(color){}
	};
	std::vector<DrawInfo> drawInfo;
	std::vector<vec2> vertices;
	void		set_region(int x1, int x2, int y1, int y2){region.set(x1, y1, x2-x1, y2-y1);}
	void		use_region()
	{
		glViewport(region.x1, region.y1, region.dx, region.dy);
		current_region=region;
	}
	void		drop_region()
	{
		glViewport(0, 0, w, h);
		current_region.set(0, 0, w, h);
	}
	void		toNDC(float xs, float ys, float &xn, float &yn)
	{
		xn=(xs+0.5f-GL2_2D::current_region.x1)*2.f/GL2_2D::current_region.dx-1;
		yn=1-(ys+0.5f-GL2_2D::current_region.y1)*2.f/GL2_2D::current_region.dy;
	}
	//void		toNDC(int xs, int ys, float &xn, float &yn)
	//{
	//	xn=(xs-current_region.x1)*2.f/current_region.dx-1;
	//	yn=1-(ys-current_region.y1)*2.f/current_region.dy;
	//}
	void		curve_begin(){GL2_2D::vertices.clear(), GL2_2D::drawInfo.clear();}
	void		curve_point(float x, float y)
	{
		vec2 ndc;
		toNDC(x, y, ndc.x, ndc.y);
		GL2_2D::vertices.push_back(ndc);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//GL2_2D::vertices.push_back(vec2(x*_2_w-1, 1-y*_2_h));
		bool increment=false;
		if(GL2_2D::drawInfo.size())
		{
			auto &di=*GL2_2D::drawInfo.rbegin();
			increment=GL2_2D::continuous&(di.color==pen_color);
		}
		if(increment)
			++GL2_2D::drawInfo.rbegin()->count;
		else
			GL2_2D::drawInfo.push_back(GL2_2D::DrawInfo(1, pen_color));
		GL2_2D::continuous=true;
	}
	void		draw_curve()
	{
		using namespace GL2_2D;
		gl_setProgram(program);									GL_CHECK();
		
		glEnableVertexAttribArray(a_vertices);					GL_CHECK();
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);													GL_CHECK();
		glBufferData(GL_ARRAY_BUFFER, (int)vertices.size()*sizeof(vec2), &GL2_2D::vertices[0], GL_STATIC_DRAW);	GL_CHECK();
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);											GL_CHECK();
		
	//	glBindBuffer(GL_ARRAY_BUFFER, 0);												GL_CHECK();
	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, &vertices[0]);		GL_CHECK();//doesn't draw

		for(int k=0, kEnd=GL2_2D::drawInfo.size(), idx=0;k<kEnd;++k)
		{
			auto &di=drawInfo[k];
			send_color(GL2_2D::u_color, di.color);								GL_CHECK();
			glDrawArrays(di.count==1?GL_POINTS:GL_LINE_STRIP, idx, di.count);	GL_CHECK();
			idx+=di.count;
			prof_add("curve", 1);//
		}
		glDisableVertexAttribArray(a_vertices);					GL_CHECK();
	}
	void		draw_line(float x1, float y1, float x2, float y2)
	{
		toNDC(x1, y1, g_fbuf[0], g_fbuf[1]);
		toNDC(x2, y2, g_fbuf[2], g_fbuf[3]);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//g_fbuf[0]=x1*_2_w-1, g_fbuf[1]=1-y1*_2_h;
		//g_fbuf[2]=x2*_2_w-1, g_fbuf[3]=1-y2*_2_h;
		gl_setProgram(program);					GL_CHECK();
		send_color(u_color, pen_color);			GL_CHECK();
		glEnableVertexAttribArray(a_vertices);	GL_CHECK();
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);			GL_CHECK();
		glBufferData(GL_ARRAY_BUFFER, 4<<2, g_fbuf, GL_STATIC_DRAW);	GL_CHECK();//use glBufferSubData
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	GL_CHECK();//use buffer

		glDrawArrays(GL_LINES, 0, 2);			GL_CHECK();
	//	prof_add("point", 1);//
	}
	void		set_pixel(float x, float y)
	{
		toNDC(x, y, g_fbuf[0], g_fbuf[1]);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//g_fbuf[0]=x*_2_w-1, g_fbuf[1]=1-y*_2_h;
		//long long broken0=broken;
		gl_setProgram(GL2_2D::program);			GL_CHECK();
		send_color(u_color, pen_color);			GL_CHECK();
		glEnableVertexAttribArray(a_vertices);	GL_CHECK();
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);			GL_CHECK();
		glBufferData(GL_ARRAY_BUFFER, 2<<2, g_fbuf, GL_STATIC_DRAW);	GL_CHECK();
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	GL_CHECK();//use buffer

		glDrawArrays(GL_POINTS, 0, 1);			GL_CHECK();
	}
	void		draw_rectangle_hollow(float x1, float x2, float y1, float y2)
	{
		float X1, X2, Y1, Y2;
		toNDC(x1, y1, X1, Y1);
		toNDC(x2, y2, X2, Y2);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//float X1=x1*_2_w-1, X2=x2*_2_w-1, Y1=1-y1*_2_h, Y2=1-y2*_2_h;
		g_fbuf[0]=X1, g_fbuf[1]=Y1;
		g_fbuf[2]=X2, g_fbuf[3]=Y1;
		g_fbuf[4]=X2, g_fbuf[5]=Y2;
		g_fbuf[6]=X1, g_fbuf[7]=Y2;
		g_fbuf[8]=X1, g_fbuf[9]=Y1;
	//	g_fbuf[10]=X1, g_fbuf[11]=Y1;
		gl_setProgram(program);						GL_CHECK();
		send_color(u_color, 0xFF000000|pen_color);	GL_CHECK();

		glBindBuffer(GL_ARRAY_BUFFER, 0);									GL_CHECK();
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);GL_CHECK();

		glEnableVertexAttribArray(a_vertices);	GL_CHECK();
		glDrawArrays(GL_LINE_STRIP, 0, 5);		GL_CHECK();
	//	glDrawArrays(GL_POINTS, 5, 1);			GL_CHECK();
	}
	void		draw_rectangle(float x1, float x2, float y1, float y2, bool opaque)
	{
		float X1, X2, Y1, Y2;
		toNDC(x1, y1, X1, Y1);
		toNDC(x2, y2, X2, Y2);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//float X1=x1*_2_w-1, X2=x2*_2_w-1, Y1=1-y1*_2_h, Y2=1-y2*_2_h;
		g_fbuf[0]=X1, g_fbuf[1]=Y1;
		g_fbuf[2]=X2, g_fbuf[3]=Y1;
		g_fbuf[4]=X2, g_fbuf[5]=Y2;
		g_fbuf[6]=X1, g_fbuf[7]=Y2;
		g_fbuf[8]=X1, g_fbuf[9]=Y1;
		gl_setProgram(program);							GL_CHECK();
		if(opaque)
			send_color(u_color, 0xFF000000|brush_color);
		else
			send_color(u_color, brush_color);
		GL_CHECK();
		glEnableVertexAttribArray(a_vertices);			GL_CHECK();
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);			GL_CHECK();
		glBufferData(GL_ARRAY_BUFFER, 10<<2, g_fbuf, GL_STATIC_DRAW);	GL_CHECK();
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	GL_CHECK();//use buffer

		glDrawArrays(GL_TRIANGLE_FAN, 0, 5);	GL_CHECK();
		if(brush_color!=pen_color)
			draw_rectangle_hollow(x1, x2, y1, y2);
	}
}
namespace		GL2_L3D//light 3D
{
	struct		TriangleInfo
	{
		int count, color;
		TriangleInfo(int count, int color):count(count), color(color){}
	};
	unsigned	program=0;
	int			a_vertices=-1, a_normals=-1, u_vpmatrix=-1, u_modelmatrix=-1, u_normalmatrix=-1, u_objectcolor=-1, u_lightcolor=-1, u_lightpos=-1, u_viewpos=-1;
	unsigned	VBO=0, EBO=0;
	std::vector<vec3> vertices;
	std::vector<int> indices;
	std::vector<TriangleInfo> drawInfo;
	void		begin()
	{
		vertices.clear();
		indices.clear();
		drawInfo.clear();
	}
	void		push_surface(vec3 const *vn, int vcount_x2, int *idx, int icount, int color)
	{
		color|=0x7F000000;
		int increment=GL2_L3D::vertices.size()>>1, istart=GL2_L3D::indices.size();
		GL2_L3D::vertices.insert(vertices.end(), vn, vn+vcount_x2);
		GL2_L3D::indices.insert(indices.end(), idx, idx+icount);
		for(int k=istart, kEnd=indices.size();k<kEnd;++k)//use glDrawElementsBaseVertex
		{
			auto &ik=indices[k];
			ik=(ik>>1)+increment;
		}
		//	indices[k]+=increment;
		GL2_L3D::drawInfo.push_back(TriangleInfo(icount, color));
	}
	void		end()
	{
	//	vbo_to_clipboard(&GL2_L3D::vertices[0], GL2_L3D::vertices.size(), &GL2_L3D::indices[0], GL2_L3D::indices.size(), 0);
		glBindBuffer(GL_ARRAY_BUFFER, GL2_L3D::VBO);															GL_CHECK();
		glBufferData(GL_ARRAY_BUFFER, (int)vertices.size()*sizeof(vec3), &GL2_L3D::vertices[0], GL_STATIC_DRAW);GL_CHECK();
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GL2_L3D::EBO);													GL_CHECK();
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (int)indices.size()<<2, &GL2_L3D::indices[0], GL_STATIC_DRAW);	GL_CHECK();
	}
	void		draw(Camera const &cam, vec3 const &lightpos)
	{
		gl_setProgram(GL2_L3D::program);
		mat4
			mView=matrixFPSViewRH(cam.p, (float)cam.a.y, (float)cam.a.x-_pi),
			mProj=perspective((float)cam.tanfov, float(w)/h, 0.1f, 1000.f),
			vp=mProj*mView;
		mat4 model(1);
		mat3 m_normal=normalMatrix(model);
		mat4 mvp=vp*model;
		glUniformMatrix4fv(GL2_L3D::u_vpmatrix, 1, GL_FALSE, mvp.data());			GL_CHECK();
		glUniformMatrix4fv(GL2_L3D::u_modelmatrix, 1, GL_FALSE, model.data());		GL_CHECK();
		glUniformMatrix3fv(GL2_L3D::u_normalmatrix, 1, GL_FALSE, m_normal.data());	GL_CHECK();
		send_color_rgb(GL2_L3D::u_lightcolor, 0xFFFFFF);							GL_CHECK();
		glUniform3fv(GL2_L3D::u_lightpos, 1, &lightpos.x);							GL_CHECK();
		vec3 cam_p=cam.p;
		glUniform3fv(GL2_L3D::u_viewpos, 1, &cam_p.x);								GL_CHECK();
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_L3D::VBO);											GL_CHECK();
		glVertexAttribPointer(GL2_L3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 6<<2, 0);				GL_CHECK();
		glEnableVertexAttribArray(GL2_L3D::a_vertices);											GL_CHECK();
		glVertexAttribPointer(GL2_L3D::a_normals, 3, GL_FLOAT, GL_FALSE, 6<<2, (void*)(3<<2));	GL_CHECK();
		glEnableVertexAttribArray(GL2_L3D::a_normals);											GL_CHECK();
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GL2_L3D::EBO);									GL_CHECK();
		glDepthMask(GL_FALSE);														GL_CHECK();//all L3D objects are transparent for now

		for(int k=0, kEnd=drawInfo.size(), start=0;k<kEnd;++k)
		{
			auto &dk=GL2_L3D::drawInfo[k];
			send_color(GL2_L3D::u_objectcolor, dk.color);										GL_CHECK();
			glDrawElements(GL_TRIANGLES, dk.count, GL_UNSIGNED_INT, (void*)(start<<2));			GL_CHECK();//use glDrawElementsBaseVertex
			start+=dk.count;
		}

		//int sum=0;
		//for(int k=0;k<drawInfo.size();++k)
		//	sum+=GL2_L3D::drawInfo[k].count;
		//send_color(GL2_L3D::u_objectcolor, drawInfo[0].color);						GL_CHECK();//
		//glDrawElements(GL_TRIANGLES, sum, GL_UNSIGNED_INT, 0);						GL_CHECK();//

		glDepthMask(GL_TRUE);														GL_CHECK();
		glDisableVertexAttribArray(GL2_L3D::a_vertices);							GL_CHECK();
		glDisableVertexAttribArray(GL2_L3D::a_normals);								GL_CHECK();
	}
	void		draw_buffer(Camera const &cam, GPUBuffer const &buffer, vec3 const &modelpos, vec3 const &lightpos, int objcolor)
	{
		gl_setProgram(GL2_L3D::program);
		mat4
			mView=matrixFPSViewRH(cam.p, (float)cam.a.y, (float)cam.a.x-_pi),
		//	mView=matrixFPSViewRH(vec3((float)cam.p.x, (float)cam.p.y, (float)cam.p.z), cam.a.y, cam.a.x-_pi),
			mProj=perspective((float)cam.tanfov, float(w)/h, 0.1f, 1000.f),
			vp=mProj*mView;
		mat4 model=translate(mat4(1), modelpos);
		mat3 m_normal=normalMatrix(model);
		mat4 mvp=vp*model;
		glUniformMatrix4fv(GL2_L3D::u_vpmatrix, 1, GL_FALSE, mvp.data());			GL_CHECK();
		glUniformMatrix4fv(GL2_L3D::u_modelmatrix, 1, GL_FALSE, model.data());		GL_CHECK();
		glUniformMatrix3fv(GL2_L3D::u_normalmatrix, 1, GL_FALSE, m_normal.data());	GL_CHECK();
		send_color(GL2_L3D::u_objectcolor, 0x7F000000|objcolor);					GL_CHECK();
	//	send_color(GL2_L3D::u_objectcolor, 0xFF0000FF);								GL_CHECK();
		send_color_rgb(GL2_L3D::u_lightcolor, 0xFFFFFF);							GL_CHECK();
		glUniform3fv(GL2_L3D::u_lightpos, 1, &lightpos.x);							GL_CHECK();
		vec3 cam_p=cam.p;
		glUniform3fv(GL2_L3D::u_viewpos, 1, &cam_p.x);								GL_CHECK();
	//	glPointSize(10);															GL_CHECK();
		
		glBindBuffer(GL_ARRAY_BUFFER, buffer.VBO);												GL_CHECK();
		glVertexAttribPointer(GL2_L3D::a_vertices, 3, GL_FLOAT, GL_FALSE, buffer.vertices_stride, (void*)buffer.vertices_start);GL_CHECK();
		glEnableVertexAttribArray(GL2_L3D::a_vertices);											GL_CHECK();
		glVertexAttribPointer(GL2_L3D::a_normals, 3, GL_FLOAT, GL_FALSE, buffer.normals_stride, (void*)buffer.normals_start);	GL_CHECK();
		glEnableVertexAttribArray(GL2_L3D::a_normals);											GL_CHECK();
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.EBO);										GL_CHECK();

		glDrawElements(GL_TRIANGLES, buffer.n_vertices, GL_UNSIGNED_INT, 0);		GL_CHECK();
	//	glDrawArrays(GL_POINTS, 0, s_vcount);										GL_CHECK();
	}
}
namespace		GL2_3D
{
	struct		VertexInfo
	{
		int count;
		unsigned type;//GL_POINTS/GL_LINES/GL_LINE_STRIP/GL_TRIANGLES/GL_TRIANGLE_FAN/GL_QUADS/...
		bool textured;
		int color, *tx, txw, txh;
		VertexInfo(int count, unsigned type, int color):count(count), type(type), textured(false), color(color), tx(nullptr), txw(0), txh(0){}
		VertexInfo(int count, unsigned type, int *tx, int txw, int txh):count(count), type(type), textured(true), color(0), tx(tx), txw(txw), txh(txh){}
	};
	unsigned	program=0;
	int			a_vertices=-1, a_texcoord=-1, u_matrix=-1, u_pointSize=-1, u_texture=-1, u_isTextured=-1, u_color=-1;
	unsigned	vertex_buffer=0, texcoord_buffer=0;//TODO: interleave vertices & texcoords
	int			transparent_start_idx=0;
	std::vector<vec3> vertices;
	std::vector<vec2> texcoord;
	std::vector<VertexInfo> drawInfo;
	void		begin(){GL2_3D::vertices.clear(), GL2_3D::texcoord.clear(), GL2_3D::drawInfo.clear();}
	void		begin_transparent(){transparent_start_idx=GL2_3D::drawInfo.size();}
	void		push_square(float x1, float x2, float y1, float y2, float z, int *tx, int txw, int txh)
	{
		vertices.push_back(vec3(x1, y1, z)), texcoord.push_back(vec2(0, 0));
		vertices.push_back(vec3(x2, y1, z)), texcoord.push_back(vec2(1, 0));
		vertices.push_back(vec3(x2, y2, z)), texcoord.push_back(vec2(1, 1));
		vertices.push_back(vec3(x1, y2, z)), texcoord.push_back(vec2(0, 1));
		drawInfo.push_back(VertexInfo(4, GL_QUADS, tx, txw, txh));
	}
	//void		push_triangle(vec3 const &a, vec3 const &b, vec3 const &c, int color, int *tx=0, int txw=0, int txh=0)
	//{
	//	vertices.push_back(a), vertices.push_back(b), vertices.push_back(c);
	//	bool increment=false;
	//	if(drawInfo.size())
	//	{
	//		auto &di=*drawInfo.rbegin();
	//		increment=di.type==GL_TRIANGLES&&di.color==color&&di.tx==tx&&di.txw==txw&&di.txh==txh;
	//	}
	//	if(increment)
	//		drawInfo.rbegin()->count+=3;
	//	else if(tx)
	//		drawInfo.push_back(VertexInfo(3, GL_TRIANGLES, tx, txw, txh));//missing texcoord
	//	else
	//		drawInfo.push_back(VertexInfo(3, GL_TRIANGLES, color));
	//}
	void		push_triangle(vec3 const &a, vec3 const &b, vec3 const &c, int color)
	{
		color|=0x7F000000;
		vertices.push_back(a), texcoord.push_back(vec2(0, 0));
		vertices.push_back(b), texcoord.push_back(vec2(0, 0));
		vertices.push_back(c), texcoord.push_back(vec2(0, 0));
		bool increment=false;
		if(drawInfo.size())
		{
			auto &di=*drawInfo.rbegin();
			increment=di.type==GL_TRIANGLES&&di.color==color;
		}
		if(increment)
			drawInfo.rbegin()->count+=3;
		else
			drawInfo.push_back(VertexInfo(3, GL_TRIANGLES, color));
	}
	void		push_line_segment(vec3 const &p1, vec3 const &p2, int color)
	{
		color|=0xFF000000;
		vertices.push_back(p1), texcoord.push_back(vec2(0, 0));
		vertices.push_back(p2), texcoord.push_back(vec2(0, 0));
		bool increment=false;
		if(drawInfo.size())
		{
			auto &di=*drawInfo.rbegin();
			increment=di.type==GL_LINES&&di.color==color;
		}
		if(increment)
			drawInfo.rbegin()->count+=2;
		else
			drawInfo.push_back(VertexInfo(2, GL_LINES, color));
	}
	void		push_point(vec3 const &p, int color)
	{
		color|=0xFF000000;
		vertices.push_back(p), texcoord.push_back(vec2(0, 0));
		bool increment=false;
		if(drawInfo.size())
		{
			auto &di=*drawInfo.rbegin();
			increment=di.type==GL_POINTS&&di.color==color;
		}
		if(increment)
			drawInfo.rbegin()->count+=1;
		else
			drawInfo.push_back(VertexInfo(1, GL_POINTS, color));
	}
	void		end()
	{
		if(vertices.size()&&texcoord.size())
		{
			glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::vertex_buffer);													GL_CHECK();
			glBufferData(GL_ARRAY_BUFFER, (int)vertices.size()*sizeof(vec3), &GL2_3D::vertices[0], GL_STATIC_DRAW);	GL_CHECK();
			glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::texcoord_buffer);													GL_CHECK();
			glBufferData(GL_ARRAY_BUFFER, (int)texcoord.size()*sizeof(vec2), &GL2_3D::texcoord[0], GL_STATIC_DRAW);	GL_CHECK();
		}
	}
	void		draw(Camera const &cam)
	{
		gl_setProgram(GL2_3D::program);		GL_CHECK();

		mat4
			mView=matrixFPSViewRH(vec3((float)cam.p.x, (float)cam.p.y, (float)cam.p.z), (float)cam.a.y, (float)cam.a.x-_pi),
			mProj=perspective((float)cam.tanfov, float(w)/h, 0.1f, 1000.f),
			mVP=mProj*mView;
		//mat4 mVP(1);
		glUniformMatrix4fv(GL2_3D::u_matrix, 1, GL_FALSE, mVP.data());					GL_CHECK();

		static unsigned tx_id=0;
		if(!tx_id)
			{glGenTextures(1, &tx_id);			GL_CHECK();}

		glEnableVertexAttribArray(GL2_3D::a_vertices);							GL_CHECK();
		glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);							GL_CHECK();
		glVertexAttribPointer(GL2_3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();
		glEnableVertexAttribArray(GL2_3D::a_texcoord);							GL_CHECK();
		glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::texcoord_buffer);					GL_CHECK();
		glVertexAttribPointer(GL2_3D::a_texcoord, 2, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();
		glPointSize(1);
		for(int k=0, p_idx=0, kEnd=drawInfo.size();k<kEnd;++k)
		{
			auto &di=GL2_3D::drawInfo[k];
			if(k==transparent_start_idx)
				glDepthMask(GL_FALSE);
			glUniform1i(GL2_3D::u_isTextured, di.textured);
			send_color(GL2_3D::u_color, di.textured?0xFFFF00FF:di.color);		GL_CHECK();//dummy color if textured
			if(di.textured)
			{
				glBindTexture(GL_TEXTURE_2D, tx_id);															GL_CHECK();
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);								GL_CHECK();
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);								GL_CHECK();
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, di.txw, di.txh, 0, GL_RGBA, GL_UNSIGNED_BYTE, di.tx);	GL_CHECK();//TODO: send textures to GPU with the vertices
				select_texture(tx_id, GL2_3D::u_texture);														GL_CHECK();

				//glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::texcoord_buffer);					GL_CHECK();
				//glVertexAttribPointer(GL2_3D::a_texcoord, 2, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();
			}
			else
				{select_texture(tx_id, GL2_3D::u_texture);		GL_CHECK();}//dummy texture if not textured
			glDrawArrays(di.type, p_idx, di.count);				GL_CHECK();
			p_idx+=di.count;
		}
		glDisableVertexAttribArray(GL2_3D::a_vertices);	GL_CHECK();
		glDisableVertexAttribArray(GL2_3D::a_texcoord);	GL_CHECK();
	//	glBindBuffer(GL_ARRAY_BUFFER, 0);				GL_CHECK();
		glDepthMask(GL_TRUE);
	}
}
namespace		GL2_Text//font & textures
{
	unsigned		program=0;
	int				a_coord2d=-1,
				//	a_texcoord=-1,
					u_mytexture=-1, u_txtColor=-1, u_bkColor=-1, u_isTexture=-1;
	unsigned		buffer=0;
	unsigned		tx_id_sf10=0, tx_id_sf8=0, tx_id_sf7=0, tx_id_sf6=0,
					font_tx_id=0;
}
float			sf10_txcoords[128<<2]={0}, sf8_txcoords[128<<2]={0}, sf7_txcoords[128<<2]={0}, sf6_txcoords[128<<2]={0},//txx1 txx2 txy1 txy2
				*text_txcoords=nullptr;
const char		*font_widths=nullptr;
int				gl_tabWidth=0;
int				gl_font_size=0, pixel_x=1, pixel_y=1, gl_fontH=resources::sf10_height,
				txtColor=0xFF000000, bkColor=0xFFFFFFFF;//WinAPI: 0xAARRGGBB, OpenGL: 0xAABBGGRR
void			gl_setTextColor(int color)
{
	gl_setProgram(GL2_Text::program);//
	send_color(GL2_Text::u_txtColor, color);
}
void			gl_setBkColor(int color)
{
	gl_setProgram(GL2_Text::program);//
	send_color(GL2_Text::u_bkColor, color);
//	GL_CHECK();
}
int				gl_setTextSize(int size)//{0, ..., 9}
{
	using namespace GL2_Text;
	using namespace resources;
	switch(size=clamp(0, size, 9))
	{
	case 0:font_tx_id=tx_id_sf6,	text_txcoords=sf6_txcoords,		font_widths=sf6_widths,		pixel_x=pixel_y=1,		gl_fontH= 8, gl_tabWidth=32;	break;
	case 1:font_tx_id=tx_id_sf7,	text_txcoords=sf7_txcoords,		font_widths=sf7_widths,		pixel_x=pixel_y=1,		gl_fontH=10, gl_tabWidth=40;	break;
	case 2:font_tx_id=tx_id_sf8,	text_txcoords=sf8_txcoords,		font_widths=sf8_widths,		pixel_x=pixel_y=1,		gl_fontH=12, gl_tabWidth=48;	break;
	case 3:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=1,		gl_fontH=16, gl_tabWidth=64;	break;
	case 4:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=2,		gl_fontH=16, gl_tabWidth=136;	break;
	case 5:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=3,		gl_fontH=16, gl_tabWidth=200;	break;
	case 6:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=4,		gl_fontH=16, gl_tabWidth=264;	break;
	case 7:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=5, pixel_y=6,	gl_fontH=16, gl_tabWidth=328;	break;
	case 8:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=5, pixel_y=7,	gl_fontH=16, gl_tabWidth=328;	break;
	case 9:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=5, pixel_y=8,	gl_fontH=16, gl_tabWidth=328;	break;
	}
	return size;
}
std::vector<float> vrtx;
int				print_array(int x, int y, const char *msg, int msg_length, int tab_origin)//glBufferData stutters on debug		TODO: send all text in one buffer per frame
{
	gl_setProgram(GL2_Text::program);				GL_CHECK();
	glUniform2f(GL2_Text::u_isTexture, 0, 1);		GL_CHECK();
//	glUniform1i(GL2_Text::u_isTexture, false);		GL_CHECK();
	int msg_width=0;
	float _2_w=2.f/w, _2_h=2.f/h;
	select_texture(GL2_Text::font_tx_id, GL2_Text::u_mytexture);
	vrtx.resize(msg_length<<4);//vx, vy, txx, txy		6.5fps
	float rect[4], *txc;
	int width, idx;
	int fontH_px=gl_fontH*pixel_y;
//	prof_add_loop(2);
	for(int k=0;k<msg_length;++k)
	{
		char c=msg[k];
		width=0;
		if(c=='\t')
			width=gl_tabWidth-(msg_width-tab_origin)%gl_tabWidth, c=' ';
		else if(c>=32&&c<0xFF&&(width=font_widths[c]))
			width*=pixel_x;
		rect[0]=(x+msg_width		)*_2_w-1, rect[1]=1- y			*_2_h;
		rect[2]=(x+msg_width+width	)*_2_w-1, rect[3]=1-(y+fontH_px)*_2_h;//y2<y1
		idx=k<<4;
		txc=text_txcoords+(c<<2);
		vrtx[idx   ]=rect[0], vrtx[idx+ 1]=rect[1],		vrtx[idx+ 2]=txc[0], vrtx[idx+ 3]=txc[2];//top left
		vrtx[idx+ 4]=rect[0], vrtx[idx+ 5]=rect[3],		vrtx[idx+ 6]=txc[0], vrtx[idx+ 7]=txc[3];//bottom left
		vrtx[idx+ 8]=rect[2], vrtx[idx+ 9]=rect[3],		vrtx[idx+10]=txc[1], vrtx[idx+11]=txc[3];//bottom right
		vrtx[idx+12]=rect[2], vrtx[idx+13]=rect[1],		vrtx[idx+14]=txc[1], vrtx[idx+15]=txc[2];//top right

		msg_width+=width;
	}
//	prof_add_loop(3);
	glBindBuffer(GL_ARRAY_BUFFER, GL2_Text::buffer);							GL_CHECK();
	glBufferData(GL_ARRAY_BUFFER, msg_length<<6, vrtx.data(), GL_STATIC_DRAW);	GL_CHECK();//set vertices & texcoords
//	prof_add_loop(4);
	glVertexAttribPointer(GL2_Text::a_coord2d, 4, GL_FLOAT, GL_FALSE, 0, 0);	GL_CHECK();
//	prof_add_loop(5);

	glEnableVertexAttribArray(GL2_Text::a_coord2d);		GL_CHECK();
	glDrawArrays(GL_QUADS, 0, msg_length<<2);			GL_CHECK();//draw the quads: 4 vertices per character quad
	glDisableVertexAttribArray(GL2_Text::a_coord2d);	GL_CHECK();
	prof_add_loop(7);
	return msg_width;
}
int				gl_print(int x, int y, int tab_origin, const char *format, ...)
{
	int msg_length=vsnprintf_s(g_buf, g_buf_size, format, (char*)(&format+1));
	return print_array(x, y, g_buf, msg_length, tab_origin);//38fps 41fps	48fps		undebugged 6.4fps
//	return print_array(x, y, g_buf, msg_length);//7.6fps
}
namespace		GL2_Texture
{//font & textures
	unsigned	program=0;
	int			a_coords=-1,
				u_mytexture=-1, u_alpha=-1;
//	unsigned	buffer=0;
}
void 			generate_glcl_texture(unsigned &tx_id, int Xplaces, int Yplaces)
{
	if(!tx_id)
		{glGenTextures(1, &tx_id);											GL_CHECK();}//generate texture id once
	glBindTexture(GL_TEXTURE_2D, tx_id);									GL_CHECK();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);			GL_CHECK();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);			GL_CHECK();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);		GL_CHECK();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);		GL_CHECK();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, Xplaces, Yplaces, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);	GL_CHECK();
}
void 			display_gl_texture(unsigned &tx_id)
{
	float _2_w=2.f/w, _2_h=2.f/h;
	float x1=0, x2=(float)w, y1=0, y2=(float)h;
	float rect[]=
	{
		x1*_2_w-1, 1-y1*_2_h,
		x2*_2_w-1, 1-y2*_2_h//y2<y1
	};
	float vrtx[]=
	{
		rect[0], rect[1],		0, 0,//top left
		rect[0], rect[3],		0, 1,//bottom left
		rect[2], rect[3],		1, 1,//bottom right

		rect[2], rect[3],		1, 1,//bottom right
		rect[2], rect[1],		1, 0,//top right
		rect[0], rect[1],		0, 0,//top left
	};
	gl_setProgram(GL2_Texture::program);
//	glUniform2f(GL2_Text::u_isTexture, 1, 0);			GL_CHECK();
	select_texture(tx_id, GL2_Texture::u_mytexture);
	glUniform1f(GL2_Texture::u_alpha, 1);		GL_CHECK();

	//if(!GL2_Texture::buffer)
	//	glGenBuffers(1, &GL2_Texture::buffer);
	//glBindBuffer(GL_ARRAY_BUFFER, GL2_Texture::buffer);	GL_CHECK();
	//glBufferData(GL_ARRAY_BUFFER, 24<<2, vrtx, GL_STATIC_DRAW);	GL_CHECK();//send vertices & texcoords
	//glVertexAttribPointer(GL2_Texture::a_coords, 4, GL_FLOAT, GL_FALSE, 4<<2, vrtx);	GL_CHECK();//select vertices & texcoord

	glBindBuffer(GL_ARRAY_BUFFER, 0);			GL_CHECK();
	glVertexAttribPointer(GL2_Texture::a_coords, 4, GL_FLOAT, GL_FALSE, 4<<2, vrtx);	GL_CHECK();//select vertices & texcoord

	glEnableVertexAttribArray(GL2_Texture::a_coords);	GL_CHECK();
	glDrawArrays(GL_TRIANGLES, 0, 6);					GL_CHECK();//draw the quad
	glDisableVertexAttribArray(GL2_Texture::a_coords);	GL_CHECK();
}
void			display_texture(int x1, int x2, int y1, int y2, int *rgb, int txw, int txh, unsigned char alpha)
{
	static unsigned tx_id=0;
	float _2_w=2.f/w, _2_h=2.f/h;
	float rect[]=
	{
		x1*_2_w-1, 1-y1*_2_h,
		x2*_2_w-1, 1-y2*_2_h//y2<y1
	};
	float vrtx[]=
	{
		rect[0], rect[1],		0, 0,//top left
		rect[0], rect[3],		0, 1,//bottom left
		rect[2], rect[3],		1, 1,//bottom right
		rect[2], rect[1],		1, 0 //top right
	};
	if(rgb)
	{
	#define	NPOT_ATIX2300_FIX
		int *rgb2, w2, h2;
#ifdef NPOT_ATIX2300_FIX
		bool expand;
		int logw=floor_log2(txw), logh=floor_log2(txh);
		if(expand=glMajorVer<3&&(txw>1<<logw||txh>1<<logh))
		{
			w2=txw>1<<logw?1<<(logw+1):txw;
			h2=txh>1<<logh?1<<(logh+1):txh;
			int size=w2*h2;
			rgb2=(int*)malloc(size<<2);
			memset(rgb2, 0, size<<2);
			for(int ky=0;ky<txh;++ky)
				memcpy(rgb2+w2*ky, rgb+txw*ky, txw<<2);
		//	memcpy(rgb2, rgb, size<<2);
			float nw=(float)txw/w2, nh=(float)txh/h2;
			vrtx[ 2]=0,		vrtx[ 3]=0;
			vrtx[ 6]=0,		vrtx[ 7]=nh;
			vrtx[10]=nw,	vrtx[11]=nh;
			vrtx[14]=nw,	vrtx[15]=0;
		}
		else
#endif
			rgb2=rgb, w2=txw, h2=txh;
		gl_setProgram(GL2_Text::program);				GL_CHECK();//select Text program
		glUniform2f(GL2_Text::u_isTexture, 1, 0);		GL_CHECK();//send isTexture
	//	glUniform1i(GL2_Text::u_isTexture, true);		GL_CHECK();
		send_color(GL2_Text::u_bkColor, alpha<<24);		GL_CHECK();//send apha
		if(!tx_id)
			{glGenTextures(1, &tx_id);															GL_CHECK();}//generate texture id once
		glBindTexture(GL_TEXTURE_2D, tx_id);													GL_CHECK();
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);						GL_CHECK();
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);						GL_CHECK();
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w2, h2, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb2);	GL_CHECK();//send bitmap to GPU

		select_texture(tx_id, GL2_Text::u_mytexture);
		glBindBuffer(GL_ARRAY_BUFFER, GL2_Text::buffer);										GL_CHECK();
		glBufferData(GL_ARRAY_BUFFER, 16<<2, vrtx, GL_STATIC_DRAW);								GL_CHECK();//send vertices & texcoords
		glVertexAttribPointer(GL2_Text::a_coord2d, 4, GL_FLOAT, GL_FALSE, 4<<2, (void*)0);		GL_CHECK();//select vertices & texcoord

		glEnableVertexAttribArray(GL2_Text::a_coord2d);		GL_CHECK();
		glDrawArrays(GL_QUADS, 0, 4);						GL_CHECK();//draw the quad
		glDisableVertexAttribArray(GL2_Text::a_coord2d);	GL_CHECK();
#ifdef NPOT_ATIX2300_FIX
		if(expand)
			free(rgb2);
#endif
	}
}
#if 0
namespace		GL2_TI2D
{
	unsigned	program=0;
	int			a_coords=-1, u_ndr=-1, u_invs=-1, u_curvecolor=-1;
	void		render(float *ndr, int curvecolor)
	{
		static unsigned tx_id=0;
		int x1=0, x2=w, y1=0, y2=h, txw=w+2, txh=h+2;
		float _2_w=2.f/w, _2_h=2.f/h;
		float rect[]=
		{
			x1*_2_w-1, 1-y1*_2_h,
			x2*_2_w-1, 1-y2*_2_h//y2<y1
		};
		float txrect[]=
		{
			1.f/txw, 1.f/txh,
			1-1.f/txw, 1-1.f/txh,
		};
		float vrtx[]=
		{
			rect[0], rect[1],		txrect[0], txrect[1],//top left
			rect[0], rect[3],		txrect[0], txrect[3],//bottom left
			rect[2], rect[3],		txrect[2], txrect[3],//bottom right
			rect[2], rect[1],		txrect[2], txrect[1] //top right
		};
		if(ndr)
		{
	#define	NPOT_ATIX2300_FIX
			float *ndr2;
			int txw2, txh2;
#ifdef NPOT_ATIX2300_FIX
			bool expand;
			int logw=floor_log2(txw), logh=floor_log2(txh);
			if(expand=glMajorVer<3&&(txw>1<<logw||txh>1<<logh))
			{
				txw2=txw>1<<logw?1<<(logw+1):txw;
				txh2=txh>1<<logh?1<<(logh+1):txh;
				int size=txw2*txh2;
				ndr2=(float*)malloc(size<<2);
				memset(ndr2, 0, size<<2);
				for(int ky=0;ky<txh;++ky)
					memcpy(ndr2+txw2*ky, ndr+txw*ky, txw<<2);
				//for(int ky=0;ky<h;++ky)
				//	for(int kx=0;kx<w;++kx)
				//		SetPixel(ghDC, kx, ky, (int)ndr2[txw2*ky+kx]&0xFFFFFF);
				float
					txx1=1.f/txw2, txy1=1.f/txh2,
					txx2=((float)txw-1)/txw2, txy2=((float)txh-1)/txh2;
				vrtx[ 2]=txx1,	vrtx[ 3]=txy1;
				vrtx[ 6]=txx1,	vrtx[ 7]=txy2;
				vrtx[10]=txx2,	vrtx[11]=txy2;
				vrtx[14]=txx2,	vrtx[15]=txy1;
			}
			else
#endif
				ndr2=ndr, txw2=txw, txh2=txh;
			gl_setProgram(GL2_TI2D::program);					GL_CHECK();//select TI2D program
			send_color(GL2_TI2D::u_curvecolor, curvecolor);		GL_CHECK();//send color
			glUniform2f(GL2_TI2D::u_invs, 1.f/txw2, 1.f/txh2);	GL_CHECK();
			if(!tx_id)
				{glGenTextures(1, &tx_id);													GL_CHECK();}//generate texture id once
			glBindTexture(GL_TEXTURE_2D, tx_id);											GL_CHECK();
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);				GL_CHECK();
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);				GL_CHECK();
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txw2, txh2, 0, GL_RED, GL_FLOAT, ndr2);	GL_CHECK();//send ndr to GPU

			select_texture(tx_id, GL2_TI2D::u_ndr);											GL_CHECK();
			glBindBuffer(GL_ARRAY_BUFFER, 0);												GL_CHECK();
			glVertexAttribPointer(GL2_TI2D::a_coords, 4, GL_FLOAT, GL_FALSE, 4<<2, vrtx);	GL_CHECK();//send vertices & texcoords

			glEnableVertexAttribArray(GL2_TI2D::a_coords);		GL_CHECK();
			glDrawArrays(GL_QUADS, 0, 4);						GL_CHECK();//draw the quad
			glDisableVertexAttribArray(GL2_TI2D::a_coords);		GL_CHECK();
#ifdef NPOT_ATIX2300_FIX
			if(expand)
				free(ndr2);
#endif
		}
	}
}
#endif
void			gl_initiate(HDC ghDC, int w, int h)
{
	tagPIXELFORMATDESCRIPTOR pfd={sizeof(tagPIXELFORMATDESCRIPTOR), 1, PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER, PFD_TYPE_RGBA, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, PFD_MAIN_PLANE, 0, 0, 0, 0};
	int PixelFormat=ChoosePixelFormat(ghDC, &pfd);
//	prof_add("ChoosePixelFormat");//either this takes 60ms
	SetPixelFormat(ghDC, PixelFormat, &pfd);
//	prof_add("SetPixelFormat");//...or this 60ms
	hRC=wglCreateContext(ghDC);
//	prof_add("wglCreateContext");//37.5ms
	wglMakeCurrent(ghDC, hRC);
//	prof_add("wglMakeCurrent");
	if(h<=0)
		h=1;
	glViewport(0, 0, w, h);
//	prof_add("glViewport");
	glShadeModel(GL_SMOOTH);
//	prof_add("glShadeModel");
	glClearColor(1, 1, 1, 1);
//	glClearColor(0, 0, 0, 1);
//	prof_add("glClearColor");
	glClearDepth(1);
//	prof_add("glClearDepth");
	glEnable(GL_DEPTH_TEST);
//	prof_add("glEnable");
	glDepthFunc(GL_LEQUAL);
//	prof_add("glDepthFunc");
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
//	prof_add("glHint");

	glGetIntegerv(GL_MAJOR_VERSION, &glMajorVer);
	glGetIntegerv(GL_MINOR_VERSION, &glMinorVer);
	GLversion=glGetString(GL_VERSION);
	if(API_not_loaded)
	{
	//	glGenVertexArrays=(decltype(glGenVertexArrays))wglGetProcAddress("glGenVertexArrays");//OpenGL 3.0
	//	glDeleteVertexArrays=(decltype(glDeleteVertexArrays))wglGetProcAddress("glDeleteVertexArrays");//OpenGL 3.0
		glBindVertexArray=(decltype(glBindVertexArray))wglGetProcAddress("glBindVertexArray");//not defined in GL ES 2
		glGenBuffers=(decltype(glGenBuffers))wglGetProcAddress("glGenBuffers");
		glBindBuffer=(decltype(glBindBuffer))wglGetProcAddress("glBindBuffer");
		glBufferData=(decltype(glBufferData))wglGetProcAddress("glBufferData");
		glBufferSubData=(decltype(glBufferSubData))wglGetProcAddress("glBufferSubData");
	//	glGetBufferSubData=(decltype(glGetBufferSubData))wglGetProcAddress("glGetBufferSubData");
		glEnableVertexAttribArray=(decltype(glEnableVertexAttribArray))wglGetProcAddress("glEnableVertexAttribArray");
		glVertexAttribPointer=(decltype(glVertexAttribPointer))wglGetProcAddress("glVertexAttribPointer");
		glDisableVertexAttribArray=(decltype(glDisableVertexAttribArray))wglGetProcAddress("glDisableVertexAttribArray");
		glCreateShader=(decltype(glCreateShader))wglGetProcAddress("glCreateShader");
		glShaderSource=(decltype(glShaderSource))wglGetProcAddress("glShaderSource");
		glCompileShader=(decltype(glCompileShader))wglGetProcAddress("glCompileShader");
		glGetShaderiv=(decltype(glGetShaderiv))wglGetProcAddress("glGetShaderiv");
		glGetShaderInfoLog=(decltype(glGetShaderInfoLog))wglGetProcAddress("glGetShaderInfoLog");
		glCreateProgram=(decltype(glCreateProgram))wglGetProcAddress("glCreateProgram");
		glAttachShader=(decltype(glAttachShader))wglGetProcAddress("glAttachShader");
		glLinkProgram=(decltype(glLinkProgram))wglGetProcAddress("glLinkProgram");
		glGetProgramiv=(decltype(glGetProgramiv))wglGetProcAddress("glGetProgramiv");
		glGetProgramInfoLog=(decltype(glGetProgramInfoLog))wglGetProcAddress("glGetProgramInfoLog");
		glDetachShader=(decltype(glDetachShader))wglGetProcAddress("glDetachShader");
		glDeleteShader=(decltype(glDeleteShader))wglGetProcAddress("glDeleteShader");
		glUseProgram=(decltype(glUseProgram))wglGetProcAddress("glUseProgram");
		glGetAttribLocation=(decltype(glGetAttribLocation))wglGetProcAddress("glGetAttribLocation");
		glDeleteProgram=(decltype(glDeleteProgram))wglGetProcAddress("glDeleteProgram");
		glDeleteBuffers=(decltype(glDeleteBuffers))wglGetProcAddress("glDeleteBuffers");
		glGetUniformLocation=(decltype(glGetUniformLocation))wglGetProcAddress("glGetUniformLocation");
		glUniform1f=(decltype(glUniform1f))wglGetProcAddress("glUniform1f");
		glUniformMatrix3fv=(decltype(glUniformMatrix3fv))wglGetProcAddress("glUniformMatrix3fv");
		glUniformMatrix4fv=(decltype(glUniformMatrix4fv))wglGetProcAddress("glUniformMatrix4fv");
		glGetBufferParameteriv=(decltype(glGetBufferParameteriv))wglGetProcAddress("glGetBufferParameteriv");
		glActiveTexture=(decltype(glActiveTexture))wglGetProcAddress("glActiveTexture");
		glUniform1i=(decltype(glUniform1i))wglGetProcAddress("glUniform1i");
		glUniform2f=(decltype(glUniform2f))wglGetProcAddress("glUniform2f");
		glUniform3f=(decltype(glUniform3f))wglGetProcAddress("glUniform3f");
		glUniform3fv=(decltype(glUniform3fv))wglGetProcAddress("glUniform3fv");
		glUniform4f=(decltype(glUniform4f))wglGetProcAddress("glUniform4f");
		API_not_loaded=false;
	}
//	prof_add("Load API");
#if 1
	GL2_2D::vertex_buffer=make_gpu_buffer(GL_ARRAY_BUFFER, 0, 2<<2);//dummy size
	
	glGenBuffers(1, &GL2_L3D::VBO);
	glGenBuffers(1, &GL2_L3D::EBO);

	GL2_3D::vertex_buffer=make_gpu_buffer(GL_ARRAY_BUFFER, 0, 3<<2);//dummy sizes
	float texcoords[]=
	{
		0, 0,
		0, 1,
		1, 1,
		1, 0
	};
	GL2_3D::texcoord_buffer=make_gpu_buffer(GL_ARRAY_BUFFER, texcoords, sizeof(texcoords));
	
	using namespace resources;
	calculate_text_txcoords(sf10_widths, sf10_height, sf10_txcoords);
	calculate_text_txcoords(sf8_widths, sf8_height, sf8_txcoords);
	calculate_text_txcoords(sf7_widths, sf7_height, sf7_txcoords);
	calculate_text_txcoords(sf6_widths, sf6_height, sf6_txcoords);
	GL2_Text::buffer=make_gpu_buffer(GL_ARRAY_BUFFER, 0, 8<<2);		GL_CHECK();//dummy size
//	prof_add("Alloc bufs");
	
	//if(one_shot)
	//{
	ShaderVar _2d_attr[]=
	{
		{&GL2_2D::a_vertices, "a_vertex", __LINE__}
	};
	ShaderVar _2d_unif[]=
	{
		{&GL2_2D::u_color, "u_color", __LINE__}
	};
	GL2_2D::program=LoadShaders(				//2D program
		"#version 120\n"
		"attribute vec2 a_vertex;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(a_vertex, 0., 1.);\n"
		"}",
		"#version 120\n"
		"uniform vec4 u_color;\n"
		"void main()\n"
		"{\n"
		"    gl_FragColor=u_color;\n"
	//	"    gl_FragColor=vec4(u_color.rgb, 0.5);\n"
		"}",
		_2d_attr, sizeof(_2d_attr)/sizeof(ShaderVar), _2d_unif, sizeof(_2d_unif)/sizeof(ShaderVar), __LINE__);
	if(!GL2_2D::program)
		GL_ERROR();
		
	ShaderVar l3d_attr[]=
	{
		{&GL2_L3D::a_vertices, "a_vertex", __LINE__},
		{&GL2_L3D::a_normals, "a_normal", __LINE__},
	//	{&GL2_L3D::a_texcoord, "a_texcoord", __LINE__}
	};
	ShaderVar l3d_unif[]=
	{
		{&GL2_L3D::u_vpmatrix, "u_vpmatrix", __LINE__},
		{&GL2_L3D::u_modelmatrix, "u_modelmatrix", __LINE__},
		{&GL2_L3D::u_normalmatrix, "u_normalmatrix", __LINE__},
	//	{&GL2_L3D::u_pointSize, "u_pointSize", __LINE__},
	//	{&GL2_L3D::u_texture, "u_texture", __LINE__},
	//	{&GL2_L3D::u_isTextured, "u_isTextured", __LINE__},
		{&GL2_L3D::u_objectcolor, "u_objectcolor", __LINE__},
		{&GL2_L3D::u_lightcolor, "u_lightcolor", __LINE__},
		{&GL2_L3D::u_lightpos, "u_lightpos", __LINE__},
		{&GL2_L3D::u_viewpos, "u_viewpos", __LINE__},
	};
	GL2_L3D::program=LoadShaders(						//3D program
		"#version 120\n"
		"uniform mat4 u_vpmatrix, u_modelmatrix;\n"
		"uniform mat3 u_normalmatrix;\n"
	//	"uniform float u_pointSize;\n"
		"attribute vec3 a_vertex;\n"
		"attribute vec3 a_normal;\n"
	//	"attribute vec2 a_texcoord;\n"
		"varying vec3 v_fragpos;\n"
		"varying vec3 v_normal;\n"
	//	"varying vec2 v_texcoord;\n"
		"varying vec4 v_glposition;\n"
	//	"varying vec4 v_noisecolor;"//CPU V5 DEBUG
	//	"float rand(vec2 co){return fract(sin(dot(co.xy, vec2(12.9898,78.233)))*43758.5453);}\n"//CPU V5 DEBUG
		"void main()\n"
		"{\n"
		"    vec4 fullpos=vec4(a_vertex, 1.);\n"
		"    gl_Position=u_vpmatrix*fullpos;\n"
		"    v_glposition=gl_Position;\n"
		"    gl_Position.z=0.;\n"
		"    v_fragpos=vec3(u_modelmatrix*fullpos);\n"
		"    v_normal=u_normalmatrix*a_normal;\n"
	//	"    v_noisecolor=vec4(rand(v_fragpos.xy), rand(v_fragpos.yz), rand(v_fragpos.xz), 1);\n"//CPU V5 DEBUG
	//	"    v_texcoord=a_texcoord;\n"
	//	"    gl_PointSize=u_pointSize;\n"
	//	"    gl_PointSize=2.;\n"
		"}",
		"#version 120\n"
		"varying vec3 v_fragpos;\n"
		"varying vec3 v_normal;\n"
		"varying vec4 v_glposition;\n"
	//	"varying vec4 v_noisecolor;"//CPU V5 DEBUG
	//	"varying vec2 v_texcoord;\n"
	//	"uniform bool u_isTextured;\n"
		"uniform vec4 u_objectcolor;\n"
		"uniform vec3 u_lightcolor, u_lightpos, u_viewpos;\n"
	//	"uniform sampler2D u_texture;\n"
		"void main()\n"
		"{\n"
		"    vec3 normal=normalize(v_normal);\n"
		"    vec3 lightdir=normalize(u_lightpos-v_fragpos);\n"
				
		"    float specularstrength=0.5;\n"
		"    vec3 viewdir=normalize(u_viewpos-v_fragpos), reflectdir=reflect(-lightdir, normal);\n"
		"    vec3 specular=specularstrength*u_lightcolor*pow(abs(dot(viewdir, reflectdir)), 32);\n"
	//	"    vec3 specular=specularstrength*u_lightcolor*pow(max(dot(viewdir, reflectdir), 0.), 32);\n"

		"    vec3 diffuse=abs(dot(normal, lightdir))*u_lightcolor;\n"
	//	"    vec3 diffuse=max(dot(normal, lightdir), 0.)*u_lightcolor;\n"
		"    gl_FragColor=vec4((0.1*u_lightcolor+diffuse+specular)*u_objectcolor.rgb, u_objectcolor.a);\n"
	//	"    gl_FragColor*=v_noisecolor;"//CPU V5 DEBUG

		"    gl_FragDepth=(-(1000.+0.1)*(-v_glposition.w)-2.*1000.*0.1)/((1000.-0.1)*v_glposition.w);\n"
		"}",
		l3d_attr, sizeof(l3d_attr)/sizeof(ShaderVar), l3d_unif, sizeof(l3d_unif)/sizeof(ShaderVar), __LINE__);
		
	glEnable(GL_PROGRAM_POINT_SIZE);
	ShaderVar _3d_attr[]=
	{
		{&GL2_3D::a_vertices, "a_vertex", __LINE__},
	//	{&GL2_3D::a_normals, "a_normal", __LINE__},
		{&GL2_3D::a_texcoord, "a_texcoord", __LINE__}
	};
	ShaderVar _3d_unif[]=
	{
		{&GL2_3D::u_matrix, "u_matrix", __LINE__},
	//	{&GL2_3D::u_pointSize, "u_pointSize", __LINE__},
		{&GL2_3D::u_texture, "u_texture", __LINE__},
		{&GL2_3D::u_isTextured, "u_isTextured", __LINE__},
		{&GL2_3D::u_color, "u_color", __LINE__}
	};
	GL2_3D::program=LoadShaders(						//3D program
		"#version 120\n"
		"uniform mat4 u_matrix;\n"
	//	"uniform float u_pointSize;\n"//can't get u_pointSize, -1
	//	"uniform u_depthcoeff;\n"//
	//	"uniform u_zfarcoeff;\n"//2/log(C*F+1)
		"attribute vec3 a_vertex;\n"
		"attribute vec2 a_texcoord;\n"
		"varying vec2 f_texcoord;\n"
	//	"noperspective varying vec4 v_fragpos;\n"
		"varying vec4 v_fragpos;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=u_matrix*vec4(a_vertex, 1.);\n"
		"    v_fragpos=gl_Position;\n"
	//	"    gl_Position.z=(1.-2.*exp(-0.01*gl_Position.w))*gl_Position.w;\n"//near clipping, no depth test
		"    gl_Position.z=0.;\n"//overlays text, no depth test
	//	"    gl_Position.z=gl_Position.w;\n"//depth=1, broken
	//	"    gl_Position.z=(log2(max(1e-6, 1.+gl_Position.w))*0.200657630123224-1.)*gl_Position.w;\n"//broken
	//	"    gl_Position.z=(2*log2(0.001*gl_Position.w+1.)/log2(0.001*1000.+1.)-1.)*gl_Position.w;\n"//broken
		"    f_texcoord=a_texcoord;\n"
	//	"    gl_PointSize=u_pointSize;\n"
	//	"    gl_PointSize=2.;\n"
		"}",
		"#version 120\n"
		"varying vec2 f_texcoord;\n"
	//	"noperspective varying vec4 v_fragpos;\n"
		"varying vec4 v_fragpos;\n"
	//	"uniform float u_isTextured;\n"
		"uniform bool u_isTextured;\n"
		"uniform vec4 u_color;\n"
		"uniform sampler2D u_texture;\n"
		"uniform int u_debug;\n"
		"void main()\n"
		"{\n"
		//"    gl_FragColor=mix(u_color, texture2D(u_texture, f_texcoord), u_isTextured);\n"
		"    if(u_isTextured)\n"
		"	     gl_FragColor=texture2D(u_texture, f_texcoord);\n"//alpha is in the texture
		"    else\n"
		"        gl_FragColor=u_color;\n"//alpha is in the color
		
	//	"    gl_FragColor.rgb=vec3(gl_FragCoord.z/gl_FragCoord.w);\n"//1 for most distance
	//	"    gl_FragColor.rgb=vec3(gl_FragCoord.z)-vec3((-(1000.+0.1)*(-v_fragpos.w)-2.*1000.*0.1)/((1000.-0.1)*v_fragpos.w));\n"//larger difference close to cam
	//	"    gl_FragColor.rgb=vec3(gl_FragCoord.z);\n"
	//	"    gl_FragColor.rgb=vec3((-(1000.+0.1)*(-v_fragpos.w)-2.*1000.*0.1)/((1000.-0.1)*v_fragpos.w));\n"

	//	"    gl_FragDepth=gl_FragCoord.z;\n"					//works
		"    gl_FragDepth=(-(1000.+0.1)*(-v_fragpos.w)-2.*1000.*0.1)/((1000.-0.1)*v_fragpos.w);\n"//works, L3D needs similar modification
	//	"    gl_FragDepth=1.-2.*exp(0.01/gl_FragCoord.w);\n"	//draws from 0 to inf, broken: no depth test
	//	"    gl_FragDepth=1.-2.*exp(-0.01/gl_FragCoord.w);\n"	//draws from 0 to inf, broken: no depth test
	//	"    gl_FragDepth=v_fragpos.w;\n"						//broken
	//	"    gl_FragDepth=v_fragpos.z;\n"						//broken
	//	"    gl_FragDepth=1.;\n"								//broken
	//	"    gl_FragDepth=-1./v_fragpos.w;\n"					//broken
	//	"    gl_FragDepth=1.-2.*exp(-0.01*v_fragpos.w);\n"		//broken
	//	"    gl_FragDepth=v_fragpos.z/v_fragpos.w;\n"			//broken
		"}",
		_3d_attr, sizeof(_3d_attr)/sizeof(ShaderVar), _3d_unif, sizeof(_3d_unif)/sizeof(ShaderVar), __LINE__);
	if(!GL2_3D::program)
		GL_ERROR();
#if 0
	ShaderVar ti2d_attr[]=
	{
		{&GL2_TI2D::a_coords, "coords", __LINE__},
	};
	ShaderVar ti2d_unif[]=
	{
		{&GL2_TI2D::u_ndr, "u_ndr", __LINE__},
		{&GL2_TI2D::u_invs, "u_invs", __LINE__},
		{&GL2_TI2D::u_curvecolor, "u_curvecolor", __LINE__}
	};
	GL2_TI2D::program=LoadShaders(					//TI2D program
		"#version 130\n"
		"attribute vec4 coords;"
		"varying vec2 f_texcoord;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(coords.xy, 0., 1.);\n"
		"    f_texcoord=coords.zw;\n"
		"}",
		"#version 130\n"
		"precision lowp float;\n"
		"varying highp vec2 f_texcoord;\n"
		"uniform vec2 u_invs;\n"//x: 1/w, y: 1/h
		"uniform sampler2D u_ndr;\n"
		"uniform vec4 u_curvecolor;\n"

		"float alpha_from_line(float x1, float y1, float x2, float y2)\n"
		"{\n"
		"	float a=x2-x1, b=y1-y2;\n"
		"	return max(0., 1.-abs(a*y1+b*x1)*inversesqrt(a*a+b*b));\n"
		"}\n"
		"float do_quadrant(float m, float R, float U, float UR)\n"//9.2fps
		"{\n"
		"	bool down=(m>0.)!=(R>0.), right=(R>0.)!=(UR>0.), up=(UR>0.)!=(U>0.), left=(U>0.)!=(m>0.);\n"
	//	"	return float(down||right||up||left);\n"
		"	float yL=m/(m-U), xD=m/(m-R), yR=R/(R-UR), xU=U/(U-UR);\n"
		"	if(left)\n"
		"	{\n"
		"		if(down)\n"//case 4 & 16
		"			return alpha_from_line(0, yL, xD, 0);\n"
		"		if(right)\n"//case 7
		"			return alpha_from_line(0, yL, 1, yR);\n"
		"		return alpha_from_line(0, yL, xU, 1);\n"//up	case 11
		"	}\n"
		"	if(down)\n"//cases 6 & 10
		"	{\n"
		"		if(right)\n"//	case 6
		"			return alpha_from_line(xD, 0, 1, yR);\n"
		"		return alpha_from_line(xD, 0, xU, 1);\n"//up	case 10
		"	}\n"
		"	if(up)\n"//&&right	case 13
		"		return alpha_from_line(xU, 1, 1, yR);\n"
		"	return 0.;\n"//case 1
		"}\n"
		"void main()\n"
		"{\n"
		"	float\n"
		"		xw=f_texcoord.x-u_invs.x, xm=f_texcoord.x, xe=f_texcoord.x+u_invs.x,\n"
		"		ys=f_texcoord.y-u_invs.y, ym=f_texcoord.y, yn=f_texcoord.y+u_invs.y;\n"
		"	float Vnw=texture2D(u_ndr, vec2(xw, yn)).r-0.5, Vnm=texture2D(u_ndr, vec2(xm, yn)).r-0.5, Vne=texture2D(u_ndr, vec2(xe, yn)).r-0.5;\n"//texture values are clamped in [0, 1]
		"	float Vmw=texture2D(u_ndr, vec2(xw, ym)).r-0.5, Vmm=texture2D(u_ndr, vec2(xm, ym)).r-0.5, Vme=texture2D(u_ndr, vec2(xe, ym)).r-0.5;\n"
		"	float Vsw=texture2D(u_ndr, vec2(xw, ys)).r-0.5, Vsm=texture2D(u_ndr, vec2(xm, ys)).r-0.5, Vse=texture2D(u_ndr, vec2(xe, ys)).r-0.5;\n"
		"	float alpha=float(Vmm==0.);\n"
		"	if(alpha!=1.)\n"
		"	{\n"
		"		alpha=do_quadrant(Vmm, Vme, Vnm, Vne);\n"
		"		alpha=max(alpha, do_quadrant(Vmm, Vnm, Vmw, Vnw));\n"
		"		alpha=max(alpha, do_quadrant(Vmm, Vmw, Vsm, Vsw));\n"
		"		alpha=max(alpha, do_quadrant(Vmm, Vsm, Vme, Vse));\n"
		"	}\n"
		"	gl_FragColor=vec4(u_curvecolor.rgb, u_curvecolor.a*alpha);\n"
		"}",
		ti2d_attr, sizeof(ti2d_attr)/sizeof(ShaderVar), ti2d_unif, sizeof(ti2d_unif)/sizeof(ShaderVar), __LINE__);
	if(!GL2_TI2D::program)
		GL_ERROR();
#endif

	ShaderVar text_attr[]=
	{
		{&GL2_Text::a_coord2d, "coords", __LINE__},
	//	{&GL2_Text::a_coord2d, "coord2d", __LINE__},
	//	{&GL2_Text::a_texcoord, "texcoord", __LINE__}
	};
	ShaderVar text_unif[]=
	{
		{&GL2_Text::u_mytexture, "mytexture", __LINE__},
		{&GL2_Text::u_isTexture, "isTexture", __LINE__},
		{&GL2_Text::u_txtColor, "txtColor", __LINE__},
		{&GL2_Text::u_bkColor, "bkColor", __LINE__}
	};
	GL2_Text::program=LoadShaders(					//Text program
		"#version 120\n"
		"attribute vec4 coords;"
		//"attribute vec2 coord2d;\n"			//coord2d, texcoord
		//"attribute vec2 texcoord;\n"
		"varying vec2 f_texcoord;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(coords.xy, 0., 1.);\n"
		"    f_texcoord=coords.zw;\n"
		//"    gl_Position=vec4(coord2d, 0., 1.);\n"
		//"    f_texcoord=texcoord;\n"
		"}",
		"#version 120\n"
		"varying vec2 f_texcoord;\n"
		"uniform sampler2D mytexture;\n"	//mytexture, isTexture, txtColor, bkColor
		"uniform vec2 isTexture;\n"
	//	"uniform bool isTexture;\n"
		"uniform vec4 txtColor;\n"
		"uniform vec4 bkColor;\n"
		"void main()\n"
		"{\n"
		"    vec4 region=texture2D(mytexture, f_texcoord);\n"
		"    gl_FragColor=mix(txtColor*region.r+bkColor*region.g, vec4(region.rgb, bkColor.a), isTexture.x);\n"//36.53fps

	//	"    gl_FragColor=mix(txtColor*region.r+vec4(bkColor.a, 0., 0., 1.)*region.g, vec4(region.rgb, bkColor.a), isTexture.x);\n"
	//	"    gl_FragColor=mix(txtColor*region.r+vec4(bkColor.rgb, 0.)*region.g, vec4(region.rgb, bkColor.a), isTexture.x);\n"
	//	"    gl_FragColor=txtColor*region.r+bkColor*region.g+0.1*isTexture.x;\n"
		//"	if(isTexture.x==0.)\n"
		//"	{\n"
		//"		if(region.r==1.)\n"//inside
		//"			gl_FragColor=txtColor;\n"
		//"		else if(region.g==1.)\n"//bk
		//"			gl_FragColor=vec4(bkColor.rgb, 0.);\n"
		//"	}\n"
		//"	else\n"
		//"		gl_FragColor=vec4(region.rgb, bkColor.a);\n"
		"}",
		text_attr, sizeof(text_attr)/sizeof(ShaderVar), text_unif, sizeof(text_unif)/sizeof(ShaderVar), __LINE__);
	if(!GL2_Text::program)
		GL_ERROR();
	
	ShaderVar texture_attr[]=								//Texture program
	{
		{&GL2_Texture::a_coords, "coords", __LINE__},
	};
	ShaderVar texture_unif[]=
	{
		{&GL2_Texture::u_mytexture, "mytexture", __LINE__},
		{&GL2_Texture::u_alpha, "alpha", __LINE__},
	};
	GL2_Texture::program=LoadShaders(
		"#version 120\n"
		"attribute vec4 coords;"			//coords
		"varying vec2 f_texcoord;\n"
		"void main()\n"
		"{\n"
		"    gl_Position=vec4(coords.xy, 0., 1.);\n"
		"    f_texcoord=coords.zw;\n"
		"}",
		"#version 120\n"
	//	"precision lowp float;\n"//not supported in GLSL 1.3
		"varying vec2 f_texcoord;\n"
		"uniform sampler2D mytexture;\n"	//mytexture, alpha
  		"uniform float alpha;\n"
		"void main()\n"
		"{\n"
		"    vec4 color=texture2D(mytexture, f_texcoord);\n"
  		"    color.a*=alpha;\n"
		"    gl_FragColor=color;\n"
		"}",
		texture_attr, sizeof(texture_attr)/sizeof(ShaderVar), texture_unif, sizeof(texture_unif)/sizeof(ShaderVar), __LINE__);
	if(!GL2_Texture::program)
		GL_ERROR();
	//	LOGERROR("Texture program not compiled.");
//	prof_add("Compile shaders");//26.7ms
	//}
	
	{
		using namespace resources;
		LoadFontTexture(sf10, sizeof(sf10), sf10_height, sf10_widths, GL2_Text::tx_id_sf10);
		LoadFontTexture(sf8, sizeof(sf8), sf8_height, sf8_widths, GL2_Text::tx_id_sf8);
		LoadFontTexture(sf7, sizeof(sf7), sf7_height, sf7_widths, GL2_Text::tx_id_sf7);
		LoadFontTexture(sf6, sizeof(sf6), sf6_height, sf6_widths, GL2_Text::tx_id_sf6);
	}
	//LoadFontTexture(compr_sf10, GL2_Text::tx_id_sf10);
	//LoadFontTexture(compr_sf8, GL2_Text::tx_id_sf8);
	//LoadFontTexture(compr_sf7, GL2_Text::tx_id_sf7);
	//LoadFontTexture(compr_sf6, GL2_Text::tx_id_sf6);
//	prof_add("Unpack textures");
	current_program=0;
	gl_setProgram(GL2_Text::program);	GL_CHECK();
	gl_font_size=gl_setTextSize(3);
	gl_setTextColor(0xFF000000);		GL_CHECK();
	gl_setBkColor(0xFFFFFFFF);			GL_CHECK();
	
		GL2_2D::set_region(0, w, 0, h);
		GL2_2D::use_region();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	GL_CHECK();
//	prof_add("Final preps");
	one_shot=false;
#endif
}
void			gl_finish()
{
	glDeleteProgram(GL2_2D::program);
	glDeleteBuffers(1, &GL2_2D::vertex_buffer);

	glDeleteProgram(GL2_L3D::program);

	glDeleteProgram(GL2_3D::program);
	glDeleteBuffers(1, &GL2_3D::vertex_buffer);
	glDeleteBuffers(1, &GL2_3D::texcoord_buffer);

//	glDeleteProgram(GL2_TI2D::program);
	
	glDeleteProgram(GL2_Text::program);
	glDeleteBuffers(1, &GL2_Text::buffer);
	
	glDeleteProgram(GL2_Texture::program);

	wglMakeCurrent(0, 0);
	wglDeleteContext(hRC);
	current_program=0;
}

//
//binder
//
double		changeFont(unsigned wParam)
{
	//if(usingOpenGL)
	//{
	//	bool mw_forward=((short*)&wParam)[1]>0;
	//	gl_font_size=gl_setTextSize(gl_font_size+mw_forward-!mw_forward);
	//	if(mw_forward)
	//		return font.larger();
	//	return font.smaller();
	//}
	return font.change(wParam);
}
double		largerFont()
{
	double A=font.larger();
	//if(usingOpenGL)
	//	gl_setTextSize(font.size());
	return A;
}
double		smallerFont()
{
	double A=font.smaller();
	//if(usingOpenGL)
	//	gl_setTextSize(font.size());
	return A;
}
double		setFont(int newFont)
{
	double A=font.assign(newFont);
	//if(usingOpenGL)
	//	gl_setTextSize(font.size());
	return A;
}
void		selectFont	()
{
	if(usingOpenGL)
		gl_setTextSize(font.size());
	else
		font.hFont=(HFONT__*)SelectObject(ghMemDC, font.hFont);
}
void		deselectFont()
{
	if(usingOpenGL)
		gl_setTextSize(3);
	else
		font.hFont=(HFONT__*)SelectObject(ghMemDC, font.hFont);
}
inline int		gl_getFontHeight(){return pixel_y*gl_fontH;}
inline int		gl_getCharWidth(char c){return pixel_x*font_widths[c];}
int			getTextWidth(const char *a, int length)
{
	if(usingOpenGL)
	{
		int width=0;
		for(int k=0;k<length;++k)
			if(a[k]>0&&a[k]<127)
				width+=gl_getCharWidth(a[k]);
		return width;
	}
	return font.getTextW(a, 0, length);
}
int			getTextWidth(const char *a, int i, int f)
{
	if(usingOpenGL)
	{
		int width=0;
		for(int k=i;k<f;++k)
			if(a[k]>0&&a[k]<127)
				width+=gl_getCharWidth(a[k]);
		return width;
	}
	return font.getTextW(a, i, f);
}
int			getBkMode()
{
	if(usingOpenGL)
		return 1+(((unsigned char*)&bkColor)[3]==0xFF);
	return GetBkMode(ghMemDC);
}
int			setBkMode(int mode)
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_Text::program);
		unsigned char alpha=((unsigned char*)&bkColor)[3];
		((unsigned char*)&bkColor)[3]=-(mode==2);
		send_color(GL2_Text::u_bkColor, bkColor);
		return 1+(alpha==0xFF);
	}
	return SetBkMode(ghMemDC, mode);
}
int			getBkColor()
{
	if(usingOpenGL)
		return bkColor&0xFFFFFF;
	return GetBkColor(ghMemDC);
}
int			setBkColor(int color)
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_Text::program);
		int oldColor=bkColor&0xFFFFFF;
		((unsigned char*)&bkColor)[0]=((unsigned char*)&color)[0];
		((unsigned char*)&bkColor)[1]=((unsigned char*)&color)[1];
		((unsigned char*)&bkColor)[2]=((unsigned char*)&color)[2];
		send_color(GL2_Text::u_bkColor, bkColor);
		return oldColor;
	}
	return SetBkColor(ghMemDC, color);
}
int			getTextColor()
{
	if(usingOpenGL)
		return txtColor&0xFFFFFF;
	return GetTextColor(ghMemDC);
}
int			setTextColor(int color)
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_Text::program);
		int oldColor=txtColor;
		txtColor=0xFF000000|0xFFFFFF&color;
		send_color(GL2_Text::u_txtColor, txtColor);
		return oldColor;
	}
	return SetTextColor(ghMemDC, color);
}
void		GUIPrint(int x, int y, const char* format, ...)
{
	int linelen=vsprintf_s(g_buf, format, (char*)(&format+1));
	if(linelen>0)
	{
		if(usingOpenGL)
			print_array(x, y, g_buf, linelen, 0);
		else
			TextOutA(ghMemDC, x, y, g_buf, linelen);
	}
}
void		GUIPrint(int x, int y, int value)
{
	int linelen=sprintf_s(g_buf, 1024, "%d", value);
	if(linelen>0)
	{
		if(usingOpenGL)
			print_array(x, y, g_buf, linelen, 0);
		else
			TextOutA(ghMemDC, x, y, g_buf, linelen);
	}
}
int			print(int x, int y, int tab_origin, const char *a, int length)
{
	if(usingOpenGL)
		return print_array(x, y, a, length, tab_origin);
	return TabbedTextOutA(ghMemDC, x, y, a, length, 0, 0, tab_origin);
}
int			print(int x, int y, const char *a, int length)
{
	if(usingOpenGL)
		return print_array(x, y, a, length, 0);
	return TextOutA(ghMemDC, x, y, a, length);
}

//2D API:
PenBrush::PenBrush(int color):pen_color(0xFF000000|color), brush_color(pen_color)
{
	hPen=CreatePen(PS_SOLID, 1, color);
	hBrush=CreateSolidBrush(color);
}
PenBrush::~PenBrush()
{
	DeleteObject(hPen);
	DeleteObject(hBrush);
}
void		PenBrush::use()
{
	if(usingOpenGL)
	{
		std::swap(pen_color, ::pen_color);
		std::swap(brush_color, ::brush_color);
	}
	else
	{
		hPen=(HPEN__*)SelectObject(ghMemDC, hPen);
		hBrush=(HBRUSH__*)SelectObject(ghMemDC, hBrush);
	}
}
void		PenBrush::drop()
{
	if(usingOpenGL)
	{
		std::swap(pen_color, ::pen_color);
		std::swap(brush_color, ::brush_color);
	}
	else
	{
		hPen=(HPEN__*)SelectObject(ghMemDC, hPen);
		hBrush=(HBRUSH__*)SelectObject(ghMemDC, hBrush);
	}
}

Pen::Pen(int color):color(0xFF000000|color){hPen=CreatePen(PS_SOLID, 1, color);}
Pen::~Pen()
{
	if(hPen)
		DeleteObject(hPen), hPen=nullptr;
}
void		Pen::set(int color){this->color=0xFF000000|color, hPen=CreatePen(PS_SOLID, 1, color);}
void		Pen::destroy(){DeleteObject(hPen), hPen=nullptr;}
void		Pen::use()
{
	if(usingOpenGL)
		std::swap(color, pen_color);
	else
		hPen=(HPEN__*)SelectObject(ghMemDC, hPen);
}
void		Pen::drop()
{
	if(usingOpenGL)
		std::swap(color, pen_color);
	else
		hPen=(HPEN__*)SelectObject(ghMemDC, hPen);
}

void		line(int x1, int y1, int x2, int y2)
{
	if(usingOpenGL)
		GL2_2D::draw_line((float)x1, (float)y1, (float)x2, (float)y2);
	else
		MoveToEx(ghMemDC, x1, y1, 0), LineTo(ghMemDC, x2, y2);
}
float		__x1, __y1;
void		moveTo(int x, int y)
{
	if(usingOpenGL)
		__x1=(float)x, __y1=(float)y;
	else
		MoveToEx(ghMemDC, x, y, 0);
}
void		lineTo(int x, int y)
{
	if(usingOpenGL)
		GL2_2D::draw_line(__x1, __y1, (float)x, (float)y);
	else
		LineTo(ghMemDC, x, y);
}
void		rectangle(int x1, int y1, int x2, int y2)
{
	if(usingOpenGL)
		GL2_2D::draw_rectangle((float)x1, (float)x2, (float)y1, (float)y2, true);
	else
		Rectangle(ghMemDC, x1, y1, x2, y2);
}
void		setPixel(int x, int y, int color)
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_2D::program);
		send_color(GL2_2D::u_color, color);
		GL2_2D::set_pixel((float)x, (float)y);
	}
	else if(x>=0&&x<w&&y>=0&&y<h)
		gBitmap.rgb[w*y+x]=color;
}

inline void	vertical_line_mul(int x, int y1, int y2, double Mr, double Mg, double Mb)
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_2D::program);
		glUniform4f(GL2_2D::u_color, (float)Mr, (float)Mg, (float)Mb, 0.5f);
		GL2_2D::draw_line((float)x, (float)y1, (float)x, (float)y2);
	}
	else
	{
		for(int y=y1;y<y2;++y)
		{
			auto px=(unsigned char*)&gBitmap.rgb[w*y+x];
			px[0]=unsigned char(Mb*px[0]);
			px[1]=unsigned char(Mg*px[1]);
			px[2]=unsigned char(Mr*px[2]);
		}
	}
}
void		vertical_line_equation(int x, int y1, int y2, double Ar, double Ag, double Ab, double a)
{
	vertical_line_mul(x, y1, y2, 1-Ar*a, 1-Ag*a, 1-Ab*a);
}
void		vertical_line_equation(int x, int y1, int y2, double a)
{
	vertical_line_mul(x, y1, y2, 1-a, 1-a, 1-a);
}
void		vertical_line_antiequation(int x, int y1, int y2, double Ar, double Ag, double Ab, double a, double Br, double Bg, double Bb)
{
	vertical_line_mul(x, y1, y2, Br+Ar*a, Bg+Ag*a, Bb+Ab*a);
}
void		vertical_line_antiequation(int x, int y1, int y2, double a)
{
	vertical_line_mul(x, y1, y2, 1-a, 1-a, 1-a);
}
void		vertical_line_logic_inequality(int x, int y1, int y2, double Ar, double Ag, double Ab)
{
	vertical_line_mul(x, y1, y2, 1-Ar, 1-Ag, 1-Ab);
}
void		vertical_line_logic_inequality(int x, int y1, int y2)
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_2D::program);
		pen_color=0x7F7F7F7F;
	//	glUniform4f(GL2_2D::u_color, 0.5f, 0.5f, 0.5f, 0.5f);
		GL2_2D::draw_line((float)x, (float)y1, (float)x, (float)y2);
	}
	else
	{
		for(int y=y1;y<y2;++y)
		{
			auto px=(unsigned char*)&gBitmap.rgb[w*y+x];
			px[0]>>=1;
			px[1]>>=1;
			px[2]>>=1;
		}
	}
//	vertical_line_mul(x, y1, y2, 0.5, 0.5, 0.5);
}
void		dim_screen()
{
	if(usingOpenGL)
	{
		gl_setProgram(GL2_2D::program);
		int pen_color2=0x7FFFFFFF, brush_color2=0x7FFFFFFF;
		std::swap(pen_color2, ::pen_color);
		std::swap(brush_color2, ::brush_color);
		GL2_2D::draw_rectangle(0, (float)w, 0, (float)h, false);
		std::swap(pen_color2, ::pen_color);
		std::swap(brush_color2, ::brush_color);
	}
	else
	{
		switch(simd_method)
		{
		case 0://ia32
			for(int k=0, kEnd=w*h;k<kEnd;++k)
			{
				auto p=(unsigned char*)(gBitmap.rgb+k);
				gBitmap.rgb[k]=~gBitmap.rgb[k];
				p[0]>>=1, p[1]>>=1, p[2]>>=1;
				gBitmap.rgb[k]=~gBitmap.rgb[k];
			}
			break;
		case 1://sse2
			lighten_sse2(gBitmap.rgb, w*h);
			break;
		case 2://avx2 required, won't compile on VS2010
			lighten_sse2(gBitmap.rgb, w*h);
			break;
		}
	}
}

//3D API:
void		newframe()
{
	if(usingOpenGL)
		glClear(GL_COLOR_BUFFER_BIT);
	else
		memset(gBitmap.rgb, 0xFF, w*h*sizeof(int));
}
void			clear_depth_buffer()
{
	if(usingOpenGL)
		glClear(GL_DEPTH_BUFFER_BIT);
	else
		memset(wbuffer, 0, w*h*sizeof(double));
	//	memset(wbuffer, 0, bw*bh*sizeof(double));
}
void			gl_line_3D(dvec3 const &p1, dvec3 const &p2, int color)
{
	GL2_3D::push_line_segment(p1, p2, color);
//	GL2_3D::push_line_segment(vec3(p1.x, p1.y, p1.z), vec3(p2.x, p2.y, p2.z), color);
}
void			point_3D(Camera const &cam, dvec3 &p_world, int color)
{
	if(usingOpenGL)
	{
		vec3 p=p_world;
		GL2_3D::push_point(p, color);
	}
	else
	{
		dvec3 cp;
		cam.world2cam(p_world, cp);
		const int bx1=0, bx2=w, by1=0, by2=h, bw=w, bh=h;
		if(cp.z>0)
		{
			double Acp=1/cp.z;
			dvec2 s;
			cam.cam2screen(cp, s);
			//double cpt=X0/(cp.z*cam.tanfov);
			//dvec2 s(X0+cp.x*cpt, Y0+cp.y*cpt);
			if((s.x>=bx1)&(s.x<bx2)&(s.y>=by1)&(s.y<by2))
			{
				int pos=w*int(s.y)+int(s.x);
				if(Acp>wbuffer[pos])
					gBitmap.rgb[pos]=color, wbuffer[pos]=Acp;
			}
		}
	}
}
void			point_3D_2x2(Camera const &cam, dvec3 &p_world, int Rcolor)
{
	if(usingOpenGL)
	{
		vec3 p=p_world;
		GL2_3D::push_point(p, Rcolor);
	}
	else
	{
		dvec3 cp;
		cam.world2cam(p_world, cp);
		if(cp.z>0)
		{
			double Acp=1/cp.z;
			dvec2 s;
			cam.cam2screen(cp, s);
		//	GUIPrint(ghMemDC, s.x, s.y, Rcolor);
			int xs1=int(s.x)-(s.x<0), ys1=int(s.y)-(s.y<0);
			_2dSet3dPoint(xs1, ys1	, Acp, Rcolor), _2dSet3dPoint(xs1+1, ys1	, Acp, Rcolor);
			_2dSet3dPoint(xs1, ys1+1, Acp, Rcolor), _2dSet3dPoint(xs1+1, ys1+1	, Acp, Rcolor);
		}
	}
}
void			point_3D_2x2(Camera const &cam, dvec3 &p_world, int Rcolor, int Icolor, int Jcolor, int Kcolor)
{
	if(usingOpenGL)
	{
		vec3 p=p_world;
		GL2_3D::push_point(p, Rcolor);
		GL2_3D::push_point(p+vec3(0.1f, 0   , 0), Icolor);
		GL2_3D::push_point(p+vec3(0   , 0.1f, 0), Jcolor);
		GL2_3D::push_point(p+vec3(0.1f, 0.1f, 0), Kcolor);
	}
	else
	{
		dvec3 cp;
		cam.world2cam(p_world, cp);
		if(cp.z>0)
		{
			double Acp=1/cp.z;
			dvec2 s;
			cam.cam2screen(cp, s);
			int xs1=int(s.x)-(s.x<0), ys1=int(s.y)-(s.y<0);
			_2dSet3dPoint(xs1		, ys1		, Acp, Rcolor);
			_2dSet3dPoint(xs1+1		, ys1		, Acp, Rcolor);
			_2dSet3dPoint(xs1+1		, ys1+1		, Acp, Rcolor);
			_2dSet3dPoint(xs1		, ys1+1		, Acp, Rcolor);

			_2dSet3dPoint(xs1+3		, ys1		, Acp, Icolor);
			_2dSet3dPoint(xs1+3+1	, ys1		, Acp, Icolor);
			_2dSet3dPoint(xs1+3+1	, ys1+1		, Acp, Icolor);
			_2dSet3dPoint(xs1+3		, ys1+1		, Acp, Icolor);

			_2dSet3dPoint(xs1		, ys1+3		, Acp, Jcolor);
			_2dSet3dPoint(xs1+1		, ys1+3		, Acp, Jcolor);
			_2dSet3dPoint(xs1+1		, ys1+3+1	, Acp, Jcolor);
			_2dSet3dPoint(xs1		, ys1+3+1	, Acp, Jcolor);

			_2dSet3dPoint(xs1+3		, ys1+3		, Acp, Kcolor);
			_2dSet3dPoint(xs1+3+1	, ys1+3		, Acp, Kcolor);
			_2dSet3dPoint(xs1+3+1	, ys1+3+1	, Acp, Kcolor);
			_2dSet3dPoint(xs1+3		, ys1+3+1	, Acp, Kcolor);
		}
	}
}
void			render_solid_transparent(Camera const &cam, dvec3 const &p1, dvec3 const &p2, dvec3 const &p3, int color)
//void			render_solid_transparent(Camera const &cam, vec3 const &p1, vec3 const &p2, vec3 const &p3, int color)
{
	if(usingOpenGL)
	{
		GL2_3D::push_triangle(p1, p2, p3, color);
		return;
	}
	vec2 s1, s2, s3;
	vec3 c1, c2, c3;
	float X0_tfov=float(X0/cam.tanfov);
	cam.world2camscreen(p1, c1, s1), cam.world2camscreen(p2, c2, s2), cam.world2camscreen(p3, c3, s3);
			if(c1.z<0){		 if(c2.z<0){		 if(c3.z<0)	return;
											else			draft_2behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_2behind(s2, s3, s1);
											else			draft_1behind(s2, s3, s1);}}
	else{					 if(c2.z<0){		 if(c3.z<0)	draft_2behind(s1, s2, s3);
											else			draft_1behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_1behind(s1, s2, s3);
											else			draft_start(s1, s2), draft(s2, s3), draft(s3, s1);}}
	vec3 upy=c2-c1,		//u12	=<12>
		upx=c3-c1,			//ux3	=<13>
		a=upy.cross(upx);	//abc	=<n>	=<12>x<13>
	float t=a.dot(c1);
	if(!t)return;
	float B7=a.x, B8=a.y, B9=a.z*X0_tfov-a.x*X0-a.y*Y0;	float cpt=t*X0_tfov;
	float A1=B7/cpt, A2=B8/cpt, A3=B9/cpt;

	double *wb_y=wbuffer;
	int lib_y, lfb_y
		, *rgb_y=gBitmap.rgb
		;
	__m128i const sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
	__m128i const sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
	const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
	const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
	__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
	__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
	__m128i m_color=_mm_set1_epi32(color);
	for(int ys=0;ys<h;++ys)
	{
	//	double *wb_y=wbuffer+w*ys;
	//	int *rgb_y=rgb+w*ys;
		lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
		float admittance=A1*lib_y+A2*ys+A3;

		int xs=lib_y;
		__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
		__m128 adm_increment=_mm_set1_ps(A1*4);
		for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)
		{
		//	__m128 wb=_mm_loadu_ps(wb_y+xs);
			__m128 wb=_mm_set_ps((float)wb_y[xs+3], (float)wb_y[xs+2], (float)wb_y[xs+1], (float)wb_y[xs]);
			__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
			__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
			__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
			c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
			c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
			if(c_or.m128i_i32[0])
			{
				c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));

				__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
				__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
				__m128i m_c2=_mm_and_si128(m_color, c_result);
				m_rgb=_mm_and_si128(m_rgb, c_result);
				
				__m128i zero=_mm_setzero_si128();
				__m128i m_c_lo=_mm_unpacklo_epi8(m_c2, zero);//{0, tx[7], ..., 0, tx[0]}
				__m128i m_c_hi=_mm_unpackhi_epi8(m_c2, zero);//{0, tx[15], ..., 0, tx[8]}
				__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
				__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
				if(lsw_transparency_multiply)
				{
					m_c_lo=_mm_mullo_epi16(m_c_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
					m_c_hi=_mm_mullo_epi16(m_c_hi, m_rgb_hi);
					m_c_lo=_mm_shuffle_epi8(m_c_lo, sh_mul_lo);
					m_c_hi=_mm_shuffle_epi8(m_c_hi, sh_mul_hi);
				}
				else
				{
					m_c_lo=_mm_add_epi16(m_c_lo, m_rgb_lo);
					m_c_hi=_mm_add_epi16(m_c_hi, m_rgb_hi);
					m_c_lo=_mm_srli_epi16(m_c_lo, 1);
					m_c_hi=_mm_srli_epi16(m_c_hi, 1);
					m_c_lo=_mm_shuffle_epi8(m_c_lo, sh_add_lo);
					m_c_hi=_mm_shuffle_epi8(m_c_hi, sh_add_hi);
				}
				m_c_lo=_mm_and_si128(m_c_lo, mask_lo);
				m_c_hi=_mm_and_si128(m_c_hi, mask_hi);
				m_c2=_mm_or_si128(m_c_lo, m_c_hi);

				m_c2=_mm_or_si128(m_c2, m_rgb2);
				_mm_storeu_si128((__m128i*)(rgb_y+xs), m_c2);
			}
			adm=_mm_add_ps(adm, adm_increment);
		}
		admittance=adm.m128_f32[0];
		for(;xs<lfb_y;++xs)
		{
			if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
				admittance>=wb_y[xs])
			{
				auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&color;
				if(lsw_transparency_multiply)
				{
					p[0]=p[0]*c[0]>>8;//b
					p[1]=p[1]*c[1]>>8;//g
					p[2]=p[2]*c[2]>>8;//r
				}
				else
				{
					p[0]=(p[0]+c[0])>>1;//b
					p[1]=(p[1]+c[1])>>1;//g
					p[2]=(p[2]+c[2])>>1;//r
				}
			}
			admittance+=A1;
		}
		wb_y=wb_y+w, rgb_y=rgb_y+w;
	}
}
void			render_solid_transparent(Camera const &cam, vec3 const &p1, vec3 const &p2, vec3 const &p3, int color)
{
	if(usingOpenGL)
		GL2_3D::push_triangle(p1, p2, p3, color);
	else
		render_solid_transparent(cam, dvec3(p1.x, p1.y, p1.z), dvec3(p2.x, p2.y, p2.z), dvec3(p3.x, p3.y, p3.z), color);
}

//basics
void			show()
{
	if(usingOpenGL)
		SwapBuffers(ghDC);
	else
		BitBlt(ghDC, 0, 0, w, h, ghMemDC, 0, 0, SRCCOPY);
}
void			resize_2D()
{
	if(usingOpenGL)
	{
		glViewport(0, 0, w, h);
		GL2_2D::set_region(0, w, 0, h);
		GL2_2D::use_region();
	}
	else
	{
		gBitmap.drop();
		gBitmap.resize(w, h);
		gBitmap.use();
	}
}

//CL-GL interop.
void			GPUBuffer::create_VN_I(float *VVVNNN, int n_floats, int *indices, int n_ints)
{
	n_vertices=n_ints;
	vertices_stride=6<<2, vertices_start=0;
	normals_stride=6<<2, normals_start=3<<2;
	if(!VBO)
		glGenBuffers(1, &VBO);
	if(!EBO)
		glGenBuffers(1, &EBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);											GL_CHECK();
	glBufferData(GL_ARRAY_BUFFER, n_floats<<2, VVVNNN, GL_STATIC_DRAW);			GL_CHECK();
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);									GL_CHECK();
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_ints<<2, indices, GL_STATIC_DRAW);	GL_CHECK();
}