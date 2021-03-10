Grapher 2: A handy graphing calculator.


How to compile:
	Dependencies on Windows:
	1) Download FFTW from:
		http://www.fftw.org/
	2) Download pre-compiled MPIR 2.7.2 and MPFR 4.0.0-dev libraries from:
		http://www.atelierweb.com/mpir-and-mpfr/
	3) Download MPFR C++ 3.6.2 from:
		http://www.holoborodko.com/pavel/mpfr/

	Option A: On an AVX-supporting compiler:
	4) Compile g2_avx.cpp & g2_common.cpp as a static library AVX.lib.
	5) Compile the other files including common.cpp and link to AVX.lib.
	This is to produce an executable that doesn't crash on older non-AVX CPUs.

	Option B: On a pre-AVX compiler:
	4) Comment out '#include"avx.h"' in g2.cpp.
	5) Compile everything excluding avx.cpp & avx.h

MSVC2013 project structure:
project: AVX
	Header Files:	g2_avx.h, g2_common.h
	Source Files:	g2_avx.cpp, g2_common.cpp
project: Grapher 2
	Header files:	g2_avx.h, g2_cl.h, g2_common.h, g2_error.h, g2_expr.h, g2_file.h, g2_graphics.h, g2_mp.h, g2_sse2.h
	Source Files:	g2.cpp, g2_cl.cpp, g2_common.cpp, g2_error.cpp, g2_expr.cpp, g2_file.cpp, g2_graphics.cpp, g2_mp.cpp, g2_sse2.cpp
	Resource Files:	g2_resources.rc


How to install:
	1) Install Microsoft Visual C++ 2013 Redistributable (x86) if not installed:
		https://www.microsoft.com/en-us/download/details.aspx?id=40784
	2) Download everything from bin folder.


How to use:
	Esc:	switch between edit and interaction modes
	F1:	context help
	F7:	toggle benchmark
	ctrl F7: switch between SIMD modes supported by CPU (IA32/SSE2/AVX).
	F11:	fullscreen


Text edit mode:
	F2/right click: parameter menu (only in edit mode)
	ctrl wheel/+/-: text size
	ctrl c:		copy selected text
	ctrl shift c:	copy with single number results


Parameter menu:
	number +/-: increment/decrement Nth parameter
	number 0: reset Nth parameter
	number enter: enter value for Nth parameter
	number space: select type of Nth parameter
	number S: select scale type for Nth parameter
	F2/right click: close parameter menu


Interaction mode:
	R:	reset view and scale
	C:	show/hide axis labels
	X/Y/Z drag/arrows:	drag x/y/z axis
	X/Y/Z wheel/+/-:	scale x/y/z axis
	Number keys:	apply operations on results (differentiate, FFT, LPF, ...)
	0:		remove all operations
	`:		show/hide contour
	TAB/shift TAB:	switch between {software, OpenGL, OpenCL-OpenGL interoperation} (OpenGL mode is best used with 3D graphs)

Numeric result:
	1:	binary
	2:	octal
	3/0:	decimal
	4:	hexadecimal
	+/-:	change numeric precision

2D:
	arrows/drag:	move
	wheel/+/-:	zoom
	F8:		reverse drag direction
	L/S:		switch between linear and log scale

3D:
	W/A/S/D/T/G:		move
	arrows/drag:		turn
	shift W/A/S/D/T/G:	move fast
	T/G:			move up/down
	wheel/+/-:		change fov
	alt wheel/+/-:		change scale
	ctrl wheel/+/-:		change speed
	L:		switch between linear and log scale




Space variables: x, y, z
Time variable: t
Complex variables: Z, W
Quaternion variable: Q




Operators:
	!	factorial, logic not
	%	modulo / percent
	~	bitwise not
	^^^	pentate
	^^	tetrate
	^ **	power, power a real integer
	* /
	@	divides
	+ -
	<< >>		bitwise shift
	< <= > >= == !=	comparison
	&, ~& 		bitwise and, nand
	#, ~# 		bitwise xor, xnor
	|, ~| 		bitwise or, nor
	&& ## ||	logical and, xor, or
	??	zero condition: x??y -> x?x:y


Functions - 1 argument:
	cos, acos, cosh, acosh, cosc,
	sec, asec, sech, asech,
	sin, asin, sinh, asinh, sinc, sinhc
	tan, atan, tanh, atanh, tanc,
	cot, acot, coth, acoth,
	exp, ln, log, sqrt, cbrt, invsqrt, sq,
	gauss, erf, fib, zeta, lngamma,
	step, sgn, rect, tent,
	ceil, floor, round, int, frac,
	abs, arg, real, imag, conjugate, polar, cartesian, rand,
	clamp		between 0 and 1


Functions - 2 arguments:
	rand, atan, log,
	beta, gamma, permutation, combination,
	bessel, neumann, hankel1
	sqwv, trwv, saw,
	mandelbrot(complex variable, number if iterations),
	clamp(x, hi)	between 0 and hi

Functions - 3 arguments
	clamp(lo, x, hi)

Functions - variadic
	min, max, av, hypot,
	norm(n, a1, ... am)		= (a1^n + ... am^n)^(1/n)

User-defined functions use C syntax (except switch and arrays are not supported)
Examples:
isprime(x)
{
	start=2;
	if(x<start)
		return false;
	if(x==start)
		return true;
	end=floor(sqrt(x));
	for(a=start;a<=end;++a)
		if(!(x%a))
			return false;
	return true;
}

gcd(a, b)
{
	while(b)
	{
		t=b;
		b=a%b;
		a=t;
	}
	return a;
}
gcd(a, b)//recursion and redefinition of the same function gcd(a, b)
{
	if(b)
		return gcd(b, a%b);
	return a;
}
gcd(a, b, c){return gcd(a, gcd(b, c));}//overloading

totient(x)//https://www.crazyforcode.com/eulers-totient-function/
{
	result=x;
	for(factor=2;factor*factor<=x;++factor)
	{
		if(!(x%factor))
		{
			while(!(x%factor))
				x/=factor;
			result-=result/factor;
		}
	}
	if(x>1)
		result-=result/x;
	return result;
}




CONSTANTS
i, j, k							Imaginary number & quaternions
_atm, atm?<		101325			Pa	Atmospheric pressure
_bbr			5.670373e-8		W/m2K4	Stefan–Boltzmann constant
_c, c<			299792458		m/s	Speed of light
_ele, permittivity	8.854187817620e-12	F/m	Electric permittivity in vacuum
_e, e			exp(1.)				Euler's	number
_g			9.80665			m/s2	Gravitational strength
_G			6.67384e-11		m3/kg/s	Gravitational constant
_h			6.62606957e-34		m2kg/s	Planck constant
_mag, permeability	1.2566370614e-6		H/m	Magnetic permeability in vacuum
_me			9.10938291e-31		kg	Electron mass
_Me			5.9736e24		kg	Earth mass
_mn			1.674927351e-27		kg	Neutron mass
_mp			1.672621777e-27		kg	Proton mass
_Ms			1.9891e30		kg	Sun mass
_Na			6.02214129e23		/mol	Avogadro constant
_phi			1.618033988749894848		Golden ratio
_pi, pi			::acos(-1.)			Pi
_q			1.602176565e-19		C	Elementary charge
_rand							Random
_R			8.3144621		J/mol.K	Gas constant
inf			1/+0				Positive infinity
ind, nan		0/0				Indeterminate quantity
























