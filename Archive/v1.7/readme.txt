How to use:

Text editor:
	Esc: switch between edit and interaction modes
	F1: context help
	F2/right click: parameter menu
	F7: toggle benchmark
	F11: fullscreen


Edit mode:
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
	`:	show/hide contour

Numeric result:
	1:	binary
	2:	octal
	3/0:	decimal
	4:	hexadecimal

2D:
	arrows:		move
	wheel/+/-:	zoom
	F8:		reverse drag direction

3D:
	W/A/S/D/T/G:		move
	arrows:			turn
	shift W/A/S/D/T/G:	move fast
	T/G:			move up/down
	wheel/+/-:		change fov
	alt wheel/+/-:		change scale
	ctrl wheel/+/-:		change speed




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
	^	power
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
	ceil, floor, round,
	abs, arg, real, imag, conjugate, polar, cartesian, rand


Functions - 2 arguments:
	rand, atan, log,
	beta, gamma, permutation, combination,
	bessel, neumann, hankel1
	sqwv, trwv, saw, min, max, hypot

User defined functions use C syntax (except switch is not supported)




CONSTANTS
i, j, k							imaginary number & quaternions
_atm, atm?<		101325			Pa	atmospheric pressure
_bbr			5.670373e-8		W/m2K4	Stefan–Boltzmann constant
_c, c<			299792458		m/s	speed of light
_ele, permittivity	8.854187817620e-12	F/m	electric permittivity in vacuum
_e, e			exp(1.)				euler's	number
_g			9.80665			m/s2	gravitational strength
_G			6.67384e-11		m3/kg/s	gravitational constant
_h			6.62606957e-34		m2kg/s	Planck constant
_mag, permeability	1.2566370614e-6		H/m	magnetic permeability in vacuum
_me			9.10938291e-31		kg	electron mass
_Me			5.9736e24		kg	earth mass
_mn			1.674927351e-27		kg	neutron mass
_mp			1.672621777e-27		kg	proton mass
_Ms			1.9891e30		kg	sun mass
_Na			6.02214129e23		/mol	Avogadro constant
_phi			1.618033988749894848		golden ratio
_pi, pi			::acos(-1.)			pi
_q			1.602176565e-19		C	elementary charge
_rand							random
_R			8.3144621		J/mol.K	gas constant
inf			1/+0				positive infinity
ind, nan		0/0				indeterminate quantity
























