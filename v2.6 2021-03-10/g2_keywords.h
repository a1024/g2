#ifdef TOKEN
	TOKEN(0, M_IGNORED),
		
	TOKEN(0, M_N),
		
	TOKEN("(", M_LPR), TOKEN(")", M_RPR),										//(		)
	TOKEN(",", M_COMMA),														//,
	TOKEN("?", M_QUESTION_MARK), TOKEN(":", M_COLON),							//?		:

		TOKEN(0, M_PROCEDURAL_START),
	TOKEN("if", M_IF), TOKEN("else", M_ELSE), TOKEN("for", M_FOR), TOKEN("do", M_DO), TOKEN("while", M_WHILE),					//if else for do while
	TOKEN("continue", M_CONTINUE), TOKEN("break", M_BREAK), TOKEN("return", M_RETURN),											//continue break return
	TOKEN("{", M_LBRACE), TOKEN("}", M_RBRACE),																					//{		}
	TOKEN(";", M_SEMICOLON),																									//;

		TOKEN(0, M_PROCEDURAL_ASSIGN_START),
	TOKEN("=", M_ASSIGN), TOKEN("*=", M_ASSIGN_MULTIPLY), TOKEN("/=", M_ASSIGN_DIVIDE), TOKEN("%=", M_ASSIGN_MOD),				//= *= /= %=			^= is confusing?
	TOKEN("+=", M_ASSIGN_PLUS), TOKEN("-=", M_ASSIGN_MINUS), TOKEN("<<=", M_ASSIGN_LEFT), TOKEN(">>=", M_ASSIGN_RIGHT),			//+= -= <<= >>=
	TOKEN("&=", M_ASSIGN_AND), TOKEN("#=", M_ASSIGN_XOR), TOKEN("|=", M_ASSIGN_OR),												//&= #= |=
		TOKEN(0, M_PROCEDURAL_ASSIGN_END),

		TOKEN(0, M_PROCEDURAL_END),

	TOKEN("++", M_INCREMENT), TOKEN("--", M_DECREMENT),							//++ --
	TOKEN("!", M_FACTORIAL_LOGIC_NOT),											//!
	TOKEN("%", M_MODULO_PERCENT),												//%
	TOKEN("~", M_BITWISE_NOT),													//~
	TOKEN("^^^", M_PENTATE),													//^^^
	TOKEN("^^", M_TETRATE),														//^^
	TOKEN("^", M_POWER), TOKEN("**", M_POWER_REAL),								//^		**
	TOKEN("*", M_MULTIPLY), TOKEN("/", M_DIVIDE), TOKEN("@", M_LOGIC_DIVIDES),	//*		/		@
	TOKEN("+", M_PLUS), TOKEN("-", M_MINUS),									//+		-
	TOKEN("<<", M_BITWISE_SHIFT_LEFT), TOKEN(">>", M_BITWISE_SHIFT_RIGHT),		//<<	>>
	TOKEN("<", M_LOGIC_LESS), TOKEN("<=", M_LOGIC_LESS_EQUAL), TOKEN(">", M_LOGIC_GREATER), TOKEN(">=", M_LOGIC_GREATER_EQUAL),	//<		<=		>		>=
	TOKEN("==", M_LOGIC_EQUAL), TOKEN("!=", M_LOGIC_NOT_EQUAL),					//==	!=
	TOKEN("&", M_BITWISE_AND), TOKEN("~&", M_BITWISE_NAND),						//&		~& |~
	TOKEN("#", M_BITWISE_XOR), TOKEN("~#", M_BITWISE_XNOR),						//#		~# #~
	TOKEN("|", M_VERTICAL_BAR), TOKEN("~|", M_BITWISE_NOR),						//|		~| &~
	TOKEN("&&", M_LOGIC_AND),													//&&
	TOKEN("##", M_LOGIC_XOR),													//##
	TOKEN("||", M_LOGIC_OR),													//||
	TOKEN("??", M_CONDITION_ZERO),												//??

	TOKEN("_=", M_S_EQUAL_ASSIGN), TOKEN("_!=", M_S_NOT_EQUAL),																	//_= _!=			not supported
	TOKEN("_<", M_S_LESS), TOKEN("_<=", M_S_LESS_EQUAL), TOKEN("_>", M_S_GREATER), TOKEN("_=>", M_S_GREATER_EQUAL),				//_< =< _> =>		not supported

		TOKEN(0, M_FSTART),//unary functions

	TOKEN("cos", M_COS), TOKEN("acos", M_ACOS), TOKEN("cosd", M_COSD), TOKEN("acosd", M_ACOSD), TOKEN("cosh", M_COSH), TOKEN("acosh", M_ACOSH), TOKEN("cosc", M_COSC),
	TOKEN("sec", M_SEC), TOKEN("asec", M_ASEC), TOKEN("secd", M_SECD), TOKEN("asecd", M_ASECD), TOKEN("sech", M_SECH), TOKEN("asech", M_ASECH),
	TOKEN("sin", M_SIN), TOKEN("asin", M_ASIN), TOKEN("sind", M_SIND), TOKEN("asind", M_ASIND), TOKEN("sinh", M_SINH), TOKEN("asinh", M_ASINH), TOKEN("sinc", M_SINC), TOKEN("sinhc", M_SINHC),
	TOKEN("csc", M_CSC), TOKEN("acsc", M_ACSC), TOKEN("cscd", M_CSCD), TOKEN("acscd", M_ACSCD), TOKEN("csch", M_CSCH), TOKEN("acsch", M_ACSCH),
	TOKEN("tan", M_TAN),						TOKEN("tand", M_TAND),							TOKEN("tanh", M_TANH), TOKEN("atanh", M_ATANH), TOKEN("tanc", M_TANC),
	TOKEN("cot", M_COT), TOKEN("acot", M_ACOT), TOKEN("cotd", M_COTD), TOKEN("acotd", M_ACOTD), TOKEN("coth", M_COTH), TOKEN("acoth", M_ACOTH),
	TOKEN("exp", M_EXP), TOKEN("ln", M_LN), TOKEN("sqrt", M_SQRT), TOKEN("cbrt", M_CBRT), TOKEN("invsqrt", M_INVSQRT), TOKEN("sq", M_SQ),
	TOKEN("gauss", M_GAUSS), TOKEN("erf", M_ERF), TOKEN("fib", M_FIB), TOKEN("zeta", M_ZETA), TOKEN("lngamma", M_LNGAMMA),
	TOKEN("step", M_STEP), TOKEN("sgn", M_SGN), TOKEN("rect", M_RECT), TOKEN("tent", M_TENT),
	TOKEN("ceil", M_CEIL), TOKEN("floor", M_FLOOR), TOKEN("round", M_ROUND), TOKEN("int", M_INT), TOKEN("frac", M_FRAC),
	TOKEN("abs", M_ABS), TOKEN("arg", M_ARG), TOKEN("re", M_REAL), TOKEN("im", M_IMAG), TOKEN("conj", M_CONJUGATE), TOKEN("polar", M_POLAR), TOKEN("cart", M_CARTESIAN),
	TOKEN("isprime", M_ISPRIME), TOKEN("totient", M_TOTIENT),
		
		TOKEN(0, M_BFSTART),//binary functions

	TOKEN("rand", M_RAND),
	TOKEN("atan", M_ATAN), TOKEN("atand", M_ATAND),
	TOKEN("log", M_LOG),
	TOKEN("beta", M_BETA), TOKEN("gamma", M_GAMMA), TOKEN("permutation", M_PERMUTATION), TOKEN("combination", M_COMBINATION),
	TOKEN("J", M_BESSEL_J), TOKEN("Y", M_BESSEL_Y), TOKEN("hankel", M_HANKEL1),
	TOKEN("sqwv", M_SQWV), TOKEN("trwv", M_TRWV), TOKEN("saw", M_SAW), TOKEN("mandelbrot", M_MANDELBROT),
	TOKEN("invmod", M_INVMOD),

		TOKEN(0, M_VFSTART),//variadic functions

	TOKEN("clamp", M_CLAMP),//up to 3 args
	TOKEN("min", M_MIN), TOKEN("max", M_MAX), TOKEN("av", M_AV), TOKEN("hypot", M_HYPOT), TOKEN("norm", M_NORM),
	TOKEN("gcd", M_GCD),

	TOKEN(0, M_USER_FUNCTION),
#endif