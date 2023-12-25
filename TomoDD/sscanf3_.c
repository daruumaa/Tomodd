#include <stdio.h>
#include "compat.h"
#include "f77types.h"

/*
 * f77-callable interface to sscanf() function
 * Klewdgy version that allows three pointer arguments
 */
F77INT
sscanf3_(s, fmt, p1, p2, p3)
const	char	*s;
const	char	*fmt;
	GENPTR	p1, p2, p3;
{

	return sscanf(s, fmt, p1, p2, p3);
}
