#include <math.h>

/* f77-callable interface to hypot function */
double
HYPOT(a, b)
	float	*a, *b;
{
	return hypot((double)*a, (double)*b);
}
