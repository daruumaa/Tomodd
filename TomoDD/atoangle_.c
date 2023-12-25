#include <string.h>
#include "compat.h"
#include "f77types.h"

	double	atoangle PARMS((char*));

/*
 * f77-callable interface to atoangle function
 */
double
atoangle_(p, l)
const	char	*p;	/* Character string		(input)	*/
	F77SLA	l;	/* Length of *p			(input)	*/
{
	char	q[l+1];	/* Buffer, with space for terminator	*/

	return atoangle(strncpy(q, p, l));
}
