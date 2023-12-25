#include <string.h>
#include "compat.h"
#include "f77types.h"

	double	atoangle PARMS((char*));

/*
 * f77-callable interface to atoangle function
 */
float
chtof_(p)
const	char	*p;	/* Character string		(input)	*/
{
	return atof(p);
}
