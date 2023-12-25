#include <string.h>
#include <stdlib.h>

#ifndef NULL
#define NULL	0
#endif NULL

/*
 * Convert string of form "degrees[:minutes[:seconds]]" to angle
 * Also useful for times
 */
double
atoangle(p)
const	char	*p;
{
	double	sign;
	double	d, m, s;

	sign = 1.0;
	if (*p == '+')
		p++;
	if (*p == '-') {
		sign = -1.0;
		p++;
	}
	d = atof(p);
	m = s = 0.0;
	if ((p = strchr(p, ':')) != NULL) {
		m = atof(++p);
		if ((p = strchr(p, ':')) != NULL)
			s = atof(++p);
	}
	return sign*(d + (m + s/60.0)/60.0);
}
