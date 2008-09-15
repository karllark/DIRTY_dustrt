#include "random_dirty.h"

/* ====================================================================== */
/*   Function to return a random number.  Function taken from "Numerical  */
/* Recipes in C, 2nd ed., 1992.  Function called ran2 in the book.  The   */
/* period is greater than 2E18.  Initialize by calling routine with a     */
/* large negative number.  Modified slightly [6/94] to make it compatible */
/* with C++.                                                              */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1 + IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

using namespace std;

//constructor
int random_dirty::random_num (long seed)

{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    //float temp;

    // setup shuffle array if not done so already
    if (seed <= 0) {
      _idum = seed;
      if (-_idum < 1) _idum = 1;
      else _idum = -_idum;
      idum2 = _idum;
      for (j = (NTAB+7); j >= 0; j--) {
	k = _idum/IQ1;
	_idum = IA1*(_idum - k*IQ1) - k*IR1;
	if (_idum < 0) _idum += IM1;
	if (j < NTAB) {
	  _iv[j] = _idum;
	}
      }
      iy = _iv[0];
    }
    return 0;
}

/* ====================================================================== */

float random_dirty::random_num ()

{

    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    float temp;

    k = _idum/IQ1;
    _idum = IA1*(_idum - k*IQ1) - k*IR1;
    if (_idum < 0) _idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j = int(iy/NDIV);
    iy = _iv[j] - idum2;
    _iv[j] = _idum;
    if (iy < 1) iy += IMM1;
    temp = AM*iy;
    if (temp > RNMX) {
//       cout << "RNMX = " << RNMX << endl;
      return RNMX;
    } else {
//       cout << "temp = " << temp << endl;
      return temp;
    }
}

/* ====================================================================== */
