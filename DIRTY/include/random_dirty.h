#define NTAB 32

#include <cassert>
#include <iostream>

class random_dirty
{
  public:
    // Constructors/destructors.
    // need some
    random_dirty (long seed);

    // get new random value
    double random_num ();

  private:
    long _idum;
    long _iv[NTAB];

    long idum2;
    long iy;
};
