#define NTAB 32 

#include <iostream>
#include <cassert>

class random_dirty {
 public:
  // Constructors/destructors. 
  // need some
  random_dirty(long seed);
  
  // get new random value
  float random_num();

 private:
  long _idum;
  long _iv[NTAB];

  long idum2;
  long iy;
};

