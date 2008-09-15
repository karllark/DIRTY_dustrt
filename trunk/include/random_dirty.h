#define NTAB 32 

#include <iostream>

class random_dirty {
 public:
  // Constructors/destructors. 
  // need some
  int random_num(long seed);
  
  // get new random value
  float random_num();

 private:
  long _idum;
  long _iv[NTAB];
};

