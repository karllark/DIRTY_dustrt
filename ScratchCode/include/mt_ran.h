#ifndef _MT_RAN_H_
#define _MT_RAN_H_

// Implement Mersenne Twister random number generator.

#include <vector> 

namespace mtran { 

  class mtrangen { 
    
  public : 

    // Default constructor with default seed
    mtrangen ( unsigned long seed = 19650218UL);
    // set the seed.
    void mt_seed ( unsigned long seed); 
    // Generate random number
    unsigned long mt_rand ( void );
    // Generate random deviates, single or series. Done INLINE
    inline float uDev( void ) { return (float)(mt_rand())*(Scale01); }
    inline std::vector <float> uDev_series ( int n ) { 
      std::vector <float> retvect; 
      for (int i=0;i<n;i++) retvect.push_back( (float)(mt_rand())*Scale01); 
      return retvect; 
    }
       
  private : 
  
    static const unsigned long N = 624; 
    static const unsigned long M = 397; 
    
    static const unsigned long A = 0x9908b0dfUL;
    static const unsigned long U = 0x80000000UL;
    static const unsigned long L = 0x7fffffffUL; 

    long next; 
    unsigned long x[N]; 

    static const double Scale01=1.0/4294967295.0;

  }; // End Class definition.

}  // End Namespace definition.

#endif
