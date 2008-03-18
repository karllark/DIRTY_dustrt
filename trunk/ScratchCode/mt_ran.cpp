#include "mt_ran.h"

namespace mtran { 

  mtrangen::mtrangen( unsigned long seed ) { mt_seed(seed); }

  void mtrangen::mt_seed ( unsigned long seed )
  {
    
    x[0] = seed & 0xffffffffUL;
    for (int i = 1; i < N; i++ )
      x[i] = 0xffffffffUL & ( 1812433253UL * ( x[i - 1] ^ ( x[i - 1] >> 30 ) ) + i );
    
  }
  
  unsigned long mtrangen::mt_rand ( void )
  {
    unsigned long y, a;
    int i;
//	std::cout << "next is: " << next << std::endl;     
    /* Refill x if exhausted */
    if ( next == N ) {
      next = 0;
      for ( i = 0; i < N - 1; i++ ) {
	y = ( x[i] & U ) | x[i + 1] & L;
	a = ( y & 0x1UL ) ? A : 0x0UL;
	x[i] = x[( i + M ) % N] ^ ( y >> 1 ) ^ a;
      }
      
      y = ( x[N - 1] & U ) | x[0] & L;
      a = ( y & 0x1UL ) ? A : 0x0UL;
      x[N - 1] = x[M - 1] ^ ( y >> 1 ) ^ a;
    }
    
    y = x[next++];
    
    /* Improve distribution */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
  
    return y;
  }
  
//   std::vector <float> mtrangen::uDev_series ( int n ) { 

//     std::vector <float> retvect; 
//     for (int i=0;i<n;i++) retvect.push_back( (float)(mt_rand())*Scale01); 
//     return retvect; 

//   }

} // End Namespace definition
