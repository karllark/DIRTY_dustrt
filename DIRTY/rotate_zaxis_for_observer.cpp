// ======================================================================
//   Procedure to classify stellar photon(s) into the output images
// and global totals.
//
// 2006 Apr/KDG - written
// ======================================================================
#include "rotate_zaxis_for_observer.h"

void
rotate_zaxis_for_observer (float transform[3][3], photon_data &photon)

{
  int i;
  int j;
  double new_position[3]; // place to store new position as it is being computed

  for (i = 0; i < 3; i++)
    {
      new_position[i] = 0.0;
      for (j = 0; j < 3; j++)
        new_position[i] += transform[i][j] * photon.position[j];
    }

#ifdef DEBUG_RZFO
  std::cout << "photon # = " << photon.number << std::endl;
  std::cout << "old photon positions = ";
  double radius = 0.0;
  for (i = 0; i < 3; i++)
    {
      std::cout << photon.position[i] << " ";
      radius += pow (photon.position[i], 2.0);
    }
  std::cout << "; radius = " << sqrt (radius) << std::endl;
  radius = 0.0;
  std::cout << "new photon positions = ";
  for (i = 0; i < 3; i++)
    {
      std::cout << new_position[i] << " ";
      radius += pow (photon.position[i], 2.0);
    }
  std::cout << "; radius = " << sqrt (radius) << std::endl;

#endif

  for (i = 0; i < 3; i++)
    photon.position[i] = new_position[i];
}
