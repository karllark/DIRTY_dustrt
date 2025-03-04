// ======================================================================
//   Procedure to check the status of the CFITSIO routines.
//
// 2006 Apr/KDG - written
// ======================================================================
#include "check_fits_io.h"

int check_fits_io(int status, const char text[100])

{
  char errtxt[FLEN_ERRMSG];

  if (status != 0) {
    fits_get_errstatus(status, errtxt);
    cout << text;
    cout << " status = " << status;
    cout << " : " << errtxt << endl;
    cout.flush();
    exit(8);
  }

  return 0;
}
