#ifndef DIRTY_COMPAT_H

#include "config.h"

#if !defined(HAVE_FINITE)
#define finite isfinite
#endif

#endif  // DIRTY_COMPAT_H
