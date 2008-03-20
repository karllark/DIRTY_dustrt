scatter_photon.o scatter_photon.d : scatter_photon.cpp include/scatter_photon.h \
  include/geometry_def.h include/NumUtils.h include/Constants.h \
  include/grid_cell.h include/photon_data.h include/random_dirty.h \
  include/roundoff_err.h include/debug.h
