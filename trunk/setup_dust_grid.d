setup_dust_grid.o setup_dust_grid.d : setup_dust_grid.cpp include/setup_dust_grid.h \
  include/ConfigFile.h include/DataFile.h include/check_input_param.h \
  include/geometry_def.h include/NumUtils.h include/Constants.h \
  include/grid_cell.h include/photon_data.h include/random_dirty.h \
  include/constants.h include/debug.h
