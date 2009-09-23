get_run_parameters.o get_run_parameters.d : get_run_parameters.cpp include/get_run_parameters.h \
  include/ConfigFile.h include/output_def.h include/NumUtils.h \
  include/Constants.h include/geometry_def.h include/grid_cell.h \
  include/runinfo_def.h include/check_input_param.h include/GrainModel.h \
  include/Grain.h include/StringManip.h include/EqTemp.h \
  include/DirtyFlags.h
