# Makefile for DIRTY (version 2)
#   Dust Radiative Transfer and Emission Code
#   Karl D. Gordon and Karl A. Misselt
#
# started KDG (Apr 2003)

DIRTY_SRCS = dirty.cpp \
	ConfigFile.cpp \
	DataFile.cpp \
	StringManip.cpp \
	GrainModel.cpp \
	Grain.cpp \
	get_run_parameters.cpp \
	get_dust_parameters.cpp \
	get_sed_parameters.cpp \
	get_wave_grid.cpp \
	setup_dust_grid.cpp \
	setup_dust_grid_sphere.cpp \
	setup_dust_grid_shell.cpp \
	setup_dust_grid_slab.cpp \
	get_dust_scat_parameters.cpp \
	radiative_transfer.cpp \
	initialize_output.cpp \
	new_photon.cpp \
	new_photon_discrete_stars.cpp \
	new_photon_diffuse_source.cpp \
	new_photon_grid_source.cpp \
	forced_first_scatter.cpp \
	stellar_weight_towards_observer.cpp \
	scattered_weight_towards_observer.cpp \
	rotate_zaxis_for_observer.cpp \
	classify_stellar_photon.cpp \
	scatter_photon.cpp \
	classify_scattered_photon.cpp \
	next_scatter.cpp \
	output_results.cpp \
	output_global_results.cpp \
	setup_absorbed_energy_grid.cpp \
	store_absorbed_energy_grid.cpp \
	get_dust_thermal_emission.cpp \
	setup_emitted_grid_for_montecarlo.cpp \
	determine_photon_position_index_initial.cpp \
	determine_photon_position_index.cpp \
	calc_delta_dist.cpp \
	calc_photon_trajectory.cpp \
	check_fits_io.cpp \
	random_dirty.cpp

DIRTY_OBJS = ${DIRTY_SRCS:.cpp=.o}

RM = rm -f

CC = g++
MAKE = make
DIFLAGS = -Y./include

#for linux
LDFLAGS = -lcfitsio
CCFLAGS = -O2 -Iinclude -Wall -Wextra #-I/home/kgordon/Bin/DHAS/MPipeline/DAT/lib/cfitsio -L/home/kgordon/Bin/DHAS/MPipeline/DAT/lib/cfitsio #-ansi -pedantic #-pg # -fpermissive
#LDFLAGS = -L../lib/CFitsIO -lcfitsio
#CCFLAGS = -Iinclude -I../src_gen/include -I../lib/CFitsIO -Wno-deprecated -O3
#for solaris
#LDFLAGS = -L../lib/CFitsIO -lcfitsio -lsocket -lnsl
#CCFLAGS = -Iinclude -I../src_gen/include -I../lib/CFitsIO -Wno-deprecated -O3 -DSOLARIS

all: dirty

dirty: $(DIRTY_OBJS)
	$(CC) -o ${@} $(DIRTY_OBJS) $(CCFLAGS) $(LDFLAGS)

depend: 
	makedepend -- $(DIFLAGS) -- $(DIRTY_SRCS)

# ----------
# cleaning targets
# ----------

clean_all: clean

clean: 
	$(RM) *.o
	$(RM) dirty

# ----------
# suffixes
# ----------

.SUFFIXES : .o .cpp

.cpp.o : $(HDS)
	 $(CC) -o ${*}.o $(CCFLAGS) -c $<

# DO NOT DELETE

dirty.o: ./include/dirty.h ./include/debug.h ./include/geometry_def.h
dirty.o: ./include/NumUtils.h ./include/Constants.h ./include/grid_cell.h
dirty.o: ./include/output_def.h ./include/runinfo_def.h
dirty.o: ./include/photon_data.h ./include/random_dirty.h
dirty.o: ./include/ConfigFile.h ./include/GrainModel.h ./include/Grain.h
dirty.o: ./include/StringManip.h
ConfigFile.o: ./include/ConfigFile.h
DataFile.o: ./include/DataFile.h
StringManip.o: ./include/StringManip.h
GrainModel.o: ./include/GrainModel.h ./include/ConfigFile.h ./include/Grain.h
GrainModel.o: ./include/Constants.h ./include/NumUtils.h
GrainModel.o: ./include/StringManip.h
Grain.o: ./include/Grain.h ./include/Constants.h ./include/NumUtils.h
Grain.o: ./include/StringManip.h
get_run_parameters.o: ./include/get_run_parameters.h ./include/ConfigFile.h
get_run_parameters.o: ./include/output_def.h ./include/NumUtils.h
get_run_parameters.o: ./include/Constants.h ./include/geometry_def.h
get_run_parameters.o: ./include/grid_cell.h ./include/runinfo_def.h
get_run_parameters.o: ./include/check_input_param.h ./include/GrainModel.h
get_run_parameters.o: ./include/Grain.h ./include/StringManip.h
get_dust_parameters.o: ./include/get_dust_parameters.h ./include/ConfigFile.h
get_dust_parameters.o: ./include/DataFile.h ./include/runinfo_def.h
get_dust_parameters.o: ./include/check_input_param.h ./include/GrainModel.h
get_dust_parameters.o: ./include/Grain.h ./include/Constants.h
get_dust_parameters.o: ./include/NumUtils.h ./include/StringManip.h
get_sed_parameters.o: ./include/get_sed_parameters.h ./include/runinfo_def.h
get_sed_parameters.o: ./include/ConfigFile.h ./include/DataFile.h
get_sed_parameters.o: ./include/NumUtils.h ./include/Constants.h
get_sed_parameters.o: ./include/check_input_param.h
get_wave_grid.o: ./include/get_wave_grid.h ./include/ConfigFile.h
get_wave_grid.o: ./include/Constants.h ./include/runinfo_def.h
get_wave_grid.o: ./include/check_input_param.h
setup_dust_grid.o: ./include/setup_dust_grid.h ./include/ConfigFile.h
setup_dust_grid.o: ./include/DataFile.h ./include/check_input_param.h
setup_dust_grid.o: ./include/geometry_def.h ./include/NumUtils.h
setup_dust_grid.o: ./include/Constants.h ./include/grid_cell.h
setup_dust_grid.o: ./include/photon_data.h ./include/random_dirty.h
setup_dust_grid.o: ./include/constants.h ./include/debug.h
setup_dust_grid_sphere.o: ./include/setup_dust_grid_sphere.h
setup_dust_grid_sphere.o: ./include/ConfigFile.h
setup_dust_grid_sphere.o: ./include/check_input_param.h
setup_dust_grid_sphere.o: ./include/geometry_def.h ./include/NumUtils.h
setup_dust_grid_sphere.o: ./include/Constants.h ./include/grid_cell.h
setup_dust_grid_sphere.o: ./include/photon_data.h ./include/random_dirty.h
setup_dust_grid_sphere.o: ./include/debug.h
setup_dust_grid_shell.o: ./include/setup_dust_grid_shell.h
setup_dust_grid_shell.o: ./include/ConfigFile.h ./include/check_input_param.h
setup_dust_grid_shell.o: ./include/geometry_def.h ./include/NumUtils.h
setup_dust_grid_shell.o: ./include/Constants.h ./include/grid_cell.h
setup_dust_grid_shell.o: ./include/photon_data.h ./include/random_dirty.h
setup_dust_grid_shell.o: ./include/debug.h
setup_dust_grid_slab.o: ./include/setup_dust_grid_slab.h
setup_dust_grid_slab.o: ./include/ConfigFile.h ./include/check_input_param.h
setup_dust_grid_slab.o: ./include/geometry_def.h ./include/NumUtils.h
setup_dust_grid_slab.o: ./include/Constants.h ./include/grid_cell.h
setup_dust_grid_slab.o: ./include/photon_data.h ./include/random_dirty.h
setup_dust_grid_slab.o: ./include/debug.h
get_dust_scat_parameters.o: ./include/get_dust_scat_parameters.h
get_dust_scat_parameters.o: ./include/geometry_def.h ./include/NumUtils.h
get_dust_scat_parameters.o: ./include/Constants.h ./include/grid_cell.h
get_dust_scat_parameters.o: ./include/runinfo_def.h
radiative_transfer.o: ./include/radiative_transfer.h ./include/geometry_def.h
radiative_transfer.o: ./include/NumUtils.h ./include/Constants.h
radiative_transfer.o: ./include/grid_cell.h ./include/runinfo_def.h
radiative_transfer.o: ./include/output_def.h ./include/photon_data.h
radiative_transfer.o: ./include/random_dirty.h ./include/debug.h
initialize_output.o: ./include/initialize_output.h ./include/output_def.h
initialize_output.o: ./include/NumUtils.h ./include/Constants.h
initialize_output.o: ./include/geometry_def.h ./include/grid_cell.h
initialize_output.o: ./include/debug.h
new_photon.o: ./include/new_photon.h ./include/geometry_def.h
new_photon.o: ./include/NumUtils.h ./include/Constants.h
new_photon.o: ./include/grid_cell.h ./include/random_dirty.h
new_photon.o: ./include/photon_data.h ./include/debug.h
new_photon_discrete_stars.o: ./include/new_photon_discrete_stars.h
new_photon_discrete_stars.o: ./include/geometry_def.h ./include/NumUtils.h
new_photon_discrete_stars.o: ./include/Constants.h ./include/grid_cell.h
new_photon_discrete_stars.o: ./include/random_dirty.h ./include/photon_data.h
new_photon_discrete_stars.o: ./include/debug.h
new_photon_diffuse_source.o: ./include/new_photon_diffuse_source.h
new_photon_diffuse_source.o: ./include/geometry_def.h ./include/NumUtils.h
new_photon_diffuse_source.o: ./include/Constants.h ./include/grid_cell.h
new_photon_diffuse_source.o: ./include/random_dirty.h ./include/photon_data.h
new_photon_diffuse_source.o: ./include/debug.h
forced_first_scatter.o: ./include/forced_first_scatter.h
forced_first_scatter.o: ./include/geometry_def.h ./include/NumUtils.h
forced_first_scatter.o: ./include/Constants.h ./include/grid_cell.h
forced_first_scatter.o: ./include/photon_data.h ./include/random_dirty.h
forced_first_scatter.o: ./include/roundoff_err.h ./include/debug.h
stellar_weight_towards_observer.o: ./include/stellar_weight_towards_observer.h
stellar_weight_towards_observer.o: ./include/geometry_def.h
stellar_weight_towards_observer.o: ./include/NumUtils.h ./include/Constants.h
stellar_weight_towards_observer.o: ./include/grid_cell.h
stellar_weight_towards_observer.o: ./include/photon_data.h ./include/debug.h
scattered_weight_towards_observer.o: ./include/scattered_weight_towards_observer.h
scattered_weight_towards_observer.o: ./include/geometry_def.h
scattered_weight_towards_observer.o: ./include/NumUtils.h
scattered_weight_towards_observer.o: ./include/Constants.h
scattered_weight_towards_observer.o: ./include/grid_cell.h
scattered_weight_towards_observer.o: ./include/photon_data.h
scattered_weight_towards_observer.o: ./include/debug.h
rotate_zaxis_for_observer.o: ./include/rotate_zaxis_for_observer.h
rotate_zaxis_for_observer.o: ./include/photon_data.h ./include/debug.h
classify_stellar_photon.o: ./include/classify_stellar_photon.h
classify_stellar_photon.o: ./include/output_def.h ./include/NumUtils.h
classify_stellar_photon.o: ./include/Constants.h ./include/geometry_def.h
classify_stellar_photon.o: ./include/grid_cell.h ./include/photon_data.h
classify_stellar_photon.o: ./include/debug.h
scatter_photon.o: ./include/scatter_photon.h ./include/geometry_def.h
scatter_photon.o: ./include/NumUtils.h ./include/Constants.h
scatter_photon.o: ./include/grid_cell.h ./include/photon_data.h
scatter_photon.o: ./include/random_dirty.h ./include/roundoff_err.h
scatter_photon.o: ./include/debug.h
classify_scattered_photon.o: ./include/classify_scattered_photon.h
classify_scattered_photon.o: ./include/output_def.h ./include/NumUtils.h
classify_scattered_photon.o: ./include/Constants.h ./include/geometry_def.h
classify_scattered_photon.o: ./include/grid_cell.h ./include/photon_data.h
classify_scattered_photon.o: ./include/debug.h
next_scatter.o: ./include/next_scatter.h ./include/geometry_def.h
next_scatter.o: ./include/NumUtils.h ./include/Constants.h
next_scatter.o: ./include/grid_cell.h ./include/photon_data.h
next_scatter.o: ./include/random_dirty.h ./include/roundoff_err.h
next_scatter.o: ./include/debug.h
output_results.o: ./include/output_results.h ./include/geometry_def.h
output_results.o: ./include/NumUtils.h ./include/Constants.h
output_results.o: ./include/grid_cell.h ./include/output_def.h
output_results.o: ./include/runinfo_def.h ./include/constants.h
output_results.o: ./include/debug.h
setup_absorbed_energy_grid.o: ./include/setup_absorbed_energy_grid.h
setup_absorbed_energy_grid.o: ./include/geometry_def.h ./include/NumUtils.h
setup_absorbed_energy_grid.o: ./include/Constants.h ./include/grid_cell.h
setup_absorbed_energy_grid.o: ./include/debug.h
store_absorbed_energy_grid.o: ./include/store_absorbed_energy_grid.h
store_absorbed_energy_grid.o: ./include/Constants.h ./include/geometry_def.h
store_absorbed_energy_grid.o: ./include/NumUtils.h ./include/grid_cell.h
store_absorbed_energy_grid.o: ./include/runinfo_def.h ./include/output_def.h
store_absorbed_energy_grid.o: ./include/debug.h
get_dust_thermal_emission.o: ./include/get_dust_thermal_emission.h
get_dust_thermal_emission.o: ./include/geometry_def.h ./include/NumUtils.h
get_dust_thermal_emission.o: ./include/Constants.h ./include/grid_cell.h
get_dust_thermal_emission.o: ./include/runinfo_def.h ./include/GrainModel.h
get_dust_thermal_emission.o: ./include/ConfigFile.h ./include/Grain.h
get_dust_thermal_emission.o: ./include/StringManip.h
determine_photon_position_index_initial.o: ./include/determine_photon_position_index_initial.h
determine_photon_position_index_initial.o: ./include/geometry_def.h
determine_photon_position_index_initial.o: ./include/NumUtils.h
determine_photon_position_index_initial.o: ./include/Constants.h
determine_photon_position_index_initial.o: ./include/grid_cell.h
determine_photon_position_index_initial.o: ./include/photon_data.h
determine_photon_position_index.o: ./include/determine_photon_position_index.h
determine_photon_position_index.o: ./include/geometry_def.h
determine_photon_position_index.o: ./include/NumUtils.h ./include/Constants.h
determine_photon_position_index.o: ./include/grid_cell.h
determine_photon_position_index.o: ./include/photon_data.h ./include/debug.h
calc_delta_dist.o: ./include/calc_delta_dist.h ./include/geometry_def.h
calc_delta_dist.o: ./include/NumUtils.h ./include/Constants.h
calc_delta_dist.o: ./include/grid_cell.h ./include/photon_data.h
calc_delta_dist.o: ./include/roundoff_err.h ./include/debug.h
calc_photon_trajectory.o: ./include/calc_photon_trajectory.h
calc_photon_trajectory.o: ./include/debug.h ./include/geometry_def.h
calc_photon_trajectory.o: ./include/NumUtils.h ./include/Constants.h
calc_photon_trajectory.o: ./include/grid_cell.h ./include/photon_data.h
check_fits_io.o: ./include/check_fits_io.h
random_dirty.o: ./include/random_dirty.h
