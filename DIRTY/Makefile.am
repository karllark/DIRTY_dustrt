AM_CPPFLAGS=-Wall -Wextra -I${srcdir}/include ${cfitsio_CFLAGS}
AM_LDFLAGS=-lm ${cfitsio_LIBS}

if DEBUG
AM_CPPFLAGS+=-O0 -g -DDEBUG -DDEBUG_CDD -DDEBUG_CDD -DDEBUG_CDEC -DDEBUG_CONT_ABS -DDEBUG_COTM -DDEBUG_DIRTY -DDEBUG_GDEE -DDEBUG_GDP -DDEBUG_GDTE -DDEBUG_GSP -DDEBUG_GWG -DDEBUG_INITOUT -DDEBUG_MRT -DDEBUG_NP1S -DDEBUG_NPDD -DDEBUG_NPDS -DDEBUG_NPGS -DDEBUG_OGR -DDEBUG_OMG -DDEBUG_OUTR -DDEBUG_OUTRANGE -DDEBUG_RZFO -DDEBUG_SAEG -DDEBUG_SDG -DDEBUG_SDGCG -DDEBUG_SDGDD -DDEBUG_SDGF -DDEBUG_SDGSOC -DDEBUG_SEGFMC -DDEBUG_SP -DDEBUG_VDG -D_DIRTY_DEBUG_
endif

if DEBUG_FIXME
AM_CPPFLAGS+=-DDEBUG_CSCP -DDEBUG_CPT -DDEBUG_CSP -DDEBUG_DPPI -DDEBUG_DPPII -DDEBUG_FFS -DDEBUG_NS -DDEBUG_RT -DDEBUG_SWTO -DDEBUG_STWTO 
endif

bin_PROGRAMS=dirty
dirty_SOURCES= \
    calc_delta_dist.cpp \
    calc_photon_trajectory.cpp \
    check_absorbed_energy_grid.cpp \
    check_de_energy_conservation.cpp \
    check_fits_io.cpp \
    classify_scattered_photon.cpp \
    classify_stellar_photon.cpp \
    ComputeDustEmission.cpp \
    compute_observer_trans_matrix.cpp \
    ConfigFile.cpp \
    continuous_absorption.cpp \
    DataFile.cpp \
    determine_photon_position_index.cpp \
    determine_photon_position_index_initial.cpp \
    dirty.cpp \
    DirtyFailure.cpp \
    fits_params_to_header.cpp \
    forced_first_scatter.cpp \
    get_dust_ere_emission.cpp \
    get_dust_parameters.cpp \
    get_dust_scat_parameters.cpp \
    get_dust_thermal_emission.cpp \
    get_run_parameters.cpp \
    get_sed_parameters.cpp \
    get_wave_grid.cpp \
    Grain.cpp \
    GrainModel.cpp \
    initialize_output.cpp \
    new_photon.cpp \
    new_photon_dexp_disk.cpp \
    new_photon_diffuse_source.cpp \
    new_photon_discrete_stars.cpp \
    new_photon_grid_source.cpp \
    new_photon_pow_sphere.cpp \
    next_scatter.cpp \
    output_global_results.cpp \
    output_model_grid.cpp \
    output_results.cpp \
    radiative_transfer.cpp \
    radiative_transfer_many_waves.cpp \
    random_dirty.cpp \
    rotate_zaxis_for_observer.cpp \
    scattered_weight_towards_observer.cpp \
    scatter_photon.cpp \
    setup_absorbed_energy_grid.cpp \
    setup_dust_grid_check_grid.cpp \
    setup_dust_grid.cpp \
    setup_dust_grid_dexp_disk.cpp \
    setup_dust_grid_file.cpp \
    setup_dust_grid_shell.cpp \
    setup_dust_grid_slab.cpp \
    setup_dust_grid_sphere.cpp \
    setup_dust_grid_subdivide_overdense_cells.cpp \
    setup_emitted_grid_for_montecarlo.cpp \
    setup_ere_dust_emission.cpp \
    setup_thermal_dust_emission_output.cpp \
    stellar_weight_towards_observer.cpp \
    StochasticHeating.cpp \
    store_absorbed_energy_grid.cpp \
    StringManip.cpp \
    verify_dust_grid.cpp

EXTRA_DIST = \
    $(dirty_SOURCES) \
    include
