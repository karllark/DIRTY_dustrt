SUBDIRS = \
    DIRTY

dirtydir=$(pkgdatadir)
dirty_exampledir=$(dirtydir)/example
dirty_example_trust_slabdir=$(dirtydir)/example/trust_slab
dirty_example_dmodel_testdir=$(dirtydir)/example/dmodel_test
dirty_licensesdir=$(dirtydir)/licenses
dirty_prodir=$(dirtydir)/pro
dirty_seddir=$(dirtydir)/SEDs
dirty_docsdir=$(mandir)

dirty_DATA = \
    README.rst

dirty_example_DATA = \
    example/cubical_clumps.param \
    example/multiple_stars_sphere.param \
    example/spherical_clumps.param \
    example/standard_shell_outputfailureadded.param \
    example/standard_shell.param \
    example/standard_slab.param \
    example/standard_sphere_diffuse.param \
    example/standard_sphere_multiwave.param \
    example/standard_sphere.param \
    example/test_sphere_tau3.param

dirty_example_dmodel_test_DATA = \
    example/dmodel_test/standard_sphere_mwave_dmodel.param

dirty_example_trust_slab_DATA = \
    example/trust_slab/basic_grid.dat \
    example/trust_slab/BB_T10000_L100000.dat \
    example/trust_slab/DIRTY_bm1_effgrain_tau_001.00_theta_090.drt \
    example/trust_slab/DIRTY_bm1_eqgrain_tau_001.00_theta_090.drt \
    example/trust_slab/DIRTY_bm1_neqgrain_tau_001.00_theta_090.drt \
    example/trust_slab/full_grid.dat \
    example/trust_slab/output/file_for_output_dir

dirty_sed_DATA = \
    SEDs/bluesg_T20000_R10.dat \
    SEDs/F98_p00_c_10.dat \
    SEDs/redsg_kycyg_T3500_R2100.dat

dirty_licenses_DATA = \
    licenses/LICENSE.rst

dirty_pro_DATA = \
    pro/blackbody_dirtyv2.pro \
    pro/close_ps_output.pro \
    pro/comp_dirty_global_seds.pro \
    pro/create_wave_grid.pro \
    pro/find_uniq_rad_fields.pro \
    pro/generate_blackbody_input.pro \
    pro/get_new_krange.pro \
    pro/klegend.pro \
    pro/koplot.pro \
    pro/kplot.pro \
    pro/krange.pro \
    pro/plot_dirty_atten.pro \
    pro/plot_dirty_global_sed.pro \
    pro/plot_dirty_image_cut.pro \
    pro/plot_model_grid_tau.pro \
    pro/plot_rad_field_info.pro \
    pro/set_symbol.pro \
    pro/setup_colors.pro \
    pro/setup_ps_output.pro

EXTRA_DIST = \
    $(dirty_DATA) \
    $(dirty_example_DATA) \
    $(dirty_example_dmodel_test_DATA) \
    $(dirty_example_trust_slab_DATA) \
    $(dirty_licenses_DATA) \
    $(dirty_pro_DATA) \
    $(dirty_sed_DATA) \
    $(top_srcdir)/docs \
    $(top_srcdir)/.build-aux/git-version-gen \
    $(top_srcdir)/bootstrap \
    $(top_srcdir)/.version

#
# Handle GIT versioning
#

BUILT_SOURCES = $(top_srcdir)/.version
$(top_srcdir)/.version:
	echo '$(VERSION)' > $@-t
	mv $@-t $@
dist-hook:
	echo '$(VERSION)' > $(distdir)/.tarball-version

if HAVE_SPHINX
all-local: man

man1_MANS = docs/_build/man/dirty.1
docs/_build/man/dirty.1:
	$(MAKE) -C $(top_srcdir)/docs man BUILDDIR=$(abs_builddir)/docs/_build

man: docs/_build/man/dirty.1

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf docs/_build
endif
