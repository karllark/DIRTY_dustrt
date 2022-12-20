#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([DIRTY_dustrt],
    m4_esyscmd([.build-aux/git-version-gen .tarball-version]),
    [https://github.com/karllark/DIRTY_dustrt/issues],
    [],
    https://github.com/karllark/DIRTY_dustrt)
AC_CONFIG_SRCDIR([DIRTY/dirty.cpp])
AC_CONFIG_HEADERS([DIRTY/include/config.h])
AC_LANG(C++)

# Not a GNU project
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
PKG_PROG_PKG_CONFIG

# Checks for libraries.
AC_ARG_WITH([cfitsio],
    AS_HELP_STRING([--with-cfitsio], [Path to prefix containing cfitsio]))

cfitsio_usage="unable to find cfitsio library.
  USE --with-cfitsio=/path/to/cfitsio
  OR set: PKG_CONFIG_PATH
  OR set: cfitsio_CFLAGS and cfitsio_LIBS"

AS_IF([test "x$with_cfitsio" != "x"],
    [test "$with_cfitsio" == "yes" && AC_SUBST(with_cfitsio, /usr/local)
        cfitsio_CFLAGS="-I${with_cfitsio}/include -I${with_cfitsio}/include/fitsio"
        cfitsio_LIBS="-L${with_cfitsio}/lib -Wl,-rpath=${with_cfitsio}/lib -lcfitsio"],
    [PKG_CHECK_MODULES([cfitsio], [cfitsio >= 3], [],
        [AC_MSG_ERROR($cfitsio_usage)])])
AC_MSG_NOTICE([cfitsio CFLAGS... $cfitsio_CFLAGS])
AC_MSG_NOTICE([cfitsio LIBS... $cfitsio_LIBS])

# Push global flags
AC_SUBST(CPPFLAGS_OLD, $CPPFLAGS)
AC_SUBST(LDFLAGS_OLD, $LDFLAGS)

# Push cfitsio flags
AC_SUBST(CPPFLAGS,$cfitsio_CFLAGS)
AC_SUBST(LDFLAGS,$cfitsio_LIBS)

# Verify cfitsio works
AC_CHECK_HEADERS([fitsio.h], [], [AC_MSG_ERROR(cfitsio headers not found! $cfitsio_usage)])
AC_CHECK_LIB(cfitsio, [ffopen], [], [AC_MSG_ERROR(cfitsio library not found! $cfitsio_usage)])
AC_CHECK_FUNC(ffopen, [], [AC_MSG_ERROR(cfitsio is not usable! $cfitsio_usage)])

# Pop global flags
AC_SUBST(CPPFLAGS, $CPPFLAGS_OLD)
AC_SUBST(LDFLAGS, $LDFLAGS_OLD)

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([atoi atof strtod strtof])
AC_CHECK_FUNCS([fabs sqrt pow finite isfinite])

# Check documentation engine
AC_CACHE_CHECK([for Sphinx], [ac_cv_path_SPHINX],
    [AC_PATH_PROGS_FEATURE_CHECK([SPHINX], [sphinx-build],
    	[AS_IF([$ac_path_SPHINX --version 2>&1 >/dev/null],
    		[ac_cv_path_SPHINX=$ac_path_SPHINX ac_path_SPHINX_found=:])],
    	[AC_MSG_WARN([Sphinx not found (man page will not be generated)])])])
AC_SUBST([SPHINX], [$ac_cv_path_SPHINX])
AM_CONDITIONAL([HAVE_SPHINX], [test -x "$SPHINX"])

# Debugging facilities
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [Turn on debugging])],
    [debug=true],
    [debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

AC_ARG_ENABLE([debug-fixme],
    [AS_HELP_STRING([--enable-debug-fixme], [Turn on debugging (manual intervention required)])],
    [debug_fixme=true; debug=true],
    [debug_fixme=false; debug=false])
AM_CONDITIONAL([DEBUG_FIXME], [test "x$debug_fixme" = "xtrue"])

AC_CONFIG_FILES([Makefile
                 DIRTY/Makefile])
AC_OUTPUT
