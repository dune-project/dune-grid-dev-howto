#! /bin/bash
# Process this file with autoconf to produce a configure script.
DUNE_AC_INIT
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/grid/identitygrid.hh])
AM_CONFIG_HEADER([config.h])

# check all dune-module stuff
DUNE_CHECK_ALL

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS) $(DUNE_LIBS)')

DUNE_SUMMARY_ALL

# write output
AC_CONFIG_FILES([Makefile
    dune/Makefile
    dune/grid/Makefile
    dune/grid/identitygrid/Makefile
    dune/grid/test/Makefile
    doc/doxygen/Makefile
    doc/doxygen/Doxyfile
    doc/Makefile
    m4/Makefile
    dune-grid-dev-howto.pc])
AC_OUTPUT
