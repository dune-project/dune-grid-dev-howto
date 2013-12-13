AC_DEFUN([DUNE_GRID_DEV_HOWTO_CHECKS],[
  AC_MSG_WARN([the module dune-grid-dev-howto is deprecated and will be removed after Dune 2.3])
])

AC_DEFUN([DUNE_GRID_DEV_HOWTO_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-grid-dev-howto...])
  DUNE_CHECK_MODULES([dune-grid-dev-howto],[grid/identitygrid.hh])
])
