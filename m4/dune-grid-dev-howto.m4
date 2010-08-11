AC_DEFUN([DUNE_GRID_DEV_HOWTO_CHECKS],[])

AC_DEFUN([DUNE_GRID_DEV_HOWTO_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-grid-dev-howto...])
  DUNE_CHECK_MODULES([dune-grid-dev-howto],[grid/identitygrid.hh])
])
