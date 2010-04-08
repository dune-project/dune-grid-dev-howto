# $Id: Makefile.am 780 2007-07-30 19:52:27Z mblatt $

SUBDIRS = doxygen

if BUILD_DOCS
  DOCFILES = grid-dev-howto.pdf
  EXTRA_DIST = grid-dev-howto.pdf
endif


# setting like in dune-web
BASEDIR=../..
CURDIR=doc/grid-dev-howto

# where to install the docs
docdir=$(datadir)/doc/dune-grid-dev-howto

include $(top_srcdir)/am/latex
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
