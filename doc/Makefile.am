# $Id: Makefile.am 780 2007-07-30 19:52:27Z mblatt $

SUBDIRS = doxygen #layout

#TEXSOURCES = istl.tex istl.bib
#EPSFILES = blockstructure.eps
EXTRA_DIST = Makefile.dist.in Makefile.dist.am
#if BUILD_DOCS
#  DOCS = istl.pdf
#endif


# setting like in dune-web
BASEDIR=../..
CURDIR=doc/subgrid
docdir=$(datadir)/doc/dune-subgrid

doc_DATA = $(DOCS)

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.ps

# files to clean only in svn tree but not in tarball
SVNCLEANFILES =  *.pdf
clean-local:
	if test -e $(top_srcdir)/doc/doxygen/Doxydep; then rm -rf $(SVNCLEANFILES); fi

dist-hook:
	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

include $(top_srcdir)/am/latex
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
