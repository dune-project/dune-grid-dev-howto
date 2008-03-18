
TESTPROGS = test-w-sgrid

# which tests to run
TESTS = $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)

AM_LDFLAGS = @AM_LDFLAGS@
AM_CXXFLAGS = -Wall -g

# define the programs
test_w_sgrid_SOURCES = test-w-sgrid.cc

include $(top_srcdir)/am/global-rules
