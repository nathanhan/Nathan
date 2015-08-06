SUBDIRS = SnowTools/bwa SnowTools/htslib SnowTools/src 

bin_PROGRAMS = nathan

nathan_SOURCES = nathan.cpp

nathan_CPPFLAGS = \
	-I$(top_srcdir)/SnowTools \
	-I$(top_srcdir)/SnowTools/src \
	-I$(top_srcdir)/SnowTools/htslib \
	-I$(top_srcdir)/SnowTools/multifast-v1.4.2

nathan_LDADD = \
	$(top_srcdir)/SnowTools/src/libsnowtools.a \
	$(top_srcdir)/SnowTools/bwa/libbwa.a \
	$(top_srcdir)/SnowTools/htslib/libhts.a
