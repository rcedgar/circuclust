UNAME_S := $(shell uname -s)

CPPFLAGS := $(CPPFLAGS) -DNDEBUG -pthread

CXX = g++
CXXFLAGS := $(CXXFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

LDFLAGS := $(LDFLAGS) -O3 -fopenmp -pthread -lpthread
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -static
endif

HDRS = \
  alndata.h \
  alpha.h \
  cmds.h \
  countsort.h \
  fastaseqsource.h \
  fastq.h \
  fastqseqsource.h \
  fileseqsource.h \
  filetype.h \
  gobuff.h \
  linereader.h \
  lockobj.h \
  lockobjs.h \
  mx.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  omplock.h \
  partition.h \
  pathinfo.h \
  seq.h \
  seqdb.h \
  seqinfo.h \
  seqsource.h \
  sfasta.h \
  sort.h \
  tracebit.h \
  xdpmem.h \

OBJS = \
  o/align.o \
  o/alignpair.o \
  o/alndata.o \
  o/alpha.o \
  o/cluster.o \
  o/codes.o \
  o/fastaseqsource.o \
  o/fastq.o \
  o/fastqseqsource.o \
  o/fileseqsource.o \
  o/filetype.o \
  o/linereader.o \
  o/lockobj.o \
  o/mx.o \
  o/objmgr.o \
  o/circuclust_main.o \
  o/myutils.o \
  o/overhangs.o \
  o/partition.o \
  o/pathinfo.o \
  o/revcompseq.o \
  o/rotate.o \
  o/rotatelongorf.o \
  o/search.o \
  o/seqdb.o \
  o/seqdbfromfasta.o \
  o/seqinfo.o \
  o/seqsource.o \
  o/sfasta.o \
  o/substmx.o \
  o/tracebackbitmem.o \
  o/translator.o \
  o/usort.o \
  o/viterbifastmem.o \
  o/writealn.o \

.PHONY: clean

o/circuclust : o/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@
	strip -d o/circuclust

o/ :
	mkdir -p o/

o/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -rf o/
