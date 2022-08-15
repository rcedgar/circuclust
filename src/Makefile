BINDIR := ../bin
OBJDIR := o
BINPATH := $(BINDIR)/circuclust

CPPFLAGS := $(CPPFLAGS) -DNDEBUG -pthread

CXX = g++
CXXFLAGS := $(CXXFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

UNAME_S := $(shell uname -s)
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
  translator.h \
  xdpmem.h \

OBJS = \
  $(OBJDIR)/align.o \
  $(OBJDIR)/alignpair.o \
  $(OBJDIR)/alndata.o \
  $(OBJDIR)/alpha.o \
  $(OBJDIR)/cluster.o \
  $(OBJDIR)/codes.o \
  $(OBJDIR)/fastaseqsource.o \
  $(OBJDIR)/fastq.o \
  $(OBJDIR)/fastqseqsource.o \
  $(OBJDIR)/fileseqsource.o \
  $(OBJDIR)/filetype.o \
  $(OBJDIR)/linereader.o \
  $(OBJDIR)/lockobj.o \
  $(OBJDIR)/mx.o \
  $(OBJDIR)/objmgr.o \
  $(OBJDIR)/circuclust_main.o \
  $(OBJDIR)/myutils.o \
  $(OBJDIR)/overhangs.o \
  $(OBJDIR)/partition.o \
  $(OBJDIR)/pathinfo.o \
  $(OBJDIR)/revcompseq.o \
  $(OBJDIR)/rotate.o \
  $(OBJDIR)/rotatelongorf.o \
  $(OBJDIR)/search.o \
  $(OBJDIR)/seqdb.o \
  $(OBJDIR)/seqdbfromfasta.o \
  $(OBJDIR)/seqinfo.o \
  $(OBJDIR)/seqsource.o \
  $(OBJDIR)/sfasta.o \
  $(OBJDIR)/substmx.o \
  $(OBJDIR)/tracebackbitmem.o \
  $(OBJDIR)/translator.o \
  $(OBJDIR)/usort.o \
  $(OBJDIR)/viterbifastmem.o \
  $(OBJDIR)/writealn.o \

.PHONY: clean

$(BINPATH) : $(BINDIR)/ $(OBJDIR)/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(BINPATH)
	strip -d $(BINPATH)

$(OBJDIR)/ :
	mkdir -p $(OBJDIR)/

$(BINDIR)/ :
	mkdir -p $(BINDIR)/

$(OBJDIR)/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR)/ $(BINPATH)
