# NCIPLOT v 3.0

# Compiler settings

include Makefile.inc

# Objects and binaries

BINS=nciplot
BINS_dbg=nciplot_dbg
OBJS=param.o tools_io.o tools_math.o reader.o props.o nciplot.o
LIBS=
INCLUDE=

%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.o: %.f
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

# General targets

all: $(BINS)

debug:
	DEBUG=1 $(MAKE) $(BINS_dbg)

clean:
	rm -f core *.mod *.o

mrproper:
	rm -rf core *.mod *.o ../bin $(BINS) $(BINS_dbg)

nciplot: $(OBJS) $(LIBS)
	$(FC) -o nciplot $(LDFLAGS) $(OBJS) $(LIBS)
	@mkdir -p ../bin; \
	ln -frs nciplot ../bin/.

nciplot_dbg: $(OBJS) $(LIBS)
	$(FC) -o nciplot_dbg $(LDFLAGS) $(OBJS) $(LIBS)
	@mkdir -p ../bin; \
	ln -frs nciplot_dbg ../bin/.

# Object dependencies

nciplot.o props.o reader.o tools_io.o : param.mod
nciplot.o props.o reader.o : tools_io.mod
nciplot.o props.o : reader.mod
nciplot.o : props.mod
nciplot.o : tools_math.o

