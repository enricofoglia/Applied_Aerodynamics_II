
SRC = ../src
# PLTOBJ = /home/codes/Xplot/libPlt.a
PLTOBJ = /Users/enricofoglia/Downloads/Xfoil/plotlib/libPlt_gDP.a # apparently not essential

#========================================
#  Default compiler flags
#FC = f77
#FLG = -O
#PLTLIB = -lX11
#FTNLIB =

#========================================
#  Uncomment for Intel Fortran Compiler
#FC = ifort
#FLG = -O -fpe0 -CB
#PLTLIB = -L/usr/X11R6/lib -lX11 # linux setting
#FTNLIB = 

#========================================
# GNU Fortran compiler
FC = gfortran
FLG = -O2 -fdefault-real-8 # double precision
PLTLIB = -L/usr/local/X11/lib -lX11 # mac setting
FTNLIB = 

#  Additional IFC stuff needed on MIT's Athena system
#FC = ifc
#FTNLIB = -Vaxlib /usr/lib/C-ctype.o /usr/lib/C_name.o /usr/lib/ctype-info.o

#========================================


OS: osgen osmap.o

clean:
	/bin/rm *.o osgen

osgen: osgen.o io.o spline.o getarg0.o
	$(FC) -o osgen osgen.o io.o spline.o getarg0.o $(FTNLIB)
osgen1: osgen1.o io.o spline.o getarg0.o
	$(FC) -o osgen1 osgen1.o io.o spline.o getarg0.o $(FTNLIB)
fscorr: fscorr.o fs.o
	$(FC) -o fscorr fscorr.o fs.o $(PLTOBJ) $(PLTLIB)
fsrun: fsrun.o fs.o
	$(FC) -o fsrun fsrun.o fs.o $(PLTOBJ) $(PLTLIB)
fs4run: fsrun.o fs4.o
	$(FC) -o fs4run fsrun.o fs4.o $(PLTOBJ) $(PLTLIB)
intai: intai.o osmap.o plutil.o
	$(FC) -o intai intai.o osmap.o plutil.o $(PLTOBJ) $(PLTLIB)
mappl1: mappl1.o ask1.o conlab.o io.o getarg0.o
	$(FC) -o mappl1 mappl1.o ask1.o conlab.o io.o getarg0.o $(PLTOBJ) $(PLTLIB)
mappl0: mappl0.o ask1.o conlab.o io.o getarg0.o
	$(FC) -o mappl0 mappl0.o ask1.o conlab.o io.o getarg0.o $(PLTOBJ) $(PLTLIB)
mappl0x: mappl0x.o ask1.o conlab.o io.o getarg0.o
	$(FC) -o mappl0x mappl0x.o ask1.o conlab.o io.o getarg0.o $(PLTOBJ) $(PLTLIB)
ncorr: ncorr.o fs.o
	$(FC) -o ncorr ncorr.o fs.o $(PLTOBJ) $(PLTLIB)
osrun: osrun.o fs.o orrs.o ospres.o plutil.o userio.o
	$(FC) -o osrun osrun.o fs.o orrs.o ospres.o \
plutil.o userio.o $(PLTOBJ) $(PLTLIB) $(FTNLIB)
osseq: osseq.o fs.o orrs.o ospres.o plutil.o userio.o
	$(FC) -o osseq osseq.o fs.o orrs.o ospres.o \
plutil.o userio.o $(PLTOBJ) $(PLTLIB) $(FTNLIB)
pfplot: pfplot.o fs.o ask1.o
	$(FC) -o pfplot pfplot.o fs.o ask1.o $(PLTOBJ) $(PLTLIB)
roll: roll.o
	$(FC) -o roll roll.o $(PLTOBJ) $(PLTLIB) $(FTNLIB)
as2bi: as2bi.o io.o getarg0.o
	$(FC) -o as2bi as2bi.o io.o getarg0.o $(FTNLIB)
bi2as: bi2as.o io.o getarg0.o
	$(FC) -o bi2as bi2as.o io.o getarg0.o $(FTNLIB)
bi2bi: bi2bi.o io.o getarg0.o
	$(FC) -o bi2bi bi2bi.o io.o getarg0.o $(FTNLIB)
otest: otest.o osmap.o
	$(FC) -o otest otest.o osmap.o
osweep: osweep.o osmap.o
	$(FC) -o osweep osweep.o osmap.o

osgen.o: $(SRC)/osgen.f
	$(FC) -c $(FLG) $(SRC)/osgen.f
osgen1.o: $(SRC)/osgen1.f
	$(FC) -c $(FLG) $(SRC)/osgen1.f
orrs.o: $(SRC)/orrs.f $(SRC)/ORRS.INC
	$(FC) -c $(FLG) $(SRC)/orrs.f
ospres.o: $(SRC)/ospres.f $(SRC)/OSPRES.INC
	$(FC) -c $(FLG) $(SRC)/ospres.f
pfplot.o: $(SRC)/pfplot.f
	$(FC) -c $(FLG) $(SRC)/pfplot.f
fscorr.o: $(SRC)/fscorr.f
	$(FC) -c $(FLG) $(SRC)/fscorr.f
fsrun.o: $(SRC)/fsrun.f
	$(FC) -c $(FLG) $(SRC)/fsrun.f
fs.o: $(SRC)/fs.f
	$(FC) -c $(FLG) $(SRC)/fs.f
fs4.o: $(SRC)/fs4.f
	$(FC) -c $(FLG) $(SRC)/fs4.f
io.o: $(SRC)/io.f
	$(FC) -c $(FLG) $(SRC)/io.f
intai.o: $(SRC)/intai.f
	$(FC) -c $(FLG) $(SRC)/intai.f
mappl1.o: $(SRC)/mappl1.f
	$(FC) -c $(FLG) $(SRC)/mappl1.f
mappl0.o: $(SRC)/mappl0.f
	$(FC) -c $(FLG) $(SRC)/mappl0.f
mappl0x.o: $(SRC)/mappl0x.f
	$(FC) -c $(FLG) $(SRC)/mappl0x.f
ncorr.o: $(SRC)/ncorr.f
	$(FC) -c $(FLG) $(SRC)/ncorr.f
osmap.o: $(SRC)/osmap.f
	$(FC) -c $(FLG) $(SRC)/osmap.f
osrun.o: $(SRC)/osrun.f
	$(FC) -c $(FLG) $(SRC)/osrun.f
osseq.o: $(SRC)/osseq.f
	$(FC) -c $(FLG) $(SRC)/osseq.f
roll.o: $(SRC)/roll.f
	$(FC) -c $(FLG) $(SRC)/roll.f
otest.o: $(SRC)/otest.f
	$(FC) -c $(FLG) $(SRC)/otest.f
osweep.o: $(SRC)/osweep.f
	$(FC) -c $(FLG) $(SRC)/osweep.f

spline.o: $(SRC)/spline.f
	$(FC) -c $(FLG) $(SRC)/spline.f
conlab.o: $(SRC)/conlab.f
	$(FC) -c $(FLG) $(SRC)/conlab.f
plutil.o: $(SRC)/plutil.f
	$(FC) -c $(FLG) $(SRC)/plutil.f
ask1.o: $(SRC)/ask1.f
	$(FC) -c $(FLG) $(SRC)/ask1.f
userio.o: $(SRC)/userio.f
	$(FC) -c $(FLG) $(SRC)/userio.f
getarg0.o: $(SRC)/getarg0.f
	$(FC) -c $(FLG) $(SRC)/getarg0.f

