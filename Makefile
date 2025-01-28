#EXP=BES3LUMI
#EXP=KLOEII
#EXP=BES3
#EXP=KLOEI
#EXP=B
#EXP=CMD
EXP=DEFAULT

DEXP=
ifneq ($(EXP),DEFAULT)
DEXP=-D$(EXP)
endif
ifeq ($(EXP),DEFAULT)
EXP=
endif

RECOLA=
## RECOLA=-DRECOLA

COLLIER=
COLLIER=-DCOLLIER
LT =
#LT =-DLOOPTOOLS

LIBFILES=
LTVER=2.16
LTDIR=LoopTools-$(LTVER)
ifeq ($(LT),-DLOOPTOOLS)
LIBFILES += $(LTDIR)/lib64/libooptools.a
LTINC=-I$(LTDIR)/include
LTSTRING =-looptools
ifdef QUAD
  LTSTRING =-looptools-quad
endif
endif

FRED=alphaQEDc19/libfred.a
FRED=

CLLDIR = collier/COLLIER-1.2.8/
ifeq ($(COLLIER),-DCOLLIER)
  LIBFILES += $(CLLDIR)/lib/libcollier.a
  CLLMOD = -I$(CLLDIR)include/
  CLLLIB = -L$(CLLDIR) -lcollier
endif

EXE = babayaga$(EXP)
PAW = -lpawlib -lpacklib -lmathlib
CERNLIB = -L/cern/pro/lib
PAW =
CERNLIB =
PAWyn = ypaw

F77 = gfortran
# tune for your processor
FFLAGS = -O3
FFLAGS = -O3 -march=native -mtune=native -ffast-math

#FFLAGS = -O3 -march=native -mtune=native -ffast-math

## C ranlux optimizations, AVX2 faster but only on recent CPUS, SSE2 slower but present on
## any neowulf node. If both are specified, -DSSE2 is ignored (and won't run on old neowulf nodes).
RLXOPT=-DSSE2 -DAVX2

CC = gcc

FPCHECK = trapfpe.c
FPCHECK =

EXTRADEPS = #Makefile

VPHLMNT=vp_hlmnt_v2_1
VPHLMNT=vp_hlmnt_v2_1_1
VPHLMNT=vp_hlmnt_v2_2

ifeq ($(RECOLA),-DRECOLA)
#        RCLDIR = recolas/recola2-collier-2.1.6
#        RCLDIR = recolas/recola-collier-1.3.7
        RCLVERSION = 1.4.0
        RCLDIR = recolas/recola-collier-$(RCLVERSION)
        RCLMOD = -I$(RCLDIR)/include/ -I$(RCLDIR)/recola-$(RCLVERSION)/include/
        RCLLIB = -L$(RCLDIR)/install/lib/ -lrecola -lcollier
# the following only for recola2
# -lmodelfile
endif

OBJECTS = main.o cuts.o sv.o matrix_model.o mapmomenta.o loops.o ffpi.o\
          routines.o sampling.o phasespace.o distributions.o $(VPHLMNT).o\
          hadr5n16.o hadr5n09.o hadr5x23.o userinterface.o intpl.o Rteubner.o $(PAWyn).o\
          hadr5n17.o hadr5n12.o c_rnlx_interface.o ranlux_common.o ranlxd.o ranlxs.o ranlux.o recola_int.o
#\
#          Acp.o Aint.o Asu3.o interface.o mtx_eeenudbb.o 

F77 += $(FFLAGS)

default: $(EXE)

SAVEDIR = release
RELEASEDIR = BabaYaga
pack: # use only to release BABAYAGA
	mkdir -p $(RELEASEDIR)/form &&\
	mkdir -p $(RELEASEDIR)/ffpi-models/ &&\
	mkdir -p $(RELEASEDIR)/c_ranlux &&\
	mkdir -p $(RELEASEDIR)/collier/ &&\
	cp -ra oneloop/ $(RELEASEDIR)/oneloop &&\
	cp -ra LoopTools-$(LTVER)-clean/ $(RELEASEDIR)/LoopTools-$(LTVER) &&\
	cp -ra collier/COLLIER-1.2.8-clean/ $(RELEASEDIR)/collier/COLLIER-1.2.8 &&\
        cp -r Makefile README input-test-pack input-gg-test shared.F trapfpe.c\
        main.F loops.F cuts.F sv.F matrix_model.F mapmomenta.F vpol_novosibirsk.dat\
        routines.F invariants.h sampling.F phasespace.F distributions.F $(VPHLMNT).F\
        hadr5n16.F hadr5n09.F ranlux.F userinterface.F intpl.F Rteubner.F strong2020common.F\
        hadr5n17.F hadr5n12.F hadr5x23.F Acp.F Aint.F Asu3.F interface.F mtx_eeenudbb.F ffpi.F\
        recola_int.F npaw.F ypaw.F dalhadshigh17.F dalhadslow17.F dalhadt17.F dalhadthigh17.F\
        $(RELEASEDIR) &&\
	cp form/*.[fF] $(RELEASEDIR)/form/ &&\
	cp ffpi-models/* $(RELEASEDIR)/ffpi-models/ &&\
	cp c_ranlux/* $(RELEASEDIR)/c_ranlux/ &&\
	mkdir $(SAVEDIR) ;\
	tar cjvf $(SAVEDIR)/babayaga.tar.bz2 $(RELEASEDIR)/ &&\
	rm -rf $(RELEASEDIR)
clean:
	rm -f $(OBJECTS) *.a
deepclean:
	rm -rf $(OBJECTS) *.a $(EXE) *~ form/*~ $(LTDIR)/lib64/ $(LTDIR)/build/ $(LTDIR)/bin $(LTDIR)/lib64/ $(LTDIR)/include/ $(CLLDIR)build/ $(CLLDIR)lib/ $(CLLDIR)collierConfig.cmake $(CLLDIR)collierConfigVersion.cmake $(CLLDIR)include/*

# C version of ranlux by Martin Luscher
# http://luscher.web.cern.ch/luscher/ranlux/index.html
c_rnlx_interface.o: c_ranlux/c_rnlx_interface.c $(EXTRADEPS) Makefile
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/c_rnlx_interface.c
ranlxd.o: c_ranlux/ranlxd.c $(EXTRADEPS) Makefile 
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlxd.c
ranlxs.o: c_ranlux/ranlxs.c $(EXTRADEPS) Makefile
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlxs.c
ranlux_common.o: c_ranlux/ranlux_common.c $(EXTRADEPS) Makefile
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlux_common.c
####


# source files
main.o: main.F $(EXTRADEPS)
	$(F77) -c main.F
cuts.o: cuts.F $(EXTRADEPS) Makefile strong2020common.F
	$(F77) $(DEXP) -c cuts.F
matrix_model.o: matrix_model.F  form/formme.F form/ee3gexact.F form/borngg.F form/formmemm.F $(EXTRADEPS)
	$(F77) -c matrix_model.F
sv.o: sv.F $(EXTRADEPS) Makefile $(CLLDIR)/lib/libcollier.a
	$(F77) -c $(COLLIER) $(CLLMOD) $(LT) $(LTINC) sv.F
2loop.o: 2loop.F $(EXTRADEPS)
	$(F77) -c 2loop.F
hplog.o: 2loop/hplog.F $(EXTRADEPS)
	$(F77) -c 2loop/hplog.F
userinterface.o: userinterface.F $(EXTRADEPS) Makefile
	$(F77) $(DEXP) $(RECOLA) -c userinterface.F
phasespace.o: phasespace.F $(EXTRADEPS) strong2020common.F Makefile
	$(F77) -c $(DEXP) phasespace.F
hadr5n16.o: hadr5n16.F $(EXTRADEPS)
	$(F77) -c hadr5n16.F
hadr5n17.o: hadr5n17.F $(EXTRADEPS)
	$(F77) -c hadr5n17.F
ffpi.o: ffpi.F $(EXTRADEPS)
	$(F77) -c ffpi.F
hadr5x23.o: hadr5x23.F $(EXTRADEPS)
	$(F77) -c hadr5x23.F
hadr5n12.o: hadr5n12.F $(EXTRADEPS)
	$(F77) -c hadr5n12.F
hadr5n09.o: hadr5n09.F $(EXTRADEPS)
	$(F77) -c hadr5n09.F
intpl.o: intpl.F $(EXTRADEPS)
	$(F77) -c intpl.F
$(VPHLMNT).o: $(VPHLMNT).F $(EXTRADEPS)
	$(F77) -c $(VPHLMNT).F
Rteubner.o: Rteubner.F $(EXTRADEPS)
	$(F77) -c Rteubner.F
mapmomenta.o: mapmomenta.F $(EXTRADEPS)
	$(F77) -c mapmomenta.F
sampling.o: sampling.F $(EXTRADEPS) Makefile strong2020common.F
	$(F77) -c $(DEXP) sampling.F
loops.o: loops.F $(EXTRADEPS) Makefile $(CLLDIR)/lib/libcollier.a
	$(F77) -c $(COLLIER) $(CLLMOD) $(LT) $(LTINC) loops.F
routines.o: routines.F $(EXTRADEPS)
	$(F77) $(COLLIER) $(CLLMOD) -c routines.F
distributions.o: distributions.F shared.F $(EXTRADEPS) Makefile strong2020common.F
	$(F77) $(DEXP) -c distributions.F
ranlux.o: ranlux.F $(EXTRADEPS)
	$(F77) -c ranlux.F
$(PAWyn).o: $(PAWyn).F $(EXTRADEPS)
	$(F77) -c $(PAWyn).F
recola_int.o: recola_int.F $(EXTRADEPS) Makefile
	$(F77) $(RECOLA) -c $(RCLMOD) recola_int.F

looptools: $(LTDIR)/lib64/libooptools.a
$(LTDIR)/lib64/libooptools.a:
	@echo "Building LoopTools"
	@echo " "
	cd $(LTDIR) && ./configure --prefix=. && make -j`nproc` && make install

collier: $(CLLDIR)/lib/libcollier.a
$(CLLDIR)/lib/libcollier.a:
	@echo "Building Collier"
	@echo " "
	mkdir -p $(CLLDIR)/build/
	cd $(CLLDIR)/build/ && cmake .. -DCMAKE_INSTALL_PREFIX=.. -Dstatic=ON && make && make install

### BEGIN ALPHA
# Acp.o: Acp.F $(EXTRADEPS)
# 	$(F77) -fno-automatic -c Acp.F
# Aint.o: Aint.F $(EXTRADEPS)
# 	$(F77) -fno-automatic -c Aint.F
# Asu3.o: Asu3.F $(EXTRADEPS)
# 	$(F77) -fno-automatic -c Asu3.F
# mtx_eeenudbb.o: mtx_eeenudbb.F $(EXTRADEPS)
# 	$(F77) -fno-automatic -c mtx_eeenudbb.F
# interface.o: interface.F $(EXTRADEPS)
# 	$(F77) -fno-automatic -c interface.F
### END ALPHA

# babayaga library
libbabayaga.a: $(OBJECTS)
	ar cr libbabayaga.a $(OBJECTS)

# executable
$(EXE): $(LIBFILES) libbabayaga.a
	$(F77) main.o $(FPCHECK) -o $(EXE) -L. -lbabayaga $(FRED) $(CLLLIB) $(RCLLIB) -L$(LTDIR)/lib64 $(LTSTRING) $(RCLEXTRA)
