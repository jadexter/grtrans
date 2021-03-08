include Makefile.top

ifeq ($(USEINTEL),1)
FC=ifort
AR=xiar
ARFLAGS=rcv
PHIFLAGS=-mmic
OMP = -qopenmp #-mkl
OMPLIB = -liomp5 -lpthread
FCNAME=intelem
ifeq ($(DEBUG),1)
OTHERFLAGS = $(OMP) -extend-source 132 -fp-model source -O3 -ipo -xHost -fPIC
else
OTHERFLAGS = $(OMP) -extend-source 132 -O3 -ipo -xHost -parallel -fp-model source -fPIC #-mkl
#OTHERFLAGS = $(OMP) -extend-source 132 -ipo -O2 -static -no-prec-div -fp-model source -fPIC -xHost
endif
endif

ifeq ($(USEGNU),1)
FC=gfortran
AR=ar
ARFLAGS=rcv
OMP = -fopenmp
OMPLIB = -lgomp
FCNAME=gnu95
PHIFLAGS= 
ifeq ($(DEBUG),1)
OTHERFLAGS = -ffixed-line-length-132 -g -pg -W -Wall -fPIC
#OMP=
#OMPLIB=
else
OTHERFLAGS = -ffixed-line-length-132 -fopenmp -O3 -fPIC
endif
endif

# compile for co-processor
ifeq ($(XEONPHI),1)
OTHERFLAGS += $(PHIFLAGS)
LIBS=$(GRTRANSDIR)/libcfitsio_phi.a
LINKFLAGS=$(OMP) $(PHIFLAGS)
else
LINKFLAGS=$(OMP)
LIBS=-L$(CFITSIODIR) -lcfitsio
endif
ifeq ($(PROFILE),1)
OTHERFLAGS += -pg
LINKFLAGS += -pg
endif

.DEFAULT:
	-echo $@ does not exist.
all: grtrans libgrtrans pgrtrans
lib: libgrtrans
public: grtrans_public
pgrtrans: pgrtrans
phi: grtrans_phi
commit: 
	python commit_grtrans.py

odepack_aux.o: ./odepack_aux.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./odepack_aux.f
interpolate_aux.o: ./interpolate_aux.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./interpolate_aux.f
geokerr_wrapper.o: ./geokerr_wrapper.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./geokerr_wrapper.f
camera.o: ./camera.f90 fits.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./camera.f90
chandra_tab24.o: ./chandra_tab24.f90 interpolate.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./chandra_tab24.f90
fits.o: fits.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fits.f90
class_four_vector.o: ./class_four_vector.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./class_four_vector.f90
emis.o: ./emis.f90 math.o polsynchemis.o chandra_tab24.o calc_maxjutt.o calc_maxcomp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./emis.f90
fluid.o: ./fluid.f90 class_four_vector.o phys_constants.o interpolate.o kerr.o fluid_model_sphacc.o fluid_model_phatdisk.o fluid_model_thindisk.o fluid_model_ffjet.o fluid_model_numdisk.o fluid_model_hotspot.o fluid_model_hotspot_schnittman.o fluid_model_harm.o fluid_model_harmpi.o fluid_model_iharm.o fluid_model_koral.o fluid_model_koral3d.o fluid_model_harm3d.o fluid_model_thickdisk.o fluid_model_mb09.o fluid_model_sariaf.o fluid_model_powerlaw.o fluid_model_toy.o calcgmin.o 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid.f90
fluid_model_sphacc.o: ./fluid_model_sphacc.f90 class_four_vector.o phys_constants.o interpolate.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_sphacc.f90
fluid_model_ffjet.o: ./fluid_model_ffjet.f90 class_four_vector.o phys_constants.o interpolate.o kerr.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_ffjet.f90
fluid_model_thindisk.o: ./fluid_model_thindisk.f90 class_four_vector.o phys_constants.o kerr.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_thindisk.f90
fluid_model_phatdisk.o: ./fluid_model_phatdisk.f90 fluid_model_thindisk.o class_four_vector.o phys_constants.o 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_phatdisk.f90
fluid_model_numdisk.o: ./fluid_model_numdisk.f90 fluid_model_thindisk.o class_four_vector.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_numdisk.f90
fluid_model_hotspot.o: ./fluid_model_hotspot.f90 fluid_model_thindisk.o class_four_vector.o kerr.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_hotspot.f90
fluid_model_hotspot_schnittman.o: ./fluid_model_hotspot_schnittman.f90 class_four_vector.o kerr.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_hotspot_schnittman.f90
fluid_model_toy.o: ./fluid_model_toy.f90 phys_constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_toy.f90

fluid_model_harm.o: ./fluid_model_harm.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_harm.f90

fluid_model_koral.o: ./fluid_model_koral.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_koral.f90

fluid_model_koral3d.o: ./fluid_model_koral3d.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_koral3d.f90

fluid_model_harm3d.o: ./fluid_model_harm3d.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_harm3d.f90

fluid_model_harmpi.o: ./fluid_model_harmpi.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_harmpi.f90

fluid_model_iharm.o: ./fluid_model_iharm.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_iharm.f90

fluid_model_thickdisk.o: ./fluid_model_thickdisk.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_thickdisk.f90
fluid_model_mb09.o: ./fluid_model_mb09.f90 class_four_vector.o interpolate.o kerr.o phys_constants.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_mb09.f90
fluid_model_sariaf.o: ./fluid_model_sariaf.f90 phys_constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_sariaf.f90

fluid_model_powerlaw.o: ./fluid_model_powerlaw.f90 phys_constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluid_model_powerlaw.f90

calc_maxjutt.o: ./calc_maxjutt.f90 polsynchemis.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./calc_maxjutt.f90

calc_maxcomp.o: ./calc_maxcomp.f90 polsynchemis.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./calc_maxcomp.f90

calcgmin.o: ./calcgmin.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./calcgmin.f90
geodesics.o: ./geodesics.f90 class_four_vector.o interpolate.o kerr.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./geodesics.f90
grtrans_driver.o: ./grtrans_driver.f90 fluid.o radtrans_integrate.o rad_trans.o geodesics.o emis.o odepack.o camera.o interpolate.o phys_constants.o kerr.o read_inputs.o chandra_tab24.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./grtrans_driver.f90
pgrtrans.o: ./pgrtrans.f90 grtrans_driver.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./pgrtrans.f90
grtrans.o: ./grtrans.f90 grtrans_driver.o pgrtrans.o read_inputs.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./grtrans.f90
grtrans_program.o: ./grtrans_program.f90 grtrans.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./grtrans_program.f90
interpolate.o: ./interpolate.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./interpolate.f90
kerr.o: ./kerr.f90 class_four_vector.o math.o phys_constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./kerr.f90
math.o: ./math.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./math.f90
odepack.o: ./odepack.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./odepack.f90
phys_constants.o: ./phys_constants.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./phys_constants.f90
bessel.o: ./bessel.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./bessel.f90
polsynchemis.o: ./polsynchemis.f90 phys_constants.o bessel.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./polsynchemis.f90
rad_trans.o: ./rad_trans.f90 read_inputs.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./rad_trans.f90
radtrans_integrate.o: ./radtrans_integrate.f90 odepack.o interpolate.o math.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./radtrans_integrate.f90
read_inputs.o: ./read_inputs.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./read_inputs.f90
SRC = ./odepack_aux.f ./interpolate_aux.f ./geokerr_wrapper.f ./interpolate.f90 ./read_inputs.f90 ./class_four_vector.f90 ./odepack.f90 ./emis.f90 ./rad_trans.f90 ./kerr.f90 ./fluid.f90 ./grtrans_driver.f90 ./calcgmin.f90 ./calc_maxjutt.f90 ./calc_maxcomp.f90 ./fluid_model_sphacc.f90 ./fluid_model_thindisk.f90 ./fluid_model_phatdisk.f90 ./fluid_model_numdisk.f90 ./fluid_model_hotspot.f90 ./fluid_model_hotspot_schnittman.f90 ./fluid_model_harm.f90 ./fluid_model_koral.f90 ./fluid_model_koral3d.f90 ./fluid_model_harmpi.f90 ./fluid_model_iharm.f90 ./fluid_model_harm3d.f90 ./fluid_model_thickdisk.f90 ./fluid_model_mb09.f90 ./fluid_model_sariaf.f90 ./fluid_model_powerlaw.f90 ./fluid_model_ffjet.f90 ./fluid_model_toy.f90          \
./polsynchemis.f90 ./geodesics.f90 ./math.f90 ./camera.f90 ./phys_constants.f90 ./bessel.f90 ./grtrans.f90 ./grtrans_program.f90
OBJ = odepack_aux.o interpolate_aux.o geokerr_wrapper.o interpolate.o read_inputs.o class_four_vector.o odepack.o emis.o rad_trans.o radtrans_integrate.o kerr.o fluid.o grtrans_driver.o calcgmin.o calc_maxjutt.o calc_maxcomp.o fluid_model_sphacc.o pgrtrans.o\
polsynchemis.o geodesics.o math.o camera.o phys_constants.o bessel.o fits.o fluid_model_thindisk.o fluid_model_phatdisk.o fluid_model_numdisk.o fluid_model_hotspot.o fluid_model_hotspot_schnittman.o fluid_model_ffjet.o fluid_model_harm.o fluid_model_koral.o fluid_model_harmpi.o fluid_model_iharm.o fluid_model_koral3d.o fluid_model_harm3d.o fluid_model_thickdisk.o fluid_model_mb09.o fluid_model_sariaf.o fluid_model_toy.o fluid_model_powerlaw.o grtrans.o grtrans_program.o chandra_tab24.o
OBJPUB = odepack_aux.o interpolate_aux.o geokerr_wrapper.o interpolate.o read_inputs.o class_four_vector.o odepack.o emis.o rad_trans.o radtrans_integrate.o kerr.o fluid.o grtrans_driver.o calcgmin.o calc_maxjutt.o calc_maxcomp.o fluid_model_sphacc.o  \
polsynchemis.o geodesics.o math.o camera.o phys_constants.o bessel.o fits.o fluid_model_thindisk.o fluid_model_phatdisk.o fluid_model_numdisk.o fluid_model_hotspot.o fluid_model_hotspot_schnittman.o fluid_model_ffjet.o fluid_model_harm.o fluid_model_koral.o fluid_model_harmpi.o fluid_model_iharm.o fluid_model_koral3d.o fluid_model_harm3d.o grtrans.o grtrans_program.o chandra_tab24.o

clean: neat
	-rm -f .grtrans.cppdefs $(OBJ) *.mod grtrans
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
grtrans_public: $(OBJPUB) 
	$(FC) -o grtrans  $(LIBS) $(OBJPUB)

grtrans: $(OBJ)
	$(FC) $(LINKFLAGS) -o grtrans  $(OBJ) $(LIBS)

libgrtrans: $(OBJ)
	$(AR) $(ARFLAGS) libgrtrans.a $(OBJ)

pgrtrans: libgrtrans.a
	f2py -c pgrtrans.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m pgrtrans $(LIBS) $(OMPLIB) $(GRTRANSDIR)/libgrtrans.a

radtrans_integrate:
	f2py -c radtrans_integrate.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m radtrans_integrate $(LIBS) $(OMPLIB) $(GRTRANSDIR)/libgrtrans.a

polsynchemis:
	f2py -c polsynchemis.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m polsynchemis $(LIBS) $(GRTRANSDIR)/libgrtrans.a

calcgmin::
	f2py -c calcgmin.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m calcgmin $(LIBS) $(OMPLIB) $(GRTRANSDIR)/libgrtrans.a

sariaf:
	f2py -c fluid_model_sariaf.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m sariaf $(LIBS) $(OMPLIB) $(GRTRANSDIR)/libgrtrans.a

maxjutt: polsynchemis.o
	f2py -c calc_maxjutt.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m calc_maxjutt $(LIBS) $(OMPLIB) polsynchemis.o

maxcomp: polsynchemis.o
	f2py -c calc_maxcomp.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m calc_maxcomp $(LIBS) $(OMPLIB) polsynchemis.o

geokerr: geokerr_wrapper.o
	f2py -c class_geokerr.f90 --fcompiler=$(FCNAME) --f90flags="$(FFLAGS) $(OTHERFLAGS)" -m geokerr $(LIBS) $(OMPLIB) geokerr_wrapper.o
