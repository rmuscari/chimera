SRC = src
OBJ = obj
MOD = $(SRC)/mod
OUT = overset3
TARGET = $(HOME)/opt/bin/$(OUT)
TECLIB = /opt/tecplot/2009/lib/tecio64.a -lstdc++

#===============================================================================
# Compilatore e opzioni di compilazione

#-----------------------
ifdef pgi
COFLAGS =
LDFLAGS =
FCFLAGS =
#
ifdef mpi
FORT = mpif90
else
FORT = pgf90
endif
#
ifdef deb
COFLAGS = -g -Mbounds
else
COFLAGS = -fast -Mcache_align -Mdalign -Mdclchk -Kieee -Ktrap=fp
endif
#
ifdef omp
COFLAGS := $(COFLAGS) -mp
endif
#
FCFLAGS = $(COFLAGS) -module $(OBJ) -c
LDFLAGS := $(COFLAGS) $(LDFLAGS)
endif

#-----------------------
ifdef intel
OBJ := $(OBJ)/intel
COFLAGS =
LDFLAGS =
FCFLAGS =
#
ifdef mpi
FORT = mpif90
else
FORT = ifort
endif
#
ifdef deb
COFLAGS = -g -check all
else
COFLAGS = -O3
endif
#
ifdef omp
COFLAGS := $(COFLAGS) -openmp
endif
#
FCFLAGS = $(COFLAGS) -fpe0 -module $(OBJ) -warn all -c
LDFLAGS := $(COFLAGS) $(LDFLAGS)
endif

#-----------------------
ifdef gnu
OBJ := $(OBJ)/gnu
COFLAGS =
LDFLAGS =
FCFLAGS =
#
ifdef mpi
FORT = mpif90
else
FORT = gfortran
endif
#
ifdef deb
COFLAGS = -g -fbounds-check
else
COFLAGS = -O2
endif
#
ifdef omp
LDFLAGS := $(LDFLAGS) -fopenmp
endif
#
FCFLAGS = $(COFLAGS) -Wall -J$(OBJ) -I$(OBJ) -fimplicit-none -c
LDFLAGS := $(COFLAGS) $(LDFLAGS)
endif

#===============================================================================
# file sorgenti
SRC1 = $(MOD)/xnavis_prec.f90
SRC2 = $(MOD)/xnavis_param.f90 $(MOD)/xnavis_chimera.f90 \
       $(MOD)/xnavis_point.f90
SRC3 = $(filter-out $(SRC1) $(SRC2), $(wildcard $(MOD)/*.f90))
SRC4 = $(wildcard $(SRC)/*.f90)

OBJ1 = $(patsubst $(MOD)/%.f90,$(OBJ)/%.o,$(SRC1))
OBJ2 = $(patsubst $(MOD)/%.f90,$(OBJ)/%.o,$(SRC2))
OBJ3 = $(patsubst $(MOD)/%.f90,$(OBJ)/%.o,$(SRC3))
OBJ4 = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SRC4))

#===============================================================================
# Rules

overset : $(OBJ4)
	@echo
	@rm -f $(TARGET)
	@echo Assemblaggio $(TARGET)
	@$(FORT) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(LDFLAGS) $(TECLIB) -o $(TARGET)
	@echo

$(OBJ4) : $(OBJ)/%.o: $(SRC)/%.f90 $(OBJ3)
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(SRC)/$(<F) -o $@

$(OBJ3) : $(OBJ)/%.o: $(MOD)/%.f90 $(OBJ2)
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

$(OBJ2) : $(OBJ)/%.o: $(MOD)/%.f90 $(OBJ1)
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

$(OBJ1) : $(OBJ)/%.o: $(MOD)/%.f90 /bin/ls
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

#
# /bin/ls è, in teoria, un qualunque file che non viene mai aggiornato
#
# la sequenza così costruita serve a eseguire "opzioni" prima di cominciare la
# compilazione
#
# se non viene inserito il passaggio per /bin/ls e si fa dipendere OBJ1
# direttamente da "opzioni", ogni volta che si esegue il make viene ricompilato
# tutto il codice
#
/bin/ls : opzioni

opzioni :
	@echo
	@echo ++++++++++++++++++++
	@echo "Sorgenti in $(SRC)"
	@echo
	@echo "Compilatore:"
	@echo "   " `type -p $(FORT)`
	@echo
	@echo "Opzioni di compilazione:"
	@echo "   $(FCFLAGS)"
	@echo
	@echo "Opzioni di linking:"
	@echo "   $(LDFLAGS)"
	@echo ++++++++++++++++++++
	@echo

backup :
	@DIR=BACKUP/`date +"%y%m%d"` ; \
      if test -d $${DIR} ; then echo "Directory di backup già esistente" ; \
      else cp -pR $(SRC) $${DIR} 1>/dev/null 2>/dev/null ; fi

clean :
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	rm -f $(TARGET)

cleanall :
	find $(OBJ) -type f -exec rm {} \;
	find $(HOME)/opt/bin -name "$(OUT)*" -exec rm {} \;
