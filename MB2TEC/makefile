SRC = src
OBJ = obj
MOD = $(SRC)
OUT = mbvtec
TARGET = $(HOME)/opt/bin/$(OUT)
TECLIB = /opt/tecplot/2009/lib/tecio64.a -lstdc++

#===============================================================================
# Compilatore e opzioni di compilazione

#-----------------------
FORT = pgf90
LDFLAGS = -fast -Mcache_align -Mdalign -Mdclchk -Kieee -Ktrap=fp
#
ifdef omp
LDFLAGS := $(LDFLAGS) -mp
endif
#
ifdef deb
LDFLAGS = -g # -Mbounds
endif
#
FCFLAGS = $(LDFLAGS) -module $(OBJ) -c

#-----------------------
ifdef gnu
OBJ := $(OBJ)/gnu
FORT = gfortran
LDFLAGS = -O2
#
ifdef omp
LDFLAGS := $(LDFLAGS) -fopenmp
endif
#
ifdef deb
LDFLAGS = -g -fbounds-check
endif
#
FCFLAGS = $(LDFLAGS) -Wall -J$(OBJ) -I$(OBJ) -fimplicit-none -c
endif

#-----------------------
ifdef intel
OBJ := $(OBJ)/intel
FORT = ifort
LDFLAGS = -O3
#
ifdef omp
LDFLAGS := $(LDFLAGS) -openmp
endif
#
ifdef deb
LDFLAGS = -g -C -check bounds
endif
#
FCFLAGS = $(LDFLAGS) -module $(OBJ) -warn all -c
endif

#===============================================================================
# file sorgenti
SRC1 = $(MOD)/prec.f90
SRC2 = $(MOD)/moddef.f90
SRC3 = $(MOD)/modpar.f90
SRC4 = $(filter-out $(SRC1) $(SRC2) $(SRC3), $(wildcard $(SRC)/*.f90))

OBJ1 = $(patsubst $(MOD)/%.f90,$(OBJ)/%.o,$(SRC1))
OBJ2 = $(patsubst $(MOD)/%.f90,$(OBJ)/%.o,$(SRC2))
OBJ3 = $(patsubst $(MOD)/%.f90,$(OBJ)/%.o,$(SRC3))
OBJ4 = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SRC4))

#===============================================================================
# Rules

mbvtec : $(OBJ4)
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
	@DIR=../BACKUP/$(OUT).`date +"%y%m%d"` ; \
      if test -d $${DIR} ; then echo "Directory di backup già esistente" ; \
      else cp -pR $(SRC) $${DIR} 1>/dev/null 2>/dev/null ; fi

clean :
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	rm -f $(TARGET)

cleanall :
	find $(OBJ) -type f -exec rm {} \;
	find $(HOME)/opt/bin -name "$(OUT)*" -exec rm {} \;
