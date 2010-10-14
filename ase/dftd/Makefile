EXE = dftd_module

all: $(EXE)

SOURCES = dftd_functions.f copyc6.f gdisp.f
INTERFACE = dftd_functions.pyf

FC = gfortran
F2 = f2py

F2FLAGS = -c --noopt --fcompiler=gnu95



# %.o: %.f
#         $(FC)  -c $< -o $@
#  575  f2py dftd_functions.f -m dftd_module -h dftd_functions.pyf --overwrite-signature
#    576  f2py -c --noopt --fcompiler=gnu95 dftd_functions.pyf dftd_functions.f
#
%.pyf: %.f
	$(F2) $< -m $(EXE) -h $@ --overwrite-signature

$(EXE): $(INTERFACE)
	$(F2) $(F2FLAGS) $(INTERFACE) $(SOURCES)

#	$(F2) -c -m $(EXE) dftd_functions.f -DF2PY_REPORT_ON_ARRAY_COPY=1 --fcompiler=gnu95



clean:
	rm -f $(EXE:=.so) $(INTERFACE) *.o *.pyc
