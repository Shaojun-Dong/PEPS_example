FC=mpif90
TensorVersion=3.5.8-4.5.0
PACKAGE= -larpack -llapack -lrefblas
LIBDIR=../../TNSP_Release/libTNSP/


TensorDIRName=Tensor-$(TensorVersion)
TENDIR=$(LIBDIR)/lib$(TensorDIRName)/


TensorDIRLink= -I$(TENDIR)
TensorlibLink= -L$(TENDIR)   -l$(TensorDIRName)    $(PACKAGE)



PEPS:= \
SimpleUpdate.o


all:   runPEPS  
runPEPS:	$(PEPS)
	$(FC)  $(openmp) -o $@ $^ $(TensorlibLink)
./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(Optionl)  $(TensorDIRLink)
clean:
	rm *.o *.mod runPEPS
