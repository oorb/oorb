include ../make.config
include ../Makefile.include

FCOPTIONS = $(FCOPTIONS_DEB_GFORTRAN)


all: liboorb.$(LIBEXT)


oorb.o: oorb.f90
	$(FC) -c $(FCOPTIONS) $(FC_INC)../$(CLASSPATH) \
	$(FC_INC)../$(MODULEPATH) \
	$(FC_INC)../main \
	oorb.f90


wrapper.o: wrapper.f90 oorb.o
	$(FC) -c $(FCOPTIONS) $(FC_INC)../$(CLASSPATH) \
	$(FC_INC)../$(MODULEPATH) \
	$(FC_INC)../main \
	wrapper.f90


liboorb.$(LIBEXT): wrapper.o oorb.o
	$(FC) $(FCOPTIONS) -shared -o liboorb.$(LIBEXT) \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	oorb.o wrapper.o \
	../main/io.o \
	../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o





# Remove library and modules:
clean:
	rm -f *~ *.mod *.o *.$(LIBEXT) *.pyc *.pyc.$(OS)

