######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
IBAMR_SRC_DIR   = $(HOME)/sfw/ibamr/IBAMR
IBAMR_BUILD_DIR = $(HOME)/sfw/ibamr/ibamr-objs-opt

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

######################################################################
## Build the application.

PDIM = 3
		
main3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) main.o
	$(CXX) main.o -o main3d $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

clean:
	$(RM) main3d 
	$(RM) *.o 