######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
IBAMR_SRC_DIR   = $(HOME)/sfw/ibamr/IBAMR
IBAMR_BUILD_DIR = $(HOME)/sfw/ibamr/ibamr-objs-dbg

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

######################################################################
## Build the application.

PDIM = 3
OBJS = boundary_condition_util.o CirculationModel.o FourierBodyForce.o main.o 
	
main3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS)
	$(CXX) -o main3d $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

clean:
	$(RM) main3d 
	$(RM) *.o 