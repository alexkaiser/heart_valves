######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
IBAMR_SRC_DIR   = $(HOME)/sfw/ibamr/IBAMR
IBAMR_BUILD_DIR = $(HOME)/sfw/ibamr/ibamr-objs-dbg_prescribed_motion

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

######################################################################
## Build the application.

PDIM = 3
OBJS = CirculationModel_with_lv.o \
       boundary_condition_util.o \
       CirculationModel.o \
       CirculationModel_RV_PA.o \
       FourierBodyForce.o \
       FeedbackForcer.o \
       pnpoly.o \
       main.o
	
main3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS)
	$(CXX) -o main3d $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

# CXXFLAGS += -DUSE_CIRC_MODEL_RV_PA
# OBJS_RV_PA = CirculationModel_with_lv.o boundary_condition_util.o CirculationModel.o FourierBodyForce.o CirculationModel_RV_PA.o pnpoly.o main.o
# main3d_rv_pa: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS_RV_PA)
# 	$(CXX) -o main3d_rv_pa -DUSE_CIRC_MODEL_RV_PA $(CXXFLAGS) $(LDFLAGS) $(OBJS_RV_PA) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

clean:
	$(RM) main3d 
	$(RM) main3d_rv_pa 
	$(RM) *.o 
