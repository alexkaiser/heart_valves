######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
IBAMR_SRC_DIR   = $(HOME)/sfw/ibamr/IBAMR
IBAMR_BUILD_DIR = $(HOME)/sfw/ibamr/ibamr-objs-dbg
# IBAMR_BUILD_DIR = $(HOME)/sfw/ibamr/ibamr-objs-dbg_prescribed_motion

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

######################################################################
## Build the application.

all: main3d main_rv_pa main_aorta

PDIM = 3

# common 
OBJS = CirculationModel_with_lv.o \
       CirculationModel_RV_PA.o \
       CirculationModel_aorta.o \
       CirculationModel.o \
       boundary_condition_util.o \
       FourierBodyForce.o \
       pnpoly.o \
       FeedbackForcer.o

MAIN = main.o

main3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) $(MAIN)
	$(CXX) -o main3d $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(MAIN) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

MAIN_RV_PA = main_rv_pa.o

main_rv_pa: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) $(MAIN_RV_PA)
	$(CXX) -o main_rv_pa $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(MAIN_RV_PA) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

MAIN_AORTA = main_aorta.o

main_aorta: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) $(MAIN_AORTA)
	$(CXX) -o main_aorta $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(MAIN_AORTA) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) 

clean:
	$(RM) main3d        
	$(RM) main_rv_pa 
	$(RM) main_aorta 
	$(RM) *.o 

