PROJECT(mysolver)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

# List your source files here 
SET(SOURCE_FILES main_aorta.cpp 
       CirculationModel_aorta.cpp
       CirculationModel_with_lv.cpp
       CirculationModel_RV_PA.cpp
       CirculationModel.cpp
       boundary_condition_util.cpp
       FourierBodyForce.cpp
       pnpoly.cpp
       FeedbackForcer.cpp)

# ADD_EXECUTABLE(main_aorta ${SOURCE_FILES})

SET(SOURCE_FILES main_aorta_fe.cpp 
       CirculationModel_aorta.cpp
       CirculationModel_with_lv.cpp
       CirculationModel_RV_PA.cpp
       CirculationModel.cpp
       boundary_condition_util.cpp
       FourierBodyForce.cpp
       pnpoly.cpp
       FeedbackForcer.cpp)

ADD_EXECUTABLE(main_aorta_fe ${SOURCE_FILES})

SET(SOURCE_FILES_PERIODIC main.cpp 
       CirculationModel_aorta.cpp
       CirculationModel_with_lv.cpp
       CirculationModel_RV_PA.cpp
       CirculationModel.cpp
       boundary_condition_util.cpp
       FourierBodyForce.cpp
       pnpoly.cpp
       FeedbackForcer.cpp)

# ADD_EXECUTABLE(main ${SOURCE_FILES_PERIODIC})


SET(SOURCE_FILES_RV_PA main_rv_pa.cpp 
       CirculationModel_aorta.cpp
       CirculationModel_with_lv.cpp
       CirculationModel_RV_PA.cpp
       CirculationModel.cpp
       boundary_condition_util.cpp
       FourierBodyForce.cpp
       pnpoly.cpp
       FeedbackForcer.cpp)

# ADD_EXECUTABLE(main_rv_pa ${SOURCE_FILES_RV_PA})


SET(SOURCE_FILES_RV_PA main_rv_pa_fe.cpp 
       CirculationModel_aorta.cpp
       CirculationModel_with_lv.cpp
       CirculationModel_RV_PA.cpp
       CirculationModel.cpp
       boundary_condition_util.cpp
       FourierBodyForce.cpp
       pnpoly.cpp
       FeedbackForcer.cpp)

# ADD_EXECUTABLE(main_rv_pa_fe ${SOURCE_FILES_RV_PA})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(IBAMR 0.16.0 REQUIRED)
# TARGET_LINK_LIBRARIES(main_aorta IBAMR::IBAMR3d)
TARGET_LINK_LIBRARIES(main_aorta_fe IBAMR::IBAMR3d)
# TARGET_LINK_LIBRARIES(main IBAMR::IBAMR3d)
# TARGET_LINK_LIBRARIES(main_rv_pa IBAMR::IBAMR3d)
# TARGET_LINK_LIBRARIES(main_rv_pa_fe IBAMR::IBAMR3d)
# IBAMR saves the flags it used to compile - you can reuse them if you want to
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})




