SELF_DIR = $(dir $(lastword $(MAKEFILE_LIST)))
C_SRC += $(wildcard $(SELF_DIR)/*.c)
CPP_SRC_ATT += $(SELF_DIR)/Propagator.cpp $(SELF_DIR)/Attitude.cpp $(SELF_DIR)/AttitudeQ.cpp
CPP_SRC_ORB += $(SELF_DIR)/Propagator.cpp $(SELF_DIR)/Orbit.cpp
CPP_SRC_EVE += $(SELF_DIR)/AnalyticalModels.cpp
H_SRC += $(wildcard $(SELF_DIR)/*.h)
INCLUDE_PATH += $ -I"$(SELF_DIR)"
LIBRARY_PATH += #-L"$(SELF_DIR)"
#include makefile fragments in subdirectories if they exist
-include $(SELF_DIR)/*/make.mk
