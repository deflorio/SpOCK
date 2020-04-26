SELF_DIR = $(dir $(lastword $(MAKEFILE_LIST)))
C_SRC += $(wildcard $(SELF_DIR)/*.c)
CPP_SRC_ATT += $(SELF_DIR)/SpaceEnvironment.cpp $(SELF_DIR)/MagneticField.cpp $(SELF_DIR)/Solarsys.cpp $(SELF_DIR)/Atmosphere.cpp $(SELF_DIR)/SolarRadiation.cpp
CPP_SRC_ORB += $(SELF_DIR)/SpaceEnvironment.cpp $(SELF_DIR)/Gravity.cpp $(SELF_DIR)/ThirdBody.cpp $(SELF_DIR)/Solarsys.cpp $(SELF_DIR)/Atmosphere.cpp $(SELF_DIR)/SolarRadiation.cpp
CPP_SRC_EVE += $(SELF_DIR)/Solarsys.cpp
H_SRC += $(wildcard $(SELF_DIR)/*.h)
INCLUDE_PATH += $ -I"$(SELF_DIR)"
LIBRARY_PATH += #-L"$(SELF_DIR)"
#include makefile fragments in subdirectories if they exist
-include $(SELF_DIR)/*/make.mk
