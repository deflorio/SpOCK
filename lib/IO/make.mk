SELF_DIR = $(dir $(lastword $(MAKEFILE_LIST)))
C_SRC += $(wildcard $(SELF_DIR)/*.c)
CPP_SRC_ATT += $(wildcard $(SELF_DIR)/*.cpp)
CPP_SRC_ORB += $(wildcard $(SELF_DIR)/*.cpp)
CPP_SRC_EVE += $(wildcard $(SELF_DIR)/*.cpp)
CPP_SRC_NAVSOLCONV += $(wildcard $(SELF_DIR)/*.cpp)
H_SRC += $(wildcard $(SELF_DIR)/*.h)
INCLUDE_PATH += $ -I"$(SELF_DIR)"
LIBRARY_PATH += #-L"$(SELF_DIR)"
#include makefile fragments in subdirectories if they exist
-include $(SELF_DIR)/*/make.mk
