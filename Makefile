###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR     - Root directory of the MOOSE project
# FRAMEWORK_DIR - Location of the MOOSE framework
#
###############################################################################
EXAMPLE_DIR        ?= $(shell dirname `pwd`)
MOOSE_DIR          ?= $(shell dirname $(EXAMPLE_DIR))
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
ADDITIONAL_CPPFLAGS += -Wall -Wextra
###############################################################################

TEST := test_ignore

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

APPLICATION_NAME := fieldmapping_test

# dep apps
APPLICATION_DIR    := $(shell pwd)
APPLICATION_NAME   := fieldmapping_test
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk
include            ../test.mk

# Include dependency files for this example
ex_srcfiles := $(shell find $(APPLICATION_DIR) -name "*.C")
ex_deps     := $(patsubst %.C, %.$(obj-suffix).d, $(ex_srcfiles))
-include $(ex_deps)
