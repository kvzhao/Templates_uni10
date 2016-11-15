PROJECT := uni10
VERSION :=
VNUM    := 0
CONFIG_FILE := Makefile.config

####### Check the version of Uni10 #######

ifeq ($(USE_ARMADILLO), 1)
	VERSION = armadillo_cpu
	VNUM = 2
else ifeq ($(USE_MAGMA), 1)	
	VERSION = magma_gpu
	VNUM = 3
else ifeq ($(USE_CUSOLVER), 1)	
	VERSION = cuda_gpu
	VNUM = 4
else
	VERSION = lapack_cpu
	VNUM = 1
endif

COMMON_FLAGS = -std=c++11

##########################################

# Explicitly check for the config file, otherwise make -k will proceed anyway.
ifeq ($(wildcard $(CONFIG_FILE)),)
$(error $(CONFIG_FILE) not found.)
endif
include $(CONFIG_FILE)

BUILD_DIR_LINK := $(BUILD_DIR)
ifeq ($(RELEASE_BUILD_DIR),)
	RELEASE_BUILD_DIR := .$(BUILD_DIR)_release
endif
ifeq ($(DEBUG_BUILD_DIR),)
	DEBUG_BUILD_DIR := .$(BUILD_DIR)_debug
endif

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	BUILD_DIR := $(DEBUG_BUILD_DIR)
	OTHER_BUILD_DIR := $(RELEASE_BUILD_DIR)
else
	BUILD_DIR := $(RELEASE_BUILD_DIR)
	OTHER_BUILD_DIR := $(DEBUG_BUILD_DIR)
endif


SRC_DIRS := $(shell find * -type d -exec bash -c "find {} -maxdepth 1 \
	\( -name '*.cpp' \) | grep -q ." \; -print)

INC_DIRS := $(shell find * -type d -exec bash -c "find {} -maxdepth 1 \
	\( -name '*.h' -o -name '*.hpp' \) | grep -q ." \; -print)

# The target shared library name
LIB_BUILD_DIR := $(BUILD_DIR)/lib
INC_BUILD_DIR := $(BUILD_DIR)/include
STATIC_NAME := $(LIB_BUILD_DIR)/lib$(PROJECT).a
DYNAMIC_NAME := $(LIB_BUILD_DIR)/lib$(PROJECT).so

##############################
# Get all source files
##############################

# CXX_SRCS are the source files excluding the test ones.
CXX_SRCS :=
ifeq ($(VNUM), 1)
	#Find the auxiliary files
	CXX_SRCS += $(shell find src/$(PROJECT) -maxdepth 1 ! -name "*_gpu.cpp" ! -name "test_*.cpp" -name "*.cpp")
	#Find UNI10_API FILES
	CXX_SRCS += $(shell find src/$(PROJECT)/uni10_api -name "*.cpp")
	#Find the others sources of the specified version.
	CXX_SRCS += $(shell find src/$(PROJECT)/uni10_$(VERSION) -name "*.cpp")
endif

# EXAMPLE_SRCS are the source files for the example binaries
EXAMPLE_SRCS := $(shell find example -name "*.cpp")

# BUILD_INCLUDE_DIR contains any generated header files we want to include.
BUILD_INCLUDE_DIR := $(BUILD_DIR)/include

##############################
# Get gtest source files
##############################

GTEST_SRC := gtest

# handle IO dependencies
USE_OPENCV ?= 1

ifeq ($(USE_OPENCV), 1)
	LIBRARIES += opencv_core opencv_highgui opencv_imgproc 

	ifeq ($(OPENCV_VERSION), 3)
		LIBRARIES += opencv_imgcodecs
	endif
		
endif

##############################
# Derive generated files
##############################

CXX_OBJS := $(addprefix $(BUILD_DIR)/, ${CXX_SRCS:.cpp=.o})

EXAMPLE_OBJS := $(addprefix $(BUILD_DIR)/, ${EXAMPLE_SRCS:.cpp=.o})

EXAMPLE_BINS := ${EXAMPLE_OBJS:.o=.bin}

##############################
# Derive compiler warning dump locations
##############################

WARNS_EXT := warnings.txt

CXX_WARNS := $(addprefix $(BUILD_DIR)/, ${CXX_SRCS:.cpp=.o.$(WARNS_EXT)})

EXAMPLE_WARNS := $(addprefix $(BUILD_DIR)/, ${EXAMPLE_SRCS:.cpp=.o.$(WARNS_EXT)})

##############################
# Derive include and lib directories
##############################

INCLUDE_DIRS += $(BUILD_INCLUDE_DIR) ./src ./include
# CPU-only configuration
ifeq ($(VNUM), 1)
	OBJS := $(CXX_OBJS)
	ALL_WARNS := $(ALL_CXX_WARNS)
	COMMON_FLAGS += -DCPU -DLAPACK
  ifeq ($(USE_BOOST), 1)
		COMMON_FLAGS += -DBOOST
	endif
else ifeq ($(VNUM), 2)
	OBJS := $(CXX_OBJS)
	COMMON_FLAGS += -DCPU -DARMADILLO
else ifeq ($(VNUM), 3)
	OBJS := $(CXX_OBJS)
	COMMON_FLAGS += -DGPU -DMAGMA
else ifeq ($(VNUM), 4)
	OBJS := $(CXX_OBJS)
	COMMON_FLAGS += -DCPU -DCUDAONLY
endif


LIBRARIES += boost_system boost_filesystem m hdf5_hl hdf5
# glog gflags protobuf

##############################
# Set build directories
##############################

ALL_BUILD_DIRS := $(sort $(BUILD_DIR) $(LIB_BUILD_DIR) \
	$(addprefix $(BUILD_DIR)/, $(SRC_DIRS)))

##############################
# Configure build
##############################

# Determine platform
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	LINUX := 1
else ifeq ($(UNAME), Darwin)
	OSX := 1
endif

# Linux
ifeq ($(LINUX), 1)
	CXX ?= /usr/bin/g++
	GCCVERSION := $(shell $(CXX) -dumpversion | cut -f1,2 -d.)
	# older versions of gcc are too dumb to build boost with -Wuninitalized
	ifeq ($(shell echo | awk '{exit $(GCCVERSION) < 4.6;}'), 1)
		WARNINGS += -Wno-uninitialized
	endif
	# boost::thread is reasonably called boost_thread (compare OS X)
	# We will also explicitly add stdc++ to the link target.
	LIBRARIES += boost_thread stdc++
endif

# OS X:
# clang++ instead of g++
# libstdc++ for NVCC compatibility on OS X >= 10.9 with CUDA < 7.0
ifeq ($(OSX), 1)
	CXX := /usr/bin/clang++
	ifneq ($(VNUM), 1)
		CUDA_VERSION := $(shell $(CUDA_DIR)/bin/nvcc -V | grep -o 'release \d' | grep -o '\d')
		ifeq ($(shell echo | awk '{exit $(CUDA_VERSION) < 7.0;}'), 1)
			CXXFLAGS += -stdlib=libstdc++
			LINKFLAGS += -stdlib=libstdc++
		endif
		# clang throws this warning for cuda headers
		WARNINGS += -Wno-unneeded-internal-declaration
	endif
	# gtest needs to use its own tuple to not conflict with clang
	COMMON_FLAGS += -DGTEST_USE_OWN_TR1_TUPLE=1
	# boost::thread is called boost_thread-mt to mark multithreading on OS X
	LIBRARIES += boost_thread-mt
	# we need to explicitly ask for the rpath to be obeyed
	DYNAMIC_FLAGS := -install_name @rpath/libuni10.so
	ORIGIN := @loader_path
else
	ORIGIN := \$$ORIGIN
endif

# Custom compiler
ifdef CUSTOM_CXX
	CXX := $(CUSTOM_CXX)
endif

# Static linking
ifneq (,$(findstring clang++,$(CXX)))
	STATIC_LINK_COMMAND := -Wl,-force_load $(STATIC_NAME)
else ifneq (,$(findstring g++,$(CXX)))
	STATIC_LINK_COMMAND := -Wl,--whole-archive $(STATIC_NAME) -Wl,--no-whole-archive
else
  # The following line must not be indented with a tab, since we are not inside a target
  $(error Cannot static link with the $(CXX) compiler)
endif

#################################
# Setting object directories
#################################

# Debugging
ifeq ($(DEBUG), 1)
	COMMON_FLAGS += -DDEBUG -g -O0
	NVCCFLAGS += -G
else
	COMMON_FLAGS += -DNDEBUG -O2
endif

# BLAS configuration (default = ATLAS)
BLAS ?= atlas
ifeq ($(BLAS), mkl)
	# MKL
	LIBRARIES += mkl_rt
	COMMON_FLAGS += -DMKL
	MKL_DIR ?= /opt/intel/mkl
	BLAS_INCLUDE ?= $(MKL_DIR)/include
	BLAS_LIB ?= $(MKL_DIR)/lib $(MKL_DIR)/lib/intel64
else ifeq ($(BLAS), open)
	# OpenBLAS
	LIBRARIES += openblas
else
	# ATLAS
	ifeq ($(LINUX), 1)
		ifeq ($(BLAS), atlas)
			# Linux simply has cblas and atlas
			#LIBRARIES += cblas atlas
			LIBRARIES += satlas tatlas
		endif
	else ifeq ($(OSX), 1)
		# OS X packages atlas as the vecLib framework
		LIBRARIES += cblas
		# 10.10 has accelerate while 10.9 has veclib
		XCODE_CLT_VER := $(shell pkgutil --pkg-info=com.apple.pkg.CLTools_Executables | grep 'version' | sed 's/[^0-9]*\([0-9]\).*/\1/')
		XCODE_CLT_GEQ_6 := $(shell [ $(XCODE_CLT_VER) -gt 5 ] && echo 1)
		ifeq ($(XCODE_CLT_GEQ_6), 1)
			BLAS_INCLUDE ?= /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
			LDFLAGS += -framework Accelerate
		else
			BLAS_INCLUDE ?= /System/Library/Frameworks/vecLib.framework/Versions/Current/Headers/
			LDFLAGS += -framework vecLib
		endif
	endif
endif

INCLUDE_DIRS += $(BLAS_INCLUDE)
LIBRARY_DIRS += $(BLAS_LIB)

LIBRARY_DIRS += $(LIB_BUILD_DIR)

# Automatic dependency generation (nvcc is handled separately)
#CXXFLAGS += -MMD -MP

# Complete build flags.
COMMON_FLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CXXFLAGS += -pthread -fPIC $(COMMON_FLAGS) $(WARNINGS)

# mex may invoke an older gcc that is too liberal with -Wuninitalized
LINKFLAGS += -pthread -fPIC $(COMMON_FLAGS) $(WARNINGS)

USE_PKG_CONFIG ?= 0

LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(PKG_CONFIG) \
		$(foreach library,$(LIBRARIES),-l$(library))

# 'superclean' target recursively* deletes all files ending with an extension
# in $(SUPERCLEAN_EXTS) below.  This may be useful if you've built older
# versions of Caffe that do not place all generated files in a location known
# to the 'clean' target.
#
# 'supercleanlist' will list the files to be deleted by make superclean.
#
# * Recursive with the exception that symbolic links are never followed, per the
# default behavior of 'find'.
SUPERCLEAN_EXTS := .so .a .o .bin .testbin .pb.cc .pb.h _pb2.py .cuo

# Set the sub-targets of the 'everything' target.
EVERYTHING_TARGETS := all py$(PROJECT) test warn lint

##############################
# Define build targets
##############################

all: lib examples

lib: $(STATIC_NAME) $(DYNAMIC_NAME)

examples: $(EXAMPLE_BINS)

tools: $(TOOL_BINS) $(TOOL_BIN_LINKS)

$(BUILD_DIR_LINK): $(BUILD_DIR)/.linked

# Create a target ".linked" in this BUILD_DIR to tell Make that the "build" link
# is currently correct, then delete the one in the OTHER_BUILD_DIR in case it
# exists and $(DEBUG) is toggled later.
$(BUILD_DIR)/.linked:
	@ mkdir -p $(BUILD_DIR)
	@ $(RM) $(OTHER_BUILD_DIR)/.linked
	@ $(RM) -r $(BUILD_DIR_LINK)
	@ cp -r include $(BUILD_DIR)/include
	@ ln -s $(BUILD_DIR) $(BUILD_DIR_LINK)
	@ touch $@

$(ALL_BUILD_DIRS): | $(BUILD_DIR_LINK)
	@ mkdir -p $@

$(DYNAMIC_NAME): $(OBJS) | $(LIB_BUILD_DIR)
	@ echo LD -o $@
	$(Q)$(CXX) -shared -o $@ $(OBJS) $(LINKFLAGS) $(LDFLAGS) $(DYNAMIC_FLAGS)

$(STATIC_NAME): $(OBJS) | $(LIB_BUILD_DIR)
	@ echo AR -o $@
	$(Q)ar rcs $@ $(OBJS)

$(BUILD_DIR)/%.o: %.cpp | $(ALL_BUILD_DIRS)
	@ echo CXX $<
	$(Q)$(CXX) $< $(CXXFLAGS) -c -o $@ 2> $@.$(WARNS_EXT) \
		|| (cat $@.$(WARNS_EXT); exit 1)
	@ cat $@.$(WARNS_EXT)

$(EXAMPLE_BINS): %.bin : %.o | $(DYNAMIC_NAME)
	@ echo CXX/LD -o $@
	$(Q)$(CXX) $< -o $@ $(LINKFLAGS) -l$(PROJECT) $(LDFLAGS) \
		-Wl,-rpath,$(ORIGIN)/../../lib

examples: $(EXAMPLE_BINS)

clean:
	@- $(RM) -rf $(ALL_BUILD_DIRS)
	@- $(RM) -rf $(OTHER_BUILD_DIR)
	@- $(RM) -rf $(BUILD_DIR_LINK)

print-%  : 
	@echo $* = $($*)

