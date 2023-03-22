# Build all flutas.* binaries at once

.PHONY: all clean

.NOTPARALLEL:

ARCH ?= nvf_gpu
USE_NVTX ?= 0

APP_LIST=two_phase_inc_ibm #basic single_phase_inc_ibm two_phase_ht_inc_ibm two_phase_inc_isot two_phase_inc_isot_turb

all: 
	@for idapp in $(APP_LIST); \
	do \
		make -C ./src clean-obj; make -C ./src ARCH=$(ARCH) APP=$${idapp} USE_NVTX=$(USE_NVTX) flutas.$${idapp}; \
	done

post: 
	make -C ./src clean-obj; make -C ./src ARCH=$(ARCH) USE_NVTX=$(USE_NVTX) DO_POSTPROC=1 flutas.post;

clean-obj:
	make -C ./src clean-obj

clean: clean
	@for idapp in $(APP_LIST); \
	do \
		make -C ./src ARCH=$(ARCH) APP=$${idapp} clean; \
	done