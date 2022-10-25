# Build all flutas.* binaries at once

.PHONY: all clean

.NOTPARALLEL:

ARCH ?= generic-gpu
USE_NVTX ?= 0

#APP_LIST=basic single_phase two_phase_ht two_phase_ibm two_phase_inc_isot two_phase_inc_isot_turb
APP_LIST=two_phase_ibm

all: 
	@for idapp in $(APP_LIST); \
	do \
		make -C ./src clean-obj; make -C ./src ARCH=$(ARCH) APP=$${idapp} USE_NVTX=$(USE_NVTX) flutas.$${idapp}; \
	done

clean-obj:
	make -C ./src clean-obj

clean: clean
	@for idapp in $(APP_LIST); \
	do \
		make -C ./src ARCH=$(ARCH) APP=$${idapp} clean; \
	done