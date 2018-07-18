CC=gcc
CXX=g++
CUDA := $(shell command -v nvcc 2> /dev/null)
CLIBS = -lgomp
CUDALIBS=-L /usr/local/cuda-8.0/targets/x86_64-linux/lib -lcudart -lcuda
LINK_OPTS=-Wl,-z,defs

SDIR=src/c
CDIR=src/cuda

_COBJ   = dpp1.o dpp2.o dpp3.o rpp1.o rpp2.o rpp3.o rpp2_multi.o
_CUDAOBJ= rpg1.o rpg2.o rpg3.o
C_OBJ   =$(patsubst %,$(SDIR)/%,$(_COBJ))
CUDA_OBJ=$(patsubst %,$(CDIR)/%,$(_CUDAOBJ))

ifdef CUDA
OBJ=$(C_OBJ) $(CUDA_OBJ)
LIBS=$(CLIBS) $(CUDALIBS)
else
OBJ = $(C_OBJ)
LIBS=$(CLIBS)
endif

SHARED_LIB_NAME = usr/lib/libnxs.so

all: gpu cpu $(SHARED_LIB_NAME)

$(SHARED_LIB_NAME): $(OBJ)
	@mkdir -p $(@D)
	$(CC) -shared -o $@ $^ $(LIBS) $(LINK_OPTS)

gpu: force_look
ifdef CUDA
	cd src/cuda; make
endif

cpu: force_look
	cd src/c; make

force_look:
	true

clean:
	rm -f $(C_OBJ) $(CUDA_OBJ) $(SHARED_LIB_NAME)
