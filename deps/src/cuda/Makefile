CUDA := $(shell command -v nvcc 2> /dev/null)
CUDAFLAGS=-Wno-deprecated-gpu-targets --ptxas-options=-v -Xcompiler '-fPIC'

OBJ= rpg1.o rpg2.o rpg3.o

all: $(OBJ)

%.o : %.cu
ifndef CUDA
	$(error "Compiling GPU library functions is only possible if the CUDA toolkit is installed.")
endif
	$(CUDA) $(CUDAFLAGS) -c $^

clean:
	rm -f $(OBJ)
