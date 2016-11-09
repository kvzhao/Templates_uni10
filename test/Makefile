TUNI10_ROOT := /home/Yun-Hsuan/GitRepo/Templates_uni10/uni10/build
#UNI10_ROOT := /Users/yunda/myLibrary/uni10
INC := $(TUNI10_ROOT)/include
LIB := $(TUNI10_ROOT)/lib
CC := g++
FLAGS:=-O3 -m64 -std=c++11 -DCPU
OBJ := $(LIB)/libuni10.a

main_cpu.e: main_cpu.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ)

clean:
	rm -f *.o *.e asd
