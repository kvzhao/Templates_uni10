TUNI10_ROOT := /Users/Yun-Hsuan/GitRepo/Templates_uni10/uni10/build
#UNI10_ROOT := /Users/yunda/myLibrary/uni10
INC := $(TUNI10_ROOT)/include
LIB := $(TUNI10_ROOT)/lib
CC := icpc
FLAGS:=-O3 -m64 -std=c++11 -DCPU -DLAPACK
OBJ := $(LIB)/libuni10.a

add.e: add.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

dot.e: matrixdot.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

trans.e: trans.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

qr.e: qr_family.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

resize.e: resize.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

svd.e: svd.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

copy_constructor.e: copy_constructor.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

main_cpu.e: main_cpu.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ)

leak_check:
	valgrind --leak-check=yes ./dot.e

clean:
	rm -f *.o *.e asd
