INC := $(TUNI10_ROOT)/include
LIB := $(TUNI10_ROOT)/lib
CC := icpc
FLAGS:= -m64 -std=c++11 -DCPU -DLAPACK
OBJ := $(LIB)/libuni10.a

tebd_1d.e: tebd_1d.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

constructor.e: constructor.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

contract.e: contract.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

eigh.e: eigh.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

exph.e: exph.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

lq.e: lq.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

qr.e: qr.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

svd.e: svd.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

dots.e: dots.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

add.e: add.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

leak_check:
	valgrind --leak-check=yes ./constructor.e

clean:
	rm -f *.o *.e asd
