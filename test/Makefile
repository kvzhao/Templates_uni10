INC := $(TUNI10_ROOT)/include
LIB := $(TUNI10_ROOT)/lib
CC := icpc
FLAGS:= -m64 -std=c++11 -DCPU -DLAPACK
OBJ := $(LIB)/libuni10.a

net.e: network.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

tract.e: contract.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

mute.e: permute.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

ten.e: tensor.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

add.e: add.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

iden.e: identity.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

dots.e: dots.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

eigh.e: eigh.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

qdr.e: qdr.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

copy_constructor.e: copy_constructor.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

qr.e: qr_family.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

dot.e: matrixdot.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

arg2.e: arg2.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

svd.e: svd.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

rng.e: rng.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

main_cpu.e: main_cpu.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

trans.e: trans.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

resize.e: resize.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -mkl=parallel

leak_check:
	valgrind --leak-check=yes ./tract.e

clean:
	rm -f *.o *.e asd
