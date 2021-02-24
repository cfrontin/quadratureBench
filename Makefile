
INCLUDES=-I/usr/include
LIB=-L/usr/lib -lquadmath

CFLAGS=-O3
# CFLAGS=-g
# CFLAGS=-g -v -ftime-report

CC=g++

TARGET=main

$(TARGET): $(TARGET).cpp math box quad
	# $(CC) $(CFLAGS) $(INCLUDES) $(TARGET).cpp -o $(TARGET) $(LIB) math.o master_box.o
	$(CC) $(CFLAGS) $(INCLUDES) $(TARGET).cpp -o $(TARGET) $(LIB) quad.o math.o master_box.o
	./$(TARGET)

build: $(TARGET).cpp math box quad
	# $(CC) $(CFLAGS) $(INCLUDES) $(TARGET).cpp -o $(TARGET) $(LIB) math.o master_box.o
	$(CC) $(CFLAGS) $(INCLUDES) $(TARGET).cpp -o $(TARGET) $(LIB) quad.o math.o master_box.o

quad: quad.cpp
	$(CC) $(CFLAGS) $(INCLUDES) quad.cpp -c $(LIB)

box: master_box.cpp
	$(CC) $(CFLAGS) $(INCLUDES) master_box.cpp -c $(LIB)

math: math.cpp
	$(CC) $(CFLAGS) $(INCLUDES) math.cpp -c $(LIB)
