
INCLUDES=-I/usr/include
LIB=-L/usr/lib -lquadmath

CFLAGS=-O3
# CFLAGS=-g

CC=g++

TARGET=main

$(TARGET): $(TARGET).cpp math box
	$(CC) $(CFLAGS) $(INCLUDES) $(TARGET).cpp -o $(TARGET) $(LIB) quad.o math.o master_box.o
	./$(TARGET)

quad: quad.cpp
	$(CC) $(CFLAGS) $(INCLUDES) quad.cpp -c $(LIB)

box: master_box.cpp
	$(CC) $(CFLAGS) $(INCLUDES) master_box.cpp -c $(LIB)

math: math.cpp
	$(CC) $(CFLAGS) $(INCLUDES) math.cpp -c $(LIB)
