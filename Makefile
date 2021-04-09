# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -std=c++14 -Wall -g
# CFLAGS = -std=c++14 -Wall -O3 -ffast-math -march=native -fwhole-program
 
# ****************************************************
# Targets needed to bring the executable up to date
 
main: main.o
	$(CC) $(CFLAGS) -o main main.o -lGL -lGLEW -lSDL2 -lGLU
 
# The main.o target can be written more simply
 
main.o: main.cpp vec3.h
	$(CC) $(CFLAGS) -c main.cpp

clean:
	rm -f main.o main
