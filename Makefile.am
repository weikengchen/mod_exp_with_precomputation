AUTOMAKE_OPTIONS = foreign

CFLAGS = -Wall -Wno-unused-result -march=native -O3 -maes -mrdrnd -lgmp
CXXFLAGS = -Wall -Wno-unused-result -march=native -O3 -maes -mrdrnd -lgmp -std=c++11
LDFLAGS =  
LIBS = -lgmp

bin_PROGRAMS = test
test_SOURCES = main.cpp libfastexp.h libfastexp.c
