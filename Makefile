CC=g++
AR=ar
CPPFLAGS= -Wall -g -Lhtslib -lhts -Lcrelib -lcre optutils/opthelper.a --std=c++17
LDFLAGS=-lz -lm -lbz2 -llzma -lpthread
LIBS=
PYTHON=python3.8
INCLUDE=/usr/include/$(PYTHON)

PROJECT_OBJS=kled.o input.o variant.o signature.o contig.o StatsManager.o StatsTracker.o
PROJECT_HEADERS=
EXAMPLE_OBJS=
HTSLIB=htslib/libhts.a
LAUNCHER=./launcher

.PHONY: all test clean

all: build

build: crelib htslib optutils kled

crelib:
	cd crelib && make

htslib:
	cd htslib && make

optutils:
	cd optutils && make

kled: $(PROJECT_OBJS)
	$(CC) $^ -o $@ -g -pg -fopenmp -Lhtslib -lhts -lz -Lcrelib -lcre optutils/opthelper.a --std=c++17

clean:
	rm *.o kled