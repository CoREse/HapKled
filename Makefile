CC=g++
AR=ar
#CPPFLAGS= -Wall -g -Lhtslib -lhts -Lcrelib -lcre optutils/opthelper.a --std=c++17
CPPFLAGS= -g --std=c++17 -fopenmp -O3
LDFLAGS=-lz -lm -lbz2 -llzma -lpthread
LIBS=
PYTHON=python3.8
INCLUDE=/usr/include/$(PYTHON)

PROJECT_OBJS=kled.o input.o signature.o contig.o StatsManager.o StatsTracker.o clustering.o report.o
PROJECT_HEADERS=kled.h input.h signature.h contig.h clustering.h report.h defines.h
EXAMPLE_OBJS=
HTSLIB=htslib/libhts.a
LAUNCHER=./launcher

HTSPREFIX=htslib/

.PHONY: all test clean

all: build

build: crelib htslib optutils kled

crelib:
	cd crelib && make

htslib:
	cd htslib && make

optutils:
	cd optutils && make

$(PROJECT_OBJS): $(PROJECT_HEADERS)

HTSLIB_LIBS = -lz -lm -lbz2 -llzma -lcurl -lpthread -lcrypto -ldeflate
kled: $(PROJECT_OBJS)
	$(CC) $^ -o $@ -g -pg -fopenmp $(HTSPREFIX)/libhts.a $(HTSLIB_LIBS) -lz -Lcrelib -lcre optutils/opthelper.a --std=c++17

clean:
	rm *.o kled
