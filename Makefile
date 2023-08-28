CC=g++
AR=ar

MAKE_PID := $(shell echo $$PPID)
JOBS := $(shell ps T | sed -n 's%.*$(MAKE_PID).*$(MAKE).* \(-j\|--jobs=\) *\([0-9][0-9]*\).*%\2%p')

#CPPFLAGS= -Wall -g -Lhtslib -lhts -Lcrelib -lcre optutils/opthelper.a --std=c++17
CPPFLAGS= --std=c++17 -O3 -fopenmp
LDFLAGS= -fopenmp
# CPPFLAGS= --std=c++17 -g -fopenmp -pg
# LDFLAGS= -fopenmp -lboost_iostreams -lboost_serialization -pg
HTSLIBDIR=htslib
HTSLIB=$(HTSLIBDIR)/libhts.a
SUBMODULES=crelib optutils $(HTSLIB)

PROJECT_OBJS=kled.o input.o signature.o contig.o StatsManager.o StatsTracker.o clustering.o report.o merge.o
PROJECT_HEADERS=kled.h input.h signature.h contig.h clustering.h report.h merge.h
EXAMPLE_OBJS=
LAUNCHER=./launcher

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

.PHONY: all clean $(SUBMODULES)

all: build

build: $(SUBMODULES) kled

$(HTSLIB): $(HTSLIBDIR)/config.status
	make -C $(HTSLIBDIR) -j $(JOBS)

crelib:
	make -C $@

optutils:
	make -C $@

# $(SUBMODULES):
# 	make -C $@

$(HTSLIBDIR)/config.status:
	cd $(HTSLIBDIR) && autoreconf -i && ./configure

$(PROJECT_OBJS): $(PROJECT_HEADERS)

HTSLIB_LIBS = -lz -lm -lbz2 -llzma -lcurl -lpthread -lcrypto -ldeflate
kled: $(PROJECT_OBJS) $(SUBMODULES) optutils/opthelper.a
	$(CC) $(PROJECT_OBJS) -o $@  $(LDFLAGS) $(HTSLIBDIR)/libhts.a $(HTSLIB_LIBS) -Lcrelib -lcre optutils/opthelper.a --std=c++17

clean:
	rm *.o kled

install: kled
	install -d $(PREFIX)/bin/
	install -m 755 kled $(PREFIX)/bin/