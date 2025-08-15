CXX = g++

ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)
NLOPT_LIBS = -lnlopt

INCLUDE_JSON = "/projects/illinois/eng/physics/chenyliu/Ryan_ciyouh2/UCNtau_Pulse_Fitting_Analysis/json"

CXXFLAGS = -Iinclude -I$(INCLUDE_JSON) -g
LDFLAGS = $(ROOT_CFLAGS) $(ROOT_LIBS) $(NLOPT_LIBS)

PROF_CXXFLAGS = $(CXXFLAGS) -pg -O2 -fno-pie
PROF_LDFLAGS  = $(LDFLAGS)  -pg -no-pie

.DEFAULT_GOAL := Pulse_Analysis

Pulse_Analysis:  src/File_Loader.cpp src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp \
			include/File_Loader.h include/Pulse_Analysis.h include/Pulse_Fitting.h
	$(CXX) -o $@ src/File_Loader.cpp src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp $(CXXFLAGS) $(LDFLAGS)

Runtime_Analysis: src/File_Loader.cpp src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp \
			include/File_Loader.h include/Pulse_Analysis.h include/Pulse_Fitting.h
	$(CXX) -o $@ src/File_Loader.cpp src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp $(PROF_CXXFLAGS) $(PROF_LDFLAGS)

Pulse_Tail: src/File_Loader.cpp src/Pulse_Tail.cpp src/Pulse_Fitting.cpp \
            include/File_Loader.h include/Pulse_Tail.h include/Pulse_Fitting.h
	$(CXX) -o $@ src/File_Loader.cpp src/Pulse_Tail.cpp src/Pulse_Fitting.cpp $(CXXFLAGS) $(LDFLAGS)

Plot_PDFs: src/File_Loader.cpp src/PDF_Lookup.cpp src/Plot_PDFs.cpp \
           include/File_Loader.h include/PDF_Lookup.h
	$(CXX) $(CXXFLAGS) src/File_Loader.cpp src/PDF_Lookup.cpp src/Plot_PDFs.cpp -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f Pulse_Analysis Runtime_Analysis_ Pulse_Tail

.PHONY: clean
