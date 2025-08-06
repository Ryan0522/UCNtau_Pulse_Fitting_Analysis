CXX = g++
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

NLOPT_LIBS = -lnlopt

INCLUDE_JSON = "/mnt/c/Users/ciyou/Desktop/UIUC/Spring 2024/Research/Prof Liu/UCNtau Counting/Ryan_Analysis"

Pulse_Analysis: src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp include/Pulse_Analysis.h include/Pulse_Fitting.h
	$(CXX) -o $@ src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp \
		-Iinclude -I$(INCLUDE_JSON) $(NLOPT_CFLAGS) -g \
		$(ROOT_CFLAGS) $(ROOT_LIBS) \
		$(NLOPT_LIBS)

Runtime_Analysis_: src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp include/Pulse_Analysis.h include/Pulse_Fitting.h
	$(CXX) -o $@ src/Pulse_Analysis.cpp src/Pulse_Fitting.cpp \
		-Iinclude -I$(INCLUDE_JSON) \
		-pg -O2 -g $(ROOT_CFLAGS) $(ROOT_LIBS) $(NLOPT_LIBS)

clean:
	rm -f Pulse_Analysis

.DEFAULT_GOAL := Pulse_Analysis
.PHONY: clean
