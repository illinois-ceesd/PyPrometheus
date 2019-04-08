include $(WORK)/cantera-2.2.1/cantera-2.2.1/include/cantera/Cantera.mak

SRCDIR = $(PROJ)/Prometheus/src
INCDIR = $(PROJ)/Prometheus/include
TSTDIR = $(PROJ)/Prometheus/test
APPDIR = ${PROJ}/Prometheus/apps
RUNDIR = $(PROJ)/Prometheus/run

CXX   = g++ -std=c++11
FLAG  = -g -O0
LIBS  = $(CANTERA_LIBS)
INCS  = -I$(INCDIR)
INCS += $(CANTERA_INCLUDES)

all: test exec

test: testSanDiego testSanDiegoFull testGRI30

testSanDiego: sanDiego.o testSanDiego.o
	$(CXX) -o testSanDiego.exe sanDiego.o testSanDiego.o $(LIBS)

testSanDiegoFull: sanDiegoFull.o testSanDiegoFull.o
	$(CXX) -o testSanDiegoFull.exe sanDiegoFull.o testSanDiegoFull.o $(LIBS)

testGRI30: gri30.o testGRI30.o
	$(CXX) -o testGRI30.exe gri30.o testGRI30.o $(LIBS)

exec: utils.o opts.o Prometheus.o driver.o
	$(CXX) -o Prometheus.exe utils.o opts.o Prometheus.o driver.o $(LIBS)

utils.o: $(SRCDIR)/utils.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

opts.o: $(SRCDIR)/opts.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

Prometheus.o: $(SRCDIR)/Prometheus.cpp utils.o
	$(CXX) -c $(FLAG) $< $(INCS)

sanDiego.o: $(SRCDIR)/sanDiego.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

sanDiegoFull.o: $(SRCDIR)/sanDiegoFull.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

gri30.o: $(SRCDIR)/gri30.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

testSanDiego.o: $(TSTDIR)/testSanDiego.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

testSanDiegoFull.o: $(TSTDIR)/testSanDiegoFull.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

testGRI30.o: $(TSTDIR)/testGRI30.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

driver.o: $(SRCDIR)/driver.cpp Prometheus.o
	$(CXX) -c $(FLAG) $< $(INCS)

clean:
	rm -f *.o *.exe

realclean: clean
	rm -f $(SRCDIR)/*~ $(INCDIR)/*~ $(TSTDIR)/*~ *~
