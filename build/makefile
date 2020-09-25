include $(WORK)/cantera-2.4.1/cantera-2.4.1/include/cantera/Cantera.mak

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

all: exec

exec: Utils.o Config.o RunPrometheus.o
	$(CXX) -o Prometheus.exe Utils.o Config.o CppPrometheus.o PyPrometheus.o RunPrometheus.o $(LIBS)

Utils.o: $(SRCDIR)/Utils.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

Config.o: $(SRCDIR)/Config.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

CppPrometheus.o: $(SRCDIR)/CppPrometheus.cpp Utils.o Config.o
	$(CXX) -c $(FLAG) $< $(INCS)

PyPrometheus.o: $(SRCDIR)/PyPrometheus.cpp Utils.o Config.o
	$(CXX) -c $(FLAG) $< $(INCS)

RunPrometheus.o: $(SRCDIR)/RunPrometheus.cpp CppPrometheus.o PyPrometheus.o
	$(CXX) -c $(FLAG) $< $(INCS)

clean:
	rm -f *.o *.exe

realclean: clean
	rm -f $(SRCDIR)/*~ $(INCDIR)/*~ $(TSTDIR)/*~ *~
