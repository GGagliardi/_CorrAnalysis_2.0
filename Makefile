IDIR=include
IPY=/usr/include/python3.8 
GCC=g++
LIBS= -lgsl -lgslcblas -lboost_system -lboost_filesystem  -lMinuit2 -fopenmp -lpython3.8 -lmpfr -lgmp
LIB_PY=/usr/lib/python3.8/config-3.8m-x86_64-linux-gnu
CPPFLAGS= -I$(IDIR) -I$(IPY) -Wall -fopenmp -std=c++17 -O3 -g -rdynamic
SOURCES=$(wildcard src/*.cpp)
OBJ_1=$(patsubst %.cpp,%.o,$(SOURCES))
OBJ=$(patsubst src/%, obj/%, $(OBJ_1))
DEP1= $(patsubst %.cpp,%.d,$(SOURCES))
DEP=$(patsubst src/%, dep/%, $(DEP1))

dep/%.d: src/%.cpp
	@set -e; rm -f $@; \
	$(GCC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$


obj/%.o: src/%.cpp
	$(GCC) -c  $(CPPFLAGS) -o $@ $<

#-include $(DEP)


stat_analysis: $(OBJ)
	$(GCC)  -o $@ $^ $(CPPFLAGS) -L $(LIB_PY) $(LIBS)



.PHONY: clean

clean :
	rm -f stat_analysis $(OBJ) $(DEP)
