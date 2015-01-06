CXX = g++
CXXFLAGS = -c -g -Wall
LIBS = -lboost_math_tr1
PROGS = ellint_ex

.PHONY: all
all: $(PROGS)

ellint_ex: ellint_ex.o 
	$(CXX) -o $@ $^ $(LIBS) 

.PHONY: clean
clean:
	$(RM) -r *.o
	$(RM) .depend

%.o: %.cc
	$(CXX) $(CXXFLAGS) $^ -o $@ 

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

