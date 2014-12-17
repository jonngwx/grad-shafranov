CXX = g++
CXXFLAGS = -c -g -Wall
LIBS = -lboost_math_tr1
PROGS = boost_test  

.PHONY: all
all: $(PROGS)

boost_test: boost_test.o 
	$(CXX) -o $@ $^ $(LIBS) 

.PHONY: clean
clean:
	$(RM) -r *.o
	$(RM) .depend

%.o: %.cc
	$(CXX) $(CXXFLAGS) $^ -o $@ 

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

