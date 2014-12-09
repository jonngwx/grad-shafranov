CXX = g++
CXXFLAGS = -c -g -Wall -std=c++11 -Iinclude
LIBS = -lboost_program_options
PROGS = gs_solver

.PHONY: all
all: $(PROGS)

gs_solver: gs_solver.o tsv_reader.o rhs_func.o grid.o field.o slow_boundary.o grad_output.o grad_output_txt.o grad_output_hdf.o create_options.o elliptic/sor.o elliptic/elliptic_solver.o
	$(CXX) -o $@ $^ $(LIBS) 

.PHONY: clean
clean:
	$(RM) -r *.o
	$(RM) .depend

%.o: %.cc
	$(CXX) $(CXXFLAGS) $^ -o $@ 

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
