CXX = g++
CXXFLAGS = -g -Wall -std=c++11
PROGS = gs_solver

.PHONY: all
all: $(PROGS)

gs_solver: gs_solver.cc tsv_reader.cc rhs_func.cc grid.cc field.cc slow_boundary.cc grad_output.cc grad_output_txt.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ -Iinclude

.PHONY: clean
clean:
	$(RM) *.o
