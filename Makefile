CXX = g++
CXXFLAGS = -g -Wall -std=c++11
LIBS = -lboost_program_options
PROGS = gs_solver

.PHONY: all
all: $(PROGS)

gs_solver: gs_solver.cc tsv_reader.cc rhs_func.cc grid.cc field.cc slow_boundary.cc grad_output.cc grad_output_txt.cc create_options.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) -Iinclude

.PHONY: clean
clean:
	$(RM) *.o

