CXX = g++
CXXFLAGS = -c -g -Wall -std=c++11 -Iinclude
LIBS = -lboost_program_options -lboost_unit_test_framework
HDF = -lhdf5_hl -lhdf5 -DHDF_MODE
PROGS = gs_solver

.PHONY: all
all: $(PROGS)

gs_solver: gs_solver.o tsv_reader.o rhs_func.o grid.o field.o slow_boundary.o grad_output.o grad_output_txt.o create_options.o elliptic/sor.o elliptic/elliptic_solver.o elliptic/gauss_seidel.o
	$(CXX) -o $@ $^ $(LIBS) 

gs_solver_hdf: gs_solver_hdf.o tsv_reader.o rhs_func.o grid.o field.o slow_boundary.o grad_output.o grad_output_txt.o grad_output_hdf.o create_options.o elliptic/sor.o elliptic/elliptic_solver.o
	$(CXX) -o $@ $^ $(LIBS) $(HDF)

gs_solver_hdf.o: gs_solver.cc
	$(CXX) $(CXXFLAGS) $(HDF) $^ -o $@

elliptic_test: elliptic/elliptic_solver.o elliptic/elliptic_test.o elliptic/sor.o grid.o field.o elliptic/gauss_seidel.o
	$(CXX) -o $@ $^ $(LIBS) 

tsv_reader_example: tsv_reader_example.o tsv_reader.o 
	$(CXX) -o $@ $^ $(LIBS)

coil_data_example: coil_data_example.o tsv_reader.o
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
