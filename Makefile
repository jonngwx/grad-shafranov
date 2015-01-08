CXX = g++
CXXFLAGS = -c -g -std=c++11 -Iinclude -Wall 
LIBS = -lboost_program_options -lboost_unit_test_framework -lboost_math_tr1
HDF = -lhdf5_hl -lhdf5 -DHDF_MODE
PROGS = gs_solver
TESTDIR = test
TESTPROGS = output_test elliptic_test $(TESTDIR)/tsv_reader_test
EXAMPLEPROGS = tsv_reader_example coil_data_example
OBJECTS = tsv_reader.o j_solver_alpha.o grid.o field.o boundary.o slow_boundary.o grad_output.o grad_output_txt.o create_options.o elliptic/sor.o elliptic/elliptic_solver.o elliptic/gauss_seidel.o elliptic/critical.o green_fcn.o elliptic/interpolate.o

.PHONY: all
all: $(PROGS) $(EXAMPLEPROGS) $(TESTPROGS)

gs_solver: gs_solver.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) 

gs_solver_hdf: gs_solver_hdf.o grad_output_hdf.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(HDF)

gs_solver_hdf.o: gs_solver.cc
	$(CXX) $(CXXFLAGS) $(HDF) $^ -o $@

elliptic_test: elliptic/elliptic_solver.o elliptic/elliptic_test.o elliptic/sor.o grid.o field.o elliptic/gauss_seidel.o boundary.o
	$(CXX) -o $@ $^ $(LIBS)

critical_test: elliptic/critical.o elliptic/critical_test.o grid.o field.o elliptic/interpolate.o
	$(CXX) -o $@ $^ $(LIBS)

interpolate_test: grid.o field.o elliptic/interpolate.o elliptic/interpolate_test.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/tsv_reader_test: tsv_reader.o $(TESTDIR)/tsv_reader_test.o
	$(CXX) -o $@ $^ $(LIBS)

output_test: test/test_output.o grad_output.o grad_output_txt.o field.o grid.o
	$(CXX) -o $@ $^ $(LIBS)

tsv_reader_example: tsv_reader_example.o tsv_reader.o 
	$(CXX) -o $@ $^ $(LIBS)

coil_data_example: coil_data_example.o tsv_reader.o
	$(CXX) -o $@ $^ $(LIBS)

.PHONY: runtests
runtests:
	# Running tests!
	test/tsv_reader_test

.PHONY: clean
clean:
	# delete .o files in subdirectories. rm -r *.o does not do this.
	find . -name '*.o' -delete
	$(RM) .depend

%.o: %.cc
	$(CXX) $(CXXFLAGS) $^ -o $@ 

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
