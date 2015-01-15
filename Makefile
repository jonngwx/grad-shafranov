CXX = g++
CXXFLAGS = -c -g -std=c++11 -Iinclude -Wall 
LIBS = -lboost_program_options -lboost_unit_test_framework -lboost_math_tr1 
HDF = -lhdf5_hl -lhdf5 -DHDF_MODE
PROGS = gs_solver
TESTDIR = test
TESTPROGS = $(TESTDIR)/output_test $(TESTDIR)/util_test $(TESTDIR)/elliptic_test $(TESTDIR)/tsv_reader_test $(TESTDIR)/grid_test $(TESTDIR)/boundary_test $(TESTDIR)/slow_boundary_test $(TESTDIR)/interpolate_test
EXDIR = exampleClassUsage
EXAMPLEPROGS = $(EXDIR)/tsv_reader_example $(EXDIR)/coil_data_example
OBJECTS = tsv_reader.o j_solver_alpha.o grid.o field.o boundary.o slow_boundary.o grad_output.o grad_output_txt.o create_options.o elliptic/sor.o elliptic/elliptic_solver.o elliptic/gauss_seidel.o elliptic/critical.o green_fcn.o elliptic/interpolate.o j_solver_nstx.o
TESTOBJECTS = $(addsuffix .o,$(TESTPROGS)) 

.PHONY: all
all: $(PROGS) $(EXAMPLEPROGS) $(TESTPROGS) $(TESTDIR)/all_test

gs_solver: gs_solver.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) 

gs_solver_hdf: gs_solver_hdf.o grad_output_hdf.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(HDF)

gs_solver_hdf.o: gs_solver.cc
	$(CXX) $(CXXFLAGS) $(HDF) $^ -o $@

$(TESTDIR)/all_test: $(TESTDIR)/test.o $(OBJECTS) $(TESTOBJECTS) 
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/elliptic_test: $(TESTDIR)/test.o elliptic/elliptic_solver.o test/elliptic_test.o elliptic/sor.o grid.o field.o elliptic/gauss_seidel.o boundary.o 
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/critical_test: $(TESTDIR)/test.o elliptic/critical.o elliptic/critical_test.o grid.o field.o elliptic/interpolate.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/interpolate_test: $(TESTDIR)/test.o grid.o field.o elliptic/interpolate.o $(TESTDIR)/interpolate_test.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/tsv_reader_test: $(TESTDIR)/test.o tsv_reader.o $(TESTDIR)/tsv_reader_test.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/output_test: $(TESTDIR)/test.o $(TESTDIR)/output_test.o grad_output.o grad_output_txt.o field.o grid.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/util_test: $(TESTDIR)/test.o $(TESTDIR)/util_test.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/grid_test: $(TESTDIR)/test.o $(TESTDIR)/grid_test.o grid.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/boundary_test: $(TESTDIR)/test.o $(TESTDIR)/boundary_test.o boundary.o grid.o field.o
	$(CXX) -o $@ $^ $(LIBS)

$(TESTDIR)/slow_boundary_test: $(TESTDIR)/test.o $(TESTDIR)/slow_boundary_test.o slow_boundary.o boundary.o grid.o field.o tsv_reader.o green_fcn.o elliptic/elliptic_solver.o elliptic/gauss_seidel.o
	$(CXX) -o $@ $^ $(LIBS)

$(EXDIR)/tsv_reader_example: $(EXDIR)/tsv_reader_example.o tsv_reader.o 
	$(CXX) -o $@ $^ $(LIBS)

$(EXDIR)/coil_data_example: $(EXDIR)/coil_data_example.o tsv_reader.o
	$(CXX) -o $@ $^ $(LIBS)

.PHONY: runtests
runtests:
	# Running tests!
	./test/all_test

.PHONY: clean
clean:
	find . -name '*.o' -delete
	$(RM) .depend

%.o: %.cc
	$(CXX) $(CXXFLAGS) $^ -o $@ 

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
