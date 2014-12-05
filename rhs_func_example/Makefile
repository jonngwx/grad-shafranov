CXX = g++
CXXFLAGS = -g -Wall -std=c++11
EXAMPLEPROGS = rhs_test

.PHONY: all
all: $(EXAMPLEPROGS)

rhs_test: rhs_test.cc tsv_reader.cc rhs_func.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) *.o
