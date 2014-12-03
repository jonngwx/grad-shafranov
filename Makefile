CXX = clang++
CXXFLAGS = -g -Wall -std=c++11

.PHONY: all
all: tsv_reader_example

tsv_reader_example: tsv_reader_example.cc tsv_reader.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ 

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) tsv_reader_example
