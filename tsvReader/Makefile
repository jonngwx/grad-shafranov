CXX = clang++
CXXFLAGS = -g -Wall -std=c++11

.PHONY: all
all: tsv_reader

csv_read: tsv_reader.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) *.o
