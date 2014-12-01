CXX = clang++
CXXFLAGS = -g -Wall -std=c++11

.PHONY: all
all: csv_reader

csv_read: csv_reader.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) *.o
