CXX = clang++
CXXFLAGS = -g -Wall -std=c++11
EXAMPLEPROGS = tsv_reader_example coil_data_example

.PHONY: all
all: $(EXAMPLEPROGS)

tsv_reader_example: tsv_reader_example.cc tsv_reader.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

coil_data_example: coil_data_example.cc tsv_reader.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) $(EXAMPLEPROGS)
