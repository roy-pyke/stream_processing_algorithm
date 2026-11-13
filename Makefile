CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
BIN = bin/bench bin/gen_synth bin/mkworkload bin/pcap2bin

all: $(BIN)

bin/%: src/%.cpp src/common.hpp src/synopses.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -rf bin
