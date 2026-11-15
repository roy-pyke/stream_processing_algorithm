CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
BIN = bin/bench bin/gen_synth bin/mkworkload bin/pcap2bin bin/jhash_test

all: $(BIN)

bin/%: src/%.cpp src/common.hpp src/synopses.hpp src/omni_java_compat.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -rf bin
