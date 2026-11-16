CXX ?= g++
# Apple clang on arm64 wants -mcpu=native
ARCHFLAG ?= $(shell uname -m | grep -q arm64 && echo -mcpu=native || echo -march=native)
CXXFLAGS ?= -O3 $(ARCHFLAG) -std=c++17 -Wall -Wextra
BIN = bin/bench bin/gen_synth bin/mkworkload bin/pcap2bin bin/jhash_test bin/bias_probe bin/mkdelstream

all: $(BIN)

bin/%: src/%.cpp src/common.hpp src/synopses.hpp src/omni_java_compat.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -rf bin
