# Simple build without CMake (Linux, macOS, MSYS2/MinGW on Windows).
#   make            -> build/libspp.a and the tools
#   make test       -> run unit tests
#   make GMP_INC=/opt/gmp/include GMP_LIB=/opt/gmp/lib   (non-standard GMP location)
CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
AR       ?= ar
GMP_INC  ?=
GMP_LIB  ?=
CPPFLAGS += -Iinclude
ifneq ($(GMP_INC),)
  CPPFLAGS += -I$(GMP_INC)
endif
ifneq ($(GMP_LIB),)
  LDFLAGS += -L$(GMP_LIB)
  ifneq ($(OS),Windows_NT)
    LDFLAGS += -Wl,-rpath,$(GMP_LIB)
  endif
endif
LDLIBS   += -lgmp
ifeq ($(OS),Windows_NT)
  LDLIBS += -lbcrypt
  EXE := .exe
endif

B := build
TOOLS := $(B)/spp_check$(EXE) $(B)/spp_gen$(EXE) $(B)/spp_bench$(EXE) $(B)/spp_example$(EXE)

all: $(B)/libspp.a $(TOOLS) $(B)/spp_tests$(EXE)

$(B):
	mkdir -p $(B)
$(B)/spp.o: src/spp.cpp include/spp/spp.hpp | $(B)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
$(B)/libspp.a: $(B)/spp.o
	$(AR) rcs $@ $^
$(B)/spp_%$(EXE): tools/spp_%.cpp tools/common.hpp $(B)/libspp.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(B)/libspp.a $(LDFLAGS) $(LDLIBS) -o $@
$(B)/spp_example$(EXE): examples/example.cpp $(B)/libspp.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(B)/libspp.a $(LDFLAGS) $(LDLIBS) -o $@
$(B)/spp_tests$(EXE): tests/spp_tests.cpp $(B)/libspp.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(B)/libspp.a $(LDFLAGS) $(LDLIBS) -o $@

test: $(B)/spp_tests$(EXE)
	./$(B)/spp_tests$(EXE)

clean:
	rm -rf $(B)

.PHONY: all test clean
