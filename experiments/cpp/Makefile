# Linux or MSYS2 UCRT64 (Windows). Requirements: g++ (C++17), GMP, OpenSSL 3 (libcrypto), make
# SPP points to the library sources (default: ../../spp).
SPP      ?= ../../spp
CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -march=x86-64
GMP_INC  ?=
GMP_LIB  ?=
CPPFLAGS += -I$(SPP)/include
ifneq ($(GMP_INC),)
  CPPFLAGS += -I$(GMP_INC)
endif
ifneq ($(GMP_LIB),)
  LDFLAGS += -L$(GMP_LIB)
  ifneq ($(OS),Windows_NT)
    LDFLAGS += -Wl,-rpath,$(GMP_LIB)
  endif
endif
ifeq ($(OS),Windows_NT)
  EXTRA := -lbcrypt
  EXE := .exe
endif
all: diss_bench$(EXE) ct_experiment$(EXE) ossl_bench$(EXE) safe_batch$(EXE) rsa_bench$(EXE) ossl_seg$(EXE) fullgen$(EXE) review_exp$(EXE)
diss_bench$(EXE): diss_bench.cpp $(SPP)/src/spp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -lgmp $(EXTRA) -o $@
ct_experiment$(EXE): ct_experiment.cpp $(SPP)/src/spp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -lgmp $(EXTRA) -o $@
ossl_bench$(EXE): ossl_bench.cpp
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -lcrypto -o $@
safe_batch$(EXE): safe_batch.cpp $(SPP)/src/spp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -lgmp $(EXTRA) -o $@
rsa_bench$(EXE): rsa_bench.cpp $(SPP)/src/spp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Wno-deprecated-declarations $^ $(LDFLAGS) -lgmp -lcrypto $(EXTRA) -o $@
ossl_seg$(EXE): ossl_seg.cpp
	$(CXX) $(CXXFLAGS) -Wno-deprecated-declarations $^ $(LDFLAGS) -lcrypto -o $@
fullgen$(EXE): fullgen.cpp $(SPP)/src/spp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -lgmp $(EXTRA) -o $@
review_exp$(EXE): review_exp.cpp $(SPP)/src/spp.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -lgmp $(EXTRA) -o $@
# the compiled-code (ctgrind) analysis is built and run by ctgrind/run_ctgrind.sh (requires Valgrind, g++, clang++)
# optional: cost of BN_gcd and BN_mod_word in OpenSSL (motivates the word-sized prefilter inside OpenSSL)
gcdcost$(EXE): gcdcost.cpp
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -lcrypto -o $@
clean:
	rm -f diss_bench$(EXE) ct_experiment$(EXE) ossl_bench$(EXE) safe_batch$(EXE) rsa_bench$(EXE) ossl_seg$(EXE) fullgen$(EXE) review_exp$(EXE) gcdcost$(EXE)
