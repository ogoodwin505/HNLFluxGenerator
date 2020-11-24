CPPFLAGS= -I $(ROOTSYS)/include\
         -I ${DK2NU}/include
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Werror -pedantic

LDFLAGS=$$(root-config --libs --cflags) \
	-l EG \
        -L${DK2NU}/lib -ldk2nuTree

%.o: srs/%.c++
	$(CXX) -o $@ $< $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS)

HNLFluxGen: HNLFluxGen.o Parameter.o HNLFlux_funcs.o
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

all: HNLFluxGen

clean:
	rm *.o HNLFluxGen

