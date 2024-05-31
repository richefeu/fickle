
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CXX = g++-14
  CXXFLAGS = -O3 -Wall -Wextra -std=c++17 -I ~/toofus
  GLLINK = `pkg-config --libs glut` -framework OpenGL
	GLFLAGS = `pkg-config --cflags glut`
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++17 -I ~/toofus
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# The list of source files
SOURCES = Particle.cpp Interaction.cpp Loading.cpp PeriodicCell.cpp PBC.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all clean

all: PBC see

clean:
	rm -f *.o
	rm -f PBC see exam DIC2conf
	rm -f libPBC.a

%.o:%.cpp
	@echo "\033[0;32m-> COMPILING OBJECT" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@

libPBC.a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING LIBRARY" $@ "\033[0m"
	ar rcs $@ $^
	
PBC: run.cpp libPBC.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o run.o
	$(CXX) -o $@ run.o libPBC.a
	
see: see.cpp libPBC.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFLAGS)
	$(CXX) -o $@ see.o libPBC.a $(GLLINK)
	
exam: FluctExamination.cpp CDF.hpp libPBC.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o FluctExamination.o
	$(CXX) -o $@ FluctExamination.o libPBC.a

DIC2conf: DIC2conf.cpp DICfile.hpp libPBC.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o DIC2conf.o
	$(CXX) -o $@ DIC2conf.o libPBC.a
	
	
	