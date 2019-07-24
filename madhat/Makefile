#
# Please modify only BOOST to match your version
#
# Boost C++ library
BOOST = ../boost_1_70_0
CXX := g++
CFLAGS = -g -O2
CPPFLAGS = -c -std=c++11
DFLAGS = 
LIBS   =
BUILD_DIR := Output

# Where the header are (no need to modify this)
INCLUDE = src

EXEC         = madhat
SOURCES      = $(wildcard src/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
OBJECTS_TEST = $(SOURCES_TEST:.cpp=.o)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)
OBJDIR := Output

$(EXEC): $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)
	$(CXX) $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_LIB)  $(LIBS) -o $(EXEC)
	mkdir -p Output

%.o: %.cpp
	$(CXX) $(CPPFLAGS)  -I src/stats -I $(INCLUDE) -I ${BOOST} $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_TEST)
