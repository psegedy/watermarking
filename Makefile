CXX=c++
CXXFLAGS=-g -pedantic -Wall -Wextra -std=c++11
SOURCES=watermarking.cpp watermarking.hpp
EXECUTABLE=watermarking
MKDIR_P=mkdir -p
OUT_DIR=attacks
MAGICK=`Magick++-config --cppflags --cxxflags --ldflags --libs`

all: directory $(EXECUTABLE)

watermarking: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@ $(MAGICK)

directory: ${OUT_DIR}

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

clean:
	rm -rf $(EXECUTABLE) ${OUT_DIR}

.PHONY: directory