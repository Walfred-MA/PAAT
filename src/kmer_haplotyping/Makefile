# Variables
CXX = g++
CXXFLAGS = -Wall -std=c++17 -O3
LDFLAGS = -lz -pthread -lstdc++fs
SOURCES = src/*.cpp
TARGET = kmer_convertor

# Default target
all: $(TARGET)

$(TARGET):
	$(CXX) $(CXXFLAGS) $(SOURCES) $(LDFLAGS) -o $@

# Clean target
clean:
	rm -f $(TARGET)
