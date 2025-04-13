CXX := g++
CXXFLAGS := -O3 -g -fno-omit-frame-pointer -Werror
LDFLAGS := -lm
TARGET := main.exe
CPP_FILES := $(wildcard *.cpp)

all:
	$(CXX) $(CXXFLAGS) $(CPP_FILES) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) perf.data flamegraph.svg perf.data.old

flamegraph: all
	perf record -F 99 -g ./$(TARGET) 1000 5000
	perf script | ~/Flamegraph/stackcollapse-perf.pl | ~/Flamegraph/flamegraph.pl > flamegraph.svg
