CXX := g++
CXXFLAGS := -O3 -g -fno-omit-frame-pointer -Werror -fopenmp
LDFLAGS := -lm
TARGET := main.exe
CPP_FILES := $(wildcard *.cpp)
PERF_REPORT := perf_report.txt

all:
	$(CXX) $(CXXFLAGS) $(CPP_FILES) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)

flamegraph: all
	perf stat -o $(PERF_REPORT) ./$(TARGET) 1000 5000
	perf record -F 99 -g ./$(TARGET) 1000 5000
	perf script | ~/Flamegraph/stackcollapse-perf.pl | ~/Flamegraph/flamegraph.pl > flamegraph.svg