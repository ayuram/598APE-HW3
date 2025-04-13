CXX := g++
CXXFLAGS := -O3 -g -fno-omit-frame-pointer -Werror -fopenmp
LDFLAGS := -lm
TARGET := main.exe
CPP_FILES := $(wildcard *.cpp)
PERF_REPORT := perf_report.txt

all:
	$(CXX) $(CXXFLAGS) $(CPP_FILES) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) perf.data flamegraph.svg perf.data.old $(PERF_REPORT)
