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

flamegraph: all
	./$(TARGET) 1000 5000 &
	sleep 1
	PID=$$(pidof $(TARGET)); \
	perf stat --per-thread -p $$PID -o $(PERF_REPORT); \
	wait $$PID
	perf record -F 99 -g ./$(TARGET) 1000 5000
	perf script | ~/Flamegraph/stackcollapse-perf.pl | ~/Flamegraph/flamegraph.pl > flamegraph.svg
