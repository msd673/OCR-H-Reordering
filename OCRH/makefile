CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra
TARGET = OCRH

SRCS = main.cpp \
       csr.cpp \
       sampling.cpp \
       kcmer_minhash.cpp \
       clustering.cpp

OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean