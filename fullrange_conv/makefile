# Compiler
CXX ?= g++

DBG_FLAGS = -g -g2 -pg -ggdb3
WARN_FLAGS = -Wall -Wextra -pedantic
OPTIMIZATION_FLAGS = -O3 -ffast-math -march=native

CXXFLAGS = -std=c++20

ifeq ($(g),1)
    CXXFLAGS += $(DBG_FLAGS)
endif

ifeq ($(w),1)
    CXXFLAGS += $(WARN_FLAGS)
endif

ifeq ($(o),1)
    CXXFLAGS += $(OPTIMIZATION_FLAGS)
endif

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = fullrange_conv.exe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

clear:
	rm -f $(OBJS)

.PHONY: all clean clear
