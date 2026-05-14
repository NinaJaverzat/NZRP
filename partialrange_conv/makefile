# Compiler
CXX ?= g++

# Flags
OPTIMIZATION_FLAGS = -O3 -march=native -flto
DBG_FLAGS = -g -g2
WARN_FLAGS = -Wall -Wextra -pedantic

# Compiler flags
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

# Source file for convolutor.exe
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = partialrange_conv.exe

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
