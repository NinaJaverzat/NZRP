# Compiler
CXX = g++ -O3 -ffast-math -march=native

DBG_FLAGS = -g -g2 -pg -ggdb3
WARN_FLAGS = -Wall -Wextra -pedantic

# Compiler flags
CXXFLAGS = -std=c++20

# Check for each flag condition
ifeq ($(g),1)
    CXXFLAGS += $(DBG_FLAGS)
endif

ifeq ($(w),1)
    CXXFLAGS += $(WARN_FLAGS)
endif

# List of source files
SRCS = main.cpp

# List of object files (derived from the source files)
OBJS = $(SRCS:.cpp=.o)

# Name of the executable
TARGET = convolutor.exe

# Default rule to build the program
all: $(TARGET)

# Rule to link the object files into the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile the source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up the object files and the executable
clean:
	rm -f $(OBJS) $(TARGET)

clear:
	rm -f $(OBJS)

# Phony targets to avoid conflicts with files named 'all' or 'clean' or 'clean'
.PHONY: all clean clear
