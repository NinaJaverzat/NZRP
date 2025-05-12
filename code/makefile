# Compiler
CXX = g++

# Flags
OPTIMIZATION_FLAGS = -O3 -march=native -flto -ffast-math
DBG_FLAGS = -g -g2
WARN_FLAGS = -Wall -Wextra -pedantic

# Compiler flags
CXXFLAGS = -std=c++20

# Enable flags based on make options
ifeq ($(g),1)
    CXXFLAGS += $(DBG_FLAGS)
endif

ifeq ($(w),1)
    CXXFLAGS += $(WARN_FLAGS)
endif

ifeq ($(o),1)
    CXXFLAGS += $(OPTIMIZATION_FLAGS)
endif

# Source files for rp.exe
SRCS_RP = main.cpp basic_functions.cpp percolation.cpp pebble_game.cpp
OBJS_RP = $(SRCS_RP:.cpp=.o)
TARGET_RP = rp.exe

# Default rule to build all executables
all: $(TARGET_RP)

# Rule to build rp.exe
$(TARGET_RP): $(OBJS_RP)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(OBJS_RP) $(TARGET_RP)

clear:
	rm -f $(OBJS_RP)

# Phony targets
.PHONY: all clean clear
