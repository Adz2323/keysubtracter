CC = gcc
CXX = g++
CFLAGS = -O3
CXXFLAGS = -O3 -std=c++11 -msha
LIBS = -lgmp -lstdc++

C_SRCS = base58/base58.c xxhash/xxhash.c
CPP_SRCS = bloom.cpp hash/sha256.cpp hash/ripemd160.cpp keysubtracter.cpp util.cpp gmpecc.cpp
OBJS = $(C_SRCS:.c=.o) $(CPP_SRCS:.cpp=.o)

# Default target
default: clean keysubtracter

# Build target for keysubtracter
keysubtracter: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Rule for compiling .c files to .o files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule for compiling .cpp files to .o files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJS) keysubtracter
