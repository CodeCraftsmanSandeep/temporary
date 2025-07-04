# Compiler and flags
CXX = g++
CXXFLAGS = -O3 -std=c++17 -I../include

# Method 1
METHOD1_SRC = method-1/main.cpp
METHOD1_OUT = method-1/exe.out

# Method 2
METHOD2_SRC = method-2/main.cpp
METHOD2_OUT = method-2/exe.out

# Method 3
METHOD3_SRC = method-3/main.cpp
METHOD3_OUT = method-3/exe.out

# Method 4
METHOD4_SRC = method-4/main.cpp
METHOD4_OUT = method-4/exe.out

# Method 5
METHOD5_SRC = method-5/main.cpp
METHOD5_OUT = method-5/exe.out

# Default: build all methods
all: method-1 method-2 method-3 method-4 method-5

# Build method-1 executable
method-1: $(METHOD1_OUT)
$(METHOD1_OUT): $(METHOD1_SRC)
	$(CXX) $(CXXFLAGS) $(METHOD1_SRC) -o $(METHOD1_OUT)

# Build method-2 executable
method-2: $(METHOD2_OUT)
$(METHOD2_OUT): $(METHOD2_SRC)
	$(CXX) $(CXXFLAGS) $(METHOD2_SRC) -o $(METHOD2_OUT)

# Build method-3 executable
method-3: $(METHOD3_OUT)
$(METHOD3_OUT): $(METHOD3_SRC)
	$(CXX) $(CXXFLAGS) $(METHOD3_SRC) -o $(METHOD3_OUT)

# Build method-4 executable
method-4: $(METHOD4_OUT)
$(METHOD4_OUT): $(METHOD4_SRC)
	$(CXX) $(CXXFLAGS) $(METHOD4_SRC) -o $(METHOD4_OUT)

# Build method-5 executable
method-5: $(METHOD5_OUT)
$(METHOD5_OUT): $(METHOD5_SRC)
	$(CXX) $(CXXFLAGS) $(METHOD5_SRC) -o $(METHOD5_OUT)

# Run targets
run-method-1: method-1
	./$(METHOD1_OUT)

run-method-2: method-2
	./$(METHOD2_OUT)

run-method-3: method-3
	./$(METHOD3_OUT)

run-method-4: method-4
	./$(METHOD4_OUT)

run-method-5: method-5
	./$(METHOD5_OUT)

# Clean up
clean:
	rm -f $(METHOD1_OUT) $(METHOD2_OUT) $(METHOD3_OUT) $(METHOD4_OUT) $(METHOD5_OUT)
