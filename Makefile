# =====================================================================
# PHASIS -- alternativa ao CMake.
#
# O CMakeLists.txt e o build canonico. Este Makefile existe porque o
# cmake nem sempre esta no PATH, e o projeto nao tem nenhuma dependencia
# que justifique exigi-lo: g++ com C++17 e a STL bastam.
#
#   make          compila
#   make test     compila e roda a suite de aceitacao
#   make clean
# =====================================================================

CXX      ?= g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -fopenmp -Iinclude -MMD -MP

BUILD := build
SRC   := src/shell.cpp src/cascade.cpp src/ode.cpp src/sweep.cpp src/integrate.cpp src/metric.cpp src/density.cpp \
         src/cross_section.cpp src/trace.cpp
OBJ   := $(SRC:%.cpp=$(BUILD)/%.o)

.PHONY: all test clean
all: $(BUILD)/trace $(BUILD)/test_phase1 $(BUILD)/test_phase2 $(BUILD)/test_phase3 $(BUILD)/test_phase4 $(BUILD)/test_phase5

$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/trace: $(OBJ) $(BUILD)/src/trace_main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD)/test_phase1: $(OBJ) $(BUILD)/tests/test_phase1.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD)/test_phase2: $(OBJ) $(BUILD)/tests/test_phase2.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD)/test_phase3: $(OBJ) $(BUILD)/tests/test_phase3.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD)/test_phase4: $(OBJ) $(BUILD)/tests/test_phase4.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD)/test_phase5: $(OBJ) $(BUILD)/tests/test_phase5.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test: $(BUILD)/test_phase1 $(BUILD)/test_phase2 $(BUILD)/test_phase3 $(BUILD)/test_phase4 $(BUILD)/test_phase5 $(BUILD)/test_phase5
	@./$(BUILD)/test_phase1
	@./$(BUILD)/test_phase2
	@./$(BUILD)/test_phase3
	@./$(BUILD)/test_phase4
	@./$(BUILD)/test_phase5

clean:
	rm -rf $(BUILD)

-include $(wildcard $(BUILD)/*/*.d)
