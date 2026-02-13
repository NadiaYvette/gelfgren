# Makefile for Gelfgren Fortran Bindings

FC = gfortran
FFLAGS = -Wall -Wextra -O2 -fPIC
LDFLAGS = -L../../target/release -lgelfgren

# Directories
SRC_DIR = src
EXAMPLE_DIR = example
TEST_DIR = test
BUILD_DIR = build
LIB_DIR = ../../target/release
INCLUDE_DIR = ../../include

# Output
MODULE_NAME = gelfgren_mod.o
EXAMPLE_EXEC = basic_usage
TEST_EXEC = test_bernstein

.PHONY: all module example test clean

all: module example

module: $(BUILD_DIR)/$(MODULE_NAME)

$(BUILD_DIR)/$(MODULE_NAME): $(SRC_DIR)/gelfgren_mod.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILD_DIR)

example: $(BUILD_DIR)/$(EXAMPLE_EXEC)

$(BUILD_DIR)/$(EXAMPLE_EXEC): $(EXAMPLE_DIR)/basic_usage.f90 $(BUILD_DIR)/$(MODULE_NAME) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -I$(BUILD_DIR) $< $(BUILD_DIR)/$(MODULE_NAME) -o $@ $(LDFLAGS)

test: $(BUILD_DIR)/$(TEST_EXEC)

$(BUILD_DIR)/$(TEST_EXEC): $(TEST_DIR)/test_bernstein.f90 $(BUILD_DIR)/$(MODULE_NAME) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -I$(BUILD_DIR) $< $(BUILD_DIR)/$(MODULE_NAME) -o $@ $(LDFLAGS)

run-example: example
	@echo "Running example..."
	@export LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH && \
	 export DYLD_LIBRARY_PATH=$(LIB_DIR):$$DYLD_LIBRARY_PATH && \
	 $(BUILD_DIR)/$(EXAMPLE_EXEC)

run-test: test
	@echo "Running tests..."
	@export LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH && \
	 export DYLD_LIBRARY_PATH=$(LIB_DIR):$$DYLD_LIBRARY_PATH && \
	 $(BUILD_DIR)/$(TEST_EXEC)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

help:
	@echo "Gelfgren Fortran Bindings Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build module and example (default)"
	@echo "  module       - Build the Fortran module"
	@echo "  example      - Build the example program"
	@echo "  test         - Build the test program"
	@echo "  run-example  - Run the example program"
	@echo "  run-test     - Run the test program"
	@echo "  clean        - Remove build artifacts"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Prerequisites:"
	@echo "  - gfortran compiler"
	@echo "  - Gelfgren C library built (cargo build --release in root)"
