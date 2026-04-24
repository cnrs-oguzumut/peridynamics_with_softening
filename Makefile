CXX ?= clang++
TARGET ?= peridynamics_with_softening
BUILD_DIR ?= build

SRC := main.cpp
OBJ := $(BUILD_DIR)/$(SRC:.cpp=.o)
DEPS := $(wildcard *.h)

CPPFLAGS ?=
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra
LDFLAGS ?=
LDLIBS ?=

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(LIBOMP_PREFIX),)
    CPPFLAGS += -I$(LIBOMP_PREFIX)/include
    CXXFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib
    LDLIBS += -lomp
  endif
else
  CXXFLAGS += -fopenmp
endif

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp $(DEPS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
