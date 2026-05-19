CXX ?= g++

TARGET := optim_par
SRC := optim_par.cpp ../quadgrid/src/particles.cpp ../timer/src/timer.cpp
OBJ := $(SRC:.cpp=.o)

INCLUDES := -I. -I../quadgrid/include -I../timer/src -I../json/single_include/nlohmann
CXXFLAGS ?= -g -std=c++20 -Wall -Wextra
LDFLAGS ?=
LDLIBS ?=
TBB_LIBDIR ?=

DEFS :=

# Enable with: make USE_DPL=1
ifeq ($(USE_DPL),1)
DEFS += -DUSE_DPL
LDLIBS += -ltbb
endif

# Enable with: make USE_TBB=1 [TBB_LIBDIR=/path/to/tbb/lib]
ifeq ($(USE_TBB),1)
LDLIBS += -ltbb
ifneq ($(strip $(TBB_LIBDIR)),)
LDFLAGS += -L$(TBB_LIBDIR) -Wl,-rpath,$(TBB_LIBDIR)
endif
endif

# Enable with: make USE_COMPRESSION=1
ifeq ($(USE_COMPRESSION),1)
DEFS += -DUSE_COMPRESSION
LDLIBS += -lboost_iostreams -lz
endif

.PHONY: all clean distclean ditclean run

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCLUDES) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)

distclean: clean
	rm -f $(OBJ)

