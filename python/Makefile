# --- Compiler and Toolchain ---
CXX := g++

# --- Configuration ---
TARGET     := ../bin/bca_reader.so
VSCRIPT    := linker.version


PYINCLUDE  := $(shell python3-config --includes)

# --- Build Flags --- 
# Reseek needs to be built with either -fPIC or -flto. LTO enables PIC to be
# emitted at link time even when -fPIC was not used when building the objects.
RESEEK_CXXFLAGS := $(shell make -C ../src -f ../src/Makefile -pn | grep '^CXXFLAGS' | tail -n 1 | cut -d= -f2-)
CXXFLAGS   := $(RESEEK_CXXFLAGS) $(PYINCLUDE) -fPIC

# Linker flags.
LDFLAGS    := -shared -Wl,--version-script=$(VSCRIPT)

# --- Source Files ---
WRAPPER_SRC := entry_point.cpp
CORE_OBJ_DIR := ../src/o
CORE_OBJS   := $(CORE_OBJ_DIR)/pdbchain.o \
               $(CORE_OBJ_DIR)/dss.o \
               $(CORE_OBJ_DIR)/dssparams.o \
               $(CORE_OBJ_DIR)/alpha.o \
               $(CORE_OBJ_DIR)/features.o \
               $(CORE_OBJ_DIR)/myutils.o \
               $(CORE_OBJ_DIR)/myss.o \
               $(CORE_OBJ_DIR)/mymalloc.o \
               $(CORE_OBJ_DIR)/valuetoint.o \
               $(CORE_OBJ_DIR)/getss.o \
               $(CORE_OBJ_DIR)/namedparams.o \
               $(CORE_OBJ_DIR)/trained_features.o \
               $(CORE_OBJ_DIR)/bcadata.o

# --- Object File Targets ---
WRAPPER_OBJ := $(patsubst %.cpp,%.o,$(WRAPPER_SRC))
ALL_OBJS    := $(CORE_OBJS) $(WRAPPER_OBJ)

# --- Build Rules ---

.PHONY: all clean check-exports

all: $(TARGET)

# Final Link Rule: Create the .so from all object files.
$(TARGET): $(ALL_OBJS) $(VSCRIPT)
	@echo "--- [LINK] Creating final shared library: $(TARGET) ---"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(ALL_OBJS)

# Compilation Rule: Compile the Python wrapper source to an LTO object file.
$(WRAPPER_OBJ): $(WRAPPER_SRC)
	@echo "--- [COMPILE] Building wrapper object: $@ ---"
	$(CXX) $(CXXFLAGS) -I../ -c $(WRAPPER_SRC) -o $@


# Clean Rule
clean:
	rm -f $(TARGET) $(WRAPPER_OBJ)