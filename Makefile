CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -I.
LDFLAGS = -flto

# Link against existing compiled reseek object files
RESEEK_OBJS = src/o/pdbchain.o \
src/o/dss.o \
src/o/dssparams.o \
src/o/alpha.o \
src/o/mx.o \
src/o/myutils.o \
src/o/myss.o \
src/o/mymalloc.o \
src/o/features.o \
src/o/valuetoint.o \
src/o/getss.o \
src/o/namedparams.o \
src/o/float_feature_bins.o \
src/o/getbins.o \
src/o/abcxyz.o \
src/o/trained_features.o \
src/o/bcadata.o

TARGET = bca_reader
SOURCE = bca_reader.cpp

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE) $(RESEEK_OBJS) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
