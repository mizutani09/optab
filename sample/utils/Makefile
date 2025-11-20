CXX = h5c++
CXXFLAGS = -std=c++11 -Wall -O2
HDF5_DIR = /usr/local/hdf5_par

# Use parallel HDF5 if available
ifneq (,$(wildcard $(HDF5_DIR)/bin/h5c++))
    CXX = $(HDF5_DIR)/bin/h5c++
else
    CXX = h5c++
    CXXFLAGS += -I$(HDF5_DIR)/include
    LDFLAGS += -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl
endif

TARGET = test_hdf5
SOURCE = test.cpp

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
