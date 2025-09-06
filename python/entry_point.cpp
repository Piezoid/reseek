#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <pybind11/numpy.h>

// Include necessary reseek headers (same as bca_reader.cpp)
#include "src/pdbchain.h"
#include "src/dss.h"
#include "src/dssparams.h"
#include "src/bcadata.h"
#include "src/features.h"

namespace py = pybind11;

// Prevent assert exiting the python process
void Die_(const char *Format, ...)
{
    // 1. Format the error message from the variadic arguments
    va_list args1;
    va_start(args1, Format);
    va_list args2;
    va_copy(args2, args1);

    // Use vsnprintf to determine the required buffer size
    int size = vsnprintf(nullptr, 0, Format, args1);
    va_end(args1);

    std::string buffer(size, '\0');
    // Now format the string into our buffer
    vsnprintf(&buffer[0], size + 1, Format, args2);
    va_end(args2);

    py::set_error(PyExc_AssertionError, buffer.c_str());
    throw py::error_already_set();
}


// Python wrapper for PDBChain
class PyPDBChain {
private:
    PDBChain m_chain;

public:
    PyPDBChain() = default;
    
    PyPDBChain(const PDBChain& chain) : m_chain(chain) {
    }
    
    std::string __str__() const {
        return m_chain.m_Label;
    }

    py::str __repr__() const {
        return py::str("<PDBChain \"{}\">").attr("format")(m_chain.m_Label);
    }
    
    py::memoryview GetSequence() {
        const std::string& seq = m_chain.m_Seq;
        // py::memoryview::from_memory creates a zero-copy view. We must mark
        // the function with py::return_value_policy::reference_internal to tie
        // the reference to PyPDBChain.
        return py::memoryview::from_memory((void*)seq.data(), seq.size(), /*readonly=*/true);

    }
    
    // Get the underlying PDBChain (for internal use by other wrapper classes)
    const PDBChain& GetChain() const {
        return m_chain;
    }
    
    
    // Python __len__ method
    int __len__() const {
        return m_chain.GetSeqLength();
    }
};

// Python wrapper for DSS
class PyDSS {
private:
    DSS m_dss;
    DSSParams m_params;

public:
    PyDSS() {
        // Initialize with default parameters
        m_params.SetDefaults();
        m_dss.SetParams(m_params);
    }
    
    py::list get_profiles(const PyPDBChain& chain) {
        // Initialize DSS with the provided chain
        const PDBChain& inner = chain.GetChain();
        m_dss.Init(inner);
    
        const uint L = inner.GetSeqLength();
        const uint FeatureCount = m_params.GetFeatureCount();
    
        // Create the Python list that we will return
        py::list result_list{FeatureCount};
    
        for (uint i = 0; i < FeatureCount; ++i) {
            // For each profile, create a numpy array
            auto profile_row_arr = py::array_t<byte>(L);
    
            // Get a direct, mutable pointer to the numpy array's buffer
            byte* ptr = profile_row_arr.mutable_data();
            FEATURE Feature = m_params.m_Features[i];
            for (uint Pos = 0; Pos < L; ++Pos) {
                uint Letter = m_dss.GetFeature(Feature, Pos);
                if (Letter == UINT_MAX) {
                    ptr[Pos] = 31;
                } else {
                    asserta(Letter < 31);
                    ptr[Pos] = byte(Letter);
                }
            }
            result_list[i] = profile_row_arr;
        }
    
        return result_list;
    }
    
    py::list get_alphabets() {
        const uint FeatureCount = m_params.GetFeatureCount();
        
        // Create the Python list that we will return
        py::list result_list{FeatureCount};
        
        for (uint i = 0; i < FeatureCount; ++i) {
            FEATURE Feature = m_params.m_Features[i];
            uint size = DSS::GetAlphaSize(Feature);
            const char* name = FeatureToStr(Feature);
            
            // Create a dictionary with the size and name
            py::dict alphabet_dict;
            alphabet_dict["size"] = size;
            alphabet_dict["name"] = std::string(name);
            result_list[i] = alphabet_dict;
        }
        
        return result_list;
    }
};

// Python wrapper for BCAData
class PyBCAFile {
private:
    BCAData m_bca;
    std::string m_filename;
    bool m_opened = false;

public:
    PyBCAFile() = default;
    
    void open(const std::string& filename) {
        if (m_opened) {
            m_bca.Close();
        }
        
        try {
            m_bca.Open(filename);
            m_filename = filename;
            m_opened = true;
        } catch (const std::exception& e) {
            throw std::runtime_error("Cannot open BCA file: " + std::string(e.what()));
        }
    }
    
    void close() {
        if (m_opened) {
            m_bca.Close();
            m_opened = false;
        }
    }
    
    // Get chain count
    uint get_chain_count() const {
        if (!m_opened) {
            throw std::runtime_error("BCA file not opened");
        }
        return m_bca.GetChainCount();
    }
    
    // Get chain by index
    PyPDBChain get_chain(uint index) {
        if (!m_opened) {
            throw std::runtime_error("BCA file not opened");
        }
        
        uint chain_count = m_bca.GetChainCount();
        if (index >= chain_count) {
            throw std::runtime_error("Chain index out of range");
        }
        
        PDBChain chain;
        m_bca.ReadChain(index, chain);
        return PyPDBChain(chain);
    }
    
    // Python __len__ method
    int __len__() const {
        return get_chain_count();
    }
    
    // Python __getitem__ method
    PyPDBChain __getitem__(int index) {
        if (index < 0) {
            index += get_chain_count();
        }
        return get_chain(index);
    }
    
    ~PyBCAFile() {
        close();
    }
};

PYBIND11_MODULE(bca_reader, m) {
    m.doc() = "Python bindings for BCA file reader";
    
    // Expose PDBChain wrapper
    py::class_<PyPDBChain>(m, "PDBChain")
        .def(py::init<>())
        .def("__len__", &PyPDBChain::__len__)
        .def("__str__", &PyPDBChain::__str__)
        .def("__repr__", &PyPDBChain::__repr__)
        .def("get_sequence",
            &PyPDBChain::GetSequence,
            py::return_value_policy::reference_internal, // ties the lifetime of the returned memoryview
            "Returns a zero-copy memoryview of the sequence.");
    
    // Expose DSS wrapper
    py::class_<PyDSS>(m, "DSS")
        .def(py::init<>())
        .def("get_profiles", &PyDSS::get_profiles, "Get profiles for a PDBChain")
        .def("get_alphabets", &PyDSS::get_alphabets, "Get alphabet sizes for each feature");
    
    // Expose BCAFile wrapper
    py::class_<PyBCAFile>(m, "BCAFile")
        .def(py::init<>())
        .def("open", &PyBCAFile::open)
        .def("close", &PyBCAFile::close)
        .def("get_chain_count", &PyBCAFile::get_chain_count)
        .def("get_chain", &PyBCAFile::get_chain)
        .def("__len__", &PyBCAFile::__len__)
        .def("__getitem__", &PyBCAFile::__getitem__);
}
