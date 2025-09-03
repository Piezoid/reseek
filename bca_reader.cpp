#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <vector>



// Include necessary reseek headers
#include "src/myutils.h"
#include "src/pdbchain.h"
#include "src/dss.h"
#include "src/dssparams.h"
#include "src/features.h"
#include "src/bcadata.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename.bca>" << std::endl;
        std::cerr << "  Loads a .bca file, extracts first chain, and prints profiles" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    
    // Check file extension
    if (filename.length() < 4 || filename.substr(filename.length() - 4) != ".bca") {
        std::cerr << "Warning: File doesn't have .bca extension" << std::endl;
    }
    
    // Use BCAData class to read the BCA file
    BCAData bca;
    try {
        bca.Open(filename);
    } catch (const std::exception& e) {
        std::cerr << "Error: Cannot open BCA file '" << filename << "': " << e.what() << std::endl;
        return 1;
    }
    
    uint64_t chainCount = bca.GetChainCount();
    if (chainCount == 0) {
        std::cerr << "Error: No chains found in BCA file" << std::endl;
        bca.Close();
        return 1;
    }
    
    // Read the first chain
    PDBChain chain;
    bca.ReadChain(0, chain);
    
    // Close the BCA file
    bca.Close();
    
    // Set up DSS parameters with default features
    DSSParams params;
    params.SetDefaults();
    
    // Create DSS object and generate profiles
    DSS dss;
    dss.SetParams(params);
    dss.Init(chain);

    // std::cout << "g_AlphaSizes2[FEATURE_AA] = " << g_AlphaSizes2[FEATURE_AA] << std::endl;
    // std::cout << "g_ScoreMxs2[FEATURE_AA] = " << (void*)g_ScoreMxs2[FEATURE_AA] << std::endl;
    
    std::vector<std::vector<unsigned char>> profiles;
    dss.GetProfile(profiles);
    
    // Display results
    std::cout << "BCA file loaded successfully!" << std::endl;
    std::cout << "Number of sequences: " << chainCount << std::endl;
    std::cout << "First chain length: " << chain.GetSeqLength() << " residues" << std::endl;
    std::cout << "First chain label: " << chain.m_Label << std::endl;
    
    // Print profiles as vectors of ints
    std::cout << "\nDSS Profiles (as vectors of ints):" << std::endl;
    std::cout << "Number of features: " << profiles.size() << std::endl;
    
    for (size_t i = 0; i < profiles.size(); i++) {
        std::cout << "Feature " << i << " (" << profiles[i].size() << " values): ";
        for (size_t j = 0; j < std::min(profiles[i].size(), size_t(20)); j++) {
            std::cout << (int)profiles[i][j] << " ";
        }
        if (profiles[i].size() > 20) {
            std::cout << "... and " << (profiles[i].size() - 20) << " more";
        }
        std::cout << std::endl;
    }
    
    return 0;
}
