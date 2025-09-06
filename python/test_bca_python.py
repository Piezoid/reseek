#!/usr/bin/env python3
"""
Test script for the BCA Python extension module.

This script tests the BCA file reader and DSS profile generation functionality.
The DSS wrapper is used to generate feature profiles from PDB chains.
"""

import sys
import os
import bca_reader

def test_bca_file(filename):
    """Test opening a BCA file and reading chains."""
    print(f"\nTesting BCA file: {filename}")
    

    # Create BCA file object
    bca = bca_reader.BCAFile()
    
    # Open the file
    bca.open(filename)
    
    # Get chain count
    chain_count = len(bca)
    print(f"Number of chains: {chain_count}")
    
    if chain_count > 0:
        # Get first chain
        chain = bca[0]
        print(f"First chain label: {chain!s}")
        print(f"First chain sequence length: {len(chain)}")
        print(f"First chain sequence: {bytes(chain.get_sequence()[:50])}...")
        
        # Get profiles using DSS wrapper (new API)
        dss = bca_reader.DSS()
        profiles = dss.get_profiles(chain)
        
        # Get alphabet information
        alphabets = dss.get_alphabets()
        
        for profile, alphabet in zip(profiles, alphabets):
            print(f"  Profile ({alphabet['name']}): {len(profile)} bytes, alphabet size: {alphabet['size']}")
            if len(profile) > 0:
                print(f"    First few values: {profile[:10]}")
    
    # Close the file
    bca.close()
    print("File closed successfully")


def main():
    """Main test function."""
    print("BCA Python Extension Module Test")
    print("=" * 40)
    
    # Check if test file exists
    test_file = "test.bca"
    if not os.path.exists(test_file):
        print(f"Test file '{test_file}' not found.")
        print("Please provide a .bca file path as an argument.")
        if len(sys.argv) > 1:
            test_file = sys.argv[1]
        else:
            print("Usage: python3 test_bca_python.py <filename.bca>")
            sys.exit(1)
    
    # Test the module
    test_bca_file(test_file)


if __name__ == "__main__":
    main()
