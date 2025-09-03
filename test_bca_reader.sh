#!/bin/bash

echo "Testing BCA Reader Program"
echo "=========================="

# Check if reseek tool exists
if [ ! -f "./bin/reseek" ]; then
    echo "Error: reseek tool not found at ./bin/reseek"
    echo "Please build the project first or run from the correct directory"
    exit 1
fi

# Check if bca_reader program exists
if [ ! -f "./bca_reader" ]; then
    echo "Error: bca_reader program not found"
    echo "Please compile it first with: make"
    exit 1
fi

echo "1. Creating a test .bca file from test data..."

# Create a small test .bca file from the existing test data
# First, let's try to convert a PDB file to .bca format
if [ -f "./test_structures/PDB_1hhs.pdb.gz" ]; then
    echo "   Converting PDB_1hhs.pdb.gz to test.bca..."
    ./bin/reseek -convert ./test_structures/PDB_1hhs.pdb.gz -bca test.bca
    
    if [ $? -eq 0 ] && [ -f "test.bca" ]; then
        echo "   Successfully created test.bca"
        
        echo ""
        echo "2. Testing bca_reader program..."
        echo "   Running: ./bca_reader test.bca"
        echo ""
        
        ./bca_reader test.bca
        
        if [ $? -eq 0 ]; then
            echo ""
            echo "✅ Test passed! bca_reader successfully read the .bca file"
        else
            echo ""
            echo "❌ Test failed! bca_reader encountered an error"
        fi
        
        # Clean up
        echo ""
        echo "3. Cleaning up test files..."
        rm -f test.bca
        echo "   Removed test.bca"
        
    else
        echo "   Failed to create test.bca file"
        echo "   This might be because the reseek tool doesn't support -convert with -bca"
        echo "   or there's an issue with the input file format"
        exit 1
    fi
    
else
    echo "   test_structures/PDB_1hhs.pdb.gz not found, skipping test"
    echo "   You can manually test with any .bca file using:"
    echo "   ./bca_reader <your_file.bca>"
fi

echo ""
echo "Test completed!"
