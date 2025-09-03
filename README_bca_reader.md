# BCA Reader - Standalone Program

This is a standalone C++ program that reads `.bca` (Binary C-alpha) PDB chain database files and displays basic information about them.

## What it does

The program loads a `.bca` file and prints:
- The total number of protein chains/sequences
- Sequence lengths for the first few chains
- Basic validation of the file format

## Compilation

```bash
make
```

This will create the `bca_reader` executable.

## Usage

```bash
./bca_reader <filename.bca>
```

## Example Output

```
BCA file loaded successfully!
Number of sequences: 1234

First few sequence lengths:
  Chain 0: 156 residues
  Chain 1: 89 residues
  Chain 2: 234 residues
  Chain 3: 67 residues
  Chain 4: 189 residues
  ... and 1229 more chains
```

## Testing

To test the program with a sample file:

```bash
./test_bca_reader.sh
```

This script will:
1. Try to create a test `.bca` file from existing test data
2. Run the `bca_reader` program on it
3. Clean up the test files

## File Format

The `.bca` format is a binary format that stores:
- File header with magic number (0xBCABCA)
- Number of chains
- Sequence lengths for each chain
- Binary data for each chain (amino acid sequences + compressed coordinates)
- Chain labels

## Dependencies

- C++11 compatible compiler (g++, clang++)
- Standard C++ libraries only
- No external dependencies

## Comparison with Full Reseek

This standalone program provides just the basic `.bca` file reading functionality. The full `reseek` project provides:

- Multiple file format support (`.bca`, `.cal`, `.pdb`, `.cif`)
- Full protein chain data extraction
- Coordinate conversion and manipulation
- Thread-safe reading
- Integration with larger bioinformatics workflows

## Error Handling

The program includes error handling for:
- File not found
- Invalid file format (wrong magic number)
- Corrupted file headers
- File read errors

## Cleanup

```bash
make clean
```

Removes the compiled executable.
