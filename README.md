# Die-Swell EVP Pendant Drop

## Usage

### Running Test Cases

*Note:* conda and mpicc can cause some issues. Deactivate conda environment before running the test cases.

To compile and run test cases, use the `runCodesInParallel.sh` script in the `testCases` directory. The script takes a filename (without the .c extension) as a required argument and an optional number of processes for MPI:

```bash
cd testCases
bash runCodesInParallel.sh <filename> [number_of_processes]
```

For example:
```bash
# Run with default 4 processes
bash runCodesInParallel.sh die-swell_ViscoElastic

# Run with 8 processes
bash runCodesInParallel.sh die-swell_ViscoElastic 8
```

This will:
1. Create a directory with the same name as the input file
2. Compile the corresponding .c file
3. Place the executable in the newly created directory
4. Run the executable using MPI with the specified number of processes (default: 4)

Available test cases:
- die-swell_ViscoElastic.c
- die-swell_Newt.c
