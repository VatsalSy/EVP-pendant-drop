# EVP-pendant-drop

This repository has moved to [comphy-lab/Pendant-Drop-EVP](https://github.com/comphy-lab/Pendant-Drop-EVP).

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![OpenMP](https://img.shields.io/badge/OpenMP-enabled-brightgreen.svg)](https://www.openmp.org/)
[![MPI](https://img.shields.io/badge/MPI-enabled-brightgreen.svg)](https://www.open-mpi.org/)

A high-performance computational framework for studying viscoelastic pendant drops using [Basilisk C](http://basilisk.fr/). This repository provides a complete suite of simulation tools for investigating two-phase viscoelastic flows with adaptive mesh refinement, focusing on the dynamics of pendant drops and die-swell phenomena.

## Key Features

### **Advanced Viscoelastic Flow Modeling**
- Log-conformation solver for viscoelastic fluids using the Oldroyd-B model
- Multiple implementation approaches:
  - Scalar-based approach for both 2D and 3D
  - Traditional tensor-based approach for 2D
- Two-phase flow capabilities with density/viscosity/elastic property interpolation
- Volume of Fluid (VOF) method for sharp interface tracking
- Accurate surface tension implementation using the Brackbill method

### **High-Performance Computing**
- Adaptive mesh refinement based on multiple criteria
- Parallel computation support:
  - OpenMP for shared-memory systems
  - MPI for distributed computing

### **Specialized Numerical Methods**
- Eigen-decomposition for 3x3 symmetric matrices
- Householder transformations for matrix reduction
- QL algorithm for eigenvalue computation
- Axisymmetric computations support

## Project Structure

```plaintext
├── basilisk/src/               # Core Basilisk C framework
│   ├── navier-stokes/         # Flow solvers
│   │   ├── centered.h         # Main centered NS solver
│   │   └── conserving.h       # Conservative form solver
│   ├── reduced.h              # Reduced gravity approach
│   ├── curvature.h           # Interface property calculations
│   ├── tension.h             # Surface tension implementation
│   ├── tracer.h              # Passive tracer advection
│   ├── diffusion.h           # Diffusion equation solver
│   ├── vof.h                 # Volume of Fluid method
│   ├── viscosity.h           # Implicit viscous stress solver
│   └── axi.h                 # Axisymmetric computations
├── src-local/                 # Project-specific implementations
│   ├── eigen_decomposition.h  # 3x3 matrix eigenvalue solver
│   ├── log-conform-viscoelastic-scalar-2D.h  # 2D scalar log-conformation
│   ├── log-conform-viscoelastic-scalar-3D.h  # 3D scalar log-conformation
│   ├── log-conform-viscoelastic.h  # Traditional log-conformation
│   └── two-phaseVE.h         # Two-phase viscoelastic solver
├── postProcess/               # Analysis and visualization tools
└── testCases/                # Validation cases
```

## Installation

### Prerequisites
- GCC compiler (version 7.0+)
- MPI implementation (OpenMPI or MPICH)
- Basilisk C installation
- Xcode Command Line Tools (for MacOS users)

### Setting up Basilisk
1. Clone this repository:
   ```bash
   git clone [repository-url]
   cd EVP-pendant-drop
   ```

2. Build the Basilisk tree shipped in `basilisk/` (no network access, no viewer):
   ```bash
   bash reset_install_requirements.sh
   source .project_config
   ```
3. Add the src-local directory to your Basilisk include path

## Usage

### Running Simulations

#### Local Development (OpenMP)
```bash
# Compile with OpenMP support
qcc -O2 -Wall -disable-dimensions -fopenmp -I$(PWD)/src-local your_case.c -o simulation -lm

# Run with desired number of threads
export OMP_NUM_THREADS=4
./simulation
```

#### HPC Deployment (MPI)
```bash
# Compile with MPI support
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions \
-I$(PWD)/src-local your_case.c -o simulation -lm

# Run with MPI
mpirun -np 4 ./simulation
```

## Contributing

We welcome contributions! Please feel free to submit pull requests or create issues using our templates:

- [🐛 Report a Bug](../../issues/new?template=bug_report.md&labels=bug)
- [✨ Request a Feature](../../issues/new?template=feature_request.md&labels=enhancement)
- [📚 Improve Documentation](../../issues/new?template=documentation.md&labels=documentation)
- [⚡ Suggest Performance Enhancement](../../issues/new?template=performance.md&labels=performance)
- [📝 Open a Blank Issue](../../issues/new?template=blank.md)

When contributing, please:
1. Fork the repository and create your branch from `main`
2. Add tests for any new features
3. Ensure your code follows our coding standards
4. Update documentation as needed
5. Submit a pull request with a comprehensive description of changes

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite our work (citation details to be added).

## Acknowledgments
- Based on the [Basilisk C](http://basilisk.fr/) framework
- Contributors and researchers who have helped develop and test the code
