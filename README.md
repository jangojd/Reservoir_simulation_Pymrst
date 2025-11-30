# PyMRST: Python Reservoir Simulation Toolbox

A comprehensive Python implementation of reservoir simulation based on MATLAB Reservoir Simulation Toolbox (MRST) principles. This project provides tools for modeling single-phase and multi-phase fluid flow in porous media.

## Result figure
![Alt text](https://github.com/jangojd/Reservoir_simulation_Pymrst/blob/b7559b247bdaed7f5cc22bc9a55fc51794ed9aff/download.png)


## 🎯 Overview

PyMRST is an open-source Python library designed to simulate fluid flow in oil and gas reservoirs. It implements finite-volume numerical methods to solve pressure and saturation equations in complex geological formations. The toolbox is ideal for students, researchers, and professionals working in petroleum engineering, hydrogeology, and energy research.

### Key Features

- **Grid Generation**: Flexible cartesian grid creation with customizable cell dimensions
- **Rock Properties**: Realistic porosity and permeability distributions
- **Fluid Modeling**: Single-phase and multi-phase flow support
- **Well Management**: Producer and injector well definitions with various control strategies
- **Pressure Solver**: Implicit finite-volume pressure equation solver
- **Visualization**: Comprehensive 2D/3D plotting of results
- **History Matching**: Built-in tools for pressure and production data analysis
- **Performance Optimization**: Sparse matrix operations for efficient computation

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Usage Examples](#usage-examples)
- [API Documentation](#api-documentation)
- [Theory](#theory)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)



## 🚀 Quick Start

### Basic Example: Single-Phase Oil Depletion

```python
from pymrst_simulator import CartesianGrid, RockProperties, FluidProperties, run_simulation, plot_results

# Create a simple 2D grid
grid = CartesianGrid(nx=10, ny=10, nz=1, Lx=1000, Ly=1000, Lz=100)

# Define properties
rock = RockProperties(grid.num_cells, porosity=0.2, perm_mean=100)
fluid = FluidProperties()

# Run simulation
results = run_simulation()

# Visualize
plot_results(results)
```

## 📁 Project Structure

```
PyMRST/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── LICENSE                      # MIT License
├── .gitignore                   # Git ignore patterns
│
├── pymrst_simulator.py          # Main simulation engine
├── grid.py                      # Grid generation module
├── properties.py                # Rock and fluid properties
├── wells.py                     # Well management
├── solvers.py                   # Pressure/saturation solvers
├── visualization.py             # Plotting utilities
├── utils.py                     # Helper functions
│
├── examples/
│   ├── simple_depletion.py     # Single-phase oil depletion
│   ├── waterflooding.py        # Two-phase flow example
│   ├── multi_well_producer.py  # Multiple well simulation
│   ├── history_matching.py     # History matching workflow
│   └── sensitivity_analysis.py # Parameter sensitivity study
│
├── tests/
│   ├── test_grid.py            # Grid tests
│   ├── test_properties.py      # Properties tests
│   ├── test_solvers.py         # Solver tests
│   └── test_integration.py     # Integration tests
│
└── docs/
    ├── theory.md               # Mathematical theory
    ├── api_reference.md        # API documentation
    ├── tutorials.md            # Detailed tutorials
    └── troubleshooting.md      # FAQ and troubleshooting
```

## 📖 Usage Examples

### Example 1: Simple Reservoir Depletion

```python
from pymrst_simulator import run_simulation, plot_results

# Run default 1000-day simulation
results = run_simulation()

# Access results
print(f"Initial Pressure: {results['pressure_avg'][0]:.1f} bar")
print(f"Final Pressure: {results['pressure_avg'][-1]:.1f} bar")
print(f"Total Production: {results['well_flux'][-1]:.3e} m³/s")

# Plot all results
plot_results(results)
```

### Example 2: Custom Grid and Properties

```python
from pymrst_simulator import CartesianGrid, RockProperties, FluidProperties

# Create larger grid
grid = CartesianGrid(nx=20, ny=20, nz=5, Lx=5000, Ly=5000, Lz=500)

# Custom rock properties
rock = RockProperties(
    num_cells=grid.num_cells,
    porosity=0.25,
    perm_mean=150,
    perm_var=50
)

# Custom fluid properties
fluid = FluidProperties()
fluid.mu = 2e-3  # 2 cp oil
fluid.rho = 850  # kg/m³
```

### Example 3: Multiple Wells

```python
from pymrst_simulator import Well, compute_transmissibilities, assemble_pressure_matrix

wells = [
    Well("Producer_1", cell_idx=50, well_type='producer', target_pressure=50e5),
    Well("Producer_2", cell_idx=150, well_type='producer', target_pressure=50e5),
    Well("Injector_1", cell_idx=10, well_type='injector', target_pressure=150e5)
]

# Transmissibilities with multiple wells
trans_matrix = compute_transmissibilities(grid, rock, fluid)
A, b = assemble_pressure_matrix(grid, rock, trans_matrix, fluid, wells)
```

See `examples/` directory for more detailed examples.

## 📚 API Documentation

### Core Classes

#### `CartesianGrid`

Creates a structured cartesian grid for the reservoir domain.

```python
grid = CartesianGrid(nx=10, ny=10, nz=1, Lx=1000, Ly=1000, Lz=100)

# Properties
grid.num_cells      # Total number of cells
grid.cell_vol       # Cell volumes (array)
grid.cell_centers   # Cell center coordinates
grid.dx, grid.dy, grid.dz  # Cell dimensions
```

#### `RockProperties`

Defines rock characteristics (porosity, permeability).

```python
rock = RockProperties(num_cells=100, porosity=0.2, perm_mean=100, perm_var=20)

# Properties
rock.porosity       # Porosity values (fraction)
rock.permeability   # Permeability values (m²)
```

#### `FluidProperties`

Defines fluid characteristics (viscosity, density, compressibility).

```python
fluid = FluidProperties()

# Properties
fluid.mu   # Viscosity (Pa·s)
fluid.rho  # Density (kg/m³)
fluid.c    # Compressibility (1/Pa)
```

#### `Well`

Represents a production or injection well.

```python
well = Well(
    name="Producer_1",
    cell_idx=50,
    well_type='producer',  # or 'injector'
    target_pressure=50e5
)

# Properties
well.name
well.cell_idx
well.type
well.target_pressure
well.productivity_index
```

### Core Functions

#### `compute_transmissibilities(grid, rock, fluid)`

Calculates flow transmissibilities between adjacent grid cells.

**Returns:** Sparse transmissibility matrix

#### `assemble_pressure_matrix(grid, rock, trans_matrix, fluid, wells)`

Assembles the pressure equation system matrix and RHS vector.

**Returns:** Tuple (A, b) - system matrix and RHS vector

#### `solve_pressure(A, b, pressure_prev, dt, grid, rock, fluid)`

Solves the implicit pressure equation.

**Returns:** Pressure array

#### `run_simulation(num_steps=50, dt_days=10, total_days=1000)`

Runs complete reservoir simulation with default parameters.

**Returns:** Dictionary with simulation results

#### `plot_results(results)`

Generates comprehensive visualization of simulation results.

**Plots generated:**
- Pressure depletion history
- Production rate vs time
- Cumulative production
- Final pressure distribution
- Permeability distribution
- Simulation summary statistics

## 🧮 Theory

### Mathematical Foundation

The simulator solves the following equations:

**Pressure Equation (Implicit):**
```
∇·(k/μ ∇p) = φc ∂p/∂t + q
```

Where:
- k = permeability (m²)
- μ = viscosity (Pa·s)
- p = pressure (Pa)
- φ = porosity (fraction)
- c = compressibility (1/Pa)
- q = source/sink term (wells)

**Discretization:** Finite Volume Method (FVM)
**Linear Solver:** Sparse matrix solver (SciPy)

### Transmissibility

```
T_ij = (k_avg * A) / (d * μ)
```

Where:
- k_avg = harmonic average permeability between cells
- A = interface area (m²)
- d = distance between cell centers (m)

## 🤝 Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Development Setup

```bash
# Clone your fork
git clone https://github.com/yourusername/PyMRST.git
cd PyMRST

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install in development mode
pip install -e .
pip install -r requirements-dev.txt

# Run tests
pytest tests/

# Run linting
flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 References

- Lie, K. A. (2019). "An Introduction to Reservoir Simulation Using MATLAB/Octave". SINTEF.
- Aziz, K., & Settari, A. (1979). "Petroleum Reservoir Simulation". Applied Science Publishers.
- Craft, B. C., & Hawkins, M. F. (2015). "Applied Petroleum Reservoir Engineering". Prentice Hall.
- MRST Documentation: https://www.sintef.no/projectweb/mrst/



## 🙏 Acknowledgments

- Inspired by MATLAB Reservoir Simulation Toolbox (MRST)
- Built for educational and research purposes
- Thanks to the open-source community

---

**Last Updated:** December 2024
**Version:** 1.0.0
**Status:** Active Development
