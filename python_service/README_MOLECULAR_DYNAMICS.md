# Molecular Dynamics Infrastructure for CarsimMD

This folder contains the molecular dynamics and simulation infrastructure for the CarsimMD application. The infrastructure is designed to enable physics-based simulations of molecular systems while keeping the UI completely separate.

## Components

1. **molecular_dynamics.py** - Core simulation framework
   - `MolecularForceField` - Base class for force field calculations
   - `SimpleMolecularForceField` - Implementation focusing on bond lengths
   - `MolecularDynamicsSimulation` - Main simulation engine using Velocity Verlet integration
   - Factory functions for creating simulations from different molecular representations

2. **bond_physics.py** - Bond length calculation and physics
   - `BondPhysicsCalculator` - Calculates accurate bond lengths based on physical properties
   - Functions to analyze bond geometries without changing the visualization

3. **bond_length_analyzer.py** - Analysis of bond lengths
   - Calculates ideal bond lengths for molecules
   - Provides statistics and analysis of bond length deviations
   - Does not modify molecule visualization

## Key Features

- **Physics-Based Force Fields**: Calculates forces between atoms based on physical principles
- **Molecular Dynamics**: Simulates molecule motion in response to forces
- **Bond Length Optimization**: Calculates accurate bond lengths based on:
  - Electronegativity differences
  - Bond order (single, double, triple)
  - Hybridization state
  - Resonance effects
  - Ring strain
- **Extensible Design**: Framework designed for adding more physical effects

## Future Extensions

The current infrastructure focuses on bond lengths but is designed to be extended for:

1. **Solvent Simulation**: Modeling molecule movement in water
2. **Torsional Dynamics**: Simulating bond rotations due to steric hindrance
3. **Electrostatic Interactions**: Effects from partial charges
4. **Energy Minimization**: Finding optimal molecular conformations

## API Endpoints

The following endpoints have been added to the Flask API:

- `/analyze-bond-lengths` - Analyzes correct bond lengths without changing visualization
- `/simulate-molecule` - Creates a new molecular dynamics simulation
- `/simulation/<id>/run` - Runs a simulation for a number of steps
- `/simulation/<id>/status` - Gets current simulation status
- `/simulation/<id>/stop` - Stops a running simulation
- `/simulation/<id>/bond-data` - Gets detailed bond data
- `/simulation/<id>/trajectory` - Gets trajectory data (if stored)
- `/simulate-and-analyze` - Convenience endpoint for quick simulation
- `/force-field-components` - Gets energy components from force field

## Usage Examples

```python
# Example 1: Analyzing bond lengths
from bond_length_analyzer import calculate_ideal_bond_lengths

# Analyze a molecule from SMILES
results = calculate_ideal_bond_lengths("CC(=O)O", data_format="smiles")
print(f"Found {len(results['bond_data'])} bonds")
print(f"Average adjustment needed: {results['statistics']['average_adjustment']:.4f} Å")

# Example 2: Running a simulation
from molecular_dynamics import create_simulation_from_mol_data

# Create simulation from PDB data
simulation = create_simulation_from_mol_data(pdb_data, data_format="pdb")

# Initialize velocities from temperature
simulation.initialize_velocities(temperature=300)  # 300K

# Run simulation for 1000 steps
simulation.run_simulation(num_steps=1000, store_trajectory=True)

# Wait for completion
import time
while simulation.state.value == "running":
    time.sleep(0.1)
    
# Get results
final_positions = simulation.positions
energies = simulation.energies
bond_data = simulation.get_bond_data()

# Print statistics
differences = [bond["difference"] for bond in bond_data]
rms_deviation = (sum(d**2 for d in differences) / len(differences)) ** 0.5
print(f"RMS deviation from ideal bond lengths: {rms_deviation:.4f} Å")
```

## Technical Notes

- All simulations run in separate threads to avoid blocking the main application
- Bond length calculations use established chemical principles and empirical data
- The Velocity Verlet integration scheme provides good energy conservation
- Energy is reported in kJ/mol, distances in Ångströms, time in picoseconds 