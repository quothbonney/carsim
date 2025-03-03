"""
Molecular Dynamics Module for CarsimMD

This module provides a framework for molecular simulation, supporting:
1. Force field calculations
2. Dynamics simulations including:
   - Movement in solvent (water)
   - Bond torsion due to steric hindrance
   - Accurate bond lengths based on electronic properties

The module is designed to be extensible through plugin-style physics scripts.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import json
import time
import threading
import logging
from enum import Enum
from typing import Dict, List, Tuple, Callable, Optional, Union, Any

# Import bond physics
from bond_physics import BondPhysicsCalculator

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ForceType(Enum):
    """Types of forces that can be applied in simulations"""
    BOND_LENGTH = "bond_length"
    BOND_ANGLE = "bond_angle"
    TORSION = "torsion"
    VAN_DER_WAALS = "van_der_waals"
    ELECTROSTATIC = "electrostatic"
    SOLVENT = "solvent"
    CUSTOM = "custom"

class SimulationState(Enum):
    """Possible states of a simulation"""
    IDLE = "idle"
    RUNNING = "running"
    PAUSED = "paused"
    COMPLETED = "completed"
    ERROR = "error"

class MolecularForceField:
    """Base class for force field calculations"""
    
    def __init__(self, name: str = "generic"):
        self.name = name
        self.energy_components: Dict[ForceType, float] = {
            force_type: 0.0 for force_type in ForceType
        }
        
    def calculate_forces(self, mol: Chem.Mol, positions: np.ndarray) -> np.ndarray:
        """
        Calculate forces on all atoms based on current positions
        
        Args:
            mol: RDKit molecule
            positions: Nx3 array of atom positions
            
        Returns:
            Nx3 array of forces (in x, y, z directions)
        """
        num_atoms = mol.GetNumAtoms()
        forces = np.zeros((num_atoms, 3))
        
        # This is a base implementation that does nothing
        # Subclasses should override this method
        
        return forces
    
    def calculate_energy(self, mol: Chem.Mol, positions: np.ndarray) -> float:
        """
        Calculate total energy of the system
        
        Args:
            mol: RDKit molecule
            positions: Nx3 array of atom positions
            
        Returns:
            Total energy
        """
        # Reset energy components
        self.energy_components = {force_type: 0.0 for force_type in ForceType}
        
        # This is a base implementation that returns zero
        # Subclasses should override this method
        
        return sum(self.energy_components.values())

class SimpleMolecularForceField(MolecularForceField):
    """
    A simple molecular force field implementation focusing on bond lengths
    """
    
    def __init__(self):
        super().__init__(name="simple_force_field")
        self.bond_calculator = BondPhysicsCalculator()
        
    def calculate_bond_forces(self, mol: Chem.Mol, positions: np.ndarray) -> np.ndarray:
        """
        Calculate forces arising from bond length constraints
        
        Args:
            mol: RDKit molecule
            positions: Nx3 array of atom positions
            
        Returns:
            Nx3 array of forces from bond constraints
        """
        num_atoms = mol.GetNumAtoms()
        forces = np.zeros((num_atoms, 3))
        
        # Force constant for bonds (can be adjusted)
        k_bond = 100.0  # N/m or equivalent
        
        # Process each bond
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            
            # Calculate current bond vector
            pos1 = positions[idx1]
            pos2 = positions[idx2]
            bond_vector = pos2 - pos1
            current_length = np.linalg.norm(bond_vector)
            
            if current_length < 1e-6:  # Avoid division by zero
                continue
                
            # Get equilibrium bond length from bond physics
            atom1 = mol.GetAtomWithIdx(idx1)
            atom2 = mol.GetAtomWithIdx(idx2)
            
            # Use the bond calculator to get accurate length
            equilibrium_length = self.bond_calculator.calculate_bond_length(bond, mol)
            
            # Calculate force magnitude (Hooke's law: F = -k(x - x0))
            force_magnitude = -k_bond * (current_length - equilibrium_length)
            
            # Normalize bond vector to get direction
            bond_direction = bond_vector / current_length
            
            # Calculate force vector
            force_vector = force_magnitude * bond_direction
            
            # Apply forces (equal and opposite)
            forces[idx1] -= force_vector  # Force on first atom
            forces[idx2] += force_vector  # Equal and opposite force on second atom
            
            # Track energy
            self.energy_components[ForceType.BOND_LENGTH] += 0.5 * k_bond * (current_length - equilibrium_length)**2
            
        return forces
    
    def calculate_forces(self, mol: Chem.Mol, positions: np.ndarray) -> np.ndarray:
        """
        Calculate all forces on atoms
        
        Args:
            mol: RDKit molecule
            positions: Nx3 array of atom positions
            
        Returns:
            Nx3 array of net forces
        """
        forces = np.zeros((mol.GetNumAtoms(), 3))
        
        # Add bond forces
        bond_forces = self.calculate_bond_forces(mol, positions)
        forces += bond_forces
        
        # In the future, we can add:
        # - angle_forces = self.calculate_angle_forces(mol, positions)
        # - torsion_forces = self.calculate_torsion_forces(mol, positions)
        # - vdw_forces = self.calculate_vdw_forces(mol, positions)
        # - electrostatic_forces = self.calculate_electrostatic_forces(mol, positions)
        # - solvent_forces = self.calculate_solvent_forces(mol, positions)
        
        return forces
    
    def calculate_energy(self, mol: Chem.Mol, positions: np.ndarray) -> float:
        """Calculate total energy by summing all components"""
        # Reset energy components
        self.energy_components = {force_type: 0.0 for force_type in ForceType}
        
        # Calculate forces - this will populate energy_components
        self.calculate_forces(mol, positions)
        
        # Return total energy
        return sum(self.energy_components.values())

class MolecularDynamicsSimulation:
    """
    Class for running molecular dynamics simulations
    """
    
    def __init__(self, mol: Chem.Mol, force_field: Optional[MolecularForceField] = None):
        """
        Initialize a molecular dynamics simulation
        
        Args:
            mol: RDKit molecule with 3D coordinates
            force_field: Force field to use (or create a simple one if None)
        """
        self.mol = mol
        
        # Make sure molecule has 3D coordinates
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol)
            if not mol.GetNumConformers():
                raise ValueError("Could not generate 3D coordinates for molecule")
        
        # Initialize properties
        self.num_atoms = mol.GetNumAtoms()
        self.atom_masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
        
        # Convert atom masses to kg (simulation units)
        self.atom_masses_kg = self.atom_masses * 1.66053886e-27  # amu to kg
        
        # Copy initial positions from molecule
        conf = mol.GetConformer()
        self.positions = np.array([
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(self.num_atoms)
        ])
        
        # Initialize velocities to zero
        self.velocities = np.zeros((self.num_atoms, 3))
        
        # Set up force field
        self.force_field = force_field if force_field else SimpleMolecularForceField()
        
        # Simulation parameters
        self.time_step = 0.001  # ps
        self.temperature = 300.0  # K
        self.total_steps = 0
        self.current_step = 0
        
        # Simulation state
        self.state = SimulationState.IDLE
        self.simulation_thread = None
        self.stop_requested = False
        
        # Trajectory storage
        self.trajectory = []
        self.energies = []
        self.store_trajectory = False
        
    def initialize_velocities(self, temperature: float = 300.0):
        """
        Initialize random velocities based on temperature
        
        Args:
            temperature: Temperature in Kelvin
        """
        # Boltzmann constant in J/K
        k_B = 1.380649e-23
        
        # Calculate standard deviation for each component based on Maxwell-Boltzmann distribution
        # σ = sqrt(k_B * T / m)
        sigma = np.sqrt(k_B * temperature / self.atom_masses_kg.reshape(-1, 1))
        
        # Generate random velocities from normal distribution
        self.velocities = np.random.normal(0, sigma, (self.num_atoms, 3))
        
        # Remove center of mass motion
        total_momentum = np.sum(self.atom_masses_kg.reshape(-1, 1) * self.velocities, axis=0)
        com_velocity = total_momentum / np.sum(self.atom_masses_kg)
        self.velocities -= com_velocity
        
    def run_simulation(self, num_steps: int, store_trajectory: bool = True):
        """
        Run simulation for specified number of steps
        
        Args:
            num_steps: Number of time steps to simulate
            store_trajectory: Whether to store trajectory data
        """
        self.total_steps = num_steps
        self.current_step = 0
        self.store_trajectory = store_trajectory
        self.trajectory = []
        self.energies = []
        self.stop_requested = False
        self.state = SimulationState.RUNNING
        
        # Run in a separate thread to not block
        self.simulation_thread = threading.Thread(target=self._simulation_loop)
        self.simulation_thread.start()
        
    def _simulation_loop(self):
        """Internal simulation loop (runs in separate thread)"""
        try:
            # Store initial state if tracking trajectory
            if self.store_trajectory:
                self.trajectory.append(self.positions.copy())
            
            # Velocity Verlet integration
            for step in range(self.total_steps):
                if self.stop_requested:
                    break
                    
                # Calculate forces at current positions
                forces = self.force_field.calculate_forces(self.mol, self.positions)
                
                # Calculate accelerations (F = ma)
                accelerations = forces / self.atom_masses_kg.reshape(-1, 1)
                
                # Update velocities (half step)
                self.velocities += 0.5 * accelerations * self.time_step
                
                # Update positions
                self.positions += self.velocities * self.time_step
                
                # Recalculate forces with new positions
                forces = self.force_field.calculate_forces(self.mol, self.positions)
                
                # Recalculate accelerations
                accelerations = forces / self.atom_masses_kg.reshape(-1, 1)
                
                # Update velocities (second half step)
                self.velocities += 0.5 * accelerations * self.time_step
                
                # Calculate energy if tracking
                if self.store_trajectory:
                    energy = self.force_field.calculate_energy(self.mol, self.positions)
                    self.energies.append(energy)
                    self.trajectory.append(self.positions.copy())
                
                # Update step counter
                self.current_step = step + 1
                
            # Mark as completed
            self.state = SimulationState.COMPLETED
            
        except Exception as e:
            logger.error(f"Error in simulation: {str(e)}")
            self.state = SimulationState.ERROR
            
    def stop_simulation(self):
        """Request simulation to stop"""
        self.stop_requested = True
        if self.simulation_thread and self.simulation_thread.is_alive():
            self.simulation_thread.join()
        self.state = SimulationState.PAUSED
        
    def update_mol_coordinates(self):
        """
        Update the molecular conformer with current simulation positions
        """
        if not self.mol.GetNumConformers():
            return
            
        conf = self.mol.GetConformer()
        for i in range(self.num_atoms):
            x, y, z = self.positions[i]
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
            
    def get_current_state(self) -> dict:
        """
        Get the current state of the simulation as a dictionary
        """
        return {
            "state": self.state.value,
            "current_step": self.current_step,
            "total_steps": self.total_steps,
            "temperature": self.temperature,
            "energy": self.energies[-1] if self.energies else None,
            "positions": self.positions.tolist()
        }
        
    def get_bond_data(self) -> list:
        """
        Get current bond data including current and equilibrium lengths
        """
        if isinstance(self.force_field, SimpleMolecularForceField):
            bond_calculator = self.force_field.bond_calculator
        else:
            bond_calculator = BondPhysicsCalculator()
            
        results = []
        
        for bond in self.mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            
            # Get current bond length from positions
            pos1 = self.positions[idx1]
            pos2 = self.positions[idx2]
            current_length = np.linalg.norm(pos2 - pos1)
            
            # Calculate equilibrium bond length
            equilibrium_length = bond_calculator.calculate_bond_length(bond, self.mol)
            
            atom1 = self.mol.GetAtomWithIdx(idx1)
            atom2 = self.mol.GetAtomWithIdx(idx2)
            
            results.append({
                'atom1': {
                    'index': idx1,
                    'symbol': bond_calculator.get_element_symbol(atom1.GetAtomicNum()),
                    'position': pos1.tolist()
                },
                'atom2': {
                    'index': idx2,
                    'symbol': bond_calculator.get_element_symbol(atom2.GetAtomicNum()),
                    'position': pos2.tolist()
                },
                'bond_order': bond_calculator.get_bond_order(bond),
                'is_aromatic': bond_calculator.is_aromatic(bond),
                'is_conjugated': bond_calculator.is_conjugated(bond),
                'is_in_ring': bond_calculator.is_in_ring(bond),
                'current_length': float(current_length),
                'equilibrium_length': float(equilibrium_length),
                'difference': float(equilibrium_length - current_length)
            })
            
        return results


# Factory function to create a simulation from molecule data
def create_simulation_from_mol_data(mol_data: str, data_format: str = 'smiles') -> MolecularDynamicsSimulation:
    """
    Create a molecular dynamics simulation from molecular data
    
    Args:
        mol_data: Molecular data (SMILES, PDB, etc.)
        data_format: Format of the molecular data ('smiles', 'pdb', etc.)
        
    Returns:
        MolecularDynamicsSimulation object
    """
    mol = None
    
    # Create RDKit molecule based on format
    if data_format.lower() == 'smiles':
        mol = Chem.MolFromSmiles(mol_data)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol)
    elif data_format.lower() == 'pdb':
        mol = Chem.MolFromPDBBlock(mol_data)
        if mol:
            mol = Chem.AddHs(mol, addCoords=True)
    else:
        raise ValueError(f"Unsupported data format: {data_format}")
        
    if not mol:
        raise ValueError(f"Could not create molecule from {data_format} data")
        
    # Create simulation
    return MolecularDynamicsSimulation(mol) 