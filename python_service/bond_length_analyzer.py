"""
Bond Length Analyzer for CarsimMD

This script calculates the correct bond lengths of molecules based on 
physical principles without modifying the visualization. It serves as a
diagnostic tool for the molecular visualization engine.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import json
import numpy as np
import logging
from typing import Dict, List, Optional, Union, Any

# Import our custom modules
from bond_physics import BondPhysicsCalculator, calculate_molecule_bond_physics
from molecular_dynamics import create_simulation_from_mol_data

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def calculate_ideal_bond_lengths(mol_data: str, data_format: str = 'pdb') -> Dict[str, Any]:
    """
    Calculate ideal bond lengths for a molecule without changing the visualization.
    
    Args:
        mol_data: Molecular data (SMILES, PDB, etc.)
        data_format: Format of the molecular data ('smiles', 'pdb', etc.)
        
    Returns:
        Dictionary containing bond length analysis and statistics
    """
    try:
        # Use the existing bond physics calculator for basic analysis
        physics_data = calculate_molecule_bond_physics(mol_data, data_format)
        
        # Check if there was an error
        if "error" in physics_data:
            logger.error(f"Error in bond physics calculation: {physics_data['error']}")
            return physics_data
        
        # Get bond data
        bond_data = physics_data.get("bond_data", [])
        
        # Calculate statistics
        if bond_data:
            # Extract adjustments
            adjustments = [bond["adjustment_needed"] for bond in bond_data]
            
            # Calculate statistics
            avg_adjustment = sum(adjustments) / len(adjustments)
            max_adjustment = max(adjustments, key=abs)
            min_adjustment = min(adjustments, key=abs)
            
            # Identify bonds with significant adjustments needed
            significant_adjustments = [
                bond for bond in bond_data 
                if abs(bond["adjustment_needed"]) > 0.1  # Threshold for significance (Angstroms)
            ]
            
            # Calculate RMS deviation
            rms_deviation = np.sqrt(sum(adj**2 for adj in adjustments) / len(adjustments))
            
            # Add statistics to results
            physics_data["statistics"] = {
                "average_adjustment": float(avg_adjustment),
                "max_adjustment": float(max_adjustment),
                "min_adjustment": float(min_adjustment),
                "rms_deviation": float(rms_deviation),
                "significant_adjustment_count": len(significant_adjustments),
                "total_bond_count": len(bond_data)
            }
            
            # Add list of significant adjustments
            physics_data["significant_adjustments"] = significant_adjustments
            
        return physics_data
        
    except Exception as e:
        logger.exception(f"Error calculating ideal bond lengths: {str(e)}")
        return {"error": str(e)}

def analyze_molecule_using_simulator(mol_data: str, data_format: str = 'pdb') -> Dict[str, Any]:
    """
    Use the molecular dynamics simulator to analyze bond lengths.
    
    This approach uses the more comprehensive molecular dynamics infrastructure
    to calculate bond properties.
    
    Args:
        mol_data: Molecular data (SMILES, PDB, etc.)
        data_format: Format of the molecular data ('smiles', 'pdb', etc.)
        
    Returns:
        Dictionary containing bond length analysis and statistics
    """
    try:
        # Create a simulation
        simulation = create_simulation_from_mol_data(mol_data, data_format)
        
        # Get bond data from the simulation
        bond_data = simulation.get_bond_data()
        
        # Calculate statistics
        if bond_data:
            # Extract differences (equilibrium - current)
            differences = [bond["difference"] for bond in bond_data]
            
            # Calculate statistics
            avg_difference = sum(differences) / len(differences)
            max_difference = max(differences, key=abs)
            min_difference = min(differences, key=abs)
            
            # Identify bonds with significant differences
            significant_differences = [
                bond for bond in bond_data 
                if abs(bond["difference"]) > 0.1  # Threshold for significance (Angstroms)
            ]
            
            # Calculate RMS deviation
            rms_deviation = np.sqrt(sum(diff**2 for diff in differences) / len(differences))
            
            # Return results
            return {
                "bond_data": bond_data,
                "statistics": {
                    "average_difference": float(avg_difference),
                    "max_difference": float(max_difference),
                    "min_difference": float(min_difference),
                    "rms_deviation": float(rms_deviation),
                    "significant_difference_count": len(significant_differences),
                    "total_bond_count": len(bond_data)
                },
                "significant_differences": significant_differences,
                "atom_count": simulation.num_atoms,
                "molecule_energy": float(simulation.force_field.calculate_energy(
                    simulation.mol, simulation.positions))
            }
            
        return {"error": "No bond data available", "bond_data": []}
        
    except Exception as e:
        logger.exception(f"Error analyzing molecule with simulator: {str(e)}")
        return {"error": str(e)}

def get_complete_molecule_analysis(mol_data: str, data_format: str = 'pdb') -> Dict[str, Any]:
    """
    Provide a complete analysis of molecule bond lengths using both methods.
    
    Args:
        mol_data: Molecular data (SMILES, PDB, etc.)
        data_format: Format of the molecular data ('smiles', 'pdb', etc.)
        
    Returns:
        Dictionary containing comprehensive bond length analysis
    """
    # Get basic bond physics analysis
    basic_analysis = calculate_ideal_bond_lengths(mol_data, data_format)
    
    # Get simulator-based analysis
    simulation_analysis = analyze_molecule_using_simulator(mol_data, data_format)
    
    # Combine results
    return {
        "basic_analysis": basic_analysis,
        "simulation_analysis": simulation_analysis,
        "data_format": data_format,
        "has_errors": "error" in basic_analysis or "error" in simulation_analysis
    }

# Main function for command-line usage
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python bond_length_analyzer.py <smiles_or_pdb_string> [format]")
        sys.exit(1)
        
    mol_data = sys.argv[1]
    data_format = sys.argv[2] if len(sys.argv) > 2 else "smiles"
    
    # Analyze molecule
    results = get_complete_molecule_analysis(mol_data, data_format)
    
    # Print results
    print(json.dumps(results, indent=2)) 