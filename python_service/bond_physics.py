"""
Bond Physics Module for CarsimMD

This module provides functions for calculating accurate bond lengths based on:
1. Electronegativity of atoms
2. Resonance effects
3. Bond order
4. Hybridization state

The calculations use standard chemistry principles and are based on both 
theoretical models and empirical adjustments.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, GetPeriodicTable
import numpy as np
import tempfile

# Pauling electronegativity values for common elements
ELECTRONEGATIVITY = {
    'H': 2.20, 'He': 0.00,
    'Li': 0.98, 'Be': 1.57, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98, 'Ne': 0.00,
    'Na': 0.93, 'Mg': 1.31, 'Al': 1.61, 'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Ar': 0.00,
    'K': 0.82, 'Ca': 1.00, 'Sc': 1.36, 'Ti': 1.54, 'V': 1.63, 'Cr': 1.66, 'Mn': 1.55, 'Fe': 1.83,
    'Co': 1.88, 'Ni': 1.91, 'Cu': 1.90, 'Zn': 1.65, 'Ga': 1.81, 'Ge': 2.01, 'As': 2.18, 'Se': 2.55,
    'Br': 2.96, 'Kr': 0.00, 'I': 2.66
}

# Standard bond lengths (in Angstroms) for common bonds
STANDARD_BOND_LENGTHS = {
    ('C', 'C', 1): 1.54,  # single bond
    ('C', 'C', 2): 1.34,  # double bond
    ('C', 'C', 3): 1.20,  # triple bond
    ('C', 'N', 1): 1.47,
    ('C', 'N', 2): 1.29,
    ('C', 'N', 3): 1.16,
    ('C', 'O', 1): 1.43,
    ('C', 'O', 2): 1.23,
    ('C', 'H', 1): 1.09,
    ('N', 'H', 1): 1.01,
    ('O', 'H', 1): 0.96,
    ('C', 'S', 1): 1.82,
    ('C', 'S', 2): 1.60,
    ('C', 'F', 1): 1.35,
    ('C', 'Cl', 1): 1.77,
    ('C', 'Br', 1): 1.94,
    ('C', 'I', 1): 2.14,
}

class BondPhysicsCalculator:
    """Class for calculating accurate bond lengths based on physical properties"""
    
    def __init__(self):
        self.periodic_table = GetPeriodicTable()
        
    def get_element_symbol(self, atomic_num):
        """Get the element symbol from atomic number"""
        return self.periodic_table.GetElementSymbol(atomic_num)
    
    def get_electronegativity(self, atom):
        """Get electronegativity value for an atom"""
        symbol = self.get_element_symbol(atom.GetAtomicNum())
        return ELECTRONEGATIVITY.get(symbol, 2.0)  # Default to 2.0 if not found
    
    def get_hybridization_factor(self, atom):
        """Calculate a factor based on atom hybridization"""
        hyb = atom.GetHybridization()
        if hyb == Chem.rdchem.HybridizationType.SP3:
            return 1.0
        elif hyb == Chem.rdchem.HybridizationType.SP2:
            return 0.95
        elif hyb == Chem.rdchem.HybridizationType.SP:
            return 0.90
        else:
            return 1.0
    
    def is_in_ring(self, bond):
        """Check if bond is part of a ring"""
        return bond.IsInRing()
    
    def is_aromatic(self, bond):
        """Check if bond is aromatic"""
        return bond.GetIsAromatic()
    
    def is_conjugated(self, bond):
        """Check if bond is conjugated"""
        return bond.GetIsConjugated()
    
    def get_bond_order(self, bond):
        """Get numeric bond order from bond type"""
        bond_type = bond.GetBondType()
        if bond_type == Chem.rdchem.BondType.SINGLE:
            return 1
        elif bond_type == Chem.rdchem.BondType.DOUBLE:
            return 2
        elif bond_type == Chem.rdchem.BondType.TRIPLE:
            return 3
        elif bond_type == Chem.rdchem.BondType.AROMATIC:
            return 1.5
        else:
            return 1
    
    def calculate_resonance_factor(self, bond, mol):
        """Calculate resonance factor based on bond and molecular environment"""
        # Base factor
        factor = 1.0
        
        # Adjust for aromatics
        if self.is_aromatic(bond):
            factor *= 0.97  # Aromatic bonds are typically shorter than single bonds
            
        # Adjust for conjugation 
        elif self.is_conjugated(bond):
            factor *= 0.98  # Conjugated bonds are slightly shorter
            
        # Ring strain factor
        if self.is_in_ring(bond):
            ring_size = self._get_smallest_ring_size(bond, mol)
            if ring_size == 3:
                factor *= 0.96  # Three-membered rings have significant strain
            elif ring_size == 4:
                factor *= 0.97  # Four-membered rings have moderate strain
                
        return factor
    
    def _get_smallest_ring_size(self, bond, mol):
        """Get the size of the smallest ring containing this bond"""
        if not bond.IsInRing():
            return 0
            
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        
        # Get all rings containing both atoms
        ri = mol.GetRingInfo()
        ring_sizes = []
        
        for ring in ri.AtomRings():
            if atom1_idx in ring and atom2_idx in ring:
                ring_sizes.append(len(ring))
                
        return min(ring_sizes) if ring_sizes else 0
    
    def calculate_bond_length(self, bond, mol):
        """Calculate accurate bond length based on physical properties"""
        atom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        atom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        
        el1 = self.get_element_symbol(atom1.GetAtomicNum())
        el2 = self.get_element_symbol(atom2.GetAtomicNum())
        bond_order = self.get_bond_order(bond)
        
        # Get base bond length from standard values
        base_length = self._get_base_bond_length(el1, el2, bond_order)
        
        # Calculate adjustment factors
        en_diff = abs(self.get_electronegativity(atom1) - self.get_electronegativity(atom2))
        # Amplify the electronegativity effect for visualization
        en_factor = 1.0 - (0.1 * en_diff)  # Increased from 0.05 to 0.1
        
        hyb_factor = (self.get_hybridization_factor(atom1) + self.get_hybridization_factor(atom2)) / 2
        resonance_factor = self.calculate_resonance_factor(bond, mol)
        
        # Make the hybridization effect stronger
        if atom1.GetHybridization() == Chem.rdchem.HybridizationType.SP2 or \
           atom2.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            hyb_factor *= 0.9  # More dramatic effect
        
        # Calculate final bond length
        adjusted_length = base_length * en_factor * hyb_factor * resonance_factor
        
        # For debugging, print significant adjustments
        current_length = self._get_current_bond_length(bond, mol)
        if abs(adjusted_length - current_length) > 0.1:
            print(f"Bond {el1}-{el2} (order {bond_order}): Current {current_length:.2f}Å, "
                  f"Adjusted to {adjusted_length:.2f}Å (diff: {adjusted_length-current_length:.2f}Å)")
        
        return adjusted_length
    
    def _get_base_bond_length(self, el1, el2, bond_order):
        """Get the base bond length from standard values"""
        # Try to find the specific bond in our standards
        for bond_key in [(el1, el2, bond_order), (el2, el1, bond_order)]:
            if bond_key in STANDARD_BOND_LENGTHS:
                return STANDARD_BOND_LENGTHS[bond_key]
        
        # If not found, estimate based on covalent radii
        cov_rad1 = self.periodic_table.GetRcovalent(self.periodic_table.GetAtomicNumber(el1))
        cov_rad2 = self.periodic_table.GetRcovalent(self.periodic_table.GetAtomicNumber(el2))
        
        # Adjust for bond order
        bo_factor = 1.0
        if bond_order == 2:
            bo_factor = 0.85
        elif bond_order == 3:
            bo_factor = 0.75
        elif bond_order == 1.5:  # Aromatic
            bo_factor = 0.9
            
        return (cov_rad1 + cov_rad2) * bo_factor
    
    def process_molecule(self, mol):
        """Process a molecule and calculate accurate bond lengths for all bonds"""
        # Ensure molecule has 3D coordinates
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol)
            
        results = []
        conf = mol.GetConformer()
        
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            
            # Get current bond length from coordinates
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            current_length = np.sqrt((pos1.x - pos2.x)**2 + 
                                    (pos1.y - pos2.y)**2 + 
                                    (pos1.z - pos2.z)**2)
            
            # Calculate accurate bond length
            accurate_length = self.calculate_bond_length(bond, mol)
            
            atom1 = mol.GetAtomWithIdx(idx1)
            atom2 = mol.GetAtomWithIdx(idx2)
            
            results.append({
                'atom1': {
                    'index': idx1,
                    'symbol': self.get_element_symbol(atom1.GetAtomicNum()),
                    'position': [pos1.x, pos1.y, pos1.z]
                },
                'atom2': {
                    'index': idx2,
                    'symbol': self.get_element_symbol(atom2.GetAtomicNum()),
                    'position': [pos2.x, pos2.y, pos2.z]
                },
                'bond_order': self.get_bond_order(bond),
                'is_aromatic': self.is_aromatic(bond),
                'is_conjugated': self.is_conjugated(bond),
                'is_in_ring': self.is_in_ring(bond),
                'current_length': float(current_length),
                'accurate_length': float(accurate_length),
                'adjustment_needed': float(accurate_length - current_length)
            })
            
        return results
        
    def adjust_bond_lengths(self, mol):
        """Generate a new molecule with adjusted bond lengths"""
        # This is a placeholder for future implementation
        # Actually adjusting 3D coordinates while maintaining molecular geometry
        # is complex and would require more advanced algorithms
        
        # For now, we'll just return the bond length data
        return self.process_molecule(mol)

    def _get_current_bond_length(self, bond, mol):
        """Get the current bond length from 3D coordinates"""
        if not mol.GetNumConformers():
            return 0.0
        
        conf = mol.GetConformer()
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        
        pos1 = conf.GetAtomPosition(idx1)
        pos2 = conf.GetAtomPosition(idx2)
        
        return np.sqrt((pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2)

def calculate_molecule_bond_physics(mol_data, data_format='smiles'):
    """Process molecule data and return bond physics information"""
    calculator = BondPhysicsCalculator()
    
    try:
        if data_format == 'smiles':
            # Parse SMILES string
            mol = Chem.MolFromSmiles(mol_data)
            if not mol:
                return {"error": "Invalid SMILES string"}
                
            # Generate 3D coordinates
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            
        elif data_format == 'pdb':
            # Parse PDB data
            print(f"Processing PDB data of length: {len(mol_data)}")
            mol = Chem.MolFromPDBBlock(mol_data)
            if not mol:
                print("Failed to parse PDB data")
                # Try a more lenient PDB parsing approach
                try:
                    # Create a temporary file
                    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
                        f.write(mol_data)
                        temp_file = f.name
                    
                    # Read with more lenient options
                    mol = Chem.MolFromPDBFile(temp_file, sanitize=False, removeHs=False)
                    if mol:
                        try:
                            Chem.SanitizeMol(mol)
                        except:
                            print("Sanitization failed, but will continue with basic structure")
                    
                    # Clean up
                    import os
                    os.unlink(temp_file)
                    
                    if not mol:
                        return {"error": "Failed to parse PDB data even with lenient options"}
                except Exception as e:
                    return {"error": f"Failed alternative PDB parsing: {str(e)}"}
            
            # Add hydrogens if they're missing and we have a workable molecule
            if mol:
                try:
                    mol = Chem.AddHs(mol, addCoords=True)
                except:
                    print("Could not add hydrogens, using molecule as is")
                
        else:
            return {"error": "Unsupported data format"}
        
        if not mol:
            return {"error": "Failed to create molecule object"}
            
        # Print molecule info for debugging
        print(f"Molecule parsed successfully: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")
        
        # Check if we have 3D coordinates
        if mol.GetNumConformers() == 0:
            print("No conformers found, embedding molecule")
            AllChem.EmbedMolecule(mol)
            if mol.GetNumConformers() == 0:
                return {"error": "Failed to generate 3D coordinates"}
        
        # Calculate bond physics
        bond_data = calculator.process_molecule(mol)
        
        # Enhanced debugging
        if bond_data:
            print(f"Processed {len(bond_data)} bonds")
            # Print a few examples
            if len(bond_data) > 0:
                print("Example bond data:")
                for i in range(min(3, len(bond_data))):
                    b = bond_data[i]
                    print(f"Bond {i}: {b['atom1']['symbol']}-{b['atom2']['symbol']}, " 
                          f"Current: {b['current_length']:.2f}Å, "
                          f"Accurate: {b['accurate_length']:.2f}Å, "
                          f"Diff: {b['adjustment_needed']:.2f}Å")
        else:
            print("No bond data generated")
        
        # Return results
        return {
            "bond_data": bond_data,
            "atom_count": mol.GetNumAtoms(),
            "bond_count": mol.GetNumBonds(),
            "debug_info": {
                "has_conformers": mol.GetNumConformers() > 0,
                "format": data_format
            }
        }
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"error": f"Error in bond physics calculation: {str(e)}"}

def create_adjusted_pdb(mol_data, data_format='pdb'):
    """
    Creates a new PDB file with atom coordinates adjusted based on calculated bond physics.
    
    This modifies the actual 3D structure to reflect accurate bond lengths based on
    electronegativity and resonance effects.
    """
    try:
        print(f"Generating physics-adjusted PDB from {data_format} data")
        
        if data_format == 'smiles':
            # Parse SMILES and generate 3D coords
            mol = Chem.MolFromSmiles(mol_data)
            if not mol:
                return {"error": "Invalid SMILES string"}
                
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            
        elif data_format == 'pdb':
            # Parse PDB data
            mol = Chem.MolFromPDBBlock(mol_data)
            if not mol:
                print("Failed to parse PDB data, trying alternative approach")
                # Try a more lenient PDB parsing approach
                try:
                    # Create a temporary file
                    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
                        f.write(mol_data)
                        temp_file = f.name
                    
                    # Read with more lenient options
                    mol = Chem.MolFromPDBFile(temp_file, sanitize=False, removeHs=False)
                    if mol:
                        try:
                            Chem.SanitizeMol(mol)
                        except:
                            print("Sanitization failed, but will continue with basic structure")
                    
                    # Clean up
                    import os
                    os.unlink(temp_file)
                    
                    if not mol:
                        return {"error": "Failed to parse PDB data even with lenient options"}
                except Exception as e:
                    return {"error": f"Failed alternative PDB parsing: {str(e)}"}
            
            # Add hydrogens if they're missing
            try:
                mol = Chem.AddHs(mol, addCoords=True)
            except:
                print("Could not add hydrogens, using molecule as is")
                
        else:
            return {"error": "Unsupported data format"}
            
        if not mol:
            return {"error": "Failed to create molecule object"}
            
        if mol.GetNumConformers() == 0:
            print("No conformers found, embedding molecule")
            AllChem.EmbedMolecule(mol)
            if mol.GetNumConformers() == 0:
                return {"error": "Failed to generate 3D coordinates"}
                
        # Create a copy of the molecule to adjust
        adjusted_mol = Chem.Mol(mol)
        conf = adjusted_mol.GetConformer()
        
        # Set up bond physics calculator
        calculator = BondPhysicsCalculator()
        
        # Calculate bond adjustments
        bond_data = []
        atoms_to_adjust = {}  # Store atom indices and their displacements
        
        # First pass: calculate all the bond physics data
        for bond in adjusted_mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            
            # Get current positions
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            
            # Calculate distance
            current_length = np.sqrt((pos1.x - pos2.x)**2 + 
                                    (pos1.y - pos2.y)**2 + 
                                    (pos1.z - pos2.z)**2)
            
            # Calculate accurate bond length
            accurate_length = calculator.calculate_bond_length(bond, adjusted_mol)
            
            # Calculate adjustment factor
            adjustment = accurate_length / current_length if current_length > 0 else 1.0
            
            # Store bond data
            atom1 = adjusted_mol.GetAtomWithIdx(idx1)
            atom2 = adjusted_mol.GetAtomWithIdx(idx2)
            
            bond_data.append({
                'atom1': {
                    'index': idx1,
                    'symbol': calculator.get_element_symbol(atom1.GetAtomicNum()),
                    'position': [pos1.x, pos1.y, pos1.z]
                },
                'atom2': {
                    'index': idx2,
                    'symbol': calculator.get_element_symbol(atom2.GetAtomicNum()),
                    'position': [pos2.x, pos2.y, pos2.z]
                },
                'bond_order': calculator.get_bond_order(bond),
                'is_aromatic': calculator.is_aromatic(bond),
                'is_conjugated': calculator.is_conjugated(bond),
                'is_in_ring': calculator.is_in_ring(bond),
                'current_length': float(current_length),
                'accurate_length': float(accurate_length),
                'adjustment_needed': float(accurate_length - current_length),
                'adjustment_factor': float(adjustment)
            })
            
            # Track adjustments for each atom
            if idx1 not in atoms_to_adjust:
                atoms_to_adjust[idx1] = []
            if idx2 not in atoms_to_adjust:
                atoms_to_adjust[idx2] = []
                
            # Direction vectors for movement
            dx = pos2.x - pos1.x
            dy = pos2.y - pos1.y
            dz = pos2.z - pos1.z
            
            # Normalize
            norm = np.sqrt(dx*dx + dy*dy + dz*dz)
            if norm > 0:
                dx /= norm
                dy /= norm
                dz /= norm
                
            # Calculate displacement for each atom
            # Move each atom half of the required adjustment
            displacement = (accurate_length - current_length) / 2
            
            # For atom1, move in negative direction of bond vector
            atoms_to_adjust[idx1].append((-dx * displacement, -dy * displacement, -dz * displacement))
            
            # For atom2, move in positive direction of bond vector
            atoms_to_adjust[idx2].append((dx * displacement, dy * displacement, dz * displacement))
            
        # Print adjustment summary
        significant_adjustments = [b for b in bond_data if abs(b['adjustment_needed']) > 0.1]
        print(f"Found {len(significant_adjustments)} bonds with significant length adjustments")
        
        # Second pass: apply the average adjustment to each atom's position
        for atom_idx, displacements in atoms_to_adjust.items():
            if not displacements:
                continue
                
            # Average all displacement vectors for this atom
            avg_dx = sum(d[0] for d in displacements) / len(displacements)
            avg_dy = sum(d[1] for d in displacements) / len(displacements)
            avg_dz = sum(d[2] for d in displacements) / len(displacements)
            
            # Get current position
            pos = conf.GetAtomPosition(atom_idx)
            
            # Apply adjustment - with a scale factor to make changes more visible
            scale_factor = 1.0  # Set to >1.0 to exaggerate changes for visibility
            new_x = pos.x + (avg_dx * scale_factor)
            new_y = pos.y + (avg_dy * scale_factor)
            new_z = pos.z + (avg_dz * scale_factor)
            
            # Update position
            conf.SetAtomPosition(atom_idx, (new_x, new_y, new_z))
            
        # Generate PDB block from adjusted structure
        adjusted_pdb = Chem.MolToPDBBlock(adjusted_mol)
        
        # Calculate bond data for the adjusted molecule for verification
        adjusted_bond_data = []
        for bond in adjusted_mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            
            # Get adjusted positions
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            
            # Calculate new distance
            new_length = np.sqrt((pos1.x - pos2.x)**2 + 
                               (pos1.y - pos2.y)**2 + 
                               (pos1.z - pos2.z)**2)
            
            atom1 = adjusted_mol.GetAtomWithIdx(idx1)
            atom2 = adjusted_mol.GetAtomWithIdx(idx2)
            
            adjusted_bond_data.append({
                'atom1_symbol': calculator.get_element_symbol(atom1.GetAtomicNum()),
                'atom2_symbol': calculator.get_element_symbol(atom2.GetAtomicNum()),
                'new_length': float(new_length)
            })
            
        return {
            "original_pdb": mol_data,
            "adjusted_pdb": adjusted_pdb,
            "bond_data": bond_data,
            "adjusted_bond_data": adjusted_bond_data,
            "atom_count": adjusted_mol.GetNumAtoms(),
            "bond_count": adjusted_mol.GetNumBonds()
        }
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"error": f"Error adjusting molecule structure: {str(e)}"} 