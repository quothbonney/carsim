"""
CarsimMD Protein Analysis Module

This module provides advanced protein structure analysis functionality using BioPython.
It includes functions for analyzing protein structures, calculating structural properties,
and preparing data for visualization in the CarsimMD web application.
"""

from flask import jsonify
import os
import tempfile
import numpy as np

try:
    from Bio.PDB import PDBParser, DSSP, NeighborSearch, PDBList
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

class ProteinAnalyzer:
    """Class for analyzing protein structures"""
    
    def __init__(self):
        """Initialize the analyzer with a PDB parser"""
        if BIOPYTHON_AVAILABLE:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None
            
    def is_available(self):
        """Check if BioPython is available"""
        return BIOPYTHON_AVAILABLE
        
    def analyze_pdb_content(self, pdb_content):
        """Analyze PDB content and return structural information"""
        if not BIOPYTHON_AVAILABLE:
            return {
                "error": "BioPython is not installed. Install it with: pip install biopython"
            }
            
        try:
            # Create a temporary file to write the PDB content
            with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', mode='w') as f:
                f.write(pdb_content)
                temp_path = f.name
                
            # Parse the PDB file
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', temp_path)
            
            # Get the first model
            model = structure[0]
            
            # Count residues and atoms
            residue_count = sum(1 for _ in model.get_residues())
            atom_count = sum(1 for _ in model.get_atoms())
            
            # Get chain information
            chains = []
            for chain in model:
                chain_residues = list(chain.get_residues())
                chains.append({
                    "chainId": chain.id,
                    "residueCount": len(chain_residues)
                })
            
            # Try to calculate secondary structure with DSSP if available
            secondary_structure = self._calculate_secondary_structure(model, temp_path)
            
            # Clean up the temporary file
            os.unlink(temp_path)
            
            # Compile the analysis results
            result = {
                "residueCount": residue_count,
                "atomCount": atom_count,
                "chains": chains
            }
            
            if secondary_structure:
                result["secondaryStructure"] = secondary_structure
                
            return result
            
        except Exception as e:
            return {"error": f"Error analyzing PDB: {str(e)}"}
    
    def _get_structure_stats(self, structure):
        """Extract basic statistics from the structure"""
        stats = {
            "chains": [],
            "totalResidues": 0,
            "totalAtoms": 0,
            "secondaryStructure": {
                "helix": 0,
                "sheet": 0,
                "loop": 0
            }
        }
        
        # Count atoms and residues
        for model in structure:
            for chain in model:
                chain_stats = {
                    "chainId": chain.id,
                    "residueCount": 0
                }
                
                for residue in chain:
                    if residue.id[0] == ' ':  # standard amino acid
                        stats["totalResidues"] += 1
                        chain_stats["residueCount"] += 1
                    
                    for atom in residue:
                        stats["totalAtoms"] += 1
                
                stats["chains"].append(chain_stats)
                
        # Try to calculate secondary structure if DSSP is available
        try:
            model = structure[0]
            dssp = DSSP(model, None, dssp='mkdssp')
            
            for key in dssp.keys():
                ss = dssp[key][2]
                if ss in ('H', 'G', 'I'):  # Helix types
                    stats["secondaryStructure"]["helix"] += 1
                elif ss in ('E', 'B'):     # Sheet types
                    stats["secondaryStructure"]["sheet"] += 1
                else:                      # Loops and others
                    stats["secondaryStructure"]["loop"] += 1
        except:
            # DSSP might not be available, or structure might not be suitable
            pass
            
        return stats
    
    def get_protein_from_pdb_id(self, pdb_id):
        """Download a protein structure from the PDB database"""
        if not BIOPYTHON_AVAILABLE:
            return {"error": "BioPython is not installed. Install it with: pip install biopython"}
            
        try:
            pdbl = PDBList()
            filename = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
            
            with open(filename, 'r') as f:
                pdb_content = f.read()
                
            # Get analysis
            stats = self.analyze_pdb_content(pdb_content)
            
            return {
                "pdbContent": pdb_content,
                "analysis": stats
            }
        except Exception as e:
            return {"error": f"Failed to retrieve PDB {pdb_id}: {str(e)}"}
            
    def find_binding_sites(self, pdb_content, distance_cutoff=5.0):
        """Find potential binding sites in the protein"""
        if not BIOPYTHON_AVAILABLE:
            return {"error": "BioPython is not installed. Install it with: pip install biopython"}
            
        try:
            # Write content to a temporary file
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                tmp.write(pdb_content.encode('utf-8'))
                tmp_path = tmp.name
                
            # Parse the PDB file
            structure = self.parser.get_structure('protein', tmp_path)
            model = structure[0]
            
            # Identify hetero atoms (typically ligands)
            hetero_atoms = []
            protein_atoms = []
            
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':  # Not a standard amino acid
                        for atom in residue:
                            hetero_atoms.append(atom)
                    else:
                        for atom in residue:
                            protein_atoms.append(atom)
            
            # If no hetero atoms, return an empty result
            if not hetero_atoms:
                os.unlink(tmp_path)
                return {"bindingSites": []}
            
            # Use neighbor search to find interactions
            ns = NeighborSearch(protein_atoms)
            binding_sites = {}
            
            for atom in hetero_atoms:
                nearby = ns.search(atom.coord, distance_cutoff)
                for nearby_atom in nearby:
                    residue = nearby_atom.get_parent()
                    chain_id = residue.get_parent().id
                    res_id = residue.id[1]
                    res_name = residue.resname
                    
                    key = f"{chain_id}:{res_id} {res_name}"
                    if key not in binding_sites:
                        binding_sites[key] = {
                            "chainId": chain_id,
                            "residueId": res_id,
                            "residueName": res_name,
                            "atoms": []
                        }
                    
                    if nearby_atom.name not in binding_sites[key]["atoms"]:
                        binding_sites[key]["atoms"].append(nearby_atom.name)
            
            # Clean up the temporary file
            os.unlink(tmp_path)
            
            return {"bindingSites": list(binding_sites.values())}
        
        except Exception as e:
            return {"error": f"Failed to analyze binding sites: {str(e)}"} 