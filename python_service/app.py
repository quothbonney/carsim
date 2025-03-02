from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
import io
import base64
import tempfile
import os
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint to verify the service is running."""
    logger.info("Health check requested")
    return jsonify({"status": "healthy", "service": "molecular-processor"})

@app.route('/process_smiles', methods=['POST'])
def process_smiles():
    """
    Process a SMILES string to generate 3D coordinates and return molecular data.
    
    Expected JSON input: {"smiles": "SMILES_STRING"}
    """
    try:
        data = request.json
        if not data or 'smiles' not in data:
            logger.error("Invalid request: missing SMILES string")
            return jsonify({"error": "Missing SMILES string"}), 400
        
        smiles = data['smiles']
        logger.info(f"Processing SMILES: {smiles}")
        
        # Create RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES string: {smiles}")
            return jsonify({"error": "Invalid SMILES string"}), 400
        
        try:
            # Add hydrogens and generate 3D coordinates
            logger.info("Adding hydrogens to molecule")
            mol = Chem.AddHs(mol)
            
            logger.info("Generating 3D coordinates")
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                logger.error("Failed to generate 3D coordinates")
                return jsonify({"error": "Failed to generate 3D coordinates"}), 500
            
            logger.info("Optimizing molecule geometry")
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Generate PDB string
            logger.info("Generating PDB string")
            pdb_string = Chem.MolToPDBBlock(mol)
            
            # Generate 2D image
            logger.info("Generating 2D image")
            img = Draw.MolToImage(mol, size=(300, 300))
            img_io = io.BytesIO()
            img.save(img_io, format='PNG')
            img_io.seek(0)
            img_data = base64.b64encode(img_io.read()).decode('utf-8')
            
            # Calculate molecular properties
            logger.info("Calculating molecular properties")
            mol_weight = Descriptors.MolWt(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            
            logger.info(f"Successfully processed SMILES: {smiles}")
            
            return jsonify({
                "pdb_string": pdb_string,
                "image_data": img_data,
                "properties": {
                    "molecular_weight": mol_weight,
                    "formula": formula,
                    "num_atoms": num_atoms,
                    "num_bonds": num_bonds
                }
            })
        except Exception as inner_e:
            logger.error(f"Error during molecule processing: {str(inner_e)}")
            return jsonify({"error": f"Error during molecule processing: {str(inner_e)}"}), 500
    
    except Exception as e:
        logger.error(f"Error processing SMILES: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/convert_to_3d', methods=['POST'])
def convert_to_3d():
    """
    Simplified endpoint that only returns the PDB string for a given SMILES.
    
    Expected JSON input: {"smiles": "SMILES_STRING"}
    """
    try:
        data = request.json
        if not data or 'smiles' not in data:
            logger.error("Invalid request: missing SMILES string")
            return jsonify({"error": "Missing SMILES string"}), 400
        
        smiles = data['smiles']
        logger.info(f"Converting SMILES to 3D: {smiles}")
        
        # Create RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES string: {smiles}")
            return jsonify({"error": "Invalid SMILES string"}), 400
        
        try:
            # Add hydrogens and generate 3D coordinates
            logger.info("Adding hydrogens to molecule")
            mol = Chem.AddHs(mol)
            
            logger.info("Generating 3D coordinates")
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                logger.error("Failed to generate 3D coordinates")
                return jsonify({"error": "Failed to generate 3D coordinates"}), 500
            
            logger.info("Optimizing molecule geometry")
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Generate PDB string
            logger.info("Generating PDB string")
            pdb_string = Chem.MolToPDBBlock(mol)
            
            return jsonify({"pdb_string": pdb_string})
        except Exception as inner_e:
            logger.error(f"Error during molecule processing: {str(inner_e)}")
            return jsonify({"error": f"Error during molecule processing: {str(inner_e)}"}), 500
    
    except Exception as e:
        logger.error(f"Error converting SMILES to 3D: {str(e)}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    logger.info(f"Starting molecular processing service on port {port}")
    app.run(host='0.0.0.0', port=port, debug=True) 