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
from protein_analysis import ProteinAnalyzer

# Import our new modules
from bond_physics import BondPhysicsCalculator, calculate_molecule_bond_physics
from bond_length_analyzer import calculate_ideal_bond_lengths, analyze_molecule_using_simulator, get_complete_molecule_analysis
from molecular_dynamics import create_simulation_from_mol_data, MolecularDynamicsSimulation
import numpy as np

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

# Initialize the protein analyzer
protein_analyzer = ProteinAnalyzer()

# Store ongoing simulations by ID
active_simulations = {}

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

@app.route('/process', methods=['POST'])
def process_molecule():
    try:
        data = request.json
        smiles = data.get('smiles', '')
        
        if not smiles:
            return jsonify({'error': 'No SMILES string provided'}), 400
        
        # Create a molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string'}), 400
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Convert to PDB format
        pdb_string = Chem.MolToPDBBlock(mol)
        
        # Return PDB string and molecule info
        return jsonify({
            'pdb': pdb_string,
            'num_atoms': mol.GetNumAtoms(),
            'num_bonds': mol.GetNumBonds(),
            'formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/analyze-pdb', methods=['POST'])
def analyze_pdb():
    """Analyze a PDB file and return structural information"""
    try:
        data = request.json
        pdb_content = data.get('pdbContent', '')
        
        if not pdb_content:
            return jsonify({'error': 'No PDB content provided'}), 400
            
        logger.info("Analyzing PDB content of length %d", len(pdb_content))
        
        # Analyze the PDB content
        result = protein_analyzer.analyze_pdb_content(pdb_content)
        
        # Check if there was an error
        if 'error' in result:
            logger.error("Error analyzing PDB: %s", result['error'])
            return jsonify(result), 400
        
        logger.info("PDB analysis complete")
        # Return the analysis results
        return jsonify(result)
        
    except Exception as e:
        logger.exception("Exception in analyze_pdb: %s", str(e))
        return jsonify({'error': str(e)}), 500

@app.route('/fetch-pdb/<pdb_id>', methods=['GET'])
def fetch_pdb(pdb_id):
    """Fetch a protein structure from the PDB database by ID"""
    try:
        if not pdb_id or len(pdb_id) != 4:
            return jsonify({'error': 'Invalid PDB ID format. Should be 4 characters.'}), 400
            
        logger.info("Fetching PDB with ID: %s", pdb_id)
        
        # Get the protein structure
        result = protein_analyzer.get_protein_from_pdb_id(pdb_id)
        
        # Check if there was an error
        if 'error' in result:
            logger.error("Error fetching PDB %s: %s", pdb_id, result['error'])
            return jsonify(result), 400
            
        logger.info("Successfully fetched PDB %s", pdb_id)
        
        # Return the structure and analysis
        return jsonify(result)
        
    except Exception as e:
        logger.exception("Exception in fetch_pdb: %s", str(e))
        return jsonify({'error': str(e)}), 500

@app.route('/binding-sites', methods=['POST'])
def find_binding_sites():
    """Find potential binding sites in a protein structure"""
    try:
        data = request.json
        pdb_content = data.get('pdbContent', '')
        distance_cutoff = data.get('distanceCutoff', 5.0)
        
        if not pdb_content:
            return jsonify({'error': 'No PDB content provided'}), 400
            
        # Find binding sites
        result = protein_analyzer.find_binding_sites(pdb_content, distance_cutoff)
        
        # Return the binding sites
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/status', methods=['GET'])
def get_status():
    """Get the status of the protein analysis service"""
    biopython_available = protein_analyzer.is_available()
    
    return jsonify({
        'status': 'ok',
        'version': '1.0.0',
        'biopython_available': biopython_available,
        'rdkit_version': Chem.__version__
    })

@app.route('/analyze-bond-lengths', methods=['POST'])
def analyze_bond_lengths():
    """
    Analyze bond lengths of a molecule without modifying it.
    
    Expected JSON input: {"pdbContent": "PDB_STRING"} or {"smiles": "SMILES_STRING"}
    """
    try:
        data = request.json
        
        if not data:
            return jsonify({"error": "No data provided"}), 400
            
        if "pdbContent" in data:
            mol_data = data["pdbContent"]
            data_format = "pdb"
            logger.info("Analyzing bond lengths from PDB data")
        elif "smiles" in data:
            mol_data = data["smiles"]
            data_format = "smiles"
            logger.info(f"Analyzing bond lengths from SMILES: {mol_data}")
        else:
            return jsonify({"error": "No molecular data found. Provide either 'pdbContent' or 'smiles'."}), 400
        
        # Use function from bond_length_analyzer module
        results = calculate_ideal_bond_lengths(mol_data, data_format)
        
        if "error" in results:
            logger.error(f"Error analyzing bond lengths: {results['error']}")
            return jsonify(results), 500
            
        logger.info(f"Bond length analysis completed: {len(results.get('bond_data', []))} bonds analyzed")
        return jsonify(results)
        
    except Exception as e:
        logger.exception(f"Error in bond length analysis: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulate-molecule', methods=['POST'])
def create_molecule_simulation():
    """
    Create a molecular dynamics simulation for a molecule.
    
    Expected JSON input: 
    {
        "pdbContent": "PDB_STRING" or "smiles": "SMILES_STRING",
        "numSteps": 1000,                 # Optional: number of simulation steps
        "temperature": 300,               # Optional: temperature in K
        "storeTrajectory": true/false     # Optional: whether to store trajectory
    }
    
    Returns: Simulation ID and initial state
    """
    try:
        data = request.json
        
        if not data:
            return jsonify({"error": "No data provided"}), 400
            
        if "pdbContent" in data:
            mol_data = data["pdbContent"]
            data_format = "pdb"
        elif "smiles" in data:
            mol_data = data["smiles"]
            data_format = "smiles"
        else:
            return jsonify({"error": "No molecular data found. Provide either 'pdbContent' or 'smiles'."}), 400
        
        # Create simulation
        try:
            logger.info(f"Creating simulation for {data_format} data")
            simulation = create_simulation_from_mol_data(mol_data, data_format)
        except ValueError as ve:
            logger.error(f"Error creating simulation: {str(ve)}")
            return jsonify({"error": str(ve)}), 400
        
        # Generate a unique ID for this simulation
        import uuid
        simulation_id = str(uuid.uuid4())
        
        # Store in active simulations
        active_simulations[simulation_id] = simulation
        
        # Set optional parameters
        if "temperature" in data:
            simulation.temperature = float(data["temperature"])
            simulation.initialize_velocities(simulation.temperature)
        
        # Return simulation ID and initial state
        initial_bond_data = simulation.get_bond_data()
        initial_state = simulation.get_current_state()
        
        return jsonify({
            "simulation_id": simulation_id,
            "initial_state": initial_state,
            "bond_data": initial_bond_data,
            "molecule_info": {
                "num_atoms": simulation.num_atoms,
                "atom_masses": simulation.atom_masses.tolist()
            }
        })
        
    except Exception as e:
        logger.exception(f"Error creating simulation: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulation/<simulation_id>/run', methods=['POST'])
def run_simulation(simulation_id):
    """
    Run a created simulation for a specified number of steps.
    
    Expected JSON input: 
    {
        "numSteps": 1000,                 # Number of simulation steps
        "storeTrajectory": true/false     # Optional: whether to store trajectory
    }
    """
    try:
        if simulation_id not in active_simulations:
            return jsonify({"error": "Simulation not found"}), 404
            
        data = request.json
        if not data or "numSteps" not in data:
            return jsonify({"error": "numSteps parameter is required"}), 400
            
        simulation = active_simulations[simulation_id]
        num_steps = int(data["numSteps"])
        store_trajectory = data.get("storeTrajectory", False)
        
        # Run simulation
        logger.info(f"Running simulation {simulation_id} for {num_steps} steps")
        simulation.run_simulation(num_steps, store_trajectory)
        
        return jsonify({
            "status": "simulation started",
            "simulation_id": simulation_id,
            "steps_requested": num_steps
        })
        
    except Exception as e:
        logger.exception(f"Error running simulation: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulation/<simulation_id>/status', methods=['GET'])
def get_simulation_status(simulation_id):
    """Get the current status of a simulation."""
    try:
        if simulation_id not in active_simulations:
            return jsonify({"error": "Simulation not found"}), 404
            
        simulation = active_simulations[simulation_id]
        state = simulation.get_current_state()
        
        # Include bond data if simulation is not running (to avoid slowing it down)
        if simulation.state.value != "running":
            bond_data = simulation.get_bond_data()
            state["bond_data"] = bond_data
        
        return jsonify(state)
        
    except Exception as e:
        logger.exception(f"Error getting simulation status: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulation/<simulation_id>/stop', methods=['POST'])
def stop_simulation(simulation_id):
    """Stop a running simulation."""
    try:
        if simulation_id not in active_simulations:
            return jsonify({"error": "Simulation not found"}), 404
            
        simulation = active_simulations[simulation_id]
        simulation.stop_simulation()
        
        return jsonify({
            "status": "simulation stopped",
            "simulation_id": simulation_id,
            "steps_completed": simulation.current_step
        })
        
    except Exception as e:
        logger.exception(f"Error stopping simulation: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulation/<simulation_id>/bond-data', methods=['GET'])
def get_simulation_bond_data(simulation_id):
    """Get detailed bond data from a simulation."""
    try:
        if simulation_id not in active_simulations:
            return jsonify({"error": "Simulation not found"}), 404
            
        simulation = active_simulations[simulation_id]
        bond_data = simulation.get_bond_data()
        
        # Calculate statistics
        if bond_data:
            differences = [bond["difference"] for bond in bond_data]
            avg_difference = sum(differences) / len(differences)
            max_difference = max(differences, key=abs)
            min_difference = min(differences, key=abs)
            rms_deviation = np.sqrt(sum(diff**2 for diff in differences) / len(differences))
            
            return jsonify({
                "bond_data": bond_data,
                "statistics": {
                    "average_difference": float(avg_difference),
                    "max_difference": float(max_difference),
                    "min_difference": float(min_difference),
                    "rms_deviation": float(rms_deviation),
                    "total_bond_count": len(bond_data)
                }
            })
        else:
            return jsonify({"error": "No bond data available", "bond_data": []})
        
    except Exception as e:
        logger.exception(f"Error getting bond data: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulation/<simulation_id>/trajectory', methods=['GET'])
def get_simulation_trajectory(simulation_id):
    """Get the trajectory data from a completed simulation."""
    try:
        if simulation_id not in active_simulations:
            return jsonify({"error": "Simulation not found"}), 404
            
        simulation = active_simulations[simulation_id]
        
        if not simulation.store_trajectory:
            return jsonify({"error": "Trajectory was not stored for this simulation"}), 400
            
        return jsonify({
            "trajectory": [positions.tolist() for positions in simulation.trajectory],
            "energies": simulation.energies,
            "steps": len(simulation.trajectory),
            "atom_count": simulation.num_atoms
        })
        
    except Exception as e:
        logger.exception(f"Error getting trajectory: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/simulate-and-analyze', methods=['POST'])
def simulate_and_analyze():
    """
    Create a simulation, run it briefly, and return analysis results.
    
    This is a convenience endpoint that combines simulation creation, 
    running, and analysis in one call.
    
    Expected JSON input: 
    {
        "pdbContent": "PDB_STRING" or "smiles": "SMILES_STRING",
        "numSteps": 100                   # Optional: number of simulation steps
    }
    """
    try:
        data = request.json
        
        if not data:
            return jsonify({"error": "No data provided"}), 400
            
        if "pdbContent" in data:
            mol_data = data["pdbContent"]
            data_format = "pdb"
        elif "smiles" in data:
            mol_data = data["smiles"]
            data_format = "smiles"
        else:
            return jsonify({"error": "No molecular data found. Provide either 'pdbContent' or 'smiles'."}), 400
        
        # Get the number of steps
        num_steps = int(data.get("numSteps", 100))
        
        # Create and run simulation
        try:
            logger.info(f"Creating simulation for {data_format} data")
            simulation = create_simulation_from_mol_data(mol_data, data_format)
            
            # Initialize velocities
            simulation.initialize_velocities()
            
            # Run simulation synchronously (not in a thread)
            logger.info(f"Running quick simulation for {num_steps} steps")
            
            # Store initial state
            initial_bond_data = simulation.get_bond_data()
            
            # Run manually instead of using run_simulation (which is async)
            for step in range(num_steps):
                # Calculate forces at current positions
                forces = simulation.force_field.calculate_forces(simulation.mol, simulation.positions)
                
                # Calculate accelerations (F = ma)
                accelerations = forces / simulation.atom_masses_kg.reshape(-1, 1)
                
                # Update velocities (half step)
                simulation.velocities += 0.5 * accelerations * simulation.time_step
                
                # Update positions
                simulation.positions += simulation.velocities * simulation.time_step
                
                # Recalculate forces with new positions
                forces = simulation.force_field.calculate_forces(simulation.mol, simulation.positions)
                
                # Recalculate accelerations
                accelerations = forces / simulation.atom_masses_kg.reshape(-1, 1)
                
                # Update velocities (second half step)
                simulation.velocities += 0.5 * accelerations * simulation.time_step
            
            # Get final bond data
            final_bond_data = simulation.get_bond_data()
            
            # Calculate statistics
            diffs_before = [bond["difference"] for bond in initial_bond_data]
            diffs_after = [bond["difference"] for bond in final_bond_data]
            
            rms_before = np.sqrt(sum(d**2 for d in diffs_before) / len(diffs_before))
            rms_after = np.sqrt(sum(d**2 for d in diffs_after) / len(diffs_after))
            
            return jsonify({
                "initial_bond_data": initial_bond_data,
                "final_bond_data": final_bond_data,
                "statistics": {
                    "initial_rms_deviation": float(rms_before),
                    "final_rms_deviation": float(rms_after),
                    "improvement": float(rms_before - rms_after),
                    "steps_run": num_steps
                },
                "energy": float(simulation.force_field.calculate_energy(
                    simulation.mol, simulation.positions))
            })
            
        except ValueError as ve:
            logger.error(f"Error in simulation: {str(ve)}")
            return jsonify({"error": str(ve)}), 400
        
    except Exception as e:
        logger.exception(f"Error in simulate-and-analyze: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/force-field-components', methods=['POST'])
def get_force_field_components():
    """
    Get the energy components from a force field for a molecule.
    
    Expected JSON input: 
    {
        "pdbContent": "PDB_STRING" or "smiles": "SMILES_STRING"
    }
    """
    try:
        data = request.json
        
        if not data:
            return jsonify({"error": "No data provided"}), 400
            
        if "pdbContent" in data:
            mol_data = data["pdbContent"]
            data_format = "pdb"
        elif "smiles" in data:
            mol_data = data["smiles"]
            data_format = "smiles"
        else:
            return jsonify({"error": "No molecular data found. Provide either 'pdbContent' or 'smiles'."}), 400
        
        # Create simulation with force field
        try:
            simulation = create_simulation_from_mol_data(mol_data, data_format)
            
            # Calculate energy components
            total_energy = simulation.force_field.calculate_energy(
                simulation.mol, simulation.positions)
            
            # Convert energy components to serializable format
            energy_components = {
                force_type.value: float(energy) 
                for force_type, energy in simulation.force_field.energy_components.items()
            }
            
            return jsonify({
                "total_energy": float(total_energy),
                "energy_components": energy_components,
                "molecule_info": {
                    "num_atoms": simulation.num_atoms,
                    "num_bonds": simulation.mol.GetNumBonds()
                }
            })
            
        except ValueError as ve:
            logger.error(f"Error in force field calculation: {str(ve)}")
            return jsonify({"error": str(ve)}), 400
        
    except Exception as e:
        logger.exception(f"Error calculating force field components: {str(e)}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    logger.info(f"Starting molecular processing service on port {port}")
    app.run(host='0.0.0.0', port=port, debug=True) 