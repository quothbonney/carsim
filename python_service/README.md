# CarsimMD Python Service

This is the Python backend service for the CarsimMD molecular visualization and simulation application.

## Features

- Molecular structure processing and 3D coordinate generation
- Protein structure analysis
- PDB file parsing and visualization support
- Molecular dynamics simulation (NEW)
- Accurate bond length calculation (NEW)
- Physics-based force fields (NEW)

## Setup

1. Make sure you have Python 3.8+ installed
2. Create a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```
4. Run the service:
   ```
   python app.py
   ```

The service will start on port 5000 by default.

## Dependencies

- Flask - Web framework
- RDKit - Cheminformatics toolkit
- NumPy - Numerical computing
- Biopython - Biological computation (optional)

## API Endpoints

### Molecular Structure

- `POST /process_smiles` - Process a SMILES string and return 3D coordinates
- `POST /convert_to_3d` - Convert a SMILES string to 3D PDB format
- `POST /process` - Process a molecule (simplified endpoint)

### Protein Analysis

- `POST /analyze-pdb` - Analyze a PDB structure
- `GET /fetch-pdb/<pdb_id>` - Fetch a protein structure from PDB database
- `POST /binding-sites` - Find potential binding sites in a protein

### Molecular Dynamics (NEW)

- `POST /analyze-bond-lengths` - Analyze bond lengths without changing visualization
- `POST /simulate-molecule` - Create a molecular dynamics simulation
- `POST /simulation/<id>/run` - Run a simulation for a number of steps
- `GET /simulation/<id>/status` - Get simulation status
- `POST /simulation/<id>/stop` - Stop a running simulation
- `GET /simulation/<id>/bond-data` - Get detailed bond data
- `GET /simulation/<id>/trajectory` - Get trajectory data
- `POST /simulate-and-analyze` - Quick simulation and analysis
- `POST /force-field-components` - Get energy components

### Utility

- `GET /health` - Health check endpoint
- `GET /status` - Service status and version info

## Molecular Dynamics Infrastructure

The service now includes a robust molecular dynamics simulation infrastructure. For details, see [README_MOLECULAR_DYNAMICS.md](./README_MOLECULAR_DYNAMICS.md).

## Development

To add new features or fix bugs:

1. Create a new branch
2. Make your changes
3. Update tests if necessary
4. Submit a pull request 