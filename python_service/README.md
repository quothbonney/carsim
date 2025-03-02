# Molecular Processing Service

This service provides molecular processing capabilities using RDKit. It exposes a REST API for converting SMILES notation to 3D molecular structures and calculating molecular properties.

## Setup

1. Create a virtual environment:
   ```
   python -m venv venv
   ```

2. Activate the virtual environment:
   - Windows: `venv\Scripts\activate`
   - macOS/Linux: `source venv/bin/activate`

3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

## Running the Service

Start the service with:
```
python app.py
```

The service will run on port 5000 by default. You can change this by setting the `PORT` environment variable.

## API Endpoints

### Health Check
```
GET /health
```
Returns the status of the service.

### Process SMILES
```
POST /process_smiles
Content-Type: application/json

{
    "smiles": "SMILES_STRING"
}
```
Converts a SMILES string to a 3D structure and returns:
- PDB string
- 2D image (base64 encoded)
- Molecular properties

### Convert to 3D
```
POST /convert_to_3d
Content-Type: application/json

{
    "smiles": "SMILES_STRING"
}
```
Simplified endpoint that only returns the PDB string for a given SMILES.

## Error Handling

All endpoints return appropriate HTTP status codes and error messages in case of failure. 