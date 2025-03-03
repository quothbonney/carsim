import { useState, useRef } from 'react';
import '../styles/MoleculeControls.css';

interface MoleculeControlsProps {
  smiles: string;
  setSmiles: (smiles: string) => void;
  onProcess: (smiles: string) => void;
  isLoading: boolean;
  onPdbUpload?: (pdbData: string) => void;
}

const MoleculeControls = ({ 
  smiles, 
  setSmiles, 
  onProcess, 
  isLoading,
  onPdbUpload 
}: MoleculeControlsProps) => {
  const [expanded, setExpanded] = useState(true);
  const [pdbExpanded, setPdbExpanded] = useState(true);
  const [pdbId, setPdbId] = useState<string>('');
  const [isFetchingPdb, setIsFetchingPdb] = useState(false);
  const [uploadError, setUploadError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onProcess(smiles);
  };

  const toggleExpanded = () => {
    setExpanded(!expanded);
  };

  const togglePdbExpanded = () => {
    setPdbExpanded(!pdbExpanded);
  };

  // Handle PDB file upload
  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      if (content && onPdbUpload) {
        console.log('PDB file loaded, length:', content.length);
        setUploadError(null);
        onPdbUpload(content);
      }
    };
    reader.onerror = () => {
      setUploadError("Failed to read the uploaded file");
    };
    reader.readAsText(file);
    
    // Reset the input value so the same file can be selected again
    if (event.target) {
      event.target.value = '';
    }
  };

  // Trigger file input click
  const triggerFileUpload = () => {
    if (fileInputRef.current) {
      fileInputRef.current.click();
    }
  };

  // Fetch PDB by ID
  const fetchPdbById = async () => {
    if (!pdbId || pdbId.length !== 4) {
      setUploadError('PDB ID must be 4 characters long');
      return;
    }
    
    try {
      setIsFetchingPdb(true);
      setUploadError(null);
      
      const response = await fetch(`http://localhost:5000/fetch-pdb/${pdbId}`);
      const data = await response.json();
      
      if (data.error) {
        setUploadError(data.error);
        return;
      }
      
      if (data.pdbContent && onPdbUpload) {
        onPdbUpload(data.pdbContent);
      } else {
        setUploadError('Failed to fetch PDB structure');
      }
    } catch (err) {
      setUploadError(`Failed to fetch PDB: ${err}`);
    } finally {
      setIsFetchingPdb(false);
    }
  };

  // Example molecules
  const exampleMolecules = [
    { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
    { name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O' },
    { name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O' },
    { name: 'Penicillin G', smiles: 'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C' },
  ];

  return (
    <div className="molecule-controls">
      <div className="panel-header" onClick={toggleExpanded}>
        <h3>SMILES Input</h3>
        <span className="panel-toggle">{expanded ? '▼' : '▶'}</span>
      </div>
      
      {expanded && (
        <div className="panel-content">
          <form onSubmit={handleSubmit}>
            <div className="form-group">
              <label htmlFor="smiles-input">SMILES String</label>
              <textarea
                id="smiles-input"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                placeholder="Enter SMILES notation..."
                rows={3}
              />
            </div>
            
            <button 
              type="submit" 
              className="process-button"
              disabled={isLoading || !smiles.trim()}
            >
              {isLoading ? 'Processing...' : 'Process Molecule'}
            </button>
          </form>
          
          <div className="examples-section">
            <h4>Example Molecules</h4>
            <div className="examples-list">
              {exampleMolecules.map((molecule) => (
                <button
                  key={molecule.name}
                  className="example-button"
                  onClick={() => setSmiles(molecule.smiles)}
                  disabled={isLoading}
                >
                  {molecule.name}
                </button>
              ))}
            </div>
          </div>
        </div>
      )}

      {/* PDB Upload Section */}
      <div className="panel-header" onClick={togglePdbExpanded}>
        <h3>PDB Structure</h3>
        <span className="panel-toggle">{pdbExpanded ? '▼' : '▶'}</span>
      </div>
      
      {pdbExpanded && (
        <div className="panel-content">
          <div className="pdb-controls">
            <button 
              onClick={triggerFileUpload} 
              className="upload-button"
              disabled={isLoading || isFetchingPdb}
            >
              <span className="icon">📁</span> Upload PDB
            </button>
            
            <input
              ref={fileInputRef}
              type="file"
              accept=".pdb"
              onChange={handleFileUpload}
              style={{ display: 'none' }}
            />
            
            <div className="pdb-id-input-container">
              <input
                type="text"
                value={pdbId}
                onChange={(e) => setPdbId(e.target.value.toUpperCase())}
                placeholder="PDB ID (e.g. 1CRN)"
                className="pdb-id-input"
                maxLength={4}
                disabled={isLoading || isFetchingPdb}
              />
              <button 
                onClick={fetchPdbById} 
                className="fetch-button"
                disabled={isLoading || isFetchingPdb || !pdbId || pdbId.length !== 4}
              >
                {isFetchingPdb ? 'Fetching...' : 'Fetch'}
              </button>
            </div>
            
            {uploadError && (
              <div className="upload-error">
                {uploadError}
              </div>
            )}
          </div>
          
          <div className="help-section">
            <h4>Help</h4>
            <p>
              Upload a PDB file or enter a 4-character PDB ID to visualize protein structures.
            </p>
            <p>
              PDB files contain atomic coordinate data for molecular structures from the 
              Protein Data Bank (rcsb.org).
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default MoleculeControls; 