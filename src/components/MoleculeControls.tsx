import { useState } from 'react';
import '../styles/MoleculeControls.css';

interface MoleculeControlsProps {
  smiles: string;
  setSmiles: (smiles: string) => void;
  onProcess: (smiles: string) => void;
  isLoading: boolean;
}

const MoleculeControls = ({ smiles, setSmiles, onProcess, isLoading }: MoleculeControlsProps) => {
  const [expanded, setExpanded] = useState(true);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onProcess(smiles);
  };

  const toggleExpanded = () => {
    setExpanded(!expanded);
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
        <h3>Molecule Input</h3>
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
          
          <div className="help-section">
            <h4>Help</h4>
            <p>
              Enter a SMILES string to visualize a molecule in 3D. SMILES (Simplified Molecular Input Line Entry System) 
              is a notation that represents molecular structures as text strings.
            </p>
            <p>
              You can also select one of the example molecules above to get started.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default MoleculeControls; 