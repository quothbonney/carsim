import React from 'react';
import '../styles/MoleculeControls.css';

interface MoleculeControlsProps {
  smiles: string;
  setSmiles: (smiles: string) => void;
  onProcess: (smiles: string) => void;
  isLoading: boolean;
}

const MoleculeControls: React.FC<MoleculeControlsProps> = ({
  smiles,
  setSmiles,
  onProcess,
  isLoading
}) => {
  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onProcess(smiles);
  };

  const exampleMolecules = [
    { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
    { name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O' },
    { name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O' },
    { name: 'Penicillin G', smiles: 'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C' }
  ];

  return (
    <div className="molecule-controls">
      <h2>Molecule Input</h2>
      <form onSubmit={handleSubmit}>
        <div className="input-group">
          <label htmlFor="smiles-input">SMILES Notation:</label>
          <input
            id="smiles-input"
            type="text"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="Enter SMILES notation"
            disabled={isLoading}
          />
        </div>
        <button type="submit" disabled={isLoading || !smiles.trim()}>
          {isLoading ? 'Processing...' : 'Process Molecule'}
        </button>
      </form>

      <div className="example-molecules">
        <h3>Example Molecules</h3>
        <div className="examples-list">
          {exampleMolecules.map((molecule) => (
            <button
              key={molecule.name}
              onClick={() => setSmiles(molecule.smiles)}
              disabled={isLoading}
              className="example-button"
            >
              {molecule.name}
            </button>
          ))}
        </div>
      </div>

      <div className="smiles-help">
        <h3>About SMILES Notation</h3>
        <p>
          SMILES (Simplified Molecular Input Line Entry System) is a notation that
          represents molecular structures as text strings. It encodes the structure
          of a molecule in a compact and human-readable format.
        </p>
        <p>
          <a
            href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system"
            target="_blank"
            rel="noopener noreferrer"
          >
            Learn more about SMILES notation
          </a>
        </p>
      </div>
    </div>
  );
};

export default MoleculeControls; 