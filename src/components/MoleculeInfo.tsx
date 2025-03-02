import React from 'react';
import '../styles/MoleculeInfo.css';

interface MolecularProperties {
  molecular_weight: number;
  formula: string;
  num_atoms: number;
  num_bonds: number;
}

interface MoleculeInfoProps {
  properties: MolecularProperties;
  imageData: string;
}

const MoleculeInfo: React.FC<MoleculeInfoProps> = ({ properties, imageData }) => {
  return (
    <div className="molecule-info">
      <h2>Molecule Information</h2>
      
      <div className="info-container">
        <div className="molecule-image">
          <h3>2D Structure</h3>
          {imageData ? (
            <img 
              src={`data:image/png;base64,${imageData}`} 
              alt="2D Molecular Structure" 
            />
          ) : (
            <div className="no-image">No image available</div>
          )}
        </div>
        
        <div className="molecule-properties">
          <h3>Properties</h3>
          <table>
            <tbody>
              <tr>
                <td>Molecular Formula:</td>
                <td>{properties.formula}</td>
              </tr>
              <tr>
                <td>Molecular Weight:</td>
                <td>{properties.molecular_weight.toFixed(2)} g/mol</td>
              </tr>
              <tr>
                <td>Number of Atoms:</td>
                <td>{properties.num_atoms}</td>
              </tr>
              <tr>
                <td>Number of Bonds:</td>
                <td>{properties.num_bonds}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
};

export default MoleculeInfo; 