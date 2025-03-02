import { useState } from 'react';
import '../styles/MoleculeInfo.css';

interface MoleculeInfoProps {
  properties: {
    molecular_weight: number;
    formula: string;
    num_atoms: number;
    num_bonds: number;
  };
  imageData: string;
}

const MoleculeInfo = ({ properties, imageData }: MoleculeInfoProps) => {
  const [propertiesExpanded, setPropertiesExpanded] = useState(true);
  const [structureExpanded, setStructureExpanded] = useState(true);

  const toggleProperties = () => {
    setPropertiesExpanded(!propertiesExpanded);
  };

  const toggleStructure = () => {
    setStructureExpanded(!structureExpanded);
  };

  return (
    <div className="molecule-info">
      <div className="panel-section">
        <div className="panel-header" onClick={toggleProperties}>
          <h4>Properties</h4>
          <span className="panel-toggle">{propertiesExpanded ? '▼' : '▶'}</span>
        </div>
        
        {propertiesExpanded && (
          <div className="panel-content">
            <table className="properties-table">
              <tbody>
                <tr>
                  <td>Formula</td>
                  <td>{properties.formula}</td>
                </tr>
                <tr>
                  <td>Molecular Weight</td>
                  <td>{properties.molecular_weight.toFixed(2)} g/mol</td>
                </tr>
                <tr>
                  <td>Atoms</td>
                  <td>{properties.num_atoms}</td>
                </tr>
                <tr>
                  <td>Bonds</td>
                  <td>{properties.num_bonds}</td>
                </tr>
              </tbody>
            </table>
          </div>
        )}
      </div>
      
      <div className="panel-section">
        <div className="panel-header" onClick={toggleStructure}>
          <h4>2D Structure</h4>
          <span className="panel-toggle">{structureExpanded ? '▼' : '▶'}</span>
        </div>
        
        {structureExpanded && (
          <div className="panel-content">
            <div className="structure-image">
              <img src={`data:image/png;base64,${imageData}`} alt="Molecular structure" />
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default MoleculeInfo; 