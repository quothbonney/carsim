import React, { useState } from 'react';
import '../styles/Subwindow.css';

interface PhysicsSubwindowProps {
  onClose: () => void;
  onApplyCorrectBondLengths: () => void;
  onResetToGenericBondLengths: () => void;
  isInitialized: boolean;
}

const PhysicsSubwindow: React.FC<PhysicsSubwindowProps> = ({
  onClose,
  onApplyCorrectBondLengths,
  onResetToGenericBondLengths,
  isInitialized
}) => {
  const [isCalculating, setIsCalculating] = useState(false);
  const [correctBondsEnabled, setCorrectBondsEnabled] = useState(false);
  
  const handleApplyCorrectBondLengths = () => {
    setIsCalculating(true);
    // Simulate calculation taking time
    setTimeout(() => {
      onApplyCorrectBondLengths();
      setIsCalculating(false);
      setCorrectBondsEnabled(true);
    }, 1000);
  };
  
  const handleResetToGeneric = () => {
    onResetToGenericBondLengths();
    setCorrectBondsEnabled(false);
  };
  
  return (
    <div className="subwindow physics-subwindow">
      <div className="subwindow-header">
        <h3>Physics Options</h3>
        <button className="close-button" onClick={onClose}>×</button>
      </div>
      
      <div className="subwindow-content">
        <div className="control-section">
          <h4>Bond Lengths</h4>
          <div className="option-description">
            <p>Molecular dynamics can calculate correct bond lengths based on atomic properties, electronegativity, and more.</p>
          </div>
          <div className="button-group">
            <button 
              className={`control-button ${correctBondsEnabled ? 'active' : ''}`}
              onClick={handleApplyCorrectBondLengths}
              disabled={!isInitialized || isCalculating || correctBondsEnabled}
            >
              <span className="icon">🔬</span>
              <span className="label">
                {isCalculating ? 'Calculating...' : 'Apply Correct Bond Lengths'}
              </span>
            </button>
            
            <button 
              className="control-button"
              onClick={handleResetToGeneric}
              disabled={!isInitialized || !correctBondsEnabled}
            >
              <span className="icon">🔄</span>
              <span className="label">Reset to Generic Bond Lengths</span>
            </button>
          </div>
        </div>
        
        <div className="info-box">
          <h4>Physics Simulation</h4>
          <p>
            Our molecular dynamics infrastructure uses advanced physics calculations
            to determine realistic bond lengths based on multiple factors:
          </p>
          <ul>
            <li>Electronegativity differences between atoms</li>
            <li>Bond order (single, double, triple)</li>
            <li>Hybridization state</li>
            <li>Resonance effects</li>
          </ul>
        </div>
      </div>
    </div>
  );
};

export default PhysicsSubwindow; 