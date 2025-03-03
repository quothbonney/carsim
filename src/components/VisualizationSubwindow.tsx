import React from 'react';
import '../styles/Subwindow.css';

interface VisualizationSubwindowProps {
  onClose: () => void;
  viewMode: string;
  onChangeViewMode: (mode: string) => void;
  showSurface: boolean;
  onToggleSurface: () => void;
  onResetView: () => void;
  isInitialized: boolean;
}

const VisualizationSubwindow: React.FC<VisualizationSubwindowProps> = ({
  onClose,
  viewMode,
  onChangeViewMode,
  showSurface,
  onToggleSurface,
  onResetView,
  isInitialized
}) => {
  return (
    <div className="subwindow visualization-subwindow">
      <div className="subwindow-header">
        <h3>Visualization Options</h3>
        <button className="close-button" onClick={onClose}>×</button>
      </div>
      
      <div className="subwindow-content">
        <div className="control-section">
          <h4>View Mode</h4>
          <div className="button-group">
            <button 
              className={`control-button ${viewMode === 'stick' ? 'active' : ''}`}
              onClick={() => onChangeViewMode('stick')}
              disabled={!isInitialized}
            >
              <span className="icon">🧩</span>
              <span className="label">Stick</span>
            </button>
            
            <button 
              className={`control-button ${viewMode === 'cartoon' ? 'active' : ''}`}
              onClick={() => onChangeViewMode('cartoon')}
              disabled={!isInitialized}
            >
              <span className="icon">🧵</span>
              <span className="label">Ribbon</span>
            </button>
            
            <button 
              className={`control-button ${viewMode === 'sphere' ? 'active' : ''}`}
              onClick={() => onChangeViewMode('sphere')}
              disabled={!isInitialized}
            >
              <span className="icon">⚪</span>
              <span className="label">Sphere</span>
            </button>
            
            <button 
              className={`control-button ${viewMode === 'line' ? 'active' : ''}`}
              onClick={() => onChangeViewMode('line')}
              disabled={!isInitialized}
            >
              <span className="icon">╱</span>
              <span className="label">Line</span>
            </button>
          </div>
        </div>
        
        <div className="control-section">
          <h4>Surface Options</h4>
          <div className="button-group">
            <button 
              className={`control-button ${showSurface ? 'active' : ''}`}
              onClick={onToggleSurface}
              disabled={!isInitialized}
            >
              <span className="icon">◯</span>
              <span className="label">Surface</span>
            </button>
          </div>
        </div>
        
        <div className="control-section">
          <h4>Camera</h4>
          <div className="button-group">
            <button 
              className="control-button"
              onClick={onResetView}
              disabled={!isInitialized}
            >
              <span className="icon">⊕</span>
              <span className="label">Reset View</span>
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};

export default VisualizationSubwindow; 