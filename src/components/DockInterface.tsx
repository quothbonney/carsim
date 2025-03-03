import React, { useState } from 'react';
import '../styles/DockInterface.css';

// Subwindow types
export enum SubwindowType {
  VISUALIZATION = 'visualization',
  PHYSICS = 'physics',
  NONE = 'none'
}

interface DockInterfaceProps {
  onOpenSubwindow: (type: SubwindowType) => void;
  activeSubwindow: SubwindowType;
  onReloadViewer: () => void;
}

const DockInterface: React.FC<DockInterfaceProps> = ({ 
  onOpenSubwindow, 
  activeSubwindow,
  onReloadViewer
}) => {
  return (
    <div className="dock-interface">
      <button 
        className={`dock-button ${activeSubwindow === SubwindowType.VISUALIZATION ? 'active' : ''}`}
        onClick={() => onOpenSubwindow(SubwindowType.VISUALIZATION)}
      >
        <span className="dock-icon">👁️</span>
        <span className="dock-label">Visualization</span>
      </button>
      
      <button 
        className={`dock-button ${activeSubwindow === SubwindowType.PHYSICS ? 'active' : ''}`}
        onClick={() => onOpenSubwindow(SubwindowType.PHYSICS)}
      >
        <span className="dock-icon">⚛️</span>
        <span className="dock-label">Physics</span>
      </button>
      
      <button 
        className="dock-button dock-reload-button"
        onClick={onReloadViewer}
      >
        <span className="dock-icon">⟳</span>
        <span className="dock-label">Reload</span>
      </button>
    </div>
  );
};

export default DockInterface; 