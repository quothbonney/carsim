import { useEffect, useRef } from 'react';
import '../styles/MoleculeViewer.css';

interface MoleculeViewerProps {
  pdbData: string;
  isLoading: boolean;
}

const MoleculeViewer = ({ pdbData, isLoading }: MoleculeViewerProps) => {
  const viewerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);

  useEffect(() => {
    // Load 3Dmol.js from CDN if it's not already loaded
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = 'https://3dmol.org/build/3Dmol-min.js';
      script.async = true;
      script.onload = () => {
        console.log('3Dmol.js loaded successfully');
        initViewer();
      };
      script.onerror = () => {
        console.error('Failed to load 3Dmol.js');
      };
      document.body.appendChild(script);
    } else {
      initViewer();
    }

    return () => {
      // Clean up viewer when component unmounts
      if (viewerInstanceRef.current) {
        try {
          viewerInstanceRef.current.clear();
        } catch (e) {
          console.error('Error cleaning up viewer:', e);
        }
      }
    };
  }, []);

  useEffect(() => {
    // Update the viewer when PDB data changes
    if (pdbData && viewerInstanceRef.current && window.$3Dmol) {
      try {
        updateMolecule(pdbData);
      } catch (e) {
        console.error('Error updating molecule:', e);
      }
    }
  }, [pdbData]);

  const initViewer = () => {
    if (!viewerRef.current || !window.$3Dmol) return;

    try {
      console.log('Initializing 3Dmol viewer');
      // Clear the container first
      viewerRef.current.innerHTML = '';
      
      // Create a new viewer instance
      const viewer = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'white',
        antialias: true,
        id: 'molecule-viewer'
      });
      
      viewerInstanceRef.current = viewer;
      
      // Set initial view settings
      viewer.setStyle({}, { stick: {} });
      viewer.zoomTo();
      viewer.render();
      
      console.log('3Dmol viewer initialized successfully');
      
      // If we already have PDB data, load it
      if (pdbData) {
        updateMolecule(pdbData);
      }
    } catch (e) {
      console.error('Error initializing viewer:', e);
    }
  };

  const updateMolecule = (pdbData: string) => {
    if (!viewerInstanceRef.current || !window.$3Dmol) return;
    
    try {
      const viewer = viewerInstanceRef.current;
      
      // Clear any existing molecules
      viewer.clear();
      
      // Add the new molecule from PDB data
      viewer.addModel(pdbData, 'pdb');
      
      // Set style and coloring
      viewer.setStyle({}, { stick: {} });
      viewer.setViewStyle({ style: 'outline' });
      
      // Add surface representation
      viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
        opacity: 0.7,
        color: 'lightblue'
      });
      
      // Center and zoom to fit the molecule
      viewer.zoomTo();
      
      // Render the scene
      viewer.render();
      
      console.log('Molecule updated successfully');
    } catch (e) {
      console.error('Error updating molecule:', e);
    }
  };

  return (
    <div className="molecule-viewer-container">
      <h2>3D Molecule Viewer</h2>
      {isLoading ? (
        <div className="loading-indicator">Loading...</div>
      ) : (
        <>
          {!pdbData ? (
            <div className="empty-state">
              Enter a SMILES string and click "Process" to view the 3D structure
            </div>
          ) : (
            <div className="viewer-wrapper">
              <div ref={viewerRef} className="viewer" id="molecule-viewer"></div>
              <div className="viewer-controls">
                <button onClick={() => viewerInstanceRef.current?.rotate(10, 'y')}>
                  Rotate Y
                </button>
                <button onClick={() => viewerInstanceRef.current?.rotate(10, 'x')}>
                  Rotate X
                </button>
                <button onClick={() => viewerInstanceRef.current?.zoomTo()}>
                  Reset View
                </button>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
};

export default MoleculeViewer; 