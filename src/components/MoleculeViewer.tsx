import { useEffect, useRef, useState } from 'react';
import '../styles/MoleculeViewer.css';

// No need to redeclare the Window interface here as it's defined in global.d.ts

interface MoleculeViewerProps {
  pdbData: string;
  isLoading: boolean;
}

const MoleculeViewer = ({ pdbData, isLoading }: MoleculeViewerProps) => {
  const viewerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);
  const [viewerError, setViewerError] = useState<string | null>(null);
  const [initAttempts, setInitAttempts] = useState(0);
  const [scriptsLoaded, setScriptsLoaded] = useState(false);

  // Ensure scripts are loaded
  useEffect(() => {
    const loadScripts = async () => {
      try {
        // Check if jQuery is loaded
        if (!window.jQuery) {
          console.log('jQuery not loaded, loading it now');
          const jqueryScript = document.createElement('script');
          jqueryScript.src = 'https://code.jquery.com/jquery-3.6.0.min.js';
          jqueryScript.async = true;
          
          // Create a promise to wait for script load
          await new Promise((resolve, reject) => {
            jqueryScript.onload = resolve;
            jqueryScript.onerror = reject;
            document.head.appendChild(jqueryScript);
          });
          
          console.log('jQuery loaded successfully');
        }
        
        // Wait a moment after jQuery loads before loading 3DMol
        await new Promise(resolve => setTimeout(resolve, 500));
        
        // Check if 3DMol is loaded
        if (!window.$3Dmol) {
          console.log('3DMol.js not loaded, loading it now');
          const threeDMolScript = document.createElement('script');
          threeDMolScript.src = 'https://3dmol.org/build/3Dmol-min.js';
          threeDMolScript.async = true;
          
          // Create a promise to wait for script load
          await new Promise((resolve, reject) => {
            threeDMolScript.onload = resolve;
            threeDMolScript.onerror = reject;
            document.head.appendChild(threeDMolScript);
          });
          
          console.log('3DMol.js loaded successfully');
        }
        
        // Wait a moment after 3DMol loads to ensure it's fully initialized
        await new Promise(resolve => setTimeout(resolve, 1000));
        setScriptsLoaded(true);
      } catch (error) {
        console.error('Error loading scripts:', error);
        setViewerError(`Error loading required scripts: ${error}`);
      }
    };
    
    loadScripts();
    
    // Cleanup function
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

  // Initialize the viewer when scripts are loaded
  useEffect(() => {
    if (!scriptsLoaded) return;
    
    console.log('Scripts loaded, attempting to initialize viewer');
    
    // Wait for DOM to be fully rendered
    const initTimer = setTimeout(() => {
      if (viewerRef.current && window.$3Dmol) {
        console.log('DOM and 3DMol.js available, initializing viewer');
        initViewer();
      } else {
        console.error('Cannot initialize viewer: DOM element or 3DMol.js not available, will retry');
        if (initAttempts < 5) {
          setInitAttempts(prev => prev + 1);
        } else {
          setViewerError('Cannot initialize viewer: DOM element or 3DMol.js not available after multiple attempts');
        }
      }
    }, 2000); // Increased delay to 2 seconds

    return () => {
      clearTimeout(initTimer);
    };
  }, [initAttempts, scriptsLoaded]);

  // This effect runs when pdbData changes
  useEffect(() => {
    if (!pdbData || !viewerInstanceRef.current || !window.$3Dmol) return;
    
    console.log('PDB data changed, updating molecule');
    try {
      updateMolecule(pdbData);
    } catch (e) {
      console.error('Error updating molecule:', e);
      setViewerError(`Error updating molecule: ${e}`);
    }
  }, [pdbData]);

  const initViewer = () => {
    if (!viewerRef.current) {
      setViewerError('Cannot initialize viewer: DOM element not available');
      return;
    }

    if (!window.$3Dmol) {
      setViewerError('Cannot initialize viewer: 3DMol.js not available');
      return;
    }

    try {
      console.log('Initializing 3DMol viewer');
      
      // Make sure the container is empty
      while (viewerRef.current.firstChild) {
        viewerRef.current.removeChild(viewerRef.current.firstChild);
      }
      
      // Set explicit dimensions to ensure the viewer is visible
      viewerRef.current.style.width = '100%';
      viewerRef.current.style.height = '400px'; // Set a fixed height
      
      // Force a reflow to ensure dimensions are applied
      void viewerRef.current.offsetHeight;
      
      // Get actual dimensions for logging
      const actualWidth = viewerRef.current.clientWidth;
      const actualHeight = viewerRef.current.clientHeight;
      console.log(`Viewer container dimensions: ${actualWidth}x${actualHeight}`);
      
      if (actualWidth === 0 || actualHeight === 0) {
        console.warn('Viewer container has zero dimensions, this may cause issues');
        // Set minimum dimensions if container has zero size
        viewerRef.current.style.width = '800px';
        viewerRef.current.style.height = '600px';
      }
      
      // Create a new viewer instance with explicit config
      const config = {
        backgroundColor: '#1e1e2e', // Dark background like Blender
        antialias: true,
        id: 'molecule-viewer-' + Date.now(), // Unique ID to avoid conflicts
      };
      
      console.log('Creating viewer with config:', config);
      
      // Create the viewer with explicit dimensions
      const viewer = window.$3Dmol.createViewer(
        window.jQuery ? window.jQuery(viewerRef.current) : viewerRef.current,
        {
          ...config,
          width: actualWidth || 800,
          height: actualHeight || 600
        }
      );
      
      if (!viewer) {
        setViewerError('Failed to create 3DMol viewer');
        return;
      }
      
      viewerInstanceRef.current = viewer;
      
      // Set initial view settings with more professional styling
      viewer.setStyle({}, { 
        stick: { radius: 0.15, colorscheme: 'cyanCarbon' },
        sphere: { scale: 0.25 } 
      });
      viewer.setViewStyle({ style: "outline" });
      viewer.zoomTo();
      viewer.render();
      
      console.log('3DMol viewer initialized successfully');
      setViewerError(null);
      
      // If we already have PDB data, load it
      if (pdbData) {
        console.log('Initial PDB data available, loading molecule');
        updateMolecule(pdbData);
      }
    } catch (e) {
      const error = `Error initializing viewer: ${e}`;
      console.error(error);
      setViewerError(error);
    }
  };

  const updateMolecule = (pdbData: string) => {
    if (!viewerInstanceRef.current) {
      console.error('Cannot update molecule: viewer not available');
      return;
    }
    
    if (!window.$3Dmol) {
      console.error('Cannot update molecule: 3DMol.js not available');
      return;
    }
    
    try {
      const viewer = viewerInstanceRef.current;
      
      console.log('Clearing existing molecules');
      viewer.clear();
      
      console.log('Adding new molecule from PDB data');
      viewer.addModel(pdbData, 'pdb');
      
      console.log('Setting molecule style');
      viewer.setStyle({}, { 
        stick: { radius: 0.15, colorscheme: 'cyanCarbon' },
        sphere: { scale: 0.25 } 
      });
      
      // Add surface representation with a semi-transparent surface
      console.log('Adding surface representation');
      if (window.$3Dmol.SurfaceType) {
        viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
          opacity: 0.5,
          color: 'lightblue'
        });
      }
      
      console.log('Zooming to fit molecule');
      viewer.zoomTo();
      
      console.log('Rendering scene');
      viewer.render();
      
      console.log('Molecule updated successfully');
    } catch (e) {
      console.error('Error updating molecule:', e);
      setViewerError(`Error updating molecule: ${e}`);
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
              <div 
                ref={viewerRef} 
                className="viewer" 
                id="molecule-viewer"
                style={{ width: '100%', height: '400px' }} // Set fixed height here too
              ></div>
              
              {viewerError && (
                <div className="error-message">
                  Error: {viewerError}
                </div>
              )}
              
              <div className="viewer-controls">
                <button onClick={() => {
                  if (viewerInstanceRef.current) {
                    viewerInstanceRef.current.rotate(10, 'y');
                    viewerInstanceRef.current.render();
                  }
                }}>
                  <span className="icon">⟲</span>
                  <span className="label">Y</span>
                </button>
                <button onClick={() => {
                  if (viewerInstanceRef.current) {
                    viewerInstanceRef.current.rotate(10, 'x');
                    viewerInstanceRef.current.render();
                  }
                }}>
                  <span className="icon">⟲</span>
                  <span className="label">X</span>
                </button>
                <button onClick={() => {
                  if (viewerInstanceRef.current) {
                    viewerInstanceRef.current.zoomTo();
                    viewerInstanceRef.current.render();
                  }
                }}>
                  <span className="icon">⊕</span>
                  <span className="label">Reset</span>
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