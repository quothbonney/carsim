import { useEffect, useRef, useState } from 'react';
import '../styles/MoleculeViewer.css';

interface MoleculeViewerProps {
  pdbData: string;
  isLoading: boolean;
}

// A robust implementation with support for protein visualization and molecular dynamics
const MoleculeViewer = ({ pdbData, isLoading }: MoleculeViewerProps) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const [error, setError] = useState<string | null>(null);
  const [isInitialized, setIsInitialized] = useState(false);
  const [forceRender, setForceRender] = useState(0);
  const viewerRef = useRef<any>(null);
  const lastPdbRef = useRef<string>('');
  const containerId = useRef(`viewer-container-${Date.now()}`);
  const [showSurface, setShowSurface] = useState(false);
  const modelRef = useRef<any>(null); // Store molecular model for future manipulations
  
  // Protein visualization options
  const [viewMode, setViewMode] = useState<'stick' | 'cartoon' | 'sphere' | 'line'>('stick');
  const [colorScheme, setColorScheme] = useState<string>('chainHetatm');
  
  // PDB ID fetching
  const [proteinAnalysis, setProteinAnalysis] = useState<any>(null);
  
  // Show protein analysis panel
  const [showAnalysis, setShowAnalysis] = useState(false);

  // Debug state to track molecule changes
  const [moleculeCount, setMoleculeCount] = useState(0);
  
  // Load 3DMol.js script once on mount
  useEffect(() => {
    const loadScript = () => {
      return new Promise<void>((resolve, reject) => {
        if (typeof window.$3Dmol !== 'undefined') {
          console.log('3DMol.js already loaded');
          resolve();
          return;
        }

        console.log('Loading 3DMol.js script');
        const script = document.createElement('script');
        script.src = 'https://3dmol.org/build/3Dmol-min.js';
        script.async = true;
        script.onload = () => {
          console.log('3DMol.js loaded successfully');
          resolve();
        };
        script.onerror = () => {
          console.error('Failed to load 3DMol.js');
          reject(new Error('Failed to load 3DMol.js'));
        };
        document.head.appendChild(script);
      });
    };

    // Setup container on mount
    if (containerRef.current) {
      // Prepare container
      containerRef.current.innerHTML = '';
      const viewerDiv = document.createElement('div');
      viewerDiv.id = containerId.current;
      viewerDiv.style.width = '100%';
      viewerDiv.style.height = '100%';
      viewerDiv.style.position = 'absolute';
      viewerDiv.style.top = '0';
      viewerDiv.style.left = '0';
      containerRef.current.appendChild(viewerDiv);
      
      console.log('Viewer container created with ID:', containerId.current);
    }

    // Load 3DMol.js and initialize viewer
    loadScript()
      .then(() => {
        if (!containerRef.current) {
          console.error('Container ref is null');
          return;
        }
        
        try {
          console.log('Initializing 3DMol viewer...');
          const viewerElement = document.getElementById(containerId.current);
          
          if (!viewerElement) {
            console.error('Viewer element not found');
            setError('Failed to initialize: viewer element not found');
            return;
          }
          
          console.log('Viewer element found:', viewerElement);
          
          const config = { backgroundColor: 'black' };
          const viewer = window.$3Dmol.createViewer(
            viewerElement,
            config
          );
          
          if (viewer) {
            viewerRef.current = viewer;
            setIsInitialized(true);
            console.log('3DMol viewer initialized successfully');
            
            // If we have data, display it
            if (pdbData) {
              console.log('Initial PDB data available, displaying molecule');
              displayMolecule(pdbData, viewer);
            }
          } else {
            console.error('Failed to create viewer object');
            setError('Failed to create 3D viewer');
          }
        } catch (err) {
          console.error('Failed to initialize 3DMol viewer:', err);
          setError('Failed to initialize molecular viewer');
        }
      })
      .catch(err => {
        console.error('Error in 3DMol.js setup:', err);
        setError('Failed to load molecular viewer library');
      });
    
    // Cleanup function
    return () => {
      if (viewerRef.current) {
        try {
          // Best effort to clean up viewer
          viewerRef.current = null;
        } catch (err) {
          console.error('Error cleaning up 3DMol viewer:', err);
        }
      }
    };
  }, []);
  
  // Re-setup the viewer if it's reloaded
  useEffect(() => {
    if (forceRender > 0 && containerRef.current && typeof window.$3Dmol !== 'undefined') {
      try {
        console.log(`Reinitializing viewer (force render: ${forceRender})`);
        const viewerElement = document.getElementById(containerId.current);
        
        if (!viewerElement) {
          console.error('Viewer element not found during reload');
          setError('Failed to reinitialize: viewer element not found');
          return;
        }
        
        const config = { backgroundColor: 'black' };
        const viewer = window.$3Dmol.createViewer(
          viewerElement,
          config
        );
        
        if (viewer) {
          viewerRef.current = viewer;
          setIsInitialized(true);
          modelRef.current = null; // Clear model reference
          console.log('3DMol viewer reinitialized');
          
          // If we have data, redisplay it
          if (pdbData) {
            console.log('Re-displaying molecule after viewer reload');
            displayMolecule(pdbData, viewer);
          }
        }
      } catch (err) {
        console.error('Failed to reinitialize 3DMol viewer:', err);
        setError('Failed to reinitialize molecular viewer');
      }
    }
  }, [forceRender, pdbData]);

  // Handle PDB data changes
  useEffect(() => {
    if (!isInitialized || !viewerRef.current) return;
    
    console.log('PDB data changed check, current length:', pdbData?.length);
    console.log('Last PDB ref length:', lastPdbRef.current?.length);
    
    // Skip if same PDB as before or no data
    if (pdbData === lastPdbRef.current) {
      console.log('Same PDB data, skipping render');
      return;
    }
    
    // Update last PDB reference
    lastPdbRef.current = pdbData;
    
    console.log('PDB data changed, updating viewer');
    setMoleculeCount(prev => prev + 1);
    
    // Reset errors
    setError(null);
    
    // Clear and reset the viewer on data change
    if (viewerRef.current) {
      try {
        console.log('Clearing viewer');
        viewerRef.current.clear();
        modelRef.current = null;
      } catch (e) {
        console.error('Error clearing viewer:', e);
      }
    }
    
    if (!pdbData) {
      console.log('No PDB data, showing empty viewer');
      if (viewerRef.current) {
        viewerRef.current.render();
      }
      return;
    }
    
    try {
      // Display the molecule
      console.log('Displaying new molecule');
      displayMolecule(pdbData, viewerRef.current);
      
      // Analyze the PDB structure
      analyzePdbData(pdbData).catch(err => {
        console.error('Error analyzing PDB data:', err);
      });
    } catch (err) {
      console.error('Error displaying PDB:', err);
      setError(`Failed to display molecule: ${err}`);
    }
  }, [pdbData, isInitialized]);

  // Re-apply style when view mode changes
  useEffect(() => {
    if (isInitialized && viewerRef.current && modelRef.current) {
      console.log('View mode changed, applying new style');
      applyCurrentStyle();
    }
  }, [viewMode, colorScheme, isInitialized]);

  // Update surface visibility when showSurface changes
  useEffect(() => {
    if (isInitialized && viewerRef.current && modelRef.current) {
      console.log('Surface visibility changed');
      updateSurfaceVisibility();
    }
  }, [showSurface, isInitialized]);

  // Analyze PDB data for additional information
  const analyzePdbData = async (pdbContent: string) => {
    try {
      const response = await fetch('http://localhost:5000/analyze-pdb', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ pdbContent }),
      });
      
      const data = await response.json();
      
      if (data.error) {
        console.error('Error from server:', data.error);
        return;
      }
      
      setProteinAnalysis(data);
      console.log('Analysis data:', data);
    } catch (err) {
      console.error('Failed to analyze PDB:', err);
    }
  };

  // Handle displaying the molecule in the viewer
  const displayMolecule = (data: string, viewer: any) => {
    try {
      console.log('Displaying molecule, data length:', data.length);
      
      // Make sure the viewer is fully cleared
      viewer.clear();
      viewer.removeAllLabels();
      viewer.removeAllShapes();
      viewer.removeAllSurfaces();
      
      // First try to use the minimal rendering for large proteins
      const didMinimalRendering = tryMinimalRendering(data, viewer);
      
      if (!didMinimalRendering) {
        // For smaller molecules, or if minimal rendering failed, use full rendering
        console.log('Using standard rendering');
        
        // Parse the data explicitly for better error handling
        try {
          const model = viewer.addModel(data, "pdb");
          if (!model) {
            console.error('Failed to add model to viewer');
            setError('Failed to create molecular model');
            return;
          }
          
          modelRef.current = model;
          
          applyCurrentStyle();
          
          // Center and zoom
          viewer.zoomTo();
          viewer.render();
          console.log('Model rendered successfully');
        } catch (modelErr) {
          console.error('Error creating model:', modelErr);
          setError(`Error creating model: ${modelErr}`);
          
          // Try one more time with a simple model
          viewer.clear();
          const simpleModel = viewer.addModel(data, "pdb", {keepH: false, parseSettings: {singleModel: true}});
          if (simpleModel) {
            modelRef.current = simpleModel;
            simpleModel.setStyle({}, {line:{}});
            viewer.zoomTo();
            viewer.render();
            console.log('Used simple fallback model');
          } else {
            throw new Error('Failed to create model after multiple attempts');
          }
        }
      }
    } catch (err) {
      console.error('Error in displayMolecule:', err);
      setError(`Failed to display molecule: ${err}`);
    }
  };
  
  // Apply the current style settings to the molecule
  const applyCurrentStyle = () => {
    if (!viewerRef.current || !modelRef.current) return;
    
    try {
      const viewer = viewerRef.current;
      const model = modelRef.current;
      
      // Clear old styles
      model.setStyle({}, {});
      
      // Apply new style based on view mode
      switch (viewMode) {
        case 'stick':
          model.setStyle({}, { stick: {} });
          break;
        case 'cartoon':
          model.setStyle({}, { cartoon: { color: colorScheme } });
          // Also show ligands with sticks if in cartoon mode
          model.setStyle({hetflag: true}, { stick: {} });
          break;
        case 'sphere':
          model.setStyle({}, { sphere: {} });
          break;
        case 'line':
          model.setStyle({}, { line: {} });
          break;
      }
      
      // Update surface if needed
      updateSurfaceVisibility();
      
      // Re-render
      viewer.render();
    } catch (err) {
      console.error('Error applying style:', err);
      setError(`Failed to apply visualization style: ${err}`);
    }
  };
  
  // For large proteins, use a more efficient rendering approach
  const tryMinimalRendering = (data: string, viewer: any) => {
    if (data.length > 200000) { // If PDB is large (lowered threshold)
      try {
        console.log('Using minimal rendering for large molecule');
        const model = viewer.addModel(data, "pdb", { keepH: false });
        if (!model) {
          console.error('Failed to add model in minimal rendering');
          return false;
        }
        
        modelRef.current = model;
        
        // For large proteins, default to cartoon representation
        model.setStyle({}, { cartoon: { color: colorScheme } });
        setViewMode('cartoon');
        
        // Add stick representation for ligands and important residues
        model.setStyle({hetflag: true}, { stick: {} });
        
        viewer.zoomTo();
        viewer.render();
        console.log('Minimal rendering successful');
        return true;
      } catch (err) {
        console.error('Minimal rendering failed, will try normal rendering:', err);
        return false;
      }
    }
    return false;
  };
  
  // Toggle surface visibility
  const updateSurfaceVisibility = () => {
    if (!viewerRef.current || !modelRef.current) return;
    
    try {
      const viewer = viewerRef.current;
      
      // Always remove all surfaces first to prevent stacking
      viewer.removeAllSurfaces();
      
      if (showSurface) {
        console.log('Adding surface');
        const model = modelRef.current;
        model.addSurface(window.$3Dmol.SurfaceType.VDW, {
          opacity: 0.7,
          color: colorScheme === 'spectrum' ? 'spectrum' : 'white'
        });
      }
      
      viewer.render();
    } catch (err) {
      console.error('Error updating surface:', err);
    }
  };
  
  // Toggle surface visibility
  const toggleSurface = () => {
    setShowSurface(!showSurface);
  };
  
  // Rotate around specified bond
  const rotateBond = (atomIdx1: number, atomIdx2: number, angle: number) => {
    if (viewerRef.current && modelRef.current) {
      modelRef.current.rotateSelected(atomIdx1, atomIdx2, angle, true);
      viewerRef.current.render();
    }
  };
  
  // Rotate around Y axis
  const rotateY = () => {
    if (viewerRef.current) {
      viewerRef.current.rotate(0, 1, 0, 0.5);
      viewerRef.current.render();
    }
  };
  
  // Rotate around X axis
  const rotateX = () => {
    if (viewerRef.current) {
      viewerRef.current.rotate(1, 0, 0, 0.5);
      viewerRef.current.render();
    }
  };
  
  // Reset the view to default
  const resetView = () => {
    if (viewerRef.current) {
      viewerRef.current.setView([]);
      viewerRef.current.zoomTo();
      viewerRef.current.render();
    }
  };
  
  // Toggle protein analysis panel
  const toggleAnalysis = () => {
    setShowAnalysis(!showAnalysis);
  };
  
  // Manual reload function
  const reloadViewer = () => {
    // Generate a new container ID to force DOM recreation
    containerId.current = `viewer-container-${Date.now()}`;
    
    if (containerRef.current) {
      containerRef.current.innerHTML = '';
      const viewerDiv = document.createElement('div');
      viewerDiv.id = containerId.current;
      viewerDiv.style.width = '100%';
      viewerDiv.style.height = '100%';
      viewerDiv.style.position = 'absolute';
      viewerDiv.style.top = '0';
      viewerDiv.style.left = '0';
      containerRef.current.appendChild(viewerDiv);
    }
    
    // Reset state for complete reinitialization
    modelRef.current = null;
    setForceRender(prev => prev + 1);
  };

  // Change view mode
  const changeViewMode = (mode: 'stick' | 'cartoon' | 'sphere' | 'line') => {
    setViewMode(mode);
  };

  // Change color scheme
  const changeColorScheme = (scheme: string) => {
    setColorScheme(scheme);
  };
  
  return (
    <div className="molecule-viewer-container">
      {proteinAnalysis && (
        <div className={`analysis-panel ${showAnalysis ? 'active' : ''}`}>
          <div className="analysis-header">
            <h3>Protein Analysis</h3>
            <button className="close-button" onClick={toggleAnalysis}>×</button>
          </div>
          <div className="analysis-content">
            {proteinAnalysis.residueCount && (
              <div className="analysis-item">
                <strong>Total residues:</strong> {proteinAnalysis.residueCount}
              </div>
            )}
            
            {proteinAnalysis.secondaryStructure && (
              <div className="analysis-item">
                <strong>Secondary Structure:</strong>
                <div className="secondary-structure">
                  <div className="ss-bar">
                    {proteinAnalysis.secondaryStructure.helix > 0 && (
                      <div 
                        className="ss-helix" 
                        style={{ 
                          width: `${(proteinAnalysis.secondaryStructure.helix / proteinAnalysis.residueCount) * 100}%`
                        }}
                        title={`Helix: ${proteinAnalysis.secondaryStructure.helix} residues`}
                      ></div>
                    )}
                    {proteinAnalysis.secondaryStructure.sheet > 0 && (
                      <div 
                        className="ss-sheet" 
                        style={{ 
                          width: `${(proteinAnalysis.secondaryStructure.sheet / proteinAnalysis.residueCount) * 100}%`
                        }}
                        title={`Sheet: ${proteinAnalysis.secondaryStructure.sheet} residues`}
                      ></div>
                    )}
                    {proteinAnalysis.secondaryStructure.loop > 0 && (
                      <div 
                        className="ss-loop" 
                        style={{ 
                          width: `${(proteinAnalysis.secondaryStructure.loop / proteinAnalysis.residueCount) * 100}%`
                        }}
                        title={`Loop: ${proteinAnalysis.secondaryStructure.loop} residues`}
                      ></div>
                    )}
                  </div>
                  <div className="ss-legend">
                    <div className="ss-legend-item">
                      <div className="ss-helix-icon"></div> Helix: {proteinAnalysis.secondaryStructure.helix}
                    </div>
                    <div className="ss-legend-item">
                      <div className="ss-sheet-icon"></div> Sheet: {proteinAnalysis.secondaryStructure.sheet}
                    </div>
                    <div className="ss-legend-item">
                      <div className="ss-loop-icon"></div> Loop: {proteinAnalysis.secondaryStructure.loop}
                    </div>
                  </div>
                </div>
              </div>
            )}
            
            {proteinAnalysis.chains && proteinAnalysis.chains.length > 0 && (
              <div className="analysis-item">
                <strong>Chains:</strong>
                <ul className="chains-list">
                  {proteinAnalysis.chains.map((chain: any, index: number) => (
                    <li key={index}>
                      Chain {chain.chainId}: {chain.residueCount} residues
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        </div>
      )}
      
      {isLoading ? (
        <div className="loading-indicator">
          Loading...
        </div>
      ) : (
        <>
          {!pdbData ? (
            <div className="empty-state">
              Enter a SMILES string or upload a PDB file to view the 3D structure
            </div>
          ) : (
            <div className="viewer-wrapper">
              <div 
                ref={containerRef} 
                className="viewer-container"
                data-molecule-count={moleculeCount}
              ></div>
              
              {error && (
                <div className="error-message">
                  Error: {error}
                  <button onClick={reloadViewer} className="reload-button">
                    Reload Viewer
                  </button>
                </div>
              )}
              
              <div className="viewer-controls fixed-bottom">
                {/* View mode controls */}
                <div className="control-group">
                  <button 
                    onClick={() => changeViewMode('stick')} 
                    disabled={!isInitialized}
                    className={viewMode === 'stick' ? "view-button active" : "view-button"}
                  >
                    <span className="icon">⊕</span>
                    <span className="label">Stick</span>
                  </button>
                  <button 
                    onClick={() => changeViewMode('cartoon')} 
                    disabled={!isInitialized}
                    className={viewMode === 'cartoon' ? "view-button active" : "view-button"}
                  >
                    <span className="icon">⊗</span>
                    <span className="label">Cartoon</span>
                  </button>
                  <button 
                    onClick={() => changeViewMode('sphere')} 
                    disabled={!isInitialized}
                    className={viewMode === 'sphere' ? "view-button active" : "view-button"}
                  >
                    <span className="icon">◯</span>
                    <span className="label">Sphere</span>
                  </button>
                </div>
                
                {/* Rotation controls */}
                <div className="control-group">
                  <button onClick={rotateY} disabled={!isInitialized}>
                    <span className="icon">⟲</span>
                    <span className="label">Y</span>
                  </button>
                  <button onClick={rotateX} disabled={!isInitialized}>
                    <span className="icon">⟲</span>
                    <span className="label">X</span>
                  </button>
                  <button onClick={resetView} disabled={!isInitialized}>
                    <span className="icon">⊕</span>
                    <span className="label">Reset</span>
                  </button>
                </div>
                
                {/* Surface and reload controls */}
                <div className="control-group">
                  <button 
                    onClick={toggleSurface} 
                    disabled={!isInitialized} 
                    className={showSurface ? "surface-button active" : "surface-button"}
                  >
                    <span className="icon">◯</span>
                    <span className="label">Surface</span>
                  </button>
                  <button onClick={reloadViewer} className="reload-button">
                    <span className="icon">⟳</span>
                    <span className="label">Reload</span>
                  </button>
                </div>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
};

export default MoleculeViewer;