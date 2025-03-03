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
  
  // Local PDB file loading
  const [localPdbData, setLocalPdbData] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // PDB ID fetching
  const [pdbId, setPdbId] = useState<string>('');
  const [isFetchingPdb, setIsFetchingPdb] = useState(false);
  const [proteinAnalysis, setProteinAnalysis] = useState<any>(null);
  
  // Show protein analysis panel
  const [showAnalysis, setShowAnalysis] = useState(false);
  
  // Use either provided pdbData or localPdbData
  const effectivePdbData = localPdbData || pdbData;

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
      viewerDiv.style.height = '400px';
      viewerDiv.style.position = 'relative';
      containerRef.current.appendChild(viewerDiv);

      loadScript()
        .then(() => {
          // Initialize after script is loaded
          console.log('Script loaded, initializing viewer');
          setForceRender(prev => prev + 1);
        })
        .catch(err => {
          setError(err.message);
        });
    }

    return () => {
      // Cleanup on unmount
      if (viewerRef.current) {
        try {
          viewerRef.current.clear();
        } catch (e) {
          console.log('Error during cleanup:', e);
        }
      }
    };
  }, []);

  // Create viewer instance
  useEffect(() => {
    if (!containerRef.current || forceRender === 0) return;

    console.log(`Initializing viewer (force render: ${forceRender})`);
    
    // Clean up any existing viewer
    if (viewerRef.current) {
      try {
        viewerRef.current.clear();
      } catch (e) {
        console.log('Error clearing existing viewer:', e);
      }
      viewerRef.current = null;
    }

    // Create new viewer with a slight delay
    setTimeout(() => {
      try {
        if (typeof window.$3Dmol === 'undefined') {
          setError('3DMol library not available');
          return;
        }

        const element = document.getElementById(containerId.current);
        if (!element) {
          setError('Could not find viewer container element');
          return;
        }

        console.log('Creating viewer instance');
        // Create viewer with optimized settings that support future dynamics
        const viewer = window.$3Dmol.createViewer(element, {
          backgroundColor: '#1e1e2e',
          antialias: true,
          disableFog: false, // Enable fog for better depth perception (important for dynamics)
          outline: false // Disable outline for performance
        });

        if (!viewer) {
          setError('Failed to create 3DMol viewer');
          return;
        }

        viewerRef.current = viewer;
        setIsInitialized(true);
        setError(null);

        // If we have molecule data, display it immediately
        if (effectivePdbData) {
          console.log('Initial display of molecule');
          displayMolecule(effectivePdbData, viewer);
        }
      } catch (err) {
        console.error('Error initializing viewer:', err);
        setError(`Error initializing viewer: ${err}`);
      }
    }, 100);
  }, [forceRender]);

  // Handle PDB data changes
  useEffect(() => {
    if (!isInitialized || !viewerRef.current || !effectivePdbData) return;
    
    // If data changed, update the display
    if (effectivePdbData !== lastPdbRef.current) {
      console.log('PDB data changed, updating display');
      lastPdbRef.current = effectivePdbData;
      displayMolecule(effectivePdbData, viewerRef.current);
      
      // If no analysis yet, try to get it
      if (!proteinAnalysis && localPdbData) {
        analyzePdbData(localPdbData);
      }
    }
  }, [effectivePdbData, isInitialized]);

  // Analyze PDB data using the server
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
      if (!data.error) {
        setProteinAnalysis(data);
      }
    } catch (err) {
      console.error('Failed to analyze PDB:', err);
    }
  };

  // Handle surface visibility changes
  useEffect(() => {
    if (!isInitialized || !viewerRef.current || !modelRef.current) return;
    
    try {
      updateSurfaceVisibility();
    } catch (e) {
      console.error('Error updating surface:', e);
    }
  }, [showSurface, isInitialized]);

  // Handle view mode changes
  useEffect(() => {
    if (!isInitialized || !viewerRef.current || !modelRef.current) return;
    
    try {
      applyCurrentStyle();
    } catch (e) {
      console.error('Error updating view mode:', e);
    }
  }, [viewMode, colorScheme, isInitialized]);

  // Main function to display the molecule
  const displayMolecule = (data: string, viewer: any) => {
    if (!viewer) return;
    
    try {
      console.log('Displaying molecule');
      
      // Clear any existing content
      viewer.clear();
      viewer.removeAllModels();
      viewer.removeAllShapes();
      viewer.removeAllSurfaces();
      
      // Detect if the PDB data is a protein (check for CA atom patterns)
      const isProtein = data.includes("CA") && (
        data.includes("ALA") || data.includes("ARG") || data.includes("ASN") || 
        data.includes("ASP") || data.includes("CYS") || data.includes("GLN") || 
        data.includes("GLU") || data.includes("GLY") || data.includes("HIS") || 
        data.includes("ILE") || data.includes("LEU") || data.includes("LYS") || 
        data.includes("MET") || data.includes("PHE") || data.includes("PRO") || 
        data.includes("SER") || data.includes("THR") || data.includes("TRP") || 
        data.includes("TYR") || data.includes("VAL")
      );
      
      // Parse and add the model - store reference for future modifications
      const model = viewer.addModel(data, 'pdb', {
        keepH: true, // Keep hydrogen atoms (important for dynamics)
        assignBonds: true // Make sure bonds are properly detected
      });
      
      // Store model reference for future modifications
      modelRef.current = model;
      
      // If it's a protein, default to cartoon view
      if (isProtein && viewMode === 'stick') {
        setViewMode('cartoon');
      }
      
      // Apply current styling
      applyCurrentStyle();
      
      // Update surface based on current state
      updateSurfaceVisibility();
      
      // Zoom to fit the molecule and render
      viewer.zoomTo();
      viewer.render();
      
      setError(null);
    } catch (err) {
      console.error('Error displaying molecule:', err);
      setError(`Error displaying molecule: ${err}`);
      
      // If display fails, try minimal rendering
      tryMinimalRendering(data, viewer);
    }
  };

  // Apply the current style based on viewMode and colorScheme
  const applyCurrentStyle = () => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    
    // Clear previous styles
    viewer.setStyle({}, {});
    
    // Apply selected style
    switch (viewMode) {
      case 'cartoon':
        viewer.setStyle({}, { cartoon: { color: colorScheme } });
        // Also show sticks for non-protein parts (ligands, etc.)
        viewer.setStyle({hetflag: true}, { stick: { radius: 0.15, colorscheme: 'greenCarbon' } });
        break;
      case 'stick':
        viewer.setStyle({}, { 
          stick: { radius: 0.15, colorscheme: colorScheme === 'chainHetatm' ? 'cyanCarbon' : colorScheme },
          sphere: { scale: 0.25 } 
        });
        break;
      case 'sphere':
        viewer.setStyle({}, { sphere: { scale: 0.6, colorscheme: colorScheme } });
        break;
      case 'line':
        viewer.setStyle({}, { line: { colorscheme: colorScheme } });
        break;
    }
    
    viewer.render();
  };

  // Try minimal rendering if normal display fails
  const tryMinimalRendering = (data: string, viewer: any) => {
    setTimeout(() => {
      try {
        if (viewer) {
          viewer.clear();
          const model = viewer.addModel(data, 'pdb');
          modelRef.current = model;
          viewer.setStyle({}, { line: { width: 1.0 } });
          viewer.zoomTo();
          viewer.render();
        }
      } catch (retryErr) {
        console.error('Retry failed:', retryErr);
        setForceRender(prev => prev + 1);
      }
    }, 100);
  };

  // Update surface visibility based on showSurface state
  const updateSurfaceVisibility = () => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    
    // Remove any existing surfaces
    viewer.removeAllSurfaces();
    
    // Add surface if enabled
    if (showSurface) {
      try {
        if (window.$3Dmol && window.$3Dmol.SurfaceType) {
          viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
            opacity: 0.5,
            color: 'lightblue'
          });
        }
      } catch (e) {
        console.warn('Error adding surface:', e);
      }
    }
    
    // Re-render to apply changes
    viewer.render();
  };

  // Toggle surface visibility
  const toggleSurface = () => {
    setShowSurface(prev => !prev);
  };

  // Example of a function that could modify molecule structure
  const rotateBond = (atomIdx1: number, atomIdx2: number, angle: number) => {
    // This would be implemented to allow bond rotation
    // 3DMol.js allows for atom position manipulation
    console.log(`Future feature: Rotate bond between atoms ${atomIdx1}-${atomIdx2} by ${angle} degrees`);
  };

  // Handler functions for viewer controls
  const rotateY = () => {
    if (viewerRef.current) {
      try {
        viewerRef.current.rotate(10, 'y');
        viewerRef.current.render();
      } catch (e) {
        console.error('Error rotating molecule:', e);
      }
    }
  };

  const rotateX = () => {
    if (viewerRef.current) {
      try {
        viewerRef.current.rotate(10, 'x');
        viewerRef.current.render();
      } catch (e) {
        console.error('Error rotating molecule:', e);
      }
    }
  };

  const resetView = () => {
    if (viewerRef.current) {
      try {
        viewerRef.current.zoomTo();
        viewerRef.current.render();
      } catch (e) {
        console.error('Error resetting view:', e);
      }
    }
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
      viewerDiv.style.height = '400px';
      containerRef.current.appendChild(viewerDiv);
    }
    
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

  // Handle PDB file upload
  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      if (content) {
        setLocalPdbData(content);
      }
    };
    reader.onerror = () => {
      setError("Failed to read the uploaded file");
    };
    reader.readAsText(file);
  };

  // Trigger file input click
  const triggerFileUpload = () => {
    if (fileInputRef.current) {
      fileInputRef.current.click();
    }
  };

  // Handle PDB ID fetch
  const fetchPdbById = async () => {
    if (!pdbId || pdbId.length !== 4) {
      setError('PDB ID must be 4 characters long');
      return;
    }
    
    try {
      setIsFetchingPdb(true);
      setError(null);
      
      const response = await fetch(`http://localhost:5000/fetch-pdb/${pdbId}`);
      const data = await response.json();
      
      if (data.error) {
        setError(data.error);
        return;
      }
      
      if (data.pdbContent) {
        setLocalPdbData(data.pdbContent);
        setProteinAnalysis(data.analysis);
      } else {
        setError('Failed to fetch PDB structure');
      }
    } catch (err) {
      setError(`Failed to fetch PDB: ${err}`);
    } finally {
      setIsFetchingPdb(false);
    }
  };

  // Handle PDB ID input
  const handlePdbIdChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setPdbId(e.target.value.toUpperCase());
  };

  // Toggle analysis panel
  const toggleAnalysis = () => {
    setShowAnalysis(prev => !prev);
  };

  return (
    <div className="molecule-viewer-container">
      <h2>CarsimMD</h2>
      <div className="file-upload-bar">
        <button onClick={triggerFileUpload} className="upload-button">
          <span className="icon">📂</span>
          <span className="label">Upload PDB</span>
        </button>
        <input 
          type="file" 
          ref={fileInputRef}
          accept=".pdb" 
          onChange={handleFileUpload}
          style={{ display: 'none' }}
        />
        
        <div className="pdb-id-input-container">
          <input 
            type="text" 
            value={pdbId} 
            onChange={handlePdbIdChange}
            placeholder="PDB ID (e.g. 1CRN)" 
            maxLength={4}
            className="pdb-id-input"
          />
          <button 
            onClick={fetchPdbById} 
            disabled={isFetchingPdb || pdbId.length !== 4}
            className="fetch-button"
          >
            <span className="icon">🔍</span>
            <span className="label">Fetch</span>
          </button>
        </div>
        
        {localPdbData && (
          <>
            <div className="file-name">
              {proteinAnalysis ? (
                `PDB loaded: ${proteinAnalysis.totalResidues} residues, ${proteinAnalysis.totalAtoms} atoms`
              ) : (
                'PDB file loaded'
              )}
            </div>
            <button 
              onClick={() => {
                setLocalPdbData(null);
                setProteinAnalysis(null);
              }} 
              className="clear-button"
            >
              <span className="icon">✖</span>
              <span className="label">Clear</span>
            </button>
            {proteinAnalysis && (
              <button 
                onClick={toggleAnalysis} 
                className={showAnalysis ? "analysis-button active" : "analysis-button"}
              >
                <span className="icon">ℹ️</span>
                <span className="label">Analysis</span>
              </button>
            )}
          </>
        )}
      </div>
      
      {showAnalysis && proteinAnalysis && (
        <div className="analysis-panel">
          <h3>Protein Analysis</h3>
          <div className="analysis-content">
            <div className="analysis-item">
              <strong>Total Residues:</strong> {proteinAnalysis.totalResidues}
            </div>
            <div className="analysis-item">
              <strong>Total Atoms:</strong> {proteinAnalysis.totalAtoms}
            </div>
            
            {proteinAnalysis.secondaryStructure && (
              <div className="analysis-item">
                <strong>Secondary Structure:</strong>
                <div className="secondary-structure">
                  <div className="ss-bar">
                    {proteinAnalysis.secondaryStructure.helix > 0 && (
                      <div 
                        className="ss-helix" 
                        style={{ 
                          width: `${(proteinAnalysis.secondaryStructure.helix / proteinAnalysis.totalResidues) * 100}%` 
                        }}
                        title={`Helix: ${proteinAnalysis.secondaryStructure.helix} residues`}
                      ></div>
                    )}
                    {proteinAnalysis.secondaryStructure.sheet > 0 && (
                      <div 
                        className="ss-sheet" 
                        style={{ 
                          width: `${(proteinAnalysis.secondaryStructure.sheet / proteinAnalysis.totalResidues) * 100}%` 
                        }}
                        title={`Sheet: ${proteinAnalysis.secondaryStructure.sheet} residues`}
                      ></div>
                    )}
                    {proteinAnalysis.secondaryStructure.loop > 0 && (
                      <div 
                        className="ss-loop" 
                        style={{ 
                          width: `${(proteinAnalysis.secondaryStructure.loop / proteinAnalysis.totalResidues) * 100}%` 
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
      
      {isLoading || isFetchingPdb ? (
        <div className="loading-indicator">
          {isFetchingPdb ? 'Fetching PDB structure...' : 'Loading...'}
        </div>
      ) : (
        <>
          {!effectivePdbData ? (
            <div className="empty-state">
              Enter a SMILES string or upload a PDB file to view the 3D structure
            </div>
          ) : (
            <div className="viewer-wrapper">
              <div 
                ref={containerRef} 
                className="viewer-container"
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