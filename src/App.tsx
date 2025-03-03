import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import "./App.css";
import MoleculeViewer from "./components/MoleculeViewer";
import MoleculeControls from "./components/MoleculeControls";
import MoleculeInfo from "./components/MoleculeInfo";

interface MolecularProperties {
  molecular_weight: number;
  formula: string;
  num_atoms: number;
  num_bonds: number;
}

interface MoleculeData {
  pdb_string: string;
  image_data: string;
  properties: MolecularProperties;
}

function App() {
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [serviceStatus, setServiceStatus] = useState<string>("Not started");
  const [moleculeData, setMoleculeData] = useState<MoleculeData | null>(null);
  const [smiles, setSmiles] = useState<string>("");
  const [sidebarOpen, setSidebarOpen] = useState(true);
  const [infoOpen, setInfoOpen] = useState(true);
  const [uploadedPdbData, setUploadedPdbData] = useState<string | null>(null);

  // Start the Python service when the app loads
  useEffect(() => {
    startPythonService();
  }, []);

  const startPythonService = async () => {
    setIsLoading(true);
    setError(null);
    try {
      const result = await invoke<string>("start_python_service");
      setServiceStatus("Running");
      console.log("Python service started:", result);
    } catch (err) {
      setError(`Failed to start Python service: ${err}`);
      setServiceStatus("Error");
      console.error("Error starting Python service:", err);
    } finally {
      setIsLoading(false);
    }
  };

  const processMolecule = async (smilesInput: string) => {
    if (!smilesInput.trim()) {
      setError("Please enter a SMILES string");
      return;
    }

    setIsLoading(true);
    setError(null);
    
    // Clear any existing PDB data
    setUploadedPdbData(null);
    setMoleculeData(null);
    
    try {
      const result = await invoke<any>("process_smiles", { smiles: smilesInput });
      
      // Set new data
      setMoleculeData(result);
      console.log("Processed molecule data:", result);
    } catch (err) {
      setError(`Failed to process molecule: ${err}`);
      console.error("Error processing molecule:", err);
    } finally {
      setIsLoading(false);
    }
  };

  const handlePdbUpload = (pdbData: string) => {
    console.log("PDB data received in App component, length:", pdbData.length);
    
    // First clear current state
    setUploadedPdbData(null);
    setMoleculeData(null);
    
    // Then set new data after a short delay to ensure proper rerender
    setTimeout(() => {
      setUploadedPdbData(pdbData);
    }, 10);
  };

  const toggleSidebar = () => {
    setSidebarOpen(!sidebarOpen);
  };

  const toggleInfo = () => {
    setInfoOpen(!infoOpen);
  };

  // Determine which PDB data to use - either from SMILES or uploaded
  const effectivePdbData = uploadedPdbData || (moleculeData?.pdb_string || "");

  console.log("Effective PDB data available:", 
    uploadedPdbData ? "Uploaded PDB" : 
    moleculeData?.pdb_string ? "SMILES-derived PDB" : 
    "None");

  return (
    <div className="app-container">
      <header className="app-header">
        <div className="header-left">
          <button className="sidebar-toggle" onClick={toggleSidebar}>
            {sidebarOpen ? "◀" : "▶"}
          </button>
          <h1>Molecular Viewer</h1>
        </div>
        <div className="service-status">
          Service: <span className={`status-${serviceStatus.toLowerCase()}`}>{serviceStatus}</span>
          {serviceStatus !== "Running" && (
            <button onClick={startPythonService} disabled={isLoading}>
              {isLoading ? "Starting..." : "Start Service"}
            </button>
          )}
        </div>
      </header>

      <main className="app-content">
        {sidebarOpen && (
          <div className="sidebar">
            <MoleculeControls
              smiles={smiles}
              setSmiles={setSmiles}
              onProcess={processMolecule}
              isLoading={isLoading}
              onPdbUpload={handlePdbUpload}
            />
          </div>
        )}

        <div className="main-viewport">
          <MoleculeViewer
            pdbData={effectivePdbData}
            isLoading={isLoading}
          />
        </div>

        {infoOpen && moleculeData && (
          <div className="info-panel">
            <div className="info-header">
              <h3>Molecule Info</h3>
              <button className="info-toggle" onClick={toggleInfo}>✕</button>
            </div>
            <MoleculeInfo
              properties={moleculeData.properties}
              imageData={moleculeData.image_data}
            />
          </div>
        )}

        {!infoOpen && moleculeData && (
          <button className="info-toggle-button" onClick={toggleInfo}>
            ℹ️
          </button>
        )}
      </main>

      {error && <div className="error-message">{error}</div>}
    </div>
  );
}

export default App;
