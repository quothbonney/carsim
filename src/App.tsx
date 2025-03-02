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
    try {
      const result = await invoke<any>("process_smiles", { smiles: smilesInput });
      setMoleculeData(result);
      console.log("Processed molecule data:", result);
    } catch (err) {
      setError(`Failed to process molecule: ${err}`);
      console.error("Error processing molecule:", err);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="app-container">
      <header className="app-header">
        <h1>Molecular Viewer</h1>
        <div className="service-status">
          Service Status: <span className={`status-${serviceStatus.toLowerCase()}`}>{serviceStatus}</span>
          {serviceStatus !== "Running" && (
            <button onClick={startPythonService} disabled={isLoading}>
              {isLoading ? "Starting..." : "Start Service"}
            </button>
          )}
        </div>
      </header>

      <main className="app-content">
        <div className="left-panel">
          <MoleculeControls
            smiles={smiles}
            setSmiles={setSmiles}
            onProcess={processMolecule}
            isLoading={isLoading}
          />
          {moleculeData && (
            <MoleculeInfo
              properties={moleculeData.properties}
              imageData={moleculeData.image_data}
            />
          )}
        </div>

        <div className="right-panel">
          <MoleculeViewer
            pdbData={moleculeData?.pdb_string || ""}
            isLoading={isLoading}
          />
        </div>
      </main>

      {error && <div className="error-message">{error}</div>}
    </div>
  );
}

export default App;
