use std::process::{Command, Stdio};
use std::io::Write;
use std::thread;
use std::sync::mpsc;
use std::time::Duration;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};

// Learn more about Tauri commands at https://tauri.app/develop/calling-rust/
#[tauri::command]
fn greet(name: &str) -> String {
    format!("Hello, {}! You've been greeted from Rust!", name)
}

#[derive(Serialize, Deserialize, Debug)]
struct PythonServiceResponse {
    pdb_string: Option<String>,
    image_data: Option<String>,
    properties: Option<MolecularProperties>,
    error: Option<String>,
    status: Option<String>,
}

#[derive(Serialize, Deserialize, Debug)]
struct MolecularProperties {
    molecular_weight: f64,
    formula: String,
    num_atoms: i32,
    num_bonds: i32,
}

#[derive(Serialize, Deserialize, Debug)]
struct SmilesRequest {
    smiles: String,
}

#[tauri::command]
fn start_python_service() -> Result<String, String> {
    // Start the Python service in a separate process
    let python_cmd = if cfg!(target_os = "windows") {
        "python"
    } else {
        "python3"
    };

    let service_path = "python_service/app.py";
    
    let child = Command::new(python_cmd)
        .arg(service_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start Python service: {}", e))?;
    
    // Wait a bit for the service to start
    thread::sleep(Duration::from_secs(2));
    
    // Check if the service is running by making a health check request
    let (tx, rx) = mpsc::channel();
    
    thread::spawn(move || {
        let client = Client::new();
        match client.get("http://localhost:5000/health").send() {
            Ok(response) => {
                if response.status().is_success() {
                    tx.send(Ok("Python service started successfully".to_string())).unwrap();
                } else {
                    tx.send(Err(format!("Service returned error status: {}", response.status()))).unwrap();
                }
            },
            Err(e) => {
                tx.send(Err(format!("Failed to connect to Python service: {}", e))).unwrap();
            }
        }
    });
    
    // Wait for the health check result with a timeout
    match rx.recv_timeout(Duration::from_secs(5)) {
        Ok(result) => result,
        Err(_) => Err("Timeout waiting for Python service to start".to_string()),
    }
}

#[tauri::command]
fn process_smiles(smiles: String) -> Result<PythonServiceResponse, String> {
    let client = Client::new();
    let request_body = SmilesRequest { smiles };
    
    match client.post("http://localhost:5000/process_smiles")
        .json(&request_body)
        .send() {
            Ok(response) => {
                if response.status().is_success() {
                    match response.json::<PythonServiceResponse>() {
                        Ok(data) => Ok(data),
                        Err(e) => Err(format!("Failed to parse response: {}", e)),
                    }
                } else {
                    Err(format!("Service returned error status: {}", response.status()))
                }
            },
            Err(e) => Err(format!("Failed to connect to Python service: {}", e)),
        }
}

#[tauri::command]
fn convert_to_3d(smiles: String) -> Result<String, String> {
    let client = Client::new();
    let request_body = SmilesRequest { smiles };
    
    match client.post("http://localhost:5000/convert_to_3d")
        .json(&request_body)
        .send() {
            Ok(response) => {
                if response.status().is_success() {
                    match response.json::<PythonServiceResponse>() {
                        Ok(data) => {
                            match data.pdb_string {
                                Some(pdb) => Ok(pdb),
                                None => Err("No PDB string returned".to_string()),
                            }
                        },
                        Err(e) => Err(format!("Failed to parse response: {}", e)),
                    }
                } else {
                    Err(format!("Service returned error status: {}", response.status()))
                }
            },
            Err(e) => Err(format!("Failed to connect to Python service: {}", e)),
        }
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            greet,
            start_python_service,
            process_smiles,
            convert_to_3d
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
