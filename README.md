# Molecular Viewer

A cross-platform desktop application for molecular visualization. This application converts SMILES notation to 3D molecular structures and renders them in an interactive viewer.

## Features

- Convert SMILES notation to 3D molecular structures
- Interactive 3D visualization of molecules
- Display molecular properties (formula, weight, atom count, etc.)
- Show 2D molecular structure images
- Built-in examples of common molecules

## Tech Stack

- **Tauri**: Framework for building desktop applications with web technologies
- **Rust**: Backend for Tauri, handling system-level operations
- **React + TypeScript**: Frontend UI framework
- **Python + Flask**: Service for molecular processing using RDKit
- **3DMol.js**: JavaScript library for molecular visualization

## Architecture

The application follows a modular architecture:

1. **Frontend (React/TypeScript)**
   - Handles UI rendering and user interactions
   - Communicates with the Python service via Tauri commands
   - Renders 3D molecules using 3DMol.js

2. **Backend (Tauri/Rust)**
   - Manages the application window and system integration
   - Provides commands to start and communicate with the Python service
   - Handles file system operations and IPC

3. **Python Service**
   - Processes SMILES strings using RDKit
   - Converts molecular representations to PDB format
   - Generates 2D molecular images
   - Calculates molecular properties

## Setup and Installation

### Prerequisites

- Node.js (v16+)
- Rust (latest stable)
- Python 3.8+ with pip
- RDKit (can be installed via conda or pip)

### Installation Steps

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/molecular-viewer.git
   cd molecular-viewer
   ```

2. Install JavaScript dependencies:
   ```
   npm install
   ```

3. Set up the Python service:
   ```
   cd python_service
   python -m venv venv
   # On Windows:
   venv\Scripts\activate
   # On macOS/Linux:
   source venv/bin/activate
   pip install -r requirements.txt
   ```

4. Run the development server:
   ```
   npm run tauri dev
   ```

## Building for Production

To build the application for production:

```
npm run tauri build
```

This will create platform-specific installers in the `src-tauri/target/release/bundle` directory.

## Project Structure

```
molecular-viewer/
├── src/                    # React frontend
│   ├── components/         # React components
│   ├── styles/             # CSS styles
│   ├── types/              # TypeScript type definitions
│   ├── App.tsx             # Main application component
│   └── main.tsx            # Entry point
├── src-tauri/              # Tauri/Rust backend
│   ├── src/                # Rust source code
│   ├── Cargo.toml          # Rust dependencies
│   └── tauri.conf.json     # Tauri configuration
├── python_service/         # Python service for molecular processing
│   ├── app.py              # Flask application
│   └── requirements.txt    # Python dependencies
└── package.json            # Node.js dependencies
```

## Extending the Application

This application is designed to be scalable and extensible. Here are some ways to extend it:

1. **Add more molecular formats**: Support additional molecular formats like MOL, SDF, etc.
2. **Implement molecular analysis**: Add features for analyzing molecular properties, simulations, etc.
3. **Enhance visualization**: Add more visualization options, coloring schemes, etc.
4. **Batch processing**: Add support for processing multiple molecules at once.

## License

[MIT License](LICENSE)

## Acknowledgements

- [RDKit](https://www.rdkit.org/) for molecular processing
- [3DMol.js](https://3dmol.csb.pitt.edu/) for molecular visualization
- [Tauri](https://tauri.app/) for the desktop application framework
- [React](https://reactjs.org/) for the UI framework
