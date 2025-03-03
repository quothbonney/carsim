import subprocess
import sys
import os

def install_biopython():
    """
    Install BioPython in the current Python environment.
    """
    print("Attempting to install BioPython...")
    
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython==1.81"])
        print("BioPython successfully installed!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error installing BioPython: {e}")
        return False

if __name__ == "__main__":
    success = install_biopython()
    
    if success:
        print("You can now use BioPython in your Python service.")
    else:
        print("\nAlternative installation methods:")
        print("1. Try installing with conda: conda install -c conda-forge biopython")
        print("2. Install manually in the python_service/venv environment:")
        
        # Check if we're on Windows or Unix
        if os.name == 'nt':  # Windows
            print("   python_service\\venv\\Scripts\\pip install biopython==1.81")
        else:  # Unix/Linux/Mac
            print("   python_service/venv/bin/pip install biopython==1.81") 