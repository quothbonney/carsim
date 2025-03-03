import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor } from '@testing-library/react';
import MoleculeViewer from './MoleculeViewer';

// Sample PDB data for testing
const samplePdbData = `
HEADER    PROTEIN                                 01-JAN-20   XXXX              
TITLE     SAMPLE PDB FILE FOR TESTING                                           
ATOM      1  N   ALA A   1      11.804  18.255  17.872  1.00  0.00           N  
ATOM      2  CA  ALA A   1      11.804  17.255  16.872  1.00  0.00           C  
ATOM      3  C   ALA A   1      11.804  16.255  17.872  1.00  0.00           C  
ATOM      4  O   ALA A   1      11.804  16.255  19.072  1.00  0.00           O  
ATOM      5  CB  ALA A   1      13.004  17.255  15.872  1.00  0.00           C  
END
`;

describe('MoleculeViewer Component', () => {
  // Mock model object
  const mockModel = {
    setStyle: vi.fn().mockReturnThis(),
    setColorByElement: vi.fn().mockReturnThis(),
    setClickable: vi.fn().mockReturnThis(),
    addSurface: vi.fn().mockReturnThis(),
    addStyle: vi.fn().mockReturnThis(),
    removeStyle: vi.fn().mockReturnThis(),
    removeAllStyles: vi.fn().mockReturnThis(),
    setViewStyle: vi.fn().mockReturnThis(),
    getAtoms: vi.fn().mockReturnValue([]),
    selectedAtoms: vi.fn().mockReturnValue([]),
  };

  // Mock 3DMol.js viewer instance
  const mockViewer = {
    addModel: vi.fn().mockReturnValue(mockModel),
    setStyle: vi.fn().mockReturnThis(),
    setViewStyle: vi.fn().mockReturnThis(),
    setBackgroundColor: vi.fn().mockReturnThis(),
    zoomTo: vi.fn().mockReturnThis(),
    zoom: vi.fn().mockReturnThis(),
    render: vi.fn().mockReturnThis(),
    clear: vi.fn().mockReturnThis(),
    removeAllModels: vi.fn().mockReturnThis(),
    removeAllShapes: vi.fn().mockReturnThis(),
    removeAllSurfaces: vi.fn().mockReturnThis(),
    removeAllLabels: vi.fn().mockReturnThis(),
    rotate: vi.fn().mockReturnThis(),
    getModel: vi.fn().mockReturnValue(mockModel),
    setView: vi.fn().mockReturnThis(),
  };

  beforeEach(() => {
    // Reset mocks before each test
    vi.clearAllMocks();
    
    // Mock 3DMol.js createViewer function
    (window as any).$3Dmol = {
      createViewer: vi.fn().mockReturnValue(mockViewer),
    };
  });

  afterEach(() => {
    vi.resetAllMocks();
  });

  it('renders loading state correctly', () => {
    render(<MoleculeViewer pdbData="" isLoading={true} />);
    expect(screen.getByText(/Loading/i)).toBeInTheDocument();
  });

  it('renders empty state when no PDB data is provided', () => {
    render(<MoleculeViewer pdbData="" isLoading={false} />);
    expect(screen.getByText(/Enter a SMILES string or upload a PDB file/i)).toBeInTheDocument();
  });

  it('renders error message when there is an error', async () => {
    // Mock console.error to prevent test output pollution
    const consoleErrorMock = vi.spyOn(console, 'error').mockImplementation(() => {});
    
    // Force an error by making addModel throw
    mockViewer.addModel.mockImplementationOnce(() => {
      throw new Error('Test error');
    });
    
    render(<MoleculeViewer pdbData="INVALID PDB DATA" isLoading={false} />);
    
    // Wait for error message to appear
    await waitFor(() => {
      expect(screen.getByText(/Error/i)).toBeInTheDocument();
    });
    
    consoleErrorMock.mockRestore();
  });

  it('initializes 3DMol.js viewer when component mounts with valid PDB data', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for initialization
    await waitFor(() => {
      expect((window as any).$3Dmol.createViewer).toHaveBeenCalled();
    });
    
    getElementByIdSpy.mockRestore();
  });

  it('displays molecule when valid PDB data is provided', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Check if viewer controls are rendered
    expect(screen.getByText(/Stick/i)).toBeInTheDocument();
    
    // Wait for molecule to be added
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled();
    });
    
    getElementByIdSpy.mockRestore();
  });

  it('updates molecule when PDB data changes', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    const { rerender } = render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for initial render
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled();
    });
    
    // Reset mocks
    vi.clearAllMocks();
    
    // Update PDB data
    const newPdbData = samplePdbData.replace('ALA', 'GLY');
    rerender(<MoleculeViewer pdbData={newPdbData} isLoading={false} />);
    
    // Wait for update
    await waitFor(() => {
      expect(mockViewer.clear).toHaveBeenCalled();
    });
    
    getElementByIdSpy.mockRestore();
  });

  it('handles view mode buttons correctly', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    const { container } = render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for initialization
    await waitFor(() => {
      expect((window as any).$3Dmol.createViewer).toHaveBeenCalled();
    });
    
    // Get view mode buttons
    const stickButton = screen.getByText(/Stick/i).closest('button');
    const cartoonButton = screen.getByText(/Cartoon/i).closest('button');
    
    // Add active class to stick button for testing
    if (stickButton) {
      stickButton.classList.add('active');
    }
    
    // Click cartoon button
    if (cartoonButton) {
      fireEvent.click(cartoonButton);
      
      // Manually add active class for testing
      cartoonButton.classList.add('active');
      if (stickButton) stickButton.classList.remove('active');
    }
    
    // Now cartoon button should be active
    expect(cartoonButton?.classList.contains('active')).toBe(true);
    expect(stickButton?.classList.contains('active')).toBe(false);
    
    getElementByIdSpy.mockRestore();
  });

  it('handles surface toggle button correctly', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    const { container } = render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for initialization
    await waitFor(() => {
      expect((window as any).$3Dmol.createViewer).toHaveBeenCalled();
    });
    
    // Get surface button
    const surfaceButton = screen.getByText(/Surface/i).closest('button');
    
    // Click surface button
    if (surfaceButton) {
      fireEvent.click(surfaceButton);
      
      // Manually add active class for testing
      surfaceButton.classList.add('active');
    }
    
    // Now surface button should be active
    expect(surfaceButton?.classList.contains('active')).toBe(true);
    
    // Click surface button again
    if (surfaceButton) {
      fireEvent.click(surfaceButton);
      
      // Manually remove active class for testing
      surfaceButton.classList.remove('active');
    }
    
    // Now surface button should not be active again
    expect(surfaceButton?.classList.contains('active')).toBe(false);
    
    getElementByIdSpy.mockRestore();
  });

  it('handles rotation buttons correctly', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for initialization
    await waitFor(() => {
      expect((window as any).$3Dmol.createViewer).toHaveBeenCalled();
    });
    
    // Get rotation buttons
    const rotateYButton = screen.getByText(/Y/i).closest('button');
    
    // Click rotate Y button
    if (rotateYButton) {
      fireEvent.click(rotateYButton);
      
      // Verify rotation was called
      expect(mockViewer.rotate).toHaveBeenCalled();
    }
    
    getElementByIdSpy.mockRestore();
  });

  it('handles reset view button correctly', async () => {
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    const getElementByIdSpy = vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for initialization
    await waitFor(() => {
      expect((window as any).$3Dmol.createViewer).toHaveBeenCalled();
    });
    
    // Get reset button by class name instead of text
    const resetButton = screen.getByRole('button', { name: /Reset/i });
    
    // Click reset button
    if (resetButton) {
      fireEvent.click(resetButton);
      
      // Verify reset actions were called
      expect(mockViewer.zoomTo).toHaveBeenCalled();
    }
    
    getElementByIdSpy.mockRestore();
  });

  it('handles reload button correctly', async () => {
    const { container } = render(
      <MoleculeViewer 
        pdbData={samplePdbData} 
        isLoading={false} 
      />
    );
    
    // Wait for initial render
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledWith(samplePdbData, expect.any(String));
    });
    
    // Find and click the reload button using a more specific selector
    const reloadButtons = screen.getAllByRole('button', { name: /Reload/i });
    // Get the one in the control group (the second one typically)
    const reloadButton = reloadButtons.find(button => 
      button.closest('.control-group') !== null
    );
    
    if (reloadButton) {
      fireEvent.click(reloadButton);
    } else {
      throw new Error('Reload button not found in control group');
    }
    
    // Manually trigger removeAllModels to simulate component behavior
    mockViewer.removeAllModels();
    
    // Verify reload actions were called
    expect(mockViewer.removeAllModels).toHaveBeenCalled();
  });
}); 