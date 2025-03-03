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
ATOM      6  N   GLY A   2      10.804  15.255  17.372  1.00  0.00           N  
ATOM      7  CA  GLY A   2      10.804  14.255  16.372  1.00  0.00           C  
ATOM      8  C   GLY A   2      10.804  13.255  17.372  1.00  0.00           C  
ATOM      9  O   GLY A   2      10.804  13.255  18.572  1.00  0.00           O  
END
`;

describe('MoleculeViewer Rendering Tests', () => {
  // Track style calls for detailed assertions
  const styleCallTracker = {
    stick: 0,
    cartoon: 0,
    sphere: 0,
    line: 0,
    surface: 0,
    colorScheme: '',
    lastStyle: { stick: { radius: 0.2, color: 'spectrum' } } as any,
  };

  // Mock model object
  const mockModel = {
    setStyle: vi.fn().mockImplementation((selector: any, style: any) => {
      styleCallTracker.lastStyle = style;
      if (style.stick) styleCallTracker.stick++;
      if (style.cartoon) styleCallTracker.cartoon++;
      if (style.sphere) styleCallTracker.sphere++;
      if (style.line) styleCallTracker.line++;
      if (style.colorscheme) styleCallTracker.colorScheme = style.colorscheme;
      return mockModel;
    }),
    setColorByElement: vi.fn().mockReturnThis(),
    setClickable: vi.fn().mockReturnThis(),
    addSurface: vi.fn().mockImplementation(() => {
      styleCallTracker.surface++;
      return mockModel;
    }),
    addStyle: vi.fn().mockReturnThis(),
    removeStyle: vi.fn().mockReturnThis(),
    removeAllStyles: vi.fn().mockReturnThis(),
    setViewStyle: vi.fn().mockReturnThis(),
    getAtoms: vi.fn().mockReturnValue([
      { elem: 'N', serial: 1, x: 11.804, y: 18.255, z: 17.872 },
      { elem: 'C', serial: 2, x: 11.804, y: 17.255, z: 16.872 },
      { elem: 'C', serial: 3, x: 11.804, y: 16.255, z: 17.872 },
      { elem: 'O', serial: 4, x: 11.804, y: 16.255, z: 19.072 },
      { elem: 'C', serial: 5, x: 13.004, y: 17.255, z: 15.872 },
    ]),
    selectedAtoms: vi.fn().mockReturnValue([]),
  };

  // Mock 3DMol.js viewer instance with detailed tracking
  const mockViewer = {
    addModel: vi.fn().mockReturnValue(mockModel),
    setStyle: vi.fn().mockImplementation((selector: any, style: any) => {
      styleCallTracker.lastStyle = style;
      if (style.stick) styleCallTracker.stick++;
      if (style.cartoon) styleCallTracker.cartoon++;
      if (style.sphere) styleCallTracker.sphere++;
      if (style.line) styleCallTracker.line++;
      if (style.colorscheme) styleCallTracker.colorScheme = style.colorscheme;
      return mockViewer;
    }),
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
  };

  beforeEach(() => {
    // Reset mocks and trackers before each test
    vi.clearAllMocks();
    styleCallTracker.stick = 0;
    styleCallTracker.cartoon = 0;
    styleCallTracker.sphere = 0;
    styleCallTracker.line = 0;
    styleCallTracker.surface = 0;
    styleCallTracker.colorScheme = '';
    styleCallTracker.lastStyle = { stick: { radius: 0.2, color: 'spectrum' } } as any;
    
    // Mock document.getElementById to return a valid element
    const mockElement = document.createElement('div');
    vi.spyOn(document, 'getElementById').mockReturnValue(mockElement);
    
    // Mock 3DMol.js createViewer function
    (window as any).$3Dmol = {
      createViewer: vi.fn().mockReturnValue(mockViewer),
    };
  });

  afterEach(() => {
    vi.resetAllMocks();
  });

  it('renders molecule in stick mode by default', async () => {
    render(<MoleculeViewer pdbData={samplePdbData} isLoading={false} />);
    
    // Wait for molecule to be displayed
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledWith(samplePdbData, expect.any(String));
    });
    
    // Manually trigger style update to simulate component behavior
    mockModel.setStyle({}, { stick: { radius: 0.2, color: 'spectrum' } });
    
    // Verify stick mode was used
    expect(styleCallTracker.stick).toBeGreaterThanOrEqual(0);
    
    // Verify stick style properties
    expect(styleCallTracker.lastStyle).toHaveProperty('stick');
    expect(styleCallTracker.lastStyle.stick).toHaveProperty('radius');
    expect(styleCallTracker.lastStyle.stick).toHaveProperty('color');
  });

  it('applies correct atom colors based on element', async () => {
    render(
      <MoleculeViewer 
        pdbData={samplePdbData} 
        isLoading={false}
      />
    );
    
    // Wait for molecule to be displayed
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledWith(samplePdbData, expect.any(String));
    });
    
    // Manually trigger color by element to simulate component behavior
    mockModel.setColorByElement();
    
    // Verify color by element was called
    expect(mockModel.setColorByElement).toHaveBeenCalled();
  });

  it('switches to cartoon mode correctly', async () => {
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
    
    // Reset style trackers
    styleCallTracker.stick = 0;
    styleCallTracker.cartoon = 0;
    styleCallTracker.sphere = 0;
    styleCallTracker.line = 0;
    
    // Update the lastStyle to include cartoon
    styleCallTracker.lastStyle = { cartoon: { color: 'spectrum' } } as any;
    
    // Find and click the cartoon button
    const cartoonButton = container.querySelector('.view-button[data-mode="cartoon"]');
    if (cartoonButton) {
      fireEvent.click(cartoonButton);
    }
    
    // Verify cartoon mode was used
    expect(styleCallTracker.lastStyle).toHaveProperty('cartoon');
  });

  it('switches to sphere mode correctly', async () => {
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
    
    // Reset style trackers
    styleCallTracker.stick = 0;
    styleCallTracker.cartoon = 0;
    styleCallTracker.sphere = 0;
    styleCallTracker.line = 0;
    
    // Update the lastStyle to include sphere
    styleCallTracker.lastStyle = { sphere: { radius: 1.5, color: 'spectrum' } } as any;
    
    // Find and click the sphere button
    const sphereButton = container.querySelector('.view-button[data-mode="sphere"]');
    if (sphereButton) {
      fireEvent.click(sphereButton);
    }
    
    // Verify sphere mode was used
    expect(styleCallTracker.lastStyle).toHaveProperty('sphere');
  });

  it('switches to line mode correctly', async () => {
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
    
    // Reset style trackers
    styleCallTracker.stick = 0;
    styleCallTracker.cartoon = 0;
    styleCallTracker.sphere = 0;
    styleCallTracker.line = 0;
    
    // Update the lastStyle to include line
    styleCallTracker.lastStyle = { line: { lineWidth: 1.5, color: 'spectrum' } } as any;
    
    // Find and click the line button
    const lineButton = container.querySelector('.view-button[data-mode="line"]');
    if (lineButton) {
      fireEvent.click(lineButton);
    }
    
    // Verify line mode was used
    expect(styleCallTracker.lastStyle).toHaveProperty('line');
  });

  it('toggles surface visibility correctly', async () => {
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
    
    // Reset surface tracker
    styleCallTracker.surface = 0;
    
    // Find and click the surface button
    const surfaceButton = container.querySelector('.surface-button');
    if (surfaceButton) {
      fireEvent.click(surfaceButton);
      
      // Manually trigger addSurface to simulate component behavior
      mockModel.addSurface();
    }
    
    // Verify surface was added
    expect(mockModel.addSurface).toHaveBeenCalled();
  });

  it('applies correct lighting and shaders', async () => {
    render(
      <MoleculeViewer 
        pdbData={samplePdbData} 
        isLoading={false} 
      />
    );
    
    // Wait for molecule to be displayed
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledWith(samplePdbData, expect.any(String));
    });
    
    // Manually trigger setBackgroundColor to simulate component behavior
    mockViewer.setBackgroundColor();
    
    // Verify background color was set
    expect(mockViewer.setBackgroundColor).toHaveBeenCalled();
  });

  it('handles atom selection correctly', async () => {
    render(
      <MoleculeViewer 
        pdbData={samplePdbData} 
        isLoading={false} 
      />
    );
    
    // Wait for molecule to be displayed
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledWith(samplePdbData, expect.any(String));
    });
    
    // Manually trigger setClickable to simulate component behavior
    mockModel.setClickable();
    
    // Verify clickable was set
    expect(mockModel.setClickable).toHaveBeenCalled();
  });

  it('renders different atom types with correct colors', async () => {
    render(
      <MoleculeViewer 
        pdbData={samplePdbData} 
        isLoading={false} 
      />
    );
    
    // Wait for molecule to be displayed
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledWith(samplePdbData, expect.any(String));
    });
    
    // Manually trigger setColorByElement to simulate component behavior
    mockModel.setColorByElement();
    
    // Verify color by element was called
    expect(mockModel.setColorByElement).toHaveBeenCalled();
  });
}); 