import { describe, it, expect, vi, beforeEach } from 'vitest';
import { render, screen, fireEvent, waitFor } from '@testing-library/react';
import MoleculeControls from './MoleculeControls';

describe('MoleculeControls Component', () => {
  // Mock props
  const mockProps = {
    smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
    setSmiles: vi.fn(),
    onProcess: vi.fn(),
    isLoading: false,
    onPdbUpload: vi.fn(),
  };

  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('renders correctly with default props', () => {
    render(<MoleculeControls {...mockProps} />);
    
    // Check for SMILES input
    expect(screen.getByPlaceholderText(/Enter SMILES notation/i)).toBeInTheDocument();
    
    // Check for Process button
    expect(screen.getByText(/Process Molecule/i)).toBeInTheDocument();
    
    // Check for PDB section
    expect(screen.getByText(/PDB Structure/i)).toBeInTheDocument();
  });

  it('handles SMILES input changes', () => {
    render(<MoleculeControls {...mockProps} />);
    
    // Get SMILES input
    const smilesInput = screen.getByPlaceholderText(/Enter SMILES notation/i);
    
    // Change input value
    fireEvent.change(smilesInput, { target: { value: 'C1=CC=CC=C1' } });
    
    // Check if setSmiles was called with new value
    expect(mockProps.setSmiles).toHaveBeenCalledWith('C1=CC=CC=C1');
  });

  it('handles form submission', () => {
    const { container } = render(<MoleculeControls {...mockProps} />);
    
    // Get form
    const form = container.querySelector('form');
    expect(form).not.toBeNull();
    
    // Submit form
    if (form) fireEvent.submit(form);
    
    // Check if onProcess was called with current SMILES
    expect(mockProps.onProcess).toHaveBeenCalledWith(mockProps.smiles);
  });

  it('disables Process button when loading', () => {
    render(<MoleculeControls {...mockProps} isLoading={true} />);
    
    // Get Process button
    const processButton = screen.getByText(/Processing/i);
    
    // Check if button is disabled
    expect(processButton).toBeDisabled();
  });

  it('toggles SMILES section visibility', () => {
    render(<MoleculeControls {...mockProps} />);
    
    // Get toggle button
    const toggleButton = screen.getByText(/SMILES Input/i).parentElement;
    
    // Initially expanded
    expect(screen.getByPlaceholderText(/Enter SMILES notation/i)).toBeVisible();
    
    // Click toggle button
    if (toggleButton) fireEvent.click(toggleButton);
    
    // Check if section is no longer in the document
    expect(screen.queryByPlaceholderText(/Enter SMILES notation/i)).not.toBeInTheDocument();
    
    // Click toggle button again
    if (toggleButton) fireEvent.click(toggleButton);
    
    // Check if section is expanded again
    expect(screen.getByPlaceholderText(/Enter SMILES notation/i)).toBeInTheDocument();
  });

  it('toggles PDB section visibility', () => {
    render(<MoleculeControls {...mockProps} />);
    
    // Get toggle button
    const toggleButton = screen.getByText(/PDB Structure/i).parentElement;
    
    // Initially expanded
    expect(screen.getByText(/Upload PDB/i)).toBeVisible();
    
    // Click toggle button
    if (toggleButton) fireEvent.click(toggleButton);
    
    // Check if section is no longer in the document
    expect(screen.queryByText(/Upload PDB/i)).not.toBeInTheDocument();
    
    // Click toggle button again
    if (toggleButton) fireEvent.click(toggleButton);
    
    // Check if section is expanded again
    expect(screen.getByText(/Upload PDB/i)).toBeInTheDocument();
  });

  it('handles PDB file upload', async () => {
    const { container } = render(<MoleculeControls {...mockProps} />);
    
    // Mock file
    const file = new File(['HEADER    PROTEIN'], 'test.pdb', { type: 'chemical/x-pdb' });
    
    // Get file input directly since it doesn't have a label
    const fileInput = container.querySelector('input[type="file"]');
    expect(fileInput).not.toBeNull();
    
    // Upload file
    if (fileInput) {
      fireEvent.change(fileInput, { target: { files: [file] } });
    }
    
    // Check if onPdbUpload was called
    await waitFor(() => {
      expect(mockProps.onPdbUpload).toHaveBeenCalledWith('HEADER    PROTEIN');
    });
  });

  it('displays error message when PDB upload fails', async () => {
    // Mock FileReader to simulate error
    const originalFileReader = window.FileReader;
    window.FileReader = class MockFileReader {
      onerror: any;
      readAsText() {
        setTimeout(() => this.onerror(new Error('File read error')), 0);
      }
    } as any;
    
    const { container } = render(<MoleculeControls {...mockProps} />);
    
    // Mock file
    const file = new File(['INVALID'], 'test.pdb', { type: 'chemical/x-pdb' });
    
    // Get file input directly
    const fileInput = container.querySelector('input[type="file"]');
    expect(fileInput).not.toBeNull();
    
    // Upload file
    if (fileInput) {
      fireEvent.change(fileInput, { target: { files: [file] } });
    }
    
    // Check for error message
    await waitFor(() => {
      expect(screen.getByText(/Failed to read/i)).toBeInTheDocument();
    });
    
    // Restore original FileReader
    window.FileReader = originalFileReader;
  });

  it('handles PDB ID fetch', async () => {
    // Mock fetch response
    window.fetch = vi.fn().mockResolvedValueOnce({
      ok: true,
      json: () => Promise.resolve({ pdbContent: 'HEADER    PROTEIN' }),
    }) as any;
    
    render(<MoleculeControls {...mockProps} />);
    
    // Get PDB ID input
    const pdbIdInput = screen.getByPlaceholderText(/PDB ID/i);
    
    // Enter PDB ID
    fireEvent.change(pdbIdInput, { target: { value: '1ABC' } });
    
    // Get fetch button
    const fetchButton = screen.getByText(/Fetch/i);
    
    // Click fetch button
    fireEvent.click(fetchButton);
    
    // Check if fetch was called
    expect(window.fetch).toHaveBeenCalledWith(expect.stringContaining('1ABC'));
    
    // Check if onPdbUpload was called
    await waitFor(() => {
      expect(mockProps.onPdbUpload).toHaveBeenCalledWith('HEADER    PROTEIN');
    });
  });

  it('handles PDB ID fetch error', async () => {
    // Mock fetch with error response
    window.fetch = vi.fn().mockRejectedValueOnce(new Error('Network error')) as any;
    
    render(<MoleculeControls {...mockProps} />);
    
    // Get PDB ID input
    const pdbIdInput = screen.getByPlaceholderText(/PDB ID/i);
    
    // Enter PDB ID
    fireEvent.change(pdbIdInput, { target: { value: '1ABC' } });
    
    // Get fetch button
    const fetchButton = screen.getByText(/Fetch/i);
    
    // Click fetch button
    fireEvent.click(fetchButton);
    
    // Check for error message
    await waitFor(() => {
      expect(screen.getByText(/Failed to fetch PDB/i)).toBeInTheDocument();
    });
  });

  it('disables fetch button when loading', () => {
    render(<MoleculeControls {...mockProps} isLoading={true} />);
    
    // Get fetch button
    const fetchButton = screen.getByText(/Fetch/i);
    
    // Check if button is disabled
    expect(fetchButton).toBeDisabled();
  });

  it('disables upload button when loading', () => {
    render(<MoleculeControls {...mockProps} isLoading={true} />);
    
    // Get upload button
    const uploadButton = screen.getByText(/Upload PDB/i);
    
    // Check if button is disabled
    expect(uploadButton).toBeDisabled();
  });

  it('handles example molecule buttons', () => {
    render(<MoleculeControls {...mockProps} />);
    
    // Get example buttons
    const aspirinButton = screen.getByText(/Aspirin/i);
    
    // Click aspirin button
    fireEvent.click(aspirinButton);
    
    // Check if setSmiles was called with aspirin SMILES
    expect(mockProps.setSmiles).toHaveBeenCalled();
    
    // Now submit the form to trigger onProcess
    const processButton = screen.getByText(/Process Molecule/i);
    const form = processButton.closest('form');
    if (form) fireEvent.submit(form);
    
    // Check if onProcess was called
    expect(mockProps.onProcess).toHaveBeenCalled();
  });
}); 