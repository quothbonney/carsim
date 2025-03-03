import '@testing-library/jest-dom';
import { vi } from 'vitest';

// Mock 3DMol.js
(window as any).$3Dmol = {
  createViewer: vi.fn().mockReturnValue({
    addModel: vi.fn().mockReturnThis(),
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
    getModel: vi.fn().mockReturnValue({
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
      addAtoms: vi.fn().mockReturnThis(),
      removeAtoms: vi.fn().mockReturnThis(),
      addBonds: vi.fn().mockReturnThis(),
      removeBonds: vi.fn().mockReturnThis(),
      vibrate: vi.fn().mockReturnThis(),
      setDynamicStyle: vi.fn().mockReturnThis(),
    }),
    addSphere: vi.fn().mockReturnThis(),
    addArrow: vi.fn().mockReturnThis(),
    addCylinder: vi.fn().mockReturnThis(),
    addBox: vi.fn().mockReturnThis(),
    addIsosurface: vi.fn().mockReturnThis(),
    addSurface: vi.fn().mockReturnThis(),
    addLabel: vi.fn().mockReturnThis(),
    addLine: vi.fn().mockReturnThis(),
    addResLabels: vi.fn().mockReturnThis(),
    addPropertyLabels: vi.fn().mockReturnThis(),
    addUnitCell: vi.fn().mockReturnThis(),
    center: vi.fn().mockReturnThis(),
    resize: vi.fn().mockReturnThis(),
    spin: vi.fn().mockReturnThis(),
    stopAnimate: vi.fn().mockReturnThis(),
    animate: vi.fn().mockReturnThis(),
    setSlabNear: vi.fn().mockReturnThis(),
    setSlabFar: vi.fn().mockReturnThis(),
    setFog: vi.fn().mockReturnThis(),
    enableFog: vi.fn().mockReturnThis(),
    disableFog: vi.fn().mockReturnThis(),
    setWidth: vi.fn().mockReturnThis(),
    setHeight: vi.fn().mockReturnThis(),
    translate: vi.fn().mockReturnThis(),
    setProjection: vi.fn().mockReturnThis(),
    setClickable: vi.fn().mockReturnThis(),
    setColorByElement: vi.fn().mockReturnThis(),
    linkViewer: vi.fn().mockReturnThis(),
    setViewerOptions: vi.fn().mockReturnThis(),
    surfaceMeshes: {},
    surfaceObj: {},
  }),
  download: vi.fn(),
  createLabelSpec: vi.fn(),
  createSphere: vi.fn(),
  createBox: vi.fn(),
  createCylinder: vi.fn(),
  createArrow: vi.fn(),
  createShape: vi.fn(),
  GLModel: vi.fn(),
  elementColors: {
    H: '#FFFFFF',
    C: '#909090',
    N: '#3050F8',
    O: '#FF0D0D',
    F: '#90E050',
    P: '#FF8000',
    S: '#FFFF30',
    Cl: '#1FF01F',
    Br: '#A62929',
    I: '#940094',
  },
};

// Mock ResizeObserver
window.ResizeObserver = vi.fn().mockImplementation(() => ({
  observe: vi.fn(),
  unobserve: vi.fn(),
  disconnect: vi.fn(),
})) as any;

// Mock fetch
window.fetch = vi.fn() as any;

// Mock console methods for testing
console.log = vi.fn();
console.error = vi.fn();
console.warn = vi.fn(); 