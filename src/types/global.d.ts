// Define the $3Dmol global variable
declare namespace $3Dmol {
  interface ViewerConfig {
    backgroundColor?: string;
    antialias?: boolean;
    id?: string;
    width?: number;
    height?: number;
  }

  interface StyleSpec {
    stick?: {
      radius?: number;
      colorscheme?: string;
    };
    sphere?: {
      scale?: number;
    };
  }

  interface ViewStyleSpec {
    style?: string;
  }

  interface SurfaceSpec {
    opacity?: number;
    color?: string;
  }

  const SurfaceType: {
    VDW: number;
    MS: number;
    SAS: number;
    SES: number;
  };

  interface Viewer {
    clear(): void;
    addModel(data: string, format: string): void;
    setStyle(sel: object, style: StyleSpec): void;
    setViewStyle(style: ViewStyleSpec): void;
    addSurface(type: number, spec: SurfaceSpec): void;
    zoomTo(): void;
    render(): void;
    rotate(angle: number, axis: string): void;
  }

  function createViewer(element: HTMLElement, config: ViewerConfig): Viewer;
}

interface Window {
  $3Dmol: typeof $3Dmol;
} 