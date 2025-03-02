declare namespace $3Dmol {
  interface Viewer {
    addModel(data: string, format: string): void;
    setStyle(sel: any, style: any): void;
    setViewStyle(style: any): void;
    zoomTo(): void;
    render(): void;
    clear(): void;
    rotate(angle: number, axis: string): void;
    addSurface(type: number, options: any): void;
  }

  const SurfaceType: {
    VDW: number;
    MS: number;
    SAS: number;
    SES: number;
  };

  function createViewer(element: HTMLElement, options: any): Viewer;
}

interface Window {
  $3Dmol: typeof $3Dmol;
} 