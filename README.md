# T-REX Tools
T-REX scripts for detector coverage, energy resolution, frame-overlap diagrams.<br/>


TREX_FrameOverlapChecks.py:<br/>
- Script looks at a central wavelength, M-chopper speed and plots your time-distance diagram between sample and detectors.<br/><br/>
- Users can set an energy loss/gain region they are interested in. Rule of thumb is 85% loss for thermal and 20% loss for cold wavelengths.<br/><br/>

TREX_DetectorCoverage.py:<br/>
- Script shows the Q-E region covered by 4 boxes vs 10 boxes of Multi-Grid detectors.<br/><br/> - Powder detector coverage only for now.<br/><br/> - Single crystal coverage is not included yet.<br/><br/>


InelasticResolution_OptimisedFluxResolution.py: <br/>
- Calculates the inelastic resolution for a given incident wavelength for T-REX.<br/><br/>
- Uses the Lechner formual for the resolution.<br/><br/>
- The opening times are optimised following the formalism of the attached paper.<br/><br/>
- McStas Monitors show the P-chopper peak shape is a top-hat function: Scaling the opening time by 0.2881 to recreate that shape.<br/><br/>
- McStas Monitors show the M-chopper peak shape is a triangular function: Scaling the opening time by 0.2041 to recreate that shape.<br/><br/>


