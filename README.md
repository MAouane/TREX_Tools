# T-REX Tools
T-REX scripts for detector coverage, energy resolution, frame-overlap diagrams


TREX_FrameOverlapChecks.py:
Script looks at a central wavelength, M-chopper speed and plots your time-distance diagram between sample and detectors.
Users can set an energy loss/gain region they are interested in. Rule of thumb is 85% loss for thermal and 20% loss for cold wavelengths

TREX_DetectorCoverage.py:
Script shows the Q-E region covered by 4 boxes vs 10 boxes of Multi-Grid detectors
For now it's only for powder samples. Will upgrade it for single crystal coverage eventually (a bit more work)


InelasticResolution_OptimisedFluxResolution.py:
Calculates the inelastic resolution for a given incident wavelength for T-REX.
Uses the Lechner formula for the resolution function
The flux/resolution ratio is optimised by matching the P-chopper and M-chopper opening times using the PWR Optimisation in this paper:
The P-chopper peak shape is a top-hat when looking at McStas: A factor of 0.2887 is applied to the opening time tau_P to recreate the top-hat shape,
The M-chopper peak shape is triangular when looking at McStas: A factor of 0.2041 is applied to the opening time tau_M to recreate the triangular shape.
