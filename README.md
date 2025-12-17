# T-REX Tools
T-REX scripts for detector coverage, energy resolution, frame-overlap time-distance diagrams.<br/>


TREX_FrameOverlapChecks.py:<br/>
- Script looks at a central wavelength, M-chopper speed and plots your time-distance diagram between sample and detectors.<br/><br/>
- All analytically calculated wavelengths and ToF's in the .py scripts have been cross checked with McStas.<br/><br/>
- Users can set an energy loss/gain region they are interested in. For this script, energy loss is set to 80% and energy gain is set to 'infinite' where elastically scattered neutrons reach the detectors instantaneously. This is for both cold and thermal wavelengths.<br/><br/>
- Plot at the end looks at how much overlap there is and how much of the gain and loss windows you can use before overlap occurs. Needs a bit of tweaking/wrap around to match the first rep of ESS pulse (n+1) with the last rep of ESS pulse (n).<br/><br/>

TREX_DetectorCoverage.py:<br/>
- Script shows the Q-E region covered by 4 boxes vs 10 boxes of Multi-Grid detectors.<br/><br/> - Powder detector coverage only for now.<br/><br/> - Single crystal coverage is not included yet.<br/><br/> - Note that the detector coverage starts from 5 degrees. This is just to account for direct beam 'splurge' covering the lowest angles.<br/><br/>


InelasticResolution_OptimisedFluxResolution.py: <br/>
- Calculates the inelastic resolution for a given incident wavelength for T-REX.<br/><br/>
- Uses the Lechner formula for the resolution. Optimises at the elastic, can be adapted for any given energy transfer. the Delta E result IS THE FWHM OF T-REX according to the Lechner paper. <br/><br/>
- The opening times are optimised at the elastic line following the formalism of the attached paper.<br/><br/>
- One can optimise in the inelastic regime. This is more relevant for the thermal wavelengths than cold ones.<br/><br/>
- McStas Monitors show the P-chopper peak shape is a top-hat function: Scaling the opening time by 0.2881 to recreate that shape.<br/><br/>
- McStas Monitors show the M-chopper peak shape is a triangular function: Scaling the opening time by 0.2041 to recreate that shape.<br/><br/>


Things that will be added in future scripts:<br/>
- Calculation of the Q resolution of T-REX for a given wavelength.<br/><br/>
