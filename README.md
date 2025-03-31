# FlickerResonanceHuman

## The code accompanying the paper "Attention differentially modulates the amplitude of resonance frequencies in the visual cortex"

### Rasa Gulbinaite, Diane H. M. Roozendaal, Rufin VanRullen, 2019, _NeuroImage_, 203:116146
##
The repository contains:

1. **MATLAB code** used to: 
    * construct frequency-specific spatiotemporal filters that isolate steady-state evoked potentials (SSEPs) from endogenous brain activity; 
    * estimate statistical reliability of SSEP responses by constructing frequency-specific spatiotemporal filters using trials, on which neither flicker frequency nor its harmonics were present.

    The code is an extension of previously published method called Rhythmic Entrainment Source Separation (RESS; Cohen and Gulbinaite, 2017 PMID: 27916666): https://github.com/mikexcohen/RESS 

2. **Sample dataset** to be used with the code can be found here: https://osf.io/7s2vp/

   Dataset is in EEGLAB format. Epochs: -1 to 9 sec, where 0 us the cue onset.
   Trial markers are organized as follows (cue location / LED state / cue validity):
   * Attend Left/ Flicker/ Valid - 193
   * Attend Left/ Flicker/ Invalid - 194
   * Attend Left/ Static/ Valid - 195
   * Attend Left/ Static/ Invalid - 196
   * Attend Right/ Flicker/ Valid - 197
   * Attend Right/ Flicker/ Invalid - 198
   * Attend Right/ Static/ Valid - 199
   * Attend Right/ Static/ Invalid - 200

