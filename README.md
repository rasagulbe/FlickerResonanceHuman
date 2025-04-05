Code accompanying the publication: 
### [Attention differentially modulates the amplitude of resonance frequencies in the visual cortex](https://www.sciencedirect.com/science/article/abs/pii/S1053811919307372)
#### Rasa Gulbinaite, Diane H. M. Roozendaal, Rufin VanRullen, 2019, _NeuroImage_, 203:116146
##### DOI: [https://doi.org/10.1016/j.neuroimage.2019.116146](https://doi.org/10.1016/j.neuroimage.2019.116146)
##
![Github_FlickerHuman](https://github.com/user-attachments/assets/51c3a0e1-ca14-47d7-9e71-1b5827a054be)
##
### The repository contains MATLAB code for:
1. Constructing frequency-specific spatiotemporal filters that isolate steady-state evoked potentials (SSEPs) from endogenous brain activity;
2. Estimating statistical reliability of SSEP responses by constructing frequency-specific spatiotemporal filters using trials, on which neither flicker frequency nor its harmonics were present.
The code is an extension of previously published method called Rhythmic Entrainment Source Separation (RESS; Cohen and Gulbinaite, 2017 PMID: 27916666): https://github.com/mikexcohen/RESS 

##
### Dataset 
**Sample dataset** to be used with the code can be found here: https://osf.io/7s2vp/

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

