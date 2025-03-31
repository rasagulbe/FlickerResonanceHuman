# FlickerResonanceHuman



The repository contains:

    MATLAB code used to: (A) construct frequency-specific spatiotemporal filters that isolate steady-state evoked potentials (SSEPs) from endogenous brain activity; (B) estimate statistical reliability of SSEP responses by constructing frequency-specific spatiotemporal filters using trials, on which neither flicker frequency nor its harmonics were present.

    The code is an extension of previously published method called Rhythmic Entrainment Source Separation (RESS; Cohen and Gulbinaite, 2017 PMID: 27916666).

    Sample dataset, which was pre-processed as described in https://www.biorxiv.org/content/10.1101/518779v1. Dataset is in EEGLAB format (https://sccn.ucsd.edu/eeglab/index.php). Epochs: -1 to 9 sec, where 0 us the cue onset. Trial markers are organized as follows (cue location / LED state / cue validity):
        Attend Left/ Flicker/ Valid - 193
        Attend Left/ Flicker/ Invalid - 194
        Attend Left/ Static/ Valid - 195
        Attend Left/ Static/ Invalid - 196
        Attend Right/ Flicker/ Valid - 197
        Attend Right/ Flicker/ Invalid - 198
        Attend Right/ Static/ Valid - 199
        Attend Right/ Static/ Invalid - 200

