# ephys_tutorials

### what you will learn:
Ephys tutorial 1 show on an entry-level how ECOG and LFP data are organized. Topics covered are:
* how explore data you recieve
* rearrange data in readable format
* plot channels
* find the sample frequency
* interpret some information (Stimulation state, medication state, stimulator model) and where to store these
* rewrite channel names
* how to do rereferencing of the signal
* save your data in Matlab, FieldTrip or SPM readable data
* convert your data to BIDS format (brainvision files)
* explore the data in Python and use mne

### what you need:
####MATLAB
For the first part, you will need Matlab2020b as we start with this script https://github.com/neuromodulation/ephys_tutorials/blob/main/Tutorial_I_Import_and_rereference.m
the dependencies are:
spm: https://github.com/spm/spm12
wjn_toolbox: https://github.com/neuromodulation/wjn_toolbox
FieldTrip: https://github.com/fieldtrip/fieldtrip.git
MatLab Signal Processing Toolbox

####Python
For the second part, you need jupyter notebook as we end with this script https://github.com/neuromodulation/ephys_tutorials/blob/main/Tutorial_1_import_and_rereferences.ipynb
mne tools https://mne.tools/stable/index.html, that can be installed using following command:
```python
import sys
!{sys.executable} -m pip  install mne
```

This tutorial is created by the Interventional and Cognitive Neuromodulation Lab: https://www.icneuromodulation.org/
