%% First pull/clone the latest version of the perceive toolbox e.g. to 
mkdir('C:\code\')
cd('C:\code\')
!git clone https://github.com/neuromodulation/perceive/

%% This is a tutorial for perceive analyses
addpath C:\code\perceive
cd('C:\code\ephys_tutorials\')
%% Ready to start
clear all, close all, clc
% First run the full perceive toolbox to test that everything works
perceive('perceive_report.json','tutorial_subject')

%% A new folder was created including a lot of derived information that were extracted from the json
% Let's have a look at the ECG detection:

open(fullfile('tutorial_subject','ses-2020011509015499','ieeg','tutorial_subject_ses-2020011509015499_run-BSTD20200115085400_ECG_LFP_Gpi_R_02.fig'))

% this looks like the ECG detection did a good job so let's take this as a
% starter and rerun it out of the perceive toolbox:

load(fullfile('tutorial_subject','ses-2020011509015499','ieeg','tutorial_subject_ses-2020011509015499_run-BSTD20200115085400.mat'))

% This has loaded the FieldTrip data file where
data.label % shows channel names, in this case two channels were recorded GPiL02 and GPiR02
data.trial % includes the rawdata for both channels:
size(data.trial{1})
data.fsample % tells us the sampling rate (always 250 Hz for Percept)
% we only deal with resting data so we only have one trial. Let's extract
% the raw data of channel GPiR02 (the second in data.label) for which the ecg detection worked so
% well:
rawdata_gpir02 = data.trial{1}(2,:);

figure,plot(rawdata_gpir02)
% additionally the cleaned signal was stored in data.ecg_cleaned
cleandata_gpir02 = data.ecg_cleaned(2,:);

figure,plot(rawdata_gpir02), hold on, plot(cleandata_gpir02)

% You can get this ecg cleaning directly for any time-series you plugin to
% perceive_ecg:
data_to_be_cleaned = rawdata_gpir02;
sampling_rate = data.fsample;
ecg = perceive_ecg(data_to_be_cleaned,sampling_rate);

% The output structure contains some information:
ecg.cleandata % is the final result fromt he process
ecg.nandata % has the QRS complex replaced with nan
ecg.hr % stores the heartrate
ecg.stats % stores more information on the artefacts

%% Now please try to run this again for the other channel in the data variable (GPiL02):


