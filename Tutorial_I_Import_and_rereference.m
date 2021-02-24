%% FIRST PART OF THE ICN INVASIVE ELECTROPHYSIOLOGY TUTORIAL 
% This tutorial requires the following repositories: 
% spm: https://github.com/spm/spm12
% wjn_toolbox: https://github.com/neuromodulation/wjn_toolbox
% FieldTrip: https://github.com/fieldtrip/fieldtrip.git
% Additionally, you will need to install the MatLab SignalProcessing
% Toolbox
% Written February 8th 2021 by the ICN group @
% https://github.com/neuromodulation/ephys_tutorials 

%% CLEAR ALL CLOSE ALL CLC 
clear all, close all, clc
%% SET FOLDERS AND PATH
% this is where you are going to work, choose any folder of your liking, we
% will create some files here:
root = 'G:\EphysTutorial'; % This needs to be adapted!
% let's create the directory you specified as root and move there
if ~exist(root,'dir')
    mkdir(root)
end

cd(root) 

% Now add the required code to your path this is where your code 
% repositories for wjn_toolbox and SPM are.
addpath('G:\EphysTutorial\code\wjn_toolbox') % This needs to be adapted
addpath('G:\EphysTutorial\code\spm12') % This needs to be adapted
% SPM needs to be initialized to include FieldTrip folders in path:
% Do not initialize FieldTrip yet, as MatLab will get confused!
spm('defaults','eeg')

%% Let's create some data first
% Let's simulate some super valuable invasive neurophyisology data with 
% Fieldtrip first to have something to play around with, just for the fun 
% of it.
% Clean slate first:
clear all, close all, clc
info.fs = 280;
info.channels = {'STNR0','L_STNR1','LFP_R2_STN','STNR3','STNL0','STNL1','LFP_STN_L2','STNL3','ECOG_SML_1','ECOG_SML_2','ECOG_L_SM_3','ECOGL_SM4','EEG_C3','EEG_Cz','FDI_R01EMG','FDIR_02'};
info.info = 'Some super valuable invasive resting datasets that were not taken care off for a long time. The datset was recorded 1985/12/16 with TMSi SAGA, maybe like 13:35 h in the afternoon or something. 58 years old tremor dominant PD patient, her preop UPDRS was 32/12 OFF/ON, disease duraation like ~8 years,  Medication and stimulation were OFF. Sampling interval was 3.57 ms. 3389. Reference was STNL3.';
writetable(struct2table(info),'trashy_info.txt')
nsignals = length(info.channels);
% Don't worry if you don't understand this paragraph, this is just to
% create some data to work with. Have a look at ft_freqsimulation if you
% are further interested in data simulation.

cfg=[];
cfg.method = 'broadband';
cfg.fsample = 280;
cfg.trllen = 60;

raw_signals=[];
for a = 1:nsignals
    cfg.n1.ampl     = 100+randi(30);
    cfg.n1.bpfreq   = [randi(3)+12 randi(5)+15];
    cfg.n2.ampl     = 50+randi(50);
    cfg.n2.bpfreq   = [randi(5)+20 randi(5)+30];
        cfg.noise.ampl = randi(50);
    data = ft_freqsimulation(cfg);
    raw_signals(a,:) = data.trial{1}(1,:);
    if a == 8 
        raw_signals(a,:) = rand(1,size(raw_signals,2));
    end
end


% so 16 raw data time-series are now stored in 
writetable(array2table(raw_signals'),'forgotten_trashy_dataset.csv')

% You can now see two files written to this folder, a csv and a text file.
% Both are human readable and can be opened with any text editor. Have a
% look at them.
%% Now this is where I would usually start importing data that I got from 
%% collaborations (except when shared by Esther).

% Always start with a clean slate
clear all, close all, clc

% We need to find the three main information required for neurophysiology
% analysis 1) raw data, 2) sampling rate and 3) channel names.

% Additional very important information is the question of the reference 
% used during recording, but this information can be lost at times. 

% first let's try and load the raw time-series data:
% note that we are converting to array for convenience
rawdata = table2array(readtable('forgotten_trashy_dataset.csv')); 
figure
plot(rawdata)

% Problem: we don't know the sampling rate yet, so we don't know how many
% samples constitute one second. We also don't know how long the overall
% recording time was. 

% I wrote a plotting function that I always use, it aligns and z-scores the
% data to make them more accessible:
% the first input is the sampling rate, the second the rawdata and the
% third optional input are the channel names.
% Let's have a look at the effect of differing sampling rates on the time
% axis:
figure
subplot(1,2,1)
% E.g. if sampling rate was believed to be 2 kHz...
false_sampling_rate1 = 2000;
wjn_plot_raw_signals(2000,rawdata)
title(['Sampling rate: ' num2str(false_sampling_rate1)])
subplot(1,2,2)
% ... this is different when compared to 50 Hz
false_sampling_rate2 = 50;
wjn_plot_raw_signals(false_sampling_rate2,rawdata)
title(['Sampling rate: ' num2str(false_sampling_rate2)])

% So let's find the real sampling rate in the info file.
% Typical abbreviations for sampling rate (also called sampling frequency
% are), fs, fsample, samplingrate,sr.
% Sometime the sampling interval is stated. This is the amount of time
% between two samples. 

info = readtable('trashy_info.txt');
% Let's have a look at the variables in the info table:
disp(info.Properties.VariableNames')

% fs looks just right so let's use that variable:
fsample = info.fs; % luckily this information is readily available. 

% If you looked closely there was also an "info" field. If you look into 
% that info text:
disp(info.info)
% it also states a sampling interval of 3.57 ms, which would translate to
% 1000/3.57 = 280.112 Hz, which is close but not 100% accurate.

figure
wjn_plot_raw_signals(fsample,rawdata)
xlabel('Time [s]')
title({['Sampling rate: ' num2str(fsample)],'Original recording time was one minute.'})

% From the  
disp(info.Properties.VariableNames')
% it seems that the channel names are in the info table depicted by
% channel_ and a number. I have written a very handy string indicator that 
% I often use. This can be used to identify the relevant fields in the 
% table:
channel_indices = ci('channels',info.Properties.VariableNames);
channel_names = table2cell(info(:,channel_indices)); % note that we are converting from Table to cell for convenience later on

% So we can look at the data now with the most important aspects in place
figure
wjn_plot_raw_signals(fsample,rawdata,channel_names)

% Let's also extract that text info field from the table:
additional_info = info.info;
% Great! Now we have all the important information extracted for conversion 
% to electrophysiology toolboxes like FieldTrip or SPM

%% CONVERTING TO FIELDTRIP:
% always start with a clean slate
clear all, close all, clc
% Let's repeat what we did up there.
rawdata = table2array(readtable('forgotten_trashy_dataset.csv'));
info = readtable('trashy_info.txt');
fsample=info.fs;
channel_names = table2cell(info(:,ci('channel',info.Properties.VariableNames)));
additional_info = info.info;
% Converting to FieldTrip is the easiest possible step
% You just need to write this data structure from the info you have
% Note that data have to be in correct dimensions:
% Number_of_channels x Number_of_samples
% Therefore we will have to transpose our data using the ' function
% It can be useful to convert single precision data to double and it
% doesn't hurt if its double already:
data.trial{1} = double(rawdata'); % if you have just 1 trial this is correct

% now given that we have the sampling rate we can construct the time axis
% using linear interpolation with linspace
data.time{1} = linspace(0,length(rawdata)/fsample,length(rawdata));

% channel names for FieldTrip go in the label field as a cell array:
data.label = channel_names;

% the sampling rate must be stored in the data.fsample field
data.fsample = fsample;

% we can save additional info in this structure, it doesn't hurt as long as
% you do not choose fieldnames required by FieldTrip. 
data.addition_info = info;

% now we can save this FieldTrip dataset to the disk and load it in 
% MNE MatLab or Python, SPM, BrainStorm and many other Toolboxes.
save('manual_fieldtrip_dataset.mat','data')

% to save me some time, I have created a one line wrapper for this process:
% data = wjn_import_rawdata(filename,idata,chanlabels,fs,format,additional_info)
% format is optional, defaults to 'spm', but can also be 'fieldtrip'
% This function also save the data to disk:
data = wjn_import_rawdata('oneline_fieldtrip_dataset.mat',rawdata,channel_names,fsample,'fieldtrip',info);

% you can now use any FieldTrip function for further processing, e.g. the
% databrowser for visual artifact inspection. 
cfg=[];
cfg.blocksize = 10;
cfg.continuous  = 'yes';
cfg.ylim = [-100 100];
cfg.viewmode ='vertical';
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 45;
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 3;
ft_databrowser(cfg,data);

%% But maybe we can improve the dataset even before importing. 
clear all, close all,clc
% Let's delete what we have saved in the previous part
delete('*fieldtrip*')
rawdata = table2array(readtable('forgotten_trashy_dataset.csv'));
info = readtable('trashy_info.txt');
fsample=info.fs;
channel_names = table2cell(info(:,ci('channel',info.Properties.VariableNames)));
source_info = info.info;

% Let's focus on the channel names first. They don't seem to follow a
% convention. Luckily we have written a script that makes sense out of what
% we expect (at least in our clinic):

chanlabels = wjn_channel_converter(channel_names);
disp(chanlabels')

% That looks much better, what is the convention then ?
% Let's have a look at the first channel
disp(chanlabels{1})
% 1) start with recording type: LFP
% 2) take any number that is given: 01
% 3) identify the hemisphere or side of the recording: L(eft)
% 4) attribute the recording location from channel name: STN
% 5) have an additional identification of the manufacturer, which can in
% our clinic sometimes be inferred from the number of LFP channels in the
% set : MT (Medtronic). Medtronic is the only manufacturer left by whom 4
% contact DBS electrodes would be implanted.
% If you look closely at the info text, there is a 3389 hidden in there.
% That is the traditional Medtronic DBS lead. The script cannot know that but it is easy 
% to replace the '_U' for unknown in the names with _MT for Medtronic.

chanlabels(1:8) = strrep(chanlabels(1:8),'_U','_MT');



% What else is there in the info?
disp(info.info)

% It seems the dataset was recorded on a very special day, December 16th
% 1985 at 1:35:23 pm, so this is nice to know, as we may want to store that
% information in case we have more than one dataset from that subject.
additional_info.recording_date = datetime(1985,12,16,13,35,23,100);
% The reference is noted as channel STNL3. That is very important
% information for the preprocessing. Let's translate that into our
% convention:
additional_info.reference = 'LFP_3_L_STN_MT';

% also it is mentioned that it is a resting dataset 
% and both medication and stimulation states are mentioned to be OFF.
additional_info.task = 'rest';
additional_info.medication = 'MedOff';
additional_info.stimulation = 'StimOff';
% since this is the first dataset we get from a series, we will simply call
% this subject 001
additional_info.subject = 'sub-001';
% My traditional way to handle this set from here was to convert to SPM
% would be to pseudostandardize some kind of internal convention for the
% filename and condition description:
additional_info.source_info = source_info;
filename = [additional_info.subject '_' additional_info.medication '_' ...
    additional_info.stimulation '_' char(datetime(additional_info.recording_date,'Format','dd-MMM-yyyy')) '.mat'];
D=wjn_import_rawdata(['spm_' filename],rawdata,chanlabels,fsample,'spm',additional_info);
data=wjn_import_rawdata(['fieldtrip_' filename],rawdata,chanlabels,fsample,'fieldtrip',additional_info);
%% Now we have a neat datafile in both SPM and FieldTrip format for further analyis.
clear all, close all,clc
D=spm_eeg_load('spm_sub-001_MedOff_StimOff_16-Dec-1985.mat')
% all the info we gathered can be found in 
D.info
% the SPM data object stores all the remaining variables and data in
% defined parts:
D
disp(D)

figure
wjn_plot_raw_signals(D.time,D(:,:),D.chanlabels)

% Ok, this looks nice but what about the reference?
% If you look closely, channel 8 LFP_3_L_STN_MT seems to be empty. 
% Now thinking about it, it was mentioned to be used as the reference. 
% So likely this channel was not connected. At least it only carries noise. 

figure
wjn_plot_raw_signals(D.time,D([4 8],:),D.chanlabels([4 8]))


% Thus, before we would want to analyze the data, we should consider the
% reference used.


% All data were referenced to an LFP contact. That means that all LFP channels are
% contaminated by LFP signal from the reference contact. So we should rereference.
%
% Optimally, rereferencing is conducted to create local bipolar
% derivatives. This should be the first preprocessing step that you do. It
% requires some understanding of what is going on. Best is to start with 
% a new data array where you add your bipolar derivations. 
% To get bipolar data, you simply subtract the data from two channels. 
% Let's do that for the adjacent channels of the right LFP electrode.
% We can use the help of the channel indicator here:
index_lfpr = ci('R_STN',D.chanlabels);
lfp_r_bp = D(index_lfpr(1:end-1),:)-D(index_lfpr(2:end),:);
% we can write out the channel names by hand:
channels_lfp_right = {'LFP_R_01_STN_MT','LFP_R_12_STN_MT','LFP_R_23_STN_MT'};
% now we have three channels from adjacent contact pairs 
% For the right side it is a little more complicated. The reference was
D.info.reference
index_lfpl = ci('L_STN',D.chanlabels);
% so that means LFP_03_L_STN_U is already hardware bipolar:
lfp_l_bp = [D(index_lfpl(1:2),:)-D(index_lfpl(2:3),:);D(index_lfpl(3),:)];
channels_lfp_left = {'LFP_L_01_STN_MT','LFP_L_12_STN_MT','LFP_L_23_STN_MT'};

% Next, we will rereference ECOG channels. ECOG is probably least affected 
% because of the high amplitude signals. Nevertheless, we can improve the 
% spatial specificity of our analyses. For ECOG, two potential referencing 
% schemes are the most useful, which we will both perform here:s
% 1) is common average referencing (CAR):
index_ecog = ci('ECOG',D.chanlabels);
ecog_car = D(index_ecog,:)-nanmean(D(index_ecog,:),1);
channels_ecog_car = strcat(D.chanlabels(index_ecog)','_CAR')';
% 2) is bipolar referencing like the LFP above:
ecog_bp = D(index_ecog(1:end-1),:)-D(index_ecog(2:end),:);
% we can write out the channel names by hand:
channels_ecog_bp = {'ECOG_L_12_SM_U','ECOG_L_23_SM_U','ECOG_L_34_SM_U'};

% next up, the EEG channels. Those are also relatively unaffected by the
% LFP because it is an order of magnitude smaller than the EEG. However, to
% rule out contamination we can create an additional more local bipolar
% derivative from the EEG channels:
index_eeg = ci('EEG',D.chanlabels);
eeg_bp = D(index_eeg(1),:)-D(index_eeg(2),:);
channels_eeg = {'EEG_L_C3Cz_U'};
% Finally, we do the same for the EMG:
index_emg = ci('EMG',D.chanlabels);
emg_bp = D(index_emg(1),:)-D(index_emg(2),:);
channels_emg = {'EMG_12_R_FDI_U'};

% All those data can now be added to the original dataset, or stored in a
% new one:
derivative_data = [lfp_r_bp;lfp_l_bp;ecog_car;ecog_bp;eeg_bp;emg_bp];
derivative_chanlabels = [channels_lfp_right,channels_lfp_left,channels_ecog_car,channels_ecog_bp,channels_eeg,channels_emg];
D=wjn_add_channels(D.fullfile,derivative_data,derivative_chanlabels)
% The above function has created a new SPM file
% 'aspm_sub-001_MedOff_StimOff_16-Dec-1985.mat' now including:
D.nchannels
D.chanlabels'

% Finally we can remove the empty and noisy monopolar STNL3 channel. 

D=wjn_remove_channels(D.fullfile, {'LFP_3_L_STN_MT'   });

% Now we can move forward, e.g. with filtering and preprocessing: 
Df=wjn_filter(D.fullfile,1,'high');
Dff=wjn_filter(Df.fullfile,[48 52],'stop');Df.delete;
D=wjn_filter(Dff.fullfile,98,'low');Dff.delete;

% You can rename it using the spm 'move' command:
fname = D.fname;
D=D.move(['reref_spm_' fname(10:end)]);
% Or you can swap that to fieldtrip:
data = wjn_raw_spm2fieldtrip(D.fullfile);
% don't forget to take your info variable:
data.info = D.info;
fname = D.fname;
save(['reref_fieldtrip_' fname(7:end)])

%% DONE! This is where I would start my data analysis with this dataset:
clear all, close all, clc
D=spm_eeg_load('reref_spm_sub-001_MedOff_StimOff_16-Dec-1985.mat');

figure
wjn_plot_raw_signals(D.time,D(:,:),D.chanlabels)

channel_of_interest=D.indchannel('LFP_R_01_STN_MT');

% skip this if you don't have the signal processing toolbox
Dtf=wjn_tf_wavelet(D.fullfile,1:100,20,channel_of_interest);

figure
subplot(4,1,1)
wjn_plot_raw_signals(D.time,D(channel_of_interest,:),D.chanlabels(channel_of_interest))
title('Raw data');
subplot(4,1,2);
mypower(Dtf.frequencies,squeeze(Dtf(1,:,:,1)))
ylabel('PSD');xlabel('Frequency [Hz]');title('Power spectrum')
subplot(4,1,[3 4])
surface(Dtf.time,Dtf.frequencies,squeeze(log(Dtf(1,:,:,1))),'edgecolor','none')
title('Time Frequency Analysis');view(-50,70);caxis([7 15]);zlim([5 20]);
xlabel('Time [s]');ylabel('Frequency [Hz]');zlabel('PSD');figone(25,20)
%% CONVERSION TO BIDS (skip this if you did not want to download fieldtrip)
% Now as part of the INF team, I will also take the opportunity to promote
% the use of the BIDS standard (https://bids-specification.readthedocs.io/en/stable/)
% and demonstrate how simple it can be.
% For that we need some additional code from the original FieldTrip toolbox
clear all, close all, clc
load('fieldtrip_sub-001_MedOff_StimOff_16-Dec-1985')
addpath G:\EphysTutorial\code\fieldtrip % This needs to be adjusted!

% let's define some general settings
cfg = [];
cfg.method    = 'convert';
cfg.datatype  = 'ieeg';
cfg.bidsroot  = 'bids';
cfg.sub       = '001'; % subject identifier
cfg.task      = 'Rest';
cfg.ses       = 'EphysMedOff';
cfg.acq       = 'StimOff';

% Info to build your participants table using some of the info provided:
disp(data.info.source_info)

cfg.participants.age = 58;
cfg.participants.sex = 'f';
cfg.participants.updrs_off = 32;
cfg.participants.updrs_on = 12;
cfg.participants.updrs_timepoint = 'preoperative';
cfg.participants.disease_duration = 8;
cfg.participants.clinical_subtype = 'tremor dominant';
cfg.participants.DBS_electrode_type = 'Medtronic 3389';
cfg.participants.ECOG_electrode_type = 'AdTech 6 contact Narrow Body';
cfg.participants.hardware_amplifier = 'TMSi SAGA';

% specify the information for the scans.tsv file
% This is optional, you can also pass other pieces of info
cfg.scans.acq_time = char(data.info.recording_date);
% specify some general information that will be added to the eeg.json file
cfg.InstitutionName             = 'Charite - Universitaetsmedizin Berlin';
cfg.InstitutionalDepartmentName = 'Sektion fuer Bewegungsstoerungen und Neuromodulation';
cfg.InstitutionAddress          = 'Chariteplatz 1, 10117 Berlin';
cfg.dataset_description.Authors  = 'Gerd-Helge Schneider, Wolf-Julian Neumann, Andrea A. Kuehn';
% provide the mnemonic and long description of the task
cfg.TaskName        = 'Rest';
cfg.TaskDescription = 'Subjects were asked not to move and remain eyes open. Recording was performeds after withdrawal of medication.';

% these are iEEG specific
cfg.ieeg.PowerLineFrequency = 50;   % since recorded in the Europe
cfg.ieeg.iEEGReference       = 'LFP_L_3_STN_MT'; % as stated in info
cfg.ieeg.SEEGChannelCount = numel(ci('SEEG',data.hdr.chantype));
cfg.ieeg.ECOGChannelCount = numel(ci('ECOG',data.hdr.chantype));
cfg.ieeg.SoftwareFilters  = 'High-pass filter 1 Hz,  Stop-Band Filter 48-52 Hz, Low-pass filter 98 Hz,';
cfg.ieeg.RecordingDuration  = [num2str(data.time{1}(end)) ' sec'];
cfg.ieeg.RecordingType = 'Continuous';
cfg.ieeg.ElectrodeManufacturer     =  'Medtronic';
cfg.ieeg.ElectrodeManufacturersModelName = '3389';
data2bids(cfg,data);

% Feel free to have a look at this beautiful BIDS dataset in the suggested 
% universally readable BrainVision format.





