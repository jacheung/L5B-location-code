%% load data 
cell_num = 1;
T_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Raw\Interneurons\TArrays\';
S_directory = 'C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays\';
ST_directory = 'C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SpikeArrays\';
load([T_directory 'trial_array_' num2str(cell_num) '.mat'])
sweep_name = ['sweepArray_' T.mouseName '_' T.sessionName '_' T.cellNum '_' T.cellCode '.mat'];
spikes_trials = ['spikes_trials_' T.mouseName '_' T.sessionName '_' T.cellNum '_' T.cellCode '.mat'];
load([S_directory sweep_name])
load([ST_directory spikes_trials])
waveforms = s.get_spike_waveforms_all;
%%
% spike filter to grab indices of filtered spikes 
st_idx = cell(1,length(s.sweeps));
for i = 1:length(s.sweeps)
    sts = spikes_trials.spikesTrials{i}.spikeTimes;
    raw_sts = s.sweeps{i}.spikeTimes;
    [~,st_idx{i}] = intersect(raw_sts,sts);
end

all_waveforms = cell2mat(cellfun(@(x,y) normalize_var(x.spikeWaveforms{2}(:,y),-1,1),waveforms,st_idx,'uniformoutput',0));
avg_waveform = mean(all_waveforms,2);

figure(10);clf
plot(all_waveforms- mean(avg_waveform),'color',[.9 .9 .9],'linewidth',.1)
hold on; plot(avg_waveform - mean(avg_waveform),'r')



