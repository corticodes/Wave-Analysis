%%%% See 
%%%% /media/sil2/Literature/Projects/corplex/progress reports/meetings/210802/script210802.m
%%%% for automations


%clear all
fileDir='/media/sil2/Literature/Projects/corplex/progress reports/meetings/210802/Simulation/';
nCh=200;
LFP_sample_ms=0.1;
spikes_sample_ms=0.01;
smooth_window_ms=10; %ms
spikeRateSmoothing=smooth_window_ms/spikes_sample_ms; %100ms window

LFP_sampling_freq=1000/LFP_sample_ms; %samples/s
spikes_sampling_freq=1000/spikes_sample_ms; %samples/s

En=reshape(1:200,10,20)';



%open files and retrive data

LFP = readmatrix([fileDir 'lfp.txt'],'NumHeaderLines',0);
lfp_T=(1:size(LFP,2))*LFP_sample_ms; %sampling times in ms

FD = lowpass(LFP',5,LFP_sampling_freq,'ImpulseResponse','iir').';

HT=hilbert(FD').';
HTabs=abs(HT);
HTangle=angle(HT);

%Find time to half max
ch_max=max(FD,[],2);
[~,time_to_half_max]=max(FD(:,1:end-1)-ch_max/2<0 & FD(:,2:end)-ch_max/2>0,[],2);
figure
[hCbar,h]=IntensityPhysicalSpacePlot(1:nCh,time_to_half_max-min(time_to_half_max),En,'plotElectrodeNumbers',0);
title('Time to half max')


spikesBin = readmatrix([fileDir 'Spikes.txt']);
spike_T=(1:size(spikesBin,2))*spikes_sample_ms; %sampling times in ms

[neurons,spike_samples]=find(spikesBin);
spike_times=spike_T(spike_samples);



% Calculate ALSA
spikingRate = binSpikes2fireRate(spikesBin,spikes_sampling_freq,'slidingWindowSize',spikeRateSmoothing);
spikingRateDataFormatSmoothed=smoothdata(spikingRate','gaussian',spikeRateSmoothing)';
averageSpikeRate=mean(spikingRateDataFormatSmoothed,1);

ALSA=calcALSA(spikingRateDataFormatSmoothed,'En',En);


for i=1:nCh
    [~,chlocs]=findpeaks(ALSA(i,:));
    if isempty(chlocs)
        ALSA_Locs(i)=NaN;
    else
        ALSA_Locs(i)=chlocs(1);
    end
end

% startEndWave=[1 size(ALSA,2)];
startEndWave=[1 5200];
chSTDs=std(ALSA,0,2);
nSTD=2;
for i=1:nCh
    chLoc=find(ALSA(i,startEndWave(1):startEndWave(2)-1)<nSTD*chSTDs(i)&ALSA(i,(startEndWave(1)+1):startEndWave(2))>=nSTD*chSTDs(i),1);
    if isempty(chLoc)
        ALSA_Locs(i)=NaN;
    else
        ALSA_Locs(i)=chLoc;
    end
end
hasALSA=find(~isnan(ALSA_Locs));


figure
plot(spike_samples,neurons,'k.')
xlabel('Time [samples]')
ylabel('Neuron #')
hold on
plot(ALSA_Locs,1:nCh,'or')

figure
[hCbar,h]=IntensityPhysicalSpacePlot(hasALSA,ALSA_Locs(hasALSA)-min(ALSA_Locs(hasALSA)),En,'plotElectrodeNumbers',0);
title('Time to ALSA Onset')

save([fileDir 'Spikes and LFP clust_con_0.03_conectivity_const_0.4_0.4_0.3_0.1.m'],'LFP','spikesBin','lfp_T','spike_T','spike_times')
