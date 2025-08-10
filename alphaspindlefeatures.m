function [sumspind,diffspind,numspindperchanall] = alphaspindlefeatures(data, par, selectedchans)
%Calculates alpha spindle features with the option to calculate all
%features or specific features (based on important channel)

%INPUTS
%data = cell array 1*number of trials; each cell is a separate trial, that
    %  contains a #channels*length of trial array (e.g. 1*25 cell vector, each cell = 16*1282)
%par = structure containing parameters necessary for computing alpha spindles (see ComputeImageryFeatures)
    %e.g.
    % par.currentsamp = 256; %Sampling rate of data
    % par.halfsamp = par.currentsamp/2; %Size of sliding window
    % par.samp4overlap = round(0.7813*par.halfsamp); %par.samp4overlap NEEDS TO BE AN INTEGER!!! %overlap of window
    % par.freqrange = [1:0.5:40]; %frequency range of interest
    % par.numstepperHz = 2; %how many steps per Hz in frequency range
    % par.SNRthreshold = 2; %minimum signal to noise ratio for alpha spindle condition
    % par.tau = 2; %downsample factor (for approximate and sample entropy to save time)
%selectedchans = optional input for only computing alpha spindles for specific channels    

%OUTPUTS
%sumspind = added number of detected alpha spindles for entire alpha range
    %(7.5-13.5 Hz) across all possible channel pairs (ideal for identifying
    %hemispheric imbalances for Rest)
%diffspind = difference of number of detected alpha spindles across each
    %channel pair for alpha subbands [1 (7.5-9.5Hz), 2
    %(8.5-10.5Hz),3(9.5-11.5Hz),4(10.5-12.5Hz) and 5(11.5-13.5Hz)]
    %(ideal for identifying hemispheric imbalances for Imagery)
%numspindperchanall = absolute number of alpha spindles per channel for
    %frequency range 0 (7.5-13.5 Hz),1 (7.5-9.5Hz), 2 (8.5-10.5Hz),3(9.5-11.5Hz),4(10.5-12.5Hz) and 5(11.5-13.5Hz)
    %this feature is computed with the function newalphaspindle

if exist('selectedchans','var')
    for trial = 1:size(data,2)
    data{trial} = data{trial}(:,unique(selectedchans));
    end
end

[newSegmentPass0, newSegmentPass1, newSegmentPass2, newSegmentPass3, newSegmentPass4, newSegmentPass5, ~] = ...
    cellfun(@(x) arrayfun(@(y)newalphaspindle(x(:,y), par.currentsamp, par.halfsamp, par.samp4overlap, par.freqrange, ...
    [7.5:13.5],par.numstepperHz, par.SNRthreshold),...
    1:size(x,2), 'un', 0), data, 'un', 0);

numspindperchan0 = cell2mat(cellfun(@(x) cell2mat(arrayfun(@(y) sum(x{y}), 1:length(x), 'un', 0)), newSegmentPass0, 'un', 0)');        
numspindperchan1 = cell2mat(cellfun(@(x) cell2mat(arrayfun(@(y) sum(x{y}), 1:length(x), 'un', 0)), newSegmentPass1, 'un', 0)');
numspindperchan2 = cell2mat(cellfun(@(x) cell2mat(arrayfun(@(y) sum(x{y}), 1:length(x), 'un', 0)), newSegmentPass2, 'un', 0)');
numspindperchan3 = cell2mat(cellfun(@(x) cell2mat(arrayfun(@(y) sum(x{y}), 1:length(x), 'un', 0)), newSegmentPass3, 'un', 0)');
numspindperchan4 = cell2mat(cellfun(@(x) cell2mat(arrayfun(@(y) sum(x{y}), 1:length(x), 'un', 0)), newSegmentPass4, 'un', 0)');
numspindperchan5 = cell2mat(cellfun(@(x) cell2mat(arrayfun(@(y) sum(x{y}), 1:length(x), 'un', 0)), newSegmentPass5, 'un', 0)');

numspindperchanall = {numspindperchan0,numspindperchan1,numspindperchan2,numspindperchan3,numspindperchan4,numspindperchan5};
numspindperchan0 = {numspindperchan0'};
numspindperchan = {numspindperchan1',numspindperchan2',numspindperchan3',numspindperchan4',numspindperchan5'};

if size(data{1},2)>1
    diffspind = cellfun(@(x) subchan(x), numspindperchan, 'un', 0);
    sumspind = cellfun(@(x) addchan(x), numspindperchan0, 'un', 0);
else
    [diffspind, sumspind] = deal(0);
end

end

