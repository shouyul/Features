function [ mscoh ] = mscoherefeatures( data, par, selectedchans )
%Calculates magnitude squared coherence features between channel pairs for
%alpha frequency and alpha subbands

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

%OUTPUT
%mscoh = alpha magnitude squared coherence calculated with the function magsquarecoherence.m

if exist('selectedchans','var')
    for trial = 1:size(data,2)
    data{trial} = data{trial}(:,(selectedchans));
    end
end
[mscohmag0, mscohmag1, mscohmag2,mscohmag3,mscohmag4,mscohmag5, ~] = cellfun(@(x) magsquarecoherence(x, par.currentsamp, par.halfsamp, par.samp4overlap), data, 'un', 0);
mscoh = {cell2mat(mscohmag0'),cell2mat(mscohmag1'),cell2mat(mscohmag2'),cell2mat(mscohmag3'),cell2mat(mscohmag4'),cell2mat(mscohmag5')};

end

