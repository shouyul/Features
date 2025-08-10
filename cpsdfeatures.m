function [ anglephase, cpsdmag,  cpsdmagt, cpsdmagb] = cpsdfeatures( data, par, minfeat, selectedchans  )
%Calculates CPSD features for alpha theta and beta frequency ranges and
%their subbands using the function crosspower.m

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
% anglephase = Cross Power Spectral Density phase (imag) is the phase delay
% of cross power spectral density between channel pairs
% cpsdmag = Cross Power Spectral Density magnitude (abs)is the magnitude of
% the cross power spectral density between channel pairs
% cpsdmagt = CPSD for theta
% cpsdmagb = CPSD for beta

thetarange = [4:0.5:7]; %in hz 
alpharange = [7.5:0.5:13.5]; %in hz
betarange = [13:0.5:30]; %in hz
if exist('selectedchans','var')
    for trial = 1:size(data,2)
    data{trial} = data{trial}(:,(selectedchans));
    end
end
[cpsdmag0, cpsdmag1, cpsdmag2,cpsdmag3,cpsdmag4,cpsdmag5, cpsdphase1,cpsdphase2,cpsdphase3,cpsdphase4,cpsdphase5,...
    anglephase1,anglephase2,anglephase3,anglephase4,anglephase5,...
    cpsdmagt0, cpsdmagt1, cpsdmagt2, cpsdmagb0,cpsdmagb1,cpsdmagb2, cpsdmagb3,cpsdmagb4, ~, ~, ~] = cellfun(@(x) ...
        crosspower(x, minfeat, par.currentsamp, par.halfsamp, par.samp4overlap, thetarange,alpharange,betarange), data, 'un', 0);

cpsdmag = {cell2mat(cpsdmag0'), cell2mat(cpsdmag1'),cell2mat(cpsdmag2'),cell2mat(cpsdmag3'),cell2mat(cpsdmag4'),cell2mat(cpsdmag5')};
cpsdphase = {cell2mat(cpsdphase1'), cell2mat(cpsdphase2'), cell2mat(cpsdphase3'), cell2mat(cpsdphase4'), cell2mat(cpsdphase5')};
anglephase = {cell2mat(anglephase1'), cell2mat(anglephase2'), cell2mat(anglephase3'), cell2mat(anglephase4'), cell2mat(anglephase5')};

%CPSD theta/beta
cpsdmagt = {cell2mat(cpsdmagt0'), cell2mat(cpsdmagt1'),cell2mat(cpsdmagt2')};
cpsdmagb = {cell2mat(cpsdmagb0'), cell2mat(cpsdmagb1'),cell2mat(cpsdmagb2'),cell2mat(cpsdmagb3'),cell2mat(cpsdmagb4')};


end

