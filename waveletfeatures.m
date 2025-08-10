function [ difflow,diffhigh,diffmid, sumlow, sumhigh,summid,sumpower,tsumlow,tsumhigh,tsummid,tsumpower,...
            bsumlow,bsumhigh,bsummid,bsumpower, ApproxEntropy, SampleEntropy,cwtpower,tcwtpower,bcwtpower] = ...
            waveletfeatures(data,par, minfeat, selectedchans, selectedFeatureGroup)
%Calculates wavelet features for trials based on continuous wavelet
%transform into scalograms

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
%selectedchans = optional input for only computing alpha spindles for
%specific channels NOTE: can compute this with function decodefeaturechans.m
%selectedFeatureGroup = only compute specific feature groups NOTE: can compute this with function decodefeaturechans.m

%OUTPUTS (computed with function contwavetran.m)
%NOTE: these features are computed for entire specified frequency bands
%(identified with the number 0 in the variable name), as well as subbands
%(identified with the numbers 1-5, for first subband, second subband, etc.)
%difflow = subtracted alpha cwt pseudopower from channel pairs for "low"-energy areas (i.e. below 1 SD of mean cwt power)
%diffhigh = subtracted alpha cwt pseudopower from channel pairs for "high"-energy areas (i.e. above 1 SD of mean cwt power)
%diffmid = subtracted alpha cwt pseudopower from channel pairs for "medium"-energy areas (i.e. within 1 SD of mean cwt power)
%sumlow =added alpha cwt pseudopower from channel pairs for "low"-energy areas (i.e. below 1 SD of mean cwt power)
%sumhigh = added alpha cwt pseudopower from channel pairs for "high"-energy areas (i.e. above 1 SD of mean cwt power)
%summid  = added alpha cwt pseudopower from channel pairs for "medium"-energy areas (i.e. within 1 SD of mean cwt power)
%sumpower  = added alpha cwt pseudopower from channel pairs for whole cwt
%tsumlow  = added theta cwt pseudopower from channel pairs for "low"-energy areas (i.e. below 1 SD of mean cwt power)
%tsumhigh  = added theta cwt pseudopower from channel pairs for "high"-energy areas (i.e. above 1 SD of mean cwt power)
%tsummid  = added theta cwt pseudopower from channel pairs for "medium"-energy areas (i.e. within 1 SD of mean cwt power)
%tsumpower  = added theta cwt pseudopower from channel pairs for whole cwt
%bsumlow  = added beta cwt pseudopower from channel pairs for "low"-energy areas (i.e. below 1 SD of mean cwt power)
%bsumhigh  = added beta cwt pseudopower from channel pairs for "high"-energy areas (i.e. above 1 SD of mean cwt power)
%bsummid  = added beta cwt pseudopower from channel pairs for "medium"-energy areas (i.e. within 1 SD of mean cwt power)
%bsumpower  = added beta cwt pseudopower from channel pairs for whole cwt
%ApproxEntropy  = approximate entropy calculated for each channel
%SampleEntropy  = sample entropy calculated for each channel
%cwtpower  = absolute alpha cwt pseudopower each channel
%tcwtpower  = absolute theta cwt pseudopower each channel
%bcwtpower  = absolute beta cwt pseudopower each channel

if exist('selectedchans','var')
    for trial = 1:size(data,2)
    data{trial} = data{trial}(:,unique(selectedchans));
    end
end

if ~exist('selectedFeatureGroup', 'var')
    selectedFeatureGroup = [];
end

[difflow1, difflow2, difflow3, difflow4, difflow5 , diffhigh1, diffhigh2, diffhigh3, diffhigh4, diffhigh5 ,diffmid1, diffmid2, diffmid3, diffmid4, diffmid5 ,...
    sumlow1, sumlow2, sumlow3, sumlow4, sumlow5 ,sumhigh1, sumhigh2, sumhigh3, sumhigh4, sumhigh5 ,summid1, summid2, summid3, summid4, summid5,  ...   
    sumpower0, sumpower1, sumpower2, sumpower3 ,sumpower4,sumpower5, cwtpower0, cwtpower1, cwtpower2, cwtpower3, cwtpower4, cwtpower5,...
    tcwtpower0, tcwtpower1, tcwtpower2, tsumlow1, tsumlow2,tsumhigh1, tsumhigh2, tsummid1, tsummid2, tsumpower0, tsumpower1, tsumpower2,...
    bcwtpower0, bcwtpower1,bcwtpower2,bcwtpower3,bcwtpower4, bsumlow1, bsumlow2, bsumlow3, bsumlow4,bsumhigh1, bsumhigh2, bsumhigh3, bsumhigh4, ...
    bsummid1, bsummid2, bsummid3, bsummid4,bsumpower0, bsumpower1, bsumpower2, bsumpower3 ,bsumpower4, ...
    ~,~, ~, AppEnt_normrec0, AppEnt_normrec1,AppEnt_normrec2,AppEnt_normrec3,AppEnt_normrec4,AppEnt_normrec5,...
    SampEnt_normrec0,SampEnt_normrec1,SampEnt_normrec2,SampEnt_normrec3,SampEnt_normrec4,SampEnt_normrec5,...
    WEnt_normrec0,WEnt_normrec1,WEnt_normrec2,WEnt_normrec3,WEnt_normrec4,WEnt_normrec5] =  cellfun(@(x) contwavetran(x,minfeat, par.currentsamp, par.tau, selectedFeatureGroup), data, 'un', 0);

    difflow = {cell2mat(difflow1'), cell2mat(difflow2'), cell2mat(difflow3'), cell2mat(difflow4'), cell2mat(difflow5')};
    diffhigh = {cell2mat(diffhigh1'), cell2mat(diffhigh2'), cell2mat(diffhigh3'), cell2mat(diffhigh4'), cell2mat(diffhigh5')};
    diffmid = {cell2mat(diffmid1'), cell2mat(diffmid2'), cell2mat(diffmid3'), cell2mat(diffmid4'), cell2mat(diffmid5')};

    sumlow = {cell2mat(sumlow1'), cell2mat(sumlow2'), cell2mat(sumlow3'), cell2mat(sumlow4'), cell2mat(sumlow5')};
    sumhigh = {cell2mat(sumhigh1'), cell2mat(sumhigh2'), cell2mat(sumhigh3'), cell2mat(sumhigh4'), cell2mat(sumhigh5')};
    summid = {cell2mat(summid1'), cell2mat(summid2'), cell2mat(summid3'), cell2mat(summid4'), cell2mat(summid5')};
    sumpower = {cell2mat(sumpower0'),cell2mat(sumpower1'),cell2mat(sumpower2'),cell2mat(sumpower3'),cell2mat(sumpower4'),cell2mat(sumpower5')};
    cwtpower = {cell2mat(cwtpower0'),cell2mat(cwtpower1'),cell2mat(cwtpower2'),cell2mat(cwtpower3'),cell2mat(cwtpower4'),cell2mat(cwtpower5')};

    %Continuous Wavelet Tranform theta/beta
    tsumlow = {cell2mat(tsumlow1'), cell2mat(tsumlow2')};
    tsumhigh = {cell2mat(tsumhigh1'), cell2mat(tsumhigh2')};
    tsummid = {cell2mat(tsummid1'), cell2mat(tsummid2')};
    tsumpower = {cell2mat(tsumpower0'),cell2mat(tsumpower1'),cell2mat(tsumpower2')};
    tcwtpower = {cell2mat(tcwtpower0'),cell2mat(tcwtpower1'),cell2mat(tcwtpower2')};

    bsumlow = {cell2mat(bsumlow1'), cell2mat(bsumlow2'), cell2mat(bsumlow3'), cell2mat(bsumlow4')};
    bsumhigh = {cell2mat(bsumhigh1'), cell2mat(bsumhigh2'), cell2mat(bsumhigh3'), cell2mat(bsumhigh4')};
    bsummid = {cell2mat(bsummid1'), cell2mat(bsummid2'), cell2mat(bsummid3'), cell2mat(bsummid4')};
    bsumpower = {cell2mat(bsumpower0'),cell2mat(bsumpower1'),cell2mat(bsumpower2'),cell2mat(bsumpower3'),cell2mat(bsumpower4')};
    bcwtpower = {cell2mat(bcwtpower0'),cell2mat(bcwtpower1'),cell2mat(bcwtpower2'),cell2mat(bcwtpower3'),cell2mat(bcwtpower4')};
    
    %Entropy
    ApproxEntropy = {cell2mat(AppEnt_normrec0'),cell2mat(AppEnt_normrec1'),cell2mat(AppEnt_normrec2'),cell2mat(AppEnt_normrec3'),...
        cell2mat(AppEnt_normrec4'),cell2mat(AppEnt_normrec5')};
    if minfeat ==0
        SampleEntropy = {cell2mat(SampEnt_normrec0'),cell2mat(SampEnt_normrec1'), cell2mat(SampEnt_normrec2'),...
        cell2mat(SampEnt_normrec3'),cell2mat(SampEnt_normrec4'),cell2mat(SampEnt_normrec5')};
    else
        SampleEntropy = {cell2mat(SampEnt_normrec0')};
    end
    WaveEntropy = {cell2mat(WEnt_normrec0'),cell2mat(WEnt_normrec1'),cell2mat(WEnt_normrec2'),cell2mat(WEnt_normrec3'),...
        cell2mat(WEnt_normrec4'),cell2mat(WEnt_normrec5')};

end

