function [segmentpass0, segmentpass1, segmentpass2, segmentpass3, segmentpass4, segmentpass5, f] = ...
    newalphaspindle(data, samplingrate, windowlength, overlap, freqrange, maxrange, numstepinHz, SNRthreshold)
%Simon et al 2011 conditions for identifying alpha spindles in each "segment" of a spectrogram
% Simon, M. et al. Eeg alpha spindle measures as indicators of driver fatigue under real traffic conditions. 
%Clin. Neurophysiol. 122, 1168–1178 (2011).

%INPUTS:
% data = vector of data points (one trial i.e. 1*1282)
% windowlength = length of the window used for the Short-time fourier
    %transform (typically half the size of the sampling rate e.g. 128 Hz)
% overlap = overlap of windows NEEDS TO BE INTEGER [ideally round(0.7813*halfsamp)]
% samplingrate = sampling rate in Hz (e.g. 256 Hz)
% freqrange = frequency range you want to use for the spectrogram (default  = [1:0.5:40] i.e. 1-40 Hz)
% maxrange = range of alpha band of interest (can use [7.5:13.5])
% numstepinHz = number of steps per Hz for your spectrogram (so 2 for [0.5:0.5:40] and 1 for [1:1:40])
% SNRthreshold = how many times the noise does the signal have to be to be considered significant (default = 2)

%OUTPUTS: 

%segmentpass = which segments of spectrogram passed all three conditions in
    %i.e. which segments are classified as alpha spindles (vector of ones and/or zeros) 
    %NOTE these are divided by identifying whether the "maximum frequency" was in specific ranges:
    % 0 = entire frequency range specified (7.5-13.5 Hz)    
    % 1 = first subband of frequency range specified (7.5-9.5 Hz)    
    % 2 = second subband of frequency range specified (8.5-10.5 Hz)
    % 3 = third subband of frequency range specified (9.5-11.5 Hz)
    % 4 = fourth subband of frequency range specified (10.5-12.5 Hz)
    % 5 = fifth subband of frequency range specified (11.5-13.5 Hz)
% f = frequencies corresponding to each step in the spectrogram

%Additional variables that are computed and could be outputted if desired
%SNR = signal to noise ratio; i.e. the oscillation index of each segment where the SNR is the signal
    %of the FWHM over the noise for that frequency range
% freqpass = whether the segment passed the frequency condition (ones and/or zeros)
% Hammbandpass = whetehr the segment passed the Hamming band window noise  condition (ones and/or zeros)
% SNRpass = whether the segment passed the SNR condition (signal >= 2* noise of that segment)
% freqmax = maximum frequency of each segment
% p = power from spectrogram of each segment across the specified frequency range
% range = the FWHM range around the maximum of each segment


%keywords: alpha, spindle, alpha spindle

if isempty(samplingrate)
    samplingrate = 256;
end

if isempty(windowlength)
    windowlength = samplingrate;
end

if isempty(overlap)
    overlap = 0.75*samplingrate;
end

if isempty(freqrange) 
    freqrange = [1:0.5:40];
end

if isempty(maxrange)
    maxrange = [7:13];
end

if isempty(numstepinHz)
    numstepinHz = 2;
end

if isempty(SNRthreshold)
    SNRthreshold = 2;
end

%Create a spectrogram using short time fourier transform
[~,f,~,p] = spectrogram(data,windowlength,overlap,freqrange, samplingrate,'yaxis');

%Find the maximum power and the index of the maximum for each time segment
[maxp,idxmaxminusfive] = max(p(3*numstepinHz:end,:)); %indicates where in the column to start (because you make the entire spectrogram, but only need the max of 3-40 Hz)
idxmax = idxmaxminusfive+(3*numstepinHz-1*numstepinHz+1);

%Find the frequency of the maximum power (i.e. peak)
freqmax = f(idxmax);

%Find which segment maxima are within the entire alpha frequency
freqpass0 = (freqmax>= maxrange(1) & freqmax<= maxrange(7))';
%Find which segment maxima are within the first alpha subband frequency
freqpass1 = (freqmax>= maxrange(1) & freqmax<= maxrange(3))';
%Find which segment maxima are within the second alpha subband frequency
freqpass2 = (freqmax>= maxrange(2) & freqmax<= maxrange(4))';
%Find which segment maxima are within the first alpha subband frequency
freqpass3 = (freqmax>= maxrange(3) & freqmax<= maxrange(5))';
%Find which segment maxima are within the first alpha subband frequency
freqpass4 = (freqmax>= maxrange(4) & freqmax<= maxrange(6))';
%Find which segment maxima are within the first alpha subband frequency
freqpass5 = (freqmax>= maxrange(5) & freqmax<= maxrange(7))';

%Find the frequency bands that are above the halfmax
maxperwindow = (repmat(maxp,[size(p,1) 1]))*0.5;
biggerthanhalfmax = p>=maxperwindow;

%Find the size of the band around the maximum that are above halfmax
FWHMrange = zeros(2, size(idxmax,2));
sizebigband = zeros(size(idxmax));
for window = 1:size(idxmax,2)
    firstZero = idxmax(window) + 1 - find(~biggerthanhalfmax(idxmax(window):-1:1, window),1);%this moves from maximum to beginnning
    nextZero = find(~biggerthanhalfmax(idxmax(window):end, window), 1)+idxmax(window)-1;%this moves from maximum to end
        if isempty(firstZero)
            firstZero = 0;
        end
        if isempty(nextZero)
            nextZero = size(biggerthanhalfmax, 1) + 1;
        end
    FWHMrange(:,window) = [firstZero+1 nextZero-1];
    sizebigband(window) = nextZero - firstZero -1; 
end

%Check to see if size of freqband is smaller than twice the noise bandwidth of the Hamming window
Hammbandpass = sizebigband<=numstepinHz*2*(enbw(hamming(windowlength),samplingrate));

%plot exponential curve to fit 1/f curve of noise across average of all segments
avgp = mean(p,2)';
ExpNoiseFit = fit(f,avgp','exp1');
ExpNoiseCurve = feval(ExpNoiseFit, f);
areaundermean = trapz(f,avgp);

%Determine area under power spectrum for each segment, find ratio to mean
%and multiply that ratio with the Noise Curve
areaofpowerspect = trapz(f,p);
ratio = areaofpowerspect/areaundermean;
SegmentCurve = (ratio.*ExpNoiseCurve)';

%Check to see if signal of alpha (according to FWHM) is twice as big as the noise
areaunderalphasignal = zeros(size(idxmax));
areaunderalphanoise = zeros(size(idxmax));

for segment = 1:size(idxmax,2)
        if sizebigband(segment) == 1
            areaunderalphasignal(segment) = p([FWHMrange(1,segment):FWHMrange(2,segment)],segment);
            areaunderalphanoise(segment) = SegmentCurve(segment,[FWHMrange(1,segment):FWHMrange(2,segment)]);
        else
            areaunderalphasignal(segment) = trapz(f([FWHMrange(1,segment):FWHMrange(2,segment)]),p([FWHMrange(1,segment):FWHMrange(2,segment)],segment));
            areaunderalphanoise(segment) = trapz(f([FWHMrange(1,segment):FWHMrange(2,segment)]),SegmentCurve(segment,[FWHMrange(1,segment):FWHMrange(2,segment)]));
        end
end
% figure;
% plot(f,avgp)
% hold on; plot(f,ExpNoiseCurve)
SNR = areaunderalphasignal./areaunderalphanoise;
SNRpass = SNR>=SNRthreshold;

%Gives logical array if segment passes all conditions
segmentpass0 = freqpass0 & Hammbandpass & SNRpass;
segmentpass1 = freqpass1 & Hammbandpass & SNRpass;
segmentpass2 = freqpass2 & Hammbandpass & SNRpass;
segmentpass3 = freqpass3 & Hammbandpass & SNRpass;
segmentpass4 = freqpass4 & Hammbandpass & SNRpass;
segmentpass5 = freqpass5 & Hammbandpass & SNRpass;