function [cpsdmag0, cpsdmag1, cpsdmag2,cpsdmag3,cpsdmag4,cpsdmag5, cpsdphase1,cpsdphase2,cpsdphase3,cpsdphase4,cpsdphase5,...
    anglephase1,anglephase2,anglephase3,anglephase4,anglephase5,...
    cpsdmagt0, cpsdmagt1, cpsdmagt2, cpsdmagb0,cpsdmagb1,cpsdmagb2, cpsdmagb3,cpsdmagb4, tcpsdf, cpsdf, bcpsdf] = ...
    crosspower(data, minfeat, currentsamp, halfsamp, samp4overlap, thetarange,alpharange,betarange)

%calculates cross-power spectral density features (magnitude and phase delay) for alpha, theta and beta
%frequencies and their subbands 

%INPUTS
% data = data points*array of channels (1282*16)
% minfeat = 1 for minimal number features (Imagery vs Imagery) or 0 for all
%   features (Rest for Imagery - note computationally expensive)
% currentsamp = sampling rate in Hz (e.g. 256 Hz)
% halfsamp = half the ampling rate in Hz (e.g. 128 Hz)
% samp4overlap = overlap of windows NEEDS TO BE INTEGER [ideally round(0.7813*halfsamp)]
% thetarange = vector specifying desired theta frequency range
%alpharange = vector specifying desired alpha frequency range
%betarange = vector specifying desired beta frequency range

%OUTPUTS (for entire specified ranges, identified as 0; as well as subbands, identified as 1-5)
% anglephase = Cross Power Spectral Density phase (imag) is the phase delay
% of cross power spectral density between channel pairs
% cpsdmag = Cross Power Spectral Density magnitude (abs)is the magnitude of
% the cross power spectral density between channel pairs
% cpsdmagt = CPSD for theta
% cpsdmagb = CPSD for beta

if isempty(thetarange)
    thetarange = [4:0.5:7];
end
if isempty(alpharange)
    alpharange = [7.5:0.5:13.5];
end
if isempty(betarange)
    betarange = [13:0.5:30];
end

%this computes CPSD for the data for specific channel pairs
[cpsdmag0, cpsdmag1, cpsdmag2,cpsdmag3,cpsdmag4,cpsdmag5, cpsdphase1,cpsdphase2,cpsdphase3,cpsdphase4,cpsdphase5,...
    anglephase1,anglephase2,anglephase3,anglephase4,anglephase5,...
    cpsdmagt0, cpsdmagt1, cpsdmagt2, cpsdmagb0,cpsdmagb1,cpsdmagb2, cpsdmagb3,cpsdmagb4, tcpsdf, cpsdf, bcpsdf] = deal(0);
if size(data,2) == 16
[A,B] = meshgrid([1:8],[9:16]); % identifies left (chans 1-8) and right (chans 9-16) channel pairs
d=cat(2,A',B');
c=cat(1,(reshape(d,[],2)),cat(2,5*ones(7,1),[1:4,6:8]'),cat(2,11*ones(7,1),[9:10,12:16]'));
else
    c = [1 2];
end

% c = nchoosek([1:16],2);
tpxy = zeros([size(thetarange,2) size(c,1)]);
pxy = zeros([size(alpharange,2) size(c,1)]);
bpxy = zeros([size(betarange,2) size(c,1)]);

    for i = 1:size(c,1)
        if minfeat == 0
             [tpxy(:,i), tf] = cpsd(data(:,c(i,1)),data(:, c(i,2)),halfsamp,samp4overlap,thetarange,currentsamp);
             [bpxy(:,i), bf] = cpsd(data(:,c(i,1)),data(:, c(i,2)),halfsamp,samp4overlap,betarange,currentsamp);
        end
    [pxy(:,i), f] = cpsd(data(:,c(i,1)),data(:, c(i,2)),halfsamp,samp4overlap,alpharange,currentsamp);
    end

cpsdf = f;
if minfeat == 0
    tcpsdf = tf;
    bcpsdf = bf;
end

%Take the magnitude of the cross power spectral density
realmag = abs(pxy); %This gives the magnitude, not the real!!!!
%     realmag = real(pxy);

if minfeat == 0
    realmagt = real(tpxy);
    clear tpxy
    realmagb = real(bpxy);
    clear bpxy
end
%Take the phase of the cross power spectral density
% imagphase = imag(pxy);
anglephase = angle(pxy);
clear pxy

%magnitude Alpha
cpsdmag0 = mean(realmag,1);
cpsdmag1 = mean(realmag(1:5,:),1);
cpsdmag2 = mean(realmag(3:7,:),1);
cpsdmag3 = mean(realmag(5:9,:),1);
cpsdmag4 = mean(realmag(7:11,:),1);
cpsdmag5 = mean(realmag(9:13,:),1);
clear realmag

if minfeat == 0
    %magnitude theta
    cpsdmagt0 = mean(realmagt,1);
    cpsdmagt1 = mean(realmagt(1:5,:),1);
    cpsdmagt2 = mean(realmagt(3:7,:),1);
    clear realmagt

    %magnitude beta
    cpsdmagb0 = mean(realmagb,1);
    cpsdmagb1 = mean(realmagb(1:11,:),1);
    cpsdmagb2 = mean(realmagb(11:21,:),1);
    cpsdmagb3 = mean(realmagb(21:27,:),1);
    cpsdmagb4 = mean(realmagb(27:35,:),1);
    clear realmagb
end

%phase Alpha
% cpsdphase1 = mean(imagphase(1:5,:),1);
% cpsdphase2 = mean(imagphase(3:7,:),1);
% cpsdphase3 = mean(imagphase(5:9,:),1);
% cpsdphase4 = mean(imagphase(7:11,:),1);
% cpsdphase5 = mean(imagphase(9:13,:),1);

%anglephase Alpha
anglephase1 = mean(anglephase(1:5,:),1);
anglephase2 = mean(anglephase(3:7,:),1);
anglephase3 = mean(anglephase(5:9,:),1);
anglephase4 = mean(anglephase(7:11,:),1);
anglephase5 = mean(anglephase(9:13,:),1);

