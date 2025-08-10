function [mscohmag0, mscohmag1, mscohmag2,mscohmag3,mscohmag4,mscohmag5, mscohf] = ...
    magsquarecoherence(data,currentsamp, halfsamp, samp4overlap)

%this computes magnitude squared coherence for the data

%INPUTS
% data = data points*array of channels (1282*16)
% currentsamp = sampling rate in Hz (e.g. 256 Hz)
% halfsamp = half the ampling rate in Hz (e.g. 128 Hz)
% samp4overlap = overlap of windows NEEDS TO BE INTEGER [ideally round(0.7813*halfsamp)]

%OUTPUTS
%mscohmag = magnitude squared coherence for entire alpha range (0) and
%subbands (1-5)
%mscohf = freuqencies computed from coherence

% [A,B] = meshgrid([1:8],[9:16]);
% d=cat(2,A',B');
% c=cat(1,(reshape(d,[],2)),cat(2,5*ones(7,1),[1:4,6:8]'),cat(2,11*ones(7,1),[9:10,12:16]'));
c = nchoosek([1:size(data,2)],2); % calculated for all possible channel combinations
mxy = zeros([13 size(c,1)]);
    for i = 1:size(c,1)
    [mxy(:,i), f] = mscohere(data(:,c(i,1)),data(:, c(i,2)),halfsamp,samp4overlap,[7.5:0.5:13.5],currentsamp);
    end
    
mscohf = f;

%magnitude Alpha
mscohmag0 = mean(mxy,1);
mscohmag1 = mean(mxy(1:5,:),1);
mscohmag2 = mean(mxy(3:7,:),1);
mscohmag3 = mean(mxy(5:9,:),1);
mscohmag4 = mean(mxy(7:11,:),1);
mscohmag5 = mean(mxy(9:13,:),1);
clear realmag
