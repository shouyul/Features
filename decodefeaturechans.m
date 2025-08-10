function [ featureGroup, channelPairing ] = decodefeaturechans( selectedFeatures, problem )
%Function to determine which channels are used for specific feature
%computations, based on final feature number (e.g. 4857 = cpsdchans 1, 14)

%Inputs: 
%selectedFeatures = vector of final feature number
%problem = Imagery or Rest 

%Chanchoose2
chanchoose2 = nchoosek([1:16],2);

%LchanvsRchan
[A,B] = meshgrid([1:8],[9:16]);
d=cat(2,A',B');
LchanvsRchan=cat(1,(reshape(d,[],2)),cat(2,5*ones(7,1),[1:4,6:8]'),cat(2,11*ones(7,1),[9:10,12:16]'));

switch problem
    case 'Imagery'
        % each line = [ (# channel combinations) (# of feature groups with those combinations)]
        featureGroupings = [120 20; 78 5]; 
    otherwise % case 'Rest'
        % each line = [ (# channel combinations) (# of feature groups with those combinations)]
        featureGroupings = [120 54; 78 19; 16 32]; 
end

featureGroup = zeros(length(selectedFeatures),1);
channelPairing = zeros(length(selectedFeatures),2);

for featureIndex = 1:length(selectedFeatures)
    thisFeature = selectedFeatures(featureIndex);
    if thisFeature <= featureGroupings(1,1)*featureGroupings(1,2)
        %one of the 120-feature groups
        featureGroup(featureIndex) = ceil(thisFeature / featureGroupings(1,1));
        channelPairing(featureIndex, :) = chanchoose2( mod(thisFeature-1, featureGroupings(1,1))+1, :);
        
    elseif thisFeature <=featureGroupings(1,1)*featureGroupings(1,2)+featureGroupings(2,1)*featureGroupings(2,2)
        %one of ths 78 feature groups
        thisFeature = thisFeature - featureGroupings(1,1)*featureGroupings(1,2);
        featureGroup(featureIndex) = ceil(thisFeature / featureGroupings(2,1));
        channelPairing(featureIndex, :) = LchanvsRchan( mod(thisFeature-1, featureGroupings(2,1))+1, :);
        featureGroup(featureIndex) = featureGroup(featureIndex) + featureGroupings(1,2);
    else % 16
        %one of the 16 feature groups
        thisFeature = thisFeature - (featureGroupings(1,1)*featureGroupings(1,2)+featureGroupings(2,1)*featureGroupings(2,2));
        featureGroup(featureIndex) = ceil(thisFeature / featureGroupings(3,1));
        channelPairing(featureIndex, :) = [mod(thisFeature-1, featureGroupings(3,1))+1 mod(thisFeature-1, featureGroupings(3,1))+1];
        featureGroup(featureIndex) = featureGroup(featureIndex) + featureGroupings(1,2) + featureGroupings(2,2);
    end
end

return;

%IF Rest
% Features = {sumspind, ... %120 (= 120)
%     sumlow, sumhigh, summid, ... %120 each *5 alpha subbands(1,2,3,4,5) *3 (=1800)
%     sumpower, ...%120 *6 alpha subbands (0,1,2,3,4,5)( = 720)
%     tsumlow,tsumhigh, tsummid,... %120 each *2 theta subbands(1,2) * 3 (=720)
%     tsumpower,... %120 * 3 theta subbands (0,1,2) (=360)
%     bsumlow, bsumhigh, bsummid, ... % 120 each*4 beta subbands(1,2,3,4) * 3 (=1440)
%     bsumpower,... % 120 * 5 beta subbands(0,1,2,3,4) (=600)
%     mscoh, ... %120 each *6 alpha subbands(0,1,2,3,4,5)( = 720)

%     anglephase, ... %78 each * 5 alpha subbbands(1,2,3,4,5) (=390)
%     cpsdmag, ... %78 each*6 alpha subbands(0,1,2,3,4,5) (=468)
%     cpsdmagt, ...%78 each *3 theta subbands(0,1,2) (=234)
%     cpsdmagb,... %78 each *5 beta subbands(0,1,2,3,4) (=390)

%     numspindperchanall,ApproxEntropy,SampleEntropy, cwtpower, ... % 16 each *6 alpha subbands(0,1,2,3,4,5) *4 (=384)
%     tcwtpower,... %16 *3 theta subbands (0,1,2)( =48)
%     bcwtpower,... %16 *5 beta subbands (0,1,2,3,4)(=80)

%IF Imagery
% Features = {diffspind, difflow, diffhigh, diffmid ... %120 *5 alpha subbands (1,2,3,4,5) *4 (= 2400)
%     anglephase... %78 *5 alpha subbands (1,2,3,4,5) (=390)



