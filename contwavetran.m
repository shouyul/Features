function [difflow1, difflow2, difflow3, difflow4, difflow5 ,...
    diffhigh1, diffhigh2, diffhigh3, diffhigh4, diffhigh5 ,...
    diffmid1, diffmid2, diffmid3, diffmid4, diffmid5 ,...
    sumlow1, sumlow2, sumlow3, sumlow4, sumlow5 ,...
    sumhigh1, sumhigh2, sumhigh3, sumhigh4, sumhigh5 ,...
    summid1, summid2, summid3, summid4, summid5 ,...
    sumpower0, sumpower1, sumpower2, sumpower3 ,sumpower4,sumpower5, ...
    cwtpower0, cwtpower1, cwtpower2, cwtpower3, cwtpower4, cwtpower5,...
    tcwtpower0, tcwtpower1, tcwtpower2, tsumlow1, tsumlow2,tsumhigh1, ...
    tsumhigh2, tsummid1, tsummid2,tsumpower0, tsumpower1, tsumpower2, ...
    bcwtpower0, bcwtpower1,bcwtpower2,bcwtpower3,bcwtpower4, bsumlow1, ...
    bsumlow2, bsumlow3, bsumlow4,bsumhigh1, bsumhigh2, bsumhigh3, bsumhigh4, ...
    bsummid1, bsummid2, bsummid3, bsummid4,bsumpower0, bsumpower1, bsumpower2, ...
    bsumpower3 , bsumpower4, tcwtf, cwtf, bcwtf, ...
    AppEnt_normrec0, AppEnt_normrec1,AppEnt_normrec2,AppEnt_normrec3,AppEnt_normrec4,AppEnt_normrec5,...
    SampEnt_normrec0,SampEnt_normrec1,SampEnt_normrec2,SampEnt_normrec3,SampEnt_normrec4,SampEnt_normrec5,...
    WEnt_normrec0,WEnt_normrec1,WEnt_normrec2,WEnt_normrec3,WEnt_normrec4,WEnt_normrec5...
    ] = contwavetran(data, minfeat, currentsamp, tau, selectedFeatureGroup)

%Function to compute continuous wavelet tranform features and entropy
%features for each trial

%INPUTS
%data = trial length * number channel array (e.g. 1282*16)
%minfeat = 1 for imagery vs imagery (to minimize number computed features)
%0 for rest vs imagery (to compute all possible features)
%currentsamp = 256; %Sampling rate of data
%tau = 2; %downsample factor (for approximate and sample entropy to save time)
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
%AppEnt_normrec0  = approximate entropy calculated for each channel
%SampEnt_normrec0  = sample entropy calculated for each channel
%cwtpower  = absolute alpha cwt pseudopower each channel
%tcwtpower  = absolute theta cwt pseudopower each channel
%bcwtpower  = absolute beta cwt pseudopower each channel
%tcwtf = frequencies for theta
%cwtf = frequencies for alpha
%bcwtf = frequencies for beta

    if ~exist('selectedFeatureGroup', 'var')
        selectedFeatureGroup = [];
    end
    
    apEnFeatureGroups = [80:85];     % approximate entropy feature groups
    sampEnFeatureGroups = [86:91];     % approximate entropy feature groups

    
    [sumlow1, sumlow2, sumlow3, sumlow4, sumlow5 ,...
    sumhigh1, sumhigh2, sumhigh3, sumhigh4, sumhigh5 ,...
    summid1, summid2, summid3, summid4, summid5 ,...
    sumpower0, sumpower1, sumpower2, sumpower3 ,sumpower4,sumpower5, ...
    cwtpower0, cwtpower1, cwtpower2, cwtpower3, cwtpower4, cwtpower5,...
    tcwtpower0, tcwtpower1, tcwtpower2, tsumlow1, tsumlow2,tsumhigh1, ...
    tsumhigh2, tsummid1, tsummid2,tsumpower0, tsumpower1, tsumpower2, ...
    bcwtpower0, bcwtpower1,bcwtpower2,bcwtpower3,bcwtpower4, bsumlow1, ...
    bsumlow2, bsumlow3, bsumlow4,bsumhigh1, bsumhigh2, bsumhigh3, bsumhigh4, ...
    bsummid1, bsummid2, bsummid3, bsummid4,bsumpower0, bsumpower1, bsumpower2, ...
    bsumpower3 , bsumpower4, tcwtf, cwtf, bcwtf, ...
    AppEnt_normrec0, AppEnt_normrec1,AppEnt_normrec2,AppEnt_normrec3,AppEnt_normrec4,AppEnt_normrec5,...
    SampEnt_normrec0,SampEnt_normrec1,SampEnt_normrec2,SampEnt_normrec3,SampEnt_normrec4,SampEnt_normrec5,...
    WEnt_normrec0,WEnt_normrec1,WEnt_normrec2,WEnt_normrec3,WEnt_normrec4,WEnt_normrec5] = deal(0);

%this computes CWT for the data, then normalizes it by dividing it by the
%mean and takes only the specified (alpha theta or beta) frequency band
tznormCWT = zeros([9 size(data)]);
znormCWT = zeros([9 size(data)]);
bznormCWT = zeros([13 size(data)]);
normrecsignal0 = zeros([size(data)]);
normrecsignal1 = zeros([size(data)]);
normrecsignal2 = zeros([size(data)]);
normrecsignal3 = zeros([size(data)]);
normrecsignal4 = zeros([size(data)]);
normrecsignal5 = zeros([size(data)]);

for chan = 1:size(data,2)
    [tempCWT, f] = cwt(data(:,chan),currentsamp,'amor');
    temp = abs(tempCWT);
    clear tempCWT
    tempmean = mean(temp(:));
    tempnormCWT = temp/tempmean;
    clear temp
    tznormCWT(:,:,chan) = tempnormCWT(42:50,:,:);
    znormCWT(:,:,chan) = tempnormCWT(33:41,:,:);
    bznormCWT(:,:,chan) = tempnormCWT(21:33,:,:);
    normrecsignal0(:,chan) = icwt(tempnormCWT,f, [7.5 13.5]);
    normrecsignal1(:,chan) = icwt(tempnormCWT,f, [7.5 9.5]);
    normrecsignal2(:,chan) = icwt(tempnormCWT,f, [8.5 10.5]);
    normrecsignal3(:,chan) = icwt(tempnormCWT,f, [9.5 11.5]);
    normrecsignal4(:,chan) = icwt(tempnormCWT,f, [10.5 12.5]);
    normrecsignal5(:,chan) = icwt(tempnormCWT,f, [11.5 13.5]);
    clear tempnormCWT
end

if minfeat == 0
%Calculate ApEn
AppEnt_normrec0 = zeros([1 size(data,2)]);
AppEnt_normrec1 = zeros([1 size(data,2)]);
AppEnt_normrec2 = zeros([1 size(data,2)]);
AppEnt_normrec3 = zeros([1 size(data,2)]);
AppEnt_normrec4 = zeros([1 size(data,2)]);
AppEnt_normrec5 = zeros([1 size(data,2)]);
if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups)) || isempty(selectedFeatureGroup))
    for chan = 1:size(data,2)
        if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups(1))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal0(:,chan);
            stdev = std(normrectemp);
            AppEnt_normrec0(chan) = ApEn(2, 0.2*stdev, normrectemp,tau);
        end
        if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups(2))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal1(:,chan);
            stdev = std(normrectemp);
            AppEnt_normrec1(chan) = ApEn(2, 0.2*stdev, normrectemp,tau);
        end
        if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups(3))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal2(:,chan);
            stdev = std(normrectemp);
            AppEnt_normrec2(chan) = ApEn(2, 0.2*stdev, normrectemp,tau);
        end
        if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups(4))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal3(:,chan);
            stdev = std(normrectemp);
            AppEnt_normrec3(chan) = ApEn(2, 0.2*stdev, normrectemp,tau);
        end
        if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups(5))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal4(:,chan);
            stdev = std(normrectemp);
            AppEnt_normrec4(chan) = ApEn(2, 0.2*stdev, normrectemp,tau);
        end
        if  (any(ismember(selectedFeatureGroup,apEnFeatureGroups(6))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal5(:,chan);
            stdev = std(normrectemp);
            AppEnt_normrec5(chan) = ApEn(2, 0.2*stdev, normrectemp,tau);
        end
    end
end
end


%Calculate SampEn
if minfeat
    SampEnt_normrec0 = zeros([1 size(data,2)]);
else
SampEnt_normrec0 = zeros([1 size(data,2)]);
SampEnt_normrec1 = zeros([1 size(data,2)]);
SampEnt_normrec2 = zeros([1 size(data,2)]);
SampEnt_normrec3 = zeros([1 size(data,2)]);
SampEnt_normrec4 = zeros([1 size(data,2)]);
SampEnt_normrec5 = zeros([1 size(data,2)]);

    if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups)) || isempty(selectedFeatureGroup))
        for chan = 1:size(data,2)
            if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups(1))) || isempty(selectedFeatureGroup))
                normrectemp = normrecsignal0(:,chan);
                stdev = std(normrectemp);
                SampEnt_normrec0(chan) = SampEn(2, 0.2*stdev, normrectemp,tau);
            end

            if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups(2))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal1(:,chan);
            stdev = std(normrectemp);
            SampEnt_normrec1(chan) = SampEn(2, 0.2*stdev, normrectemp,tau);
            end

            if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups(3))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal2(:,chan);
            stdev = std(normrectemp);
            SampEnt_normrec2(chan) = SampEn(2, 0.2*stdev, normrectemp,tau);
            end

            if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups(4))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal3(:,chan);
            stdev = std(normrectemp);
            SampEnt_normrec3(chan) = SampEn(2, 0.2*stdev, normrectemp,tau);
            end

            if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups(5))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal4(:,chan);
            stdev = std(normrectemp);
            SampEnt_normrec4(chan) = SampEn(2, 0.2*stdev, normrectemp,tau);
            end

            if  (any(ismember(selectedFeatureGroup,sampEnFeatureGroups(6))) || isempty(selectedFeatureGroup))
            normrectemp = normrecsignal5(:,chan);
            stdev = std(normrectemp);
            SampEnt_normrec5(chan) = SampEn(2, 0.2*stdev, normrectemp,tau);
            end
        end
    end

%Calculate Wentropy
WEnt_normrec0 = zeros([1 size(data,2)]);
WEnt_normrec1 = zeros([1 size(data,2)]);
WEnt_normrec2 = zeros([1 size(data,2)]);
WEnt_normrec3 = zeros([1 size(data,2)]);
WEnt_normrec4 = zeros([1 size(data,2)]);
WEnt_normrec5 = zeros([1 size(data,2)]);

% for chan = 1:size(data,2)
%     normrectemp = normrecsignal0(:,:,chan);
%     WEnt_normrec0(chan) = wentropy(normrectemp,'shannon');
%     normrectemp = normrecsignal1(:,:,chan);
%     WEnt_normrec1(chan) = wentropy(normrectemp,'shannon');
%     normrectemp = normrecsignal2(:,:,chan);
%     WEnt_normrec2(chan) = wentropy(normrectemp,'shannon');
%     normrectemp = normrecsignal3(:,:,chan);
%     WEnt_normrec3(chan) = wentropy(normrectemp,'shannon');
%     normrectemp = normrecsignal4(:,:,chan);
%     WEnt_normrec4(chan) = wentropy(normrectemp,'shannon');
%     normrectemp = normrecsignal5(:,:,chan);
%     WEnt_normrec5(chan) = wentropy(normrectemp,'shannon');       
% end
% 
% clear normrecsignal0
% clear normrecsignal1
% clear normrecsignal2
% clear normrecsignal3
% clear normrecsignal4
% clear normrecsignal5
% clear normrecsignal
% clear normrectemp
end

tcwtf = f(42:50);
cwtf = f(33:41);
bcwtf = f(21:33);

%Compute theta power
tcwtpoweralpha = squeeze(sum(tznormCWT,2));
tcwtpower0 = sum(tcwtpoweralpha(:,:),1);
tcwtpower1 = sum(tcwtpoweralpha(3:9,:),1);
tcwtpower2 = sum(tcwtpoweralpha(1:6,:),1);

%Compute Alpha power
cwtpoweralpha = squeeze(sum(znormCWT,2));
cwtpower0 = sum(cwtpoweralpha(:,:),1);
cwtpower1 = sum(cwtpoweralpha(6:9,:),1);
cwtpower2 = sum(cwtpoweralpha(4:7,:),1);
cwtpower3 = sum(cwtpoweralpha(3:6,:),1);
cwtpower4 = sum(cwtpoweralpha(2:4,:),1);
cwtpower5 = sum(cwtpoweralpha(1:3,:),1);

%Computer Beta Power
bcwtpoweralpha = squeeze(sum(bznormCWT,2));
bcwtpower0 = sum(bcwtpoweralpha(:,:),1);
bcwtpower1 = sum(bcwtpoweralpha(9:13,:),1);
bcwtpower2 = sum(bcwtpoweralpha(5:9,:),1);
bcwtpower3 = sum(bcwtpoweralpha(3:5,:),1);
bcwtpower4 = sum(bcwtpoweralpha(1:3,:),1);

if size(data,2) > 1
    %this subtracts left and right channels - only relevant for imagery versus imagery
    c = nchoosek([1:size(data,2)],2);
    diffCWT = zeros([size(znormCWT,1) size(znormCWT,2) size(c,1)]);

        for i = 1:size(c,1)
            diffCWT(:,:,i) = znormCWT(:,:,c(i,1))-znormCWT(:,:,c(i,2));
        end

    %this adds left and right channels - only relevant for rest vs imagery
    c = nchoosek([1:size(data,2)],2);
    tsumCWT = zeros([size(tznormCWT,1) size(tznormCWT,2) size(c,1)]);
    sumCWT = zeros([size(znormCWT,1) size(znormCWT,2) size(c,1)]);
    bsumCWT = zeros([size(bznormCWT,1) size(bznormCWT,2) size(c,1)]);

        for i = 1:size(c,1)
            tsumCWT(:,:,i) = tznormCWT(:,:,c(i,1))+tznormCWT(:,:,c(i,2));
            sumCWT(:,:,i) = znormCWT(:,:,c(i,1))+znormCWT(:,:,c(i,2));
            bsumCWT(:,:,i) = bznormCWT(:,:,c(i,1))+bznormCWT(:,:,c(i,2));
        end
else
    diffCWT = zeros([size(znormCWT,1) size(znormCWT,2) 1]);
    tsumCWT = zeros([size(tznormCWT,1) size(tznormCWT,2) 1]);
    sumCWT = zeros([size(znormCWT,1) size(znormCWT,2) 1]);
    bsumCWT = zeros([size(bznormCWT,1) size(bznormCWT,2) 1]);
end
clear tznormCWT
clear znormCWT
clear bznormCWT

%need to do the sum power for the sumCWT for each theta subband
tsumtime = squeeze(sum(tsumCWT,2));
tsumalpha0 = tsumtime;
tsumalpha1 = tsumtime(3:9,:);
tsumalpha2 = tsumtime(1:6,:);

%need to do the sum power for the sumCWT for each alpha subband
sumtime = squeeze(sum(sumCWT,2));
sumalpha0 = sumtime;
sumalpha1 = sumtime(6:9,:);
sumalpha2 = sumtime(4:7,:);
sumalpha3 = sumtime(3:6,:);
sumalpha4 = sumtime(2:4,:);
sumalpha5 = sumtime(1:3,:);

%need to do the sum power for the sumCWT for each beta subband
bsumtime = squeeze(sum(bsumCWT,2));
bsumalpha0 = bsumtime;
bsumalpha1 = bsumtime(9:13,:);
bsumalpha2 = bsumtime(5:9,:);
bsumalpha3 = bsumtime(3:5,:);
bsumalpha4 = bsumtime(1:3,:);

%Create sumCWT for threshold detections for theta above and below thresholds
tsumCWTalpha1 = tsumCWT(3:9,:,:);
tsumCWTalpha2 = tsumCWT(1:6,:,:);

%Create sumCWT for threshold detections for alpha above and below thresholds
sumCWTalpha1 = sumCWT(6:9,:,:);
sumCWTalpha2 = sumCWT(4:7,:,:);
sumCWTalpha3 = sumCWT(3:6,:,:);
sumCWTalpha4 = sumCWT(2:4,:,:);
sumCWTalpha5 = sumCWT(1:3,:,:);

%Create sumCWT for threshold detections for beta above and below thresholds
bsumCWTalpha1 = bsumCWT(9:13,:,:);
bsumCWTalpha2 = bsumCWT(5:9,:,:);
bsumCWTalpha3 = bsumCWT(3:5,:,:);
bsumCWTalpha4 = bsumCWT(1:3,:,:);

% find above, below and between power using thresholds for theta
[tsumlow1, tsumhigh1, tsummid1] = thresholdmagforcwt(tsumCWTalpha1);
[tsumlow2, tsumhigh2, tsummid2] = thresholdmagforcwt(tsumCWTalpha2);

% find above, below and between power using thresholds for alpha
[sumlow1, sumhigh1, summid1] = thresholdmagforcwt(sumCWTalpha1);
[sumlow2, sumhigh2, summid2] = thresholdmagforcwt(sumCWTalpha2);
[sumlow3, sumhigh3, summid3] = thresholdmagforcwt(sumCWTalpha3);
[sumlow4, sumhigh4, summid4] = thresholdmagforcwt(sumCWTalpha4);
[sumlow5, sumhigh5, summid5] = thresholdmagforcwt(sumCWTalpha5);

% find above, below and between power using thresholds for beta
[bsumlow1, bsumhigh1, bsummid1] = thresholdmagforcwt(bsumCWTalpha1);
[bsumlow2, bsumhigh2, bsummid2] = thresholdmagforcwt(bsumCWTalpha2);
[bsumlow3, bsumhigh3, bsummid3] = thresholdmagforcwt(bsumCWTalpha3);
[bsumlow4, bsumhigh4, bsummid4] = thresholdmagforcwt(bsumCWTalpha4);

% find sum of theta power
tsumpower0 = sum(tsumalpha0,1);
tsumpower1 = sum(tsumalpha1,1);
tsumpower2 = sum(tsumalpha2,1);

% find sum of alpha power
sumpower0 = sum(sumalpha0,1);
sumpower1 = sum(sumalpha1,1);
sumpower2 = sum(sumalpha2,1);
sumpower3 = sum(sumalpha3,1);
sumpower4 = sum(sumalpha4,1);
sumpower5 = sum(sumalpha5,1);

% find sum of beta power
bsumpower0 = sum(bsumalpha0,1);
bsumpower1 = sum(bsumalpha1,1);
bsumpower2 = sum(bsumalpha2,1);
bsumpower3 = sum(bsumalpha3,1);
bsumpower4 = sum(bsumalpha4,1);

%need to do the above, below and mid for each alpha subband using threshold magforcwt
diffalpha1 = diffCWT(6:9,:,:);
diffalpha2 = diffCWT(4:7,:,:);
diffalpha3 = diffCWT(3:6,:,:);
diffalpha4 = diffCWT(2:4,:,:);
diffalpha5 = diffCWT(1:3,:,:);

[difflow1, diffhigh1, diffmid1] = thresholdmagforcwt(diffalpha1);
[difflow2, diffhigh2, diffmid2] = thresholdmagforcwt(diffalpha2);
[difflow3, diffhigh3, diffmid3] = thresholdmagforcwt(diffalpha3);
[difflow4, diffhigh4, diffmid4] = thresholdmagforcwt(diffalpha4);
[difflow5, diffhigh5, diffmid5] = thresholdmagforcwt(diffalpha5);

