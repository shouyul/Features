function [Features, FeatureNames, FeatureMat] = ComputeImageryFeatures(data, problem, minfeat, ...
    channelPairing,featureGroup, selectedFeatures)

%Function computes Imagery features based options (e.g. sampling rate) and input data (X).

%INPUTS:
%data = 3D matrix of channels*datapoints in each trial*number of trials (e.g. 16x1282x25)
%problem = string input: 'Imagery' (Computes features ideal for imagery vs imagery classification problem)
    %or "Rest" (computes features ideal for imagery vs rest problem)
%minfeat = 0 for all possible features to compute (not recommended for
    %imagery problem) or 1 for minimizing number of features (recommended for imagery problem)
%channelPairing = can specify vector of channel pairs to compute features
    %for, otherwise don't include as input argument if you want features  computed for all
    %channels NOTE: can compute this with function decodefeaturechans.m
%featureGroup = specify which feature groups you want computed (see below),
    %otherwise don't include in input NOTE: can compute this with function decodefeaturechans.m
%selectedFeatures = specific feature IDs if you know exactly which features
    %you need computed (in order to minimize computation time) otherwise leave empty

%OUTPUTS:
%Features = 1*#features cell with each cell containing groups of features
%FeatureNames = 1*#features cell with each cell containing name of featuregroups
%FeatureMat = matrix version of Features output

%Features include:
% for 78 comparisons of left vs right channels
% for 5 alpha sub-bands between 7-13, depending on the features
% for sum/mean for all but alpha spindles

%NOTE: the par structure (mostly used in alpha spindle computation) can be edited, although the below is default values
%Preparing parameter structure
par.currentsamp = 256; %Sampling rate of data
par.halfsamp = par.currentsamp/2; %Size of sliding window
par.samp4overlap = round(0.7813*par.halfsamp); %par.samp4overlap NEEDS TO BE AN INTEGER!!! %overlap of window
par.freqrange = [1:0.5:40]; %frequency range of interest
par.numstepperHz = 2; %how many steps per Hz in frequency range
par.SNRthreshold = 2; %minimum signal to noise ratio for alpha spindle condition
par.tau = 2; %downsample factor (for approximate and sample entropy to save time)

if nargin<4 %If there is no channelPairing variable, then you want to compute all the features, so make channelPairing blank
    channelPairing = [];
    %     switch problem
    %         case 'Imagery'
    %         featureGroup = [1:25];
    %         otherwise %Rest
    %         featureGroup = [1:105];
    %     end
end
if ~exist('featureGroup', 'var') %If this is not an input, i.e. you want to compute all features, 
    %create a variable featureGroup depending on the problem
    switch problem
        case 'Imagery'
            featureGroup = [1:25]';
        otherwise %Rest
            featureGroup = [1:105]';
    end
end

% reshape the incoming data
perdata = permute(data, [2 1 3]);
[m, n, t] = size(perdata);
data = squeeze(mat2cell(perdata,m,n,ones(1,t)))';

switch problem
    case 'Imagery'
        %Alpha Spindle [1:5]
        if isempty(channelPairing) %Compute all features
            [~,diffspind,~] = alphaspindlefeatures(data, par);
        else
            alphaspindlemember = ismember(featureGroup,[1:5]);
            alphaspindlefeatureGroup = featureGroup(alphaspindlemember);
            if any(alphaspindlemember)
                alphaspindlechanpairs = channelPairing(alphaspindlemember,:);
                alphaFeatureMatIdx = selectedFeatures(alphaspindlemember);
                %ImageryFeatureNumbers - {sumspind [NaN], diffspind [1:5], numspindperchanall(0-5) [NaN]}
                for chanpairidx = 1:size(alphaspindlechanpairs,1)
                    [~,diffspind,~] = alphaspindlefeatures(data, par,alphaspindlechanpairs(chanpairidx,:)); %cell = channel/col. = featureGroup
                    FeatureMatAll{alphaFeatureMatIdx(chanpairidx)} = diffspind{alphaspindlefeatureGroup(chanpairidx)}; %concatenating all the calculated channel pairs, selecting the specific featureGroups
                end
                clear diffspindtemp;
            end
        end
        
        %Continuous Wavelet Transform [6:20]
        if isempty(channelPairing) %Compute all features
            [difflow,diffhigh,diffmid, ~, ~,~,~,~,~,~,~, ~,~,~,~, ~, ~,~,~,~] = waveletfeatures(data,par, minfeat);
        else
            cwtmember = ismember(featureGroup,[6:20]);
            cwtfeatureGroup = featureGroup(cwtmember);
            if any(cwtmember)
                cwtchanpairs = channelPairing(cwtmember,:);
                cwtFeatureMatIdx = selectedFeatures(cwtmember);
                %ImageryFeatureNumbers- {difflow(a1-5) [6:10],diffhigh(a1-5)[11:15],diffmid(a1-5)[16:20], ...
                %                     sumlow(a1-5)[NaN], sumhigh(a1-5)[NaN],summid(a1-5)[NaN],sumpower(a0-5)[NaN],...
                %                     tsumlow(t1-2)[NaN],tsumhigh(t1-2)[NaN],tsummid(t1-2)[NaN],tsumpower(t0-2)[NaN],...
                %                     bsumlow(b1-4)[NaN],bsumhigh(b1-4)[NaN],bsummid(b1-4)[NaN],bsumpower(b0-4)[NaN], ...
                %                     ApproxEntropy(a0-5)[NaN], SampleEntropy(a0-5)[NaN],...
                %                     cwtpower(a0-5)[NaN], tcwtpower(t0-2)[NaN],bcwtpower(b0-4)[NaN]}
                for chanpairidx = 1:size(cwtchanpairs,1)
                    [difflow,diffhigh,diffmid, ~, ~,~,~,~,~,~,~, ~,~,~,~, ~, ~,~,~,~] = waveletfeatures(data,par, minfeat, cwtchanpairs(chanpairidx,:));
                    switch cwtfeatureGroup(chanpairidx)
                        case num2cell(6:10)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = difflow{cwtfeatureGroup(chanpairidx)-5};
                        case num2cell(11:15)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = diffhigh{cwtfeatureGroup(chanpairidx)-10};
                        case num2cell(16:20)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = diffmid{cwtfeatureGroup(chanpairidx)-15};
                    end
                end
            end
        end
        
        %CPSD cross-power spectral density features [21:25] 
        if isempty(channelPairing) %Compute all features
            [ anglephase, ~,  ~, ~,]  = cpsdfeatures( data, par, minfeat );
        else
            cpsdmember = ismember(featureGroup,[21:25]);
            cpsdfeatureGroup = featureGroup(cpsdmember);
            if any(cpsdmember)
                cpsdchanpairs = channelPairing(cpsdmember,:);
                cpsdFeatureMatIdx = selectedFeatures(cpsdmember);
                %ImageryFeatureNumbers- {anglephase(a1-5) [21:25], cpsdmag(a0-5) [NaN],   cpsdmagt(t0-2)[NaN],   cpsdmagb(t0-2)[NaN]}
                for chanpairidx = 1:size(cpsdchanpairs,1)
                    [ anglephase, ~,  ~, ~,]  = cpsdfeatures( data, par, minfeat, cpsdchanpairs(chanpairidx,:));
                    FeatureMatAll{cpsdFeatureMatIdx(chanpairidx)} = anglephase{cpsdfeatureGroup(chanpairidx)-20};
                end
            end
        end
        if isempty(channelPairing)
            Features = {diffspind, difflow, diffhigh, diffmid, anglephase};
        else
            Features = [];
        end
        
    otherwise %Rest
        
        %Alpha Spindle [1 74:79]
        if isempty(channelPairing) %Compute all features
            [sumspind,~,numspindperchanall] = alphaspindlefeatures(data, par);
        else
            alphaspindlemember = ismember(featureGroup,[1 74:79]);
            alphaspindlefeatureGroup = featureGroup(alphaspindlemember);
            if any(alphaspindlemember)
                alphaspindlechanpairs = channelPairing(alphaspindlemember,:);
                alphaFeatureMatIdx = selectedFeatures(alphaspindlemember);
                %RestFeatureNumbers -    {sumspind [1],   diffspind [NaN], numspindperchanall(0-5) [74:79]}
                for chanpairidx = 1:size(alphaspindlechanpairs,1)
                    [sumspind,~,numspindperchanall] = alphaspindlefeatures(data, par,alphaspindlechanpairs(chanpairidx,:)); %cell = channel/col. = featureGroup
                    switch alphaspindlefeatureGroup(chanpairidx)
                        case num2cell(1)
                            FeatureMatAll{alphaFeatureMatIdx(chanpairidx)} = sumspind{alphaspindlefeatureGroup(chanpairidx)};
                        case num2cell(74:79)
                            FeatureMatAll{alphaFeatureMatIdx(chanpairidx)} = numspindperchanall{alphaspindlefeatureGroup(chanpairidx)-73};
                    end
                end
            end
        end
        
        %Continuous Wavelet Transform [2:48 80:105]
        if isempty(channelPairing) %Compute all features
            [~,~,~, sumlow, sumhigh,summid,sumpower,tsumlow,tsumhigh,tsummid,tsumpower,...
                bsumlow,bsumhigh,bsummid,bsumpower, ApproxEntropy, SampleEntropy,cwtpower,tcwtpower,bcwtpower] = ...
                waveletfeatures(data,par, minfeat);
        else
            cwtmember = ismember(featureGroup,[2:48 80:105]);
            cwtfeatureGroup = featureGroup(cwtmember);
            if any(cwtmember)
                cwtchanpairs = channelPairing(cwtmember,:);
                cwtFeatureMatIdx = selectedFeatures(cwtmember);
                %RestFeatureNumbers- {difflow(a1-5) [NaN],diffhigh(a1-5)[NaN],diffmid(a1-5)[NaN], ...
                %                     sumlow(a1-5)[2:6], sumhigh(a1-5)[7:11],summid(a1-5)[12:16],sumpower(a0-5)[17:22],...
                %                     tsumlow(t1-2)[23:24],tsumhigh(t1-2)[25:26],tsummid(t1-2)[27:28],tsumpower(t0-2)[29:31],...
                %                     bsumlow(b1-4)[32:35],bsumhigh(b1-4)[36:39],bsummid(b1-4)[40:43],bsumpower(b0-4)[44:48], ...
                %                     ApproxEntropy(a0-5)[80:85],SampleEntropy(a0-5)[86:91],...
                %                     cwtpower(a0-5)[92:97],tcwtpower(t0-2)[98:100],bcwtpower(b0-4)[101:105]}
                for chanpairidx = 1:size(cwtchanpairs,1)
                    [~,~,~, sumlow, sumhigh,summid,sumpower,tsumlow,tsumhigh,tsummid,tsumpower,...
                        bsumlow,bsumhigh,bsummid,bsumpower, ApproxEntropy, SampleEntropy,cwtpower,tcwtpower,bcwtpower] = ...
                        waveletfeatures(data,par, minfeat, cwtchanpairs(chanpairidx,:), cwtfeatureGroup(chanpairidx));
                    switch cwtfeatureGroup(chanpairidx)
                        case num2cell(2:6)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = sumlow{cwtfeatureGroup(chanpairidx)-1};
                        case num2cell(7:11)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = sumhigh{cwtfeatureGroup(chanpairidx)-6};
                        case num2cell(12:16)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = summid{cwtfeatureGroup(chanpairidx)-11};
                        case num2cell(17:22)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = sumpower{cwtfeatureGroup(chanpairidx)-16};
                        case num2cell(23:24)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = tsumlow{cwtfeatureGroup(chanpairidx)-22};
                        case num2cell(25:26)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = tsumhigh{cwtfeatureGroup(chanpairidx)-24};
                        case num2cell(27:28)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = tsummid{cwtfeatureGroup(chanpairidx)-26};
                        case num2cell(29:31)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = tsumpower{cwtfeatureGroup(chanpairidx)-28};
                        case num2cell(32:35)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = bsumlow{cwtfeatureGroup(chanpairidx)-31};
                        case num2cell(36:39)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = bsumhigh{cwtfeatureGroup(chanpairidx)-35};
                        case num2cell(40:43)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = bsummid{cwtfeatureGroup(chanpairidx)-39};
                        case num2cell(44:48)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = bsumpower{cwtfeatureGroup(chanpairidx)-43};
                        case num2cell(80:85)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = ApproxEntropy{cwtfeatureGroup(chanpairidx)-79};
                        case num2cell(86:91)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = SampleEntropy{cwtfeatureGroup(chanpairidx)-85};
                        case num2cell(92:97)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = cwtpower{cwtfeatureGroup(chanpairidx)-91};
                        case num2cell(98:100)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = tcwtpower{cwtfeatureGroup(chanpairidx)-97};
                        case num2cell(101:105)
                            FeatureMatAll{cwtFeatureMatIdx(chanpairidx)} = bcwtpower{cwtfeatureGroup(chanpairidx)-100};
                    end
                end
            end
        end
        
        %mscoherence- magnitude squared coherence features [49:54]
        if isempty(channelPairing) %Compute all features
            [ mscoh ] = mscoherefeatures(data, par);
        else
            mscohmember = ismember(featureGroup,[49:54]);
            mscohfeatureGroup = featureGroup(mscohmember);
            if any(mscohmember)
                mscohchanpairs = channelPairing(mscohmember,:);
                mscohFeatureMatIdx = selectedFeatures(mscohmember);
                %RestFeatureNumbers-    {mscoh(a0-5) [49:54]}
                for chanpairidx = 1:size(mscohchanpairs,1)
                    
                    [ mscoh ] = mscoherefeatures( data, par, mscohchanpairs(chanpairidx,:) );
                    FeatureMatAll{mscohFeatureMatIdx(chanpairidx)} = mscoh{mscohfeatureGroup(chanpairidx)-48};
                    
                end
            end
        end
        
        %CPSD - cross power spectral density features [55:73]
        if isempty(channelPairing) %Compute all features
            [ anglephase, cpsdmag,  cpsdmagt, cpsdmagb]  = cpsdfeatures( data, par, minfeat);
        else
            cpsdmember = ismember(featureGroup,[55:73]);
            cpsdfeatureGroup = featureGroup(cpsdmember);
            if any(cpsdmember)
                cpsdchanpairs = channelPairing(cpsdmember,:);
                cpsdFeatureMatIdx = selectedFeatures(cpsdmember);
                %RestFeatureNumbers-    {anglephase(a1-5) [55:59], cpsdmag(a0-5) [60:65], cpsdmagt(t0-2)[66:68], cpsdmagb(t0-2)[69:73]}
                for chanpairidx = 1:size(cpsdchanpairs,1)
                    
                    [ anglephase, cpsdmag,  cpsdmagt, cpsdmagb]  = cpsdfeatures( data, par, minfeat,cpsdchanpairs(chanpairidx,:));
                    switch cpsdfeatureGroup(chanpairidx)
                        case num2cell(55:59)
                            FeatureMatAll{cpsdFeatureMatIdx(chanpairidx)} = anglephase{cpsdfeatureGroup(chanpairidx)-54};
                        case num2cell(60:65)
                            FeatureMatAll{cpsdFeatureMatIdx(chanpairidx)} = cpsdmag{cpsdfeatureGroup(chanpairidx)-59};
                        case num2cell(66:68)
                            FeatureMatAll{cpsdFeatureMatIdx(chanpairidx)} = cpsdmagt{cpsdfeatureGroup(chanpairidx)-65};
                        case num2cell(69:73)
                            FeatureMatAll{cpsdFeatureMatIdx(chanpairidx)} = cpsdmagb{cpsdfeatureGroup(chanpairidx)-68};
                    end
                    
                end
            end
        end
        
        if isempty(channelPairing)
            Features = { sumspind, sumlow, sumhigh, summid, sumpower, tsumlow,tsumhigh, tsummid, tsumpower, bsumlow, bsumhigh, bsummid, bsumpower,mscoh, ...
                anglephase, cpsdmag, cpsdmagt, cpsdmagb, numspindperchanall,ApproxEntropy,SampleEntropy, cwtpower, tcwtpower,bcwtpower}; %For Rest
        else
            Features = [];
        end
end

if isempty(channelPairing) %Compute all features
    FeatureMat = cell2mat(cellfun(@(x)cell2mat(x), Features, 'un',0)); %Turn cell array into matrix
else %specific features
    FeatureMat = cell2mat(FeatureMatAll); %For Rest
end
FeatureNames = featureGroup;