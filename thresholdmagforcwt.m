function [sumlow, sumhigh, summid] = thresholdmagforcwt(data)

numdims = 3; %third dimension will always be number of comparisons (regardless of number of comparisons) i.e. = 1 when only one comparison

sumlow = zeros(1,size(data,numdims));
sumhigh = zeros(1,size(data,numdims));
summid = zeros(1,size(data,numdims));

 for comp = 1:size(data,numdims)
    spectro = data(:,:,comp);
    allspectro = spectro(:);
    meanspectro = mean(allspectro);
%     sumspectro = sum(allspectro);
    stdspectro = std(allspectro);
    upperthresh =  meanspectro+1*stdspectro;
    lowerthresh =  meanspectro-1*stdspectro;
    above = spectro;
    above(above<=upperthresh) = 0;
    below = spectro;
    below(below>=lowerthresh) = 0;
    between = spectro;
    between(between<lowerthresh|between>upperthresh)=0;
    alllow = below(below ~=0);
    allhigh = above(above~=0);
    allmid = between(between~=0);
    if isempty(alllow)
        alllow = 0;
    end
    if isempty(allhigh)
        allhigh = 0;
    end
    if isempty(allmid)
        allmid = 0;
    end
    sumlow(comp) = sum(alllow);
%     abssumlow(comp) = abs(sumlow(comp));
    sumhigh(comp) = sum(allhigh);
    summid(comp) = sum(allmid);
%     meanlow(comp) = mean(alllow);
%     absmeanlow(comp) = abs(meanlow(comp));
%     meanhigh(comp) = mean(allhigh);
%     meanmid(comp) = mean(allmid);
end