function [diffspind, c] = subchan(data);

c = nchoosek([1:size(data,1)],2);
diffspind = zeros(size(data,2),size(c,1)); %preallocating
for i = 1:size(c,1)
        diffspind(:,i) = data(c(i,1),:)-data(c(i,2),:);
%         chancomp(i,1) = left;
%         chancomp(i,2) = right;
%         i=i+1;
end

