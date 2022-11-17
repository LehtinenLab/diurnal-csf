function out = stdshade2(amatrix,alpha,acolor)

F=1:size(amatrix,2);
amean = nanmean(amatrix,1); %get man over first dimension
astd = nanstd(amatrix,[],1)/sqrt(size(amatrix,1)); % to get sem shading

fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'FaceAlpha', alpha,'linestyle','none');

hold on;
plot(F,amean, 'color', acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
%hold off;

end