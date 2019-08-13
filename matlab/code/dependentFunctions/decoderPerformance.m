function resolution = decoderPerformance(mdl)

chance = 1/length(mdl.gof.confusionMatrix); 

figure;
subplot(1,2,1)
predProb = mdl.gof.confusionMatrix ./ sum(mdl.gof.confusionMatrix); 
imagesc(predProb)
% caxis([0 max(predProb(:))])
caxis([0 .75])
set(gca,'xtick',[],'ytick',[])
xlabel('predicted');ylabel('true')
colorbar
axis square

subplot(1,2,2);
hold on;
binedges = max(mdl.io.Y.normal)-.5; 
plotedges = max(mdl.io.Y.normal)-1;
n = histcounts(mdl.io.trueXpredicted(:,2)-mdl.io.trueXpredicted(:,1),[-binedges:1:binedges]);
x = -plotedges:1:plotedges;
y = n/sum(n); 
plot(x,y,'k');
hold on; plot(x,ones(1,length(x)).*chance,'-.k')
axis square
set(gca,'ylim',[0 .4])
xlabel('distance of prediction from true');
ylabel('proportion of trials')
xdata = y(find(x==0):end);
ydata = x(find(x==0):end);
resolution = spline(xdata(1:9),ydata(1:9),chance);

 



% DEPRECATED 190524 since no longer doing angles 
% sampledSpanPerCell = cellfun(@(x) (max(x.allTouches.theta.range) - min(x.allTouches.theta.range)),selectedArray);
% binResolutionMean = (mean(sampledSpanPerCell)./size(predProb,1)); 
% binResolutionStDev = (std(sampledSpanPerCell)./size(predProb,1));
% 
% resolution = nan(size(predProb,1),1); 
% for i = 1:size(predProb,1)
%     if ~isempty(find(predProb(i,:)<chance))
%     resolution(i) = min(abs(i - find(predProb(i,:)<chance)));
%     end
% end
% statsResolution = [mean(resolution.*binResolutionMean) std(resolution.*binResolutionMean)];
% subplot(1,2,1)
% title(['chance = ' num2str(chance) '. Resolution = ' num2str(statsResolution(1)) '+/-' num2str(statsResolution(2))])
% 
% 
