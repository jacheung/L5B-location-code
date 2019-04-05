function decoderPerformance(mdl)

figure;
subplot(1,2,1)
predProb = mdl.gof.confusionMatrix ./ sum(mdl.gof.confusionMatrix); 
imagesc(predProb)
caxis([0 max(predProb(:))])
set(gca,'xtick',[],'ytick',[])
xlabel('predicted');ylabel('true')
colorbar
axis square

subplot(1,2,2);
hold on;
binedges = max(mdl.io.Y.normal)-.5; 
plotedges = max(mdl.io.Y.normal)-1;
n = histcounts(mdl.io.trueXpredicted(:,2)-mdl.io.trueXpredicted(:,1),[-binedges:1:binedges]);
plot([-plotedges:1:plotedges],n/sum(n),'k');
axis square
xlabel('distance of prediction from true');
ylabel('proportion of trials')