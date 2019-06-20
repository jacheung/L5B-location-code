%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%% Top level parameters and definitions
touchWindow = [-25:50]; %window for analyses around touch

touchCells = touchCell(U,'off');
selectedCells = find(touchCells==1);


% Structure for quantifying tuning and evaluating decoding
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99,'off');
whisking = whisking_general(U,'off');
hilbertWhisking = whisking_hilbert(U,popV,'off');

%%
[rc]= numSubplots(length(U))

fieldsToCompare = hilbertWhisking.cellVarNames;
pctSamp = .01; 
for i = 1:4
    figure(230);clf
    featureMean = cellfun(@(x) nanmean(x),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    featureSEM = cellfun(@(x) nanstd(x) ./ sqrt(sum(~isnan(x))),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    sampleSize = cellfun(@(x) sum(~isnan(x)),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    tscore = cellfun(@(x) tinv(.0975,x-1), sampleSize,'uniformoutput',0);
    CI = cellfun(@(x,y) x.*y,featureSEM,tscore,'uniformoutput',0);
    

    selBins= cellfun(@(x) sum(~isnan(x))> round(sum(~isnan(x(:)))*pctSamp),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    x = hilbertWhisking.S_ctk.(fieldsToCompare{i});
    for d = 1:length(U)
        [WHilbertP{i}(d),~,~] = anova1(hilbertWhisking.R_ntk.(fieldsToCompare{i}){d},[],'off');

        subplot(rc(1),rc(2),d);
        if strcmp(fieldsToCompare{i},'phase')
            bar(x,featureMean{d}.*1000,'facecolor',[.8 .8 .8])
            hold on; shadedErrorBar(x,smooth(featureMean{d}.*1000,3),smooth(CI{d}(selBins{d}).*1000,3),'g')
            set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
        else
            bar(1:sum(selBins{d}),featureMean{d}(selBins{d}).*1000,'facecolor',[.8 .8 .8])
            hold on; shadedErrorBar(1:sum(selBins{d}),smooth(featureMean{d}(selBins{d}).*1000,3),smooth(CI{d}(selBins{d}).*1000,3),'g')
            xtv = x(selBins{d}); 
            set(gca,'xtick',1:3:length(xtv),'xticklabel',xtv(1:3:end),'ylim',[0 (max(featureMean{d}(selBins{d}).*1000)+.01).*1.2])
        end
    end
    
    suptitle(fieldsToCompare{i})
    
   % print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} '_whisking'],'-dpng')
   % print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} '_whisking'],'-depsc')
end

%%
wOn = [find(whisking.matrix(1,:)==1) find(whisking.matrix(1,:)==-1)];
%wOn = 1:length(U); 
hilbTune = cell2mat(WHilbertP')<.01;


ampOnly = intersect(find(hilbTune(1,:)),find(sum(hilbTune)==1));
ampT = intersect(wOn,ampOnly);

mpOnly = intersect(find(hilbTune(2,:)),find(sum(hilbTune)==1));
mpT = intersect(wOn,mpOnly);

pOnly = intersect(find(hilbTune(3,:)),find(sum(hilbTune)==1));
pT = intersect(wOn,pOnly);

angleOnly = intersect(find(hilbTune(4,:)),find(sum(hilbTune)==1));
aT = intersect(wOn,angleOnly);

for d = 0:4
    tuned{d+1} = intersect(wOn,find(sum(hilbTune)==d))
    tuneNum(d+1) = sum(numel(intersect(wOn,find(sum(hilbTune)==d))));
end

%%
pTune = cell2mat(WHilbertP');
chosenPTune = pTune(:,wOn)<.01; %looking only at cells that are whisking excited and inhibited. 

ampOnly = intersect(find(chosenPTune(1,:)),find(sum(chosenPTune)==1));


mpOnly = intersect(find(chosenPTune(2,:)),find(sum(chosenPTune)==1));

pOnly = intersect(find(chosenPTune(3,:)),find(sum(chosenPTune)==1));


angleOnly = intersect(find(chosenPTune(4,:)),find(sum(chosenPTune)==1));

twofeats = find(sum(chosenPTune)==2);
threefeats = find(sum(chosenPTune)==3);
allfeats = find(sum(chosenPTune)==4);
unmod = find(sum(chosenPTune)==0);

plotSort = [unmod ampOnly mpOnly pOnly angleOnly twofeats threefeats allfeats];

sortedMap = chosenPTune(:,plotSort);
figure(9);clf
imagesc(sortedMap')
set(gca,'xtick',1:4,'xticklabel',{'phase','amp','midpoint','angle'})



%% bayesian adaptive regression spline fitting 
[rc]= numSubplots(length(U))

fieldsToCompare = hilbertWhisking.cellVarNames;
pctSamp = .01; 
for i = 3
    featureMean = cellfun(@(x) nansum(x),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    sampleSize = cellfun(@(x) sum(~isnan(x)),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    selBins= cellfun(@(x) sum(~isnan(x))> round(sum(~isnan(x(:)))*pctSamp),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    x = hilbertWhisking.S_ctk.(fieldsToCompare{i});
     

    for g = 1:length(featureMean)
        counts = featureMean{g}(selBins{g});
        stim = x(selBins{g});
        
        if sum(counts)>numel(stim)
            fit1 = barsP(counts,[stim(1) stim(end)],round(median(sampleSize{g}(selBins{g}))));
            subplot(rc(1),rc(2),g)
            
            shadedErrorBar(stim,fit1.mean,abs(fit1.confBands(:,1) - fit1.mean),'g')
        else
            disp(['skipping cell ' num2str(g) ' b/c too few spikes'])
        end
    end
end
        
        

%% slow x fast tuning 
for g = 1:length(U)
    phase = hilbertWhisking.S_ctk.raw{g}(:,3);
    amp = hilbertWhisking.S_ctk.raw{g}(:,1);
    midpoint = hilbertWhisking.S_ctk.raw{g}(:,2); 
    angle = hilbertWhisking.S_ctk.raw{g}(:,4); 
    spikes = hilbertWhisking.R_ntk.raw{g}(:);
    spksT = find(spikes == 1);
    noSpksT = find(spikes == 0); 
    
    featX = phase; 
    featY = angle; 
    featZ = phase; 
    
    if ~isempty(spksT)
    figure(230);clf
     hold on; scatter(featX(noSpksT),featY(noSpksT),'ko','markeredgecolor',[.9 .9 .9])
    scatter(featX(spksT),featY(spksT),'ko','filled')
    title(num2str(corr(featX(spksT),featY(spksT),'rows','complete')));
    end

    
    pause
    
    
end
    








