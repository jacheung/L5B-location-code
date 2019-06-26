
%% Build feature at touch and response in touch response window 

% Defining touch response
U = defTouchResponse(U,.99,'off');

% viewing window around touch 
touchWindow = [-25:50];

preDecisionTouches = preDecisionTouchMat(U);
for rec = 1:length(U)
 [tVar] = atTouch_sorter(U{rec},touchWindow,preDecisionTouches{rec});

 if isfield(U{rec}.meta,'responseWindow')
     tResponse = U{rec}.meta.responseWindow+find(touchWindow==0) ;
%      tResponse = [-25 -5] + find(touchWindow==0); %control to find off touch responses
 else
     tResponse = [5 35] + find(touchWindow==0);
 end
  response = mean(tVar.allTouches.R_ntk(:,tResponse(1):tResponse(2)),2);
 
  hilbertTouch.S_ctk{rec} = tVar.allTouches.S_ctk(:,[3, 4, 5, 1]);
  hilbertTouch.R_ntk{rec} = response; 
  
end

%% 1 dimensional hilbert features
[rc]= numSubplots(length(U))

fieldsToCompare = fields(hilbertTouch.R_ntk.allTouches);

THilbertP = cell(1,numel(fieldsToCompare)); 
peakResponse = cell(1,numel(fieldsToCompare)); 

pctSamp = .01; 
smoothParam = .7;
for i = 1:4
    figure(230);clf
    featureMean = cellfun(@(x) nanmean(x),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    featureSEM = cellfun(@(x) nanstd(x) ./ sqrt(sum(~isnan(x))),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    sampleSize = cellfun(@(x) sum(~isnan(x)),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    tscore = cellfun(@(x) tinv(.0975,x-1), sampleSize,'uniformoutput',0);
    CI = cellfun(@(x,y) x.*y,featureSEM,tscore,'uniformoutput',0);
    

    selBins= cellfun(@(x) sum(~isnan(x))> round(sum(~isnan(x(:)))*pctSamp),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    x = hilbertTouch.S_ctk.(fieldsToCompare{i});
    for d = 1:length(U)
        [THilbertP{i}(d),~,~] = anova1(hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}){d},[],'off');

        subplot(rc(1),rc(2),d);
        if strcmp(fieldsToCompare{i},'phase')
            f = fit((1:numel(x))',(featureMean{d}.*1000)','smoothingspline','SmoothingParam',smoothParam);
            fCI = fit((1:numel(x))',(CI{d}.*1000)','smoothingspline','SmoothingParam',smoothParam);
            fitline = f(1:numel(x));
            fCIline = fCI(1:numel(x)); 
            
            bar(x,featureMean{d}.*1000,'facecolor',[.8 .8 .8])
            hold on; shadedErrorBar(x,fitline,fCIline,'g')
            set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
            
            [~,idx] = max(fitline);
            peakResponse{i}(d) = x(idx); 
        else
            xtv = x(selBins{d}); 
            
            bar(xtv,featureMean{d}(selBins{d}).*1000,'facecolor',[.8 .8 .8])
            f = fit(xtv',(featureMean{d}(selBins{d}).*1000)','smoothingspline','SmoothingParam',smoothParam);
            fCI = fit(xtv',(CI{d}(selBins{d}).*1000)','smoothingspline','SmoothingParam',smoothParam);
            fitline = f(xtv(1):xtv(end));
            fCIline = fCI(xtv(1):xtv(end)); 
            hold on; shadedErrorBar(xtv(1):xtv(end),fitline,fCIline,'g')
           
            set(gca,'xtick',xtv(1):3:xtv(end),'ylim',[0 (max(featureMean{d}(selBins{d}).*1000)+.01).*1.2])
            
            [~,idx] = max(fitline);
            ranges = xtv(1):xtv(end);
            peakResponse{i}(d) = ranges(idx); 
        end
    end
    
    suptitle(fieldsToCompare{i})
    
%    print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} '_whisking'],'-dpng')
%    print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} '_whisking'],'-depsc')
end


whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
hilbertTouch = tuningQuantification(U,popV,selectedCells,fieldsList([ 1 3 4 5]),whichTouches,touchWindow,'off');

touchFields = fields(hilbertTouch.R_ntk.allTouches);
for d = 1:lenght(touchFields)
    c_feat = hilbertTouch.R_ntk.allTouches.(touchFields{d});
    TouchR = cellfun(@(x) nanmean(x,2),c_feat,'uniformoutput',0) * 1000;
    responseSEM = cellfun(@(x) nanstd(x,[],2) ./ nansum(~isnan(x),2),c_feat,'uniformoutput',0) * 1000;
    
    for b = 1:length(c_feat)
        c_cell=c_feat{b};
        
        if ~isempty(c_cell)
        p_val{d}(b) = anova1(c_feat{b},[],'off');
        
        
        
            
           

%% %2dimensional or 3dimensional scatter of features at touch
figure(2332);clf
for k = 1:length(U)
    phase = hilbertTouch.S_ctk{k}(:,3);
    amp = hilbertTouch.S_ctk{k}(:,1);
    midpoint = hilbertTouch.S_ctk{k}(:,2);
    angle = hilbertTouch.S_ctk{k}(:,4);
    
    responses = normalize_var(hilbertTouch.R_ntk{k},.05,.95);
    [~,idx] = sort(responses);
    
    featX = phase;
    featY = angle;
    featZ = angle;
    
    colors = repmat(1 - responses,1,3);
    
    if isfield(U{k}.meta,'responseWindow')
        figure(2332);subplot(6,10,k)
        scatter(featX(idx),featY(idx),50,colors(idx,:),'filled');
        set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
    end
end

%%
for k = 1
    dmatX = hilbertTouch.S_ctk{k};
    nDmatX = (dmatX-nanmean(dmatX))./nanstd(dmatX);
    
    y = hilbertTouch.R_ntk{k};
    
    fitlm(nDmatX,y)
end
    