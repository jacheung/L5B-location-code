
%Build feature at touch and response in touch response window 

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

%%
figure(2332);clf
for k = 1:length(U)
    phase = hilbertTouch.S_ctk{k}(:,3); 
    amp = hilbertTouch.S_ctk{k}(:,1); 
    midpoint = hilbertTouch.S_ctk{k}(:,2); 
    angle = hilbertTouch.S_ctk{k}(:,4); 
    
    
    responses = normalize_var(hilbertTouch.R_ntk{k},.05,.95);
    [~,idx] = sort(responses); 
    
    featX = phase;
    featY = amp
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
    