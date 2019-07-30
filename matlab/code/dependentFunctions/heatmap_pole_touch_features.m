%How does each stimulus vary w/ pole position?
%dependent function to id all touches and pre/post decision touches
preDecisionTouches = preDecisionTouchMat(U);
viewWindow = [-25:50];

gc_vals = cell(length(U),1);
%stimulus and response variables definitions
for rec = 1:length(U)
    array = U{rec};
    [tVar] = atTouch_sorter(array,viewWindow,preDecisionTouches{rec});
    pole = tVar.allTouches.S_ctk(:,end);
    it = tVar.allTouches.S_ctk(:,1:end-1);
    dt = tVar.allTouches.dtS_ctk;
    
    %heat of just each feature correlated to pole
    corr_values = corr([pole it dt]);
    gc_vals{rec} = corr_values;
    
    
    %imagesc of heat for all feature correlation
    %     figure(42);clf
    %     imagesc(corr([pole it dt]))
    %     caxis([-1 1])
    %     colorbar
    %     colormap bone
    %     set(gca,'xticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames(1:end-1)]...
    %         ,'yticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames(1:end-1)])
    %     axis square
end

%% Population feature correlation
vect = abs(cell2mat(cellfun(@(x) x(:)',gc_vals,'uniformoutput',0)));
figure(67);clf
imagesc(reshape(nanmedian(vect),size(gc_vals{1})))

caxis([0 1])
colorbar
colormap bone
set(gca,'xticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames]...
    ,'yticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames])
axis square