clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

Type = {'Chat';'PV';'SOM';'VGLUT2';'Ctrl';};

BehaviorEvent = {'NovelExp';'Grooming';'QuietAwakeTime';'Speed';};
Excel = fullfile(WholePath,'BehaviorResults.xlsx');
Excel0 = fullfile(WholePath,'NMF_weight_perc_new.xlsx');

NCK = nchoosek(1:6,3);

ColorVec = {[240,150,000];[000,145,058];[003,110,183];[145,007,130];[0 0 0]};
OutputExcel = fullfile(WholePath,'Encode_Decode_NMF.xlsx');
for bl=1:numel(BehaviorEvent)
    
    data = xlsread(Excel,BehaviorEvent{bl});
    X_ = reshape(data,6,[],5);
    Qx = squeeze(nanmean(X_,2));

    for nl = 1:6

        data = xlsread(Excel0,['NMF_C',num2str(nl)]);
        Y_ = reshape(data,6,[],5);
        Qy = squeeze(nanmean(Y_,2));
        eval(['Qy',num2str(nl),'=Qy;']);
    end
        
    BetaALL=zeros([size(NCK,1),6]);
    PredCCALL=zeros([size(NCK,1),1]);
    for nloop=1:size(NCK,1)
        EncData = NCK(nloop,:);
        DecData = setdiff(1:6,EncData);
        
        Wx = Qx(EncData,:);
        Wy1 = Qy1(EncData,:);
        Wy2 = Qy2(EncData,:);
        Wy3 = Qy3(EncData,:);
        Wy4 = Qy4(EncData,:);
        Wy5 = Qy5(EncData,:);
        Wy6 = Qy6(EncData,:);
        [Beta0] = ridge(Wx(:),[ones(numel(Wy1),1),Wy1(:),Wy2(:),Wy3(:),Wy4(:),Wy5(:),Wy6(:)],0:0.01:0.1);
        Beta = Beta0(:,end);
        
        Dx = Qx(DecData,:);
        Dy1 = Qy1(DecData,:);
        Dy2 = Qy2(DecData,:);
        Dy3 = Qy3(DecData,:);
        Dy4 = Qy4(DecData,:);
        Dy5 = Qy5(DecData,:);
        Dy6 = Qy6(DecData,:);
        y = [ones(numel(Dy1),1),Dy1(:),Dy2(:),Dy3(:),Dy4(:),Dy5(:),Dy6(:)]*Beta;
        Dx_ = reshape(y,size(Dx));
        
        PredCC = corr(Dx(:),Dx_(:));
        
        BetaALL(nloop,:) = Beta(2:end);
        PredCCALL(nloop,:) = PredCC;
    end
    Enc = [mat2cell(NCK,ones(size(NCK,1),1),size(NCK,2))];
    TXT = {'EncodeScan','Beta1','Beta2','Beta3','Beta4','Beta5','Beta6','PredictedCC'};
    xlswrite(OutputExcel,TXT,BehaviorEvent{bl},'A1');
    xlswrite(OutputExcel,[BetaALL,PredCCALL],BehaviorEvent{bl},'B2');
end


% EigVec Map of cell-type OG
CellTypeVec = [mean(Qy1,1);mean(Qy2,1);mean(Qy3,1);mean(Qy4,1);mean(Qy5,1);mean(Qy6,1)];
for cl=1:size(CellTypeVec,2)
    CellTypeVec(:,cl) = CellTypeVec(:,cl)/norm(CellTypeVec(:,cl),2)*1.5;
end
Vec = mat2cell(CellTypeVec,size(CellTypeVec,1),ones(size(CellTypeVec,2),1));

ColorVec = {[240,150,000];[000,145,058];[003,110,183];[145,007,130];[0 0 0]};
F = MY_3D_sphere_figure(Vec(1:4),ColorVec(1:4));
saveas(F,fullfile('F:\BF_optogentics\SubMit\V5\add\',['NMF_sphere_vector.tif']));
for cl=1:size(CellTypeVec,2)
    F = MY_3D_sphere_figure(Vec(cl),ColorVec(cl));
    x=mean(Vec{cl},2);
    %x=x/norm(x,2);
    x([1 3 4 2 5 6])'

    saveas(F,fullfile('F:\BF_optogentics\SubMit\V5\add\',['NMF_sphere_vector_',(Type{cl}),'.tif']));
end


% 
TIFtype={'NovelExp';'Grooming';'QW';'Locomotion'};
Excel = fullfile(WholePath,'Encode_Decode_NMF.xlsx');

for bl=1:3%numel(BehaviorEvent)
    data = xlsread(Excel,BehaviorEvent{bl});
    X = mean(data,1);
    
    BehaVec = zeros(6,size(data,1));
    for dl=1:size(data,1)
        X = data(dl,1:6);
        X = X/norm(X,2);
        BehaVec(:,dl) = X(1:6)*1.8;
    end
    x=mean(BehaVec,2);
    x=x/norm(x,2);x([1 3 4 2 5 6])'
    Vec = mat2cell(BehaVec,size(BehaVec,1),ones(size(BehaVec,2),1));
    ColorVec = {[240,150,000];[000,145,058];[003,110,183];[145,007,130];[0 0 0]};
    F = MY_3D_sphere_figure_v2(Vec,repmat(ColorVec(bl),1,numel(Vec)));
    saveas(F,fullfile('F:\BF_optogentics\SubMit\V5\add\',[TIFtype{bl},'_sphere_vector.tif']));
    close all
end




%% OG & Behavior angle
CellTypeVec = [mean(Qy1,1);mean(Qy2,1);mean(Qy3,1);mean(Qy4,1);mean(Qy5,1);mean(Qy6,1)];
for cl=1:size(CellTypeVec,2)
    CellTypeVec(:,cl) = CellTypeVec(:,cl)/norm(CellTypeVec(:,cl),2)*1.5;
end


TIFtype={'NovelExp';'Grooming';'QW';'Locomotion'};
Excel = fullfile(WholePath,'Encode_Decode_NMF.xlsx');
BehaVec = zeros(6,numel(BehaviorEvent));
for bl=1:numel(BehaviorEvent)
    data = xlsread(Excel,BehaviorEvent{bl});
    X = median(data,1);
    X = X/norm(X,2);
    BehaVec(:,bl) = X(1:6);
end

clear ang
F = figure; a=0;
for ix=1:size(CellTypeVec,2)
    for iy=1:size(BehaVec,2)
        OA=CellTypeVec(:,ix);OB=BehaVec(:,iy);
        x=acos(dot(OA,OB)/(norm(OA,2)*norm(OB,2)));
        ang(ix,iy) = x/pi*180;
        
        a=a+1;
        subplot(size(CellTypeVec,2),size(BehaVec,2),a);
        pie([360-ang(ix,iy),ang(ix,iy)])
    end
end


%% angle CC
TIFtype={'NovelExp';'Grooming';'QW';'Locomotion'};
Excel = fullfile(WholePath,'BehaviorResults.xlsx');
for bl=1:numel(BehaviorEvent)
    
    data = xlsread(Excel,BehaviorEvent{bl});
    X_ = reshape(data,6,[],5);
    Qx = squeeze(nanmean(X_,2));

    for nl = 1:6

        data = xlsread(Excel0,['NMF_C',num2str(nl)]);
        Y_ = reshape(data,6,[],5);
        Qy = squeeze(nanmean(Y_,2));
        eval(['Qy',num2str(nl),'=Qy;']);
    end
    
    
    clear ang BehaviorEmpirical
    for nloop=1:size(NCK,1)
        EncData = NCK(nloop,:);
        DecData = setdiff(1:6,EncData);
        
        Wx = Qx(EncData,:);
        Wy1 = Qy1(EncData,:);
        Wy2 = Qy2(EncData,:);
        Wy3 = Qy3(EncData,:);
        Wy4 = Qy4(EncData,:);
        Wy5 = Qy5(EncData,:);
        Wy6 = Qy6(EncData,:);

        [Beta0] = ridge(Wx(:),[ones(numel(Wy1),1),Wy1(:),Wy2(:),Wy3(:),Wy4(:),Wy5(:),Wy6(:)],0:0.01:0.1);
        Beta = Beta0(:,end);
       
        Dx = Qx(DecData,:);
        Dy1 = Qy1(DecData,:);
        Dy2 = Qy2(DecData,:);
        Dy3 = Qy3(DecData,:);
        Dy4 = Qy4(DecData,:);
        Dy5 = Qy5(DecData,:);
        Dy6 = Qy6(DecData,:);
        y = [ones(numel(Dy1),1),Dy1(:),Dy2(:),Dy3(:),Dy4(:),Dy5(:),Dy6(:)]*Beta;
        Dx_ = reshape(y,size(Dx));
        
        for cl=1:5
            OA = CellTypeVec(:,cl);
            OA = OA/norm(OA,2);
            OB = Beta(2:end);
            x=acos(dot(OA,OB)/(norm(OA)*norm(OB)));
            ang(nloop,cl) = x/pi*180;
        end
        BehaviorEmpirical(nloop,:) = mean(Dx,1);
    end
    
    ColorVec = {[240,150,000];[000,145,058];[003,110,183];[145,007,130];[100 100 100]};
     F = figure('Position', [680 647 357 331]);
     for ql=1:size(Qx,2)
        plot(ang(:,ql),BehaviorEmpirical(:,ql),'o',...
            'MarkerFaceColor',ColorVec{ql}/255,'MarkerEdgeColor','none');
        hold on
     end
    p = polyfit(ang(:),BehaviorEmpirical(:),1);
    y_ = polyval(p,[min(ang(:)),max(ang(:))]);
    hold on;plot([min(ang(:)),max(ang(:))],y_,'k-','linewidth',1.5)
    set(gca,'tickdir','out','box','off','ticklength',[0.04 0.1]);
    saveas(F,fullfile('F:\BF_optogentics\SubMit\V5\add\',[TIFtype{bl},'_angle_CC.emf']));
    close all

    
    [CC,p] = corr(ang(:),BehaviorEmpirical(:))

end


%% NMF predicted spatial map
NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
lmask = spm_read_vols(spm_vol(fullfile(codepath,'lmask_Mouse_v38.nii')));
for bl=1:numel(BehaviorEvent)
    weight = BehaVec(:,bl);
    dest = fullfile(WholePath,'DimReduc_SC_log_voxwise','Montaged','N6');
    w=6;Vf=[];
    for wl=1:w
        C0 = spm_read_vols(spm_vol(fullfile(dest,['NMF.nii,',num2str(wl)]))); 
        C = imresize3(C0,size(lmask));
        C = flip(C,3);
        VC=fmask(C,lmask);
        %VC = VC/max(VC);
        Vf=[Vf,VC];
    end

    
    Vp = Vf*weight;

    V3D = funmask(Vp,lmask);
    
    hdr = ihdr; hdr.dt=[16 0];
    hdr.fname = fullfile(WholePath,['N_',BehaviorEvent{bl},'_1.nii']);
    spm_write_vol(hdr,V3D);
    
    dest = fullfile(WholePath);
    template = [codepath,'\meanT2.nii'];
    spmT_file = fullfile(WholePath,['N_',BehaviorEvent{bl},'_1.nii']);
    lmask0 = [codepath,'\lmask_Mouse_v38.nii'];
    Colormap(   'statfile',spmT_file,...
        'bgmfile',template,...
        'slice',9:2:35-4*2,...
        'bar_value',[-10 -2 2 10]/3,...
        'dest',dest,...
        'mapname',['N_',BehaviorEvent{bl},'_1'],...
        'denoi_profile',lmask0,...
        'cluster',10);

end



