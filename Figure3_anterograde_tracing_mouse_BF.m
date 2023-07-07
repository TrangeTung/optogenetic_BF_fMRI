clc;clear
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

WholePath = 'F:\BF_optogentics\\';
SC_Excel = fullfile(codepath,'Allen_mouse_SC_213.xlsx');
SC = xlsread(SC_Excel);


%% Anterograde tracing
%

Type = {'chat';'vglut';'pv';'sst'};
for tloop=1:numel(Type)

    filename = fullfile(WholePath,'BF_structure_connectivity',['nR_',Type{tloop},'_with_mask.nii']);
    template = fullfile(codepath,'meanT2.nii');
    
    
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Img_RGB = MY_display_function_map_3D_nothreshold(Func_Img_3D/(255/2),[0 1],10:2:36);
    
    I = cat(1,Img_RGB(:,1:size(Img_RGB,2)/2,:),Img_RGB(:,1+size(Img_RGB,2)/2:end,:));
    imwrite(uint8(I),fullfile(WholePath,['tmap_',Type{tloop},'.tiff']));
1;


barlim = [-10 -2.89 +1 50];
Func3D = filename;
Temp3D = fullfile(codepath,'Template_Mouse_v38.nii');
[f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,1,0.8);
[f2_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,barlim,codepath,1,0.6);

I = cat(1,f1_rgb.cdata,f2_rgb.cdata);

filename = fullfile(WholePath,['3D_tmap_',Type{tloop},'.tiff']);
imwrite(uint8(I),filename);
close all;
    
end




%}


%% Quantitative comparison
%
Excel = fullfile(codepath,'Label_213_v38.xlsx');
[~,~,CellData] = xlsread(Excel);
ExpTable_213 = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);

Type = {'chat';'vglut';'pv';'sst'};
for tloop=1:numel(Type)
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    Vs = zeros(213,1);
    for loop=1:213
        lmask = Labels==loop;
        Vs(loop) = nanmedian(fmask(Trace,lmask))/255;
    end
    
    filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'beta_0001.nii');
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Vf = zeros(213,1);
    VoxNum = zeros(213,1);
    for loop=1:213
        lmask = Labels==loop;
        Vf(loop) = nanmean(fmask(Func_Img_3D,lmask));
        VoxNum(loop) = numel(find(lmask==1));
    end
    ResizeVox = (3*3*9.47)/(2*2*8);
    
    
    dest = fullfile(WholePath,'BF_structure_connectivity');
    VoxThr = 100;
    Vbi = VoxNum/ResizeVox>VoxThr&~isnan(Vf);
    x_ = Vs(Vbi)*100;
    y_ = Vf(Vbi);
    ROIs = ExpTable_213.abbre(Vbi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VB = VoxNum/ResizeVox>0&~isnan(Vf)&Vs~=0; % 1st
    SCx = SC(VB,:);
    VC = mean(SCx,1)'>(mean(mean(SCx,1))+0*std(mean(SCx,1))) & ~VB;
    VN = mean(SCx,1)'<=(mean(mean(SCx,1))+0*std(mean(SCx,1))) & ~VB;
    
    F = figure('Position',[280 261 656 424]);
    subplot('Position',[0.4 0.1 0.55 .8]);
    semilogx(Vs(VB)*100,Vf(VB),'o','MarkerEdgecolor','none',...
        'MarkerFaceColor',[1 0 0]); hold on;
    scatter(Vs(VN)*100+10^-1,Vf(VN),'o','MarkerEdgecolor','none',...
        'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.2); hold on;
    scatter(Vs(VC)*100+10^-1,Vf(VC),'o','MarkerEdgecolor','none',...
        'MarkerFaceColor',[1 0 1],'MarkerFaceAlpha',0.2); hold on;
    %plot(x_(x_==0)+10^-4,y_(x_==0),'o','MarkerEdgecolor','none',...
    %    'MarkerFaceColor',[.5 .5 .5]); hold on;
    xlim([0+10^-1.1 100]);ylim([floor(min(y_)) ceil(max(y_))]);
    xlabel('Normalized output density (%)'); ylabel('BOLD amplitude (β)');
    set(gca,'fontsize',13,'linewidth',1.5,'box','off',...
        'tickdir','out');%,'xaxislocation','left'
    
    f = fittype('a*log10(x)+b');
    [fit1,gof] = fit(x_(x_~=0),y_(x_~=0),f,'StartPoint',[x_(1) y_(1)]);
    fdata = feval(fit1,x_(x_~=0));
    co = confint(fit1);
    plot(x_(x_~=0),fdata,'r','linewidth',2);
    ups=co(1,1)*log10(sort(x_(x_~=0)))+co(1,2);
    dws=co(2,1)*log10(sort(x_(x_~=0)))+co(2,2);
    title({['y=',num2str(fit1.a),'*log10(x)+',num2str(fit1.b)];...
        ['R2 = ',num2str(gof.rsquare)]},'fontsize',10);    
    
    subplot('Position',[0 0.1 0.05 .8]);
    h=histogram(Vf(VB),40,'Edgecolor','none','FaceColor',[1 0 0],'FaceAlpha',0.2);
    hold on;
    y=medfilt1(h.Values,5);[bx,ax]=butter(3,0.2);yf=filtfilt(bx,ax,y);yf(yf<0)=0;
    plot(h.BinEdges(1:end-1)+h.BinWidth/2,yf,'color',[1 0 0]);
    xlim([floor(min(y_)) ceil(max(y_))]);ylim([0 15])
    view([-90 90]);axis off;
    subplot('Position',[0.05 0.1 0.05 .8]);
    boxplot(Vf(VB),'color',[1 0 0],'width',2)
    ylim([floor(min(y_)) ceil(max(y_))]);axis off;
    
    subplot('Position',[0.2 0.1 0.05 .8]);
    h=histogram(Vf(VN),40,'Edgecolor','none','FaceColor',[.5 .5 .5],'FaceAlpha',0.2);
    hold on;
    y=medfilt1(h.Values,5);[bx,ax]=butter(3,0.2);yf=filtfilt(bx,ax,y);yf(yf<0)=0;
    plot(h.BinEdges(1:end-1)+h.BinWidth/2,yf,'color',[.5 .5 .5]);
    xlim([floor(min(y_)) ceil(max(y_))]);ylim([0 15])
    view([-90 90]);axis off;
    subplot('Position',[0.25 0.1 0.05 .8]);
    boxplot(Vf(VN),'color',[.5 .5 .5],'width',2)
    ylim([floor(min(y_)) ceil(max(y_))]);axis off;

    subplot('Position',[0.1 0.1 0.05 .8]);
    h=histogram(Vf(VC),40,'Edgecolor','none','FaceColor',[1 0 1],'FaceAlpha',0.2);
    hold on;
    y=medfilt1(h.Values,5);[bx,ax]=butter(3,0.2);yf=filtfilt(bx,ax,y);yf(yf<0)=0;
    plot(h.BinEdges(1:end-1)+h.BinWidth/2,yf,'color',[1 0 1]);
    xlim([floor(min(y_)) ceil(max(y_))]);ylim([0 15])
    view([-90 90]);axis off;
    subplot('Position',[0.15 0.1 0.05 .8]);
    boxplot(Vf(VC),'color',[1 0 1],'width',2)
    ylim([floor(min(y_)) ceil(max(y_))]);axis off;
    
    [~,p1]=ttest2(Vf(VB),Vf(VC),'tail','right');
    [~,p2]=ttest2(Vf(VB),Vf(VN),'tail','right');
    [~,p3]=ttest2(Vf(VC),Vf(VN),'tail','right');
    [p1,p2,p3]
    saveas(F,fullfile(dest,['SC_FC_',Type{tloop},'.emf']))

    1;
end
%}


%% Spatial Overlap 3D

Excel = fullfile(codepath,'Label_213_v38.xlsx');
[~,~,CellData] = xlsread(Excel);
ExpTable_213 = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Colormap_3Dviewer','R_Label_Mouse_213.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);
NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels_v38 = spm_read_vols(ihdr);

Type = {'chat';'vglut';'pv';'sst'};
for tloop=1:numel(Type)
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    Vs = zeros(213,1);
    for loop=1:213
        lmask = Labels_v38==loop;
        Vs(loop) = nanmedian(fmask(Trace,lmask))/255;
    end
    
    filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'beta_0001.nii');
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Pos = Func_Img_3D>1.7;
    Neg = Func_Img_3D<-0.7;
    FuncMap = double(Pos)+(-double(Neg));
    
    
    Vf = zeros(213,1);
    VoxNum = zeros(213,1);
    for loop=1:213
        lmask = Labels_v38==loop;
        Vf(loop) = nanmedian(fmask(FuncMap,lmask));
    end
    
    
    dest = fullfile(WholePath,'BF_structure_connectivity');
    F = figure('Color',[1 1 1],'Position', [680 52 799 926]);
    % all brain
    ALL = Labels~=0; ALL = smooth3(ALL,'box',11);
    p0 = patch(isosurface(ALL,0.5));
    isonormals(ALL,p0);
    p0.FaceColor = 'black';
    p0.EdgeColor = 'none';
    p0.FaceAlpha=0.05;
    axis equal; axis tight
    axis off;
    hold on;
    
    
    %  only positive Func 
    Lv = Vf>0; Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask1 = smooth3(Rmask,'box',3);
    p = patch(isosurface(Rmask1,0.5));
    isonormals(Rmask1,p);
    p.FaceColor = [1,0,0];p.EdgeColor = 'none';
    p.FaceAlpha = 0.10;
    axis equal; axis tight; axis off;
    
    %  only negative Func 
    Lv = Vf<0; Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask2 = smooth3(Rmask,'box',3);
    p = patch(isosurface(Rmask2,0.5));
    isonormals(Rmask2,p);
    p.FaceColor = [0,0,1];p.EdgeColor = 'none';
    p.FaceAlpha = 0.10;
    axis equal; axis tight; axis off;
    
    %  only projection color
    Lv = Vs>0; Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask3 = smooth3(Rmask,'box',3);
    p = patch(isosurface(Rmask3,0.5));
    isonormals(Rmask3,p);
    p.FaceColor = [0,1,0];p.EdgeColor = 'none';
    p.FaceAlpha = 0.10;
    axis equal; axis tight; axis off;
    
    Rmask = smooth3(Rmask1&Rmask3,'box',3);
    p = patch(isosurface(Rmask,0.5));
    isonormals(Rmask,p);
    p.FaceColor = [1,1,0];p.EdgeColor = 'none';
    p.FaceAlpha = 0.050;
    axis equal; axis tight; axis off;
    
    Rmask = smooth3(Rmask2&Rmask3,'box',3);
    p = patch(isosurface(Rmask,0.5));
    isonormals(Rmask,p);
    p.FaceColor = [0,1,1];p.EdgeColor = 'none';
    p.FaceAlpha = 0.050;
    axis equal; axis tight; axis off;
    
    view([-90 180]);
    saveas(F,fullfile(dest,['ColorMap1_3DViewer_',Type{tloop},'_dorsal.tif']));
    view([0 180]);
    saveas(F,fullfile(dest,['ColorMap1_3DViewer_',Type{tloop},'_lateral.tif']));
    
    A = imread(fullfile(dest,['ColorMap1_3DViewer_',Type{tloop},'_dorsal.tif']));
    B = imread(fullfile(dest,['ColorMap1_3DViewer_',Type{tloop},'_lateral.tif']));
    C = cat(2,A(1:end,161:1200,:),B(1:end,231:1030,:));
    imwrite(C,fullfile(dest,['ColorMap1_3DViewer_',Type{tloop},'_merge.tif']));
    close all
    1;
end



function Img_RGB = MY_display_function_map_3D_nothreshold(Func_Img_3D,bar_value,slice)

map_nothre_reshape = reshape(permute(flip(Func_Img_3D(:,:,slice),1),[2 1 3]),[size(Func_Img_3D,2) size(Func_Img_3D,1)*numel(slice)]);

MyMap = [ones(1,64);(63:-1:0)/64;(63:-1:0)/64]';

defaultMap = MyMap;
nothrebar = [min(bar_value(:)) max(bar_value(:))];
tmap = ones([size(map_nothre_reshape),3]);
map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
map_normalize(map_normalize<=0) = 0.000000001;
map_normalize(map_normalize>=1) = 1;
nan_mask = isnan(map_normalize);
map_normalize(isnan(map_normalize)) = 0.000000001;
tmap(:,:,1) = reshape(defaultMap(ceil(map_normalize*64),1)*256,size(map_normalize));
tmap(:,:,2) = reshape(defaultMap(ceil(map_normalize*64),2)*256,size(map_normalize));
tmap(:,:,3) = reshape(defaultMap(ceil(map_normalize*64),3)*256,size(map_normalize));

mask_index = find(double(nan_mask)==1);
tmap(mask_index+numel(map_nothre_reshape)*0) = 0;
tmap(mask_index+numel(map_nothre_reshape)*1) = 0;
tmap(mask_index+numel(map_nothre_reshape)*2) = 0;


tmap = flip(tmap,1);
Img_RGB = tmap;
end



