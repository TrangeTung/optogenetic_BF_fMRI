clc;clear
WholePath = 'F:\BF_optogentics\\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

SC_Excel = fullfile(codepath,'Allen_mouse_SC_213.xlsx');
SC = xlsread(SC_Excel);

NII_v213 = fullfile(codepath,'Colormap_3Dviewer','R_Label_Mouse_213.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);

Excel = fullfile(codepath,'Label_213_v38.xlsx');
[~,~,CellData] = xlsread(Excel);
ExpTable_213 = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels_v38 = spm_read_vols(ihdr);
%% pipeline figure
%{
dest = fullfile(WholePath,'DimReduc_SC');
Type = {'chat';'vglut';'pv';'sst'};
for tloop=1:numel(Type)
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    Vs = zeros(213,1);
    VoxNum = zeros(213,1);
    for loop=1:213
        lmask = Labels_v38==loop;
        VoxNum(loop) = numel(find(lmask==1));
        Vs(loop) = nanmedian(fmask(Trace,lmask))/255*100;
    end
    ResizeVox = (3*3*9.47)/(2*2*8);
    
    VoxThr = 0;

    VB = VoxNum/ResizeVox>0&Vs~=0; % 1st
    SCx = SC;SCx(~VB,:)=nan;
    
    
    A = repmat(Vs,1,10);
    MyMap = [ones(1,64);(63:-1:0)/64;(63:-1:0)/64]';
    Img_1st = MY_display_color_nothreshold(log(A),[0 .5],MyMap);
    MyMap = [(63:-1:0)/64;(63:-1:0)/64;ones(1,64)]';
    Img_2nd = MY_display_color_nothreshold(SCx,[0 .10],MyMap);
    
    FF = cat(2,uint8(Img_1st),uint8(Img_2nd));
    FF = imresize3(FF,[size(FF,1)*5 size(FF,2)*5 3],'nearest');
    imwrite(uint8(FF),fullfile(dest,['Projection_',Type{tloop},'.tif']));
    
    
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
    
    
    %  1st
    Lv = Vs>0; Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask1 = smooth3(Rmask,'box',5);
    p = patch(isosurface(Rmask1,0.5));
    isonormals(Rmask1,p);
    p.FaceColor = [1,0,0];p.EdgeColor = 'none';
    p.FaceAlpha = 0.10;
    axis equal; axis tight; axis off;
    
    %  2nd
    Lv = (nanmean(SCx,1)>mean(nanmean(SCx,1)))'&Vs<=0; Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask2 = smooth3(Rmask,'box',5);
    p = patch(isosurface(Rmask2,0.5));
    isonormals(Rmask2,p);
    p.FaceColor = [1,0,1];p.EdgeColor = 'none';
    p.FaceAlpha = 0.05;
    axis equal; axis tight; axis off;
    
    Lv = (nanmean(SCx,1)>mean(nanmean(SCx,1)))'&Vs>0; Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask2 = smooth3(Rmask,'box',5);
    p = patch(isosurface(Rmask2,0.5));
    isonormals(Rmask2,p);
    p.FaceColor = [1,0,1];p.EdgeColor = 'none';
    p.FaceAlpha = 0.05;
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
    
end



function Img_RGB = MY_display_color_nothreshold(V2D,bar_value,MyMap)

V2D(isnan(V2D))=0;
map_nothre_reshape=V2D;
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

Img_RGB = tmap;
end

%}


%% CC b.w. SC& Func

Excel = fullfile(codepath,'Label_213_v38.xlsx');
[~,~,CellData] = xlsread(Excel);
ExpTable_213 = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Colormap_3Dviewer','R_Label_Mouse_213.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);
NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels_v38 = spm_read_vols(ihdr);
lmask = spm_read_vols(spm_vol(fullfile(codepath,'lmask_Mouse_v38.nii')));

Type = {'chat';'pv';'sst';'vglut'};
for tloop=1:numel(Type)
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Func_Img_3D(isnan(Func_Img_3D))=0;
    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    V  = fmask(Trace,lmask);
    
    
    if tloop==1;Vall = zeros([size(V,1),numel(Type)*2]);end
    Vall(:,tloop*2-1) = V;
    
    filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'beta_0001.nii');
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Func_Img_3D(isnan(Func_Img_3D))=0;
    V  = fmask(Func_Img_3D,lmask);
    Vall(:,tloop*2) = V;
end
CC = corr(Vall(:,[1:2:end,2:2:end]));
CCv = flip(tril(CC,0),2);


F = figure;
for xl=1:size(CCv,1)
    for yl=1:size(CCv,2)
        X = xl; Y=10-yl;
        Radius = CCv(xl,yl);
        if Radius~=1&Radius~=0
        x = X+cos(0:pi/20:2*pi)*Radius*.6;
        y = Y+sin(0:pi/20:2*pi)*Radius*.6;
        patch(x,y,ones(size(y)),'edgecolor','none',...
            'facecolor',[.5 .5 .5],'facealpha',min([Radius+0.3,1]));
        hold on;
        text(X,Y,num2str(Radius,2),'Horizontalalignment','center')
        end
        if Radius~=0
        patch([X-.5 X-.5 X+.5 X+.5],[Y-.5 Y+.5 Y+.5 Y-.5],[1 1 1 1],...
            'edgecolor',[.8 .8 .8],'facecolor','none');
        
        
        end
    end
end
axis equal; axis off



%% NMF dimension reduction
%
dest = fullfile(WholePath,'DimReduc_SC');
Type = {'chat';'vglut';'pv';'sst'};
DMx = [];
for tloop=1:numel(Type)
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    Vs = zeros(213,1);
    VoxNum = zeros(213,1);
    for loop=1:213
        lmask = Labels_v38==loop;
        VoxNum(loop) = numel(find(lmask==1));
        Vs(loop) = nanmedian(fmask(Trace,lmask))/255*100;
    end
    ResizeVox = (3*3*9.47)/(2*2*8);
    
    VoxThr = 0;
    
    VB = VoxNum/ResizeVox>0&Vs~=0; % 1st
    SCx = SC(VB,:)';
    
    
    DM = cat(2,Vs/max(Vs),SCx);
    DMx = [DMx,DM];
end
DMx(isnan(DMx))=0;
DM=DMx;


1;
%{
R2 = nan(64,1000);
Wall = cell(64,1000);
Hall = cell(64,1000);
DM(DM>.3)=.3;
for st=1:64
    
    fprintf(['BF No. ',num2str(st),'\n']);
    parfor sh=1:1000
        [W,H,errs,loss] = nmf_euc(DM, st);
        ySC = W*H;
        
        Vec = DMx(:)'; VecP = ySC(:)';
        mu = mean(VecP,2).*mean(Vec,2) - mean(Vec.*VecP,2);
        md = mean(VecP,2).^2 - mean(VecP.^2,2);
        b = mean(Vec,2)-mu./md.*mean(VecP,2);
        Vecf = mu./md.*VecP + b;
        SSR = sum((Vecf-mean(Vec,2)).^2,2);
        SST = sum((Vec-mean(Vec,2)).^2,2);
        R2(st,sh) = SSR./SST;
        Wall{st,sh}=W;
        Hall{st,sh}=H;
    end
end
cd(fullfile(WholePath,'DimReduc_SC'));
save('Celltype_R2_64W_1000times_Oct25.mat','R2');
save('Celltype_W_64W_1000times_Oct25.mat','Wall','-v7.3');
save('Celltype_H_64W_1000times_Oct25.mat','Hall','-v7.3');

%}

%% Vector align

cd(fullfile(WholePath,'DimReduc_SC'));
load('Celltype_R2_64W_1000times_Oct25.mat');
load('Celltype_W_64W_1000times_Oct25.mat');
load('Celltype_H_64W_1000times_Oct25.mat');
% Explained variance
%{
R = R2(1:64,:)./R2(64,:);R(2,:)=R(2,:)+0.02;R(3,:)=R(3,:)+0.02;
%R = cellfun(@(x)sum(x),R2aligned);
F = figure;
semilogx(1:size(R,1),R,'color',[.8 .8 .8]);
hold on;plot(1:size(R,1),mean(R,2),'k-o','linewidth',1)
xlabel('No. of vectors');ylim([0 1]);xlim([1 100])
ylabel('Explained variance');
set(gca,'fontsize',13,'linewidth',.5,'box','off','tickdir','out');
saveas(F,'Explained_variance.emf');
%}

%
% Tree
Waligned = cell(64,1);
Haligned = cell(64,1);
R2aligned = cell(64,1);
for w=1:64
    Wx = Wall(w,:); Hx = Hall(w,:);
    R2collec = zeros(w,1000);
    for loop=1:1000
        Wxx = Wx{loop};
        Hxx = Hx{loop};
        
        R2 = zeros([w,1]);
        for xl=1:w
            ySC = Wxx(:,xl)*Hxx(xl,:);
            
            Vec = DM(:)'; VecP = ySC(:)';
            mu = mean(VecP,2).*mean(Vec,2) - mean(Vec.*VecP,2);
            md = mean(VecP,2).^2 - mean(VecP.^2,2);
            b = mean(Vec,2)-mu./md.*mean(VecP,2);
            Vecf = mu./md.*VecP + b;
            SSR = sum((Vecf-mean(Vec,2)).^2,2);
            SST = sum((Vec-mean(Vec,2)).^2,2);
            R2(xl) = SSR./SST;
        end
        [sR2,s] = sort(R2,'descend');
        R2collec(:,loop) = sR2;
        
        Wxy = Wxx(:,s);
        Wx{loop} = Wxy;
        Hxy = Hxx(s,:);
        Hx{loop} = Hxy;
    end
    Wr = cell2mat(Wx);
    Wmean = zeros(size(Wxx));
    for xl=1:w
        Wmean(:,xl) = mean(Wr(:,xl:w:end),2);
    end
    Waligned{w} = Wmean;
    Hr = cat(1,Hx{:});
    Hmean = zeros(size(Hxx));
    for xl=1:w
        Hmean(xl,:) = mean(Hr(xl:w:end,:),1);
    end
    Haligned{w} = Hmean;    
    R2aligned{w} = mean(R2collec,2);
end

cd(fullfile(WholePath,'DimReduc_SC'));
save('Celltype_Waligned_Oct25.mat','Waligned');
save('Celltype_Haligned_Oct25.mat','Haligned');
save('Celltype_R2aligned_Oct25.mat','R2aligned');
%}
%% R2 matrix
%{
F = figure;
for w = 1:64
    Ww = Waligned{w};
    Rw = R2aligned{w};
    
    for xl = 1:numel(Rw)
        X = 10-xl; Y=10-w;
        Radius = Rw(xl);
        
        x = X+cos(0:pi/20:2*pi)*Radius*2;
        y = Y+sin(0:pi/20:2*pi)*Radius*2;
        patch(x,y,ones(size(y)),'edgecolor','none',...
            'facecolor',[1 0 0],'facealpha',min([Radius+0.3,1]));
        hold on;
        patch([X-.5 X-.5 X+.5 X+.5],[Y-.5 Y+.5 Y+.5 Y-.5],[1 1 1 1],...
            'edgecolor',[.8 .8 .8],'facecolor','none');
        text(X,Y,num2str(Radius,2),'Horizontalalignment','center')

    end
    
    %{
    if w~=1;
        CC = pdist2(Wold',Ww');
        %[CC,p] = corr(Wold,Ww);
        for xl = 1:size(CC,1)
            for yl = 1:size(CC,2)
                Xx = 10-xl; Yx=10+1-w;
                Xy = 10-yl; Yy=10-w;
                
                plot([Xx,Xy],[Yx,Yy],'color',[.9 .9 .9],'linewidth',2-CC(xl,yl));
                
            end
        end
    end
    %}
    Wold = Ww;
end
axis equal; axis off;
cd(fullfile(WholePath,'DimReduc_SC'));
saveas(F,fullfile('R2_eachVector.emf'));
%}


%% 3D map
NII_v213 = fullfile(codepath,'Colormap_3Dviewer','R_Label_Mouse_213.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);

dest = fullfile(WholePath,'DimReduc_SC');
load(fullfile(dest,'Celltype_Waligned_Oct25.mat'))
w=5;Ww = Waligned{w};

ColorVec = {[0 0 1];[1 0 0];[1 1 0];}; % Cerebrum % Brain stem % Cerebellum
ColorBin = {[1:87];[88:201];[202:213];};

% exclude the distortion regions
Cimg = spm_read_vols(spm_vol(fullfile(codepath,'DistortionMask.nii')));
SVox = zeros(213,1);
for il=1:213
    X = fmask(Cimg,Labels_v38==il);
    SVox(il) = mean(X);
end
SB = SVox<.4;

VoxNum = zeros(213,1);
for loop=1:213
    lmask = Labels_v38==loop;
    VoxNum(loop) = numel(find(lmask==1));
end
ResizeVox = (3*3*9.47)/(2*2*8);
VoxThr = 10;
VB = VoxNum/ResizeVox>VoxThr; % 1st


for wl=1:w

    Lv = Ww(:,wl); [s,~] = sort(Lv,'descend');
    Lv(Lv<s(42))=0;
    Lv(SB)=0; Lv(~VB) = 0;
    n=0;
    F = figure('Position', [73 506 673 187]);
    for cl=1:numel(ColorBin)
       x = Lv(ColorBin{cl});
       ROIs = ExpTable_213.abbre(ColorBin{cl});
       [sx,s] = sort(x,'descend');
       n(cl+1)=numel(find(sx~=0));
       b=bar(sum(n(1:cl))+(1:n(cl+1)),sx(sx~=0),'edgecolor','none',...
          'facecolor',ColorVec{cl}); hold on;
       ROIv = ROIs(s(sx~=0))
       text(b(1).XEndPoints,b(1).YEndPoints,ROIv,'rotation',90)
    end
    ylim([0 0.5])
    set(gca,'box','off','linewidth',.5,'tickdir','out');
    saveas(F,fullfile(dest,['Bar_',num2str(wl),'.emf']));
    
    

    F = figure('Color',[1 1 1],'Position', [680 52 799 926]);

    ALL = Labels~=0; ALL = smooth3(ALL,'box',11);
    p0 = patch(isosurface(ALL,0.5));
    isonormals(ALL,p0);
    p0.FaceColor = 'black';
    p0.EdgeColor = 'none';
    p0.FaceAlpha=0.05;
    axis equal; axis tight
    axis off;
    hold on;


    Lv = Ww(:,wl)>mean(Ww(:,wl)); Rmask=Labels*0;
    for rl = 1:213;Rmask(Labels==rl)=Lv(rl);end
    Rmask1 = smooth3(Rmask,'box',5);
    p = patch(isosurface(Rmask1,0.5));
    isonormals(Rmask1,p);
    p.FaceColor = [1,0,0];p.EdgeColor = 'none';
    p.FaceAlpha = 0.20;
    axis equal; axis tight; axis off;

    view([-90 180]);
    saveas(F,fullfile(dest,['ColorMap_NMF_3DViewer_',num2str(wl),'_dorsal.tif']));
    view([0 180]);
    saveas(F,fullfile(dest,['ColorMap_NMF_3DViewer_',num2str(wl),'_lateral.tif']));
    
    A = imread(fullfile(dest,['ColorMap_NMF_3DViewer_',num2str(wl),'_dorsal.tif']));
    B = imread(fullfile(dest,['ColorMap_NMF_3DViewer_',num2str(wl),'_lateral.tif']));
    C = cat(2,A(1:end,161:1200,:),B(1:end,231:1030,:));
    imwrite(C,fullfile(dest,['ColorMap_NMF_3DViewer_',num2str(wl),'_merge.tif']));
    close all

    
    
    
    
end





%% cca

Excel = fullfile(codepath,'Label_213_v38.xlsx');
[~,~,CellData] = xlsread(Excel);
ExpTable_213 = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Colormap_3Dviewer','R_Label_Mouse_213.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);
NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels_v38 = spm_read_vols(ihdr);

dest = fullfile(WholePath,'DimReduc_SC');
load(fullfile(dest,'Celltype_Waligned.mat'))
w=3;Ww = Waligned{w};

Type = {'chat';'pv';'sst';'vglut';};
for tloop=1:numel(Type)
    
    VoxNum = zeros(213,1);
    for loop=1:213
        lmask = Labels==loop;
        VoxNum(loop) = numel(find(lmask==1));
    end
    ResizeVox = (3*3*9.47)/(2*2*8);
    VoxThr = 400;
    Vbi = VoxNum/ResizeVox>VoxThr;
    
    
    
    filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'beta_0001.nii');
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    
    Vf = zeros(213,1);
    VoxNum = zeros(213,1);
    for loop=1:213
        lmask = Labels_v38==loop;
        Vf(loop) = nanmedian(fmask(Func_Img_3D,lmask));
    end
    Vf(isnan(Vf))=0;
    
    [A,B,r,U,V,stats] = canoncorr(log(100*Ww(Vbi,:)+1),Vf(Vbi,:));
    r
    A
    
end


