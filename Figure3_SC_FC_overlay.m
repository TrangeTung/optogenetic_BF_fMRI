clc;clear
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
WholePath = 'F:\BF_optogentics\\';




Excel = fullfile(codepath,'Label_213_v38.xlsx');
[~,~,CellData] = xlsread(Excel);
ExpTable_213 = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels = spm_read_vols(ihdr);

ProjectType = {'chat';'vglut';'pv';'sst'};
Type = {'CHAT';'PV';'SOM';'Vglut2';};
for tloop = 1:numel(Type)
    
    blckshift=160;
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',ProjectType{tloop},'_with_mask.nii']);
    Func_Img_3D = spm_read_vols(spm_vol(filename));

    Trace = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;
    Trace = Trace/max(Trace(:))*100;
    Trace = flip(Trace,2);

    
    
    I = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii')));
    I0x = imresize3(Trace,size(I),'linear');
    I0x = flip(flip(I0x,2),3);
    I0x = smooth3(I0x,'box',5);
    Ip = permute(I,[3 2 1]);
    I0xp = permute(I0x,[3 2 1]);
    lmask = Ip<10^5;
    barlim = [1 50];
    IMG=[];
    for sl=119:30:200
        Func_Img_3D = I0xp;
        bar_value = barlim;
        slice=sl-5;
        map_nothre_reshape = reshape(permute(flip(Func_Img_3D(:,:,slice),1),[2 1 3]),[size(Func_Img_3D,2) size(Func_Img_3D,1)*numel(slice)]);
        MyMap = [(64:-1:1)',64*ones(64,1),(64:-1:1)']/64;
        defaultMap =  MyMap;
        nothrebar = [min(bar_value(:)) max(bar_value(:))];
        tmap = ones([size(map_nothre_reshape),3]);
        map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
        map_normalize(map_normalize<=0) = 0.000000001;
        map_normalize(map_normalize>=1) = 1;
        nan_mask = isnan(map_normalize);
        map_normalize(isnan(map_normalize)) = 0.000000001;
        tmap(:,:,1) = reshape(defaultMap(ceil(map_normalize*64),1)*255,size(map_normalize));
        tmap(:,:,2) = reshape(defaultMap(ceil(map_normalize*64),2)*255,size(map_normalize));
        tmap(:,:,3) = reshape(defaultMap(ceil(map_normalize*64),3)*255,size(map_normalize));
        
        mask_index = find(double(nan_mask)==1);
        tmap(mask_index+numel(map_nothre_reshape)*0) = 0;
        tmap(mask_index+numel(map_nothre_reshape)*1) = 0;
        tmap(mask_index+numel(map_nothre_reshape)*2) = 0;
        Img_RGB = tmap;
        
        Img_RGB = flip(Img_RGB,1);
        Lx = flip(lmask(:,:,sl)',2);
        SE = strel('disk',2); Lx = imdilate(Lx,SE);
        Lx = 1-Lx;
        Lx3 = repmat(Lx,[1 1 3]);
        Img_RGB = Img_RGB.*double(Lx3);
        Img_RGB(Lx3~=1)=255;
        
        if sl==119+30*2;Img_RGB(:,end+(1:blckshift),:)=255;end
        if sl==119+30*1;Img_RGB=cat(2,255*ones(size(Img_RGB,1),blckshift/2,size(Img_RGB,3)),Img_RGB);Img_RGB(:,end+(1:blckshift/2),:)=255;end
        if sl==119+30*0;Img_RGB=cat(2,255*ones(size(Img_RGB,1),blckshift,size(Img_RGB,3)),Img_RGB);end
        
        IMG = cat(1,IMG,Img_RGB);
    end
    IMG(:,end-4:end,:)=255;
    IS = IMG;
    
    
    
    
    filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'spmT_0001.nii');
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    Func_Img_3D = Func_Img_3D/2+flip(Func_Img_3D,1)/2;

    bar=[-10 -2.3 2.3 10];
    Func_Img_3D(Func_Img_3D>bar(2)&Func_Img_3D<bar(3))=0;
    Func_Img_3D = flip(Func_Img_3D,2);
    
    I = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii')));
    I0x = imresize3(Func_Img_3D,size(I),'linear');
    I0x = flip(flip(I0x,2),3);
    I0x = smooth3(I0x,'box',1);
    Ip = permute(I,[3 2 1]);
    I0xp = permute(I0x,[3 2 1]);
    lmask = Ip<10^5;
    barlim = [-10 10]*0.8;
    IMG=[];
    for sl=119:30:200
        Func_Img_3D = I0xp;
        bar_value = barlim;
        slice=sl-5;
        map_nothre_reshape = reshape(permute(flip(Func_Img_3D(:,:,slice),1),[2 1 3]),[size(Func_Img_3D,2) size(Func_Img_3D,1)*numel(slice)]);
        MyMap1 = [64*ones(64,1),(64:-1:1)',(64:-1:1)']/64;
        MyMap2 = [(64:-1:1)',(64:-1:1)',64*ones(64,1)]/64;
        defaultMap =  [MyMap2(end:-2:1,:); MyMap1(1:2:end,:)];
        nothrebar = [min(bar_value(:)) max(bar_value(:))];
        tmap = ones([size(map_nothre_reshape),3]);
        map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
        map_normalize(map_normalize<=0) = 0.000000001;
        map_normalize(map_normalize>=1) = 1;
        nan_mask = isnan(map_normalize);
        map_normalize(isnan(map_normalize)) = 0.000000001;
        tmap(:,:,1) = reshape(defaultMap(ceil(map_normalize*64),1)*255,size(map_normalize));
        tmap(:,:,2) = reshape(defaultMap(ceil(map_normalize*64),2)*255,size(map_normalize));
        tmap(:,:,3) = reshape(defaultMap(ceil(map_normalize*64),3)*255,size(map_normalize));
        
        mask_index = find(double(nan_mask)==1);
        tmap(mask_index+numel(map_nothre_reshape)*0) = 0;
        tmap(mask_index+numel(map_nothre_reshape)*1) = 0;
        tmap(mask_index+numel(map_nothre_reshape)*2) = 0;
        Img_RGB = tmap;
        
        Img_RGB = flip(Img_RGB,1);
        Lx = flip(lmask(:,:,sl)',2);
        SE = strel('disk',2); Lx = imdilate(Lx,SE);
        Lx = 1-Lx;
        Lx3 = repmat(Lx,[1 1 3]);
        Img_RGB = Img_RGB.*double(Lx3);
        Img_RGB(Lx3~=1)=255;
        
        if sl==119+30*2;Img_RGB(:,end+(1:blckshift),:)=255;end
        if sl==119+30*1;Img_RGB=cat(2,255*ones(size(Img_RGB,1),blckshift/2,size(Img_RGB,3)),Img_RGB);Img_RGB(:,end+(1:blckshift/2),:)=255;end
        if sl==119+30*0;Img_RGB=cat(2,255*ones(size(Img_RGB,1),blckshift,size(Img_RGB,3)),Img_RGB);end
        
        IMG = cat(1,IMG,Img_RGB);
    end
    IMG(:,end-4:end,:)=255;
    IF = IMG;
    
    
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_OVERLAY.tiff']);
    imwrite(uint8((IS+IF)/2),filename);
    close all;
    
    
    filename = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'spmT_0001.nii');
    Func_Img_3D = spm_read_vols(spm_vol(filename));
    X = Trace;
    Y = Func_Img_3D;
    lmask = Labels~=0;
    
    NumSC = numel(find(fmask(X,lmask)>2));
    NumfP = numel(find(fmask(Y,lmask)>bar(3)));
    NumfN = numel(find(fmask(Y,lmask)<bar(2)));
    N_SC_fP = numel(find(fmask(Y,lmask)>bar(3) & fmask(X,lmask)>2));
    N_SC_fN = numel(find(fmask(Y,lmask)<bar(2) & fmask(X,lmask)>2));
    N_null = numel(find( fmask(lmask,lmask)==1 & ~(fmask(Y,lmask)<bar(2)) & ~(fmask(X,lmask)>2) & ~(fmask(Y,lmask)<bar(2))));
    A = [NumfP NumSC NumfN]; I = [N_SC_fP  0 N_SC_fN 0];
    F = figure;
    venn(A,I,'ErrMinMode','None','FaceAlpha', 0.6);
    axis image
    filename = fullfile(WholePath,'BF_structure_connectivity',...
        ['nR_',Type{tloop},'_veen.emf']);
    xlim([-300 300]);ylim([-300 300]);
    saveas(F,filename);
    A
    I

end