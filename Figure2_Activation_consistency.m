clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));


Mask3D = [codepath,'\lmask_Mouse_v38.nii'];
lmask=spm_read_vols(spm_vol(Mask3D));

Type = {'chat';'pv';'som';'vglut2';};
clear Xmask;
for tloop=1:numel(Type)
    
    Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
        upper(Type{tloop}),'spmT_0001.nii');
    Ya = spm_read_vols(spm_vol(Func3D));

    M = Ya>=1.0 | Ya<=-2.0;
    M(~lmask)=0;
    if tloop==1; Xmask=M*0;end
    Xmask=Xmask+double(M);
end


Img_RGB = MY_display_function_map_3D_nothreshold(Xmask,[0 5],9:2:35);
Img_RGB = cat(1,Img_RGB,zeros([10,size(Img_RGB,2),3])+255);

overlay = Img_RGB;
slice = 1:numel([9:2:35]);
[SizeX,SizeY,~] = size(overlay);
slice_per_row = 5;
pump = ceil(length(slice)/slice_per_row);
blank = pump*slice_per_row-length(slice);
overlay(:,SizeY+1:SizeY+SizeY/length(slice)*blank,:) = repmat(overlay(2,2,:),[SizeX,SizeY/length(slice)*blank,1]);
for j =1:pump
    overlay2((j-1)*SizeX+1:j*SizeX,:,:) = overlay(:,(j-1)*(size(overlay,2)/pump)+1:j*(size(overlay,2)/pump),:);
end
Img_RGB = overlay2;

imwrite(uint8(Img_RGB),fullfile(WholePath,'2nd_TwoSampleTTest_results','activation_overlay.tiff'));



% In DMN and Out DMN
Func3D = fullfile(WholePath,'Zerbi','ICA_zerbi15.nii');
ICA = spm_read_vols(spm_vol(Func3D));
I = ICA(:,:,:,3)>3 & double(lmask)==1 ;
for lp=[4 3 2 1 0]
    X= Xmask==lp & double(lmask)==1;
    Xin = numel(find(X&I==1));
    Xout = numel(find( X&~I==1));
    perc(lp+1) = Xin/(Xin+Xout);
    perc(lp+1) = Xin/numel(find(I==1));
end
perc



%% Activation in DMN
CellType = {'Vglut2';'CHAT';'PV';'SOM';'Ctrl'};
Duration = {'05';'2'};

for cl = 1:numel(CellType)
    
    
    clear matrix matrix1 matrix2 scans
    fclose('all');
    spm('defaults', 'FMRI');
    set(spm('CreateIntWin','off'),'Visible','on');
    
    dest = fullfile(WholePath,'2nd_TwoSampleTTest_results',CellType{cl});
    
    %%%%%%%%%%%%    Exp    %%%%%%%%%%%%
    Excel = fullfile(codepath,'Exp_recording.xlsx');
    sheet = [CellType{cl},'_2nd'];
    [~,~,CellData] = xlsread(Excel,sheet);
    ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
    subj_num = size(ExpTable,1);
    clear scans1
    a=0;
    for idx = 1:size(ExpTable)
        path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
        GEEPI = [ExpTable.A_05(idx),ExpTable.B_05(idx),ExpTable.C_05(idx),...
            ExpTable.A_2(idx),ExpTable.B_2(idx),ExpTable.C_2(idx)];
        for gl=1:numel(GEEPI)
            a=a+1;
            scans1(a,1) = {fullfile(path,'Functions\tsfMRI\',num2str(GEEPI(gl)),'\spmT_0001.nii')};
        end
    end
    
    for al=1:a
        I = spm_read_vols(spm_vol(scans1{al}));
        if al==1;Iall=zeros([size(I),a]);end
        Iall(:,:,:,al)=I;
    end
    
    
    for iclp=1:6
        ICAmask = ICA(:,:,:,iclp)>3 & double(lmask)==1 ;
        V = fmask(Iall,ICAmask);
        Vm = nanmedian(V,1);
        eval(['V_',CellType{cl},'_ICA_',num2str(iclp),'=Vm(:);']);
    end
end

ICAname = {'Somatosensory';'Sensory(VC+AUD+Mot)';'DMN';'HPF';'BasalGanglia';'Olfactory'};
Excel = fullfile(WholePath,'Zerbi','NetworkActivation.xlsx');
for iclp=1:numel(ICAname)
    TxT = CellType(:)';
    xlswrite(Excel,TxT,ICAname{iclp},'A1');
    cl=1;  eval(['Vm=V_',CellType{cl},'_ICA_',num2str(iclp),';']);
    xlswrite(Excel,Vm,ICAname{iclp},'A2');
    cl=2;  eval(['Vm=V_',CellType{cl},'_ICA_',num2str(iclp),';']);
    xlswrite(Excel,Vm,ICAname{iclp},'B2');
    cl=3;  eval(['Vm=V_',CellType{cl},'_ICA_',num2str(iclp),';']);
    xlswrite(Excel,Vm,ICAname{iclp},'C2');
    cl=4;  eval(['Vm=V_',CellType{cl},'_ICA_',num2str(iclp),';']);
    xlswrite(Excel,Vm,ICAname{iclp},'D2');
    cl=5;  eval(['Vm=V_',CellType{cl},'_ICA_',num2str(iclp),';']);
    xlswrite(Excel,Vm,ICAname{iclp},'E2');
end


