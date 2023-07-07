clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

for ix=1:3
    ICAname = fullfile(WholePath,'Harris',...
        ['ica_all_05_icasso_iter_1000_comp_',num2str(ix-1),'_mask_z_1_allen_masked_sym_thresh.nii.gz']);
    X = spm_read_vols(spm_vol(ICAname));
    X0(:,:,:,ix) = imresize3(X,[114 80 38],'nearest');
end


Type = {'chat';'pv';'som';'vglut2';};
ProjType = {'chat';'pv';'sst';'vglut'};

WsT=[];
for tloop=1:numel(Type)

Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'spmT_0001.nii');
Ya = spm_read_vols(spm_vol(Func3D));


filename = fullfile(WholePath,'BF_structure_connectivity',...
    ['nR_',ProjType{tloop},'_with_mask.nii']);
Func_Img_3D = spm_read_vols(spm_vol(filename));
Yp = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;

F = figure('Position', [680 387 307 591]);
subplot(2,1,1)
for ic=1:3

    lmask = X0(:,:,:,ic);
    Xa = fmask(Ya,lmask);
    bar(ic*1.2,mean(Xa)); hold on;
    errorbar(ic*1.2,mean(Xa),std(Xa)/sqrt(sqrt(numel(Xa))));

end

subplot(2,1,2)
for ic=1:3
    lmask = X0(:,:,:,ic);%X==ic;
    Xp = fmask(Yp,lmask);
    Xp(Xp>100)=100;
    bar(ic*1.2,mean(Xp)); hold on;
    errorbar(ic*1.2,mean(Xp),std(Xp)/sqrt(sqrt(numel(Xp))));
end
% saveas(F,fullfile('F:\BF_optogentics\SubMit\V3\add',['ICA_',ProjType{tloop},'.emf']));
% close all
end




Type = {'chat';'pv';'som';'vglut2';};
ProjType = {'chat';'pv';'sst';'vglut'};

WsT=[];
for tloop=1:numel(Type)

Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
    upper(Type{tloop}),'spmT_0001.nii');
Ya = spm_read_vols(spm_vol(Func3D));


filename = fullfile(WholePath,'BF_structure_connectivity',...
    ['nR_',ProjType{tloop},'_with_mask.nii']);
Func_Img_3D = spm_read_vols(spm_vol(filename));
Yp = ( Func_Img_3D + flip(Func_Img_3D,1) )/2;

F = figure('Position', [680 387 307 591]);
subplot(2,1,1)
for ic=1:3

    lmask = X0(:,:,:,ic);
    
    Xa = fmask(Ya,lmask);
    bar(ic*1.2,mean(Xa)); hold on;
    errorbar(ic*1.2,mean(Xa),std(Xa)/sqrt(sqrt(numel(Xa))));

end

subplot(2,1,2)
for ic=1:3
    lmask = X0(:,:,:,ic);%X==ic;
    Xp = fmask(Yp,lmask);
    Xp(Xp>100)=100;
    bar(ic*1.2,mean(Xp)); hold on;
    errorbar(ic*1.2,mean(Xp),std(Xp)/sqrt(sqrt(numel(Xp))));
end
% saveas(F,fullfile('F:\BF_optogentics\SubMit\V3\add',['ICA_',ProjType{tloop},'.emf']));
% close all
end


