clc;clear
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

Excel = fullfile(codepath,'Exp_recording.xlsx');
[~,~,CellData] = xlsread(Excel,'RawData');
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

NII_v213 = fullfile(codepath,'Label_Mouse_213_v38.nii');
ihdr = spm_vol(NII_v213);
Labels_v38 = spm_read_vols(ihdr);
lmask = spm_read_vols(spm_vol(fullfile(codepath,'lmask_Mouse_v38.nii')));

%% regressor choices
Reg_choices = {'rp';'rp"';'10PCs'};


CellType = {'CHAT';'PV';'SOM';'Vglut2';'Ctrl'};
Duration = {'05';'2'};




%
for idx = 4:4:51
    idx
    
    path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
    RARE =  ExpTable.T2(idx);
    EPI = [ExpTable.A(idx),ExpTable.B(idx),ExpTable.C(idx),...
        ExpTable.D(idx),ExpTable.E(idx),ExpTable.F(idx)];
    GEEPI = EPI(~isnan(EPI));
    EPItp = [ExpTable.TOPUP1(idx),ExpTable.TOPUP2(idx)];

    for kk = GEEPI
        
        GEPath = fullfile(path,'Results',num2str(kk));
        filename = fullfile(GEPath,'snrmUw2dseq.nii');
%         hdr = spm_vol(filename);
%         fMRI_4D = spm_read_vols_4D(hdr);
%         
%         reg = load(fullfile(GEPath,'Multi_Regessor.txt'));
%         RegBasFuc = [ones(length(reg),1),reg];
%         
%         data_ready_regress = fmask(fMRI_4D,lmask);
%         data_ready_regress(isnan(data_ready_regress)) = 0;
%         for iii = 1:size(data_ready_regress,1)
%             [Beta,~,Residual] = regress(squeeze(data_ready_regress(iii,:))',RegBasFuc);
%             data_ready_regress(iii,:) = Residual + Beta(1);
%         end
%         data_regressout = funmask(data_ready_regress,lmask);
%         for h=1:numel(hdr);hdr(h).fname=fullfile(GEPath,'rsnrmUw2dseq.nii');end
%         spm_write_vol_4D(hdr,data_regressout);
%         clear data_* fMRI_4D
        
        filename = fullfile(GEPath,'R_snrmUw2dseq.nii');
        hdr = spm_vol(filename);
        fMRI_4D = spm_read_vols_4D(hdr);
        
        Smean = nan([213,numel(hdr)]);
        Smedian = nan([213,numel(hdr)]);
        for loop=1:213
            Xmask = Labels_v38==loop;
            Smean(loop,:) = nanmean(fmask(fMRI_4D,Xmask),1);
            Smedian(loop,:) = nanmedian(fmask(fMRI_4D,Xmask),1);
        end
        Excel = fullfile(GEPath,'TimeSeries.xlsx');
        xlswrite(Excel,Smean,'Mean','A1');
        xlswrite(Excel,Smedian,'Median','A1');
        clear fMRI_4D
    end
end