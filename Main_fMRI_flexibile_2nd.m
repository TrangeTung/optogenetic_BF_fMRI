%% preparation
% According to your experimental parameters and certain file direction, Change these bellow.
clc;clear
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

CellType = {'PV';'SOM';'CHAT';'Vglut2';'Ctrl'};
Duration = {'05';'2'};

for cl = 1:4%numel(CellType)
        
        clear matrix matrix1 matrix2 scans
        fclose('all');
        spm('defaults', 'FMRI');
        set(spm('CreateIntWin','off'),'Visible','on');
        
        %%%%%%%%exp group%%%%%%%%%%%%%%%%
        Excel = fullfile(codepath,'Exp_recording.xlsx');
        sheet = [CellType{cl},'_2nd'];
        [~,~,CellData] = xlsread(Excel,sheet);
        ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
        
        dest = fullfile(WholePath,'2nd_results',sheet);
        rmdir(dest,'s');
        if ~exist(dest,'dir');mkdir(dest);end
        
        
        subj_num = size(ExpTable,1);
        clear scans1
        for idx = 1:size(ExpTable)
            path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
            eval(['a=ExpTable.A_',Duration{1},'(',num2str(idx),');']);
            eval(['b=ExpTable.B_',Duration{1},'(',num2str(idx),');']);
            eval(['c=ExpTable.C_',Duration{1},'(',num2str(idx),');']);
            eval(['d=ExpTable.A_',Duration{2},'(',num2str(idx),');']);
            eval(['e=ExpTable.B_',Duration{2},'(',num2str(idx),');']);
            eval(['f=ExpTable.C_',Duration{2},'(',num2str(idx),');']);
            GEEPI = [a,b,c,d,e,f];
            for gl=1:numel(GEEPI)
                scans1((idx-1)*numel(GEEPI)+gl,1) = {fullfile(path,'Functions','tsfMRI',num2str(GEEPI(gl)),'con_0001.nii')};
            end
        end
        session_num = numel(GEEPI);
        a=session_num;
        n=subj_num;
        for m=1:n
            matrix1((a*m-a+1):a*m,1)=ones(a,1);
            matrix1((a*m-a+1):a*m,2)=ones(a,1);
            matrix1((a*m-a+1):a*m,3)=m*ones(a,1);
            matrix1((a*m-a+1):a*m,4)= [1:a]';
        end
        l=n;

        
        %%%%%%%%CTRL group%%%%%%%%%%%%%%%%
        Excel = fullfile(codepath,'Exp_recording.xlsx');
        sheet = [CellType{5},'_2nd'];
        [~,~,CellData] = xlsread(Excel,sheet);
        ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

        subj_num = size(ExpTable,1);
        clear scans2
        for idx = 1:size(ExpTable)
            path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
            eval(['a=ExpTable.A_',Duration{1},'(',num2str(idx),');']);
            eval(['b=ExpTable.B_',Duration{1},'(',num2str(idx),');']);
            eval(['c=ExpTable.C_',Duration{1},'(',num2str(idx),');']);
            eval(['d=ExpTable.A_',Duration{2},'(',num2str(idx),');']);
            eval(['e=ExpTable.B_',Duration{2},'(',num2str(idx),');']);
            eval(['f=ExpTable.C_',Duration{2},'(',num2str(idx),');']);
            GEEPI = [a,b,c,d,e,f];
            for gl=1:numel(GEEPI)
                scans2((idx-1)*numel(GEEPI)+gl,1) = {fullfile(path,'Functions','tsfMRI',num2str(GEEPI(gl)),'con_0001.nii')};
            end
        end
        session_num = numel(GEEPI);
        a=session_num;
        n=subj_num;
        for m=1:n
            matrix2((a*m-a+1):a*m,1)=ones(a,1);
            matrix2((a*m-a+1):a*m,2)=ones(a,1)*2;
            matrix2((a*m-a+1):a*m,3)=(m+l)*ones(a,1);
            matrix2((a*m-a+1):a*m,4)= [1:a]';
        end
        
        
        matlabbatch{1}.spm.stats.factorial_design.dir ={dest};
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = [scans1;scans2];
        matrix=[matrix1; matrix2];

        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix = matrix;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'group';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'subj';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'session';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'cre-ctrl';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights(1:l) = 1/l;
          matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights(l+1:l+n) = -1/n;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
        F = spm_figure('GetWin');
        spm_jobman('run',matlabbatch);
        hgexport(figure(F), fullfile(dest,strcat('block design')), hgexport('factorystyle'), 'Format', 'tiff');
        clear  matlabbatch;

        

        % colormap
        template = [codepath,'\meanT2.nii'];
        spmT_file = [dest '\spmT_0001.nii'];
        lmask = [codepath,'\lmask_Mouse_v38.nii'];
        Colormap(   'statfile',spmT_file,...
            'bgmfile',template,...
            'slice',9:2:35,...
            'bar_value',[-15 -0.7 0.7 15],...
            'dest',dest,...
            'mapname','tmap',...
            'denoi_profile',lmask,...
            'cluster',10,...
            'adjust_method','FDR',...
            'corrected_p',0.05);%,...
    
    
end