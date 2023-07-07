%% preparation
% According to your experimental parameters and certain file direction, Change these bellow.
clc;clear
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

Excel = fullfile(codepath,'Exp_recording.xlsx');
[~,~,CellData] = xlsread(Excel,'RawData');
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
%% regressor choices
Reg_choices = {'rp';'rp"';'10PCs'};

%
for idx = 03:3:51
    idx
    
    path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
    RARE =  ExpTable.T2(idx);
    EPI = [ExpTable.A(idx),ExpTable.B(idx),ExpTable.C(idx),...
        ExpTable.D(idx),ExpTable.E(idx),ExpTable.F(idx)];
    GEEPI = EPI(~isnan(EPI));
    EPItp = [ExpTable.TOPUP1(idx),ExpTable.TOPUP2(idx)];
    
    fclose('all');
    spm('defaults', 'FMRI');
    set(spm('CreateIntWin','off'),'Visible','on');
    
    
    for flag_stage = 10
        if flag_stage == 1
            %% Bruker2nifiti
            Bruker2nifti_multislice(path,RARE,'mouse');
            Bruker2nifti_multislice(path,GEEPI,'mouse');
        end
         if flag_stage == 2
            %% GEPI //// Pre-set the FSL topup parameters
            MY_TOPUP_pars_SETUP(WholePath,ExpTable.path{idx},EPItp,'reference','mouse');
            MY_TOPUP_pars_SETUP(WholePath,ExpTable.path{idx},GEEPI,'GEEPI','mouse');
             
            folder_txt = fullfile(WholePath,'FSL_folder','folder.txt');
            fid=fopen(folder_txt,'a+');
            fprintf(fid,'%s\n',['Path_',num2str(idx),'="',ExpTable.path{idx},'";']);
            fprintf(fid,'%s\n',['Folder_',num2str(idx),'=(',num2str(GEEPI),');']);
            fclose(fid);
        end
        if flag_stage == 3
            %% Unzip .gz file and rename Uw2dseq.nii
            FSL_TopUp_Path = fullfile(WholePath,'FSL_folder',ExpTable.path{idx});
            Unzip_path = fullfile(path,'Results');
            % GEEPI
            for kk = GEEPI
                Unzip_path = fullfile(path,'Results',num2str(kk));
                delete(fullfile(Unzip_path,'Uw2dseq.nii'));
                gunzip(fullfile(FSL_TopUp_Path,['Trange_2dseq_',num2str(kk),'.nii.gz']),Unzip_path);
                cd(Unzip_path);
                eval(['!rename,',['Trange_2dseq_',num2str(kk),'.nii'],',Uw2dseq.nii']);

                ihdr=spm_vol('Uw2dseq.nii');
                img = spm_read_vols_4D(ihdr);
                hdr=ihdr(1);hdr.fname='meanUw2dseq.nii';
                IMG=mean(img,4);
                spm_write_vol(hdr,IMG);
                
            end
        end
       if flag_stage == 4
            %% generate mask
            rootpath = fullfile(path,'\Results\');
            EXE_path = fullfile(codepath,'\GrayMattter_segment\');
            MY_braintissue_segment(rootpath,GEEPI(4),'meanUw2dseq.nii',EXE_path,'EPI')
            MY_braintissue_segment(rootpath,'T2','2dseq.nii',EXE_path,'T2')
            
            ihdr = spm_vol(fullfile(rootpath,'EPI_mask','meanUw2dseq.niimeanUw2dseq.nii'));
            img = spm_read_vols(ihdr);ihdr.fname=fullfile(path,'Results','EPI_mask.nii');
            spm_write_vol(ihdr,img);
            ihdr = spm_vol(fullfile(rootpath,'T2_mask','2dseq.nii2dseq.nii'));
            img = spm_read_vols(ihdr);ihdr.fname=fullfile(path,'Results','T2_mask.nii');
            spm_write_vol(ihdr,img);
        end
        if flag_stage ==5
            %% mask all EPI and RARE
            cd(fullfile(path,'Results'));
            EPI_mask = spm_read_vols(spm_vol('EPI_mask.nii'));
            T2_mask = spm_read_vols(spm_vol('T2_mask.nii'));
            MY_mask_images(path,'Results',GEEPI,'Uw2dseq.nii',EPI_mask,'mUw2dseq.nii','EPI');
            MY_mask_images(path,'Results',RARE,'2dseq.nii',T2_mask,'T2_m.nii','T2');
            
        end
        if flag_stage == 6
            %% Realignment
            for kk=GEEPI
                all_func = MY_find_images_in_all_scans(path,'Results',{kk},'^mUw2dseq','.nii',[1 Inf],'separate_cells');
                realign_mlb = MY_get_default_realign_batch_struct(all_func);
                F = spm_figure('GetWin');
                disp('Start to process realignment !')
                spm_jobman('run',realign_mlb);
                hgexport(figure(F), fullfile([path,'Results\',num2str(kk)],strcat('realign')), hgexport('factorystyle'), 'Format', 'tiff');
                clear realign_mlb all_func;
            end
        end
        if flag_stage == 7
            for kk = GEEPI
                copyfile(fullfile(path,'Results','T2','T2_m.nii'),...
                    fullfile(path,'Results',num2str(kk)));
                ref{1,1} = fullfile(path,'Results',num2str(kk),'meanmUw2dseq.nii,1');
                source{1,1} = fullfile(path,'Results',num2str(kk),'T2_m.nii,1');
                coreg_mlb = MY_get_default_coreg_batch_struct(ref, source, {''});
                disp('Start to process coreg!');
                F = spm_figure('GetWin');
                spm_jobman('run',coreg_mlb);
                hgexport(figure(F), fullfile([path,'Results\'],strcat('coreg')), hgexport('factorystyle'), 'Format', 'tiff');
            end
            
        end
        
        
        if flag_stage == 8          
            %% T22Template coregistration
            for kk = GEEPI
                ref{1,1} = [codepath '\Template_Mouse_v38.nii,1'];
                source{1,1} = [path 'Results\',num2str(kk),'\cT2_m.nii,1'];
                all_func = MY_find_images_in_all_scans(path,'Results',{kk},'^rmUw2dseq','.nii',[1 Inf],'all_mixed');
                all_func = [all_func;[path 'Results\',num2str(kk),'\cT2_m.nii,1']];
                OldNormalize_mlb = MY_get_default_oldnormalize_batch_struct(ref, source, all_func);
                disp('Start to process OldNormalize!');
                F = spm_figure('GetWin');
                spm_jobman('run',OldNormalize_mlb);
                hgexport(figure(F), fullfile([path 'Results\' num2str(kk)], strcat('oldnormalize')), hgexport('factorystyle'), 'Format', 'tiff');
            end
           
        end
        if flag_stage == 9
            %% smooth_space
            all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^nrmUw2dseq','.nii',[1 Inf],'all_mixed');
            Smooth_mlb = MY_get_default_smooth_batch_struct(all_func);
            disp('Start to process Smooth!');
            spm_jobman('run',Smooth_mlb);
            clear Smooth_mlb;
            
        end
        if flag_stage == 10
            
            for kk=GEEPI
                Repetitions = MY_search_bruker_method('Repetitions',kk,path);
                if Repetitions==898;duration=0.5;else;duration=2;end

                block = [15,19.5,10.5,12,13.5,19.5,12,13.5,19.5,10.5,13.5,10.5,...
                    12,16.5,10.5,16.5,18,16.5,16.5,15,18,15,19.5,15,18,12,13.5,18];
                onset=15;  colorbar = [-5 -2 2 5];          
                for i=2:numel(block);onset(i)=onset(i-1)+duration+block(i-1);end

                template = [codepath,'\Template_Mouse_v38.nii'];
                result_1st = struct('weights',1,'slice',1:30,'template',template,'FDR_pvalue',0.05,'colorbar',colorbar);
                %% individual
                defined_1st = struct('Nscans','individual','filename','snrmUw2dseq','duration',duration,'onset',onset);
                MY_task_state_statistics(path,'Results',{kk},[1 Inf],Reg_choices,defined_1st,result_1st);
            end
        end
        if flag_stage ==11
            %% Group Results
            dest = [path 'Functions\tsfMRI\Allscans'];
            delete(fullfile(dest,'SPM.mat'));
            if ~exist(dest,'dir');mkdir(dest);end
            design_unit = 'secs';
            template = [codepath,'\Template_Mouse_v38.nii'];
            colorbar = [-5 -1 1 5]/2;   
            slice = 1:35;
            for kk = 1:numel(GEEPI)
                Repetitions = MY_search_bruker_method('Repetitions',GEEPI(kk),path);
                if Repetitions==898;duration=0.5;else;duration=2;end
                block = [15,19.5,10.5,12,13.5,19.5,12,13.5,19.5,10.5,13.5,10.5,...
                    12,16.5,10.5,16.5,18,16.5,16.5,15,18,15,19.5,15,18,12,13.5,18];
                onset=15;       
                for i=2:numel(block);onset(i)=onset(i-1)+duration+block(i-1);end

                Segments = MY_search_bruker_method('Segments',GEEPI(kk),path);
                EPI_TR = MY_search_bruker_method('EPI_TR',GEEPI(kk),path)/1000*Segments;
                %% regressors
                varargin = {1;{'rp_mUw2dseq.txt';'Uw2dseq.nii';['snrmUw2dseq','.nii']};...
                    {'EPI_mask.nii';'WM_mask.nii';'CSF_mask.nii';'GS_mask.nii'}};
                RegBasFuc = MY_find_regressors_in_all_scans(path,'Results',{GEEPI(kk)},[1 Inf],Reg_choices,varargin);
                all_epi = MY_find_images_in_all_scans(path,'Results',{GEEPI(kk)},'snrmUw2dseq','.nii',[1 Inf],'separate_cells');

                sess(kk) = struct('scans',all_epi,'multi_reg',{{RegBasFuc}},...
                    'cond',struct('name',['tsfMRI_',num2str(GEEPI(kk))],'onset',onset,...
                    'duration',duration,'tmod',0,'orth',1,...
                    'pmod',struct('name', {}, 'param', {}, 'poly', {})),...
                    'regress',struct('name', {}, 'val', {}),...
                    'multi',{{''}},'hpf',128);
            end
            %% 1st level analysis
            first_level_analysis_mlb = MY_1st_level_analysis_1rodentNscan_get_default_batch_struct({dest},design_unit,EPI_TR,sess);
            F = spm_figure('GetWin');
            cd(dest);
            spm_jobman('run',first_level_analysis_mlb);
            hgexport(figure(F), fullfile(dest, strcat('1st_level_analysis')), hgexport('factorystyle'), 'Format', 'tiff');
            %% estimate
            estimate_mlb = MY_1st_level_analysis_estimate_batch_struct([dest '\SPM.mat']);
            spm_jobman('run',estimate_mlb);
            %% results
            results_mlb = MY_1st_level_analysis_results_batch_struct([dest '\SPM.mat'],'tsfMRI',1,'replsc');
            spm_jobman('run',results_mlb);
            %% display
            spmT_file = [dest '\spmT_0001.nii'];
            Colormap(   'statfile',spmT_file,...
                        'bgmfile',template,...
                        'slice',slice,...
                        'bar_value',colorbar,...
                        'dest',dest,...
                        'mapname','tmap',...
                        'cluster',1);%,...
            %             'adjust_method','FDR',...
            %             'corrected_p',FDR_pvalue);

            epi_Tem_name = [path,'\','Results', '\' num2str(GEEPI(end)) '\' 'snrmUw2dseq' '.nii,1'];
            Colormap(   'statfile',spmT_file,...
                        'bgmfile',epi_Tem_name,...
                        'slice',slice,...
                        'bar_value',colorbar,...
                        'dest',dest,...
                        'mapname','tmap_epi',...
                        'cluster',1);%,...
            %             'adjust_method','FDR',...
            %             'corrected_p',FDR_pvalue);

            fclose('all');
            
        end
    end
end