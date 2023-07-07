clc;clear
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

Excel = fullfile(codepath,'Exp_recording.xlsx');
[~,~,CellData] = xlsread(Excel,'RawData');
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
Reg_choices = {'rp';'rp"';'10PCs'};

CellType = {'PV';'SOM';'CHAT';'Vglut2';'Ctrl'};
Duration = {'05';'2'};

% all FD
%{
F =  figure;
for cl = 1:numel(CellType)
        
        clear matrix matrix1 matrix2 scans
        
        Excel = fullfile(codepath,'Exp_recording.xlsx');
        sheet = [CellType{cl},'_2nd'];
        [~,~,CellData] = xlsread(Excel,sheet);
        ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
        
        FDall = [];
        
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
                filepath = fullfile(path,'Results',num2str(GEEPI(gl)));
                rp = load(fullfile(filepath,'rp_mUw2dseq.txt'));
                drp = rp - [rp(1,:);rp(1:end-1,:)];
                drp(:,1:3)=drp(:,1:3)/20;
                drp(:,4:6)=drp(:,4:6)*5;
                FD=sqrt( drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2  );
                
                FDall = cat(1,FDall,FD(:)*1000); % μm
            end
        end
        
        histogram(FDall,200,'Normalization','count',...
            'DisplayStyle','stairs');
        y = [numel(find(FDall>75/2))/numel(FDall),numel(find(FDall>75/4))/numel(FDall)]*100
        hold on;
end
legend(CellType);
xlim([1 100])
set(gca,'xscale','log','tickdir','out');

%}


% PCs
%{
for idx = 3:3:51
    idx
    
    path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
    RARE =  ExpTable.T2(idx);
    EPI = [ExpTable.A(idx),ExpTable.B(idx),ExpTable.C(idx),...
        ExpTable.D(idx),ExpTable.E(idx),ExpTable.F(idx)];
    GEEPI = EPI(~isnan(EPI));
    EPItp = [ExpTable.TOPUP1(idx),ExpTable.TOPUP2(idx)];
    
    for gl=1:numel(GEEPI)
        
%{
        filepath = fullfile(path,'Results',num2str(GEEPI(gl)));
        ihdr = spm_vol(fullfile(filepath,'Uw2dseq.nii'));
        fMRI_4D = spm_read_vols_4D(ihdr);
        
        lmask = spm_read_vols(spm_vol(fullfile(path,'Results','EPI_mask.nii')));
        V = fmask(fMRI_4D,lmask);
        dd = V - [V(:,1),V(:,1:end-1)];
        DVARS = rms(dd); DVARS(1)=median(DVARS);
        save(fullfile(filepath,'DVARS.txt'),'DVARS','-ascii');
        
        MusMask = spm_read_vols(spm_vol(fullfile(path,'Results','Mus_mask.nii')));
        fMRI_RAM = fmask(fMRI_4D,MusMask);
        fMRI_noise_ready_PCA = (fMRI_RAM-mean(fMRI_RAM,2))./std(fMRI_RAM,0,2);
        fMRI_noise_ready_PCA(isnan(fMRI_noise_ready_PCA))=0;
        [~, score, latent, ~, explained, ~] = pca(fMRI_noise_ready_PCA','NumComponents',150,'algorithm','svd');
        save(fullfile(filepath,'PCs.txt'),'score','-ascii');
        ram = [latent explained];
        save(fullfile(filepath,'score and explained.txt'),'ram','-ascii');
        
        data_ready_regress = V;
        X = load(fullfile(filepath,'Multi_Regessor.txt'));
        RegFunc = [ones(length(X),1),X];
        for iii=1:size(data_ready_regress,1)
           [Beta,~,Residual] = regress(double(data_ready_regress(iii,:)'),RegFunc);
           data_ready_regress(iii,:) = Residual + Beta(1);
        end
        V0 = data_ready_regress;
        dd = V0 - [V0(:,1),V0(:,1:end-1)];
        DVARS = rms(dd); DVARS(1)=median(DVARS);
        save(fullfile(filepath,'DVARS_denoised.txt'),'DVARS','-ascii');
        
        
        
        filepath = fullfile(path,'Results',num2str(GEEPI(gl)));
        ihdr = spm_vol(fullfile(filepath,'Uw2dseq.nii'));
        fMRI_4D = spm_read_vols_4D(ihdr);
        lmask = spm_read_vols(spm_vol(fullfile(path,'Results','EPI_mask.nii')));
        V = fmask(fMRI_4D,lmask);
        data_ready_regress = V;
        X = load(fullfile(filepath,'Multi_Regessor.txt'));
        RegFunc = [ones(length(X),1),X];
        for iii=1:size(data_ready_regress,1)
           [Beta,~,Residual] = regress(double(data_ready_regress(iii,:)'),RegFunc);
           data_ready_regress(iii,:) = Residual + Beta(1);
        end
        V1 = data_ready_regress;
        fMRI_4D_NEW = funmask(V1,lmask);
        for ix=1:numel(ihdr)
            ihdr(ix).fname = fullfile(filepath,'R_Uw2dseq.nii');
        end
        spm_write_vol_4D(ihdr,fMRI_4D_NEW);
        
        clear V V0 V1 fMRI* data_ready_regress
        
        
        filepath = fullfile(path,'Results',num2str(GEEPI(gl)));
        ihdr = spm_vol(fullfile(filepath,'snrmUw2dseq.nii'));
        fMRI_4D = spm_read_vols_4D(ihdr);
        lmask = spm_read_vols(spm_vol(fullfile(codepath,'lmask_Mouse_v38.nii')));
        V = fmask(fMRI_4D,lmask);
        data_ready_regress = V;
        X = load(fullfile(filepath,'Multi_Regessor.txt'));
        RegFunc = [ones(length(X),1),X];
        for iii=1:size(data_ready_regress,1)
           [Beta,~,Residual] = regress(double(data_ready_regress(iii,:)'),RegFunc);
           data_ready_regress(iii,:) = Residual + Beta(1);
        end
        fMRI_4D_NEW = funmask(data_ready_regress,lmask);
        for ix=1:numel(ihdr)
            ihdr(ix).fname = fullfile(filepath,'R_snrmUw2dseq.nii');
        end
        spm_write_vol_4D(ihdr,fMRI_4D_NEW);
        
        clear V V0 V1 fMRI* data_ready_regress

        all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(gl)},'^R_snrmUw2dseq','.nii',[1 Inf],'separate_cells');
        realign_mlb = MY_get_default_realign_batch_struct(all_func);
        F = spm_figure('GetWin');
        disp('Start to process realignment !')
        spm_jobman('run',realign_mlb);
        hgexport(figure(F), fullfile([path,'Results\',num2str(GEEPI(gl))],strcat('R_realign')), hgexport('factorystyle'), 'Format', 'tiff');
        clear realign_mlb all_func;

%}
        
        

        
    end
    
end
%}


% FD before & after regression
FDX0=[];CCX1=[];
for cl = 1:numel(CellType)
    
    a_=0;
    
    clear matrix matrix1 matrix2 scans
    
    Excel = fullfile(codepath,'Exp_recording.xlsx');
    sheet = [CellType{cl},'_2nd'];
    [~,~,CellData] = xlsread(Excel,sheet);
    ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
    
    FDall = [];
    
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
            filepath = fullfile(path,'Results',num2str(GEEPI(gl)));
            rp = load(fullfile(filepath,'rp_mUw2dseq.txt'));
            drp = rp - [rp(1,:);rp(1:end-1,:)];
            drp(:,1:3)=drp(:,1:3)/20;
            drp(:,4:6)=drp(:,4:6)*5;
            FD=sqrt( drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2  );
            
            rp = load(fullfile(filepath,'rp_R_snrmUw2dseq.txt'));
            drp = rp - [rp(1,:);rp(1:end-1,:)];
            drp(:,1:3)=drp(:,1:3)/20;
            drp(:,4:6)=drp(:,4:6)*5;
            dFD=sqrt( drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2  );
        
            a_=a_+1;
            FDX0(a_,cl*2-1) = max(FD);
            FDX0(a_,cl*2-0) = max(dFD);

        end
    end
end


% C.C. b.w. PC and DVARS
%{
CCall_PC_DVARS=[];
a=0;
for idx = 1:51
    idx
    
    path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
    RARE =  ExpTable.T2(idx);
    EPI = [ExpTable.A(idx),ExpTable.B(idx),ExpTable.C(idx),...
        ExpTable.D(idx),ExpTable.E(idx),ExpTable.F(idx)];
    GEEPI = EPI(~isnan(EPI));
    EPItp = [ExpTable.TOPUP1(idx),ExpTable.TOPUP2(idx)];
    
    for kk=GEEPI
        filepath = fullfile(path,'Results',num2str(kk));
        X = load(fullfile(filepath,'DVARS.txt'));
        PCs = load(fullfile(filepath,'PCs.txt'));
        CCx = corr(PCs,X');
        
        a=a+1;
        CCall_PC_DVARS(:,a)=CCx;
        
    end
end
Cx = abs(CCall_PC_DVARS);
Ms = mean(Cx,2);
SEMs = std(Cx,0,2);
x = 1:size(Cx,1);

F = figure;
plot(x,Ms,'r','linewidth',2); hold on
patch('XData',[x,fliplr(x)],'YData',[Ms-SEMs;flipud(Ms+SEMs)],...
    'facecolor',[1,0,0],...
    'edgecolor','none','facealpha',0.2);
xlim([1 100]);
xlabel('PCs');ylabel('abs(C.C.) b.w. PCs & DVARS');
set(gca,'xscale','log','tickdir','out');

%}

% C.C. b.w. FD and DVARS
%{
CCX0=[];CCX1=[];
F =  figure;
for cl = 1:numel(CellType)
    
    a_=0;
    
    clear matrix matrix1 matrix2 scans
    
    Excel = fullfile(codepath,'Exp_recording.xlsx');
    sheet = [CellType{cl},'_2nd'];
    [~,~,CellData] = xlsread(Excel,sheet);
    ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
    
    FDall = [];
    
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
            filepath = fullfile(path,'Results',num2str(GEEPI(gl)));
            rp = load(fullfile(filepath,'rp_mUw2dseq.txt'));
            drp = rp - [rp(1,:);rp(1:end-1,:)];
            drp(:,1:3)=drp(:,1:3)/20;
            drp(:,4:6)=drp(:,4:6)*5;
            FD=sqrt( drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2  );
            
            X0 = load(fullfile(filepath,'DVARS.txt'));
            X1 = load(fullfile(filepath,'DVARS_denoised.txt'));
            a_=a_+1;
            CCX0(a_,cl*2-1) = corr(X0(:),FD(:));
            CCX0(a_,cl*2-0) = corr(X1(:),FD(:));
        end
    end
    
end
%}

% Example of regression
%{

for idx = 1
    idx
    
    path = fullfile(WholePath,[ExpTable.path{idx},filesep]);
    RARE =  ExpTable.T2(idx);
    EPI = [ExpTable.A(idx),ExpTable.B(idx),ExpTable.C(idx),...
        ExpTable.D(idx),ExpTable.E(idx),ExpTable.F(idx)];
    GEEPI = EPI(~isnan(EPI));
    EPItp = [ExpTable.TOPUP1(idx),ExpTable.TOPUP2(idx)];
    
    for kk=GEEPI(1)
        filepath = fullfile(path,'Results',num2str(kk));
        rp = load(fullfile(filepath,'rp_mUw2dseq.txt'));
        drp = rp - [rp(1,:);rp(1:end-1,:)];
        drp(:,1:3)=drp(:,1:3)/20;
        drp(:,4:6)=drp(:,4:6)*5;
        FD=sqrt( drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2  );

        filepath = fullfile(path,'Results',num2str(kk));
        ihdr = spm_vol(fullfile(filepath,'Uw2dseq.nii'));
        fMRI_4D = spm_read_vols_4D(ihdr);

        lmask = spm_read_vols(spm_vol(fullfile(path,'Results','EPI_mask.nii')));
        V0 = fmask(fMRI_4D,lmask);
        dd = V0 - [V0(:,1),V0(:,1:end-1)];
        DVARS = rms(dd); DVARS(1)=median(DVARS);
        
        X = load(fullfile(filepath,'Multi_Regessor.txt'));
        RegFunc1 = [ones(length(X),1),X(:,1:12)];
        RegFunc2 = [ones(length(X),1),X];
        V1=V0*0;V2=V0*0;
        for iii=1:size(V0,1)
           [Beta,~,Residual] = regress(double(V0(iii,:)'),RegFunc1);
           V1(iii,:) = Residual + Beta(1);
           [Beta,~,Residual] = regress(double(V0(iii,:)'),RegFunc2);
           V2(iii,:) = Residual + Beta(1);
        end
        
        
        dd = V0 - [V0(:,1),V0(:,1:end-1)];
        DVARS0 = rms(dd); DVARS0(1)=median(DVARS0);
        dd = V1 - [V1(:,1),V1(:,1:end-1)];
        DVARS1 = rms(dd); DVARS1(1)=median(DVARS1);
        dd = V2 - [V2(:,1),V2(:,1:end-1)];
        DVARS2 = rms(dd); DVARS2(1)=median(DVARS2);
        V0_ = (V0-mean(V0,2))./std(V0,0,2); V0_(isnan(V0_))=0; GS0=mean(V0_,1);
        V1_ = (V1-mean(V1,2))./std(V1,0,2); V1_(isnan(V1_))=0; GS1=mean(V1_,1);
        V2_ = (V2-mean(V2,2))./std(V2,0,2); V2_(isnan(V2_))=0; GS2=mean(V2_,1);
        
        F = figure;
        subplot(7,1,1);
        plot(FD); xlim([1 898]);
        
        subplot(7,1,2);
        plot(DVARS0); xlim([1 898]); 
        yyaxis right;plot(GS0)
        subplot(7,1,3);
        Cx=corr(V0_',GS0');[~,s]=sort(Cx);
        imagesc(V0_(s(1:1000),:)); caxis([-2 2]); colormap('gray')
        
        subplot(7,1,4);
        plot(DVARS1); xlim([1 898]); 
        yyaxis right;plot(GS1)
        subplot(7,1,5);
        Cx=corr(V1_',GS1');[~,s]=sort(Cx);
        imagesc(V1_(s(1:1000),:)); caxis([-2 2]); colormap('gray')

        subplot(7,1,6);
        plot(DVARS2); xlim([1 898]); 
        yyaxis right;plot(GS2)
        subplot(7,1,7);
        Cx=corr(V2_',GS2');[~,s]=sort(Cx);
        imagesc(V2_(s(1:1000),:)); caxis([-2 2]); colormap('gray')
        
        
        
        F = figure;
        fx=abs(fft(GS0)); fxx=fx(1:numel(fx)/2);
        p0 = 10*log10(fxx);p0(1)=[];
        fx=abs(fft(GS1)); fxx=fx(1:numel(fx)/2);
        p1 = 10*log10(fxx);p1(1)=[];
        fx=abs(fft(GS2)); fxx=fx(1:numel(fx)/2);
        p2 = 10*log10(fxx);p2(1)=[];

        x = (1:numel(p0))/numel(p0)/0.5/2;
        plot(x,p0); hold on;plot(x,p1); plot(x,p2); 
        set(gca,'xscale','log','tickdir','out');
        xlim([-0.001 1])
        
    end
end
%}

