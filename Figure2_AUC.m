
clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

Mask3D = [codepath,'\lmask_Mouse_v38.nii'];
lmask=spm_read_vols(spm_vol(Mask3D));

Func3D = fullfile(WholePath,'Zerbi','ICA_zerbi15.nii');
ICA = spm_read_vols(spm_vol(Func3D));
Vica = fmask(ICA,lmask);


ICAname = {'Somatosensory';'Sensory(VC+AUD+Mot)';'DMN';'HPF';'BasalGanglia';'Olfactory'};
Type = {'vglut2';'chat';'pv';'som'};
for tloop=1:numel(Type)
Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
upper(Type{tloop}),'spmT_0001.nii');
Ya = spm_read_vols(spm_vol(Func3D));
Vy = fmask(Ya,lmask);
F=figure;
for lp=1:6
    
    Y = Vica(:,lp)>4;
    RAM = sort((Vy),'ascend');
    ct=1000;
    Q = RAM(round((1:ct)/ct*numel(RAM)));%floor(numel(RAM)/1000):end);
    %Q(1:numel(Q)-ct)=[];
    PX = zeros(ct,1);
    PY = zeros(ct,1);
    for cl=1:ct
       X = (Vy)>Q(cl); 
       TP = numel(find(X==1&Y==1));
       FN = numel(find(X==0&Y==1));
       FP = numel(find(X==1&Y==0));
       TN = numel(find(X==0&Y==0));
       Sensi = TP/(TP+FN);
       Speci = TN/(TN+FP);
       PY(cl)=Sensi;
       PX(cl)=1-Speci;
    end
    Auc = nansum(PY)/numel(PX);
    
    AucNull = zeros(1000,1);
    parfor shf=1:1000
        Vyn = Vy(randperm(numel(Vy)));
        RAM = sort((Vyn),'ascend');
        ct=1000;
        Q = RAM(round((1:ct)/ct*numel(RAM)));%floor(numel(RAM)/1000):end);
        %Q(1:numel(Q)-ct)=[];
        PX = zeros(ct,1);
        PY = zeros(ct,1);
        for cl=1:ct
           X = (Vyn)>Q(cl); 
           TP = numel(find(X==1&Y==1));
           FN = numel(find(X==0&Y==1));
           FP = numel(find(X==1&Y==0));
           TN = numel(find(X==0&Y==0));
           Sensi = TP/(TP+FN);
           Speci = TN/(TN+FP);
           PY(cl)=Sensi;
           PX(cl)=1-Speci;
        end
        AucNull(shf) = nansum(PY)/numel(PX);
    end
    pvalue(tloop,lp) = numel(find(AucNull<Auc))/numel(AucNull);
    
    plot(PX,PY); hold on;
    AUC(tloop,lp)=Auc;
end
axis equal
xlim([0 1]);ylim([0 1])
legend(ICAname,'Location','eastoutside');
set(gca,'tickdir','out','ticklength',[0.04 1],'linewidth',2,'fontsize',16);
set(gca,'box','off');
saveas(F,fullfile(WholePath,'Zerbi',[Type{tloop},'_ROC.emf']));
end



ICAname = {'Somatosensory';'Sensory(VC+AUD+Mot)';'DMN';'HPF';'BasalGanglia';'Olfactory'};
Type = {'vglut2';'chat';'pv';'som'};
for tloop=1:numel(Type)
Func3D = fullfile(WholePath,'2nd_TwoSampleTTest_results',...
upper(Type{tloop}),'spmT_0001.nii');
Ya = spm_read_vols(spm_vol(Func3D));
Vy = fmask(Ya,lmask);
for lp=1:6
    
    Y = Vica(:,lp)>3;
    X = Vy>3 | Vy<-5;
    N = numel(find(X==1&Y==1));
    pec(tloop,lp) = N;%/numel(find(X==1));
end
end
