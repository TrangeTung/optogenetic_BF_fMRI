% clc;clear
close all
WholePath = 'F:\BF_optogentics\';
codepath = 'F:\BF_optogentics\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

%
filepath = 'F:\VIDEO_BF\Package\Data\';
Type = {'Chat';'PV';'SOM';'VGLUT2';'Ctrl';};
Excel = fullfile(WholePath,'BehaviorResults.xlsx');
delete(Excel);

for tl = 1:numel(Type)
    fs=6;
    Vdir = dir(fullfile(filepath,[Type{tl},'_*']));
    
    Speed = zeros(numel(Vdir),6);
    Grooming = zeros(numel(Vdir),6);
    NovelExp = zeros(numel(Vdir),6);
    FamiliarExp = zeros(numel(Vdir),6);
    TotalExp = zeros(numel(Vdir),6);
    NovelRatio = zeros(numel(Vdir),6);
    FamiliarRatio = zeros(numel(Vdir),6);
    QWTime = zeros(numel(Vdir),6);
    QWplusGrooming = zeros(numel(Vdir),6);
    
    for vl=1:numel(Vdir)
        vl
        
        Vname = Vdir(vl).name;
        
        CSVname = fullfile(filepath,Vname,'videos-processed',...
            ['ch',num2str(5,'%02d'),'_',Vname,'_processed.csv']);
        [~,~,CellData] = xlsread(CSVname);
        CSVTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
        
        ts = 30*60*fs+1+(1:60*49*fs);%30*60*fs+1:27000;
        S = CSVTable.Speed(ts);
        S(S>100)=nan;S=fillmissing(S,'linear');
        Speed(vl,:) = mean(reshape(S,[],6),1);
        Grooming(vl,:) = sum(reshape(CSVTable.Grooming(ts),[],6),1)/fs;
        x = abs(gradient(CSVTable.NovelObject_Exploration(ts)));
        NovelExp(vl,:) = sum(reshape(x,[],6),1);
        x = abs(gradient(CSVTable.FamiliarObject_Exploration(ts)));
        FamiliarExp(vl,:) = sum(reshape(x,[],6),1);
        TotalExp(vl,:) = NovelExp(vl,:) + FamiliarExp(vl,:);
        NovelRatio(vl,:) = NovelExp(vl,:) ./ (NovelExp(vl,:) + FamiliarExp(vl,:));
        NovelRatio(isnan(NovelRatio))=0;
        FamiliarRatio(vl,:) = FamiliarExp(vl,:) ./ (NovelExp(vl,:) + FamiliarExp(vl,:));
        FamiliarRatio(isnan(FamiliarRatio))=0;        
        QWTime(vl,:) = sum(reshape(CSVTable.QuietAwake(ts),[],6),1)/fs;
        QWplusGrooming(vl,:) = QWTime(vl,:)+Grooming(vl,:);
        
    end
    
    Excel = fullfile(WholePath,'BehaviorResults.xlsx');
    
    xlswrite(Excel,{Type{tl}},'Speed',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,Speed(:),'Speed',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'Grooming',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,Grooming(:),'Grooming',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'NovelExp',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,NovelExp(:),'NovelExp',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'FamiliarExp',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,FamiliarExp(:),'FamiliarExp',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'TotalExp',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,TotalExp(:),'TotalExp',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'NovelRatio',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,NovelRatio(:),'NovelRatio',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'FamiliarRatio',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,FamiliarRatio(:),'FamiliarRatio',[char(double('A')+tl-1),'2']);

    xlswrite(Excel,{Type{tl}},'QuietAwakeTime',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,QWTime(:),'QuietAwakeTime',[char(double('A')+tl-1),'2']);
    
    xlswrite(Excel,{Type{tl}},'QWplusGrooming',[char(double('A')+tl-1),'1']);
    xlswrite(Excel,QWplusGrooming(:),'QWplusGrooming',[char(double('A')+tl-1),'2']);
end
%}






