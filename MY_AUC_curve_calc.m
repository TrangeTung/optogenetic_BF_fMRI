function Auc = MY_AUC_curve_calc(Vy,Y);

RAM = Vy;%sort((Vy),'ascend');
ct=1000;
K = round((1:ct)/ct*numel(RAM));
Q = RAM(K);%floor(numel(RAM)/1000):end);
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

end