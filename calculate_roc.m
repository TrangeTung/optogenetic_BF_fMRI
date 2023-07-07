function [PX,PY,Auc] = calculate_roc(predict, ground_truth)
 
pos_num = sum(ground_truth == 1);
neg_num = sum(ground_truth == 0);
 
m = length(ground_truth);
[~, index] = sort(predict);
ground_truth = ground_truth(index);
PX = zeros(m+1,1);
PY = zeros(m+1,1);
Auc = 0;
PX(1) = 1; PY(1) = 1;
 
for i = 2:m
    TP = sum(ground_truth(i:m)==1);
    FP = sum(ground_truth(i:m)==0);
    PX(i) = FP/neg_num;
    PY(i) = TP/pos_num;
    Auc = Auc + (PY(i)+PY(i-1))*(PX(i-1)-PX(i))/2;   
end
PX(m+1) = 0;
PY(m+1) = 0;
Auc = Auc + PY(m)*PX(m)/2;