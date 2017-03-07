Path = 'H:/WS/Analysis/';
TrFeatures = readtable(strcat(Path,'Tr_Features_8651.csv'));
TeFeatures = readtable(strcat(Path,'Te_Features_4351.csv'));
[TrSize,NoOfFeatures] = size(TrFeatures);
[TeSize,NoOfFeatures] = size(TeFeatures);
x1 = table2array(TrFeatures(:,1:NoOfFeatures-1));
y1 = table2array(TrFeatures(:,NoOfFeatures));
x2 = table2array(TeFeatures(:,1:NoOfFeatures-1));
y2 = table2array(TeFeatures(:,NoOfFeatures));

x1 = x1(:,all(~isnan(x1)));
x2 = x2(:,all(~isnan(x1)));
x1 = x1(:,all(~isnan(x2)));
x2 = x2(:,all(~isnan(x2)));

x3 = Normalize(x1);
x4 = Normalize(x2);

sp1 = sparse(x1);
sp2 = sparse(x2);
sp3 = sparse(x3);
sp4 = sparse(x4);

model = train(y1,sp3,'-s 1 -c 7   -B 1.5 -wi 5.5');
model = train(y1,sp3,'-s 1 -c 6.5 -B 1.5 -wi 4.5');
[lab,acc,esti] = predict(y2,sp4,model);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab, 0);
[TPR, PPV, FM, Acc, FDR] = ComputePerformanceMetrics(TP, FN, TN, FP, 100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

TrP = 'H:/Datasets/Contests/MITOS_2012/Aperio/Training/Features_Patch_256/';
TrFeatures = readtable(strcat(TrP,'Training_Patch_256.csv'));
[TrSize,NoOfFeatures] = size(TrFeatures);
x1 = table2array(TrFeatures(:,3:NoOfFeatures-1));
labels = table2array(TrFeatures(:,NoOfFeatures));
y1 = zeros(TrSize,1);
for i=1:TrSize
    if(strcmp(labels(i),'Mitosis'))
        y1(i,1) = 0;
    else
        y1(i,1) = 1;
    end
end

TeP = 'H:/Datasets/Contests/MITOS_2012/Aperio/Testing/Features_Patch_256/';
TeFeatures = readtable(strcat(TeP,'Testing_Patch_256.csv'));
[TeSize,NoOfFeatures] = size(TeFeatures);
x2 = table2array(TeFeatures(:,3:NoOfFeatures-1));
labels = table2array(TeFeatures(:,NoOfFeatures));
y2 = zeros(TeSize,1);
for i=1:TeSize
    if(strcmp(labels(i),'Mitosis'))
        y2(i,1) = 0;
    else
        y2(i,1) = 1;
    end
end

x1 = x1(:,all(~isnan(x1)));
x2 = x2(:,all(~isnan(x1)));
x1 = x1(:,all(~isnan(x2)));
x2 = x2(:,all(~isnan(x2)));

x3 = Normalize(x1);
x4 = Normalize(x2);

sp1 = sparse(x1);
sp2 = sparse(x2);
sp3 = sparse(x3);
sp4 = sparse(x4);

model = train(y1,sp3,'-s 1 -c 7   -B 1.5 -wi 5.5');
model = train(y1,sp3,'-s 1 -c 6.5 -B 1.5 -wi 4.5');
[lab,acc,esti] = predict(y2,sp4,model);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab, 0);
[TPR, PPV, FM, Acc, FDR] = ComputePerformanceMetrics(TP, FN, TN, FP, 100);
