TrPath = 'H:/WS/MITOS_A_TrainingSet/Features/';
%TrPath = '/groups/becklab/MITOS/A_TrainingSet/Features/';
EvPath = 'H:/WS/MITOS_A_TestingSet/Features/';
%EvPath = '/groups/becklab/MITOS/A_TestingSet/Features/';

TrFeatures = readtable(strcat(TrPath,'Tr_FeaturesWithLabels_256.csv'));
EvFeatures = readtable(strcat(EvPath,'Te_FeaturesWithLabels_256.csv'));
[TrSize,NoOfFeatures] = size(TrFeatures);
[EvSize,NoOfFeatures] = size(EvFeatures);

addpath('/home/hi41/WS/MatCode/Lib/glmnet_matlab');
addpath('/home/hi41/Libraries/liblinear-2.01/matlab');

x1 = table2array(TrFeatures(:,1:NoOfFeatures-4));
y1 = table2array(TrFeatures(:,NoOfFeatures-3));
x2 = table2array(EvFeatures(:,1:NoOfFeatures-4));
y2 = table2array(EvFeatures(:,NoOfFeatures-3));

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

%Liblinear Model
model = train(y1,sp1);
[lab,acc,esti] = predict(y2,sp2,model);

TN = 0; TP = 0; FP = 0; FN = 0;
for i=1:length(lab)
    if(lab(i) == 1 && y2(i) == 1)
        TP = TP + 1;
    elseif(lab(i) == 1 && y2(i) == 0)
        FP = FP + 1;
    elseif(lab(i) == 0 && y2(i) == 1)
        FN = FN + 1;
    elseif(lab(i) == 0 && y2(i) == 0)
        TN = TN + 1;
    end
end

SVMStruct = svmtrain(x1,logical(y1),'Kernel_Function','rbf','Method','LS','RBF_Sigma',34); 
% SVMStruct = svmtrain(ftrs(tr,:),labls(tr),'Kernel_Function','linear','Method','LS');
% SVMStruct = svmtrain(ftrs(tr,:),labls(tr),'Kernel_Function','polynomial','polyorder',3,'Method','LS');
% SVMStruct = svmtrain(ftrs(tr,:),labls(tr),'Kernel_Function','quadratic','Method','LS');
[classes,dist_f] = svmclassify2(SVMStruct,x2);

cp = classperf(labls);
classperf(cp,classes,y2);
cr = cp.CorrectRate;

classes_array_F(type,:) = classes;
distances_array_F(type,:) = dist_f;  
Correct_Rate(1,type) = cr*100;


fit = glmnet(x1(1:500,1:50),Tr_L(1:500),'binomial')

rng(1);  % For reproducibility
CVSVMModel = fitcsvm(Tr,Tr_L,'Standardize',true,'ClassNames',{'1','0'},'CrossVal','on')
genError = kfoldLoss(CVSVMModel)

[B,FitInfo] = lassoglm(Tr,Tr_L,'binomial', 'NumLambda',25,'CV',5);
lassoPlot(B,FitInfo,'PlotType','CV');

indx = FitInfo.Index1SE;
B0 = B(:,indx);
nonzeros = sum(B0 ~= 0)

cnst = FitInfo.Intercept(indx);
B1 = [cnst;B0];

preds = glmval(B1,X,'logit');
histogram(Ybool - preds) % plot residuals
title('Residuals from lassoglm model')

plotResiduals(mdl)



