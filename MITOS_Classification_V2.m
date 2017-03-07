
addpath('/home/hi41/WS/MatCode/Lib/glmnet_matlab');
addpath('/home/hi41/Libraries/liblinear-2.01/matlab');

TrPath = 'H:/WS/MITOS_A_TrainingSet/Features/';
%TrPath = '/groups/becklab/MITOS/A_TrainingSet/Features/';
EvPath = 'H:/WS/MITOS_A_TestingSet/Features/';
%EvPath = '/groups/becklab/MITOS/A_TestingSet/Features/';

TrFeatures = readtable(strcat(TrPath,'Tr_FeaturesWithLabels_64.csv'));
TeFeatures = readtable(strcat(EvPath,'Te_FeaturesWithLabels_64.csv'));
[TrSize,NoOfFeatures] = size(TrFeatures);
[TeSize,NoOfFeatures] = size(TeFeatures);

x1 = table2array(TrFeatures(:,1:NoOfFeatures-4));
y1 = table2array(TrFeatures(:,NoOfFeatures-3));
x2 = table2array(TeFeatures(:,1:NoOfFeatures-4));
y2 = table2array(TeFeatures(:,NoOfFeatures-3));

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

%Simple Liblinear 2.0.1
 %model = train(training_label_vector, training_instance_matrix [,'liblinear_options', 'col']);
%[predicted_label, accuracy, decision_values/prob_estimates] = predict(testing_label_vector, testing_instance_matrix, model [, 'liblinear_options', 'col']);
model = train(y1,sp3,'-s 1 -c 7   -B 1.5 -wi 5.5');
model = train(y1,sp3,'-s 1 -c 6.5 -B 1.5 -wi 4.5');
[lab,acc,esti] = predict(y2,sp4,model);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab, 0);
[TPR, PPV, FM, Acc, FDR] = ComputePerformanceMetrics(TP, FN, TN, FP, 100)

%SVM Classifier fitcsvm

%'Solver','ISDA','L1QP','SMO'
%'Weights',
%wt = ones(size(x1,1),1)

SVMModel = fitcsvm(x3,y1,'Standardize',true,'KernelFunction','polynomial');
[lab,score] = predict(SVMModel,x4);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab);
disp('Polynomial - Normalization');
[TPR, PPV, FM, Acc, FDR] = ComputePerformanceMetrics(TP, FN, TN, FP, 100);

SVMModel = fitcsvm(x3,y1,'Standardize',true,'KernelFunction','linear');
[lab,score] = predict(SVMModel,x4);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab);
disp('Linear - Normalization');
[TPR, PPV, FM, Acc, FDR] = ComputePerformanceMetrics(TP, FN, TN, FP, 100)

SVMModel = fitcsvm(x3,y1,'KernelFunction','linear','Standardize',true,'Cost',6.5,'Bias',1.5,'autoscale',true);
[lab,score] = predict(SVMModel,x4);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab);

ScoreSVMModel = fitPosterior(SVMModel,x1,y1);
[lab,score] = predict(ScoreSVMModel,x2);
[ TP, FN, TN, FP ] = CompareLabelWithGT(y2, lab);

SVMModel = fitcsvm(x1,y1,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
SVMModel = fitcsvm(x1,y1,'KernelFunction','rbf','method','LS','rbf_sigma',34); 
SVMModel = fitcsvm(x1,y1,'KernelFunction','linear','method','LS');
SVMModel = fitcsvm(x1,y1,'KernelFunction','polynomial','polyorder',3,'method','LS');
SVMModel = fitcsvm(x1,y1,'KernelFunction','quadratic','method','LS');

%SVMModel = svmtrain(x1,y1,'Kernel_Function','rbf','Method','LS','RBF_Sigma',34);
%[classes,dist_f] = svmclassify2(SVMModel,x2);

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



