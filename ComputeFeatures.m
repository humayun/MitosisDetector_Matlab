%%%
% Author: Humayun Irshad
% Beck lab, BIDMC, Harvard Medical School.
% 2014_05_20
%%%

function [ Features ] = ComputeFeatures(RGB,BW,Centroids,PatchSize,GrayLevels)
%ComputeFeatures_MITOS compute Morphological, Intensity, Co-occurrence and 
% Run-length features in seven color channels (Red, Green, Blue, HSV(V), 
% Lab(L), H&E(H) and BR
%

    if(nargin < 5)
        PatchSize = 70;
    end
    if(nargin < 6)
        GrayLevels = 256;
    end
    
    myFcn=@(x)graycomatrix(x,'Offset',[1 0;0 1;1 1;-1 1],'Symmetric',...
        true, 'NumLevels', GrayLevels);

    Features = DataStructureDefinition();
    
    %% Morphological Features
    morphFeatures = regionprops(logical(BW), 'Area', ...
        'Perimeter', 'MajorAxisLength', 'MinorAxisLength', ...
        'Eccentricity', 'ConvexArea', 'Orientation', ...
        'EquivDiameter', 'Solidity', 'Extent', 'Perimeter');

    Features.Area = cat(1,morphFeatures.Area);
    Features.Perimeter = cat(1,morphFeatures.Perimeter);
    Features.MajorAxisLength = cat(1,morphFeatures.MajorAxisLength);
    Features.MinorAxisLength = cat(1,morphFeatures.MinorAxisLength);
    Features.Eccentricity = cat(1,morphFeatures.Eccentricity);
    Features.ConvexArea = cat(1,morphFeatures.ConvexArea);
    Features.Orientation = cat(1,morphFeatures.Orientation);
    Features.EquivDiameter = cat(1,morphFeatures.EquivDiameter);
    Features.Solidity = cat(1,morphFeatures.Solidity);
    Features.Extent = cat(1,morphFeatures.Extent);
    Features.Compactness = 4*pi*Features.Area./(Features.Perimeter).^2;
    
    Red = RGB(:,:,1);
    Green = RGB(:,:,2);
    Blue = RGB(:,:,3);
    cs = rgb2hsv(RGB);
    HSV_V = cs(:,:,3);
    HSV_V2 = uint8(GrayLevels*mat2gray(HSV_V));
    cs = rgb2lab(RGB);
    Lab_L = cs(:,:,1);
    Lab_L2 = uint8(GrayLevels*mat2gray(Lab_L));
    %[ A1, A2, A3, A4, A5, A6, A7, A8 ] = colordeconv( RGB, 2 );
    %[ DCH, HE_H, HE_E, R, M] = Deconvolve( RGB ); % Ruifrok Method (Adnan)
    [HE_H, ~, ~] = colordeconv2( RGB, 'PSL' ); % IPAL Method
    HE_H2 = uint8(GrayLevels*mat2gray(HE_H));
    BR = ComputeBlueRatio(RGB);
    BR2 = uint8(GrayLevels*mat2gray(BR));

    [L,num] = bwlabel(BW);
    
    for i = 1:num
    %% Intensity Features
    
    % 1 - Red Channel
        obj=Red((L==i));
        Features.Mean_R(i,1) = mean(obj);
        Features.Median_R(i,1) = median(obj);
        Features.MAD_R(i,1) = mad(double(obj));
        Features.SD_R(i,1) = std(double(obj));
        Features.IQR_R(i,1) = iqr(double(obj));
        Features.Skewness_R(i,1) = skewness(double(obj));
        Features.Kurtosis_R(i,1) = kurtosis(double(obj));

    % 2 - Greeen Channel
        obj=Green((L==i));
        Features.Mean_G(i,1) = mean(obj);
        Features.Median_G(i,1) = median(obj);
        Features.MAD_G(i,1) = mad(double(obj));
        Features.SD_G(i,1) = std(double(obj));
        Features.IQR_G(i,1) = iqr(double(obj));
        Features.Skewness_G(i,1) = skewness(double(obj));
        Features.Kurtosis_G(i,1) = kurtosis(double(obj));

    % 3 - Blue Channel
        obj=Blue((L==i));
        Features.Mean_B(i,1) = mean(obj);
        Features.Median_B(i,1) = median(obj);
        Features.MAD_B(i,1) = mad(double(obj));
        Features.SD_B(i,1) = std(double(obj));
        Features.IQR_B(i,1) = iqr(double(obj));
        Features.Skewness_B(i,1) = skewness(double(obj));
        Features.Kurtosis_B(i,1) = kurtosis(double(obj));

    % 4 - V (HSV) Color Channel
        obj=HSV_V((L==i));
        Features.Mean_V(i,1) = mean(obj);
        Features.Median_V(i,1) = median(obj);
        Features.MAD_V(i,1) = mad(double(obj));
        Features.SD_V(i,1) = std(double(obj));
        Features.IQR_V(i,1) = iqr(double(obj));
        Features.Skewness_V(i,1) = skewness(double(obj));
        Features.Kurtosis_V(i,1) = kurtosis(double(obj));

    % 5 - L (Lab) Color Channel
        obj=Lab_L((L==i));
        Features.Mean_L(i,1) = mean(obj);
        Features.Median_L(i,1) = median(obj);
        Features.MAD_L(i,1) = mad(double(obj));
        Features.SD_L(i,1) = std(double(obj));
        Features.IQR_L(i,1) = iqr(double(obj));
        Features.Skewness_L(i,1) = skewness(double(obj));
        Features.Kurtosis_L(i,1) = kurtosis(double(obj));

    % 6 - H (H&E) Color Deconvolution
        obj=HE_H((L==i));
        Features.Mean_H(i,1) = mean(obj);
        Features.Median_H(i,1) = median(obj);
        Features.MAD_H(i,1) = mad(double(obj));
        Features.SD_H(i,1) = std(double(obj));
        Features.IQR_H(i,1) = iqr(double(obj));
        Features.Skewness_H(i,1) = skewness(double(obj));
        Features.Kurtosis_H(i,1) = kurtosis(double(obj));

    % 7 - BlueRatio Image Channel
        obj=BR((L==i));
        Features.Mean_Br(i,1) = mean(obj);
        Features.Median_Br(i,1) = median(obj);
        Features.MAD_Br(i,1) = mad(double(obj));
        Features.SD_Br(i,1) = std(double(obj));
        Features.IQR_Br(i,1) = iqr(double(obj));
        Features.Skewness_Br(i,1) = skewness(double(obj));
        Features.Kurtosis_Br(i,1) = kurtosis(double(obj));
    
    %% Texture Features    
    
    % 1 - Red Channel
        patch = PatchExtraction(Red,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch);
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_R(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_R(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_R(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_R(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_R(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_R(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_R(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_R(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_R(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_R(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_R(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_R(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_R(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_R(i,1) = mean(glcmFeatures.inf2h);
        
    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_R(i,1) = mean(stats(:,1));
        Features.LRE_R(i,1) = mean(stats(:,2));
        Features.GLN_R(i,1) = mean(stats(:,3));
        Features.RLN_R(i,1) = mean(stats(:,4));
        Features.RP_R(i,1) = mean(stats(:,5));
        Features.LGRE_R(i,1) = mean(stats(:,6));
        Features.HGRE_R(i,1)= mean(stats(:,7));
        Features.SRLGE_R(i,1) = mean(stats(:,8));
        Features.SRHGE_R(i,1) = mean(stats(:,9));
        Features.LRLGE_R(i,1) = mean(stats(:,10));
        Features.LRHGE_R(i,1) = mean(stats(:,11));

    % 2 - Green Channel
        patch = PatchExtraction(Green,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch); 
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_G(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_G(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_G(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_G(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_G(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_G(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_G(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_G(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_G(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_G(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_G(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_G(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_G(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_G(i,1) = mean(glcmFeatures.inf2h);

    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_G(i,1) = mean(stats(:,1));
        Features.LRE_G(i,1) = mean(stats(:,2));
        Features.GLN_G(i,1) = mean(stats(:,3));
        Features.RLN_G(i,1) = mean(stats(:,4));
        Features.RP_G(i,1) = mean(stats(:,5));
        Features.LGRE_G(i,1) = mean(stats(:,6));
        Features.HGRE_G(i,1)= mean(stats(:,7));
        Features.SRLGE_G(i,1) = mean(stats(:,8));
        Features.SRHGE_G(i,1) = mean(stats(:,9));
        Features.LRLGE_G(i,1) = mean(stats(:,10));
        Features.LRHGE_G(i,1) = mean(stats(:,11));

    % 3 - Blue Channel
        patch = PatchExtraction(Blue,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch); 
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_B(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_B(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_B(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_B(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_B(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_B(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_B(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_B(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_B(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_B(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_B(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_B(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_B(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_B(i,1) = mean(glcmFeatures.inf2h);

    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_B(i,1) = mean(stats(:,1));
        Features.LRE_B(i,1) = mean(stats(:,2));
        Features.GLN_B(i,1) = mean(stats(:,3));
        Features.RLN_B(i,1) = mean(stats(:,4));
        Features.RP_B(i,1) = mean(stats(:,5));
        Features.LGRE_B(i,1) = mean(stats(:,6));
        Features.HGRE_B(i,1)= mean(stats(:,7));
        Features.SRLGE_B(i,1) = mean(stats(:,8));
        Features.SRHGE_B(i,1) = mean(stats(:,9));
        Features.LRLGE_B(i,1) = mean(stats(:,10));
        Features.LRHGE_B(i,1) = mean(stats(:,11));
    
    % 4 - V (HSV) Channel
        patch = PatchExtraction(HSV_V2,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch); 
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_V(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_V(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_V(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_V(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_V(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_V(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_V(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_V(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_V(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_V(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_V(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_V(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_V(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_V(i,1) = mean(glcmFeatures.inf2h);

    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_V(i,1) = mean(stats(:,1));
        Features.LRE_V(i,1) = mean(stats(:,2));
        Features.GLN_V(i,1) = mean(stats(:,3));
        Features.RLN_V(i,1) = mean(stats(:,4));
        Features.RP_V(i,1) = mean(stats(:,5));
        Features.LGRE_V(i,1) = mean(stats(:,6));
        Features.HGRE_V(i,1)= mean(stats(:,7));
        Features.SRLGE_V(i,1) = mean(stats(:,8));
        Features.SRHGE_V(i,1) = mean(stats(:,9));
        Features.LRLGE_V(i,1) = mean(stats(:,10));
        Features.LRHGE_V(i,1) = mean(stats(:,11));

    % 5 - L (Lab) Channel
        patch = PatchExtraction(Lab_L2,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch); 
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_L(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_L(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_L(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_L(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_L(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_L(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_L(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_L(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_L(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_L(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_L(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_L(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_L(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_L(i,1) = mean(glcmFeatures.inf2h);

    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_L(i,1) = mean(stats(:,1));
        Features.LRE_L(i,1) = mean(stats(:,2));
        Features.GLN_L(i,1) = mean(stats(:,3));
        Features.RLN_L(i,1) = mean(stats(:,4));
        Features.RP_L(i,1) = mean(stats(:,5));
        Features.LGRE_L(i,1) = mean(stats(:,6));
        Features.HGRE_L(i,1)= mean(stats(:,7));
        Features.SRLGE_L(i,1) = mean(stats(:,8));
        Features.SRHGE_L(i,1) = mean(stats(:,9));
        Features.LRLGE_L(i,1) = mean(stats(:,10));
        Features.LRHGE_L(i,1) = mean(stats(:,11));

    % 6 - H (H&E) Channel
        patch = PatchExtraction(HE_H2,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch); 
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_H(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_H(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_H(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_H(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_H(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_H(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_H(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_H(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_H(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_H(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_H(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_H(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_H(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_H(i,1) = mean(glcmFeatures.inf2h);
        
    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_H(i,1) = mean(stats(:,1));
        Features.LRE_H(i,1) = mean(stats(:,2));
        Features.GLN_H(i,1) = mean(stats(:,3));
        Features.RLN_H(i,1) = mean(stats(:,4));
        Features.RP_H(i,1) = mean(stats(:,5));
        Features.LGRE_H(i,1) = mean(stats(:,6));
        Features.HGRE_H(i,1)= mean(stats(:,7));
        Features.SRLGE_H(i,1) = mean(stats(:,8));
        Features.SRHGE_H(i,1) = mean(stats(:,9));
        Features.LRLGE_H(i,1) = mean(stats(:,10));
        Features.LRHGE_H(i,1) = mean(stats(:,11));
        

    % 7 - BlueRatio Image
        patch = PatchExtraction(BR2,PatchSize,Centroids(i,:));
    % CM
        glcm = myFcn(patch); 
        glcmFeatures = GLCM_Features(glcm,0);
        Features.Autocorrelation_Br(i,1) = mean(glcmFeatures.autoc);
        Features.CorrelationP_Br(i,1) = mean(glcmFeatures.corrp);
        Features.Contrast_Br(i,1) = mean(glcmFeatures.contr);
        Features.ClusterShade_Br(i,1) = mean(glcmFeatures.cshad);
        Features.ClusterProminence_Br(i,1) = mean(glcmFeatures.cprom);
        Features.Energy_Br(i,1) = mean(glcmFeatures.energ);
        Features.Entropy_Br(i,1) = mean(glcmFeatures.entro);
        Features.HomogeneityP_Br(i,1) = mean(glcmFeatures.homop);
        Features.InverseDiffNormalized_Br(i,1) = mean(glcmFeatures.indnc);
        Features.InverseDiffMomentNormalized_Br(i,1) = mean(glcmFeatures.idmnc);
        Features.Dissimilarity_Br(i,1) = mean(glcmFeatures.dissi);
        Features.MaxProbability_Br(i,1) = mean(glcmFeatures.maxpr);
        Features.InfoMeasureCorr1_Br(i,1) = mean(glcmFeatures.inf1h);
        Features.InformeasureCorr2_Br(i,1) = mean(glcmFeatures.inf2h);
        
    % RL
        [RLM,~]= grayrlmatrix(patch,'OFFSET',[1;2;3;4],'NumLevels', ...
            GrayLevels,'GrayLimits',[]); 
        stats = grayrlprops(RLM,{'SRE','LRE','GLN','RLN','RP','LGRE',...
            'HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
        Features.SRE_Br(i,1) = mean(stats(:,1));
        Features.LRE_Br(i,1) = mean(stats(:,2));
        Features.GLN_Br(i,1) = mean(stats(:,3));
        Features.RLN_Br(i,1) = mean(stats(:,4));
        Features.RP_Br(i,1) = mean(stats(:,5));
        Features.LGRE_Br(i,1) = mean(stats(:,6));
        Features.HGRE_Br(i,1)= mean(stats(:,7));
        Features.SRLGE_Br(i,1) = mean(stats(:,8));
        Features.SRHGE_Br(i,1) = mean(stats(:,9));
        Features.LRLGE_Br(i,1) = mean(stats(:,10));
        Features.LRHGE_Br(i,1) = mean(stats(:,11));
    end
    Features = struct2table(Features);
end

