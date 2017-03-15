%%%
% Author: Humayun Irshad
% Beck lab, BIDMC, Harvard Medical School.
% 2014_05_20
%%%
function MITOS_SegFeature()

addpath('/home/hi41/WS/MatCode/Lib');
addpath('/home/hi41/WS/MatCode/Features');

% Specify the Path for Training Data set
RGBPath = 'H:\Datasets\Contests\MITOS_2012\Aperio\';

% Specify the Nuclei Segmentation Method
SegMethod = 2;

% Creating folder to save segmentation and overlay images of candidates
SegPath = strcat(RGBPath,'Segmentation_M',num2str(SegMethod),'/');
mkdir(SegPath);

% Creating folder to save features and centroids of candidates
FeaturesPath = strcat(RGBPath,'MITOS_A_Training_M',num2str(SegMethod),'/');
mkdir(FeaturesPath);

% Image Extension
ImageExt = '.bmp';
% Read all the images in training set folder with specific image extension
srcFiles = dir(strcat(RGBPath,'*',ImageExt));

% Split the touching and overlapping nuclei
SplitNuclei = 1;
% Pixel Sizes
ResolutionX = 0.2456;
ResolutionY = 0.2273;
% Radius of Mitotic Nuclei in pixels
Radius = 8;
% Number of gray levels for computing texture features
GrayLevels = 256;
% Size of patch for computing texture features
PatchSize = 70;
% candidate size ranges (min, max)
MinPixel = round(12.5 / (ResolutionX * ResolutionY));
MaxPixel = round(120 / (ResolutionX * ResolutionY));
% Feature strcture definition
Features256 = struct2table(DataStructureDefinition());
Features256.Labels = 0;
Features256(1,:) = [];
% result table definition
DetectionResults = table;

for i = 1:length(srcFiles)
    [~, ImageName, ImageExt] = fileparts(srcFiles(i).name);
    disp(strcat('Iteration:',num2str(i),'    Image Name:',ImageName));
    RGB = imread( strcat( RGBPath, ImageName, ImageExt ));

    % Detect and Segment the candidate in image 
    [~, BW] = NucleiDetection_HE(RGB,MinPixel,MaxPixel,SegMethod,SplitNuclei);
    imwrite(BW, strcat(SegPath,ImageName,'_Binary.png'));
    %BW = imread( strcat(SegPath,ImageName,'_Binary.png'));
    
    % Overlay the candidate segmentation on H&E Images
    overlay = superimpose( RGB, BW, [0 1 0]); 
    imwrite(overlay, strcat(SegPath,ImageName, '_Overlay.png'));

    % Extract the centroids for each candidate and find it is mitosis or
    % not and save in CSV file
    [L,~] = bwlabel(BW);
    centroids = regionprops(L,'Centroid'); 
    centroids = floor(cat(1,centroids.Centroid));
    GTCentroids = csvread(strcat(RGBPath,ImageName,'_Centroid.csv'));
    [vTP, vFP, vFN, CentroidsLabel ] = ComputeLabels ...
            (GTCentroids,centroids,Radius,ResolutionX,ResolutionY);
    csvwrite(strcat(FeaturesPath,ImageName,'_CentroidsLabels.csv'),CentroidsLabel);
    st = struct;
    st.ImageName =  str2double(ImageName);
    [st.GT,~] = size(GTCentroids);
    [st.TP,~] = size(vTP);
    [st.FN,~] = size(vFN);
    [st.FP,~] = size(vFP);
    DetectionResults(i,:) = struct2table(st);
    
    % Compute the features for each detected candidate
    Features = ComputeFeatures_MITOS(RGB,BW,centroids,PatchSize,GrayLevels);
    % Copy the label of the candidate at the end of feature 
    Features.Labels = CentroidsLabel(:,3);
    % Write the features and label of each candidate of this image
    writetable(Features,strcat(FeaturesPath,ImageName,'_Features_', ...
                                            num2str(GrayLevels),'.csv'));
    % Concatenate the features of all images into one feature table
    Features256 = vertcat(Features256,Features);

end
% Write the features with labels of all candidates in all the images into a
% csv file
writetable(Features256,strcat(FeaturesPath,'MITOS_A_Training_M', ...
                    num2str(SegMethod),'_', num2str(GrayLevels),'.csv'));

%% Evalution Set

% Specify the Path for Testing Data set
RGBPath = 'H:\Datasets\Contests\MITOS_2012\Aperio\';

% Creating folder to save segmentation and overlay images of candidates
SegPath = strcat(RGBPath,'Segmentation_M',num2str(SegMethod),'/');
mkdir(SegPath);
% Creating folder to save features and centroids of candidates
FeaturesPath = strcat(RGBPath,'MITOS_A_Testing_M',num2str(SegMethod),'/');
mkdir(FeaturesPath);

% Read all the images in training set folder with specific image extension
srcFiles = dir(strcat(RGBPath,'*',ImageExt));

% Feature strcture definition
Features256 = struct2table(DataStructureDefinition());
Features256.Labels = 0;
Features256(1,:) = [];
% result table definition
DetectionResults = table;

for i = 1:length(srcFiles)
    [~, ImageName, ImageExt] = fileparts(srcFiles(i).name);
    disp(strcat('Iteration:',num2str(i),'    Image Name:',ImageName));
    RGB = imread( strcat( RGBPath, ImageName, ImageExt ));

    % Detect and Segment the candidate in image 
    [~, BW] = NucleiDetection_HE(RGB,MinPixel,MaxPixel,SegMethod,SplitNuclei);
    imwrite(BW, strcat(SegPath,ImageName,'_Binary.png'));
    %BW = imread( strcat(SegPath,ImageName,'_Binary.png'));
    
    % Overlay the candidate segmentation on H&E Images
    overlay = superimpose( RGB, BW, [0 1 0]); 
    imwrite(overlay, strcat(SegPath,ImageName, '_Overlay.png'));

    % Extract the centroids for each candidate and find it is mitosis or
    % not and save in CSV file
    [L,~] = bwlabel(BW);
    centroids = regionprops(L,'Centroid'); 
    centroids = floor(cat(1,centroids.Centroid));
    GTCentroids = csvread(strcat(RGBPath,ImageName,'_Centroid.csv'));
    [vTP, vFP, vFN, CentroidsLabel ] = ComputeLabels ...
                    (GTCentroids,centroids,Radius,ResolutionX,ResolutionY);
    csvwrite(strcat(FeaturesPath,ImageName,'_CentroidsLabels.csv'),CentroidsLabel);
    st = struct;
    st.ImageName =  str2double(ImageName);
    [st.GT,~] = size(GTCentroids);
    [st.TP,~] = size(vTP);
    [st.FN,~] = size(vFN);
    [st.FP,~] = size(vFP);
    DetectionResults(i,:) = struct2table(st);
    
    % Compute the features for each detected candidate
    Features = ComputeFeatures_MITOS(RGB,BW,centroids,PatchSize,GrayLevels);
    % Copy the label of the candidate at the end of feature 
    Features.Labels = CentroidsLabel(:,3);
    % Write the features and label of each candidate of this image
    writetable(Features,strcat(FeaturesPath,ImageName,'_Features_', ...
                                            num2str(GrayLevels),'.csv'));
    % Concatenate the features of all images into one feature table
    Features256 = vertcat(Features256,Features);

end
% Write the features with labels of all candidates in all the images into a
% csv file
writetable(Features256,strcat(FeaturesPath,'MITOS_A_Testing_M', ...
                    num2str(SegMethod),'_', num2str(GrayLevels),'.csv'));
