%%%
% Author: Humayun Irshad
% Beck lab, BIDMC, Harvard Medical School.
% 2014_05_20
%%%
function [ vTP, vFP, vFN, Centroids ] = ComputeLabels( GTCentroids, Centroids, Radius, ResolutionX, ResolutionY )
%ComputeLabels Summary of this function goes here
%   Detailed explanation goes here

    GTCentroids = double(GTCentroids);
    GTCentroids(:,3) = 0;
    Centroids = double(Centroids);
    TP=0;FP=0;FN=0;
    vTP=[];
    vFN=[];
    vFP=[];
    [sizeGT,~] = size(GTCentroids);
    [sizeCen,~] = size(Centroids);
    
    for j = 1:sizeCen
        found = false;
        for k = 1:sizeGT
            dif_x = ((GTCentroids(k,1)-Centroids(j,1))*ResolutionX)^2;
            dif_y = ((GTCentroids(k,2)-Centroids(j,2))*ResolutionY)^2;
            if(sqrt(double(dif_x + dif_y)) < Radius)
                GTCentroids(k,3) = 1;
                found = true;
                break;
            end
        end
        if(found)
            TP = TP + 1;
            vTP(TP,1) = Centroids(j,1);
            vTP(TP,2) = Centroids(j,2);
            Centroids(j,3) = 1;
        else
            FP=FP+1;
            vFP(FP,1) = Centroids(j,1);
            vFP(FP,2) = Centroids(j,2);
            Centroids(j,3) = 0;
        end
    end
    for j = 1:sizeGT
        if(~GTCentroids(j,3))
            FN = FN + 1;
            vFN(FN,:) = GTCentroids(j,:);
        end
    end
end

