function [ Patch ] = PatchExtraction( Image, PatchSize, Centroid )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    [Width,Height] = size(Image);
    
    if(Centroid(1) <= PatchSize/2)
        top_x = 1;
    else
        top_x = Centroid(1) - PatchSize/2;
    end
    if( (top_x + PatchSize) > Width)
        top_x = top_x - (top_x + PatchSize - Width);
    end
    
    if(Centroid(2) <= PatchSize/2)
        top_y = 1;
    else
        top_y = Centroid(2) - PatchSize/2;
    end
    if( (top_y + PatchSize) > Height)
        top_y = top_y - (top_y + PatchSize - Height);
    end
    if(ndims(Image) == 3)
        Patch = Image(top_x:top_x+PatchSize,top_y:top_y+PatchSize,:);
    else
        Patch = Image(top_x:top_x+PatchSize,top_y:top_y+PatchSize);
    end
end

