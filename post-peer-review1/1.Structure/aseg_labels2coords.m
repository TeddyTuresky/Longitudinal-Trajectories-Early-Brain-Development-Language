function coord2 = aseg_labels2coords(aseg, labels)

% Generates coordinates from list of aseg labels
% assumes 1 mm^3 voxel dimensions

% index voxels 
for i = 1:size(labels,2)
    [x, y, z] = ind2sub(size(aseg), find(aseg == labels(i)));
    labs = (zeros(size(x,1),1)) + labels(i);
    if i == 1
        coord2 = [x y z labs];
    else
        coord2 = vertcat(coord2, [x y z labs]);
    end
    clearvars x y z labs
end




