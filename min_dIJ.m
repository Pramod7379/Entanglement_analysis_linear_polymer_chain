function A = min_dIJ(Mat_dist)
% Finds the minimum distance pair of atoms on two different chains
[C, D] = find(Mat_dist == min(min(Mat_dist)));
A = [C, D];
end

