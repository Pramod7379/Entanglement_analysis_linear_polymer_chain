function atom_indices = segmenting(m1_ind, m2_ind,m1,m2, r_cut, abs_lnk_lim,n_c)

% Finds out which segments are close and entangled in given chains m1 & m2.
% m1 and m2 are arrays consisting the coordinates of atoms on the
% respective chains (Nx3). m1_ind and m2_ind are the respective chain
% indices.
% r_cut is the minimum cut-off distance for the atoms to be described as
% close.
% seg_min - (minimum no. of atoms required to form a proper segment - 1);
% should be an even no.
% increm - increment in no. of atoms
% abs_lnk_lim - cut-off limit on the linking number of a pair of segments to be
% classified as an entanglement.

% Initialization
atom_indices = [];


% Matrix of distance between atoms of respective chains
Mat_dist = all_atoms_dist(m1, m2);

% Check for minimum value in Mat_dist less than r_cut
if min(min(Mat_dist)) >= r_cut
    return;
end

% Chain sizes
N_m1 = length(m1);
N_m2 = length(m2);

while min(min(Mat_dist)) < r_cut
    
    % cl_ind = [C, D]; C = atom on chain 1, D = atom on chain 2, in the closest pair
    cl_ind = min_dIJ(Mat_dist);
    
    c1 = max((cl_ind(1)-n_c), 1):min((cl_ind(1)+n_c), N_m1);
    c2 = max((cl_ind(2)-n_c), 1):min((cl_ind(2)+n_c), N_m2);
    
    lnk = 0;
    length_0 = length(c1);
    length_1 = length(c2);

    for i = 1:length_0-1
        for j = 1:length_1-1
            % a, b, c and d are four vectors. Solid angle of the prism bounded by
            % these four vectors are needed in the calculation of Gaussian
            % linking number.
            
            
            
            a = m2(c2(j), :) + 0 - m1(c1(i), :);
            b = m2(c2(j), :) - m1(c1(i+1), :);
            c = m2(c2(j+1), :) - m1(c1(i+1), :);
            d = m2(c2(j+1), :) + 0 - m1(c1(i), :);
            lnk = lnk + solid_angle_quadrilateral(a', b', c', d');
        end
    end

    lnk = lnk/(4 *pi); 
    
    
    if abs(lnk) > abs_lnk_lim
        atom_indices = [atom_indices; m1_ind, m2_ind, ...
            cl_ind(1), cl_ind(2), c1(1), c1(end), c2(1), c2(end), Mat_dist(cl_ind(1), cl_ind(2)), lnk];
    end
    
    % Revised Mat_dist
    Mat_dist(c1, c2) = r_cut + 1;
end
 






end

