function Mat_dist = all_atoms_dist(m1, m2)

% m1 and m2 are arrays consisting the coordinates of atoms on the
% respective chains (Nx3)

N_m1 = size(m1);
N_m2 = size(m2);

% x-coordinate
m1_x = m1(:, 1);
m2_x = m2(:, 1);
X1 = repmat(m1_x, 1, N_m2(1));
X2_dum = repmat(m2_x, 1, N_m1(1));
X2 = (X2_dum)';
A_x = X1-X2;

% y-coordinate
m1_y = m1(:, 2);
m2_y = m2(:, 2);
Y1 = repmat(m1_y, 1, N_m2(1));
Y2_dum = repmat(m2_y, 1, N_m1(1));
Y2 = (Y2_dum)';
A_y = Y1-Y2;

% z-coordinate
m1_z = m1(:, 3);
m2_z = m2(:, 3);
Z1 = repmat(m1_z, 1, N_m2(1));
Z2_dum = repmat(m2_z, 1, N_m1(1));
Z2 = (Z2_dum)';
A_z = Z1-Z2;

% Distance calculation
Mat_dist = sqrt((A_x).^2 + (A_y).^2 + (A_z).^2);


end

