function omega = solid_angle_triangle(a, b, c)
%This function takes the three vectors a, b and c, and calculates the solid
% using formula given by Van Oosterom and Strackee.

det = sum([b(2)*c(3) - b(3)*c(2); b(3)*c(1) - b(1)*c(3); b(1)*c(2) - b(2)*c(1)] .* a);

al = norm(a);
bl = norm(b);
cl = norm(c);
 
div = al * bl * cl + a' * b * cl + c' * a * bl + b' * c *al;
at=atan2(det,div);

omega = 2 * at;
