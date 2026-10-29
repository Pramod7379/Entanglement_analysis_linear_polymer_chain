function omega = solid_angle_quadrilateral(a, b, c, d)
% This function takes four vectors a, b, c and d, and produces solid angle 
% made by these four vectors.  

omega = solid_angle_triangle(a, b, c) + solid_angle_triangle(c, d, a);
