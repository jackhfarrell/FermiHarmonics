SetFactory("OpenCASCADE");

length_x = 1.0;
half_height = 0.25;

nx = 33;
ny = 9;

Point(1) = {-length_x / 2, -half_height, 0, 0.08};
Point(2) = { length_x / 2, -half_height, 0, 0.08};
Point(3) = { length_x / 2,  half_height, 0, 0.08};
Point(4) = {-length_x / 2,  half_height, 0, 0.08};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1, 3} = nx Using Progression 1;
Transfinite Curve {2, 4} = ny Using Progression 1;
Transfinite Surface {1};

Mesh.RecombineAll = 1;
Recombine Surface {1};

Physical Surface("domain") = {1};
Physical Curve("inlet") = {4};
Physical Curve("outlet") = {2};
Physical Curve("walls") = {1, 3};
