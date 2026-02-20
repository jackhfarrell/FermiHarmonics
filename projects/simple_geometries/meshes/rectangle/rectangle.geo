SetFactory("OpenCASCADE");

// Characteristic length
lc = 0.02;

// Rectangle vertices (intersection.geo with the right arm removed)
Point(1) = {-0.075, -0.5, 0, lc};
Point(2) = { 0.075, -0.5, 0, lc};
Point(3) = { 0.075,  0.5, 0, lc};
Point(4) = {-0.075,  0.5, 0, lc};

// Boundary lines
Line(1) = {1,2};   // bottom contact
Line(2) = {2,3};   // right wall
Line(3) = {3,4};   // top contact
Line(4) = {4,1};   // left wall

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("domain") = {1};
Physical Curve("contact_bottom") = {1};
Physical Curve("contact_top")    = {3};
Physical Curve("walls")          = {2,4};
