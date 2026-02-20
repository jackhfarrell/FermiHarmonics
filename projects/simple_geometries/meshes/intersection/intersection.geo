SetFactory("OpenCASCADE");

// Characteristic length
lc = 0.02;

// Polygon vertices (junction.geo with left arm removed)
Point(1)  = {-0.075, -0.5,    0, lc};
Point(2)  = { 0.075, -0.5,    0, lc};
Point(3)  = { 0.075, -0.075,  0, lc};
Point(4)  = { 0.5,   -0.075,  0, lc};
Point(5)  = { 0.5,    0.075,  0, lc};
Point(6)  = { 0.075,  0.075,  0, lc};
Point(7)  = { 0.075,  0.5,    0, lc};
Point(8)  = {-0.075,  0.5,    0, lc};

// Boundary lines
Line(1)  = {1,2};   // bottom contact
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};   // right contact
Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};   // top contact
Line(8)  = {8,1};

Curve Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("domain") = {1};
Physical Curve("contact_bottom") = {1};
Physical Curve("contact_right")  = {4};
Physical Curve("contact_top")    = {7};
Physical Curve("walls")          = {2,3,5,6,8};
