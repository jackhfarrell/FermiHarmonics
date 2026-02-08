SetFactory("OpenCASCADE");

// Characteristic length
lc = 0.02;

// Polygon vertices
Point(1)  = {-0.075, -0.5,    0, lc};
Point(2)  = { 0.075, -0.5,    0, lc};
Point(3)  = { 0.075, -0.075,  0, lc};
Point(4)  = { 0.5,   -0.075,  0, lc};
Point(5)  = { 0.5,    0.075,  0, lc};
Point(6)  = { 0.075,  0.075,  0, lc};
Point(7)  = { 0.075,  0.5,    0, lc};
Point(8)  = {-0.075,  0.5,    0, lc};
Point(9)  = {-0.075,  0.075,  0, lc};
Point(10) = {-0.075,  0.0375, 0, lc};
Point(11) = {-0.5,    0.0375, 0, lc};
Point(12) = {-0.5,   -0.0375, 0, lc};
Point(13) = {-0.075, -0.0375, 0, lc};

// Boundary lines
Line(1)  = {1,2};   // bottom
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};   // right
Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};   // middle
Line(8)  = {8,9};
Line(9)  = {9,10};
Line(10) = {10,11};
Line(11) = {11,12}; // left
Line(12) = {12,13};
Line(13) = {13,1};

Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("domain") = {1};
Physical Curve("bottom") = {1};
Physical Curve("right")  = {4};
Physical Curve("middle") = {7};
Physical Curve("left")   = {11};
Physical Curve("walls")  = {2,3,5,6,8,9,10,12,13};
