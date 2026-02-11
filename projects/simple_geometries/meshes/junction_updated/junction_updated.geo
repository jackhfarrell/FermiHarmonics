SetFactory("OpenCASCADE");

// Characteristic length
lc = 0.02;

// Geometry parameters (in your current normalized units)
w_top   = 0.15;    // width of the top prong (and the vertical stem)
L_top   = 0.425;   // length of the top prong measured from the junction (y=+w_top/2) to the top contact

// Side prongs: half the width and half the length of the top prong
w_side  = 0.5*w_top;     // 0.075
L_side  = 0.5*L_top;     // 0.2125

// Convenience
x0 = 0.5*w_top;          // 0.075  (half-width of vertical stem)
y0 = 0.5*w_top;          // 0.075  (half-width of vertical stem)
ys = 0.5*w_side;         // 0.0375 (half-width of side prongs)
xr = x0 + L_side;        // 0.2875 (right prong end x)
xl = -xr;                // -0.2875 (left prong end x)
ybot = -0.5;
ytop =  y0 + L_top;      // 0.5

// Polygon vertices (counterclockwise)
Point(1)  = {-x0, ybot, 0, lc};
Point(2)  = { x0, ybot, 0, lc};
Point(3)  = { x0,-y0,   0, lc};

// Right prong (narrowed + shortened)
Point(4)  = { x0,-ys,   0, lc};
Point(5)  = { xr,-ys,   0, lc};
Point(6)  = { xr, ys,   0, lc};
Point(7)  = { x0, ys,   0, lc};

// Back to the vertical stem / top prong
Point(8)  = { x0, y0,   0, lc};
Point(9)  = { x0, ytop, 0, lc};
Point(10) = {-x0, ytop, 0, lc};
Point(11) = {-x0, y0,   0, lc};

// Left prong (matched to right: same width, half top prong width; same length, half top prong length)
Point(12) = {-x0, ys,   0, lc};
Point(13) = { xl, ys,   0, lc};
Point(14) = { xl,-ys,   0, lc};
Point(15) = {-x0,-ys,   0, lc};

// Close the loop back down the left side to the bottom
// (no extra corner needed; straight down from y=-ys to y=ybot)

// Boundary lines
Line(1)  = {1,2};    // bottom contact
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};
Line(5)  = {5,6};    // right contact
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,9};
Line(9)  = {9,10};   // top/middle contact
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};  // left contact
Line(14) = {14,15};
Line(15) = {15,1};

Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("domain") = {1};

Physical Curve("bottom") = {1};
Physical Curve("right")  = {5};
Physical Curve("middle") = {9};
Physical Curve("left")   = {13};

// Everything else is a wall
Physical Curve("walls")  = {2,3,4,6,7,8,10,11,12,14,15};
