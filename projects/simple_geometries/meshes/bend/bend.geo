SetFactory("OpenCASCADE");

// Mesh size
lc = 0.02;

// Normalized vertical extent (match rectangle.geo convention)
y0 = -0.5;
y1 =  0.5;
L  =  y1 - y0;

// Channel width (same at inlet/outlet contacts)
W  = 0.12;
h  = W/2;

// Lateral offset amplitude for the bend (set to 0 for straight channel)
A  = 0.18;

// Helper: smooth "single-bend" centerline offset with zero slope at contacts
// x_c(y) = A * sin^2(pi * (y - y0) / L)

// Sample y-levels for wall splines
ya = y0 + 0.20*L;
yb = y0 + 0.40*L;
ym = y0 + 0.50*L;
yc = y0 + 0.60*L;
yd = y0 + 0.80*L;

// Centerline offsets at sample points
xa = A * Sin(Pi*(ya - y0)/L)^2;
xb = A * Sin(Pi*(yb - y0)/L)^2;
xm = A * Sin(Pi*(ym - y0)/L)^2;
xc = A * Sin(Pi*(yc - y0)/L)^2;
xd = A * Sin(Pi*(yd - y0)/L)^2;

// --- Contact corner points (keep contacts straight & horizontal) ---
// Bottom contact is centered at x=0; top contact also centered at x=0
Point(1) = {-h, y0, 0, lc};   // bottom-left
Point(2) = { h, y0, 0, lc};   // bottom-right
Point(3) = { h, y1, 0, lc};   // top-right
Point(4) = {-h, y1, 0, lc};   // top-left

// --- Right wall control points (bottom -> top) ---
Point(11) = {xa + h, ya, 0, lc};
Point(12) = {xb + h, yb, 0, lc};
Point(13) = {xm + h, ym, 0, lc};
Point(14) = {xc + h, yc, 0, lc};
Point(15) = {xd + h, yd, 0, lc};

// --- Left wall control points (top -> bottom), same centerline but minus half-width ---
Point(21) = {xd - h, yd, 0, lc};
Point(22) = {xc - h, yc, 0, lc};
Point(23) = {xm - h, ym, 0, lc};
Point(24) = {xb - h, yb, 0, lc};
Point(25) = {xa - h, ya, 0, lc};

// --- Boundary curves (match rectangle.geo numbering & names) ---
Line(1)    = {1,2};                        // bottom contact
BSpline(2) = {2,11,12,13,14,15,3};          // right wall (bottom -> top)
Line(3)    = {3,4};                        // top contact
BSpline(4) = {4,21,22,23,24,25,1};          // left wall (top -> bottom)

// Surface
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Physical groups (exactly like rectangle.geo)
Physical Surface("domain") = {1};
Physical Curve("contact_bottom") = {1};
Physical Curve("contact_top")    = {3};
Physical Curve("walls")          = {2,4};
