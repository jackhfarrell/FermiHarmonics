SetFactory("OpenCASCADE");

// Mesh size
lc = 0.02;

// Total height = 1 (normalized like rectangle.geo)
y0 = -0.5;
y1 =  0.5;

// Contact width (same at both contacts) and constriction (min) width in the middle
Wt = 0.12;     // contact width
Wc = 0.06;     // minimum width at the constriction (set < Wt)

// Half widths
ht = Wt/2;
hc = Wc/2;

// Control point y-positions (symmetric about 0)
ya = -0.28;
ym =  0.00;
yb =  0.28;

// Optional shaping: s=1.0 gives a gentle constriction, <1 makes it sharper
s = 1.00;
hmid = s*hc;

// --- Corner points (contacts) ---
Point(1) = {-ht, y0, 0, lc};   // bottom-left
Point(2) = { ht, y0, 0, lc};   // bottom-right
Point(3) = { ht, y1, 0, lc};   // top-right
Point(4) = {-ht, y1, 0, lc};   // top-left

// --- Right wall control points (converge then diverge) ---
// Endpoints are Point(2) -> Point(3), matching rectangle.geo topology
Point(11) = { ht - 0.35*(ht - hmid), ya, 0, lc}; // easing in
Point(12) = { hmid,                ym, 0, lc};  // tightest point
Point(13) = { ht - 0.35*(ht - hmid), yb, 0, lc}; // easing out

// --- Left wall control points (mirror) ---
// Endpoints are Point(4) -> Point(1), matching rectangle.geo topology
Point(21) = {-ht + 0.35*(ht - hmid), yb, 0, lc};
Point(22) = {-hmid,                ym, 0, lc};
Point(23) = {-ht + 0.35*(ht - hmid), ya, 0, lc};

// --- Boundary curves (match rectangle.geo numbering & names) ---
Line(1)    = {1,2};              // bottom contact
BSpline(2) = {2,11,12,13,3};      // right wall (bottom -> top)
Line(3)    = {3,4};              // top contact
BSpline(4) = {4,21,22,23,1};      // left wall (top -> bottom)

// Surface
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Physical groups (exactly like rectangle.geo)
Physical Surface("domain") = {1};
Physical Curve("contact_bottom") = {1};
Physical Curve("contact_top")    = {3};
Physical Curve("walls")          = {2,4};
