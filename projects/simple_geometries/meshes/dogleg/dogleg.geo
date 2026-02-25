// dogleg_rectstyle.geo
// L-shaped (dog-leg) channel with orthogonal contacts of equal width.
// Boundary naming: source, drain, walls, domain.
//
// Conventions:
//  - Physical Curve("source") is the inlet contact (horizontal segment)
//  - Physical Curve("drain")  is the outlet contact (vertical segment)
//  - Physical Curve("walls")  are all remaining boundaries
//  - Physical Surface("domain") is the fluid region

SetFactory("OpenCASCADE");

// ---------- Parameters ----------
W  = 0.30;   // channel / contact width
Lv = 1.20;   // vertical arm length (from source up to top wall)
Lh = 1.20;   // horizontal arm length (from left wall out to drain)
lc = 0.03;   // target mesh size (can be overridden by your meshing setup)

// Derived
y0 = 0.0;
y1 = Lv - W; // bottom of horizontal arm
y2 = Lv;     // top wall
x0 = 0.0;
x1 = W;      // right wall of vertical arm
x2 = Lh;     // drain wall

// ---------- Points ----------
Point(1) = {x0, y0, 0, lc}; // source left
Point(2) = {x1, y0, 0, lc}; // source right

Point(3) = {x1, y1, 0, lc}; // inner corner (concave)
Point(4) = {x2, y1, 0, lc}; // drain bottom

Point(5) = {x2, y2, 0, lc}; // drain top
Point(6) = {x0, y2, 0, lc}; // top-left

// ---------- Boundary curves ----------
// Contacts
Line(1) = {1, 2}; // source (horizontal)
Line(4) = {4, 5}; // drain  (vertical)

// Walls (all remaining boundaries)
Line(2) = {2, 3}; // right wall of vertical arm
Line(3) = {3, 4}; // bottom wall of horizontal arm
Line(5) = {5, 6}; // top wall
Line(6) = {6, 1}; // left wall

// ---------- Surface ----------
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// ---------- Physical groups ----------
Physical Curve("source") = {1};
Physical Curve("drain")  = {4};
Physical Curve("walls")  = {2, 3, 5, 6};
Physical Surface("domain") = {1};
