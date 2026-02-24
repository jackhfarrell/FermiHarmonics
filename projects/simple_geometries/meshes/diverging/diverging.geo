SetFactory("OpenCASCADE");

// Mesh size
lc = 0.02;

// Total height = 1 (original style)
y0 = -0.5;
y1 =  0.5;

// Contact (throat) width and chamber (max) width
Wt = 0.12;     // same width at both contacts
Wc = 0.55;     // BIG divergence in the middle (tune: 0.45–0.70)

// Half widths
ht = Wt/2;
hc = Wc/2;

// Control point y-positions (symmetric about 0)
ya = -0.28;
ym =  0.00;
yb =  0.28;

// Optional shaping: push the bulge outward a bit more around midpoints
// (set s = 1.0 for "normal", >1 for stronger bow)
s = 1.10;
hmid = s*hc;

// --- Points ---
// Contacts
Point(1) = {-ht, y0, 0, lc};   // bottom-left
Point(2) = { ht, y0, 0, lc};   // bottom-right
Point(3) = { ht, y1, 0, lc};   // top-right
Point(4) = {-ht, y1, 0, lc};   // top-left

// Right wall control points (diverge then converge)
Point(10) = { ht,   y0, 0, lc};      // start (same as Point 2 position, but separate id ok)
Point(11) = { 0.35*hmid, ya, 0, lc}; // easing out
Point(12) = { hmid, ym, 0, lc};      // max bulge
Point(13) = { 0.35*hmid, yb, 0, lc}; // easing back
Point(14) = { ht,   y1, 0, lc};      // end (same as Point 3 position)

// Left wall control points (mirror)
Point(20) = {-ht,   y1, 0, lc};        // start (top)
Point(21) = {-0.35*hmid, yb, 0, lc};
Point(22) = {-hmid, ym, 0, lc};
Point(23) = {-0.35*hmid, ya, 0, lc};
Point(24) = {-ht,   y0, 0, lc};        // end (bottom)

// --- Curves ---
// Contacts (straight)
Line(1) = {1,2};   // contact_bottom
Line(2) = {3,4};   // contact_top (note orientation; physical group doesn't care)

// Curved walls (BSplines)
BSpline(3) = {10,11,12,13,14};  // right wall: bottom -> top
BSpline(4) = {20,21,22,23,24};  // left wall:  top -> bottom

// Close loop: go bottom contact -> right wall -> top contact -> left wall
Curve Loop(1) = {1, 3, 2, 4};
Plane Surface(1) = {1};

// Physical groups (match your rectangle.geo naming style)
Physical Surface("domain")        = {1};
Physical Curve("contact_bottom")  = {1};
Physical Curve("contact_top")     = {2};
Physical Curve("walls")           = {3,4};
