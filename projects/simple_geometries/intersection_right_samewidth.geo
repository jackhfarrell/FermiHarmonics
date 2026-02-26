// intersection_right_samewidth.geo
// 3-terminal T-intersection: main vertical channel width Wmain, side arm to the RIGHT with the SAME width.
// Equal arm length L from junction to each contact (top/bottom/side).

SetFactory("Built-in");

// -------------------- Parameters --------------------
Wmain = 0.4;     // width of straight channel
Warm  = Wmain;   // width of side channel (set equal to main)
L     = 1.0;     // arm length from junction to each contact
lc    = 0.03;    // target mesh size (tune as needed)

// Derived
xm    = Wmain/2.0;            // half-width of main channel
ya    = Warm/2.0;             // half-width of side arm
xside = xm + L;               // x-location of side contact plane (gives side length = L)

// -------------------- Geometry (single outer boundary polygon) --------------------
// Polygon vertices (counterclockwise):
// 1 (-xm,  -L)
// 2 ( xm,  -L)
// 3 ( xm, -ya)
// 4 (xside,-ya)
// 5 (xside, ya)
// 6 ( xm,  ya)
// 7 ( xm,   L)
// 8 (-xm,   L)

Point(1) = {-xm,   -L, 0, lc};
Point(2) = { xm,   -L, 0, lc};
Point(3) = { xm,  -ya, 0, lc};
Point(4) = {xside,-ya, 0, lc};
Point(5) = {xside, ya, 0, lc};
Point(6) = { xm,   ya, 0, lc};
Point(7) = { xm,    L, 0, lc};
Point(8) = {-xm,    L, 0, lc};

Line(1) = {1,2}; // bottom contact edge
Line(2) = {2,3}; // right wall (bottom segment)
Line(3) = {3,4}; // arm bottom wall
Line(4) = {4,5}; // side contact edge (right)
Line(5) = {5,6}; // arm top wall
Line(6) = {6,7}; // right wall (top segment)
Line(7) = {7,8}; // top contact edge
Line(8) = {8,1}; // left wall

Curve Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};

// -------------------- Physical groups --------------------
Physical Surface("domain") = {1};

Physical Curve("contact_bottom") = {1};
Physical Curve("contact_side")   = {4};
Physical Curve("contact_top")    = {7};

// Everything else is wall
Physical Curve("walls") = {2,3,5,6,8};
