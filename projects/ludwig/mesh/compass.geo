SetFactory("OpenCASCADE");

// Characteristic length
lc = 0.2;

// ============================================================
// 8-contact Hall-bar geometry (10 µm central square)
//
// Central square spans (-5, -5) to (5, 5).
// Cardinal contacts (N, E, S, W): narrow rectangular leads.
// Diagonal contacts (NE, SE, SW, NW): narrow rectangular leads
// at 45°, cutting into the corners.
// ============================================================

// --- Lead dimensions ---
cw     = 0.75;   // Cardinal lead half-width
cl     = 2.5;    // Cardinal lead length
dw     = cw;     // Diagonal lead half-width (same as cardinal)
dl     = 2.5;    // Diagonal lead length along 45° direction
hs     = 5.0;    // Half-size of central square

// Derived quantities
d45    = dl / Sqrt(2);           // lead length projected onto x or y
dw_e   = dw * Sqrt(2);          // half-width projected onto square edge

// ============================================================
// Points — CW traversal starting at S contact
//
// For each diagonal contact (e.g. SE at corner (hs,-hs)):
//   Inner pt on bottom edge : (hs - dw_e, -hs)
//   Inner pt on right edge  : (hs, -hs + dw_e)
//   Outer pt (from bottom)  : (hs - dw_e + d45, -hs - d45)
//   Outer pt (from right)   : (hs + d45, -hs + dw_e - d45)
// This gives a constant-width rectangular channel at 45°.
// ============================================================

// S contact (bottom, extends downward)
Point(1)  = {-cw,  -hs,      0, lc};
Point(2)  = {-cw,  -hs - cl, 0, lc};
Point(3)  = { cw,  -hs - cl, 0, lc};
Point(4)  = { cw,  -hs,      0, lc};

// Bottom edge to SE contact
Point(5)  = { hs - dw_e, -hs, 0, lc};

// SE contact (rectangular channel at 45° toward bottom-right)
Point(6)  = { hs - dw_e + d45, -hs - d45,        0, lc};
Point(7)  = { hs + d45,        -hs + dw_e - d45,  0, lc};
Point(8)  = { hs,  -hs + dw_e, 0, lc};

// Right edge to E contact
Point(9)  = { hs, -cw, 0, lc};

// E contact (extends rightward)
Point(10) = { hs + cl, -cw, 0, lc};
Point(11) = { hs + cl,  cw, 0, lc};
Point(12) = { hs,       cw, 0, lc};

// Right edge to NE contact
Point(13) = { hs,  hs - dw_e, 0, lc};

// NE contact (rectangular channel at 45° toward top-right)
Point(14) = { hs + d45,        hs - dw_e + d45, 0, lc};
Point(15) = { hs - dw_e + d45, hs + d45,        0, lc};
Point(16) = { hs - dw_e, hs, 0, lc};

// Top edge to N contact
Point(17) = { cw, hs, 0, lc};

// N contact (extends upward)
Point(18) = { cw, hs + cl, 0, lc};
Point(19) = {-cw, hs + cl, 0, lc};
Point(20) = {-cw, hs,      0, lc};

// Top edge to NW contact
Point(21) = {-hs + dw_e, hs, 0, lc};

// NW contact (rectangular channel at 45° toward top-left)
Point(22) = {-hs + dw_e - d45, hs + d45,        0, lc};
Point(23) = {-hs - d45,        hs - dw_e + d45, 0, lc};
Point(24) = {-hs, hs - dw_e, 0, lc};

// Left edge to W contact
Point(25) = {-hs, cw, 0, lc};

// W contact (extends leftward)
Point(26) = {-hs - cl,  cw, 0, lc};
Point(27) = {-hs - cl, -cw, 0, lc};
Point(28) = {-hs,      -cw, 0, lc};

// Left edge to SW contact
Point(29) = {-hs, -hs + dw_e, 0, lc};

// SW contact (rectangular channel at 45° toward bottom-left)
Point(30) = {-hs - d45,        -hs + dw_e - d45, 0, lc};
Point(31) = {-hs + dw_e - d45, -hs - d45,        0, lc};
Point(32) = {-hs + dw_e, -hs, 0, lc};

// ============================================================
// Boundary lines
// ============================================================

// S contact
Line(1)  = {1, 2};
Line(2)  = {2, 3};     // S contact edge
Line(3)  = {3, 4};

// Bottom edge: S -> SE
Line(4)  = {4, 5};

// SE contact
Line(5)  = {5, 6};
Line(6)  = {6, 7};     // SE contact edge
Line(7)  = {7, 8};

// Right edge: SE -> E
Line(8)  = {8, 9};

// E contact
Line(9)  = {9, 10};
Line(10) = {10, 11};   // E contact edge
Line(11) = {11, 12};

// Right edge: E -> NE
Line(12) = {12, 13};

// NE contact
Line(13) = {13, 14};
Line(14) = {14, 15};   // NE contact edge
Line(15) = {15, 16};

// Top edge: NE -> N
Line(16) = {16, 17};

// N contact
Line(17) = {17, 18};
Line(18) = {18, 19};   // N contact edge
Line(19) = {19, 20};

// Top edge: N -> NW
Line(20) = {20, 21};

// NW contact
Line(21) = {21, 22};
Line(22) = {22, 23};   // NW contact edge
Line(23) = {23, 24};

// Left edge: NW -> W
Line(24) = {24, 25};

// W contact
Line(25) = {25, 26};
Line(26) = {26, 27};   // W contact edge
Line(27) = {27, 28};

// Left edge: W -> SW
Line(28) = {28, 29};

// SW contact
Line(29) = {29, 30};
Line(30) = {30, 31};   // SW contact edge
Line(31) = {31, 32};

// Bottom edge: SW -> S
Line(32) = {32, 1};

// ============================================================
// Curve loop and surface
// ============================================================
Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
Plane Surface(1) = {1};

// ============================================================
// Physical groups
// ============================================================
Physical Surface("domain") = {1};

// Contacts (outer edges of each lead)
Physical Curve("S")  = {2};
Physical Curve("SE") = {6};
Physical Curve("E")  = {10};
Physical Curve("NE") = {14};
Physical Curve("N")  = {18};
Physical Curve("NW") = {22};
Physical Curve("W")  = {26};
Physical Curve("SW") = {30};

// All remaining boundary segments are walls
Physical Curve("walls") = {1,3,4,5,7,8,9,11,12,13,15,16,
                           17,19,20,21,23,24,25,27,28,29,31,32};
