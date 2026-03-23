// Tesla valve (single-stage, 2D)
// Generated from tesla_valve.control
// Forward (easy) flow direction: lower-left (inlet) → right (outlet)
//
// Channel width d = 0.12 as true perpendicular width throughout.
//
// Main channel:
//   RIGHT portion: horizontal, y ∈ [-0.06, +0.06]
//   LEFT  portion: 45° diagonal arm
//     Top wall:    y = x              (through (0.06, 0.06))
//     Bottom wall: y = x − 0.12√2    (through (0.1097, -0.06))
//     Perpendicular distance = 0.12 = d ✓
//   Inlet face: orthogonal to arm, from (-0.3400, -0.3400) to (-0.2551, -0.4249)
//
// Bypass loop: tilted ellipse arcs, apex up and to the right.
//   Ellipse centre C = (0.4382, 0.1051)
//   Major axis angle ≈ 34.4°, a = 0.4271, b = 0.2880
//   Outer arc: left (0.06, 0.06) → apex (0.5752, 0.4437) → right (0.78, 0.06)
//   Inner arc: left (0.1825, 0.06) → apex (0.5770, 0.3271) → right (0.6375, 0.06)
//   Each arc split into two segments through apex (avoids >180° arc in GMSH)
//   Outer major axis point: (0.7907, 0.3463)
//   Inner major axis point: (0.6917, 0.2785)
//
//   Left  junction gap: x ∈ [0.06, 0.1825], width ≈ 0.12 = d ✓
//   Right junction gap: x ∈ [0.6375, 0.78], width ≈ 0.12 = d ✓
//
// Both Curve Loops CCW; OpenCASCADE resolves fin as hole by containment.

SetFactory("OpenCASCADE");
lc = 0.03; // characteristic length; adjust as desired

// ── Ellipse reference points (not on any boundary curve) ─────────────────────
Point(20) = { 0.4382,  0.1051, 0, lc};  // shared ellipse centre
Point(21) = { 0.7907,  0.3463, 0, lc};  // outer ellipse major axis point
Point(22) = { 0.6917,  0.2785, 0, lc};  // inner ellipse major axis point

// ── Outer boundary points (CCW) ───────────────────────────────────────────────
Point(1) = {-0.2551, -0.4249, 0, lc};  // inlet face bottom
Point(2) = { 0.1097, -0.0600, 0, lc};  // 45° bottom wall meets horizontal bottom
Point(3) = { 0.7000, -0.0600, 0, lc};  // outlet bottom  (shifted right for ellipse)
Point(4) = { 0.9000, -0.0600, 0, lc};  // outlet bottom far right
Point(5) = { 0.9000,  0.0600, 0, lc};  // outlet top
Point(6) = { 0.7800,  0.0600, 0, lc};  // outer arc right endpoint
Point(7) = { 0.5752,  0.4437, 0, lc};  // outer arc apex
Point(8) = { 0.0600,  0.0600, 0, lc};  // outer arc left endpoint = 45° arm top corner
Point(9) = {-0.3400, -0.3400, 0, lc};  // inlet face top

// ── Fin (island) points (CCW) ─────────────────────────────────────────────────
Point(11) = {0.1825, 0.0600, 0, lc};  // inner arc left endpoint
Point(12) = {0.6375, 0.0600, 0, lc};  // inner arc right endpoint
Point(13) = {0.5770, 0.3271, 0, lc};  // inner arc apex

// ── Outer boundary curves (CCW) ───────────────────────────────────────────────
Line(1)      = {1, 2};         // 45° arm bottom wall
Line(2)      = {2, 4};         // horizontal bottom wall
Line(3)      = {4, 5};         // outlet (right wall)
Line(4)      = {5, 6};         // horizontal top wall, right of right junction
Ellipse(5)   = {6, 20, 21, 7}; // outer arc, right segment  (CCW: right → apex)
Ellipse(6)   = {7, 20, 21, 8}; // outer arc, left  segment  (CCW: apex → left)
Line(7)      = {8, 9};         // 45° arm top wall
Line(8)      = {9, 1};         // inlet face (orthogonal to arm)

// ── Fin boundary curves (CCW) ─────────────────────────────────────────────────
Line(9)      = {11, 12};        // fin bottom (left → right along y = 0.06)
Ellipse(10)  = {12, 20, 22, 13};// inner arc, right segment  (CCW: right → apex)
Ellipse(11)  = {13, 20, 22, 11};// inner arc, left  segment  (CCW: apex → left)

// ── Curve loops & surface ─────────────────────────────────────────────────────
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};  // outer boundary (CCW)
Curve Loop(2) = {9, 10, 11};                 // fin / hole     (CCW, contained in Loop 1)
Plane Surface(1) = {1, 2};

// ── Physical groups ───────────────────────────────────────────────────────────
Physical Surface("domain") = {1};
Physical Curve("inlet")  = {8};
Physical Curve("outlet") = {3};
Physical Curve("walls")  = {1, 2, 4, 5, 6, 7, 9, 10, 11};
