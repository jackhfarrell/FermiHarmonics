// horseshoe.geo
// Fat single-bend (horseshoe/U-turn) channel with straight inlet/outlet contacts.
// Inlet: bottom of left leg. Outlet: bottom of right leg.

SetFactory("Built-in");

// Geometry parameters (edit freely)
W   = 0.22;    // channel width (fat)
Rc  = 0.28;    // centerline bend radius
L   = 0.55;    // straight-leg length below the bend
lc  = 0.03;    // mesh size hint

Rout = Rc + W/2;
Rin  = Rc - W/2;
If (Rin <= 0)
  Error("Choose Rc > W/2 so inner radius is positive.");
EndIf

// Key points
// Left leg (x < 0)
Point(1) = {-Rc + W/2, -L, 0, lc}; // inner left bottom
Point(2) = {-Rc - W/2, -L, 0, lc}; // outer left bottom
Point(3) = {-Rc - W/2,  0, 0, lc}; // outer left at bend

// Right leg (x > 0)
Point(4) = { Rc + W/2,  0, 0, lc}; // outer right at bend
Point(5) = { Rc + W/2, -L, 0, lc}; // outer right bottom
Point(6) = { Rc - W/2, -L, 0, lc}; // inner right bottom
Point(7) = { Rc - W/2,  0, 0, lc}; // inner right at bend
Point(8) = {-Rc + W/2,  0, 0, lc}; // inner left at bend

// Bend center
Point(9) = {0, 0, 0, lc};

// Boundary curves (counterclockwise loop)
Line(1)  = {1, 2};           // contact_in (bottom left)
Line(2)  = {2, 3};           // outer left wall
Circle(3)= {4, 9, 3};        // outer arc defined right->left (CCW goes through TOP)
Line(4)  = {4, 5};           // outer right wall
Line(5)  = {5, 6};           // contact_out (bottom right)
Line(6)  = {6, 7};           // inner right wall
Circle(7)= {7, 9, 8};        // inner arc right->left (CCW through TOP)
Line(8)  = {8, 1};           // inner left wall

// Surface
Curve Loop(1) = {1, 2, -3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Physical groups for boundary conditions
Physical Surface("domain")     = {1};
Physical Curve("contact_in")  = {1};
Physical Curve("contact_out") = {5};
Physical Curve("walls")       = {2, 3, 4, 6, 7, 8};
