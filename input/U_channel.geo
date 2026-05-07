
// global dimensions

H  = 0.5;  // height of channel
A  = 1.0;  // size of the (square) hole
Ll = 1.2;  // length of the left part
Lr = 1.5;  // length of the right part

// mesh parameters

lc1 = .025;

// outer rectangle

Point(1) = {-Ll,0,0,lc1};
Point(2) = {-Ll/2,0,0,lc1};
Point(3) = {0,0,0,lc1};
Point(4) = {0,-A,0,lc1};
Point(5) = {A,-A,0,lc1};
Point(6) = {A,0,0,lc1};
Point(7) = {A+Lr/2,0,0,lc1};
Point(8) = {A+Lr,0,0,lc1};
Point(9) = {A+Lr,H,0,lc1};
Point(10) = {-Ll,H,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,1};
Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};

// surface

Plane Surface(1) = {11};

// physical entities

Physical Line(1) = {1,7};         // free slip
Physical Line(2) = {2,3,4,5,6};   // no slip
Physical Line(3) = {8};           // outflow
Physical Line(4) = {9};           // symmetry
Physical Line(5) = {10};          // inflow
Physical Surface(6) = {1};        

// creates second order mesh and saves

Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'U_channel.msh';


