
// global dimensions

L=2.2; // length of channel
H=.41; // height of channel
D=0.2;  // distance of obstacle
R=0.05;  // radius of obstacle

// mesh parameters

lc1=.025;
lc2=.01;

// lc1=.05;
// lc2=.02;

// outer rectangle

Point(1) = {0,0,0,lc1};
Point(2) = {L,0,0,lc1};
Point(3) = {L,H,0,lc1};
Point(4) = {0,H,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

// inner circel

Point(5) = {D,D,0,lc2};
Point(6) = {D+R,D,0,lc2};
Point(7) = {D,D+R,0,lc2};
Point(8) = {D-R,D,0,lc2};
Point(9) = {D,D-R,0,lc2};

Circle(6) = {6,5,7};
Circle(7) = {7,5,8};
Circle(8) = {8,5,9};
Circle(9) = {9,5,6};
Line Loop(10) = {-6,-7,-8,-9};

// surface

Plane Surface(1) = {5,10};

// physical entities

Physical Line(1) = {1,3,6,7,8,9}; // zero vel
Physical Line(2) = {4};   // inflow
Physical Surface(3) = {1};
Physical Line(4) = {2};   // outflow

// creates second order mesh and saves

Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'cylinder.msh';
