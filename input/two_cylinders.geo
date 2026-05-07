
// global dimensions

L = 40.0; // length of channel
H = 15.0; // height of channel
D = 5.0;  // horizontal coord. of the center of the cylinder
R = 0.5;  // radius of the cylinder
g = 0.74; // gap

// mesh parameters

lc1=.5;
lc2=.1;

// outer rectangle

Point(1) = {0,-0.5*H,0,lc1};
Point(2) = {L,-0.5*H,0,lc1};
Point(3) = {L,0.5*H,0,lc1};
Point(4) = {0,0.5*H,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

// down cylinder

Point(5) = {D,-0.5*g-R,0,lc2};
Point(6) = {D+R,-0.5*g-R,0,lc2};
Point(7) = {D,-0.5*g,0,lc2};
Point(8) = {D-R,-0.5*g-R,0,lc2};
Point(9) = {D,-0.5*g-2*R,0,lc2};

Circle(6) = {6,5,7};
Circle(7) = {7,5,8};
Circle(8) = {8,5,9};
Circle(9) = {9,5,6};
Line Loop(10) = {-6,-7,-8,-9};

// up cylinder

Point(10) = {D,0.5*g+R,0,lc2};
Point(11) = {D+R,0.5*g+R,0,lc2};
Point(12) = {D,0.5*g+2*R,0,lc2};
Point(13) = {D-R,0.5*g+R,0,lc2};
Point(14) = {D,0.5*g,0,lc2};

Circle(11) = {11,10,12};
Circle(12) = {12,10,13};
Circle(13) = {13,10,14};
Circle(14) = {14,10,11};
Line Loop(15) = {-11,-12,-13,-14};

// surface

Plane Surface(1) = {5,10,15};

// physical entities

Physical Line(1) = {1};            // down part
Physical Line(2) = {2};            // outflow
Physical Line(3) = {3};            // up part
Physical Line(4) = {4};            // inflow
Physical Line(5) = {6,7,8,9};      // down cylinder
Physical Line(6) = {11,12,13,14};  // up cylinder
Physical Surface(7) = {1};         // fluid

// creates second order mesh and saves

Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'two_cylinders.msh';
