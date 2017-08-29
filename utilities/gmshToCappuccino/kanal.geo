cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {2, 0, 0, 1};
Point(3) = {2, 0.02, 0, 1};
Point(4) = {0, 0.02, 0, 1};
Point(5) = {0, 0.02, 0.02, 1};
Point(6) = {0, 0, 0.02, 1};
Point(7) = {2, 0, 0.02, 1};
Point(8) = {2, 0.02, 0.02, 1};
//connect them with lines
Line(1) = {6, 5};
Line(2) = {6, 1};
Line(3) = {1, 4};
Line(4) = {4, 5};
Line(5) = {7, 8};
Line(6) = {8, 3};
Line(7) = {3, 2};
Line(8) = {2, 7};
Line(9) = {7, 6};
Line(10) = {5, 8};
Line(11) = {3, 4};
Line(12) = {1, 2};
//
Line Loop(13) = {1, -4, -3, -2};
Plane Surface(14) = {13};
Line Loop(15) = {9, 1, 10, -5};
Plane Surface(16) = {15};
Line Loop(17) = {8, 5, 6, 7};
Plane Surface(18) = {17};
Line Loop(19) = {6, 11, 4, 10};
Plane Surface(20) = {19};
Line Loop(21) = {2, 12, 8, 9};
Plane Surface(22) = {21};
Line Loop(23) = {7, -12, 3, -11};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};
//And the next lines will convert default generation strategy from tetrahedral meshes to hexagonal
Transfinite Line {1,2,3,4,5,6,7,8} = 41 Using Bump 0.25;
Transfinite Line {9,10,11,12} = 100;
//Transfinite Surface "*";
//Recombine Surface "*";
//Transfinite Volume "*";
//To use generated mesh in Cappuccino, one has to define physical groups also:
Physical Surface("inlet") = {14};
Physical Surface("outlet") = {18};
Physical Surface("wall") = {16, 20, 22, 24};
Physical Volume("kanal") = {26};
