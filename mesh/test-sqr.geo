cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {10, 0, 0, 1};
Point(3) = {10, 10, 0, 1};
Point(4) = {0, 10, 0, 1};
Point(5) = {5, 5, 0, 1};
Point(6) = {6.7725, 6.7725, 0, 1};
Point(7) = {6.7725, 3.2275, 0, 1};
Point(8) = {3.2275, 3.2275, 0, 1};
Point(9) = {3.2275, 6.7725, 0, 1};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {8, 9};
Line(6) = {9, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line Loop(11) = {2, 3, 4, 1, -5, -8, -7, -6};
Plane Surface(11) = {11};
Physical Line("dirichlet 0") = {1, 2, 3, 4};
Physical Line("dirichlet 1") = {5, 6, 7, 8};
Physical Surface(14) = {11};
