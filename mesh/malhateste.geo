cl__1 = 1;
Point(1) = {0, 0, 0, 1.5};
Point(2) = {5, 0, 0, 1.5};
Point(3) = {5, 5, 0, 1.5};
Point(4) = {0, 5, 0, 1.5};
Point(5) = {2, 2, 0, 0.4};
Point(6) = {3, 2, 0, 0.4};
Point(7) = {3, 3, 0, 0.4};
Point(8) = {2, 3, 0, 0.4};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {5, 8};
Line(6) = {8, 7};
Line(7) = {7, 6};
Line(8) = {6, 5};
Line Loop(11) = {1, 2, 3, 4, -8, -7, -6, -5};
Plane Surface(11) = {11};
Physical Line("neumann 0") = {1};
Physical Line("neumann 1") = {2};
Physical Line("neumann 2") = {3};
Physical Line("neumann 3") = {4};
Physical Line("aa") = {5, 6, 7, 8};
Physical Surface(17) = {11};
