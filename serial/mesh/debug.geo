cl__1 = 1;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {-0.5, 0, 0, cl__1};
Point(3) = {-0.25, 0, 0, cl__1};
Point(4) = {0, 0.25, 0, cl__1};
Point(5) = {0, 0.5, 0, cl__1};
Point(6) = {0.25, 0, 0, cl__1};
Point(7) = {0.5, 0, 0, cl__1};
Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 6};
Line(3) = {6, 7};
Circle(4) = {7, 1, 5};
Circle(5) = {5, 1, 2};
Line(6) = {2, 3};
Line Loop(7) = {1,2,3,4,5,6};
Plane Surface(7) = {7};
//+
Transfinite Surface {7} = {3, 6, 7, 2};
//+
Transfinite Curve {1, 2, 5, 4} = 6 Using Progression 1;
//+
Transfinite Curve {6, 3} = 11 Using Progression 1;
