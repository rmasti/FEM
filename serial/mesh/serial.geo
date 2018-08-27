cl__1 = 1;
Point(1) = {-0.5, 0, 0, 1};
Point(2) = {-0.25, 0, 0, 1};
Point(3) = {0.25, 0, 0, 1};
Point(4) = {0.5, 0, 0, 1};
Point(5) = {0, 0, 0, 1};
Circle(1) = {3, 5, 2};
Circle(2) = {4, 5, 1};
Line(3) = {1, 2};
Line(4) = {3, 4};
Line Loop(6) = {2, 3, -1, 4};
Plane Surface(6) = {6};
Transfinite Surface {6} = {3, 4, 1, 2};
Transfinite Line {3, 4} = 65 Using Progression 1;
Transfinite Line {1, 2} = 384 Using Progression 1;
Delete {
  Point{5};
}
Delete {
  Point{5};
}
Delete {
  Point{5};
}
