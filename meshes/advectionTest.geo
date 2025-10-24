cl__1 = 0.4;
cl__2 = 0.6;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {5, 0, 0, cl__1};
Point(3) = {5, 2, 0, cl__1};
Point(4) = {0, 2, 0, cl__1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Curve {4, 2} = 10 Using Progression 1;
Transfinite Curve {3, 1} = 10 Using Progression 1;

Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Transfinite Surface {6};
Physical Line("Outflow",2) = {2};
Physical Line("Wall",1) = {1,3, 4};
Physical Surface("Domain") = {6};




