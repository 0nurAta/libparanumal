cl__1 = 0.6;
cl__2 = 0.6;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {4, 0, 0, cl__1};

Point(3) = {0, 2, 0, cl__1};
Point(4) = {4, 2, 0, cl__1};

Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {3, 1};
Line(4) = {4, 2};

Point(5) = {8, 0, 0, cl__1};
Point(6) = {8, 2, 0, cl__1};

Line(5) = {2, 5};
Line(6) = {4, 6};
Line(7) = {5, 6};
//Transfinite Curve {4, 2,7} = 20 Using Progression 1;
Transfinite Curve {2, 1} = 20 Using Progression 1;

Transfinite Curve {3} = 10 Using Progression 1;
Transfinite Curve {4} = 8 Using Progression 1;
Transfinite Curve {5,6} = 10 Using Progression 1;
Transfinite Curve {7} = 5 Using Progression 1;

//+
Curve Loop(1) = {3, 1, -4, -2};
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, 7, -6};
Plane Surface(2) = {2};

//Transfinite Surface {1};
Physical Line("Outflow",2) = {7};
Physical Line("Wall",1) = {1,2,3,5,6};
Physical Surface("Domain") = {1,2};


