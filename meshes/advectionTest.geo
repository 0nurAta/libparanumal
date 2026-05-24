cl__1 = 0.6;
cl__2 = 0.6;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {4, 0, 0, cl__1};

Point(4) = {0, 1.5, 0, cl__1};

Point(5) = {8, 0, 0, cl__1};
Point(6) = {8, 1.5, 0, cl__1};

Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 4};
Line(4) = {4, 1};


Line(6) = {2, 5};
Line(7) = {6, 5};

Transfinite Curve {4, 2,7} = 20 Using Progression 1;
Transfinite Curve {3, 1} = 10 Using Progression 1;

Transfinite Curve {5, 6} = 50 Using Progression 1;


//Line Loop(6) = {4, 1, 2, 3};
//Plane Surface(6) = {6};

//Curve Loop(7) = {5, 7, -6, 2};
//Surface(7) = {7};





//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {6, -7, -2};
Surface(2) = {2};


Transfinite Surface {6,7};
Physical Line("Outflow",2) = {7};
Physical Line("Wall",1) = {1,3, 4,5,6};
Physical Surface("Domain") = {1,2};
