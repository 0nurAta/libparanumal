//+
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Recombine Surface {1};
//+
Sphere(2) = {0, 0.540, -0, 0.045, -Pi/2, Pi/2, 2*Pi};
//+

//+

//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
