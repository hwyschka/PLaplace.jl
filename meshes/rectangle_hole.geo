// Gmsh project created on Mon Sep 12 10:52:29 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {2, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {1, 0.5, 0, 0.15, 0, 2*Pi};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("bottom", 1001) = {1};
//+
Physical Curve("top", 1002) = {3};
//+
Physical Curve("left", 1003) = {4};
//+
Physical Curve("right", 1004) = {2};
//+
Physical Curve("hole", 1005) = {5};
//+
Physical Surface("domain", 10001) = {1};
