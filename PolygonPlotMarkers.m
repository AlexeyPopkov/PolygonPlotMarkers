BeginPackage["PolygonPlotMarkers`"];

ClearAll[PolygonMarker];
PolygonMarker::usage="\!\(\*RowBox[{\"PolygonMarker\", \"[\", RowBox[{StyleBox[\"shape\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) returns Polygon of \!\(\*StyleBox[\"shape\", \"TI\"]\) with centroid at {0,0} and area \!\(\*SuperscriptBox[StyleBox[\"size\", \"TI\"], StyleBox[\"2\", \"TR\"]]\).";
SyntaxInformation[PolygonMarker]={"ArgumentsPattern"->{_,_.,_.}};

Begin["`Private`"];

ClearAll[PolygonArea,PolygonCentroid,LineIntersectionPoint,ngon,nstar,ncross,scale,coords];
(* The shoelace method for computing the area of polygon
http://mathematica.stackexchange.com/a/22587/280 *)
PolygonArea[pts_?MatrixQ]:=Abs@Total[Det/@Partition[pts,2,1,1]]/2;
(* http://mathematica.stackexchange.com/a/7715/280 *)
PolygonCentroid[pts_?MatrixQ]:=With[{dif=Map[Det,Partition[pts,2,1,{1,1}]]},ListConvolve[{{1,1}},Transpose[pts],{-1,-1}].dif/(3 Total[dif])];
(* http://mathematica.stackexchange.com/a/51399/280 *)
LineIntersectionPoint[{a_,b_},{c_,d_}]:=(Det[{a,b}] (c-d)-Det[{c,d}] (a-b))/Det[{a-b,c-d}];

ngon[n_,phase_:0]:=Table[{0,1}.RotationMatrix[2k Pi/n+phase],{k,0,n-1}];
(* 
  nn - number of vertices in related polygram
step - step at which vertices in the polygram are connected (must be lesser than nn/2)
n - number of points in the final star (must be divisor of nn) 
an illustration: http://en.wikipedia.org/wiki/Star_polygon# Simple _isotoxal _star _polygons
*)
nstar[n_/;n>=5,phase_:0]:=nstar[n,2,n,phase];
nstar[nn_,step_,n_,phase_:0]/;Divisible[nn,n]&&nn/2>step>nn/n:=Module[{a1,a2,b1,b2,ab},
{a1,a2,b1,b2}=ngon[nn][[{1,1+step,1+nn/n,nn/n-step}]];
ab=LineIntersectionPoint[{a1,a2},{b1,b2}];
Flatten[Table[{a1,ab}.RotationMatrix[2k Pi/n+phase],{k,0,n-1}],1]];
(* a - semiwidths of the crossing stripes *)
ncross[n_,phase_:0,a_:1/10]:=Flatten[NestList[#.RotationMatrix[2Pi/n]&,{{-a,1},{a,1},{a,a Cot[Pi/n]}}.RotationMatrix[phase],n-1],1];

(* Unitizes the area of the polygon *)
scale[coords_]:=Chop[#/Sqrt@PolygonArea@#]&@N[coords,{18,18}];

coords["UpTriangle"|"Triangle"]=ngon[3]//scale;
coords["DownTriangle"]=ngon[3,Pi/3]//scale;
coords["LeftTriangle"]=ngon[3,Pi/6]//scale;
coords["RightTriangle"]=ngon[3,-Pi/6]//scale;
coords["ThreePointedStar"]=nstar[12,5,3]//scale;
coords["DiagonalSquare"|"Diamond"]=ngon[4,0]//scale;
coords["Square"]=ngon[4,Pi/4]//scale;
coords["FourPointedStar"]=nstar[8,3,4]//scale;
coords["DiagonalFourPointedStar"]=nstar[8,3,4,Pi/4]//scale;
coords["Pentagon"]=ngon[5]//scale;
coords["FivePointedStar"]=nstar[5]//scale;
coords["FivePointedStarThick"]=nstar[20,7,5]//scale;
coords["Hexagon"]=ngon[6]//scale;
coords["SixPointedStar"]=nstar[6]//scale;
coords["SixPointedStarSlim"]=nstar[12,5,6]//scale;
coords["SevenPointedStar"]=nstar[7]//scale;
coords["SevenPointedStarNeat"]=nstar[14,5,7]//scale;
coords["SevenPointedStarSlim"]=nstar[14,6,7]//scale;
coords["Cross"]=ncross[4]//scale;
coords["DiagonalCross"|"CrossDiagonal"]=ncross[4,Pi/4]//scale;
coords["TripleCross"|"TripleCrossUp"]=ncross[3]//scale;
coords["TripleCrossDown"|"Y"]=ncross[3,Pi/3]//scale;
coords["FivefoldCross"]=ncross[5]//scale;
coords["SixfoldCross"]=ncross[6]//scale;
coords["SevenfoldCross"]=ncross[7]//scale;
coords["EightfoldCross"]=ncross[8]//scale;
(* The truncated triangle shape originates from the Cross's Theorem
http://demonstrations.wolfram.com/CrosssTheorem/
 *)
coords["UpTriangleTruncated"|"TriangleTruncated"|"TruncatedTriangle"]=Flatten[{{-3,6+Sqrt[3]},{3,6+Sqrt[3]}}.RotationMatrix[# Pi/3]&/@{0,2,4},1]//scale;
coords["DownTriangleTruncated"]=coords["UpTriangleTruncated"].ReflectionMatrix[{0,1}];
coords["LeftTriangleTruncated"]=coords["UpTriangleTruncated"].RotationMatrix[Pi/6];
coords["RightTriangleTruncated"]=coords["UpTriangleTruncated"].RotationMatrix[-Pi/6];
(* Disk approximated by 24-gon *)
coords["Disk"|"Circle"]=ngon[24]//scale;

(* Plotting symbols recommended in [Cleveland W.S. The Elements of Graphing Data (1985)] *)
(* Symmetric symbol "H" *)
coords["H"]=Join[#,-#]&@Join[#,Reverse@#.{{1,0},{0,-1}}]&@{{333,108},{333,630},{585,630}}//scale;
(* Symmetric symbol "I" *)
coords["I"]=Join[#,-#]&@{{-20,-68},{-64,-68},{-64,-104},{64,-104},{64,-68},{20,-68}}//scale;
(* Antisymmetric symbol "N" *)
coords["N"]=Join[#,-#]&@{{18,-32},{30,-32},{30,32},{17,32},{17,-12}}//scale;
(* Antisymmetric symbol "Z" *)
coords["Z"]=Join[#,-#]&@{{-567,-432},{-567,-630},{567,-630},{567,-414},{-234,-414}}//scale;
(* Antisymmetric symbol "S" (simple) *)
coords["S"]=Join[#,-#]&@{{-176,-54},{116,-54},{167,-100},{167,-170},{116,-216},{-284,-216},{-284,-324},{176,-324},{293,-216},{293,-54}}//scale;
(* Antisymmetric symbol "S" (curved, long) *)
coords["LongS"|"SLong"|"Sl"]=Join[#,-#]&@{-(29/32)},{-(49/16),-(3/11)},{-(425/91),23/28},{-(141/26),31/12},{-(165/32),88/19},{-(167/45),106/17},{-(24/17),149/21},{121/69,233/33},{130/27,31/5},{130/27,118/29},{127/47,199/39},{7/20,233/42},{-(12/7),139/26},{-(65/21),139/31},{-(395/113),114/35},{-(157/52),77/39},{-(83/44),56/41},{9/22,39/43}}//scale;

PolygonMarker[name_String,size_?NumericQ]:=Polygon[size coords[name]];
PolygonMarker[name_String,(h:Scaled|Offset)[size_?NumericQ]]:=Polygon[h[size #,{0,0}]&/@coords[name]];
PolygonMarker[coords:{{_?NumericQ,_?NumericQ}..},size_?NumericQ]:=Polygon[size N[scale[Transpose[Transpose[coords]-PolygonCentroid[coords]]],{16,16}]];
PolygonMarker[coords:{{_?NumericQ,_?NumericQ}..},Scaled[size_?NumericQ]]:=Polygon[Scaled[size #,{0,0}]&/@N[scale[Transpose[Transpose[coords]-PolygonCentroid[coords]]],{16,16}]];
PolygonMarker[arg:_String|{{_?NumericQ,_?NumericQ}..},size:_?NumericQ|(Scaled|Offset)[_?NumericQ],positions:{_?NumericQ,_?NumericQ}|{{_?NumericQ,_?NumericQ}..}]:=Translate[PolygonMarker[arg,size],positions];
PolygonMarker[]=PolygonMarker[All]={"TripleCross","Y","UpTriangle","UpTriangleTruncated","DownTriangle","DownTriangleTruncated","LeftTriangle","LeftTriangleTruncated","RightTriangle","RightTriangleTruncated","ThreePointedStar","Cross","DiagonalCross","Diamond","Square","FourPointedStar","DiagonalFourPointedStar","FivefoldCross","Pentagon","FivePointedStar","FivePointedStarThick","SixfoldCross","Hexagon","SixPointedStar","SixPointedStarSlim","SevenfoldCross","SevenPointedStar","SevenPointedStarNeat","SevenPointedStarSlim","EightfoldCross","Disk","H","I","N","Z","S","Sl"};
(* A subset of plot markers suitable for use when plotting symbols on the plot significantly overlap. *)
PolygonMarker["Overlap"]={"TripleCross","Y","UpTriangle","DownTriangle","LeftTriangle","RightTriangle","ThreePointedStar","Cross","DiagonalCross","Diamond","Square","FourPointedStar","DiagonalFourPointedStar","FivefoldCross","FivePointedStar","FivePointedStarThick","Disk","H","I","N","Z","S","Sl"};

End[];

EndPackage[];
