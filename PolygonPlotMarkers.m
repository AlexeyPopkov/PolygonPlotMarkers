BeginPackage["PolygonPlotMarkers`"];

ClearAll[PolygonMarker];
PolygonMarker::usage=ToString[Row[{"PolygonMarker","[",Style["shape","TI"],",",Style["size","TI"],"]"," returns ","Polygon"," describing ", Style["shape","TI"]," with centroid at ","{0,0}"," and area ",Superscript[Style["size","TI"],2],".\n","PolygonMarker[All]"," returns the list of supported ",Style["shape","TI"]," names."}],StandardForm];
SyntaxInformation[PolygonMarker]={"ArgumentsPattern"->{_,_,_.}};

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
  
  an illustration: http://en.wikipedia.org/wiki/Star_polygon# Simple_isotoxal _star _polygons
*)
nstar[n_/;n>=5,phase_:0]:=nstar[n,2,n,phase];
nstar[nn_,step_,n_,phase_:0]/;Divisible[nn,n]&&nn/2>step>nn/n:=Module[{a1,a2,b1,b2,ab},
{a1,a2,b1,b2}=ngon[nn][[{1,1+step,1+nn/n,nn/n-step}]];
ab=LineIntersectionPoint[{a1,a2},{b1,b2}];
Flatten[Table[{a1,ab}.RotationMatrix[2k Pi/n+phase],{k,0,n-1}],1]];
(* a - semiwidths of the crossing stripes *)
ncross[n_,phase_:0,a_:1/10]:=Flatten[NestList[#.RotationMatrix[2Pi/n]&,{{-a,1},{a,1},{a,a Cot[Pi/n]}}.RotationMatrix[phase],n-1],1];

(* Unitizes the area of the polygon *)
scale[coords_]:=Chop[#/Sqrt@PolygonArea@#]&@N[coords,{16,16}];

coords["UpTriangle"|"Triangle"]=ngon[3]//scale;
coords["DownTriangle"]=ngon[3,Pi/3]//scale;
coords["LeftTriangle"]=ngon[3,Pi/6]//scale;
coords["RightTriangle"]=ngon[3,-Pi/6]//scale;
coords["ThreePointedStar"]=nstar[12,5,3]//scale;
coords["DiagonalSquare"|"Diamond"]=ngon[4,0]//scale;
coords["Square"]=ngon[4,Pi/4]//scale;
coords["FourPointedStar"]=nstar[8,3,4]//scale;
coords["DiagonalFourPointedStar"] = nstar[8, 3, 4, Pi/4] // scale;
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
coords["TripleCross"]=ncross[3]//scale;
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
(* Circle approximated by 24-gon *)
coords["Circle" | "Disk"] = ngon[24] // scale;

PolygonMarker[name_String,size_?NumericQ]:=Polygon[size coords[name]];
PolygonMarker[name_String,(h:Scaled|Offset)[size_?NumericQ]]:=Polygon[h[size #,{0,0}]&/@coords[name]];
PolygonMarker[coords:{{_?NumericQ,_?NumericQ}..},size_?NumericQ]:=Polygon[size N[scale[Transpose[Transpose[coords]-PolygonCentroid[coords]]],{16,16}]];
PolygonMarker[coords:{{_?NumericQ,_?NumericQ}..},Scaled[size_?NumericQ]]:=Polygon[Scaled[size #,{0,0}]&/@N[scale[Transpose[Transpose[coords]-PolygonCentroid[coords]]],{16,16}]];
PolygonMarker[arg:_String|{{_?NumericQ,_?NumericQ}..},size:_?NumericQ|(Scaled|Offset)[_?NumericQ],positions:{_?NumericQ,_?NumericQ}|{{_?NumericQ,_?NumericQ}..}]:=Translate[PolygonMarker[arg,size],positions];
PolygonMarker[]=PolygonMarker[All]={"TripleCross","UpTriangle","UpTriangleTruncated","DownTriangle","DownTriangleTruncated","LeftTriangle","LeftTriangleTruncated","RightTriangle","RightTriangleTruncated","ThreePointedStar","Cross","DiagonalCross","Diamond","Square","FourPointedStar","DiagonalFourPointedStar","FivefoldCross","Pentagon","FivePointedStar","FivePointedStarThick","SixfoldCross","Hexagon","SixPointedStar","SixPointedStarSlim","SevenfoldCross","SevenPointedStar","SevenPointedStarNeat","SevenPointedStarSlim","EightfoldCross", "Circle"};

End[];

EndPackage[];
