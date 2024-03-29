BeginPackage["PolygonPlotMarkers`"];

ClearAll[PolygonMarker];
PolygonMarker::usage="\!\(\*RowBox[{\"PolygonMarker\", \"[\",StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"name\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True], \"]\"}]\) returns a unit area Polygon describing the shape \!\(\*StyleBox[\"\\\"\\!\\(\\*StyleBox[\\\"name\\\",\\\"TI\\\"]\\)\\\"\", ShowStringCharacters->True]\) with centroid at {0,0}.\n\!\(\*RowBox[{\"PolygonMarker\", \"[\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"p\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", \ StyleBox[\"\[Ellipsis]\", \"TR\"], \",\", SubscriptBox[StyleBox[\"p\", \"TI\"], StyleBox[\"n\", \"TI\"]]}], \"}\"}], \"]\"}]\) returns a unit area Polygon with shape described by points \!\(\*SubscriptBox[StyleBox[\"p\", \"TI\"], StyleBox[\"i\", \"TI\"]]\) and centroid at {0,0}.\n\!\(\*RowBox[{\"PolygonMarker\", \"[\", RowBox[{StyleBox[\"shape\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) returns Polygon of \!\(\*StyleBox[\"shape\", \"TI\"]\) with centroid at {0,0} and area \!\(\*SuperscriptBox[StyleBox[\"size\", \"TI\"], StyleBox[\"2\", \"TR\"]]\).\n\!\(\*RowBox[{\"PolygonMarker\", \"[\", RowBox[{StyleBox[\"shape\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"style\", \"TI\"]}], \"]\"}]\) returns a Graphics object which can be used as a marker for PlotMarkers where the style of \!\(\*StyleBox[\"shape\", \"TI\"]\) is defined by \!\(\*StyleBox[\"style\", \"TI\"]\).\n\!\(\*RowBox[{\"PolygonMarker\", \"[\",\"All\", \"]\"}]\) returns the list of names of predefined shapes.";
SyntaxInformation[PolygonMarker]={"ArgumentsPattern"->{_,_.,_.,OptionsPattern[]}};
Options[PolygonMarker] = {AlignmentPoint -> {0,0}, BaselinePosition -> Axis, AspectRatio -> Automatic, Axes -> False, AxesLabel -> None, AxesOrigin -> {0,0}, AxesStyle -> {}, Background -> None, BaseStyle -> {},  ContentSelectable -> Automatic, CoordinatesToolOptions -> Automatic, DisplayFunction :> Identity, Epilog -> {}, FormatType :> TraditionalForm, Frame -> False, FrameLabel -> None, FrameStyle -> {}, FrameTicks -> Automatic, FrameTicksStyle -> {}, GridLines -> None, GridLinesStyle -> {}, ImageMargins -> 0., ImagePadding -> All, ImageSize -> Automatic, ImageSizeRaw -> Automatic, LabelStyle -> {}, Method -> Automatic, PlotLabel -> None, PlotRange -> All, PlotRangeClipping -> False, PlotRangePadding -> Automatic, PlotRegion -> Automatic, PreserveImageOptions -> Automatic, Prolog -> {}, RotateLabel -> True, Ticks -> Automatic, TicksStyle -> {}};

Begin["`Private`"];

ClearAll[PolygonArea,PolygonCentroid,LineIntersectionPoint,ngon,nstar,ncross,scale,coords,dancingStar];
(*The shoelace method for computing the area of polygon http://mathematica.stackexchange.com/a/22587/280*)
PolygonArea[pts_?MatrixQ]:=Abs@Total[Det/@Partition[pts,2,1,1]]/2;
(*http://mathematica.stackexchange.com/a/7715/280*)
PolygonCentroid[pts_?MatrixQ]:=With[{dif=Map[Det,Partition[pts,2,1,{1,1}]]},ListConvolve[{{1,1}},Transpose[pts],{-1,-1}] . dif/(3 Total[dif])];
(*http://mathematica.stackexchange.com/a/51399/280*)
LineIntersectionPoint[{a_,b_},{c_,d_}]:=(Det[{a,b}] (c-d)-Det[{c,d}] (a-b))/Det[{a-b,c-d}];

ngon[n_,phase_:0]:=Table[{0,1} . RotationMatrix[2k Pi/n+phase],{k,0,n-1}];
(* 
   nn - number of vertices in related polygram
   step - step at which vertices in the polygram are connected (must be lesser than nn/2)
   n - number of points in the final star (must be divisor of nn)  an illustration: http://en.wikipedia.org/wiki/Star_polygon# Simple_isotoxal _star _polygons
*)
nstar[n_/;n>=5,phase_:0]:=nstar[n,2,n,phase];
nstar[nn_,step_,n_,phase_:0]/;Divisible[nn,n]&&nn/2>step>nn/n:=Module[{a1,a2,b1,b2,ab},
   {a1,a2,b1,b2}=ngon[nn][[{1,1+step,1+nn/n,nn/n-step}]];
   ab=LineIntersectionPoint[{a1,a2},{b1,b2}];
   Flatten[Table[{a1,ab} . RotationMatrix[2k Pi/n+phase],{k,0,n-1}],1]];
(*a-semiwidths of the crossing stripes*)
ncross[n_,phase_:0,a_:1/10]:=Flatten[NestList[# . RotationMatrix[2Pi/n]&,{{-a,1},{a,1},{a,a Cot[Pi/n]}} . RotationMatrix[phase],n-1],1];

dancingStarRight[n_,k_]:=Join@@Table[RotationTransform[2 i Pi/n,{0,0}][{{0,Cos[(2 \[Pi])/k]},{0,1}}],{i,0,k-1}];
dancingStarLeft[n_,k_]:=Join@@Table[RotationTransform[2 i Pi/n,{0,0}][{{0,1},{0,Cos[(2 \[Pi])/k]}}],{i,0,k-1}];

(*Unitizes the area of the polygon*)
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
coords["FivePointedStarSlim"]=nstar[40,18,5]//scale;
coords["Hexagon"]=ngon[6]//scale;
coords["SixPointedStar"]=nstar[6]//scale;
coords["SixPointedStarSlim"]=nstar[12,5,6]//scale;
coords["SevenPointedStar"]=nstar[7]//scale;
coords["SevenPointedStarNeat"]=nstar[14,5,7]//scale;
coords["SevenPointedStarSlim"]=nstar[14,6,7]//scale;
coords["Cross"|"+"]=ncross[4]//scale;
coords["DiagonalCross"|"CrossDiagonal"|"X"|"x"]=ncross[4,Pi/4]//scale;
coords["TripleCross"|"TripleCrossUp"]=ncross[3]//scale;
coords["TripleCrossDown"|"Y"|"y"]=ncross[3,Pi/3]//scale;
coords["FivefoldCross"]=ncross[5]//scale;
coords["SixfoldCross"]=ncross[6]//scale;
coords["SevenfoldCross"]=ncross[7]//scale;
coords["EightfoldCross"]=ncross[8]//scale;
(*The truncated triangle shape originates from the Cross's Theorem http://demonstrations.wolfram.com/CrosssTheorem/*)
coords["UpTriangleTruncated"|"TriangleTruncated"|"TruncatedTriangle"]=Flatten[{{-3,6+Sqrt[3]},{3,6+Sqrt[3]}} . RotationMatrix[# Pi/3]&/@{0,2,4},1]//scale;
coords["DownTriangleTruncated"]=coords["UpTriangleTruncated"] . ReflectionMatrix[{0,1}];
coords["LeftTriangleTruncated"]=coords["UpTriangleTruncated"] . RotationMatrix[Pi/6];
coords["RightTriangleTruncated"]=coords["UpTriangleTruncated"] . RotationMatrix[-Pi/6];
(*Disk approximated by 24-gon*)
coords["Disk"|"Circle"]=ngon[24]//scale;
(*Dancing stars*)
coords["DancingStar"|"DancingStarLeft"|"FivefoldPinwheel"|"FivefoldPinwheelLeft"]=dancingStarLeft[5,5]//scale;
coords["DancingStarRight"|"FivefoldPinwheel"|"FivefoldPinwheelRight"]=dancingStarRight[5,5]//scale;
coords["DancingStarThick"|"DancingStarThickLeft"]=dancingStarLeft[5,6]//scale;
coords["DancingStarThickRight"]=dancingStarRight[5,6]//scale;
coords["SixfoldPinwheel"|"SixfoldPinwheelLeft"]=dancingStarLeft[6,6]//scale;
coords["SixfoldPinwheelRight"]=dancingStarRight[6,6]//scale;
coords["SevenfoldPinwheel"|"SevenfoldPinwheelLeft"]=dancingStarLeft[7,7]//scale;
coords["SevenfoldPinwheelRight"]=dancingStarRight[7,7]//scale;

(*Plotting symbols recommended in[Cleveland W.S.The Elements of Graphing Data (1985)]*)
(*Symmetric symbol "H"*)
coords["H"]=Join[#,-#]&@Join[#,Reverse@# . {{1,0},{0,-1}}]&@{{333,108},{333,630},{585,630}}//scale;
(*Symmetric symbol "I"*)
coords["I"]=Join[#,-#]&@{{-20,-68},{-64,-68},{-64,-104},{64,-104},{64,-68},{20,-68}}//scale;
(*Antisymmetric symbol "N"*)
coords["N"]=Join[#,-#]&@{{18,-32},{30,-32},{30,32},{17,32},{17,-12}}//scale;
(*Antisymmetric symbol "Z"*)
coords["Z"]=Join[#,-#]&@{{-567,-432},{-567,-630},{567,-630},{567,-414},{-234,-414}}//scale;
(*Antisymmetric symbol "S" (simple)*)
coords["S"]=Join[#,-#]&@{{-176,-54},{116,-54},{167,-100},{167,-170},{116,-216},{-284,-216},{-284,-324},{176,-324},{293,-216},{293,-54}}//scale;
(*Antisymmetric symbol "S" (curved,long)*)
coords["LongS"|"SLong"|"Sl"]=Join[#,-#]&@{{-(49/16),-(3/11)},{-(425/91),23/28},{-(141/26),31/12},{-(165/32),88/19},{-(167/45),106/17},{-(24/17),149/21},{121/69,233/33},{130/27,31/5},{130/27,118/29},{127/47,199/39},{7/20,233/42},{-(12/7),139/26},{-(65/21),139/31},{-(395/113),114/35},{-(157/52),77/39},{-(83/44),56/41},{9/22,39/43}}//scale;
(*Antisymmetric symbol "S" (curved,wide)*)
coords["WideS"|"SWide"|"Sw"]=Join[#,-#]&@{{80/11,-(3/5)},{49/6,-(9/4)},{97/12,-(41/11)},{39/5,-(35/8)},{88/13,-(65/12)},{51/10,-(49/8)},{2,-(13/2)},{-(20/11),-(13/2)},{-(37/8),-(81/13)},{-(81/13),-(40/7)},{-(59/8),-(54/11)},{-(81/10),-(26/7)},{-(70/11),-(29/9)},{-(57/11),-(46/11)},{-(11/4),-(33/7)},{11/7,-(19/4)},{16/3,-(37/9)},{31/5,-(38/11)},{32/5,-(38/13)},{37/6,-(49/24)},{61/13,-(6/5)},{23/7,-(13/14)},{-(25/9),-(4/5)},{-(23/4),-(3/13)}}//scale;

PolygonMarker[All]=PolygonMarker[]={"TripleCross", "Y", "UpTriangle", "UpTriangleTruncated", "DownTriangle", "DownTriangleTruncated", "LeftTriangle", "LeftTriangleTruncated", "RightTriangle", "RightTriangleTruncated", "ThreePointedStar", "Cross", "DiagonalCross", "Diamond", "Square", "FourPointedStar", "DiagonalFourPointedStar", "FivefoldCross", "Pentagon", "FivePointedStar","FivePointedStarSlim", "FivePointedStarThick","DancingStar","DancingStarRight", "DancingStarThick","DancingStarThickRight","SixfoldCross", "Hexagon", "SixPointedStar", "SixPointedStarSlim","SixfoldPinwheel","SixfoldPinwheelRight", "SevenfoldCross", "SevenPointedStar", "SevenPointedStarNeat", "SevenPointedStarSlim","SevenfoldPinwheel","SevenfoldPinwheelRight", "EightfoldCross", "Disk", "H", "I", "N", "Z", "S", "Sw", "Sl"};

PolygonMarker["Overlap"] = {"TripleCross", "Y", "UpTriangle", "DownTriangle", "LeftTriangle", "RightTriangle", "ThreePointedStar", "Cross", "DiagonalCross", "Diamond", "Square", "FourPointedStar", "DiagonalFourPointedStar", "FivefoldCross", "FivePointedStar","FivePointedStarSlim", "FivePointedStarThick","DancingStar","DancingStarRight", "Disk", "H", "I", "N", "Z", "S", "Sw", "Sl"};
 
PolygonMarker[name_String] := Polygon[coords[name]];
 
PolygonMarker[name_String, size_?NumericQ] := Polygon[size*coords[name]];

PolygonMarker[name_String, {size_?NumericQ,angle_?NumericQ}] := Polygon[size*RotationTransform[angle][coords[name]]];
 
PolygonMarker[name_String, (h:Scaled | Offset)[size_?NumericQ]] := Polygon[(h[size*#1, {0, 0}] & ) /@ coords[name]];

PolygonMarker[name_String, {(h:Scaled | Offset)[size_?NumericQ],angle_?NumericQ}] := Polygon[(h[size*#1, {0, 0}] & ) /@ RotationTransform[angle][coords[name]]];
 
PolygonMarker[coords:{{_?NumericQ, _?NumericQ}..}] := PolygonMarker[coords, 1];
 
PolygonMarker[coords:{{_?NumericQ, _?NumericQ}..}, size_?NumericQ] := Polygon[size*N[scale[Transpose[Transpose[coords] - PolygonCentroid[coords]]], {16, 16}]];

PolygonMarker[coords:{{_?NumericQ, _?NumericQ}..}, {size_?NumericQ,angle_?NumericQ}]:=Polygon[size*RotationTransform[angle][N[scale[Transpose[Transpose[coords] - PolygonCentroid[coords]]], {16, 16}]]];
 
PolygonMarker[coords:{{_?NumericQ, _?NumericQ}..}, (h:Scaled | Offset)[size_?NumericQ]] := Polygon[(h[size*#1, {0, 0}] & ) /@ N[scale[Transpose[Transpose[coords] - PolygonCentroid[coords]]], {16, 16}]];

PolygonMarker[coords:{{_?NumericQ, _?NumericQ}..}, {(h:Scaled | Offset)[size_?NumericQ],angle_?NumericQ}]:=Polygon[(h[size*#1, {0, 0}] & ) /@RotationTransform[angle][N[scale[Transpose[Transpose[coords] - PolygonCentroid[coords]]], {16, 16}]]];
 
PolygonMarker[arg:_String | {{_?NumericQ, _?NumericQ}..}, size:_?NumericQ , position:{_?NumericQ, _?NumericQ} ] := TranslationTransform[position]/@PolygonMarker[arg, size];

PolygonMarker[arg:_String | {{_?NumericQ, _?NumericQ}..}, size:_?NumericQ , positions:{{_?NumericQ, _?NumericQ}..} ] := Table[PolygonMarker[arg, size,position],{position,positions}];

PolygonMarker[arg:_String | {{_?NumericQ, _?NumericQ}..}, size:(Scaled | Offset)[_?NumericQ], positions:{_?NumericQ, _?NumericQ} | {{_?NumericQ, _?NumericQ}..}] := Translate[PolygonMarker[arg, size], positions];

PolygonMarker[arg:_String | {{_?NumericQ, _?NumericQ}..},{size:_?NumericQ ,angle_?NumericQ}, position:{_?NumericQ, _?NumericQ}] := TranslationTransform[position]/@PolygonMarker[arg, {size,angle}];

PolygonMarker[arg:_String | {{_?NumericQ, _?NumericQ}..}, {size:_?NumericQ ,angle_?NumericQ} , positions:{{_?NumericQ, _?NumericQ}..} ] := Table[PolygonMarker[arg, {size,angle},position],{position,positions}];

PolygonMarker[arg:_String | {{_?NumericQ, _?NumericQ}..},{size:(Scaled | Offset)[_?NumericQ],angle_?NumericQ}, positions:{_?NumericQ, _?NumericQ} | {{_?NumericQ, _?NumericQ}..}] := Translate[PolygonMarker[arg, {size,angle}], positions];
 
PolygonMarker[shape_, size:_?NumericQ | (Scaled | Offset)[_?NumericQ], {{g___}, {primitives___}}] := Block[{p = PolygonMarker[shape, size]}, Graphics[{{g, p}, {primitives}}, AlignmentPoint -> {0, 0}, ImagePadding -> All, PlotRange -> All] /; Head[p] === Polygon];
(* This form allows to construct composite plot markers containing additional graphics primitives *)
PolygonMarker[shape_, {size:_?NumericQ | (Scaled | Offset)[_?NumericQ],angle_?NumericQ}, {{g___}, {primitives___}}] := Block[{p = PolygonMarker[shape, {size,angle}]}, Graphics[{{g, p}, {primitives}}, AlignmentPoint -> {0, 0}, ImagePadding -> All, PlotRange -> All] /; Head[p] === Polygon];
 PolygonMarker[shape_, size:_?NumericQ | (Scaled | Offset)[_?NumericQ], {g___}] := Block[{p = PolygonMarker[shape, size]}, Graphics[{g, p}, AlignmentPoint -> {0, 0}, ImagePadding -> All, PlotRange -> All] /; Head[p] === Polygon];
PolygonMarker[shape_, {size:_?NumericQ | (Scaled | Offset)[_?NumericQ],angle_?NumericQ}, {g___}] := Block[{p = PolygonMarker[shape, {size,angle}]}, Graphics[{g, p}, AlignmentPoint -> {0, 0}, ImagePadding -> All, PlotRange -> All] /; Head[p] === Polygon];
 PolygonMarker[shape_, size:_?NumericQ | (Scaled | Offset)[_?NumericQ], g_]:=PolygonMarker[shape, size, {g}];

PolygonMarker[shape_, {size:_?NumericQ | (Scaled | Offset)[_?NumericQ],angle_?NumericQ}, g_]:=PolygonMarker[shape, {size,angle}, {g}];
 (* This form allows to pass any Graphics options as an argument of PolygonMarker *)
PolygonMarker[shape_, size:_?NumericQ | (Scaled | Offset)[_?NumericQ], style_, opts:OptionsPattern[]] := Block[{gr = PolygonMarker[shape, size, style]}, Show[gr, opts] /; Head[gr] === Graphics];
PolygonMarker[shape_, {size:_?NumericQ | (Scaled | Offset)[_?NumericQ],angle_?NumericQ}, style_, opts:OptionsPattern[]] := Block[{gr = PolygonMarker[shape, {size,angle}, style]}, Show[gr, opts] /; Head[gr] === Graphics];

End[];

EndPackage[];
