Do[
 Print[ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"]];
 Print[ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"]];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 PairsIndices = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 BoundaryList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 NonBoundaryList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 bwVert = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 rholist = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 CornerVert = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 OrigLine = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 DualLine = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 FirstSecondVerts = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 AdjList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 DualList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 OrigList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 PairToCrossingVerts = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 AdjToInfty = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 AdjToAdjToInftyList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 KitePairs = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 Face2VertexList = ToExpression[a];
 a = ReadLine["C:\\Users\\14846\\Desktop\\9_250math.txt"];
 Vertex2FaceList = ToExpression[a];
 
 RhoIndex2VertexDict = <||>;
 For[i = 1, i < Length[bwVert] + 1, i++, 
  AssociateTo[RhoIndex2VertexDict, i -> bwVert[[i]]]];
 AdjDict = <| |>;
 For[i = 1, i < Length[AdjList], i += 2, 
  AssociateTo[AdjDict, AdjList[[i]] -> AdjList[[i + 1]]]];
 PairToCVDict = <||>;
 For[i = 1, i < Length[PairToCrossingVerts], i += 2, 
  AssociateTo[PairToCVDict, 
   PairToCrossingVerts[[i]] -> PairToCrossingVerts[[i + 1]]]];
 AdjToAdjToInftyDict = <||>;
 For[i = 1, i < Length[AdjToAdjToInftyList], i += 2, 
  AssociateTo[AdjToAdjToInftyDict, 
   AdjToAdjToInftyList[[i]] -> AdjToAdjToInftyList[[i + 1]]]];
 Vertex2FaceDict = <||>;
 For[i = 1, i < Length[Vertex2FaceList], i += 2, 
  AssociateTo[Vertex2FaceDict, 
   Vertex2FaceList[[i]] -> Vertex2FaceList[[i + 1]]]];
 Face2VertexDict = <||>;
 For[i = 1, i < Length[Face2VertexList], i += 2, 
  AssociateTo[Face2VertexDict, 
   Face2VertexList[[i]] -> Face2VertexList[[i + 1]]]];
 
 rho = Table[Subscript[r, i], {i, 1, Length[rholist]}];
 F[x_] := (1/2) I (PolyLog[2, -I E^x] - PolyLog[2, I E^x]);
 G[a_, b_] := F[a - b] + F[b - a] - (\[Pi]/2) (a + b);
 BSNew = Sum[
    G[Subscript[r, PairsIndices[[i]][[1]]], Subscript[r, 
     PairsIndices[[i]][[2]]]], {i, 1, Length[PairsIndices]}] + 
   Sum[Subscript[r, v]*If[MemberQ[BoundaryList, v], Pi, 2*Pi], {v, 1, 
     Length[rho]}];
 NM = NMinimize[{Re[BSNew], 
    Sum[Subscript[r, j], {j, 1, Length[rho]}] == 0}, rho, 
   AccuracyGoal -> 35, PrecisionGoal -> 35, WorkingPrecision -> 45];
 RADII = Exp[Re[rho /. (NM[[2]])]];
 RADII = RADII/
   RADII[[1]];(* normalize r1 \[Rule] 1, and drop some decimals \
(inaccurate?) *)
 radii = <||>;
 For[i = 1, i <= Length[rho], i++,
  AssociateTo[radii, Lookup[RhoIndex2VertexDict, i] -> RADII[[i]]]]; 
 
 SendToInt[\[Epsilon]_ : 0][x_] := 
  If[Abs[x - Floor[x]] < \[Epsilon] , Floor[x],
   If[Abs[x - Ceiling[x]] < \[Epsilon] , Ceiling[x], x], x];
 
 getCartesian[invCoordList_] := Module[{bz1, bz2, b, z2, z1},
   bz1 = invCoordList[[3]];
   bz2 = invCoordList[[4]];
   b = invCoordList[[2]];
   z1 = bz1/b;
   z2 = bz2/b;
   {z1, z2}];
 
 getCartesianAndRad[invCoordList_] := Module[{bz1, bz2, b, z2, z1},
   bz1 = invCoordList[[3]];
   bz2 = invCoordList[[4]];
   b = invCoordList[[2]];
   z1 = bz1/b;
   z2 = bz2/b;
   {{z1, z2}, 1/b}];
 
 getRemainingKiteVertex[vertex1_, vertex2_, vertex3_] := 
  Module[{vert1 = vertex1, vert2 = vertex2, vert3 = vertex3},
   list1 = Lookup[Vertex2FaceDict, vert1];
   list2 = Lookup[Vertex2FaceDict, vert2];
   list3 = Lookup[Vertex2FaceDict, vert3];
   commonlist = Intersection[list1, list2, list3];
   VertexSet = Lookup[Face2VertexDict, commonlist[[1]]];
   VertexSet = DeleteCases[VertexSet, vert1];
   VertexSet = DeleteCases[VertexSet, vert2];
   VertexSet = DeleteCases[VertexSet, vert3];
   VertexSet[[1]]];
 (** This function reflects a point over a line defined by two points \
**)
 reflection[PointToReflect_, Point1OfLineToReflectOver_, 
   Point2OfLineToReflectOver_] := 
  Module[{p = PointToReflect, lp1 = Point1OfLineToReflectOver, 
    lp2 = Point2OfLineToReflectOver, x1, x2, y1, y2, slope, deltax, 
    deltay}, 
   x1 = lp1[[1]];
   x2 = lp2[[1]];
   y1 = lp1[[2]];
   y2 = lp2[[2]];
   If[ (** IF the line is NOT vertical **) x1 != x2,
    If[(** AND IF the line is NOT horizontal **) y1 != y2,
     slope = (y2 - y1)/(x2 - x1);
     orthogonalSlope = -1/slope;
     normalVectorToLine = {1, orthogonalSlope};
     rt = ReflectionTransform[normalVectorToLine, lp1];
     answer = rt[p],
     (** ELSE the line is not vertical, but horizontal, 
     so compute this **)
     deltay = 2*Abs[p[[2]] - y1];
     answer = {p[[1]], p[[2]] + deltay}, (** 
     ELSE the line is vertical so compute this **)
     deltax = 2*Abs[p[[1]] - x1];
     answer = {p[[1]] + deltax, p[[2]]}
     ]
    ]; answer
   ];
 Q = ( {
    {0, 1/2, 0, 0},
    {1/2, 0, 0, 0},
    {0, 0, -1, 0},
    {0, 0, 0, -1}
   } );
 InvCoordDict = <||>;
 AssociateTo[InvCoordDict, OrigLine -> {0, 0, 0, -1}];
 AssociateTo[InvCoordDict, DualLine -> {0, 0, -1, 0}];
 FirstSecondVerts = Lookup[AdjDict, CornerVert];
 CVDone = {CornerVert};
 CoordKnown = {OrigLine, DualLine};
 
 
 (** For each vertex in FirstSecondVerts--These are the vertices \
adjacent to CornerVert, and the centers of the first circles we are \
placing **)
 For[i = 1, i < Length[FirstSecondVerts] + 1, i++,
  coor = {x, y, z, w};
  eq1 = coor.Q.coor == -1;
  eq2 = y == 1/Lookup[radii, FirstSecondVerts[[i]]];
  If[ (** If the vertex of FirstSecondVerts is an original vertex**) 
   MemberQ[OrigList, FirstSecondVerts[[i]]], (** 
   THEN these equations hold**) 
   eq3 = coor.Q.Lookup[InvCoordDict, DualLine] == 0;
   eq4 = coor.Q.Lookup[InvCoordDict, OrigLine] == 1; 
   (**We solve these equations and put the solution in the \
InvCoordDict**)
   AssociateTo[InvCoordDict, 
    FirstSecondVerts[[i]] -> {x, y, z, w} /. 
     Solve[eq1 && eq2 && eq3 && eq4, {x, y, z, w}][[1]]];
   firstwhitevert = FirstSecondVerts[[i]];
   AppendTo[CoordKnown, FirstSecondVerts[[i]]],(**
   OTHERWISE--the vertex is a dual vertex and these equations hold**)
\

   	eq3 = coor.Q.Lookup[InvCoordDict, DualLine] == 1;
   	eq4 = coor.Q.Lookup[InvCoordDict, OrigLine] == 0; 
   	(**We solve these equations and put the solution in the \
InvCoordDict**)
   	AssociateTo[InvCoordDict, 
    FirstSecondVerts[[i]] -> {x, y, z, w} /. 
     Solve[eq1 && eq2 && eq3 && eq4, {x, y, z, w}][[1]]];
   	firstblackvert = FirstSecondVerts[[i]] ;
   	AppendTo[CoordKnown, FirstSecondVerts[[i]]]
   ];
  ];
 CV = {CornerVert};
 CartesianDict = <||>;
 InCartesianList = {CornerVert};
 AssociateTo[CartesianDict, CornerVert -> {0, 0}];
 
 kitelist = 
  Lookup[PairToCVDict, {{firstwhitevert, firstblackvert}}][[1]];
 Length[kitelist];
 branchinglist = {};
 For[i = 1, i < Length[kitelist] + 1, i++, 
  AppendTo[branchinglist, kitelist[[i]]]];
 
 (** getKite() Function **)
 getKite[cv_] := Module[{CV = cv, Radius},
   knownB = 1000;
   knownW = 1000;
   (** Print["We are working with:"]; **)
   (** Print[CV]; **)
   If[(** If this crossing vertex is "done"**) 
    MemberQ[CV, CVDone], (** THEN delete it from the branching list**) 
    DeleteCases[branchinglist, CV],  (** 
    ELSE proceed with the function**)
    adlist = Lookup[AdjDict, CV];
    (** Print[adlist]; **)
    knownWList = {};
    knownBList = {};
    For[(** For each vertex in the adjacency list **) i = 1, 
     i < Length[adlist] + 1, i++, 
     (** Print[adlist[[i]]];
     Print["Adjacency List above"]; **)
     If[ (** If the vertex already has coordinates **) 
      MemberQ[CoordKnown, adlist[[i]]],
      (**THEN delete the vertex from the adjacency list and check \
whether the vertex is original or dual **)
      DeleteCases[adlist, adlist[[i]]];
      (** Print["Vertex already has coordinates"] **)
      If[(** If the vertex is original **) 
       MemberQ[OrigList, adlist[[i]]], (** 
       THEN label this vertex as knownW **)
       knownW = adlist[[i]];
       AppendTo[knownWList, knownW],
       (** Print["This is a known White vert"] **) (** 
       ELSE label this vertex as known B **)
       knownB = adlist[[i]];
       AppendTo[knownBList, knownB]
       (** Print["This a known Black vert"] **)
       ]
      ]
     ];
    (** Note that knownWList and knownBList will have maximum length \
2 **)
    If[(** Only proceed if there is a knownW and knownB **) 
     knownW != 1000 && knownB != 1000,
     (** Print[
     "We are proceeding because we have enough information"]; **)
     pair1 = {knownWList[[1]], knownBList[[1]]};
     If[(** If this is a legitimate pair **) 
      MemberQ[KitePairs, pair1],(** 
      then look up the known cross vert **) 
      knownCrossVert = 
       getRemainingKiteVertex[CV, pair1[[1]], pair1[[2]]];
      If[(** If we have visited this known cross vert, i.e. 
       it has coordinates **) MemberQ[CVDone, knownCrossVert],
       CVCartesianCoord = 
         reflection[Lookup[CartesianDict, knownCrossVert], 
          getCartesian[Lookup[InvCoordDict, pair1[[2]]]], 
          getCartesian[Lookup[InvCoordDict, pair1[[1]]]]];
       ]; 
      
      If[Length[knownBList == 2],
       pair2 = {knownWList[[1]], knownBList[[2]]};
       If[MemberQ[KitePairs, pair2], 
        knownCrossVert = 
         getRemainingKiteVertex[CV, pair2[[1]], pair2[[2]]];
        If[MemberQ[CVDone, knownCrossVert], 
         CVCartesianCoord = 
           reflection[Lookup[CartesianDict, knownCrossVert], 
            getCartesian[Lookup[InvCoordDict, pair2[[2]]]], 
            getCartesian[Lookup[InvCoordDict, pair2[[1]]]]];
         ]]]; 
      
      If[Length[knownWList] == 2,
       pair3 = {knownWList[[2]], knownBList[[1]]};
       If[MemberQ[KitePairs, pair3], 
        knownCrossVert = 
         getRemainingKiteVertex[CV, pair3[[1]], pair3[[2]]];
        If[MemberQ[CVDone, knownCrossVert], 
         CVCartesianCoord = 
           reflection[Lookup[CartesianDict, knownCrossVert], 
            getCartesian[Lookup[InvCoordDict, pair3[[2]]]], 
            getCartesian[Lookup[InvCoordDict, pair3[[1]]]]];
         ]]]; 
      
      
      If[Length[knownWList] == 2 && Length[knownBList] == 2,
       pair4 = {knownWList[[2]], knownBList[[2]]};
       If[MemberQ[KitePairs, pair4], 
        knownCrossVert = 
         getRemainingKiteVertex[CV, pair4[[1]], pair4[[2]]];
        If[MemberQ[KitePairs, pair4], 
         knownCrossVert = 
          getRemainingKiteVertex[CV, pair4[[1]], pair4[[2]]];
         If[MemberQ[CVDone, knownCrossVert], 
          CVCartesianCoord = 
            reflection[Lookup[CartesianDict, knownCrossVert], 
             getCartesian[Lookup[InvCoordDict, pair4[[2]]]], 
             getCartesian[Lookup[InvCoordDict, pair4[[1]]]]];
          
          
          ]]]];
      
      
      If[(** If CV is not the corner vert **)CV != CornerVert,
       (** Print["knownW,", knownWList];
       Print["knownB,", knownBList];
       Print["knownCrossVert", knownCrossVert]; **)
       AssociateTo[CartesianDict, CV -> CVCartesianCoord]];
      ;
      
      
      For[(** AGAIN for each vertex in the adjacency list **) i = 1, 
       i < Length[adlist] + 1, i++,
       If [(** 
        If the vertex does NOT have coordinates **) ! 
         MemberQ [CoordKnown, adlist[[i]]],
        (**THEN proceed to get coordinates **)
        (** Print["We are getting coordinates for a vertex"]; **)
        (** Print[adlist[[i]]]; **)
        coor = {x, y, z, w};
        eq1 = coor.Q.coor == -1;
        Radius = Lookup[radii, adlist[[i]]];
        eq2 = y == 1/Lookup[radii, adlist[[i]]];
        (** Print[Lookup[radii, adlist[[i]]]]; **)
        If[(** IF the vertex is original, use these equations **) 
         MemberQ[OrigList, adlist[[i]]],
         eq3 = coor.Q.Lookup[InvCoordDict, knownW] == 1;
         (** Print["InversiveCoordinates"]; **)
         (** Print[Lookup[InvCoordDict, knownW]]; **)
         eq4 = coor. Q.Lookup[InvCoordDict, knownB] == 0;
         eq5 = coor.Q.Lookup[InvCoordDict, OrigLine] >= .95,
         (** ELSE the vertex is dual and use these equations **)
         eq3 = coor.Q.Lookup[InvCoordDict, knownW] == 0;
         (** Print["InversiveCoordinates"]; **)
         (** Print[Lookup[InvCoordDict, knownB]]; **)
         eq4 = coor. Q.Lookup[InvCoordDict, knownB] == 1;
         eq5 = coor.Q.Lookup[InvCoordDict, DualLine] >= .95]; (** 
        I've put .95 instead of 1 to give the computations some \
leeway **)
        (** Solve the equations **)
        solutionList = 
         Solve[eq1 && eq2 && eq3 && eq4 && eq5, {x, y, z, w}];
        (** Print[solutionList]; **)
        If[(** IF there is one solution **) Length[solutionList] == 1, 
         (** Print["There is one solution"]; **)
         (** 
         Then assign the solution of inversive coordinates to the \
vertex **)
         AssociateTo[InvCoordDict, 
          adlist[[i]] -> {x, y, z, w} /. solutionList[[1]]];
         AppendTo[CoordKnown, adlist[[i]]],
         (** Print["There are two solutions"]; **)
         solution1 = {x, y, z, w} /. solutionList[[1]];
         center1 = getCartesian[solution1];
         (** Print["solution 1"]; **)
         (** Print[N[EuclideanDistance[center1, 
         CVCartesianCoord]]]; **)
         (** Print["radius is ", N[Radius]]; **)
         (** Print["the crossing verts is", CVCartesianCoord]; **)
         (** Print["solution 2"]; **)
         solution2 = {x, y, z, w} /. solutionList[[2]];
         center2 = getCartesian[solution2];
         (** 
         IF the center circle has radius approx equal to the distance \
from the CV to the circle center (a for respective solutions), 
         then choose that solution **)
         If[(N[EuclideanDistance[center1, CVCartesianCoord], 3]) == 
           N[Radius, 3], 
          AssociateTo[InvCoordDict, 
           adlist[[i]] -> {x, y, z, w} /. solutionList[[1]]];
          AppendTo[CoordKnown, adlist[[i]]],
          (** Else **)
          
          AssociateTo[InvCoordDict, 
           adlist[[i]] -> {x, y, z, w} /. solutionList[[2]]];
          AppendTo[CoordKnown, adlist[[i]]]
          ]];
        
        potentialcrossingverts = Lookup[AdjDict, adlist[[i]]];
        (** 
        We only add a vertex from potential crossing verts to the \
branching list if it is not CV and if it not in the CVDone list, i.e. 
        if we have not run the code on the vertex already **)
        For[j = 1, j < Length[potentialcrossingverts] + 1, j++, 
         If[potentialcrossingverts[[j]] != CV, 
          If[(! MemberQ[CVDone, potentialcrossingverts[[j]] ]), 
           AppendTo[branchinglist, potentialcrossingverts[[j]]]]]
         (** If[(!MemberQ[CVDone, potentialcrossingverts[[j]] ]), 
         AppendTo[branchinglist,potentialcrossingverts[[j]]];
         Print["We appended something to the branching list"];
         Print[potentialcrossingverts[[j]]];
         Print[branchinglist] **)
         ]
        ]
       ]; AppendTo[CVDone, CV];
      branchinglist = DeleteCases[branchinglist, CV];
      branchinglist = DeleteDuplicates[branchinglist];
      (** Print[branchinglist]; **)
      previousCV = CV
      ]
     ]
    ]
   ];
 getKite[CornerVert];
 (** Given the inversive coordinates for a circle of finite radius, \
this function returns the center of the circle **) 
 While[Length[branchinglist] != 0 (** 
  While there are vertices left in the branchinglist **),
  Scan[getKite, branchinglist] 
  ] ;
 
 Module[{eq1, eq2, eq3, eq4}, 
  For[(** For each line **) i = 1, i < Length[AdjToInfty] + 1, i++, 
   If[(** IF coordinates of the line are not knonwn, 
     prooceed **) ! MemberQ[CoordKnown, AdjToInfty[[i]]], 
     (** Print["we are proceeding"]; **)
     ln = {x, y, z, w};
     eq1 = ln.Q.ln == -1;
     eq2 = y == 0;
     crossinglist = Lookup[AdjToAdjToInftyDict, AdjToInfty[[i]]];
     whitecenterverts = {};
     blackcenterverts = {};
     For[k = 1, k < Length[crossinglist] + 1, k++, 
      centerverts = Lookup[AdjDict, crossinglist[[k]]];
      (** Print[centerverts]; **)
      For[j = 1, j < Length[centerverts] + 1, j++, 
       If[MemberQ[OrigList, centerverts[[j]]], 
        AppendTo[whitecenterverts, centerverts[[j]]];
        (** Print["adding to whitecenterverts"];,
        Print["adding to blackcenterverts"]; **)
        AppendTo[blackcenterverts, centerverts[[j]]]]
       ];
      
      ];
     whitecenterverts = DeleteDuplicates[whitecenterverts];
     blackcenterverts = DeleteDuplicates[blackcenterverts];
     If[(** IF the line is an original line **) 
      MemberQ[OrigList, AdjToInfty[[i]]],
      (** THEN proceed here**) 
      
      w1c = Lookup[InvCoordDict, whitecenterverts[[1]]];
      
      eq3 = w1c.Q.ln == 1;
      eq4 = Lookup[InvCoordDict, OrigLine].Q.ln == 1,
      
      (** ELSE the line is dual and proceed here **) 
      
      w0c = Lookup[InvCoordDict, whitecenterverts[[1]]];
      
      eq3 = w0c.Q.ln == 0;
      eq4 = Lookup[InvCoordDict, DualLine].Q.ln == 1];
     
     solutionList = Solve[eq1 && eq2 && eq3 && eq4, {x, y, z, w}];
     equationlist = {eq1, eq2, eq3, eq4};
     AssociateTo[InvCoordDict, 
      AdjToInfty[[i]] -> {x, y, z, w} /. solutionList[[1]]];
     AppendTo[CoordKnown, AdjToInfty[[i]]]
     ]; 
   ]
  ];
 InvCoordVec = {};
 For[i = 1, i < Length[CoordKnown] + 1, i++, 
  vec = Simplify[
    RootApproximant[
     N[SendToInt[1/100] /@ Lookup[InvCoordDict, CoordKnown[[i]]], 
      25], 3]];
  AssociateTo[InvCoordDict, CoordKnown[[i]] -> vec];
  AppendTo[InvCoordVec, vec]
  ];
 ICVPlace2Vert = <||>;
 AllVerts = Join[OrigList, DualList];
 For[i = 1, i < Length[CoordKnown] + 1, i++, 
  AssociateTo[ICVPlace2Vert, CoordKnown[[i]] -> i]];
 AllVerts = Join[OrigList, DualList];
 CheckCluster[List_] := Module[{adlist, adlist2, intersectionlist},
   For[i = 1, i < Length[List], i++,
    (** For every circle in the cluster **)
    (** get the vertices adjacent to it **)
    adlist = Lookup[AdjDict, List[[i]]];
    (** Now for every other vertex, also get the adjacency list **)
    For[j = 1, j < Length[AllVerts], j++,
     adlist2 = Lookup[AdjDict, AllVerts[[j]]];
     intersectionlist = Intersection[adlist, adlist2];
     Cond1 = Length[intersectionlist] != 0;
     Cond2 = AllVerts[[j]] != List[[i]];
     (** If the circle vertex shares a crossing vertex with another \
vertex **)
     (** Check if the second vertex is original **)
     If[Cond1 && Cond2,
      If[MemberQ[OrigList, AllVerts[[j]]],
       (** Check if there is a 1 in the Gram  matrix **)
       If[
        GramMatrix[[(** **) 
           Lookup[ICVPlace2Vert, AllVerts[[j]]] (** **)]][[(** **) 
          Lookup[ICVPlace2Vert, OrigList[[i]]] (** **)]] != 1,
        AppendTo[CheckingIfGramIsCorrect, 5]]];
      (** Check if the second vertex is dual **)
      If[MemberQ[DualList, AllVerts[[j]]],
       (** Check if there is a 0 in the Gram matrix **)
       If[
        GramMatrix[[(** **) 
           Lookup[ICVPlace2Vert, AllVerts[[j]]](** **)]][[(** **)
          Lookup[ICVPlace2Vert, OrigList[[i]]](** **)]] != 0,
        Print[
         GramMatrix[[(** **) 
           Lookup[ICVPlace2Vert, AllVerts[[j]]](** **)]][[(** **)
          Lookup[ICVPlace2Vert, OrigList[[i]]](** **)]]];
        AppendTo[CheckingIfGramIsCorrect, 5]]]]]
    ]];
 GramMatrix = {};
 Gram[InvCoordVec_] := Module[{},
   For[j = 1, j < Length[InvCoordVec] + 1, j++, 
    vector1 = InvCoordVec[[j]];
    gramRow = {}; 
    For[k = 1, k < Length[InvCoordVec] + 1, k++, 
     vector2 = InvCoordVec[[k]]; 
     AppendTo[gramRow, Simplify[vector1.Q.vector2]]];
    AppendTo[GramMatrix, gramRow]];
   GramMatrix];
 DrawCircle[{bhat_, b_, bx_, by_}] := 
  If[b != 0, Circle[{bx, by}/b, 1/b], 
   InfiniteLine[{bx, by} bhat/2, {-by, bx}]];
 Render[DualList_, OrigList_] := Graphics[{Blue, Thick,
    Table[
     DrawCircle[Lookup[InvCoordDict, OrigList[[i]]]], {i, 1, 
      Length[OrigList]}],
    Red, Dashed,
    Table[
     DrawCircle[Lookup[InvCoordDict, DualList[[i]]]], {i, 1, 
      Length[DualList]}]
    }];
 (** Render[DualList, OrigList] **)
 cluster = {};
 cocluster = {};
 For[i = 1, i < Length[OrigList] + 1, i++,
  AppendTo[cluster, Lookup[InvCoordDict, OrigList[[i]]]]
  ];
 For[i = 1, i < Length[DualList] + 1, i++,
  AppendTo[cocluster, Lookup[InvCoordDict, DualList[[i]]]]
  ];
 Vout = MatrixForm[InvCoordVec];
 n = Length[OrigList];
 R[vhat_] := IdentityMatrix[4] + 2*Q.Transpose[{vhat}].{vhat};
 BMat[n_] := Table[b[i, j], {i, 1, n}, {j, 1, n}];
 
 Bend3[V_, n_, cocluster_, cluster_, RING_ : Integers][i_] :=(* 
  V is a supercluster matrix. iso is the indices of the cluster. 
  1\[LessEqual]i\[LessEqual]|cocluster| *)
  Module[{sol, assignments, rest, j, l, others1, others2, 
    co = cocluster, cl = cluster},
   sol = Solve[BMat[n].cl == cl.R[co[[i]]], RING]; (** Solve BV=
   VR **)
   (** Print["sol", sol]; **)
   If[sol == {}, Null,
    sol = sol[[1]];
    assignments = Association[sol]; (** 
    make a dictionary of the variables, i.e. elements of BMat, 
    and the solutions **) 
    
    (** Print["assignments", assignments]; **)
    rest = Complement[Flatten[BMat[n]], Keys[assignments]];  (** 
    rest is a list of variables that do not have a solution **)
    (** Print["rest", rest]; **)
    l = Length[rest];
    (**
    If[l==0,MatrixForm[Simplify[BMat[n]/.assignments]],
    others1=Solve[rest==Table[0,{j,1,l}]][[1]];
    others2=<||>;
    For[j=1,j\[LessEqual]l,j++,others2[rest[[j]]]=Subscript[a, j]];
    others2=Normal[others2];
    MatrixForm[Simplify[(BMat[n]/.(assignments/.others1))/.others1]]
    ]
    **)
    If[l == 0, (** IF all variables have solutions **)
     (** THEN BMat[n] /. 
     assignments puts all of the solutions into the matrix B **)
     (** /. replaces any free variable c with 0?? **)
     (** output the matrix form **)
     MatrixForm[
      Simplify[
       BMat[n] /. assignments /. {C[1] -> 0, C[2] -> 0, C[3] -> 0, 
         C[4] -> 0, C[5] -> 0, C[6] -> 0, C[7] -> 0, C[8] -> 0, 
         C[9] -> 0, C[10] -> 0, C[11] -> 0, C[12] -> 0, C[13] -> 0, 
         C[14] -> 0, C[15] -> 0, C[16] -> 0, C[17] -> 0, C[18] -> 0, 
         C[19] -> 0, C[20] -> 0, C[21] -> 0}]], (** 
     ELSE we make rest a table of 0's of length of rest?? **)
     
     others1 = Solve[rest == Table[0, {j, 1, l}]][[1]];
     
     others2 = <||>;
     (** for j from 1 to length of rest, 
     make a key in others2 equal to the index of each element of rest \
and make the value a subscripted to that index **)
     For[j = 1, j <= l, j++, others2[rest[[j]]] = Subscript[a, j]];
     (** Normal[] converts the association to a list of rules **)
     others2 = Normal[others2];
     MatrixForm[
      Simplify[(BMat[n] /. (assignments /. others1)) /. others1]]]
    ]];
 
 Bend4[V_, n_][i_] := 
  With[{B = Bend3[V, n, cocluster, cluster, Integers][i]},
   If[Length[B] > 0, (** 
    IF the system of equations has integer solutions **) 
    {B, "integral"}, (** 
    THEN return the bend matrix and "integral" **)
    {Bend3[V, n, cocluster, cluster, Reals][i], "nonintegral"} (** 
    Otherwise, find the real solutions and return "nonintegral" **)
    ]
   ];
 bends = {};
 CheckCluster[List_] := Module[{adlist, adlist2, intersectionlist},
   For[i = 1, i < Length[List], i++,
    (** For every circle in the cluster **)
    (** get the vertices adjacent to it **)
    adlist = Lookup[AdjDict, List[[i]]];
    (** Now for every other vertex, also get the adjacency list **)
    For[j = 1, j < Length[AllVerts], j++,
     adlist2 = Lookup[AdjDict, AllVerts[[j]]];
     intersectionlist = Intersection[adlist, adlist2];
     Cond1 = Length[intersectionlist] != 0;
     Cond2 = AllVerts[[j]] != List[[i]];
     (** If the circle vertex shares a crossing vertex with another \
vertex **)
     (** Check if the second vertex is original **)
     If[Cond1 && Cond2,
      If[MemberQ[OrigList, AllVerts[[j]]],
       (** Check if there is a 1 in the Gram  matrix **)
       If[
        GramMatrix[[(** **) 
           Lookup[ICVPlace2Vert, AllVerts[[j]]] (** **)]][[(** **) 
          Lookup[ICVPlace2Vert, OrigList[[i]]] (** **)]] != 1,
        AppendTo[CheckingIfGramIsCorrect, 5]]];
      (** Check if the second vertex is dual **)
      If[MemberQ[DualList, AllVerts[[j]]],
       (** Check if there is a 0 in the Gram matrix **)
       If[
        GramMatrix[[(** **) 
           Lookup[ICVPlace2Vert, AllVerts[[j]]](** **)]][[(** **)
          Lookup[ICVPlace2Vert, OrigList[[i]]](** **)]] != 0,
        AppendTo[CheckingIfGramIsCorrect, 5]]]]]
    ]];
 CheckGram[OrigList_, DualList_, Matrix_] := Module[{},
   CheckingIfGramIsCorrect = {};
   CheckCluster[OrigList];
   For[i = 1, i < Length[InvCoordVec], i++, 
    num = FullSimplify[GramMatrix[[i]][[i]]]; 
    If[num != -1, AppendTo[CheckingIfGramIsCorrect, 7]]];
   Length[CheckingIfGramIsCorrect] == 0];
 
 Everything[InvCoordVec_, DualList_, OrigList_] := Module[{},
   Column[{GramMatrix = Gram[InvCoordVec];
     a = CheckGram[OrigList, DualList, GramMatrix];
     Print["Gram matrix is verified, ", a]
     }]
   ];
 Everything[InvCoordVec, DualList, OrigList];
 If[a, Print[MatrixForm[GramMatrix]]; 
  bends = Bend4[Vout, Length[OrigList]] /@ Range@Length[DualList];
  Print[bends]; Render[DualList, OrigList]], 250]
