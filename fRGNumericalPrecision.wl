(* ::Package:: *)

(* ::Title:: *)
(*Numerical Accuracy of the Derivative-Expansion-Based Functional Renormalization Group*)


(* ::Text:: *)
(*	Program analyzing the precision of a numerical implementation *)
(*	of the derivative-expansion-based functional renormalization group *)
(*	in the three-dimensional O(N) models.*)
(*	This program serves as an appendix to the article under the same name. *)
(*		Article status: Preprint*)
(*		ArXiv Link: https://arxiv.org/abs/2404.18707*)
(*		Article DOI: Unassigned*)
(**)
(*	Copyright (C) 2024 Andrzej Chlebicki*)
(**)
(*	This program is free software: you can redistribute it and/or modify*)
(*	it under the terms of the GNU General Public License as published by*)
(*	the Free Software Foundation, either version 3 of the License, or*)
(*	(at your option) any later version.*)
(**)
(*	This program is distributed in the hope that it will be useful,*)
(*	but WITHOUT ANY WARRANTY; without even the implied warranty of*)
(*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*)
(*	GNU General Public License for more details.*)
(**)
(*	You should have received a copy of the GNU General Public License*)
(*	along with this program.  If not, see https://www.gnu.org/licenses/.*)


(* ::Section::Closed:: *)
(*Initialization*)


$HistoryLength = 10;


arguments = $CommandLine;
If[MemberQ[$CommandLine, "-script"], notebookMode = False, notebookMode = True];
Protect[notebookMode];

If[notebookMode, 
	SetDirectory[NotebookDirectory[]];, 
	
	notebookPath = arguments[[3]];
	slashes = StringPosition[notebookPath, "/"];
	If[Length[slashes]!=0, 
		projectPath = StringTake[notebookPath, slashes[[-1, 1]]];
		SetDirectory[projectPath];
	];
];
debug = False; (* In debug mode many logs about the calculations are printed *)

articleFiguresDirectory = Directory[] <> "/Article figures/";


Print["
Numerical Accuracy of the derivative-expansion-based functional renormalization group  Copyright (C) 2024 Andrzej Chlebicki
	Program investigating the precision of a numerical implementation 
	of the Derivative-Expansion-Based Functional Renormalization Group 
	in the three-dimensional O(N) models.
    
	This program comes with ABSOLUTELY NO WARRANTY.
	This is free software, and you are welcome to redistribute it under certain conditions.
"];


(* ::Section::Closed:: *)
(*Numerical Library*)


(* ::Subsubsection::Closed:: *)
(*PrintLog*)


(* ::Input:: *)
(*PrintLog[message_]:= If[!TrueQ[debug == False], Print[message]]*)


(* ::Subsection::Closed:: *)
(*Transforming associations*)


(* ::Subsubsection::Closed:: *)
(*TransposeAssociation*)


(* ::Text:: *)
(*Exchanges primary and secondary keys in a nested dictionary*)


TransposeAssociation[association_] := 
	Block[{i, j, transposed, primaryKeys, primaryKey, secondaryKeys, secondaryKey}, 
	transposed = Association[];
	secondaryKeys = Keys[association];
	primaryKeys = Sort[DeleteDuplicates[Flatten[Values[Map[Keys, association]], 1]], #1>#2&];
	For[i=1, i<=Length[primaryKeys], i++, 
		primaryKey = primaryKeys[[i]];
		transposed[primaryKey] = Association[];
		For[j=1, j<=Length[secondaryKeys], j++, 
			secondaryKey = secondaryKeys[[j]];
			If[MemberQ[Keys[association[secondaryKey]], primaryKey], 
				transposed[primaryKey][secondaryKey] = association[secondaryKey][primaryKey]];
			];
		];
	Return[transposed];
];


(* ::Subsubsection::Closed:: *)
(*OneLevelNest*)


(* ::Text:: *)
(*Transforms a dictionary with iterable keys into a nested dictionary. One level nest - extracts last value from a list-like key and organizes it into a nested dictionary*)


OneLevelNest[association_] := 
	Block[{i, newEntry, nested, key, keyLastElement, keyRemainder}, 
	nested = Association[];
	If[Length[Keys[association]] == 0 || Length[Keys[association][[1]]] == 0, 
		Return[association]]; 
	For[i=1, i<=Length[association], i++, 
		key = Keys[association][[i]];
		keyLastElement = key[[-1]];
		keyRemainder = key[[1;;-2]];
		If[Length[keyRemainder] == 1, keyRemainder = keyRemainder[[1]]];
		newEntry = Association[keyLastElement -> association[[i]]];
		If[MemberQ[Keys[nested], keyRemainder], 
			nested[keyRemainder] = KeySort[Join[nested[keyRemainder], newEntry], #1>#2&], 
			nested[keyRemainder] = newEntry
		];
	];

	Return[nested];
];


(* ::Subsubsection::Closed:: *)
(*NestAssociation*)


(* ::Text:: *)
(*Transforms a dictionary with iterable keys into a nested dictionary*)


NestAssociation[association_] := 
	Block[{newEntry, nested}, 

	nested = association;
	While[Length[Keys[nested]] > 0 && Keys[nested][[1]][[0]] == List, 
		nested = KeySort[OneLevelNest[nested], #1>#2&];
	];

	Return[nested];
];


(* ::Subsubsection::Closed:: *)
(*UnnestAssociation*)


(* ::Text:: *)
(*Combines keys of a nested dictionary into first-level dictionary*)


UnnestAssociation[association_] := 
	Block[{i, j, keysLevels, level, newlevel, newkey, primaryKey, secondaryKey}, 
	keysLevels = {};

	level = association;
	While[AssociationQ[level[[1]]], 
		newlevel = Association[];
		For[i=1, i<=Length[level], i++, 
			primaryKey = Keys[level][[i]];
			For[j=1, j<=Length[level[[i]]], j++, 
				secondaryKey = Keys[level[primaryKey]][[j]];
				newkey = Flatten[{primaryKey, secondaryKey}];
				newlevel[newkey] = level[primaryKey][secondaryKey];
			];
		];
		level = newlevel;
	];
	Return[level];
];


(* ::Subsubsection::Closed:: *)
(*Normal2*)


(* ::Text:: *)
(*Transform an association into a list of pairs {Key, Value}*)


Normal2[assoc_] := Transpose[{Keys[assoc], Values[assoc]}]


(* ::Subsection::Closed:: *)
(*Finite difference approximation for derivatives*)


derivativePointsCount = 5; (* Number of points used in derivative approximations on the \[Rho] grid, can be changed later *)
maxDerivativeOrder = 2; (* Maximal order of derivative appearing in the \[Beta] functions *)
derivativePointsOneSide := (derivativePointsCount-1)/2; 
Protect[derivativePointsOneSide];
maximalNumericalPrecision = 60;


(* ::Subsubsection::Closed:: *)
(*Finite difference coefficients*)


GetDerCoefs[n_, d_, shift_]:= Block[{rep, stencil, mat, vec, MyCForm, pos, declaration, half}, 
	rep = s[j_] -> j-shift;
	stencil = Table[s[j], {j, 0, n-1}] /. rep;
	mat = Table[s[j]^i, {i, 0, n-1}, {j, 0, n-1}] /. rep;
	vec = d! Table[KroneckerDelta[d, i], {i, 0, n-1}];
	Return[Inverse[mat] . vec];
]


EvaluationPoints[f_, x0_] := Block[{start, end}, 
	start = x0 - derivativePointsOneSide;
	end = x0 + derivativePointsOneSide;
	If[start < 0, start = 0; end = derivativePointsCount-1];
	If[end > gridSize, start = gridSize - derivativePointsCount+1; end = gridSize];
	Return[Thread[f[ Range[start, end]]]];
]


(* ::Subsubsection::Closed:: *)
(*Function for replacing derivative with finite differences*)


NumericDerivatives::noGrid := "gridSize not specified as a big enough integer"
NumericDerivatives::derivativeError := "Finite difference approximation parameters are not set up properly"
NumericDerivatives[expression_] := Block[{lowBorderReps, highBorderReps, centralReps, fullReps}, 
	If[!IntegerQ[derivativePointsCount] || !IntegerQ[maxDerivativeOrder] || derivativePointsCount < maxDerivativeOrder+1, 
		Message[NumericDerivatives::derivativeError];
		Return[Null];
	];

	If[!IntegerQ[gridSize] || gridSize <= derivativePointsCount, 
		Message[NumericDerivatives::noGrid];
		Return[Null];
	];

	lowBorderReps =
		Outer[Derivative[#1][ff_][#2] /; !TrueQ[ff == r] -> 
			GetDerCoefs[derivativePointsCount, #1, #2] . EvaluationPoints[ff, #2]&, 
			Range[maxDerivativeOrder], Range[0, (derivativePointsCount-1)/2-1]];

	highBorderReps = 
		Outer[Derivative[#1][ff_][gridSize+#2] /; !TrueQ[ff == r] -> 
			GetDerCoefs[derivativePointsCount, #1, derivativePointsCount-1+#2] . EvaluationPoints[ff, gridSize+#2]&, 
			Range[maxDerivativeOrder], Range[-(derivativePointsCount-1)/2+1, 0]];

	centralReps = Map[Derivative[#][ff_][x_]/; !TrueQ[ff == r] -> 
		GetDerCoefs[derivativePointsCount, #1, derivativePointsOneSide] . EvaluationPoints[ff, x]&, 
		Range[maxDerivativeOrder]];

	fullReps = Flatten[{lowBorderReps, highBorderReps, centralReps}];

	Return[expression /. fullReps];
]


(* ::Subsubsection::Closed:: *)
(*3-point first-derivative approximation on unevenly spaced grid*)


FirstDerivativeUnevenGrid[pointData_] := 
	Block[{normalData, xs, values, spacings, derivativeCoefficients, d1, d2, derivative, a0, a1, a2}, 
	
	If[TrueQ[pointData[[0]] == Association], 
		normalData = Normal2[pointData], 
		normalData = pointData;
	];
	
	xs = normalData[[All, 1]];
	spacings = xs[[2;;-1]] - xs [[1;;-2]];
	d1 = spacings[[1;;-2]];
	d2 = spacings[[2;;-1]];

	derivativeCoefficients = Solve[Table[ 
		KroneckerDelta[i, 1] f'[x] == SeriesCoefficient[a1 f[x - \[Delta]1 h] + a0 f[x] + a2 f[x + \[Delta]2 h], {h, 0, i}], 
		{i, 0, 2}], {a0, a1, a2}][[1]];

	values = normalData[[All, 2]];

	{a0, a1, a2} = {a0, a1, a2} /. derivativeCoefficients /. {\[Delta]1 -> d1, \[Delta]2 -> d2};

	derivative = a1 values[[1;;-3]] + a0 values[[2;;-2]] + a2 values[[3;;-1]] ;

	Return[Association[Map[xs[[#+1]] -> derivative[[#]]&, Range[Length[derivative]]]]];
]


(* ::Subsubsection::Closed:: *)
(*Numerical Central Derivative*)


NumericCentralDerivative[baseFunction_, parameter_, oneSideDerivativePoints_, epsilon_] := Block[{derivativePoints, coefficients, evaluations}, 
	If[FreeQ[baseFunction, parameter], Return[0]];
	derivativePoints = 2oneSideDerivativePoints + 1;
	coefficients = GetDerCoefs[derivativePoints, 1, oneSideDerivativePoints][[oneSideDerivativePoints+2;;-1]];

	evaluations = Map[ ((baseFunction /. parameter -> (parameter + # epsilon)) - 
					(baseFunction /. parameter -> (parameter - # epsilon)))&, 
		Range[1, oneSideDerivativePoints]];

	If[ListQ[baseFunction], 
		Return[Transpose[evaluations] . coefficients / epsilon], 
		Return[evaluations . coefficients / epsilon];
	];
]


(* ::Subsubsection::Closed:: *)
(*Numerical Gradient*)


NumericGradient[baseFunction_List, parameters_, oneSideDerivativePoints_, epsilon_] := Block[{}, 
	Return[Map[NumericGradient[#, parameters, oneSideDerivativePoints, epsilon]&, baseFunction]];
]


NumericGradient[baseFunction_, parameters_, oneSideDerivativePoints_, epsilon_] := Block[{}, 
	Return[Map[NumericCentralDerivative[baseFunction, #, oneSideDerivativePoints, epsilon]&, parameters]];
]


(* ::Subsection::Closed:: *)
(*Gauss-Legendre (GL) integral*)


(* ::Text:: *)
(*Gauss-Legendre integral is performed with the function GaussLegendreIntegrate over the variable 'y'.*)
(*The integral is performed on on interval y \[Element] [GLIntegralLowerBound, GLIntegralUpperBound] with GLIntegralPointCount integrand evaluations.*)
(*The GL-integral parameters mentioned above are 'protected' to secure the integrity of the program. They can only be altered via the SetGLIntegralParameters function.*)


GLIntegralReplacementExpressions = {r[y^2], r'[y^2], r''[y^2], y^2, y};
regulator = {};


SetGLIntegralParameters[parameterSet_] := SetGLIntegralParameters[parameterSet[[1]], parameterSet[[2]], parameterSet[[3]]]

SetGLIntegralParameters[lowerBound_, upperBound_, pointCount_] := Block[{}, 
	Unprotect[GLIntegralLowerBound, GLIntegralUpperBound, GLIntegralPointCount];
	GLIntegralLowerBound = lowerBound;
	GLIntegralUpperBound = upperBound;
	GLIntegralPointCount = pointCount;

	Protect[GLIntegralLowerBound, GLIntegralUpperBound, GLIntegralPointCount];

	RecalculateGLIntegralParameters[];
]


RecalculateGLIntegralParameters[] := Block[{legendrePolynomial, legendrePolynomialDerivative, legendreRoots, legendreWeights}, 
	Unprotect[GLIntegralIntervalHalfSpan, GLIntegralIntervalMidpoint, GLIntegralEvaluationPoints, 
		GLIntegralReplacementRules, GLIntegralWeights];
	
	Block[{$MaxPrecision = Infinity}, 
		legendrePolynomial = LegendreP[GLIntegralPointCount, x];
		legendrePolynomialDerivative = D[legendrePolynomial, x];
		legendreRoots = NSolve[legendrePolynomial == 0, x, WorkingPrecision -> maximalNumericalPrecision];
		legendreWeights = Map[2/((1-x^2)legendrePolynomialDerivative^2) /. # &, legendreRoots];
	];

	GLIntegralIntervalHalfSpan = (GLIntegralUpperBound-GLIntegralLowerBound) / 2;
	GLIntegralIntervalMidpoint = (GLIntegralUpperBound+GLIntegralLowerBound) / 2;
	GLIntegralEvaluationPoints = GLIntegralIntervalMidpoint + GLIntegralIntervalHalfSpan x /. legendreRoots;
	GLIntegralReplacementRules = Outer[#2 -> (Limit[#2 /. regulator, y -> #1])&, GLIntegralEvaluationPoints, GLIntegralReplacementExpressions];
	GLIntegralWeights = GLIntegralIntervalHalfSpan * legendreWeights;

	Protect[GLIntegralIntervalHalfSpan, GLIntegralIntervalMidpoint, GLIntegralEvaluationPoints, 
		GLIntegralReplacementRules, GLIntegralWeights];
	PrintLog["Gauss-Legendre integral parameters recalculated"];
]


GaussLegendreIntegrate::invalidBounds := "The loop integral calculated with the Litim regulator has bounds `1` and `2`"
GaussLegendreIntegrate[integrand_] := Block[{yvalues, fvalues}, 
	If[TrueQ[regulatorLabel == "Litim"] && (GLIntegralLowerBound != 0 || GLIntegralUpperBound != 1), 
		Message[GaussLegendreIntegrate::invalidBounds, GLIntegralLowerBound, GLIntegralUpperBound]
	];
	If[Length[GLIntegralEvaluationPoints] != GLIntegralPointCount || 
		Length[GLIntegralWeights] != GLIntegralPointCount, 
		
		RecalculateGLIntegralParameters[];
	];
	
	Return[Total[MapThread[(#1*integrand /. #2&), {GLIntegralWeights, GLIntegralReplacementRules}]]];

];


(* ::Text:: *)
(*Standard sets of parameters for Gauss-Legendre integrals, precision is tested at the end of the notebook.*)


GLIntegralLowPrecisionParameters = {0, 5, 20}; (* Low precision, maximal error: 10^-5 *)
GLIntegralStandardParameters = {0, 5, 35}; (* Standard precision, maximal error: 10^-10 *)
GLIntegralReferenceParameters = {0, 5.5, 50}; (* High precision, maximal error: 10^-13 *)
SetGLIntegralParameters[GLIntegralStandardParameters];


(* ::Subsection::Closed:: *)
(*Simpson integral (3/8 rule)*)


regulator = {};
simsponIntegralReplacementExpressions = {r[y^2], r'[y^2], r''[y^2], y^2, y};


RecalculateSimpsonIntegralParameters[] := Block[{}, 
	Unprotect[simpsonIntegralSmallParameter, simpsonIntegralEvaluationPoints, simpsonIntegralReplacementRules];
	simpsonIntegralSmallParameter = (simpsonIntegralUpperBound-simpsonIntegralLowerBound)/simpsonIntegralPointCount;
	simpsonIntegralEvaluationPoints = N[Range[0, simpsonIntegralPointCount] * simpsonIntegralSmallParameter, $MinPrecision];

	simpsonIntegralReplacementRules = Outer[#2 -> (Limit[#2 /. regulator, y -> #1])&, 
		simpsonIntegralEvaluationPoints, simsponIntegralReplacementExpressions];

	Protect[simpsonIntegralSmallParameter, simpsonIntegralEvaluationPoints, simpsonIntegralReplacementRules];
	PrintLog["Simpson integral parameters recalculated"];
]


SetSimpsonIntegralParameters[lowerBound_:simpsonIntegralLowerBound, upperBound_:simpsonIntegralUpperBound, pointCount_:simpsonIntegralPointCount] := Block[{}, 
	Unprotect[simpsonIntegralLowerBound, simpsonIntegralUpperBound, simpsonIntegralPointCount, 
		simpsonIntegralSmallParameter, simpsonIntegralEvaluationPoints];
	
	simpsonIntegralLowerBound = lowerBound;
	simpsonIntegralUpperBound = upperBound;
	simpsonIntegralPointCount = Ceiling[pointCount/3]*3;
	RecalculateSimpsonIntegralParameters[];
	
	Protect[simpsonIntegralLowerBound, simpsonIntegralUpperBound, simpsonIntegralPointCount, 
		simpsonIntegralSmallParameter, simpsonIntegralEvaluationPoints];
]


SetSimpsonIntegralParameters[0, 8, 60]


SimpsonIntegrate::invalidBounds := "The loop integral calculated with the Litim regulator has bounds `1` and `2`"
SimpsonIntegrate[integrand_] := Block[{ coefficients, yvalues, fvalues}, 
	If[TrueQ[regulatorLabel == "Litim" || regulatorLabel == "Litim2"] && 
		(simpsonIntegralLowerBound != 0 || simpsonIntegralUpperBound != 1), 
		Message[SimpsonIntegrate::invalidBounds, simpsonIntegralLowerBound, simpsonIntegralUpperBound]
	];
	If[Length[simpsonIntegralEvaluationPoints] != 1+simpsonIntegralPointCount || 
		simpsonIntegralLowerBound != simpsonIntegralEvaluationPoints[[1]] || 
		simpsonIntegralUpperBound != simpsonIntegralEvaluationPoints[[-1]] ||
		simpsonIntegralEvaluationPoints[[1]] != simpsonIntegralEvaluationPoints[[2]] - simpsonIntegralSmallParameter, 
		
		RecalculateSimpsonIntegralParameters[];
	];

	coefficients = Map[2 + Boole[Mod[ #, 3]!=0]&, Range[0, simpsonIntegralPointCount]];
	coefficients[[1]] = 1;
	coefficients[[-1]] = 1;
	fvalues= MapThread[(#1*integrand /. #2&), {coefficients, simpsonIntegralReplacementRules}];

	Return[simpsonIntegralSmallParameter 3/8 * Total[fvalues]];
];


(* ::Subsection::Closed:: *)
(*Working precision manipulation*)


(* ::Subsubsection::Closed:: *)
(*Standard C++ floating point types*)


floatPrecision = 24 Log10[2];
doublePrecision = MachinePrecision; (* 53 Log10[2]; *)
longDoublePrecision = 64 Log10[2];

floatLabel = "Float (32 bit)";
doubleLabel = "Double (64 bit)";
longDoubleLabel = "Long double (80 bit)";

floatingPointTypes = {{floatLabel, floatPrecision}, {doubleLabel, doublePrecision}, {longDoubleLabel, longDoublePrecision}};


(* ::Subsubsection::Closed:: *)
(*Recasting floating point numbers to high precision*)


(* ::Text:: *)
(*Takes a finite-precision number and recasts them as an exact number with the number of decimal digits equal to `precision`*)


RecastFloatingPointReal[number_, precision_:maximalNumericalPrecision]:= 
	Block[{digits, extendedNumber, currentPrecision}, 
	
	currentPrecision = Min[maximalNumericalPrecision, Precision[number]];
	extendedNumber = SetPrecision[number, maximalNumericalPrecision];
	digits = RealDigits[extendedNumber, 10, maximalNumericalPrecision];
	digits[[1]] = digits[[1, 1;;Ceiling[currentPrecision]]];
	Return[N[FromDigits[digits, 10], precision]];
]


RecastFloatingPoint[number_, precision_:maximalNumericalPrecision]:= Block[{digits, extendedNumber}, 
	Return[Sign[Re[number]]RecastFloatingPointReal[Abs[Re[number]], precision] 
		+ I Sign[Im[number]] RecastFloatingPointReal[Abs[Im[number]], precision]];
]


SetAttributes[RecastFloatingPointReal, Listable]
SetAttributes[RecastFloatingPoint, Listable]


(* ::Subsubsection::Closed:: *)
(*Setting fixed floating point precision*)


(* ::Text:: *)
(*Switches from abitrary-precision calculations to fixed-precision calculations*)


SetFixedFloatPrecision::PrecisionLimit := "Precision of `1` digits exceeds the maximal value of `2` digits.";
SetFixedFloatPrecision[nPrecision_]:= Block[{}, 
	If[nPrecision > maximalNumericalPrecision, 
		Message[SetFixedFloatPrecision::PrecisionLimit, nPrecision, maximalNumericalPrecision]
	];

	$MaxExtraPrecision =0;
	If[nPrecision <= $MinPrecision, 
		$MinPrecision = nPrecision;
		$MaxPrecision = nPrecision, 

		$MaxPrecision = nPrecision;
		$MinPrecision = nPrecision;
	];

	$MaxExtraPrecision = 0;
	(*RecalculateSimpsonIntegralParameters[];*);
	RecalculateGLIntegralParameters[];
]


(* ::Subsubsection::Closed:: *)
(*Restoring arbitrary floating point precision*)


(* ::Text:: *)
(*Switches back to abitrary-precision calculations from fixed-precision calculations*)


RestoreArbitraryPrecision[] := Block[{}, 
	$MaxExtraPrecision = 50;
	$MaxPrecision = Infinity;
	$MinPrecision = 0;
]


(* ::Section::Closed:: *)
(*Visualisation Library*)


(* ::Subsection::Closed:: *)
(*Data visualisation styles*)


DefaultPlotColors = ColorData[97, "ColorList"];


(* ::Subsubsection::Closed:: *)
(*Lines*)


DefaultThicknessValue = 4;
DefaultThickness = AbsoluteThickness[DefaultThicknessValue];
DefaultDashing = {Null, 
AbsoluteDashing[{4 DefaultThicknessValue, 1.5 DefaultThicknessValue}], (* Dashed *)
AbsoluteDashing[{1.5 DefaultThicknessValue, 1.5 DefaultThicknessValue}], (* Dotted *)
AbsoluteDashing[{4 DefaultThicknessValue, 1.5 DefaultThicknessValue, 1.5 DefaultThicknessValue, 1.5 DefaultThicknessValue}] (* Dot-dashed *)} ;


LineDefaultStyle = Table[Directive[DefaultPlotColors[[i]], DefaultThickness, DefaultDashing[[i]]], {i, 1, Min[Length[DefaultDashing], Length[DefaultPlotColors]]}];
LineFunctionalFitStyle = Map[Directive[Black, DefaultThickness, #]&, DefaultDashing];


(* ::Subsubsection::Closed:: *)
(*Markers*)


DefaultMarkerSize = 26;
DefaultLegendMarkerSize = 32;


FilledMarkers= {"\[FilledCircle]", "\[FilledSquare]", "\[FilledUpTriangle]", "\[FilledDownTriangle]"};
OpenMarkers= {"\[EmptyCircle]", "\[EmptySquare]", "\[EmptyUpTriangle]", "\[EmptyDownTriangle]"};
StarMarkers= {"\[FivePointedStar]", "\[SixPointedStar]", "\[FilledDiamond]"};


ScatterDefaultStyle = Table[Directive[DefaultPlotColors[[i]], DefaultThickness, DefaultDashing[[i]]], 
	{i, 1, Min[Length[DefaultDashing], Length[DefaultPlotColors]]}];
ScatterDefaultMarkers := Map[{#, DefaultMarkerSize}&, FilledMarkers];
ScatterDefaultLegendMarkers := Map[{#, DefaultLegendMarkerSize}&, FilledMarkers];


(* ::Subsection::Closed:: *)
(*Labels*)


DefaultLabel = FontSize -> 32;


(* ::Subsection::Closed:: *)
(*Legends*)


CommonLegendStyles = {};
ScatterLegendStyle =Join[{LegendMarkers -> ScatterDefaultLegendMarkers, LegendMarkerSize -> {24, 24}}, CommonLegendStyles];
LineLegendStyle = Join[{LegendMarkerSize -> {48, 24}}, CommonLegendStyles];
LineScatterLegendStyle = Join[{LegendMarkers -> ScatterDefaultLegendMarkers, LegendMarkerSize -> {48, 24}}, CommonLegendStyles];


(* ::Subsection::Closed:: *)
(*Default plot settings*)


CommonSettings = {LabelStyle -> DefaultLabel, 
	GridLinesStyle -> LightGray, ImageSize -> Large};
LinePlotSettings = Join[{PlotStyle -> LineDefaultStyle}, CommonSettings];
FitPlotSettings = Join[{PlotStyle -> LineFunctionalFitStyle}, CommonSettings];
ScatterPlotSettings = Join[{PlotStyle -> ScatterDefaultStyle, PlotMarkers -> ScatterDefaultMarkers}, CommonSettings];
LineScatterPlotSettings = Join[{PlotStyle -> ScatterDefaultStyle, 
	PlotMarkers -> ScatterDefaultMarkers, Joined -> True}, CommonSettings];


(* ::Subsection::Closed:: *)
(*Plot functions*)


MakeLinePlot[functions_, argument_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	Plot[functions, argument, Evaluate[LinePlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[LineLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeLineLogLogPlot[functions_, argument_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	LogPlot[functions, argument, Evaluate[LinePlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[LineLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeLineLogLogPlot[functions_, argument_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	LogPlot[functions, argument, Evaluate[LinePlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[LineLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeFitPlot[functions_, argument_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	Plot[functions, argument, Evaluate[FitPlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[LineLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeFitLogPlot[functions_, argument_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	LogPlot[functions, argument, Evaluate[FitPlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[LineLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeFitLogLogPlot[functions_, argument_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	LogLogPlot[functions, argument, Evaluate[FitPlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[LineLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeScatterPlot[data_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	ListPlot[data, Evaluate[ScatterPlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[ScatterLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeScatterLogPlot[data_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	ListLogPlot[data, Evaluate[ScatterPlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[ScatterLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


MakeScatterLogLogPlot[data_, legend_, legendPosition_:Bottom, extra_:{}, legendExtra_:{}] := 
	ListLogLogPlot[data, Evaluate[ScatterPlotSettings], 
	PlotLegends -> Placed[LineLegend[legend, Evaluate[ScatterLegendStyle], Evaluate[legendExtra]], legendPosition], 
	Evaluate[extra]]


(* ::Subsection::Closed:: *)
(*Plot testing*)


(* ::Input:: *)
(*testLines = {x, x^2, x^3, x^4};*)
(*testData = Map[Transpose[{x, #} /. {x -> Range[0, 20]/10}]&, testLines];*)
(*testLegend = {x, Superscript[x, 2], Superscript[x, 3], Superscript[x, 4]};*)


(* ::Input:: *)
(*ListPlot[testData, Evaluate[ScatterPlotSettings], PlotLegends -> PointLegend[testLegend, Evaluate[ScatterLegendStyle]] ]*)


(* ::Input:: *)
(*ListPlot[testData, Evaluate[LineScatterPlotSettings], PlotLegends -> PointLegend[testLegend, Evaluate[LineScatterLegendStyle]] ]*)


(* ::Input:: *)
(*Plot[testLines, {x, 0, 2}, Evaluate[LinePlotSettings], PlotLegends -> Placed[LineLegend[testLegend, Evaluate[LineLegendStyle]], Below]]*)


(* ::Subsection::Closed:: *)
(*Common complex symbols*)


\[Rho]0 = ToString[Subscript[OverTilde[\[Rho]], 0], StandardForm];
N\[Rho] = ToString[Subscript[N, \[Rho]], StandardForm];
h\[Rho] = ToString[Subscript[h, \[Rho]], StandardForm];


(* ::Section::Closed:: *)
(*NPRG library*)


(* ::Subsection::Closed:: *)
(*Protecting important symbols*)


(* ::Text:: *)
(*Various symbols should not be overwritten as it might impede the program's function. They are marked as protected below.*)


(* ::Subsubsection::Closed:: *)
(*Flow equation constants and variables*)


Protect[d, vd, n, \[Epsilon], y];


$Assumptions = d>0;


(* ::Subsubsection::Closed:: *)
(*Effective action parameters*)


EffectiveActionParameters = {V, Zs, \[Delta]Zs, Zp, \[Delta]Zp, \[Eta]};


Map[Protect, EffectiveActionParameters];


(* ::Subsubsection::Closed:: *)
(*Regulator function*)


Protect[regulator, regulatorLabel, regulatorReplacement];


(* ::Subsubsection::Closed:: *)
(*Reserving symbols for eigenvalues*)


maxEigenvalueOrder = 20;
eigenvalueSymbols = Map[Subscript[e, #]&, Range[maxEigenvalueOrder]];
eigenvalueShortSymbols = Map[ToExpression["e" <> ToString[#]]&, Range[maxEigenvalueOrder]];
MapThread[(#1 = #2)&, {eigenvalueShortSymbols, eigenvalueSymbols}];
Protect[e];


(* ::Subsection::Closed:: *)
(*Flow equations*)


(* ::Subsubsection::Closed:: *)
(*DE2 guess data*)


(* ::Text:: *)
(*Guess values for the DE2 fixed-point solution*)


guessDE2v = {{0, -0.47770721436870717`}, {1/4, -0.4610871742592099`}, {1/2, -0.44410007864702195`}, {3/4, -0.42673862279138874`}, {1, -0.4089955148467867`}, {5/4, -0.3908634861318007`}, {3/2, -0.3723353014780051`}, {7/4, -0.3534037696042524`}, {2, -0.3340617534532111`}, {9/4, -0.31430218042626223`}, {5/2, -0.29411805245261874`}, {11/4, -0.2735024558290309`}, {3, -0.25244857076781646`}, {13/4, -0.23094968059318086`}, {7/2, -0.2089991805288659`}, {15/4, -0.18659058602403805`}, {4, -0.16371754056892654`}, {17/4, -0.14037382295698364`}, {9/2, -0.11655335395614569`}, {19/4, -0.09225020235803426`}, {5, -0.06745859038051583`}, {21/4, -0.04217289840582259`}, {11/2, -0.016387669043288557`}, {23/4, 0.009902389487445153`}, {6, 0.03670240065027947`}, {25/4, 0.06401731755475144`}, {13/2, 0.09185192157162007`}, {27/4, 0.12021082158172536`}, {7, 0.14909845383298767`}, {29/4, 0.17851908237610428`}, {15/2, 0.20847680004583466`}, {31/4, 0.2389755299517656`}, {8, 0.270019027440115`}, {33/4, 0.30161088248646356`}, {17/2, 0.3337545224782643`}, {35/4, 0.36645321534554615`}, {9, 0.39971007299834294`}, {37/4, 0.4335280550299986`}, {19/2, 0.46790997264656253`}, {39/4, 0.5028584927839342`}, {10, 0.5383761423761783`}, {41/4, 0.5744653127404583`}, {21/2, 0.6111282640462551`}, {43/4, 0.648367129838899`}, {11, 0.6861839215898979`}, {45/4, 0.7245805332490369`}, {23/2, 0.7635587457757155`}, {47/4, 0.8031202316294499`}, {12, 0.843266559201851`}, {49/4, 0.8839991971746866`}, {25/2, 0.9253195187908207`}, {51/4, 0.9672288060268744`}, {13, 1.0097282536583685`}, {53/4, 1.0528189732098754`}, {27/2, 1.0965019967843415`}, {55/4, 1.1407782807671842`}, {14, 1.18564870940213`}, {57/4, 1.231114098236908`}, {29/2, 1.2771751974379826`}, {59/4, 1.323832694974411`}, {15, 1.3710872196717077`}, {61/4, 1.4189393441372888`}, {31/2, 1.4673895875596297`}, {63/4, 1.516438418383765`}, {16, 1.5660862568661456`}, {65/4, 1.6163334775121918`}, {33/2, 1.6671804114001207`}, {67/4, 1.7186273483948225`}, {17, 1.7706745392556795`}, {69/4, 1.8233221976423164`}, {35/2, 1.876570502022309`}, {71/4, 1.9304195974848841`}, {18, 1.9848695974646249`}, {73/4, 2.039920585379146`}, {37/2, 2.09557261618464`}, {75/4, 2.151825717853095`}, {19, 2.2086798927749123`}, {77/4, 2.266135119090512`}, {39/2, 2.3241913519544144`}, {79/4, 2.3828485247351665`}, {20, 2.4421065501543295`}, {81/4, 2.5019653213676394`}, {41/2, 2.562424712991313`}, {83/4, 2.6234845820763275`}, {21, 2.685144769033394`}, {85/4, 2.7474050985112024`}, {43/2, 2.810265380230392`}, {87/4, 2.8737254097755915`}, {22, 2.9377849693477334`}, {89/4, 3.0024438284787616`}, {45/2, 3.067701744710708`}, {91/4, 3.1335584642410472`}, {23, 3.2000137225360983`}, {93/4, 3.26706724491417`}, {47/2, 3.334718747100058`}, {95/4, 3.4029679357523923`}, {24, 3.471814508965268`}, {97/4, 3.541258156745506`}, {49/2, 3.611298561466826`}, {99/4, 3.6819353983021177`}, {25, 3.753168335634965`}, {101/4, 3.8249970354514806`}, {51/2, 3.8974211537134655`}, {103/4, 3.9704403407138638`}, {26, 4.044054241415385`}, {105/4, 4.118262495773178`}, {53/2, 4.19306473904234`}, {107/4, 4.2684606020710225`}, {27, 4.344449711579865`}, {109/4, 4.421031690428417`}, {55/2, 4.498206157869195`}, {111/4, 4.575972729789992`}, {28, 4.654331018944996`}, {113/4, 4.733280635175253`}, {57/2, 4.8128211856190335`}, {115/4, 4.8929522749125045`}, {29, 4.973673505381257`}, {117/4, 5.054984477223045`}, {59/2, 5.136884788682229`}, {119/4, 5.219374036216201`}, {30, 5.302451814654306`}};
guessDE2zs = {{0, 1}, {1/4, 1.0032517770406784`}, {1/2, 1.006474314140323`}, {3/4, 1.0096648269307462`}, {1, 1.0128205752776709`}, {5/4, 1.0159388763030224`}, {3/2, 1.0190171172909137`}, {7/4, 1.02205276834321`}, {2, 1.025043394643287`}, {9/4, 1.0279866681881524`}, {5/2, 1.0308803788531502`}, {11/4, 1.0337224446601616`}, {3, 1.0365109211295382`}, {13/4, 1.0392440096078714`}, {7/2, 1.0419200644778257`}, {15/4, 1.0445375991724084`}, {4, 1.047095290933794`}, {17/4, 1.0495919842757737`}, {9/2, 1.052026693128592`}, {19/4, 1.054398601664873`}, {5, 1.0567070638250513`}, {21/4, 1.0589516015796931`}, {11/2, 1.061131901983944`}, {23/4, 1.0632478130955598`}, {6, 1.0652993388423435`}, {25/4, 1.0672866329369286`}, {13/2, 1.069209991946613`}, {27/4, 1.071069847633165`}, {7, 1.072866758682211`}, {29/4, 1.0746014019439398`}, {15/2, 1.0762745633066144`}, {31/4, 1.0778871283218474`}, {8, 1.0794400726960653`}, {33/4, 1.0809344527562657`}, {17/2, 1.0823713959903818`}, {35/4, 1.0837520917536145`}, {9, 1.085077782222236`}, {37/4, 1.086349753665996`}, {19/2, 1.0875693280995673`}, {39/4, 1.0887378553627793`}, {10, 1.089856705668914`}, {41/4, 1.090927262650277`}, {21/2, 1.0919509169207908`}, {43/4, 1.0929290601665964`}, {11, 1.0938630797677142`}, {45/4, 1.094754353946741`}, {23/2, 1.095604247434407`}, {47/4, 1.0964141076365699`}, {12, 1.0971852612828843`}, {49/4, 1.0979190115338826`}, {25/2, 1.0986166355205542`}, {51/4, 1.0992793822885654`}, {13, 1.099908471118019`}, {53/4, 1.1005050901890137`}, {27/2, 1.10107039556314`}, {55/4, 1.1016055104513776`}, {14, 1.1021115247395772`}, {57/4, 1.1025894947437183`}, {29/2, 1.1030404431683887`}, {59/4, 1.1034653592433785`}, {15, 1.1038651990148416`}, {61/4, 1.1042408857691348`}, {31/2, 1.1045933105691326`}, {63/4, 1.104923332884523`}, {16, 1.105231781299265`}, {65/4, 1.1055194542810263`}, {33/2, 1.1057871209990027`}, {67/4, 1.106035522178003`}, {17, 1.1062653709781163`}, {69/4, 1.1064773538905845`}, {35/2, 1.106672131641733`}, {71/4, 1.1068503400979388`}, {18, 1.1070125911656448`}, {73/4, 1.1071594736813533`}, {37/2, 1.1072915542873938`}, {75/4, 1.1074093782899865`}, {19, 1.1075134704968153`}, {77/4, 1.107604336031903`}, {39/2, 1.1076824611261118`}, {79/4, 1.1077483138820394`}, {20, 1.1078023450124683`}, {81/4, 1.107844988551888`}, {41/2, 1.10787666254086`}, {83/4, 1.1078977696832775`}, {21, 1.107908697976746`}, {85/4, 1.1079098213164884`}, {43/2, 1.1079015000733254`}, {87/4, 1.1078840816463673`}, {22, 1.1078579009911598`}, {89/4, 1.1078232811240767`}, {45/2, 1.1077805336038102`}, {91/4, 1.1077299589908212`}, {23, 1.1076718472856568`}, {93/4, 1.1076064783470245`}, {47/2, 1.1075341222905213`}, {95/4, 1.1074550398689111`}, {24, 1.1073694828348153`}, {97/4, 1.1072776942866713`}, {49/2, 1.1071799089987886`}, {99/4, 1.1070763537362922`}, {25, 1.1069672475557366`}, {101/4, 1.1068528020921187`}, {51/2, 1.1067332218330024`}, {103/4, 1.106608704380428`}, {26, 1.106479440701251`}, {105/4, 1.1063456153665192`}, {53/2, 1.1062074067804681`}, {107/4, 1.1060649873996835`}, {27, 1.1059185239429523`}, {109/4, 1.105768177592284`}, {55/2, 1.1056141041855758`}, {111/4, 1.1054564544013492`}, {28, 1.1052953739359677`}, {113/4, 1.1051310036737334`}, {57/2, 1.1049634798502057`}, {115/4, 1.1047929342091087`}, {29, 1.1046194941531209`}, {117/4, 1.104443282888877`}, {59/2, 1.1042644195664464`}, {119/4, 1.1040830194135514`}, {30, 1.1038991938648042`}};
guessDE2zp = {{0, 1}, {1/4, 1.0003400209628095`}, {1/2, 1.00066517380496`}, {3/4, 1.0009752194036339`}, {1, 1.0012699357785124`}, {5/4, 1.0015491191355352`}, {3/2, 1.0018125848553683`}, {7/4, 1.0020601684163568`}, {2, 1.0022917262425615`}, {9/4, 1.0025071364681186`}, {5/2, 1.0027062996097946`}, {11/4, 1.0028891391403245`}, {3, 1.0030556019559878`}, {13/4, 1.0032056587328315`}, {7/2, 1.0033393041670426`}, {15/4, 1.0034565570961336`}, {4, 1.003557460498852`}, {17/4, 1.0036420813730016`}, {9/2, 1.003710510491697`}, {19/4, 1.0037628620398646`}, {5, 1.0037992731341199`}, {21/4, 1.0038199032303754`}, {11/2, 1.003824933424716`}, {23/4, 1.0038145656541608`}, {6, 1.003789021804892`}, {25/4, 1.0037485427363837`}, {13/2, 1.0036933872305618`}, {27/4, 1.0036238308756822`}, {7, 1.0035401648950162`}, {29/4, 1.003442694930674`}, {15/2, 1.0033317397929966`}, {31/4, 1.0032076301858814`}, {8, 1.0030707074182197`}, {33/4, 1.0029213221112858`}, {17/2, 1.0027598329114953`}, {35/4, 1.0025866052173904`}, {9, 1.002402009929094`}, {37/4, 1.0022064222277771`}, {19/2, 1.0020002203919323`}, {39/4, 1.0017837846564719`}, {10, 1.0015574961198586`}, {41/4, 1.0013217357036732`}, {21/2, 1.001076883168232`}, {43/4, 1.000823316187075`}, {11, 1.0005614094824113`}, {45/4, 1.000291534022896`}, {23/2, 1.0000140562844464`}, {47/4, 0.9997293375742118`}, {12, 0.9994377334172387`}, {49/4, 0.9991395930048924`}, {25/2, 0.9988352587036525`}, {51/4, 0.9985250656225257`}, {13, 0.9982093412369945`}, {53/4, 0.9978884050671609`}, {27/2, 0.9975625684075297`}, {55/4, 0.9972321341057094`}, {14, 0.9968973963871998`}, {57/4, 0.9965586407233495`}, {29/2, 0.9962161437395409`}, {59/4, 0.9958701731606371`}, {15, 0.9955209877907715`}, {61/4, 0.9951688375245932`}, {31/2, 0.9948139633871593`}, {63/4, 0.9944565975997542`}, {16, 0.9940969636690172`}, {65/4, 0.9937352764968768`}, {33/2, 0.9933717425089098`}, {67/4, 0.9930065597988763`}, {17, 0.9926399182873124`}, {69/4, 0.9922719998921982`}, {35/2, 0.9919029787098524`}, {71/4, 0.9915330212043418`}, {18, 0.9911622864038221`}, {73/4, 0.9907909261023565`}, {37/2, 0.9904190850658793`}, {75/4, 0.9900469012410953`}, {19, 0.9896745059662156`}, {77/4, 0.9893020241825351`}, {39/2, 0.9889295746459674`}, {79/4, 0.9885572701377398`}, {20, 0.9881852176735435`}, {81/4, 0.9878135187105214`}, {41/2, 0.9874422693515446`}, {83/4, 0.987071560546309`}, {21, 0.9867014782888449`}, {85/4, 0.9863321038110927`}, {43/2, 0.9859635137722564`}, {87/4, 0.9855957804436927`}, {22, 0.985228971889144`}, {89/4, 0.9848631521401598`}, {45/2, 0.9844983813665962`}, {91/4, 0.9841347160421048`}, {23, 0.9837722091045661`}, {93/4, 0.9834109101114389`}, {47/2, 0.9830508653900236`}, {95/4, 0.9826921181826631`}, {24, 0.9823347087869139`}, {97/4, 0.9819786746907442`}, {49/2, 0.9816240507028273`}, {99/4, 0.9812708690780061`}, {25, 0.9809191596380208`}, {101/4, 0.9805689498875971`}, {51/2, 0.9802202651259989`}, {103/4, 0.9798731285541545`}, {26, 0.9795275613774754`}, {105/4, 0.9791835829044805`}, {53/2, 0.97884121064135`}, {107/4, 0.978500460382531`}, {27, 0.978161346297516`}, {109/4, 0.9778238810139204`}, {55/2, 0.9774880756969808`}, {111/4, 0.9771539401255974`}, {28, 0.9768214827650414`}, {113/4, 0.9764907108364475`}, {57/2, 0.9761616303832092`}, {115/4, 0.9758342463343939`}, {29, 0.9755085625652894`}, {117/4, 0.9751845819551914`}, {59/2, 0.974862306442551`}, {119/4, 0.9745417370775632`}, {30, 0.9742228740723354`}};


(* ::Subsubsection::Closed:: *)
(*FP guesses*)


(* ::Text:: *)
(*Constructing guesses for LPA and DE2 fixed-point solutions*)


guessMinimum = 4.5;
guessLPA := Table[{V[i], -0.2 (1 - 2 i/gridSize \[Rho]Max/guessMinimum)}, {i, 0, gridSize}];


guessDE2vInterpolation = Interpolation[guessDE2v];
guessDE2zsInterpolation = Interpolation[guessDE2zs];
guessDE2zpInterpolation = Interpolation[guessDE2zp];


guessDE2 := Quiet[Join[{{\[Eta], 45/1000}}, 
Table[{V[i], guessDE2vInterpolation[i*\[Rho]Max/gridSize]}, {i, 0, gridSize}], 
Table[{Zs[i], guessDE2zsInterpolation[i*\[Rho]Max/gridSize]}, {i, 1, gridSize}], 
Table[{Zp[i], guessDE2zpInterpolation[i*\[Rho]Max/gridSize]}, {i, 1, gridSize}]]];


guessDE2\[Delta] := Join[{{\[Eta], 45/1000}}, 
Table[{V[i], guessDE2vInterpolation[i*\[Rho]Max/gridSize]}, {i, 0, gridSize}], 
Table[{\[Delta]Zs[i], guessDE2zsInterpolation[i*\[Rho]Max/gridSize]-1}, {i, 1, gridSize}], 
Table[{\[Delta]Zp[i], guessDE2zpInterpolation[i*\[Rho]Max/gridSize]-1}, {i, 1, gridSize}]];


(* ::Subsubsection::Closed:: *)
(*LPA flow equations*)


(* ::Text:: *)
(*Flow equations at the LPA level. Flow equations form a list of 4-element lists; each entry corresponds to a separate parametrizing function.*)
(*The entries for each function denotes:*)
(*1. The symbol for the function, e.g. V[\[Rho]/\[Epsilon]] (\[Rho]/\[Epsilon] is an integer i denoting a grid point)*)
(*2. The rescaling contribution to the flow of the function*)
(*3. The \[Beta] function contribution to the flow*)
(*4. The \[Rho]=0 evaluation of the \[Beta] function (this element surpasses a problem with apparent nonanalyticities of form A/\[Rho], which have to be resolved analytically before the numerical treatment).*)


flowEquationsLPA ={{V[\[Rho]/\[Epsilon]], 2 V[\[Rho]/\[Epsilon]]-((-2+d) \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon], 4 vd y^(-1+d) (((2 \[Alpha] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) Derivative[1][V][\[Rho]/\[Epsilon]](n-1))/(2 \[Epsilon] (y^2+\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]])^2)+((2 \[Alpha] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/(2 (y^2+\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^2)), (2 vd y^(-1+d) (2 \[Alpha] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) Derivative[1][V][0] (n+2))/(\[Epsilon] (y^2+\[Alpha] r[y^2]+V[0])^2)}};


(* ::Subsubsection::Closed:: *)
(*DE2 flow equations*)


(* ::Text:: *)
(*Flow equations at the DE2 level. Flow equations form a list of 4-element lists; each entry corresponds to a separate parametrizing function.*)
(*The entries for each function denotes:*)
(*1. The symbol for the function, e.g. V[\[Rho]/\[Epsilon]] (\[Rho]/\[Epsilon] is an integer i denoting a grid point)*)
(*2. The rescaling contribution to the flow of the function*)
(*3. The \[Beta] function contribution to the flow*)
(*4. The \[Rho]=0 evaluation of the \[Beta] function (this element surpasses a problem with apparent nonanalyticities of form A/\[Rho], which have to be resolved analytically before the numerical treatment).*)


flowEquationsDE2 = {{V[\[Rho]/\[Epsilon]], -((-2+\[Eta]) V[\[Rho]/\[Epsilon]])-((-2+d+\[Eta]) \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon], -4 vd y^(-1+d) (-(((-1+n) (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]))/(2 (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^2))-((2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/(2 (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^2)), -((4 vd y^(-1+d) (1/2 (1-n) (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) (Derivative[1][V][0]/\[Epsilon]+(y^2 Derivative[1][Zp][0])/\[Epsilon])-1/2 (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) ((3 Derivative[1][V][0])/\[Epsilon]+(y^2 Derivative[1][Zs][0])/\[Epsilon])))/(\[Alpha] r[y^2]+V[0]+y^2 Zs[0])^2)}, {Zs[\[Rho]/\[Epsilon]], -\[Eta] Zs[\[Rho]/\[Epsilon]]-((-2+d+\[Eta]) \[Rho] Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon], -(1/\[Rho])2 vd y^(-1+d) (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) (((4 (-1+n) y^2 \[Rho]^2 Derivative[1][Zp][\[Rho]/\[Epsilon]]^2)/(d \[Epsilon]^2)+4 (-1+n) \[Rho] (-Zp[\[Rho]/\[Epsilon]]+Zs[\[Rho]/\[Epsilon]]) (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]))/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^3+(16 (-1+n) y^2 \[Rho]^2 Zp[\[Rho]/\[Epsilon]]^2 (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])^2+32 (-1+n) y^2 \[Alpha] \[Rho]^2 Zp[\[Rho]/\[Epsilon]] Derivative[1][r][y^2] (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])^2+16 (-1+n) y^2 \[Alpha]^2 \[Rho]^2 Derivative[1][r][y^2]^2 (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])^2)/(d (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^5)-((-1+n) (Zp[\[Rho]/\[Epsilon]]-Zs[\[Rho]/\[Epsilon]]+(\[Rho] Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]))/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^2+1/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^4 (-4 (-1+n) \[Rho]^2 Zp[\[Rho]/\[Epsilon]] (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])^2-4 (-1+n) \[Alpha] \[Rho]^2 Derivative[1][r][y^2] (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])^2+(-((16 (-1+n) y^2 \[Rho]^2 Zp[\[Rho]/\[Epsilon]] Derivative[1][Zp][\[Rho]/\[Epsilon]] (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]))/\[Epsilon])-(16 (-1+n) y^2 \[Alpha] \[Rho]^2 Derivative[1][r][y^2] Derivative[1][Zp][\[Rho]/\[Epsilon]] (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]))/\[Epsilon]-8 (-1+n) y^2 \[Alpha] \[Rho]^2 (Derivative[1][V][\[Rho]/\[Epsilon]]/\[Epsilon]+(y^2 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])^2 Derivative[2][r][y^2])/d)+((4 y^2 \[Rho]^2 Derivative[1][Zs][\[Rho]/\[Epsilon]]^2)/(d \[Epsilon]^2)+(8 \[Rho]^2 Derivative[1][Zs][\[Rho]/\[Epsilon]] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/\[Epsilon])/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^3+(16 y^2 \[Rho]^2 Zs[\[Rho]/\[Epsilon]]^2 ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2)^2+32 y^2 \[Alpha] \[Rho]^2 Zs[\[Rho]/\[Epsilon]] Derivative[1][r][y^2] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2)^2+16 y^2 \[Alpha]^2 \[Rho]^2 Derivative[1][r][y^2]^2 ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2)^2)/(d (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^5)+1/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^4 (-4 \[Rho]^2 Zs[\[Rho]/\[Epsilon]] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2)^2-4 \[Alpha] \[Rho]^2 Derivative[1][r][y^2] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2)^2+(-((16 y^2 \[Rho]^2 Zs[\[Rho]/\[Epsilon]] Derivative[1][Zs][\[Rho]/\[Epsilon]] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/\[Epsilon])-(16 y^2 \[Alpha] \[Rho]^2 Derivative[1][r][y^2] Derivative[1][Zs][\[Rho]/\[Epsilon]] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/\[Epsilon]-8 y^2 \[Alpha] \[Rho]^2 Derivative[2][r][y^2] ((3 Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon]+(y^2 Derivative[1][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][V][\[Rho]/\[Epsilon]])/\[Epsilon]^2)^2)/d)-(\[Rho] (Derivative[1][Zs][\[Rho]/\[Epsilon]]/\[Epsilon]+(2 \[Rho] Derivative[2][Zs][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^2), (2 vd y^(-1+d) (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) (-(Derivative[1][Zp][0]/\[Epsilon])+(n Derivative[1][Zp][0])/\[Epsilon]+Derivative[1][Zs][0]/\[Epsilon]))/(\[Alpha] r[y^2]+V[0]+y^2 Zs[0])^2}, {Zp[\[Rho]/\[Epsilon]], -\[Eta] Zp[\[Rho]/\[Epsilon]]-((-2+d+\[Eta]) \[Rho] Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon], -(1/\[Rho])2 vd y^(-1+d) (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) ((4 y^2 Zp[\[Rho]/\[Epsilon]]^2+8 y^2 \[Alpha] Zp[\[Rho]/\[Epsilon]] Derivative[1][r][y^2]+4 y^2 \[Alpha]^2 Derivative[1][r][y^2]^2)/(d (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^3)+(4 y^2 Zs[\[Rho]/\[Epsilon]]^2+8 y^2 \[Alpha] Zs[\[Rho]/\[Epsilon]] Derivative[1][r][y^2]+4 y^2 \[Alpha]^2 Derivative[1][r][y^2]^2)/(d (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^3)+(-2 Zp[\[Rho]/\[Epsilon]]-2 \[Alpha] Derivative[1][r][y^2]-((-1+n) \[Rho] Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]-(4 y^2 \[Alpha] Derivative[2][r][y^2])/d)/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^2+((-8 y^2 \[Alpha] Zp[\[Rho]/\[Epsilon]] Derivative[1][r][y^2]-4 y^2 \[Alpha]^2 Derivative[1][r][y^2]^2-4 y^2 (Zp[\[Rho]/\[Epsilon]]^2-(\[Rho]^2 Derivative[1][Zp][\[Rho]/\[Epsilon]]^2)/\[Epsilon]^2))/(d (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]])^2)+(4 \[Alpha] Derivative[1][r][y^2]+4 (Zp[\[Rho]/\[Epsilon]]+(\[Rho] Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon])+(8 y^2 \[Alpha] Derivative[2][r][y^2])/d)/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]]))/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])+(-Zp[\[Rho]/\[Epsilon]]-Zs[\[Rho]/\[Epsilon]]-2 \[Alpha] Derivative[1][r][y^2]+(-8 y^2 \[Alpha] Zs[\[Rho]/\[Epsilon]] Derivative[1][r][y^2]-4 y^2 \[Alpha]^2 Derivative[1][r][y^2]^2+4 y^2 (Zp[\[Rho]/\[Epsilon]]+(\[Rho] Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]) (Zp[\[Rho]/\[Epsilon]]-2 Zs[\[Rho]/\[Epsilon]]+(\[Rho] Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]))/(d (\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zp[\[Rho]/\[Epsilon]]))-(4 y^2 \[Alpha] Derivative[2][r][y^2])/d-\[Rho] ((5 Derivative[1][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]+(2 \[Rho] Derivative[2][Zp][\[Rho]/\[Epsilon]])/\[Epsilon]^2))/(\[Alpha] r[y^2]+V[\[Rho]/\[Epsilon]]+y^2 Zs[\[Rho]/\[Epsilon]]+(2 \[Rho] Derivative[1][V][\[Rho]/\[Epsilon]])/\[Epsilon])^2), (2 vd y^(-1+d) (2 \[Alpha] r[y^2]-\[Alpha] \[Eta] r[y^2]-2 y^2 \[Alpha] Derivative[1][r][y^2]) (-(Derivative[1][Zp][0]/\[Epsilon])+(n Derivative[1][Zp][0])/\[Epsilon]+Derivative[1][Zs][0]/\[Epsilon]))/(\[Alpha] r[y^2]+V[0]+y^2 Zs[0])^2}};


(* ::Text:: *)
(*flowEquationsDE2\[Delta] parameterize Zs and Zp via their difference with unity. This slightly improves the numerical accuracy.*)


flowEquationsDE2\[Delta] = flowEquationsDE2;
flowEquationsDE2\[Delta][[All, 2;;-1]] = flowEquationsDE2\[Delta][[All, 2;;-1]] /. {Zs -> (1+\[Delta]Zs[#]&), Zp -> (1+\[Delta]Zp[#]&)};
flowEquationsDE2\[Delta][[All, 1]] = flowEquationsDE2\[Delta][[All, 1]] /. {Zs -> (\[Delta]Zs[#]&), Zp -> (\[Delta]Zp[#]&)};


(* ::Subsubsection::Closed:: *)
(*Replacing Constants*)


(* ::Text:: *)
(*Function for replacing standard constants appearing in the flow equations*)


ReplaceConstants[expression_] := expression /. {\[Alpha] -> alpha, d -> dimension, vd -> 1, 
	\[Epsilon] -> \[Rho]Max/gridSize, Zs[0] -> 1, Zp[0] -> 1, \[Delta]Zs[0] -> 0, \[Delta]Zp[0] -> 0, n -> NValue};


(* ::Subsection::Closed:: *)
(*Infrared regulator functions*)


(* ::Subsubsection::Closed:: *)
(*Defining replacement rules for regulators*)


(* ::Text:: *)
(*Standard regulator functions. In this program the dimensionfull regulator takes the form R[k, q^2] = Z k^2 \[Alpha] r[q^2/k^2], note that unlike in the typical convention \[Alpha] is moved outside of the dimensionless function r.*)


litimRegulator = {\[Alpha] -> 1, r -> ((1-#)&)};
litim2Regulator = {\[Alpha] -> 1, r -> ((1-#)^2&)};
exponentialRegulator = {r -> (Exp[-#]&)};
wetterichRegulator = {r -> (#/(Exp[#]-1)&)};


(* ::Subsubsection::Closed:: *)
(*Selecting regulator function *)


(* ::Text:: *)
(*Similarly to the integrating functions, the regulator is protected and can only be changed via the SelectRegulator function.*)


SelectRegulator::InvalidArgument := "`1` is not a valid regulator name, no regulator function has been selected.";
SelectRegulator[regulatorName_]:= Block[{litimNames, exponentialNames, wetterichNames, litim2Names}, 
	litimNames = {"litim", "theta"};
	exponentialNames = {"exponential"};
	wetterichNames = {"wetterich", "smooth"};
	litim2Names = {"litim2"};

	Unprotect[regulator, regulatorLabel, regulatorReplacement];
	If[!MemberQ[Join[litimNames, exponentialNames, wetterichNames, litim2Names], ToLowerCase[regulatorName]], 
		Message[SelectRegulator::InvalidArgument, regulatorName];
		regulator = Null;
		regulatorLabel = Null;
	];

	If[MemberQ[litimNames, ToLowerCase[regulatorName]], 
		regulator = litimRegulator;
		regulatorLabel = "Litim";
	];
	
	If[MemberQ[exponentialNames, ToLowerCase[regulatorName]], 
		regulator = exponentialRegulator;
		regulatorLabel = "Exponential"; 
	];
	If[MemberQ[litim2Names, ToLowerCase[regulatorName]], 
		regulator = litim2Regulator;
		regulatorLabel = "Litim2"; 
	];
	
	If[MemberQ[wetterichNames, ToLowerCase[regulatorName]], 
		regulator = wetterichRegulator;
		regulatorLabel = "Wetterich"; 
	];

	If[ValueQ[regulator] && !TrueQ[regulator == Null], 
		PrintLog[regulatorLabel <> " regulator selected"];
	];
	
	regulatorReplacement = {r[\[Xi]_] -> (Limit[r[\[Zeta]] /. regulator, \[Zeta] -> \[Xi]]), 
						Derivative[n_][r][\[Xi]_] -> (Limit[Derivative[n][r][\[Zeta]] /. regulator, \[Zeta] -> \[Xi]])};

	Protect[regulator, regulatorLabel, regulatorReplacement];
]


(* ::Subsection::Closed:: *)
(*Processing flow equations*)


(* ::Subsubsection::Closed:: *)
(*Discretization of flow equations*)


(* ::Text:: *)
(*Converts a list of functional flow into their discrete grid representation.*)
(*Returns a three-element list:*)
(*	1. List of effective-action parameters in the grid representation*)
(*	2. List of rescaling terms in time derivative of parameters in the grid representation*)
(*	3. List of \[Beta] functions in the grid representation (note that \[Beta] functions can be integrated or not at this point)*)


DiscretizeFlowEquations::noRegulator := "Regulator function has not been properly selected"
DiscretizeFlowEquations::noGrid := "gridSize is not properly specified"
DiscretizeFlowEquations::integralError := "Parameters of numerical integrals are not properly specified"
DiscretizeFlowEquations[flowEquations_, include\[Eta]_:False] := Block[{ ProcessOneParameter, discreteEquations, allParameters, log, allScalingTerms, allIntegrands}, 

	If[!ValueQ[regulatorLabel] || TrueQ[regulatorLabel == Null], 
		Message[DiscretizeFlowEquations::noRegulator];
		Return[Null];
	];

	If[!IntegerQ[gridSize] || gridSize < 5, 
		Message[DiscretizeFlowEquations::noGrid];
		Return[Null];
	];

	If[regulatorLabel != "Litim" && 
		(!NumberQ[GLIntegralUpperBound] || 
		!IntegerQ[GLIntegralPointCount] || GLIntegralPointCount < 0), 

		Message[DiscretizeFlowEquations::integralError];
		Return[Null];
	];

	ProcessOneParameter[flowEq_] := Block[{parameter, gridPoints, discreteParameters, scaling, scalingTerms, integrand, 
		discreteIntegrands, fullFlow, rho0Integrand}, 
		
		{parameter, scaling, integrand, rho0Integrand} = flowEq;
		
		If[MemberQ[{Zs[\[Rho]/\[Epsilon]], Zp[\[Rho]/\[Epsilon]], \[Delta]Zs[\[Rho]/\[Epsilon]], \[Delta]Zp[\[Rho]/\[Epsilon]]}, parameter], 
			gridPoints = Range[1, gridSize], 
			gridPoints = Range[0, gridSize];
		];

		discreteParameters = Map[(parameter /. \[Rho] -> # \[Epsilon])&, gridPoints] ;
		scalingTerms = Map[NumericDerivatives[(scaling /. \[Rho] -> # \[Epsilon])]&, gridPoints];
		discreteIntegrands = Map[NumericDerivatives[(integrand /. \[Rho] -> # \[Epsilon])]&, gridPoints];

		If[include\[Eta] && (parameter == Zs[\[Rho]/\[Epsilon]] || parameter == \[Delta]Zs[\[Rho]/\[Epsilon]]), 
			discreteParameters = Join[{\[Eta]}, discreteParameters];
			scalingTerms = Join[{NumericDerivatives[(scaling /. \[Rho] -> 0)]}, scalingTerms];
			discreteIntegrands = Join[{NumericDerivatives[rho0Integrand]}, discreteIntegrands];
		];

		Return[{discreteParameters, scalingTerms, discreteIntegrands}];
	];
 
	discreteEquations = Map[ProcessOneParameter, flowEquations];

	allParameters = Flatten[discreteEquations[[All, 1]]];
	allScalingTerms = Flatten[discreteEquations[[All, 2]]];
	allIntegrands = Flatten[discreteEquations[[All, 3]]];

	log = "Flow equations discretized on a " <> ToString[gridSize+1] <> "-point grid, with " <> 
		ToString[derivativePointsCount] <> "-point derivatives.";

	PrintLog[log];

	Return[{allParameters, allScalingTerms, allIntegrands}];
]


(* ::Subsubsection::Closed:: *)
(*Performing loop integrals*)


(* ::Text:: *)
(*Performs the loop integrals in the \[Beta] functions. *)
(*Important: \[Beta] functions can be integrated before or after discretization. If they're integrated before discretization the integral has to be performed separately for the general \[Beta] functions ('flowEquationsLPA[[All, 3]]', 'flowEquationsDE2[[All, 3]]') and their \[Rho]=0 forms ('flowEquationsLPA[[All, 4]]', 'flowEquationsDE2[[All, 4]]').*)
(*The integrals can are performed analytically when the Litim regulator is selected and the argument 'analyticalIntegrals' is True.*)


IntegrateLitim[integrand_] := Integrate[Simplify[integrand /. litimRegulator, Assumptions -> 0<y<1], {y, 0, 1}]


PerformMomentumIntegrals[integrands_, analyticalIntegrals_:False] := Block[{integratedLoop}, 
	
	If[regulatorLabel == "Litim" && analyticalIntegrals, 
		integratedLoop = IntegrateLitim[integrands];
		PrintLog["Momentum integrals performed analytically."], 

		integratedLoop = GaussLegendreIntegrate[integrands];
		(*integratedLoop = SimpsonIntegrate[integrands];*)
		PrintLog["Momentum integrals performed via the " <> 
			ToString[GLIntegralPointCount] <> "-point Gauss-Legendre quadrature on the interval y \[Element] [0, " <> 
			ToString[GLIntegralUpperBound] <> "]."];
	];
	Return[integratedLoop];
]


(* ::Subsection::Closed:: *)
(*Fixed point equation solving*)


(* ::Subsubsection::Closed:: *)
(*Find fixed point*)


(* ::Text:: *)
(*Solves the discretized fixed-point equations starting from the 'guess'.*)


FindFixedPoint[discreteEquations_, guess_] := Block[{numericalEquations, solution, equationPrecision, precision, guessList}, 
	guessList = Transpose[{guess[[All, 1]], guess[[All, 2]]}];
	precision = $MinPrecision;
	If[precision <= 0, 
		precision = MachinePrecision];
		
	Map[Unprotect, EffectiveActionParameters];
	numericalEquations = ReplaceConstants[discreteEquations];
	solution = Quiet[FindRoot[numericalEquations, guessList, 
			WorkingPrecision -> precision, AccuracyGoal -> precision-2, PrecisionGoal -> precision-2, MaxIterations -> 1000, DampingFactor -> .1], 
		{FindRoot::precw, N::meprec, FindRoot::bddir, FindRoot::lstol}];
	Map[Protect, EffectiveActionParameters];

	Return[solution];
]


(* ::Subsubsection::Closed:: *)
(*Plot effective action*)


(* ::Text:: *)
(*Creates a simple plot visualizing the effective-action parameters.*)


PlotEffectiveAction[action_] := Block[{actionReplacementList, parameters, symbols = {V, Zs, Zp}, SerializeParameter, reparametrization, LPAReplacement}, 
	LPAReplacement = {\[Eta] -> 0, Zs[i_] -> 1, Zp[i_] -> 1, \[Delta]Zs[i_] -> 0, \[Delta]Zp[i_] -> 0};
	reparametrization = {Zs[i_] -> 1 + \[Delta]Zs[i], Zp[i_] -> 1 + \[Delta]Zp[i]};
	actionReplacementList =Join[Map[#[[1]] -> #[[2]]&, action], {\[Epsilon] -> \[Rho]Max/gridSize}];
	SerializeParameter[parameterSymbol_] := Map[{\[Epsilon] #, parameterSymbol[#]}&, Range[0, gridSize]] /. actionReplacementList 
		 /. reparametrization /. actionReplacementList /. LPAReplacement;
	parameters = Map[SerializeParameter, symbols];

	Return[ListPlot[parameters, PlotLegends -> symbols, AxesLabel -> {"\[Rho]", ""}, 
		PlotLabel -> "\[Eta]=" <> ToString[NumberForm[\[Eta] /. actionReplacementList /. \[Eta] -> 0, {4, 3}]], PlotRange -> All]];
]


(* ::Subsubsection::Closed:: *)
(*Get potential minimum*)


(* ::Text:: *)
(*Returns the minimum of the effective action*)


GetPotentialMinimum[action_] := Block[{actionReplacementList, VInterpolation, solution}, 
	actionReplacementList = Join[Map[#[[1]] -> #[[2]]&, action], {\[Epsilon] -> \[Rho]Max/gridSize}];
	VInterpolation = Interpolation[Table[{i \[Epsilon], V[i]}, {i, 0, gridSize}] /. actionReplacementList];
	solution = FindRoot[VInterpolation[\[Rho]] == 0, {\[Rho], \[Rho]Max/2, 0, \[Rho]Max}, DampingFactor -> 2];
	Return[\[Rho] /. solution[[1]]];
]


(* ::Subsubsection::Closed:: *)
(*Find potential minimum*)


(* ::Text:: *)
(*Solves the fixed point equations with low accuracy to quickly determine an approximation for the potential minimum \[Rho]0.*)


FindPotentialMinimum[nn_, DEOrder_, regulatorLabel_, analyticalIntegrals_] := Block[{
	\[Rho]Max = 16, gridSize = 40, NValue = nn, dimension = 3, alpha = 1, derivativePointsCount = 7, 
	regulator, GLIntegralEvaluationPoints, GLIntegralSmallParameter, GLIntegralReplacementRules, 
	parameters, scalingTerms, discreteEquations, integratedEquations, fpSolution, potentialMinimum, 
	fpGuess, flowEquations
	}, 
	If[DEOrder == "LPA", 
		flowEquations = flowEquationsLPA;
		fpGuess = guessLPA, 
		flowEquations = flowEquationsDE2;
		fpGuess = guessDE2;
	];

	SelectRegulator[regulatorLabel];
	SetGLIntegralParameters[GLIntegralLowPrecisionParameters];
	SetFixedFloatPrecision[doublePrecision];
	
	integratedEquations = flowEquations;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], analyticalIntegrals];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], analyticalIntegrals];
	
	{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];
	fpSolution = FindFixedPoint[scalingTerms + discreteEquations, fpGuess];
	potentialMinimum = GetPotentialMinimum[fpSolution];
	
	RestoreArbitraryPrecision[];
	Return[RecastFloatingPoint[potentialMinimum]];
]


(* ::Subsection::Closed:: *)
(*Stability matrix*)


(* ::Subsubsection::Closed:: *)
(*Calculated analytically*)


(* ::Text:: *)
(*Calculates the stability matrix from 'flowEquations' presented in the functional form around the FP solution 'fixedPoint'.*)
(*The integrals can are performed analytically when the Litim regulator is selected and the argument 'analyticalIntegrals' is True.*)
(*The derivatives are performed analytically if 'epsilon' == 0 or 'oneSideDerivativePointsCount' <= 0. Otherwise, the derivatives are approximated with finite differences.*)


StabilityMatrix[flowEquations_, fixedPoint_, analyticalIntegrals_:False, oneSideDerivativePointsCount_:0, epsilon_:0] := 
	Block[{parameters, scalingTerms, discreteEquations, integratedLoop, flowEquationsIntegrated, loopMatrix, scalingMatrix, gradientFunction}, 
	gradientFunction = Grad;
	If[oneSideDerivativePointsCount > 0 && epsilon != 0, 
		gradientFunction = NumericGradient[#1, #2, oneSideDerivativePointsCount, epsilon]&;
	];
	
	If[TrueQ[regulatorLabel == "Litim"] && analyticalIntegrals && Length[flowEquations] == 1, 
		Return[StabilityMatrixAnalyticalIntegrals[flowEquations, fixedPoint, gradientFunction]], 
		Return[StabilityMatrixNumericalIntegrals[flowEquations, fixedPoint, gradientFunction]]
	];
]


StabilityMatrixAnalyticalIntegrals[flowEquations_, fixedPoint_, gradientFunction_]:= 
	Block[{parameters, scalingTerms, discreteEquations, integratedLoop, flowEquationsIntegrated, loopMatrix, scalingMatrix}, 
	flowEquationsIntegrated = flowEquations;
	flowEquationsIntegrated[[All, 3]] = PerformMomentumIntegrals[flowEquationsIntegrated[[All, 3]], True];
	flowEquationsIntegrated[[All, 4]] = PerformMomentumIntegrals[flowEquationsIntegrated[[All, 4]], True];

	{parameters, scalingTerms, discreteEquations} = ReplaceConstants[DiscretizeFlowEquations[flowEquationsIntegrated, False]];

	Return[gradientFunction[scalingTerms + discreteEquations, parameters] /. fixedPoint];
];


StabilityMatrixNumericalIntegrals[flowEquations_, fixedPoint_, gradientFunction_] := 
	Block[{parameters, scalingTerms, discreteEquations, integratedLoop, flowEquationsIntegrated, 
		loopMatrix, scalingMatrix, loopMatrixIntegrated, 
		\[Eta]Search, \[Eta]Found, \[Eta]Index, \[Eta]Equation, \[Eta]Solution, \[Eta]Gradient, \[Eta]ScalingMatrix, \[Eta]LoopMatrix, \[Eta]Matrix}, 
	
	{parameters, scalingTerms, discreteEquations} = ReplaceConstants[DiscretizeFlowEquations[flowEquations, True]];

	\[Eta]Search = Position[parameters, \[Eta]];
	\[Eta]Found = Length[\[Eta]Search] > 0;
	If[\[Eta]Found, 
		\[Eta]Index = Flatten[\[Eta]Search][[1]];
		\[Eta]Equation = scalingTerms[[\[Eta]Index]] + PerformMomentumIntegrals[discreteEquations[[\[Eta]Index]], False];
		\[Eta]Solution = \[Eta] /. Flatten[Solve[0 == \[Eta]Equation, \[Eta]]][[1]];
		parameters = Drop[parameters, {\[Eta]Index}];
		scalingTerms = Drop[scalingTerms, {\[Eta]Index}];
		discreteEquations = Drop[discreteEquations, {\[Eta]Index}];

		\[Eta]Gradient = gradientFunction[\[Eta]Solution, parameters] /. fixedPoint;
		
		\[Eta]ScalingMatrix = TensorProduct[D[scalingTerms, \[Eta]], \[Eta]Gradient] /. fixedPoint;
		\[Eta]LoopMatrix = TensorProduct[PerformMomentumIntegrals[D[discreteEquations, \[Eta]] /. fixedPoint, False], \[Eta]Gradient];
		\[Eta]Matrix = \[Eta]LoopMatrix + \[Eta]ScalingMatrix, 

		\[Eta]Matrix = Table[0, {i, 1, Length[parameters]}, {j, 1, Length[parameters]}];
	];

	loopMatrix = gradientFunction[discreteEquations, parameters] /. fixedPoint;
	loopMatrixIntegrated = PerformMomentumIntegrals[loopMatrix, False];

	scalingMatrix = gradientFunction[scalingTerms, parameters] /. fixedPoint;

	Return[scalingMatrix + loopMatrixIntegrated + \[Eta]Matrix];
]


(* ::Subsubsection::Closed:: *)
(*Selecting eigenvalues*)


(* ::Text:: *)
(*Generates a list of 'count' of leading eigenvalues of 'matrix' sorted in the descending order with respect to the real part and subsequently to the imaginary part.*)


LeadingEigenvalues[matrix_, count_:20] := Sort[Eigenvalues[matrix, -count], Re[#1]>Re[#2] || (Re[#1] == Re[#2] && Im[#1]>Im[#2]) &];


(* ::Section:: *)
(*Calculations*)


(* ::Subsection::Closed:: *)
(*Tracking the compactification error (\[Rho]Max dependence)*)


(* ::Subsubsection::Closed:: *)
(*LPA*)


Track\[Rho]MaxErrorDependenceLPA[nn_]:=Block[{NValue = nn, dimension = 3, alpha = 1, derivativePointsCount = 7, 
	gridSpacing = 1/20, optimal\[Rho]Max = 35/10, 
	integratedEquations, Scan\[Rho]MaxDependenceLPA, 
	rhoMaxPoints, baseGrid, baseEpsilon, arguments, ScanIntegralErrorLPA, 
	timeElapsed, exponentsLPA, exponentsLPANested, errorsLPA, logerrorsLPA, 
	SubtractReferenceValue, LogSubtractReferenceValue, 
	fitdata, fits, averageFit, exponent, constant, 
	maxval, minval, xmargin, ymargin, plotRange, 
	GetLabel, standardPrecisionPoints, plotLPA, filename, 
	potentialMinimum, referenceValuesNested, referenceValues, 
	subdirectory, directoryName, regulator,
	startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking error of grid compactification at the LPA level for N=" <> ToString[NValue] <> "."];
	
	subdirectory = "Rho Max Dependence/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	potentialMinimum = FindPotentialMinimum[NValue, "LPA", "litim", True];
	
	SelectRegulator["litim"];
	SetFixedFloatPrecision[doublePrecision];
	integratedEquations = flowEquationsLPA;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], True];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], True];
		
	Scan\[Rho]MaxDependenceLPA[rhoMax_, nPoints_, nPrecision_:doublePrecision]:=Block[{
		\[Rho]Max = rhoMax, gridSize = nPoints, 
		$MaxPrecision, $MinPrecision, $MaxExtraPrecision, 
		integratedLoop, fpSolution, stabilityMatrix, 
		eigenvalues, integralError, parameters, scalingTerms, discreteEquations}, 

		SetFixedFloatPrecision[nPrecision];
		
		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];

		fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessLPA];
		stabilityMatrix = StabilityMatrix[flowEquationsLPA, fpSolution, True];

		eigenvalues = LeadingEigenvalues[stabilityMatrix, 3];
		RestoreArbitraryPrecision[];
		eigenvalues = RecastFloatingPoint[eigenvalues];

		Return[<|{e1, \[Rho]Max / potentialMinimum} -> eigenvalues[[1]], 
				{e2, \[Rho]Max / potentialMinimum} -> eigenvalues[[2]], 
				{e3, \[Rho]Max / potentialMinimum} -> eigenvalues[[3]]|>];
	];
	
	rhoMaxPoints = (3/2 + 1/2 * Range[0, 13]);
	arguments = 
		Sort[Map[{potentialMinimum #, Ceiling[# / gridSpacing], doublePrecision}&, rhoMaxPoints], #1[[2]] > #2[[2]]&];
	
	referenceValues = Scan\[Rho]MaxDependenceLPA[optimal\[Rho]Max potentialMinimum, 
		Ceiling[optimal\[Rho]Max / gridSpacing], longDoublePrecision];
	
	exponentsLPA = ParallelMap[Scan\[Rho]MaxDependenceLPA[#[[1]], #[[2]], #[[3]]]&, arguments, Method -> "FinestGrained"];
	exponentsLPANested = KeySort[NestAssociation[Association[exponentsLPA]]];
	
	referenceValuesNested = KeyMap[#[[1]]&, referenceValues];
	errorsLPA = Map[Normal2, Abs[exponentsLPANested-referenceValuesNested]];
	
	xticks = Map[{#, If[Mod[#, 2] == 0, Round[#], ""]}&, rhoMaxPoints];
	yticks = Map[{10^-#, If[Mod[#, 3] == 1, Superscript[10, -#], ""]}&, Range[1, 15]];
	ticks = {xticks, yticks};
	gridLines = {{2, 3, 4, 5, 6, 7, 8}, {10^-13, 10^-10, 10^-7, 10^-4, 10^-1}};
	
	plotLPA = MakeScatterLogPlot[Values[errorsLPA], 
			Map[Subscript[\[CapitalDelta], #]&, Keys[errorsLPA]], Bottom, 
			{AxesLabel -> {Subscript[OverTilde[\[Rho]], "Max"]/Subscript[OverTilde[\[Rho]], 0], ""}, 
			 Ticks -> ticks, GridLines -> gridLines}];

	filename = "LPA N=" <> ToString[NValue];
	Export[directoryName <> filename <> ".wdx", exponentsLPA];
	Export[directoryName <> filename <> ".pdf", plotLPA];
		
	If[NValue == 2, 
		CopyFile[directoryName <> filename <> ".pdf", articleFiguresDirectory <> "compactification LPA.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
]


Map[Track\[Rho]MaxErrorDependenceLPA, {1, 2, 3}];


(* ::Subsubsection::Closed:: *)
(*DE2*)


Track\[Rho]MaxErrorDependenceDE2[nn_]:=Block[{NValue = nn, dimension = 3, alpha = 1, derivativePointsCount = 7, 
	gridSpacing = 1/20, optimal\[Rho]Max = 35/10, 
	integratedEquations, Scan\[Rho]MaxDependenceDE2, 
	rhoMaxPoints, baseGrid, baseEpsilon, arguments, ScanIntegralErrorDE2, 
	timeElapsed, exponentsDE2, exponentsDE2Nested, errorsDE2, logerrorsDE2, 
	SubtractReferenceValue, LogSubtractReferenceValue, 
	fitdata, fits, averageFit, exponent, constant, 
	maxval, minval, xmargin, ymargin, plotRange, 
	GetLabel, standardPrecisionPoints, plotDE2, filename, referenceValues, 
	potentialMinimum, referenceValuesNested, 
	subdirectory, directoryName, regulator,
	startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking error of grid compactification at the DE2 level for N=" <> ToString[NValue] <> "."];
		
	subdirectory = "Rho Max Dependence/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	potentialMinimum = FindPotentialMinimum[NValue, "DE2", "exponential", False];
	
	SelectRegulator["exponential"];
	SetFixedFloatPrecision[doublePrecision];
	SetGLIntegralParameters[GLIntegralStandardParameters];
	
	integratedEquations = flowEquationsDE2;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], True];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], True];
	
	Scan\[Rho]MaxDependenceDE2[rhoMax_, nPoints_, nPrecision_:doublePrecision] := Block[{
		\[Rho]Max = rhoMax, gridSize = nPoints, 
		integratedLoop, fpSolution, stabilityMatrix, 
		eigenvalues, integralError, parameters, scalingTerms, discreteEquations, 
		$MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
		
		SetFixedFloatPrecision[nPrecision];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];

		fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessDE2];		
		stabilityMatrix = StabilityMatrix[flowEquationsDE2, fpSolution, True];

		eigenvalues = LeadingEigenvalues[stabilityMatrix];
		RestoreArbitraryPrecision[];
		eigenvalues = RecastFloatingPoint[eigenvalues];

		Return[<|{e1, \[Rho]Max / potentialMinimum} -> eigenvalues[[1]], 
				{e2, \[Rho]Max / potentialMinimum} -> eigenvalues[[2]], 
				{e3, \[Rho]Max / potentialMinimum} -> eigenvalues[[3]], 
				{\[Eta], \[Rho]Max / potentialMinimum} -> (\[Eta] /. fpSolution)|>];
	];
	
	rhoMaxPoints = (3/2 + 1/2 Range[0, 13]);
	arguments = 
		Sort[Map[{potentialMinimum #, Ceiling[# / gridSpacing]}&, rhoMaxPoints], #1[[2]] > #2[[2]]&];
		
	referenceValues = Scan\[Rho]MaxDependenceDE2[optimal\[Rho]Max potentialMinimum, 
		Ceiling[optimal\[Rho]Max / gridSpacing], longDoublePrecision];
	
	exponentsDE2 = ParallelMap[Scan\[Rho]MaxDependenceDE2[#[[1]], #[[2]]]&, arguments, Method -> "FinestGrained"];
	exponentsDE2Nested = NestAssociation[Association[exponentsDE2]];
	
	referenceValuesNested = KeyMap[#[[1]]&, referenceValues];
	errorsDE2 = Map[Normal2, Abs[exponentsDE2Nested-referenceValuesNested]];
		
	xticks = Map[{#, If[Mod[#, 2] == 0, Round[#], ""]}&, rhoMaxPoints];
	yticks = Map[{10^-#, If[Mod[#, 3] == 1, Superscript[10, -#], ""]}&, Range[1, 15]];
	ticks = {xticks, yticks};
	gridLines = {{2, 3, 4, 5, 6, 7, 8}, {10^-13, 10^-10, 10^-7, 10^-4, 10^-1}};
	
	plotDE2 = MakeScatterLogPlot[Values[errorsDE2], 
			Map[Subscript[\[CapitalDelta], #]&, Keys[errorsDE2]], Bottom, 
			{AxesLabel -> {Subscript[OverTilde[\[Rho]], "Max"]/Subscript[OverTilde[\[Rho]], 0], ""}, 
			 Ticks -> ticks, GridLines -> gridLines}];
	
	filename = "DE2 N=" <> ToString[NValue];
	Export[directoryName <> filename <> ".wdx", exponentsDE2];
	Export[directoryName <> filename <> ".pdf", plotDE2];
	
	If[NValue == 2, 
		CopyFile[directoryName <> filename <> ".pdf", articleFiguresDirectory <> "compactification DE2.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
]


Map[Track\[Rho]MaxErrorDependenceDE2, {1, 2, 3}];


(* ::Subsection:: *)
(*Propagation of the discretization error*)


(* ::Subsubsection:: *)
(*LPA*)


TrackDiscretizationErrorLPA[nn_, \[Rho]MaxTo\[Rho]0_:35/10] := Block[{NValue = nn, dimension = 3, alpha = 1, 
		minimalSpacing = 1/5, referenceSpacing = 1/40, \[Rho]Max, 
		integratedEquations, GetDiscretizationErrorLPA, minimalGrid = 10, referenceGrid, 
		gridSizes, derivativePoints, standardPrecision, referencePrecision, 
		referenceArguments, referenceValues, arguments, potentialMinimum, 
		timeElapsed, exponents, nearExactValues, approximateValues, errors, PlotOneExponent, 
		maxval, minval, xmargin, ymargin, plotRange, dataFileName, low\[Rho]Max = (\[Rho]MaxTo\[Rho]0 < 3), 
		subdirectory, directoryName, regulator,
		startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking discretization-error propagation at the LPA level for N=" <> ToString[NValue] <> "."];
		
	If[!low\[Rho]Max,
		subdirectory = "Discretization Precision/", 
		subdirectory = "Discretization Precision Low Rho Max/";
	];
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	referenceGrid = Floor[\[Rho]MaxTo\[Rho]0 / referenceSpacing];
	potentialMinimum = FindPotentialMinimum[NValue, "LPA", "litim", True];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;
	
	SelectRegulator["litim"];
	integratedEquations = flowEquationsLPA;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], True];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], True];
	
	GetDiscretizationErrorLPA[gSize_, derivativePoints_, nPrecision_]:=Block[{parameters, scalingTerms, discreteEquations, 
		integratedLoop, fpSolution, stabilityMatrix, eigenvalues, 
		gridSize=gSize, derivativePointsCount=derivativePoints, 
		$MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
		SetFixedFloatPrecision[nPrecision];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];
		fpSolution = FindFixedPoint[scalingTerms+discreteEquations, guessLPA];
		stabilityMatrix = StabilityMatrix[flowEquationsLPA, fpSolution, True];

		eigenvalues = LeadingEigenvalues[stabilityMatrix, minimalGrid];
		RestoreArbitraryPrecision[];
		Return[<|{e1, derivativePoints, \[Rho]MaxTo\[Rho]0/gridSize, nPrecision} -> RecastFloatingPoint[eigenvalues[[1]]], 
				{e2, derivativePoints, \[Rho]MaxTo\[Rho]0/gridSize, nPrecision} -> RecastFloatingPoint[eigenvalues[[2]]], 
				{e3, derivativePoints, \[Rho]MaxTo\[Rho]0/gridSize, nPrecision} -> RecastFloatingPoint[eigenvalues[[3]]]|>];
	];
		
	gridSizes = Floor[\[Rho]MaxTo\[Rho]0/minimalSpacing * (5/4)^Range[0, 18]];
	gridSizes = Select[gridSizes, # >= minimalGrid&];
	derivativePoints = {3, 5, 7};
	standardPrecision = doublePrecision;
	referencePrecision = longDoublePrecision;
	referenceArguments = {{referenceGrid, 9, referencePrecision}};
	arguments = DeleteDuplicates[Join[referenceArguments, 
					Tuples[{gridSizes, derivativePoints, {standardPrecision}}]]];
	arguments = Reverse[Sort[arguments]];
	
	exponents = ParallelMap[GetDiscretizationErrorLPA[#[[1]], #[[2]], #[[3]]]&, arguments, Method -> "FinestGrained"];
	
	referenceValues = KeySort[KeyMap[#[[1;;2]]&, KeySelect[Association[exponents], #[[2]] == 9&]]];
	approximateValues = KeySort[OneLevelNest[KeyMap[#[[1;;3]]&, KeySelect[Association[exponents], #[[2]] != 9&]]]];
	errors = Abs[NestAssociation[approximateValues] - KeyMap[#[[1]]&, referenceValues]];
	
	maxval = Max[Map[Max, errors]];
	minval = Min[Map[Min, errors]];
	
	xmargin = (Max[gridSizes] / Min[gridSizes])^0.05;
	ymargin = (maxval / minval)^0.05;
	plotRange = {{\[Rho]MaxTo\[Rho]0 / Max[gridSizes] / xmargin, \[Rho]MaxTo\[Rho]0 / Min[gridSizes] * xmargin}, 
		{minval / ymargin, maxval * ymargin}};
	
	
	dataFileName = "LPA N=" <> ToString[NValue];
	Export[directoryName <> dataFileName <> ".wdx", exponents];
	
	PlotOneExponent[exponent_] := Block[{plotData = KeySort[errors[exponent]], fits, powerLaw, plot}, 
		plotData = Map[KeySort, plotData];

		If[low\[Rho]Max, 
			fits = Map[Normal[NonlinearModelFit[Log[Normal2[#1][[-6;;-1]]], A x + B, {A, B}, x]]&, Values[plotData]], 
			fits = MapThread[Normal[NonlinearModelFit[Log[Normal2[#1][[-6;;-1]]], (#2-1) x + B, B, x]]&, {Values[plotData], Keys[plotData]}];
		];
		
		xticks = Sort[Map[{#[[1]] 10^-#[[2]], If[#[[1]] == 1, Superscript[10, -#[[2]]], ""]} &, 
			Tuples[{Range[1, 9], Range[1, 3]}]]];
		yticks = Sort[Map[{10^(-2#), If[Mod[#, 2]==1, Superscript[10, -2#], ""]}&, Range[1, 7]]];
		ticks = {xticks, yticks};
		
		gridLines = {xticks[[1;;-1;;9, 1]], yticks[[1;;-1;;2]]};
		
		powerLaw = Map[Superscript[Subscript[h, \[Rho]], NumberForm[D[#, x], {2, 1}]]&, fits];
	
		plot = Show[
			MakeScatterLogLogPlot[Values[plotData], Map[ToString[#] <> "-points"&, Keys[plotData]], Bottom, 
				{AxesLabel -> {h\[Rho] / \[Rho]0, Subscript[\[CapitalDelta], exponent]}, PlotRange -> plotRange, 
				Ticks -> ticks, GridLines -> gridLines}, {LegendLayout -> {"Row", 2}}], 
			MakeFitLogLogPlot[Exp[fits] /. x -> Log[\[CapitalDelta]], {\[CapitalDelta], plotRange[[1, 1]], plotRange[[1, 2]]}, 
				powerLaw, Bottom, {}, {LegendLayout -> {"Row", 2}}]];
	
		Export[directoryName <> dataFileName <> " " <> ToString[exponent[[1]]] <> ToString[exponent[[2]]] <> ".pdf", plot];
	];
	Map[PlotOneExponent, {e1, e2, e3}];
	
	If[NValue == 2, 
		If[low\[Rho]Max, 
			CopyFile[directoryName <> dataFileName <> " e1.pdf", articleFiguresDirectory <> "discretization LPA e1 lowrho.pdf", OverwriteTarget -> True];
			CopyFile[directoryName <> dataFileName <> " e2.pdf", articleFiguresDirectory <> "discretization LPA e2 lowrho.pdf", OverwriteTarget -> True], 
		
			CopyFile[directoryName <> dataFileName <> " e1.pdf", articleFiguresDirectory <> "discretization LPA e1.pdf", OverwriteTarget -> True];
			CopyFile[directoryName <> dataFileName <> " e2.pdf", articleFiguresDirectory <> "discretization LPA e2.pdf", OverwriteTarget -> True];
		];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
];


Map[TrackDiscretizationErrorLPA, {1, 2, 3}];


Map[TrackDiscretizationErrorLPA[#, 2]&, {1, 2, 3}];


(* ::Subsubsection::Closed:: *)
(*DE2*)


TrackDiscretizationErrorDE2[nn_, \[Rho]MaxTo\[Rho]0_:35/10] := Block[{NValue = nn, dimension = 3, \[Rho]Max, alpha = 1, 
		minimalSpacing = 1/5, referenceSpacing = 1/30, minimalGrid = 10, 
		GetDiscretizationErrorDE2, gridSizes, derivativePoints, 
		standardPrecision, referenceGrid, referencePrecision, referenceArguments, arguments, 
		timeElapsed, exponents, referenceValues, approximateValues, 
		errors, flattenedErrors, minarg, maxarg, xmargin, 
		maxval, minval, ymargin, plotRange, dataFileName, PlotOneExponent, 
		potentialMinimum, low\[Rho]Max = (\[Rho]MaxTo\[Rho]0 < 3), 
		subdirectory, directoryName, regulator,
		startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking discretization-error propagation at the DE2 level for N=" <> ToString[NValue] <> "."];
		
	If[low\[Rho]Max, 
		subdirectory = "Discretization Precision Low Rho Max/", 	
		subdirectory = "Discretization Precision/";
	];
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	referenceGrid = Floor[\[Rho]MaxTo\[Rho]0 / referenceSpacing];
	potentialMinimum = FindPotentialMinimum[NValue, "DE2", "exponential", False];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;
	
	SelectRegulator["exponential"];
	
	GetDiscretizationErrorDE2[gSize_, derivativePointsC_, nPrecision_] := Block[{parameters, scalingTerms, discreteEquations, 
		GLIntegralIntervalHalfSpan, GLIntegralIntervalMidpoint, GLIntegralEvaluationPoints, 
		GLIntegralReplacementRules, GLIntegralWeights, 
		integratedEquations, fpSolution, stabilityMatrix, eigenvalues, gridEpsilon, 
		gridSize=gSize, derivativePointsCount=derivativePointsC, $MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
		
		SetGLIntegralParameters[GLIntegralStandardParameters];
		SetFixedFloatPrecision[nPrecision];
		
		integratedEquations = flowEquationsDE2;
		integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], True];
		integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], True];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];
		fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessDE2];
		gridEpsilon = potentialMinimum / gridSize;
		
		stabilityMatrix = StabilityMatrix[flowEquationsDE2, fpSolution];
		eigenvalues = LeadingEigenvalues[stabilityMatrix];
		RestoreArbitraryPrecision[];
		eigenvalues = RecastFloatingPoint[eigenvalues];
		
		Return[<|{e1, derivativePointsCount, \[Rho]MaxTo\[Rho]0 / gridSize, nPrecision} -> eigenvalues[[1]], 
				{e2, derivativePointsCount, \[Rho]MaxTo\[Rho]0 / gridSize, nPrecision} -> eigenvalues[[2]], 
				{e3, derivativePointsCount, \[Rho]MaxTo\[Rho]0 / gridSize, nPrecision} -> eigenvalues[[3]], 
				{\[Eta], derivativePointsCount, \[Rho]MaxTo\[Rho]0 / gridSize, nPrecision} -> RecastFloatingPoint[\[Eta] /. fpSolution]|>];
	];
	
	gridSizes = Floor[\[Rho]MaxTo\[Rho]0/minimalSpacing * (5/4)^Range[0, 14]];
	gridSizes = Select[gridSizes, # >= minimalGrid&];
	derivativePoints = {3, 5, 7};
	standardPrecision = doublePrecision;
	referencePrecision = longDoublePrecision;
	referenceArguments = {{referenceGrid, 9, referencePrecision}};
	arguments = DeleteDuplicates[Join[referenceArguments, 
					Tuples[{gridSizes, derivativePoints, {standardPrecision}}]]];
	arguments = Reverse[Sort[arguments]];
	
	exponents = Map[GetDiscretizationErrorDE2[#[[1]], #[[2]], #[[3]]]&, arguments];
	(* exponents = ParallelMap[GetDiscretizationErrorDE2[#[[1]], #[[2]], #[[3]]]&, arguments, Method->"FinestGrained"];*);
	
	referenceValues = KeySort[KeyMap[#[[1;;2]]&, KeySelect[Association[exponents], #[[2]] == 9&]]];
	approximateValues = KeySort[OneLevelNest[KeyMap[#[[1;;3]]&, KeySelect[Association[exponents], #[[2]] != 9&]]]];
	errors = Abs[NestAssociation[approximateValues]- KeyMap[#[[1]]&, referenceValues]];
	
	
	maxval = Max[errors];
	minval = Min[errors];
	flattenedErrors = UnnestAssociation[errors];
	minarg = Min[Map[#[[3]]&, Keys[flattenedErrors]]];
	maxarg = Max[Map[#[[3]]&, Keys[flattenedErrors]]];
	xmargin = (maxarg / minarg)^0.05;
	ymargin = (maxval / minval)^0.05;
	plotRange = {{minarg / xmargin, maxarg * xmargin}, {minval / ymargin, maxval * ymargin}};
	
	
	dataFileName = "DE2 N=" <> ToString[NValue];
	Export[directoryName <> dataFileName <> ".wdx", exponents];
	
	PlotOneExponent[exponent_] := Block[{plotData = KeySort[errors[exponent]], fits, powerLaw, plot, exponentLabel, theoreticalExponents}, 
		plotData = Map[KeySort, plotData];

		If[low\[Rho]Max, 
			fits = Map[Normal[NonlinearModelFit[Log[Normal2[#1][[-7;;-2]]], A x + B, {A, B}, x]]&, Values[plotData]], 

			theoreticalExponents = Keys[plotData] - 1;
			fits = MapThread[Normal[NonlinearModelFit[Log[Normal2[#1][[-7;;-2]]], #2 x + B, B, x]]&, {Values[plotData], theoreticalExponents}];
		];
		
		powerLaw = Map[Superscript[Subscript[h, \[Rho]], NumberForm[D[#, x], {2, 1}]]&, fits];
		
		
		xticks = Sort[Map[{#[[1]] 10^-#[[2]], If[#[[1]] == 1, Superscript[10, -#[[2]]], ""]} &, 
			Tuples[{Range[1, 9], Range[1, 3]}]]];
		yticks = Sort[Map[{10^(-2#), If[Mod[#, 2]==1, Superscript[10, -2#], ""]}&, Range[1, 7]]];
		ticks = {xticks, yticks};
		
		gridLines = {xticks[[1;;-1;;9, 1]], yticks[[1;;-1;;2]]};
	
		plot = Show[
			MakeScatterLogLogPlot[Values[plotData], Map[ToString[#] <> "-points"&, Keys[plotData]], Bottom, 
				{AxesLabel -> {h\[Rho] / \[Rho]0, Subscript[\[CapitalDelta], exponent]}, PlotRange -> plotRange, 
				Ticks -> ticks, GridLines -> gridLines}, 
				{LegendLayout -> {"Row", 2}}], 
			MakeFitLogLogPlot[Exp[fits] /. x -> Log[\[CapitalDelta]], {\[CapitalDelta], minarg, maxarg}, powerLaw, 
				Bottom, {}, {LegendLayout -> {"Row", 2}}]
			];		
		
		If[Length[exponent] > 0, 
			exponentLabel = ToString[exponent[[1]]] <> ToString[exponent[[2]]], 
			exponentLabel = ToString[exponent];
		];
	
		Export[directoryName <> dataFileName <> " " <> exponentLabel <> ".pdf", plot];
	];
	
	Map[PlotOneExponent, {e1, e2, e3, \[Eta]}];
		
	If[NValue == 2, 
		If[low\[Rho]Max, 
			CopyFile[directoryName <> dataFileName <> " \[Eta].pdf", articleFiguresDirectory <> "discretization DE2 eta lowrho.pdf", OverwriteTarget -> True];
			CopyFile[directoryName <> dataFileName <> " e1.pdf", articleFiguresDirectory <> "discretization DE2 e1 lowrho.pdf", OverwriteTarget -> True];
			CopyFile[directoryName <> dataFileName <> " e2.pdf", articleFiguresDirectory <> "discretization DE2 e2 lowrho.pdf", OverwriteTarget -> True], 
		
			CopyFile[directoryName <> dataFileName <> " \[Eta].pdf", articleFiguresDirectory <> "discretization DE2 eta.pdf", OverwriteTarget -> True];
			CopyFile[directoryName <> dataFileName <> " e1.pdf", articleFiguresDirectory <> "discretization DE2 e1.pdf", OverwriteTarget -> True];
			CopyFile[directoryName <> dataFileName <> " e2.pdf", articleFiguresDirectory <> "discretization DE2 e2.pdf", OverwriteTarget -> True];
		];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
];


ParallelMap[TrackDiscretizationErrorDE2[#[[1]], #[[2]]]&, 
	{{1, 3.5}, {2, 3.5}, {3, 3.5}, 
	 {1, 2}, {2, 2}, {3, 2}}];


(* ::Subsection::Closed:: *)
(*Comparison of the discretization error between subsequent eigenvalues*)


(* ::Subsubsection::Closed:: *)
(*LPA*)


TrackEigenvalueErrorLPA[nn_] := Block[{NValue = nn, dimension = 3, \[Rho]Max = 15, alpha = 1, derivativePointsCount = 7, 
	\[Rho]0To\[Rho]Max = 35/10, minGridSpacing = 1/10, minGridSize, integratedEquations, 
	potentialMinimum, gridSizes, standardPrecision, referencePrecision, 
	referenceArguments, arguments, timeElapsed, exponentsLPA, 
	referenceValuesLPA, approximateValuesLPA, errors, plot, dataFileName, 
	maxval, minval, xmargin, ymargin, minx, maxx, plotRange, errorThreshold=0.01, 
	subdirectory, directoryName, regulator, GetEigenvaluesLPA,
	startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Comparing errors of subsequent eigenvalues at the LPA level for N=" <> ToString[NValue] <> "."];
	
	subdirectory = "Eigenvalue Precision/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
		
	minGridSize = \[Rho]0To\[Rho]Max / minGridSpacing;
	potentialMinimum = FindPotentialMinimum[NValue, "LPA", "litim", True];
	\[Rho]Max = \[Rho]0To\[Rho]Max*potentialMinimum;
	
	SelectRegulator["litim"];
	
	integratedEquations = flowEquationsLPA;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], True];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], True];
	
	GetEigenvaluesLPA[gSize_, nPrecision_]:=Block[{parameters, scalingTerms, discreteEquations, 
	integratedLoop, fpSolution, stabilityMatrix, eigenvalues, 
		gridSize=gSize, $MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
		SetFixedFloatPrecision[nPrecision];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];
		fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessLPA];
		stabilityMatrix = StabilityMatrix[flowEquationsLPA, fpSolution, True];

		eigenvalues = LeadingEigenvalues[stabilityMatrix, 20];

		RestoreArbitraryPrecision[];
		Return[<|{gridSize, nPrecision} -> RecastFloatingPoint[eigenvalues]|>];
	];
		
	gridSizes = minGridSize {1, 2, 4};
	standardPrecision = doublePrecision;
	referencePrecision = longDoublePrecision;
	referenceArguments = {{6 minGridSize, referencePrecision}};
	arguments = DeleteDuplicates[Join[referenceArguments, 
					Tuples[{gridSizes, {standardPrecision}}]]];
	arguments = Reverse[Sort[arguments]];
	
	exponentsLPA = ParallelMap[GetEigenvaluesLPA[#[[1]], #[[2]]]&, arguments, Method -> "FinestGrained"];
	
	referenceValuesLPA = Association[exponentsLPA][referenceArguments[[1]]];
	approximateValuesLPA = KeySort[KeyMap[#[[1]]&, KeySelect[Association[exponentsLPA], # != referenceArguments[[1]]&]]];
	
	errors = Map[Select[Transpose[{Range[Length[#]], Abs[# - referenceValuesLPA]}], #[[2]] < errorThreshold&]&, approximateValuesLPA];

	maxval = Max[errors[[All, All, 2]]];
	minval = Min[errors[[All, All, 2]]];
	xmargin = 0.75;
	ymargin = (maxval / minval)^0.075;
	minx = Min[errors[[All, All, 1]]];
	maxx = Max[errors[[All, All, 1]]];
	plotRange = {{minx - xmargin, maxx + xmargin}, {minval / ymargin, maxval * ymargin}};
	
	xticks = Map[{#, If[Mod[#, 5]==0, #, ""]}&, Range[1, 20]];
	yticks = Map[{10^-#, If[Mod[#, 3]==1, Superscript[10,-#], ""]}&, Range[1, 15]];
	ticks = {xticks, yticks};
	gridLines = {5 Range[4], 10^(-1-3 Range[5])};
		
	plot = MakeScatterLogPlot[Values[errors], 
			Map[h\[Rho] <> "=" <> ToString[\[Rho]0 / Rationalize[# / \[Rho]0To\[Rho]Max], StandardForm]&, Keys[errors]], Right, 
			{PlotRange -> plotRange, AxesLabel -> {i, Subscript[\[CapitalDelta], Subscript[e, i]]},
			 Ticks -> ticks, GridLines -> gridLines}, LegendLayout -> "Column"];
	
	dataFileName = "LPA N=" <> ToString[NValue];
	Export[directoryName <> dataFileName <> ".wdx", exponentsLPA];
	Export[directoryName <> dataFileName <> ".pdf", plot];
	
	If[NValue == 2, 
		CopyFile[directoryName <> dataFileName <> ".pdf", articleFiguresDirectory <> "eigenvalue LPA.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
]


Map[TrackEigenvalueErrorLPA, {1, 2, 3}];


(* ::Subsubsection::Closed:: *)
(*DE2*)


TrackEigenvalueErrorDE2[nn_] := Block[{NValue = nn, dimension = 3, \[Rho]Max, alpha = 1, derivativePointsCount = 7, 
	\[Rho]0To\[Rho]Max = 35/10, minGridSpacing = 1/10, minGridSize, 
	integratedEquations, potentialMinimum, gridSizes, standardPrecision, referencePrecision, 
	referenceArguments, arguments, timeElapsed, exponentsDE2, 
	referenceValuesDE2, approximateValuesDE2, approximateValuesDE2Selected, referenceValuesDE2Selected, 
	errorsDE2, plotValues, plotStyle, plotMarkers, plot, plots, 
	dataFileName, maxval, minval, xmargin, ymargin, minx, maxx, plotRange, errorThreshold=0.01, 
	subdirectory, directoryName, regulator, GetEigenvaluesDE2,
	startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Comparing errors of subsequent eigenvalues at the DE2 level for N=" <> ToString[NValue] <> "."];
	
	subdirectory = "Eigenvalue Precision/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	minGridSize = \[Rho]0To\[Rho]Max / minGridSpacing;
	potentialMinimum = FindPotentialMinimum[NValue, "DE2", "exponential", False];
	\[Rho]Max = \[Rho]0To\[Rho]Max*potentialMinimum;
	
	SelectRegulator["exponential"];
	SetGLIntegralParameters[GLIntegralStandardParameters];
	
	integratedEquations = flowEquationsDE2;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]]];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]]];
	
	GetEigenvaluesDE2[ gSize_, nPrecision_]:=Block[{parameters, scalingTerms, discreteEquations, 
		integratedLoop, fpSolution, stabilityMatrix, eigenvalues, 
		gridSize=gSize, $MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
		SetFixedFloatPrecision[nPrecision];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];
		fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessDE2];
		stabilityMatrix = StabilityMatrix[flowEquationsDE2, fpSolution, True];

		eigenvalues = LeadingEigenvalues[stabilityMatrix, 30];
		RestoreArbitraryPrecision[];
		Return[<|{gridSize, nPrecision} -> RecastFloatingPoint[eigenvalues]|>];
	];
		
	gridSizes = minGridSize {1, 2};
	standardPrecision = doublePrecision;
	referencePrecision = longDoublePrecision;
	referenceArguments = {{3 minGridSize, referencePrecision}};
	arguments = DeleteDuplicates[Join[referenceArguments, 
					Tuples[{gridSizes, {standardPrecision}}]]];
	arguments = Reverse[Sort[arguments]];
	
	exponentsDE2 = ParallelMap[GetEigenvaluesDE2[#[[1]], #[[2]]]&, arguments, Method -> "FinestGrained"];
	
	referenceValuesDE2 = Association[exponentsDE2][referenceArguments[[1]]];
	approximateValuesDE2 = KeySort[KeyMap[#[[1]]&, KeySelect[Association[exponentsDE2], # != referenceArguments[[1]]&]]];
	approximateValuesDE2Selected = Map[Select[#, Im[#] >= 0&]&, approximateValuesDE2];
	referenceValuesDE2Selected = Select[referenceValuesDE2, Im[#]>=0&];
	errorsDE2 = Map[Transpose[{Range[Length[#]], # - referenceValuesDE2Selected}]&, approximateValuesDE2Selected];

	plotValues = Map[{Abs[Select[#, Im[#[[2]]] == 0 && Abs[#[[2]]] < errorThreshold &]], 
					 Abs[Select[#, Im[#[[2]]] != 0 && Abs[#[[2]]] < errorThreshold &]]}&, errorsDE2];
	plotStyle = Map[{#, #}&, DefaultPlotColors][[1;;Length[plotValues]]];
	plotMarkers = MapThread[{{#1, DefaultMarkerSize}, {#2, DefaultMarkerSize}}&, {FilledMarkers, OpenMarkers}][[1;;Length[plotValues]]];
	
	maxval = Max[plotValues[[All, All, All, 2]]];
	minval = Min[plotValues[[All, All, All, 2]]];
	xmargin = 0.75;
	ymargin = (maxval / minval)^0.075;
	minx = Min[plotValues[[All, All, All, 1]]];
	maxx = Max[plotValues[[All, All, All, 1]]];
	plotRange = {{minx - xmargin, maxx + xmargin}, {minval / ymargin, maxval * ymargin}};
		
	xticks = Map[{#, If[Mod[#, 5]==0, #, ""]}&, Range[1, 20]];
	yticks = Map[{10^-#, If[Mod[#, 3]==1, Superscript[10,-#], ""]}&, Range[1, 15]];
	ticks = {xticks, yticks};
	gridLines = {5 Range[4], 10^(-1-3 Range[5])};
	
	plots = MapThread[ListLogPlot[#1, Evaluate[CommonSettings], PlotMarkers -> #3, PlotStyle -> #4, 
		PlotRange -> plotRange, AxesLabel -> {i, Subscript[\[CapitalDelta], Subscript[e, i]]}, 
		PlotLegends -> PointLegend[{"Real", "Complex"}, 
			LegendLabel -> h\[Rho] <> "=" <> ToString[\[Rho]0 / Rationalize[#2 / \[Rho]0To\[Rho]Max], StandardForm]],
		Ticks -> ticks, GridLines -> gridLines]&, 
		{Values[plotValues], Keys[plotValues], plotMarkers, plotStyle}];
	
	plot = Show[plots];
	
	dataFileName = "DE2 N=" <> ToString[NValue];
	Export[directoryName <> dataFileName <> ".wdx", exponentsDE2];
	Export[directoryName <> dataFileName <> ".pdf", plot];
		
	If[NValue == 2, 
		CopyFile[directoryName <> dataFileName <> ".pdf", articleFiguresDirectory <> "eigenvalue DE2.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
	
]


Map[TrackEigenvalueErrorDE2, {1, 2, 3}];


(* ::Subsection::Closed:: *)
(*Propagation of the loop-integral error*)


(* ::Subsubsection::Closed:: *)
(*LPA*)


TrackIntegralPrecisionErrorLPA[nn_]:=Block[{NValue = nn, gridSize, \[Rho]Max, dimension = 3, alpha = 1, derivativePointsCount = 7, 
	\[Rho]MaxTo\[Rho]0 = 35/10, gridSpacing = 1/20, integralErrorThreshold = 10^-13, 
	referenceParameters, referenceScalingTerms, referenceIntegratedEquations, referenceDiscreteEquations, 
	referenceFPSolution, referenceStabilityMatrix, referenceEigenvalues, 
	parametersReference, scalingTermsReference, discreteEquationsReference, 
	fpSolutionReference, potentialMinimum, 
	ScanIntegralErrorLPA, timeElapsed, scanPoints, exponentsLPA, exponentsLPANested, 
	SubtractReferenceValue, LogSubtractReferenceValue, errorsLPA, logerrorsLPA, 
	fits, averageFit, plotLPA, filename, subdirectory, directoryName, regulator,
	startTime, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking integral-error propagation at the LPA level for N=" <> ToString[NValue] <> "."];
	
	gridSize = Ceiling[\[Rho]MaxTo\[Rho]0 / gridSpacing];
	
	subdirectory = "Integral Precision/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
	CreateDirectory[directoryName]];
	
	potentialMinimum = FindPotentialMinimum[NValue, "LPA", "exponential", False];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;
	
	SelectRegulator["exponential"];
	SetFixedFloatPrecision[longDoublePrecision];
	
	SetGLIntegralParameters[GLIntegralReferenceParameters];
	referenceIntegratedEquations = flowEquationsLPA;
	referenceIntegratedEquations[[All, 3]] = PerformMomentumIntegrals[referenceIntegratedEquations[[All, 3]]];
	referenceIntegratedEquations[[All, 4]] = PerformMomentumIntegrals[referenceIntegratedEquations[[All, 4]]];
	{referenceParameters, referenceScalingTerms, referenceDiscreteEquations} = ReplaceConstants[DiscretizeFlowEquations[referenceIntegratedEquations, True]];
	
	referenceFPSolution = FindFixedPoint[referenceScalingTerms + referenceDiscreteEquations, guessLPA];
	referenceStabilityMatrix = StabilityMatrix[flowEquationsLPA, referenceFPSolution, False];

	referenceEigenvalues = LeadingEigenvalues[referenceStabilityMatrix];
	RestoreArbitraryPrecision[];
	referenceEigenvalues = RecastFloatingPoint[referenceEigenvalues];
	referenceEigenvalues = <|e1 -> referenceEigenvalues[[1]], e2 -> referenceEigenvalues[[2]], e3 -> referenceEigenvalues[[3]]|>;
	
	ScanIntegralErrorLPA[integralPointCount_] := Block[{integratedLoop, fpSolution, stabilityMatrix, GLIntegralPointCount, 
		eigenvalues, integralError, parameters, scalingTerms, discreteEquations, 
		GLIntegralEvaluationPoints, GLIntegralSmallParameter, GLIntegralReplacementRules, 
		$MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 

		SetGLIntegralParameters[0, GLIntegralReferenceParameters[[2]], integralPointCount];
		SetFixedFloatPrecision[doublePrecision];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[flowEquationsLPA, True];
		integratedLoop = PerformMomentumIntegrals[discreteEquations];

		fpSolution = FindFixedPoint[scalingTerms+integratedLoop, referenceFPSolution];
		integralError = Max[Abs[(referenceDiscreteEquations - ReplaceConstants[integratedLoop]) /. fpSolution]];
		stabilityMatrix = StabilityMatrix[flowEquationsLPA, fpSolution, False];

		eigenvalues = LeadingEigenvalues[stabilityMatrix];
		RestoreArbitraryPrecision[];
		eigenvalues = RecastFloatingPoint[eigenvalues];

		Return[<|{integralError, e1} -> eigenvalues[[1]], 
				{integralError, e2} -> eigenvalues[[2]], 
				{integralError, e3} -> eigenvalues[[3]]|>];
	];
	
	scanPoints = 5 Reverse[Range[1, 8]];
	exponentsLPA = ParallelMap[ScanIntegralErrorLPA, scanPoints, Method -> "FinestGrained"];
	exponentsLPANested = TransposeAssociation[KeySort[NestAssociation[Association[exponentsLPA]]]];
	
	SubtractReferenceValue[assoc_] := Map[{#, Abs[assoc[#] - referenceEigenvalues]}&, Keys[KeySort[assoc]]];
	LogSubtractReferenceValue[assoc_] := Map[{Log[#], Log[Abs[assoc[#] - referenceEigenvalues]]}&, Keys[KeySort[assoc]]];
	errorsLPA = Abs[exponentsLPANested - referenceEigenvalues];
	logerrorsLPA = Map[Log[Normal2[#]]&, errorsLPA];
	
	fits = Map[Normal[NonlinearModelFit[#[[-8;;-1]], x + b, b, x]]&, Values[logerrorsLPA]];
	averageFit = Simplify[Exp[Mean[fits]] /. x -> Log[\[CapitalDelta]]];
	
	ticks = {{{10^-12, ""}, {10^-10, Superscript[10, -10]}, {10^-8, ""}, 
		{10^-6, Superscript[10, -6]}, {10^-4, ""}, {10^-2, Superscript[10, -2]}}, 
		{{10^-10, Superscript[10, -10]}, {10^-6, Superscript[10, -6]}, {10^-2, Superscript[10, -2]}}}; 
	gridLines = ticks[[All, All, 1]];
	
	plotLPA = Show[
		MakeScatterLogLogPlot[Values[errorsLPA], 
			Map[Subscript[\[CapitalDelta], #]&, Keys[errorsLPA]], Bottom, 
			{AxesLabel -> {Subscript[\[CapitalDelta], "I"], ""}, Ticks -> ticks, GridLines -> gridLines}], 
		MakeFitLogLogPlot[averageFit, {\[CapitalDelta], 10^-16, 1}, {Superscript[Subscript[\[CapitalDelta], "I"], 1]}, Bottom, PlotRange -> All]
	];

	filename = "LPA N=" <> ToString[NValue];
	Export[directoryName <> filename <> ".wdx", exponentsLPA];
	Export[directoryName <> filename <> ".pdf", plotLPA];
	
	If[NValue == 2, 
		CopyFile[directoryName <> filename <> ".pdf", articleFiguresDirectory <> "integral LPA.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
]


Map[TrackIntegralPrecisionErrorLPA, {1, 2, 3}];


(* ::Subsubsection::Closed:: *)
(*DE2*)


TrackIntegralPrecisionErrorDE2[nn_]:=Block[{NValue = nn, gridSize, dimension = 3, \[Rho]Max, alpha = 1, 
	derivativePointsCount = 7, \[Rho]MaxTo\[Rho]0 = 35/10, gridSpacing = 1/20, integralErrorThreshold = 10^-13, 
	referenceIntegratedEquations, referenceParameters, referenceScalingTerms, referenceDiscreteEquations, 
	referenceFPSolution, referenceStabilityMatrix, referenceEigenvalues, referenceExponents, 
	ScanIntegralErrorDE2, timeElapsed, scanPoints, exponentsDE2, exponentsDE2Nested, 
	parametersReference, scalingTermsReference, discreteEquationsReference, fpSolutionReference, potentialMinimum, 
	errorsDE2, logerrorsDE2, fits, averageFit, 
	GetLabel, standardPrecisionPoints, plotDE2, filename, subdirectory, directoryName, regulator,
	startTime, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking integral-error propagation at the DE2 level for N=" <> ToString[NValue] <> "."];
		
	subdirectory = "Integral Precision/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
	CreateDirectory[directoryName]];
	
	potentialMinimum = FindPotentialMinimum[NValue, "DE2", "exponential", False];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;
	
	gridSize = \[Rho]MaxTo\[Rho]0 / gridSpacing;
	SelectRegulator["exponential"];
	SetFixedFloatPrecision[longDoublePrecision];
	SetGLIntegralParameters[GLIntegralReferenceParameters];

	referenceIntegratedEquations = flowEquationsDE2;
	referenceIntegratedEquations[[All, 3]] = PerformMomentumIntegrals[referenceIntegratedEquations[[All, 3]]];
	referenceIntegratedEquations[[All, 4]] = PerformMomentumIntegrals[referenceIntegratedEquations[[All, 4]]];
	{referenceParameters, referenceScalingTerms, referenceDiscreteEquations} = 
		ReplaceConstants[DiscretizeFlowEquations[referenceIntegratedEquations, True]];
	
	referenceFPSolution = FindFixedPoint[referenceScalingTerms + referenceDiscreteEquations, guessDE2];
	referenceStabilityMatrix = StabilityMatrix[flowEquationsDE2, referenceFPSolution, False];

	referenceEigenvalues = LeadingEigenvalues[referenceStabilityMatrix];
	RestoreArbitraryPrecision[];
	referenceEigenvalues = RecastFloatingPoint[referenceEigenvalues];
	
	referenceExponents = <|e1 -> referenceEigenvalues[[1]], 
		e2 -> referenceEigenvalues[[2]], 
		e3 -> referenceEigenvalues[[3]], 
		\[Eta] -> RecastFloatingPoint[\[Eta] /. referenceFPSolution]
		|>;
	
	ScanIntegralErrorDE2[integralPointCount_] := Block[{integratedLoop, fpSolution, stabilityMatrix, GLIntegralPointCount, 
		eigenvalues, integralError, parameters, scalingTerms, discreteEquations, 
		GLIntegralEvaluationPoints, GLIntegralSmallParameter, GLIntegralReplacementRules, 
		$MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 

		SetGLIntegralParameters[0, GLIntegralReferenceParameters[[2]], integralPointCount];
		SetFixedFloatPrecision[doublePrecision];

		{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[flowEquationsDE2, True];
		integratedLoop = PerformMomentumIntegrals[discreteEquations];

		fpSolution = FindFixedPoint[scalingTerms+integratedLoop, guessDE2];
		integralError = Max[Abs[(referenceDiscreteEquations - ReplaceConstants[integratedLoop]) /. fpSolution]];
		stabilityMatrix = StabilityMatrix[flowEquationsDE2, fpSolution, False];

		eigenvalues = LeadingEigenvalues[stabilityMatrix];
		RestoreArbitraryPrecision[];
		eigenvalues = RecastFloatingPoint[eigenvalues];

		Return[<|{integralError, e1} -> eigenvalues[[1]], 
				{integralError, e2} -> eigenvalues[[2]], 
				{integralError, e3} -> eigenvalues[[3]], 
				{integralError, \[Eta]} -> RecastFloatingPoint[\[Eta] /. fpSolution]|>];
	];
	
	scanPoints = 5 Reverse[Range[1, 8]];
	exponentsDE2 = ParallelMap[ScanIntegralErrorDE2, scanPoints, Method -> "FinestGrained"];
	exponentsDE2Nested = TransposeAssociation[KeySort[NestAssociation[Association[exponentsDE2]]]];
	
	errorsDE2 = Abs[exponentsDE2Nested - referenceExponents];
	logerrorsDE2 = Map[Log[Normal2[#]]&, errorsDE2];
	
	fits = Map[Normal[NonlinearModelFit[#[[-8;;-1]], x + b, b, x]]&, Values[logerrorsDE2]];
	averageFit = Simplify[Exp[Mean[fits]] /. x -> Log[\[CapitalDelta]]];
		
	ticks = {{{10^-12, ""}, {10^-10, Superscript[10, -10]}, {10^-8, ""}, 
		{10^-6, Superscript[10, -6]}, {10^-4, ""}, {10^-2, Superscript[10, -2]}}, 
		{{10^-10, Superscript[10, -10]}, {10^-6, Superscript[10, -6]}, {10^-2, Superscript[10, -2]}}}; 
	gridLines = ticks[[All, All, 1]];
	
	plotDE2 = Show[
		MakeScatterLogLogPlot[Values[errorsDE2], 
			Map[Subscript[\[CapitalDelta], #]&, Keys[errorsDE2]], Bottom, 
			{AxesLabel -> {Subscript[\[CapitalDelta], "I"], ""}, Ticks -> ticks, GridLines -> gridLines}], 
		MakeFitLogLogPlot[averageFit, {\[CapitalDelta], 10^-16, 1}, {Superscript[Subscript[\[CapitalDelta], "I"], 1]}, Bottom, PlotRange -> All]
	];

	filename = "DE2 N=" <> ToString[NValue];
	Export[directoryName <> filename <> ".wdx", exponentsDE2];
	Export[directoryName <> filename <> ".pdf", plotDE2];
	
	If[NValue == 2, 
		CopyFile[directoryName <> filename <> ".pdf", articleFiguresDirectory <> "integral DE2.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
]


Map[TrackIntegralPrecisionErrorDE2, {1, 2, 3}];


(* ::Subsection::Closed:: *)
(*Error of the finite-difference stability-matrix approximation*)


(* ::Subsubsection::Closed:: *)
(*LPA*)


TrackStabilityMatrixApproximationErrorLPA[nn_] := Block[{NValue = nn, derivativePointsCount = 7, gridSize, 
	dimension = 3, \[Rho]Max, alpha = 1, \[Rho]MaxTo\[Rho]0 = 35/10, gridSpacing = 1/20, 
	integratedEquations, parameters, scalingTerms, discreteEquations, fpSolution, 
	GetStabilityMatrixErrorLPA, epsilonValues, jacobianPointCount, arguments, 
	timeElapsed, eigenvaluesLPA, eigenvaluesLPANested, 
	SubtractExactValue, SubtractExactValueLogLog, 
	errorsLPA, logerrorsLPA, 
	flattenedErrors, minarg, maxarg, xmargin, 
	minval, maxval, ymargin, plotRange, 
	dataFileName, PlotOneExponent, 
	parametersReference, scalingTermsReference, discreteEquationsReference, 
	fpSolutionReference, potentialMinimum, FilterHighError, 
	subdirectory, directoryName, regulator,
	startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking stability-matrix-approximation-error propagation at the LPA level for N=" <> ToString[NValue] <> "."];
	
	subdirectory = "Stability Matrix Approximation Precision/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	gridSize = Ceiling[\[Rho]MaxTo\[Rho]0 / gridSpacing];
	potentialMinimum = FindPotentialMinimum[NValue, "LPA", "litim", True];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;

	SelectRegulator["litim"];
	SetFixedFloatPrecision[doublePrecision];
	integratedEquations = flowEquationsLPA;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]], True];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]], True];	
	{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];

	fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessLPA];
	
	RestoreArbitraryPrecision[];
	
	GetStabilityMatrixErrorLPA[epsilon_, nPointsStabilityMatrix_] := 
		Block[{stabilityMatrix, eigenvalues, $MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
			
			SetFixedFloatPrecision[doublePrecision];
			stabilityMatrix = StabilityMatrix[flowEquationsLPA, fpSolution, True, nPointsStabilityMatrix, epsilon];
			eigenvalues = LeadingEigenvalues[stabilityMatrix];
			RestoreArbitraryPrecision[];
			
			Return[<|{e1, nPointsStabilityMatrix, epsilon} -> RecastFloatingPoint[eigenvalues[[1]]], 
				{e2, nPointsStabilityMatrix, epsilon} -> RecastFloatingPoint[eigenvalues[[2]]], 
				{e3, nPointsStabilityMatrix, epsilon} -> RecastFloatingPoint[eigenvalues[[3]]]|>];
		];
	
	epsilonValues = N[Join[{0}, 10^(-2-Range[0, 24]/3)]];
	jacobianPointCount = {1, 2, 3};
	arguments = Tuples[{epsilonValues, jacobianPointCount}];
	
	eigenvaluesLPA = ParallelMap[GetStabilityMatrixErrorLPA[#[[1]], #[[2]]]&, arguments, Method -> "FinestGrained"];
	eigenvaluesLPANested = OneLevelNest[Association[eigenvaluesLPA]];
	
	SubtractExactValue[exponent_] := Association[Map[# -> Abs[eigenvaluesLPANested[exponent][#]-eigenvaluesLPANested[exponent][0.]]&, 
		Sort[Keys[eigenvaluesLPANested[exponent]]][[2;;-1]]]];
	SubtractExactValueLogLog[exponent_] := Association[Map[Log[#] -> Log[Abs[eigenvaluesLPANested[exponent][#]-eigenvaluesLPANested[exponent][0.]]]&, 
		Sort[Keys[eigenvaluesLPANested[exponent]]][[2;;-1]]]];

	errorsLPA = NestAssociation[Association[Map[# -> SubtractExactValue[#]&, Keys[eigenvaluesLPANested]]]];
	logerrorsLPA = NestAssociation[Association[Map[# -> SubtractExactValueLogLog[#]&, Keys[eigenvaluesLPANested]]]];
	
	FilterHighError[assoc_, threshold_:0.05] := Map[Select[#, #<threshold&]&, assoc];
	errorsLPA = Map[FilterHighError, errorsLPA];
	
	flattenedErrors = UnnestAssociation[errorsLPA];
	minarg = Min[Map[#[[3]]&, Keys[flattenedErrors]]];
	maxarg = Max[Map[#[[3]]&, Keys[flattenedErrors]]];
	xmargin = (maxarg / minarg)^0.05;
	maxval = Min[Max[Map[Max, errorsLPA]], 10^-1];
	minval = Min[Map[Min, errorsLPA]];
	ymargin = (maxval / minval)^0.05;
	plotRange = {{minarg / xmargin, maxarg * xmargin}, {minval / ymargin, maxval * ymargin}};
	
	dataFileName = "LPA N=" <> ToString[NValue];
	Export[directoryName <> dataFileName <> ".wdx", eigenvaluesLPA];
	
	PlotOneExponent[exponent_] := Block[{plotData = KeySort[errorsLPA[exponent]], fits, powerLaw, plot}, 
		plotData = Map[KeySort, plotData];

		fits = MapThread[Normal[NonlinearModelFit[Log[Normal2[#1][[-4;;-1]]], 2 #2 x + B, B, x]]&, 
			{Values[plotData], Keys[plotData]}];
		powerLaw = Map[Superscript[\[Epsilon], 2 #]&, Keys[plotData]];
		
		xticks = Sort[Map[{10^-#, If[Mod[#, 2] == 0, Superscript[10, -#], ""]} &, Range[1, 10]]];
		yticks = Sort[Map[{10^-(2#), If[Mod[#, 2] == 1, Superscript[10, -2#], ""]} &, Range[1, 7]]];
		ticks = {xticks, yticks};
		
		gridLines = {xticks[[1;;-1;;2, 1]], yticks[[1;;-1;;1, 1]]};
	
		plot = Show[
			MakeScatterLogLogPlot[Values[plotData], Map[ToString[2#+1] <> "-points"&, Keys[plotData]], Bottom, 
				{AxesLabel -> {\[Epsilon], Subscript[\[CapitalDelta], exponent]}, PlotRange -> plotRange, 
				 Ticks -> ticks, GridLines -> gridLines}, 
				{LegendLayout -> {"Row", 2}}], 
			MakeFitLogLogPlot[Exp[fits] /. x -> Log[\[CapitalDelta]], {\[CapitalDelta], 10^-10, 0.1}, powerLaw, Bottom, {}, {LegendLayout -> {"Row", 2}}]];
	
		Export[directoryName <> dataFileName <> " " <> ToString[exponent[[1]]] <> ToString[exponent[[2]]] <> ".pdf", plot];
	];
	Map[PlotOneExponent, {e1, e2, e3}];
	
	
	If[NValue == 2, 
		CopyFile[directoryName <> dataFileName <> " e1.pdf", articleFiguresDirectory <> "stability matrix LPA e1.pdf", OverwriteTarget -> True];
		CopyFile[directoryName <> dataFileName <> " e2.pdf", articleFiguresDirectory <> "stability matrix LPA e2.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
];


Map[TrackStabilityMatrixApproximationErrorLPA, {1, 2, 3}];


(* ::Subsubsection::Closed:: *)
(*DE2*)


TrackStabilityMatrixApproximationErrorDE2[nn_] := Block[{NValue = nn, derivativePointsCount = 7, gridSize, dimension = 3, 
	\[Rho]Max, alpha = 1, \[Rho]MaxTo\[Rho]0 = 35/10, gridSpacing = 1/20, 
	integratedEquations, parameters, scalingTerms, discreteEquations, fpSolution, 
	GetStabilityMatrixErrorDE2, epsilonValues, jacobianPointCount, arguments, 
	timeElapsed, eigenvaluesDE2, eigenvaluesDE2Nested, 
	SubtractExactValue, SubtractExactValueLogLog, 
	errorsDE2, logerrorsDE2, 
	flattenedErrors, minarg, maxarg, xmargin, 
	minval, maxval, ymargin, plotRange, 
	dataFileName, PlotOneExponent, potentialMinimum, FilterHighError, 
	subdirectory, directoryName, regulator,
	startTime, xticks, yticks, ticks, gridLines
	}, 
	
	startTime = AbsoluteTime[];
	Print["Tracking stability-matrix-approximation-error propagation at the DE2 level for N=" <> ToString[NValue] <> "."];

	subdirectory = "Stability Matrix Approximation Precision/";
	directoryName = Directory[] <> "/" <> subdirectory;
	If[!DirectoryQ[directoryName], 
		CreateDirectory[directoryName]];
	
	gridSize = Ceiling[\[Rho]MaxTo\[Rho]0 / gridSpacing];
	potentialMinimum = FindPotentialMinimum[NValue, "DE2", "exponential", False];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;

	SelectRegulator["exponential"];
	SetFixedFloatPrecision[doublePrecision];
	SetGLIntegralParameters[GLIntegralStandardParameters];
	integratedEquations = flowEquationsDE2;
	integratedEquations[[All, 3]] = PerformMomentumIntegrals[integratedEquations[[All, 3]]];
	integratedEquations[[All, 4]] = PerformMomentumIntegrals[integratedEquations[[All, 4]]];
	{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[integratedEquations, True];

	fpSolution = FindFixedPoint[scalingTerms + discreteEquations, guessDE2];
	
	GetStabilityMatrixErrorDE2[epsilon_, nPointsStabilityMatrix_] := 
		Block[{stabilityMatrix, eigenvalues, $MaxPrecision, $MinPrecision, $MaxExtraPrecision}, 
		
			stabilityMatrix = StabilityMatrix[flowEquationsDE2, fpSolution, True, nPointsStabilityMatrix, epsilon];

			eigenvalues = LeadingEigenvalues[stabilityMatrix];
			RestoreArbitraryPrecision[];
			Return[<|{e1, nPointsStabilityMatrix, epsilon} -> RecastFloatingPoint[eigenvalues[[1]]], 
				{e2, nPointsStabilityMatrix, epsilon} -> RecastFloatingPoint[eigenvalues[[2]]], 
				{e3, nPointsStabilityMatrix, epsilon} -> RecastFloatingPoint[eigenvalues[[3]]]|>];
		];
	
	epsilonValues = N[Join[{0}, 10^(-2-Range[0, 24]/3)]];
	jacobianPointCount = {1, 2, 3};
	arguments = Tuples[{epsilonValues, jacobianPointCount}];
	
	eigenvaluesDE2 = ParallelMap[GetStabilityMatrixErrorDE2[#[[1]], #[[2]]]&, arguments, Method -> "FinestGrained"];
	eigenvaluesDE2Nested = OneLevelNest[Association[eigenvaluesDE2]];
	
	SubtractExactValue[exponent_] := Association[Map[# -> Abs[eigenvaluesDE2Nested[exponent][#]-eigenvaluesDE2Nested[exponent][0.]]&, 
		Sort[Keys[eigenvaluesDE2Nested[exponent]]][[2;;-1]]]];
	SubtractExactValueLogLog[exponent_] := Association[Map[Log[#] -> Log[Abs[eigenvaluesDE2Nested[exponent][#]-eigenvaluesDE2Nested[exponent][0.]]]&, 
		Sort[Keys[eigenvaluesDE2Nested[exponent]]][[2;;-1]]]];

	errorsDE2 = NestAssociation[Association[Map[# -> SubtractExactValue[#]&, Keys[eigenvaluesDE2Nested]]]];
	logerrorsDE2 = NestAssociation[Association[Map[# -> SubtractExactValueLogLog[#]&, Keys[eigenvaluesDE2Nested]]]];
	
	FilterHighError[assoc_, threshold_:0.05] := Map[Select[#, #<threshold&]&, assoc];
	errorsDE2 = Map[FilterHighError, errorsDE2];
		
	flattenedErrors = UnnestAssociation[errorsDE2];
	minarg = Min[Map[#[[3]]&, Keys[flattenedErrors]]];
	maxarg = Max[Map[#[[3]]&, Keys[flattenedErrors]]];
	xmargin = (maxarg / minarg)^0.05;
	maxval = Min[Max[Map[Max, errorsDE2]], 10^-1];
	minval = Min[Map[Min, errorsDE2]];
	ymargin = (maxval / minval)^0.05;
	plotRange = {{minarg / xmargin, maxarg * xmargin}, {minval / ymargin, maxval * ymargin}};
	
	dataFileName = "DE2 N=" <> ToString[NValue];
	Export[directoryName <> dataFileName <> ".wdx", eigenvaluesDE2];
	
	PlotOneExponent[exponent_] := Block[{plotData = KeySort[errorsDE2[exponent]], fits, powerLaw, plot}, 
		plotData = Map[KeySort, plotData];

		fits = MapThread[Normal[NonlinearModelFit[Log[Normal2[#1][[-5;;-1]]], 2 #2 x + B, B, x]]&, {Values[plotData], Keys[plotData]}];
		powerLaw = Map[Superscript[\[Epsilon], 2 #]&, Keys[plotData]];
			
		xticks = Sort[Map[{10^-#, If[Mod[#, 2] == 0, Superscript[10, -#], ""]} &, 
			Range[1, 10]]];
		yticks = Sort[Map[{10^-(2#), If[Mod[#, 2] == 1, Superscript[10, -2#], ""]} &, 
			Range[1, 7]]];
		ticks = {xticks, yticks};
		
		gridLines = {xticks[[1;;-1;;2, 1]], yticks[[1;;-1;;1, 1]]};
	
		plot = Show[
			MakeScatterLogLogPlot[Values[plotData], Map[ToString[2#+1] <> "-points"&, Keys[plotData]], Bottom, 
				{AxesLabel -> {\[Epsilon], Subscript[\[CapitalDelta], exponent]}, PlotRange -> plotRange,
				Ticks -> ticks, GridLines -> gridLines}, 
				{LegendLayout -> {"Row", 2}}], 
			MakeFitLogLogPlot[Exp[fits] /. x -> Log[\[CapitalDelta]], {\[CapitalDelta], 10^-10, 0.1}, powerLaw, Bottom, {}, {LegendLayout -> {"Row", 2}}]];
	
		Export[directoryName <> dataFileName <> " " <> ToString[exponent[[1]]] <> ToString[exponent[[2]]] <> ".pdf", plot];
	];
	Map[PlotOneExponent, {e1, e2, e3}];
	
	If[NValue == 2, 
		CopyFile[directoryName <> dataFileName <> " e1.pdf", articleFiguresDirectory <> "stability matrix DE2 e1.pdf", OverwriteTarget -> True];
		CopyFile[directoryName <> dataFileName <> " e2.pdf", articleFiguresDirectory <> "stability matrix DE2 e2.pdf", OverwriteTarget -> True];
	];
	
	timeElapsed = AbsoluteTime[] - startTime;
	Print["Calculation completed in " <> ToString[Floor[timeElapsed/60]] <> " minutes and " <> 
		ToString[Mod[Floor[timeElapsed], 60]] <> " seconds. \n"];
];


Map[TrackStabilityMatrixApproximationErrorDE2, {1, 2, 3}];


(* ::Section::Closed:: *)
(*Extras*)


(* ::Subsection::Closed:: *)
(*Testing the accuracy of the Gauss-Legendre integrals*)


TestIntegrals[DEOrder_, nn_, integralParameters_] := Block[{NValue = nn, dimension = 3, alpha = 1, 
	\[Rho]MaxTo\[Rho]0 = 35/10, gridSpacing = 1/20, gridSize, \[Rho]Max, 
	regulator, $MaxPrecision, $MinPrecision, $MaxExtraPrecision, 
	integratedEquations, flowEquations, fpGuess, fpSolution, 
	parameters, scalingTerms, discreteEquations, potentialMinimum, betaFunctions, 
	exactIntegrals, approximateIntegrals, error
	}, 
	
	gridSize = Ceiling[\[Rho]MaxTo\[Rho]0 / gridSpacing];
	potentialMinimum = FindPotentialMinimum[NValue, DEOrder, "exponential", False];
	\[Rho]Max = \[Rho]MaxTo\[Rho]0*potentialMinimum;
	
	If[DEOrder == "LPA", 
		flowEquations = flowEquationsLPA;
		fpGuess = guessLPA, 
		
		flowEquations = flowEquationsDE2;
		fpGuess = guessDE2;
	];

	SelectRegulator["exponential"];
	SetFixedFloatPrecision[doublePrecision];
	SetGLIntegralParameters[integralParameters];
	
	{parameters, scalingTerms, discreteEquations} = DiscretizeFlowEquations[flowEquations, True];
	integratedEquations = PerformMomentumIntegrals[discreteEquations];
	
	fpSolution = FindFixedPoint[scalingTerms + integratedEquations, fpGuess];
	RestoreArbitraryPrecision[];
	
	fpSolution = Map[#[[1]] -> RecastFloatingPoint[#[[2]]]&, fpSolution];
	
	betaFunctions = ReplaceConstants[discreteEquations /. regulator /. fpSolution];
	
	exactIntegrals = NIntegrate[betaFunctions, {y, 0, 12}, PrecisionGoal -> doublePrecision+1, 
		MaxRecursion -> 20, WorkingPrecision -> longDoublePrecision];
	approximateIntegrals = PerformMomentumIntegrals[betaFunctions];
	error = Max[Abs[RecastFloatingPoint[exactIntegrals] - RecastFloatingPoint[approximateIntegrals]]];
	
	Print["The integration error for the O(" <> ToString[NValue] <> ") " <> DEOrder <> " \[Beta] functions on the interval |q|\[Element]["
		 <> ToString[integralParameters[[1]]] <> ", " <> ToString[integralParameters[[2]]] <> "] with "
		 <> ToString[integralParameters[[3]]] <> " integrand-evaluation points reached the value " 
		 <> ToString[ScientificForm[error, 3], StandardForm]];
]


(* ::Text:: *)
(*Testing low-precision parameters at the LPA level*)


Map[TestIntegrals["LPA", #, GLIntegralLowPrecisionParameters]&, {1, 2, 3}];


(* ::Text:: *)
(*Testing standard-precision parameters at the LPA level*)


Map[TestIntegrals["LPA", #, GLIntegralStandardParameters]&, {1, 2, 3}];


(* ::Text:: *)
(*Testing reference-precision parameters at the LPA level*)


Map[TestIntegrals["LPA", #, GLIntegralReferenceParameters]&, {1, 2, 3}];


(* ::Text:: *)
(*Testing low-precision parameters at the order DE2*)


Map[TestIntegrals["DE2", #, GLIntegralLowPrecisionParameters]&, {1, 2, 3}];


(* ::Text:: *)
(*Testing standard-precision parameters at the order DE2*)


Map[TestIntegrals["DE2", #, GLIntegralStandardParameters]&, {1, 2, 3}];


(* ::Text:: *)
(*Testing reference-precision parameters at the order DE2*)


Map[TestIntegrals["DE2", #, GLIntegralReferenceParameters]&, {1, 2, 3}];


(* ::Subsection::Closed:: *)
(*Finite-difference derivative operators presented in the article*)


(* ::Text:: *)
(*3-point derivative operators on a 5-point grid*)


derivativePointsCount=3;
gridSize=4;
Print[ToString[Subsuperscript[D, 3, 1], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff'[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];
Print[ToString[Subsuperscript[D, 3, 2], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff''[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];


(* ::Text:: *)
(*5-point derivative operators on a 7-point grid*)


derivativePointsCount=5;
gridSize=6;
Print[ToString[Subsuperscript[D, 3, 1], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff'[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];
Print[ToString[Subsuperscript[D, 3, 2], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff''[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];


(* ::Text:: *)
(*7-point derivative operators on a 9-point grid*)


derivativePointsCount=7;
gridSize=8;
Print[ToString[Subsuperscript[D, 3, 1], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff'[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];
Print[ToString[Subsuperscript[D, 3, 2], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff''[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];


(* ::Text:: *)
(*9-point derivative operators on a 11-point grid*)


derivativePointsCount=9;
gridSize=10;
Print[ToString[Subsuperscript[D, 4, 1], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff'[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];
Print[ToString[Subsuperscript[D, 4, 2], StandardForm] <> " = " <> 
	ToString[MatrixForm[Table[D[NumericDerivatives[ff''[i]], ff[j]], {i, 0, gridSize}, {j, 0, gridSize}]], StandardForm]];
