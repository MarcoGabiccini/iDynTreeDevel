(* ::Package:: *)

(* 
 *        iDynTree.m - Lie Geometric Algorithms for the Kinematics, Dynamics and
                       Control of Multibody Systems (MBS)
                       
          by Marco Gabiccini and Silvia Manara
          
          Dipartimento di Ingegneria Civile ed Industriale (DICI)
          Scuola di Ingegneria
          Universita' di Pisa
          56122 Pisa PI - Italy
          
          Copyright (c) 2014, Marco Gabiccini             
 *
 *)
 
Needs["NumericalCalculus`"];
BeginPackage["iDynTreeDevel`"]
 
 (* Error Messages *)

AxisToSkew::wrongD           = "The first argument is not a 3D vector.";
SkewToAxis::notskewsymmetric = "The first argument is not a skew-symmetric matrix";
Screws::wrongDimensions      = "`1 argument`: Dimensions of input matrices incorrect.";
Screws::notSquare            = "`1 argument`: Input matrix is not square.";
Screws::notVector            = "`1 argument`: Input is not a vector.";
 
 (* Matrix Utility Functions *)
 
AppendColumns::usage=
"AppendColumns[mat1, mat2, ...] gives a new matrix composed of the
submatrices mat1, mat2, ... , by joining their columns. The submatrices must
all have the same number of columns."

AppendRows::usage=
"AppendRows[mat1, mat2, ...] gives a new matrix composed of the submatrices
mat1, mat2, ..., by joining their rows. The submatrices must all have the same
number of rows."

StackCols::usage=
"StakCols[mat1, mat2, ...] gives a new matrix obtained by stacking the columns 
of the submatrices mat1, mat2, ... . The submatrices must all have the same number
of rows."

StackRows::usage=
"StackRows[mat1, mat2, ...] gives a new matrix obtained by stacking the rows
of the submatrices mat1, mat2, ... . The submatrices must all have the same number 
of columns."

BlockMatrix::usage=
"BlockMatrix[block] gives a matrix composed of a matrix of matrices."

BlockDiag::usage=
"BlockDiag[list] gives a block-diagonal matrix composed of the matrices listed in list."

TakeRows::usage=
"TakeRows[mat, part] takes the rows of matrix mat given by the partspec,
part. This uses the same notation as Take."

TakeColumns::usage=
"TakeColumns[mat, part] takes the columns of matrix mat given by the
partspec, part. This uses the same notation as Take."

TakeMatrix::usage=
"TakeMatrix[mat, start, end] gives the submatrix starting at position start
and ending at position end, where each position is expressed as a 
{row, column} pair, e.g., {1, 2}."

SubMatrix::usage=
"SubMatrix[mat, start, dim] gives the submatrix of dimension dim = {m, n} starting
at position start = {i, j} in mat."

ZeroMatrix::usage=
"ZeroMatrix[n] gives the nxn zero matrix."

Eye::usage=
"Eye[n] gives the nxn identity matrix."

SelectionMatrix::usage=
"SelectionMatrix[q_List, indexlist_List] or SelectionMatrix[n_Integer, indexlist_List] returns the selection matrix with n=Length[q] columns that returns the
indexlist components of q."

SelectionMatrixColumns::usage=
"SelectionMatrixColumns[n, indexlist] returns the matrix that selects the indexlist columns from a matrix with n columns."

Magnitude::usage=
 "Magnitude[v] gives length of v."

 (* Rotation Matrices - forward kin map *)
 
RotX::usage=
"RotX[alpha] gives the rotation matrix about the X axis.";

RotY::usage=
"RotY[alpha] gives the rotation matrix about the Y axis.";

RotZ::usage=
"RotZ[alpha] gives the rotation matrix about the Z axis.";

ATan2::usage =
"ATan2[y, x] gives the angle with tangent given by y/x in the right quadrant."

EulZXZToMat::usage=
"EulZXZToMat[alpha, beta, gamma] returns the rotation matrix corresponding to the ZXZ Euler angles.
This is the one of the standard parametrizations employed in Multibody Dynamics."

EulZYZToMat::usage=
"EulZYZToMat[alpha, beta, gamma] returns the rotation matrix corresponding to the ZYZ Euler angles.
This is the one of the standard parametrizations employed in Robotics."

EulZYXToMat::usage = RPYToMat::usage=
"EulZYXToMat[alpha, beta, gamma] returns the rotation matrix corresponding to the ZYX Euler angles,
also known as Roll, Pitch and Yaw.
This is one of the standard parametrizations employed in Robotics."

EulXYZToMat::usage=
"EulXYZToMat[alpha, beta, gamma] returns the rotation matrix corresponding to the XYZ Euler angles.
This is useful, accompanied by its corresponding spatial Jacobian for small rotations."

RodriguezToMat::usage=
"RodriguezToMat[g1, g2, g3] returns the rotation matrix corresponding to the Rodriguez parameters g1, g2, g3.\
Recall that {g1, g2, g3} = r tan(theta/2)." 

QuatToMat::usage=
"QuatToMat[q] returns the rotation matrix corresponding to the unit quaternion q (1-vector)."

 (* Rotation Matrices - forward differential kin map *)

EulZXZToSpatialJac::usage=
"EulZXZToSpatialJac[alpha, beta, gamma] returns the Spatial Jacobian Js corresponding to the ZXZ Euler angles.
This is such that w_s = Js dot{alpha, beta, gamma}, with w_s spatial components."

EulZXZToBodyJac::usage=
"EulZXZToBodyJac[alpha, beta, gamma] returns the Body Jacobian Jb corresponding to the ZXZ Euler angles.
This is such that w_b = Jb dot{alpha, beta, gamma}, with w_b body fixed components."

EulZYZToSpatialJac::usage=
"EulZYZToSpatialJac[alpha, beta, gamma] returns the Spatial Jacobian Js corresponding to the ZYZ Euler angles.
This is such that w_s = Js dot{alpha, beta, gamma}, with w_s spatial components."

EulZYZToBodyJac::usage=
"EulZYZToBodyJac[alpha, beta, gamma] returns the Body Jacobian Jb corresponding to the ZYZ Euler angles.
This is such that w_b = Jb dot{alpha, beta, gamma}, with w_b body fixed components."

EulZYXToSpatialJac::usage=
"EulZYXToSpatialJac[alpha, beta, gamma] returns the Spatial Jacobian Js corresponding to the ZYX Euler angles.
This is such that w_s = Js dot{alpha, beta, gamma}, with w_s spatial components."

EulZYXToBodyJac::usage=
"EulZYXToBodyJac[alpha, beta, gamma] returns the Body Jacobian Jb corresponding to the ZYX Euler angles.
This is such that w_b = Jb dot{alpha, beta, gamma}, with w_b body fixed components."

EulXYZToSpatialJac::usage=
"EulXYZToSpatialJac[alpha, beta, gamma] returns the Spatial Jacobian Js corresponding to the XYZ Euler angles.
This is such that w_s = Js dot{alpha, beta, gamma}, with w_s spatial components."

QuatToSpatialJac::usage=
"QuatToSpatialJac[q] returns the Spatial Jacobian Js corresponding to the unit quaternion q (1-vector).
If q = (q0, q1, q2, q3), w_s = Js dot{q}."

QuatToBodyJac::usage=
"QuatToBodyJac[q] returns the Body Jacobian Jb corresponding to the unit quaternion q (1-vector).
If q = (q0, q1, q2, q3), w_b = Jb dot{q}."

RodriguezToSpatialJac::usage=
"RodriguezToSpatialJac[g1, g2, g3] returns the Spatial Jacobian Js corresponding to Rodriguez parameters.\
This is such that w_s = Js dot{g1, g2, g3}, with w_s spatial components."

RodriguezToBodyJac::usage=
"RodriguezToBodyJac[g1, g2, g3] returns the Body Jacobian Jb corresponding to Rodriguez parameters.\
This is such that w_b = Jb dot{g1, g2, g3}, with w_b body fixed components."

 (* Utilities for Testing Matrix Properties *)
 
RotationQ::usage=
"RotationQ[m] tests whether matrix m is a rotation matrix.";

skewQ::usage=
"skewQ[m] tests whether matrix m is a skew-symmetrix matrix."; 

 (* Rotation Matrices - inverse kin map *)

MatToEulZYZ::usage=
"MatToEulZYZ[R] returns the Euler angles for the parametrization ZYZ corresponding to the rotation matrix R."

MatToEulZYX::usage=MatToRPY::usage=
"MatToEulZYX[R] returns the Euler angles for the parametrization ZYX (also called RPY) corresponding to the rotation matrix R."

MatToVect::usage=
"MatToVect[R] or MatToVect[ {R, Rd} ] returns the orientation error e_O = r sin(theta), where r and theta are the components of the axis angle parametrization.
e_O represents also the axis form of the skew-symm matrix R_SS = (1/2) (R - R^T)."

FramesToVect::usage=
"FramesTovect[{ R, Rd }] returns the orientation error e_O = r sin(theta) of frame Rd (desired) w.r.t. R (actual) in the spatial frame\
where both R and Rd are expressed."

 (* Homogeneous Representations *)

HomogeneousToVector::usage=
"HomogeneousToVector[q] extracts the cartesian coordinates from homogeneous components q.";

PointToHomogeneous::usage=
"PointToHomogeneous[q] gives the homogeneous representation of a point.";

VectorToHomogeneous::usage=
"VectorToHomogeneous[q] gives the homogeneous representation of a vector.";

RPToHomogeneous::usage=
"RPToHomogeneous[R,p] forms homogeneous matrix from rotation matrix R \
 and position vector p.";
  
RigidOrientation::usage=
"RigidOrientation[g] extracts rotation matrix R from g.";

RigidPosition::usage=
"RigidPosition[g] extracts position vector p from g.";
  
RigidInverse::usage=
"RigidInverse[g] gives the inverse transformation of g.";  

RotationAxis::usage=
"RotationAxis[R] returns the rotation axis of R in SO(3)."  
  
RotationParam::usage=
 "RotationParam[R] returns the rotation axis the and amount of rotation of R in SO(3)."  

 (* Lie Algebra Operations *)

AxisToSkew::usage=Skew::usage=Hat::usage=
  "AxisToSkew[w] generates skew-symmetric matrix given 3 vector w";

SkewToAxis::usage=UnSkew::usage=HatInv::usage=Vee::usage=
  "SkewToAxis[S] extracts vector w from skew-symmetric matrix S";

SkewExp::usage=
 "SkewExp[w,(theta)] gives the matrix exponential of an axis w.
  Default value of Theta is 1.";

xitow::usage=
"xitow[xi] gives the angular part of twist xi."

xitov::usage=
"xitov[xi] gives the translational part of twist xi."

(* Operations on se(e), the Lie algebra of SE(3) *)

TwistExp::usage=
 "TwistExp[xi,(Theta)] gives the matrix exponential of a twist xi.
  Default value of Theta is 1.";

TwistToHomogeneous::usage=
  "TwistToHomogeneous[xi] converts xi from a 6 vector to a 4x4 matrix.";

HomogeneousToTwist::usage=
  "HomogeneousToTwist[xi] converts xi from a 4x4 matrix to a 6 vector.";

RigidTwist::usage=
  "RigidTwist[g] extracts 6 vector xi from g and the angle theta that generates g.";

TwistPitch::usage=
  "TwistPitch[xi] gives pitch of screw corresponding to a twist.";

TwistAxis::usage=
  "TwistAxis[xi] gives axis of screw corresponding to a twist as q (point on axis), w (axis)";

TwistMagnitude::usage=
  "TwistMagnitude[xi] gives magnitude of screw corresponding to a twist.";
  
ScrewToTwist::usage=
  "ScrewToTwist[h, q, w] builds a twist from its screw coordinates h (pitch), q (point on axis), w (axis).";  

(* Adjoint matrix calculations *)

RigidAdjoint::usage = Ad::usage = 
  "RigidAdjoint[g] / Ad[g] gives the Adjoint matrix corresponding to g.";
  
InverseRigidAdjoint::usage = AdInv::usage =
  "InverseRigidAdjoint[g] / AdInv[g] gives the inverse Adjoint matrix corresponding to g.";  

TransposeRigidAdjoint::usage = AdTr::usage =
  "TransposeRigidAdjoint[g] / AdTr[g] gives the transpose Adjoint matrix corresponding to g.";  

InverseTransposeRigidAdjoint::usage = AdStar::usage =
  "InverseTransposeRigidAdjoint[g] / AdStar[g] gives the inverse transpose Adjoint matrix corresponding to g."

(* adjoint (small) matrix calculations - aka Lie Bracket *)

ad::usage = 
  "ad[xi] gives the adjoint matrix associated to twist xi."
  
adTr::usage =
  "adTr[xi] gives the transpose adjoint matrix associated to xi."
  
adStar::usage =
  "adStar[xi] gives the negative transpose adjoint matrix associated to xi."    

LieBracketDull::usage = 
  "LieBracket[xia, xib] / LieBracket[xia^, xib^] gives Vee[(xia^ xib^ - xib^ xia^)]."

LieBracket::usage = 
  "LieBracket[xia, xib] / LieBracket[xia^, xib^] gives Vee[(xia^ xib^ - xib^ xia^)]."

(* Operations on SO(n), the Lie group of rotation matrices *)

LocalTranTwist::usage =
"LocalTranTwist[gStart, gEnd] returns the twist and theta such that: gStart.TwistExp[xi,theta] = gEnd.
 The xi coordinate are expressed in the local (gStart) reference frame."

GlobalTranTwist::usage =
"GlobalTranTwist[gStart, gEnd] returns the twist and theta such that: TwistExp[xi,theta].gStart = gEnd.
The xi coordinate are expressed in the global (identity) reference frame."

 (* Systematic Methods to Build Forward Kinematics and Spatial/Body Jacobians *)

RevoluteTwist::usage=
 "RevoluteTwist[q, w] builds the 6-vector corresponding to point q on the axis with unit vector w for a revolute joint.";

PrismaticTwist::usage=
 "RevoluteTwist[q, w] builds the 6-vector corresponding to point q on the axis with unit vector w for a prismatic joint.";

ForwardKinematics::usage=
 "ForwardKinematics[{xi1,th1}, ... ,{xiN,thN}, g0] computes the forward kinematics via the product of exponentials formula.
 g0 is the initial affine tranformation if any.";

ForwardKinematicsLeft::usage=
"ForwardKinematicsLeft[g0, {xi1,th1}, ... ,{xiN,thN}] computes the forward kinematics via the product of exponentials formula, premultiplicating the offset matrix.
 g0 is the initial affine tranformation if any.";

ForwardKinematicsLocalPOE::usage=
 "ForwardKinematicsLocalPOE[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the forward kinematics via the Local Product of Exponentials formula.
 The outcome is g = M1.exp(X1*q1) ... Mn.exp(Xn*qn).";

ForwardKinematicsLocalPOERight::usage=
"ForwardKinematicsLocalPOERight[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the forward kinematics via the Local Product of Exponentials formula.
 Each offset matrix is post-multiplicated.
 The outcome is g = exp(X1*q1).M1 ... exp(Xn*qn).Mn .";

ForwardKinematicsLocalPOETranslated::usage=
"ForwardKinematicsLocalPOETranslated[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the forward kinematics via the Local Product of Exponentials formula.
 The outcome is g = M1.exp(X1*q1) ... Mn.exp(Xn*qn). Please note that this result is computed by initially calculating the Global version of the twists
 and then applying a formula similar to the one used in ForwardKinematics.";

ForwardKinematicsDH::usage=
"ForwardKinematicsDH[{theta1,d1,alpha1,a1}, ... ,{thetaN,dN,alphaN,aN}] computes the forward kinematics via the Denavit-Hartenberg formula.";

SpatialJacobianGlobalPOE::usage=
"SpatialJacobianGlobalPOE[{xi1,th1}, ... ,{xiN, thN} ,g0] computes the Spatial Manipulator Jacobian of a robot defined by the given twists
described in the Global POE form.";

BodyJacobianGlobalPOE::usage=
"BodyJacobianGlobalPOE[{xi1,th1}, ..., {xiN, thN}, g0] computes the Body Manipulator Jacobian of a robot defined by the given twists.";

BodyJacobianGlobalPOELeft::usage=
"BodyJacobianGlobalPOELeft[g0, {xi1,th1}, ..., {xiN, thN}] computes the Body Manipulator Jacobian of a robot defined by the Global POE Left form.";

SpatialJacobianGlobalPOELeft::usage=
"SpatialJacobianGlobalPOELeft[g0, {xi1,th1}, ..., {xiN, thN}] computes the Spatial Manipulator Jacobian of a robot defined by the Global POE Left form.";

SpatialJacobianLocalPOERight::usage=
"SpatialJacobianLocalPOERight[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the Spatial Manipulator Jacobian of a robot defined by the given twists
described in the Local POE Right form.";

BodyJacobianLocalPOERight::usage=
"BodyJacobianLocalPOERight[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the Body Manipulator Jacobian of a robot defined by the given twists
described in the Local POE Right form.";

SpatialJacobianLocalPOE::usage=
"SpatialJacobianLocalPOE[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the Spatial Manipulator Jacobian of a robot defined by the given twists
described in the Local Product of Exponentials form.";

BodyJacobianLocalPOE::usage=
"BodyJacobianLocalPOE[{M1, X1, q1}, ... ,{Mn, Xn, qn}] computes the Body Manipulator Jacobian of a robot defined by the given twists
described in the Local Product of Exponentials form.";

SpatialTwistDerivative::usage=
"SpatialTwistDerivative[{X1, q1, q1p, q1pp}, ... {Xn, qn, qnp, qnpp}, gst0] computes the derivative of a Spatial Twist given the twists
expressed in the Global frame, the joint angles, joint velocities and joint accelerations."

SpatialTwistDerivativeLocal::usage=
"SpatialTwistDerivativeLocal[{M1, X1, q1, q1p, q1pp}, ... {Mn, Xn, qn, qnp, qnpp}] computes the spatial derivative of a Twist given the twists
expressed in the Local frame, the joint angles, joint velocities and joint accelerations."

BodyTwistDerivative::usage=
"BodyTwistDerivative[{X1, q1, q1p, q1pp}, ... {Xn, qn, qnp, qnpp}, gst0] computes the derivative, expressed in the Body fixed reference frame, of a Spatial Twist given the twists
expressed in the Global frame, the joint angles, joint velocities and joint accelerations."

BodyTwistDerivativeLocal::usage=
"BodyTwistDerivativeLocal[{M1, X1, q1, q1p, q1pp}, ... {Mn, Xn, qn, qnp, qnpp}] computes the derivative, expressed in the Body fixed reference frame, of a Twist given the twists
expressed in the Local frame, the joint angles, joint velocities and joint accelerations."

FBSpatialJacobianGlobalPOE::usage=
"FBSpatialJacobianGlobalPOE[{{x,quat},{Y2,q2}, ... {Yn,qn}}, g1n(0), JointChainList, Nj] computes the Spatial Jacobian of a rigid body inside a Floating based, Branched kinematic tree,
given the description of the joints (position vector and quaternion for the 6DoF joint, twists in the Global frame for the following), the list of the joints which belong to the kinematic chain
between the floating base and the body we are considering and the total number of joints belonging to the kinematic tree."

FBSpatialJacobianLocalPOE::usage=
"FBSpatialJacobianLocalPOE[{{x,quat},{M2,X2,q2}, ... {Mn,Xn,qn}}, JointChainList, Nj] computes the Spatial Jacobian of a rigid body inside a Floating based, Branched kinematic tree,
given the description of the joints position vector and quaternion for the 6DoF joint, twists in the Local frame for the following), the list of joints which belong to the kinematic 
chain between the floating base and the body we are considering and the total number of joints belonging to the kinematic tree."

FBBodyJacobianGlobalPOE::usage=
"FBBodyJacobianGlobalPOE[{{x,quat},{Y2,q2}, ... {Yn,qn}}, g1n(0), JointChainList, Nj] computes the Body Jacobian of a rigid body inside a Floating based, Branched kinematic tree,
given the description of the joints (position vector and quaternion for the 6DoF joint, twists in the Global frame for the following), the list of the joints which belong to the kinematic chain
between the floating base and the body we are considering and the total number of joints belonging to the kinematic tree."

FBBodyJacobianLocalPOE::usage=
"FBBodyJacobianLocalPOE[{{x,quat},{M2,X2,q2}, ... {Mn,Xn,qn}}, JointChainList, Nj] computes the Body Jacobian of a rigid body inside a Floating based, Branched kinematic tree, given
the description of the joints (position vector and quaternion for the 6DoF joint, twists in the Local frame for the following), the list of the joint which belong to the kinematic
chain between the floating base and the body we are considering and the total number of joints belonging to the kinematic tree." 

SpatialJacobianDerivative::usage=
"SpatialJacobianDerivative[Js,{q1dot, ... , qndot}] computes the Spatial derivative of the Jacobian of a kinematic chain, given the Spatial Jacobian and the velocities 
of the joint variables."

BodyJacobianDerivative::usage=
"BodyJacobianDerivative[Jb,{q1dot, ... , qndot}] computes che Body derivative of the Jacobian of a kinemati chain, given the Body Jacobian and the velocities of the joint
variables."

FBSpatialJacobianDerivative::usage=
"FBSpatialJacobianDerivative[FBJs, {xdot, quatdot, q2dot, ... , qndot}] computes the Spatial derivative of the Jacobian of a rigid body inside a Floating based, Branched kinematic
tree, given its Spatial Jacobian and the velocities of the joint variables."

FBBodyJacobianDerivative::usage=
"FBBodyJacobianDerivative[FBJb, {xdot, quatdot, q2dot, ... , qndot}] computes the Body derivative of the Jacobian of a rigid body inside a Floating based, Branched kinematic
tree, given its Body Jacobian and the velocities of the joint variables."


(******************************************************************************
 ******************************************************************************
 ******************** PRIVATE SECTION *****************************************
 ******************************************************************************
 ******************************************************************************
*) 
 
 
Begin["Private`"]

(* Matrix Dimensions and Type Test Functions *)

SameColumnSize[l_List] := (SameQ @@ (Dimensions[#][[2]]& /@ l) )

SameRowSize[l_List] := (SameQ @@ (Dimensions[#][[1]]& /@ l) )

RotationQ[mat_] :=
  Module[
    {nr, nc, zmat},

    (* First check to make sure that this is square matrix *)
    If[Not[MatrixQ[mat]] || ({nr, nc} = Dimensions[mat]; nr != nc),
	Message[Screws::notSquare];    
        Return[False]];

    (* Check to see if R^T R = Identity *)
    zmat = Simplify[mat . Transpose[mat]] - IdentityMatrix[nr];
    Return[ And @@ Map[TrueQ[Simplify[#] == 0]&, Flatten[zmat]]]
  ];

skewQ[mat_] :=
  Module[
    {nr, nc, zmat},

    (* First check to make sure that this is square matrix *)
    If[Not[MatrixQ[mat]] || ({nr, nc} = Dimensions[mat]; nr != nc),
	Message[Screws::notSquare];    
        Return[False]];

    (* Check to see if A = -A^T *)
    zmat = mat + Transpose[mat];
    Return[ And @@ Map[TrueQ[Simplify[#] == 0]&, Flatten[zmat]]]
];


(* Matrix Utility Functions *)
 
AppendColumns[l__?MatrixQ] := Join[l] /; SameColumnSize[{l}] 

AppendRows[l__?MatrixQ] := MapThread[Join, {l}] /; 
   SameRowSize[{l}]
   
StackCols[l__?MatrixQ] := AppendRows[l];   

StackRows[l__?MatrixQ] := AppendColumns[l];

BlockMatrix[block_] :=
	AppendColumns @@ Apply[AppendRows, block, {1}]; 

BlockDiag[list_] :=
	Module[{r, mlist, nlist, m, n, i, res},
		r = Length[list];
		mlist = Dimensions[#][[1]]& /@ list;
		nlist = Dimensions[#][[2]]& /@ list;
		m = Plus @@ mlist;
		n = Plus @@ nlist;
		res = ZeroMatrix[m, n];
		
		res[[ 1;;mlist[[1]], 1;;nlist[[1]] ]]= list[[1]];
		
		
		For[ i=2, i <= r, ++i,
			res[[ Plus@@Take[mlist, i-1]+1 ;; Plus@@Take[mlist, i], Plus@@Take[nlist, i-1]+1 ;; Plus@@Take[nlist, i] ]] = list[[i]]
		];
		
		Return[res];
	];

TakeRows[mat_?MatrixQ, part_] := 
	Take[mat, part];

TakeColumns[mat_?MatrixQ, part_] := 
	Take[mat, All, part];

TakeMatrix[mat_?MatrixQ, start:{startR_Integer, startC_Integer},
end:{endR_Integer, endC_Integer}] :=
	Take[mat, {startR, endR}, {startC, endC}] /;
	And @@ Thread[Dimensions[mat] >= start] && 
	And @@ Thread[Dimensions[mat] >= end]

SubMatrix[mat_List, start:{_Integer, _Integer}, dim:{_Integer,_Integer}] :=
	TakeMatrix[mat, start, start+dim-1]; 

ZeroMatrix[0, ___] := {};
ZeroMatrix[m_Integer, 0] := Table[{}, {m}];

ZeroMatrix[m_Integer,n_Integer] := 
	Normal[SparseArray[{},{m, n}]] /; m >= 0 && n>=0

ZeroMatrix[m_Integer] := ZeroMatrix[m, m] /; m >= 0 

Eye[m_Integer] := IdentityMatrix[m] /; m >= 1

SelectionMatrix[q_List, indexlist_List] :=
	Module[{m, eye},
		m = Length[q];
		eye = IdentityMatrix[m];
		res = eye[[indexlist]];
		Return[res];
	];

SelectionMatrix[m_Integer, indexlist_List] :=
	Module[{eye},
		eye = IdentityMatrix[m];
		res = eye[[indexlist]];
		Return[res];
	];

SelectionMatrixColumns[n_Integer, indexlist_List] :=	
	Module[{eye},
		eye = IdentityMatrix[n];
		res = Transpose[ eye[[indexlist]] ];
		Return[res];
	];

Magnitude[v_] :=
 Module[
  {},
  
  If[Not[VectorQ[v]],
    Message[Screws::wrongDimensions, "Vector"];
    Return[Null];
  ];

  Sqrt[v.v]
];

 (* Rotation Matrices - forward kin map *)
 
RotX[alpha_] :=
  Module[ {ca = Cos[alpha], sa = Sin[alpha]},
         {{1,  0, 0},
          {0, ca, -sa}, 
          {0, sa,  ca}
          }
  ];

RotY[alpha_] :=
  Module[ {ca = Cos[alpha], sa = Sin[alpha]},
         {{ca,  0, sa},
          {0,   1,  0}, 
          {-sa, 0,  ca}
          }
  ];

RotZ[alpha_] :=
  Module[ {ca = Cos[alpha], sa = Sin[alpha]},
         {{ca, -sa, 0},
          {sa,  ca, 0}, 
          {0,    0, 1}
          }
  ];
     
ATan2[y_, x_] := ArcTan[x, y];

EulZXZToMat[alpha_, beta_, gamma_]:=
	Module[ {},
		
	RotZ[alpha].RotX[beta].RotZ[gamma]
	
	];

EulZYZToMat[alpha_, beta_, gamma_]:=
	Module[ {},
		
	RotZ[alpha].RotY[beta].RotZ[gamma]
	
	];

EulZYXToMat[alpha_, beta_, gamma_]:=
	Module[ {},
		
	RotZ[alpha].RotY[beta].RotX[gamma]
	
	];
	
RPYToMat[alpha_, beta_, gamma_]:= EulZYXToMat[alpha, beta, gamma];

EulXYZToMat[alpha_, beta_, gamma_]:=
	Module[ {},
		
	RotX[alpha].RotY[beta].RotZ[gamma]
	
	];

(* Rodriguez parameters gamma = r tan(theta/2) *)
RodriguezToMat[gamma1_, gamma2_, gamma3_] :=
	Module[ {gamma, hatgamma, modulusgammasquared},
	
			gamma = {gamma1, gamma2, gamma3};
			hatgamma = Hat[gamma];
			modulusgammasquared = gamma.gamma;
			
			IdentityMatrix[3] + 2/(1 + modulusgammasquared) (hatgamma + hatgamma.hatgamma)
	      ]; 

QuatToMat[qList_]:=
Module[ {b0=qList[[1]], b1=qList[[2]], b2=qList[[3]], b3=qList[[4]]},
{{b0^2+b1^2-b2^2-b3^2, 2*(b1*b2-b0*b3), 2*(b0*b2+b1*b3)},{2*(b1*b2-b0*b3),b0^2-b1^2+b2^2-b3^2,2*(b2*b3-b0*b1)},{2*(b1*b3-b0*b2),2*(b0*b1+b2*b3),b0^2-b1^2-b2^2+b3^2}}
];

 (* Rotation Matrices - forward differential kin map *)

EulZXZToSpatialJac[phi_, theta_, psi_]:=
	Module[{},
		{{ 0, Cos[phi],  Sin[phi] Sin[theta] },
		 { 0, Sin[phi], -Cos[phi] Sin[theta] },
		 { 1,        0,           Cos[theta] }
		}
		
	];
	
EulZXZToBodyJac[phi_, theta_, psi_]:=
	Module[{},
		{{ Sin[theta] Sin[psi],  Cos[psi], 0 },
		 { Sin[theta] Cos[psi], -Sin[psi], 0 },
		 {          Cos[theta],         0, 1 }
		}
		
	];	
	
EulZYZToSpatialJac[phi_, theta_, psi_]:=
	Module[{},
		{{ 0, -Sin[phi], Cos[phi] Sin[theta]  },
		 { 0,  Cos[phi], Sin[phi] Sin[theta]  },
		 { 1,        0,           Cos[theta]  }
		}
		
	];	

EulZYZToBodyJac[phi_, theta_, psi_]:=
	Module[{},
		{{ -Cos[psi] Sin[theta], Sin[psi], 0 },
		 {  Sin[psi] Sin[theta], Cos[psi], 0 },
		 {           Cos[theta],        0, 1 }
		}
		
	];	
	
EulZYXToSpatialJac[phi_, theta_, psi_]:=
	Module[{},
		{{ 0, -Sin[phi], Cos[theta] Cos[phi]  },
		 { 0,  Cos[phi], Cos[theta] Sin[phi]  },
		 { 1,        0,          -Sin[theta]  }
		}
		
	];		
	
EulZYXToBodyJac[phi_, theta_, psi_]:=
	Module[{},
		{{          -Sin[theta],         0,  1 },         
		 {  Cos[theta] Sin[psi],  Cos[psi],  0 },
		 {  Cos[theta] Cos[psi], -Sin[psi],  0 }
		}
		
	];	
	
EulXYZToSpatialJac[phi_, theta_, psi_]:=
	Module[{},
		
		{{1, 	0,	  		 Sin[theta] 	     },
		 {0,	Cos[phi],	-Cos[theta] Sin[phi] },
		 {0,	Sin[phi],    Cos[theta] Cos[phi]}
		}
	];		

QuatToSpatialJac[ bList_] :=
	Module[{b0, b1, b2, b3},
	b0 = bList[[1]];
	b1 = bList[[2]];
	b2 = bList[[3]];
	b3 = bList[[4]];
	
	2 * { { -b1,  b0, -b3,  b2  },
		  { -b2,  b3,  b0, -b1 },
		  { -b3, -b2,  b1,  b0 }
	    }
	
	];

QuatToBodyJac[ bList_ ] :=
	Module[{b0, b1, b2, b3},
	b0 = bList[[1]];
	b1 = bList[[2]];
	b2 = bList[[3]];
	b3 = bList[[4]];
	
	2 * { { -b1,  b0,  b3, -b2  },
		  { -b2, -b3,  b0,  b1 },
		  { -b3,  b2, -b1,  b0 }
	    }
	
	];	
 
RodriguezToSpatialJac[gamma1_, gamma2_, gamma3_] :=
	Module[ {gamma, modulusgammasquared},
		gamma = {gamma1, gamma2, gamma3};
		modulusgammasquared = gamma.gamma;
		2/(1 + modulusgammasquared) { {      1,   -gamma3,    gamma2},
									  { gamma3,         1,   -gamma1},
									  {-gamma2,    gamma1,         1}
									}  	
	
	];
	
RodriguezToBodyJac[gamma1_, gamma2_, gamma3_] :=
	Module[ {gamma, modulusgammasquared},
		gamma = {gamma1, gamma2, gamma3};
		modulusgammasquared = gamma.gamma;
		2/(1 + modulusgammasquared) { {      1,   gamma3,    -gamma2},
									  { -gamma3,         1,   gamma1},
									  {  gamma2,    -gamma1,         1}
									}  	
	
	];
	
 (* Rotation Matrices - inverse kin map *)		 

MatToEulZYZ[R_] :=
	Module[{phi, theta, psi, thetaiszero, thetaisPi},
	        
	(* Check the singularity of representation ZYZ *)        
	thetaiszero = Abs[ R[[3,3]] - 1] < 10^(-10);
	thetaisPi = Abs[R[[3,3]] + 1] < 10^(-10);     
	        
	(* In cases of singularity we set arbitrarily psi = 0 *)        
	        
	If[ thetaiszero,  
	     phi = ATan2[ R[[2,1]] , R[[1,1]] ];
	     theta = 0;
	     psi = 0;
	     ];
	  
	If[ thetaisPi,
	     phi = ATan2[ R[[2,1]], R[[1,1]] ];
	     theta = Pi;
	     psi = 0;
	     ];  
	   
	If[ !(thetaiszero || thetaisPi),
	     phi = ATan2[R[[2,3]], R[[1,3]]];
	     theta = ATan2[Sqrt[ R[[1,3]]^2 + R[[2,3]]^2 ], R[[3,3]] ];
	     psi = ATan2[R[[3,2]], -R[[3,1]]];
	     ];
	     
	     Return[{phi, theta, psi}];
	     
	];

(* with theta \in (-Pi/2, Pi/2) *)
MatToEulZYX[R_] :=
	Module[{phi, theta, psi, thetaisplushalfPi, thetaisminushalfPi},
	        
	(* Check the singularity of representation ZYX *)        
	thetaisplushalfPi = Abs[ R[[3,1]] + 1] < 10^(-10);
	thetaisminushalfPi = Abs[R[[3,1]] - 1] < 10^(-10);     
	        
	(* In cases of singularity we set arbitrarily psi = 0 *)        
	        
	If[ thetaisplushalfPi,  
	     phi = ATan2[ R[[2,3]] , R[[1,3]] ];
	     theta = Pi/2;
	     psi = 0;
	     ];
	  
	If[ thetaisminushalfPi,
	     phi = ATan2[ -R[[2,3]], -R[[1,3]] ];
	     theta = -Pi/2;
	     psi = 0;
	     ];  
	   
	If[ !(thetaisplushalfPi || thetaisminushalfPi),
	     phi = ATan2[R[[2,1]], R[[1,1]]];
	     theta = ATan2[ -R[[3,1]], Sqrt[ R[[3,2]]^2 + R[[3,3]]^2 ] ];
	     psi = ATan2[ R[[3,2]], R[[3,3]] ];
	     ];
	     
	     Return[{phi, theta, psi}];
	     
	];
	
MatToRPY[R_] := MatToEulZYX[R];	

FramesToVect[list_] :=
	Module[{ R  = list[[1]],
			 Rd = list[[2]],
			 hn,  hs,    ha,
			  nd,  sd,    ad,
			 eO },
			
			 hn = Hat[R[[All, 1]]];    
			 hs = Hat[R[[All, 2]]];   
			 ha = Hat[R[[All, 3]]];
			
		   	 nd = Rd[[All, 1]]; 
		   	 sd = Rd[[All, 2]]; 
		   	 ad = Rd[[All, 3]];
			
			 eO = (1/2) (hn.nd + hs.sd + ha.ad)
	];

MatToVect[R_] := 
	Module[{axis, theta},
	{axis, theta} = RotationParam[R];
	axis*Sin[theta]
	]/; Length[R]==3  
	
MatToVect[list_]	:=
	Module[{ R  = list[[1]],
			 Rd = list[[2]],
			 hn,  hs,    ha,
			  nd,  sd,    ad,
			 eO },
			
			 hn = Hat[R[[All, 1]]];    
			 hs = Hat[R[[All, 2]]];   
			 ha = Hat[R[[All, 3]]];
			
		   	 nd = Rd[[All, 1]]; 
		   	 sd = Rd[[All, 2]]; 
		   	 ad = Rd[[All, 3]];
			
			 eO = (1/2) (hn.nd + hs.sd + ha.ad)
	]/; Length[list]==2
 
 (* Homogeneous Representations *)
 
HomogeneousToVector[p_] :=
	Block[{},
		Take[p, 3]
	];

(* Convert a point into homogeneous coordinates *)
PointToHomogeneous[p_] :=
  Block[{},
    (* Check to make sure the dimensions of the args make sense *)
    (* If[Not[VectorQ[p]], Message[Screws::notVector, "PointToHomogeneous"]]; *)

    (* Now put everything together into a homogeneous vector *)
    Append[p, 1]
  ];  

(* Convert a vector into homogeneous coordinates *)
VectorToHomogeneous[p_] :=
  Block[{},
    (* Check to make sure the dimensions of the args make sense *)
    (* If[Not[VectorQ[p]], Message[Screws::notVector, "VectorToHomogeneous"]]; *)

    (* Now put everything together into a homogeneous vector *)
    Append[p, 0]
  ];

(* Convert a rotation + a translation to a homogeneous matrix *)
RPToHomogeneous[R_, p_] :=
  Module[
    {n},

    (* Check to make sure the dimensions of the args make sense *)
    If[Not[VectorQ[p]] || Not[MatrixQ[R]] ||
       (n = Length[p]; Dimensions[R] != {n, n}),
	Message[Screws::wrongDimensions, "RPToHomogeneous:"];
    ];

    (* Now put everything together into a homogeneous transformation *)
    
    BlockMatrix[{{R,       Transpose[{p}]},
                 {ZeroMatrix[1,3],  {{1}}}
                 }
    ]
    
  ];  
  
(* Calculate the inverse rigid transformation *)
RigidInverse[g_?MatrixQ] := 
  Module[
    {R = RigidOrientation[g], p = RigidPosition[g]},
    RPToHomogeneous[Transpose[R], -Transpose[R].p]
  ];  

(* Extract the orientation portion from a homogeneous transformation *)
RigidOrientation[g_?MatrixQ]:=
  Module[
    {nr, nc},

    (* Check to make sure that we were passed a square matrix *)
    If[Not[MatrixQ[g]] || ({nr, nc} = Dimensions[g]; nr != nc) || nr < 3,
        Message[Screws::wrongDimensions, "RigidOrientation"];
	Return Null;
    ];

    (* Extract the 3x3 upper left corner *)
    SubMatrix[g, {1,1}, {nc-1,nc-1}]
  ];

(* Extract the orientation portion from a homogeneous transformation *)
RigidPosition[g_?MatrixQ]:=
  Module[
    {nr, nc},

    (* Check to make sure that we were passed a square matrix *)
    If[Not[MatrixQ[g]] || ({nr, nc} = Dimensions[g]; nr != nc) || nr < 3,
        Message[Screws::wrongDimensions, "RigidPosition"];
	Return Null;
    ];

    (* Extract the upper left column *)
    Flatten[SubMatrix[g, {1, nc}, {nc - 1 ,1}]]
  ];

(* Find the axis of a rotation matrix *)
RotationAxis[R_, theta_] :=
  Module[
    {nr, nc},

    (* Check to make sure that our input makes sense *)
    If[Not[MatrixQ[R]] || ({nr, nc} = Dimensions[R]; nr != nc),
        Message[Screws::wrongDimensions, "RotationAxis"];
	Return[Null];
    ];

    If[theta<0 || theta>Pi,
        Message[Screws::wrongTheta, "RotationAxis"];
	Return[Null];
    ];
 
    If[theta==Pi,
      axis=NullSpace[R-IdentityMatrix[3]][[1]];
      axis=axis/Magnitude[axis];
    ,
      axis={R[[3,2]]-R[[2,3]],R[[1,3]]-R[[3,1]],R[[2,1]]-R[[1,2]]}/(2*Sin[theta]);
    ];
    Return[axis];
];


RotationParam[R_] :=
  Module[
    {nr, nc},

    (* Check to make sure that our input makes sense *)
    If[Not[MatrixQ[R]] || ({nr, nc} = Dimensions[R]; nr != nc) || nr != 3,
        Message[Screws::wrongDimensions, "RotationAxis"];
	Return[Null];
    ];
     
     
    t = (Sum[R[[i,i]],{i, nr}]-1)/2;
    If[t<-1, t=-1;];
    If[t>1,t=1;];

    theta = ArcCos[t];

    If[theta != 0,
       axis=RotationAxis[R, theta];,
       axis = Table[0, {nc}];
       theta=0;
    ];

    Return[{axis, theta}]; 
 ];
 
 (* Lie Algebra Operations *) 
 
(* Generate a skew symmetric matrix from an axis*)
Hat[w_] := AxisToSkew[w];  (* synonym *) 

Skew[w_] := AxisToSkew[w]; (* synonym *)

AxisToSkew[omega_?VectorQ]:=
  Module[
    {},
    (* Check to make sure the dimensions are okay *)
    If[Not[VectorQ[omega]] || Length[omega] != 3,
      Message[Screws::wrongDimension];
      Return Null;
    ];

    (* Return the appropriate matrix *)
    {{ 0,          -omega[[3]],  omega[[2]]},
     { omega[[3]], 0,           -omega[[1]]},
     {-omega[[2]], omega[[1]],  0          }}
  ];

(* Generate an axis from a skew symmetric matrix *)
HatInv[S_] := SkewToAxis[S];  (* synonyms *)

UnSkew[S_] := SkewToAxis[S];  (* synonym *) 

Vee[S_] := SkewToAxis[S];     (* synonym *)

SkewToAxis[S_]:=
  Module[
    {},
    (* First check to make sure we have a skew symmetric matrix *)
    If[Not[skewQ[S]] || Dimensions[S] != {3,3},
      Message[Screws::wrongDimension];
      Return Null
    ];

    (* Now extract the appropriate component *)
    {S[[3,2]], S[[1,3]], S[[2,1]]}
  ];

(* Matrix exponential for a skew symmetric matrix or a vector *)

SkewExp[v_?VectorQ,theta_:1] := SkewExp[AxisToSkew[v],theta];

SkewExp[S_?skewQ,theta_:1]:=
  Module[
    {},
    (* Use Rodrigues's formula *)
    IdentityMatrix[3] + Sin[theta] S + (1 - Cos[theta]) S.S
  ];
  
(*
*  Rigid transformation matrices
*  Operations on se(3), the Lie algebra of rigid motions SE(3)
 *)    
 
 (* Figure out the dimension of a twist [private] *)
xidim[xi_?VectorQ] :=
  Module[
    {l = Length[xi], n},

    (* Check the dimensions of the vector to make sure everything is okay *)
    n = (Sqrt[1 + 8l] - 1)/2;
    If[Not[IntegerQ[n]],
      Message[Screws::wrongDimensions, "xidim"];
      Return 0;
    ];
    n
];

(* Extract the angular portion of a twist [private] *)
xitow[xi_?VectorQ] :=
  Module[
    {n = xidim[xi]},

    (* Make sure that the vector had a reasonable length *)   
    If[n == 0, Return Null];

    (* Extract the angular portion of the twist *)
    (* SetPrecision[Take[xi, -(n (n-1) / 2)],PRECISION] *)
    Take[xi, -(n (n-1) / 2)]
  ];

(* Extract the linear portion of a twist [private] *)
xitov[xi_?VectorQ] :=
  Module[
    {n = xidim[xi]},

    (* Make sure that the vector had a reasonable length *)   
    If[n == 0, Return Null];

    (* Extract the linear portion of the twist *)
    (* SetPrecision[Take[xi, n],PRECISION] *)
    Take[xi, n]
  ];

(* Check to see if a matrix is a twist matrix *)
(*! Not implemented !*)
TwistMatrixQ[A_] := MatrixQ[A];

(* Convert a homogeneous matrix to a twist *)
(*! This only works in dimensions 2 and 3 for now !*)
HomogeneousToTwist[A_] :=
  Module[
    {nr, nc},

    (* Check to make sure that our input makes sense *)
    If[Not[MatrixQ[A]] || ({nr, nc} = Dimensions[A]; nr != nc),
        Message[Screws::wrongDimensions, "HomogeneousToTwist"];
	Return Null;
    ];

    (* Make sure that we have a twist and not a rigid motion *)
    If[A[[nr,nc]] != 0,
        Message[Screws::notTwistMatrix, "HomogeneousToTwist"];
	Return Null;
    ];

    (* Extract the skew part and the vector part and make a vector *)
    Join[
      Flatten[SubMatrix[A, {1, nr}, {nr - 1, 1}]],
      SkewToAxis[ SubMatrix[A, {1, 1}, {nr - 1 ,nr - 1}] ]
    ]
  ];

(* Convert a twist to homogeneous coordinates *)
TwistToHomogeneous[xi_?VectorQ] :=
  Module[
    {w = xitow[xi], v = xitov[xi]},
    
    (* Make sure that we got a real twist *)
    If[w == Null || v == NULL, Return Null];

    (* Now put everything together into a homogeneous transformation *)
    BlockMatrix[{{AxisToSkew[w],   Transpose[{v}]},
                 {ZeroMatrix[1,3], {{0}} }
                 }
    ]
  ];  

(* Take the exponential of a twist *)
(*! This only works in dimension 3 for now !*)
TwistExp[xi_?MatrixQ, theta_:1]:=TwistExp[HomogeneousToTwist[xi], theta]; 

TwistExp[xi_?VectorQ, theta_:1] :=
  Module[
    {w = xitow[xi], ws, v = xitov[xi], R, p},
      
    (* Make sure that we got a real twist *)
    If[w == Null || v == NULL, Return Null];

    (* Use the exponential formula from MLS *)
    If [(MatchQ[w,{0,0,0}] || MatchQ[w, {{0},{0},{0}}]),
      R = IdentityMatrix[3];
      p = v * theta;,
     (* else *)
      ws=Skew[w];
      R = SkewExp[ws, theta];
      p = (IdentityMatrix[3] - R) . (ws . v) + w (w.v) theta;
    ];
    RPToHomogeneous[R, p]
  ];

(* Find the twist which generates a rigid motion - NEW VERSION *)
RigidTwist[g_?MatrixQ] :=
  Module[
    {R, p, v, theta, w, hatw, Ainv},

    (* Make sure the dimensions are okay *)
    (*! Missing !*)

    (* Extract the appropriate pieces of the homogeneous transformation *)
    R = RigidOrientation[g];
    p = RigidPosition[g];

    (* Now find the axis from the rotation *)    
    (*w = RotationAxis[R];
    *theta = RotationAngle[R];*)

    {w,theta}=RotationParam[R];
    hatw = Hat[w];
     
    (* Split into cases depending on whether theta is zero *)
    If[theta == 0,
      theta = Magnitude[p];
      If[theta == 0,  
        Return[{{0,0,0,0,0,0},0}];
        ];
      v = p/theta;,
    (* else *)
      (* Solve a linear equation to figure out what v is *)   
      Ainv = IdentityMatrix[3]/theta - (1/2) hatw + (1/theta - (1/2) Cot[theta/2]) MatrixPower[hatw, 2];
      
      v = Ainv.p;
    ];
    
	Return[{Flatten[{v, w}],theta}];
 
  ];

(*
 * Geometric attributes of twists and wrenches.
 *
 * For twists in R^3, find the attributes of that twist.
 *
 * Wrench attributes are defined by switching the role of linear
 * and angular portions
 *)

(* Build a twist from a screw *)
ScrewToTwist[Infinity, q_, w_] := Join[w, {0,0,0}];

ScrewToTwist[h_, q_, w_] := Join[-AxisToSkew[w] . q + h w, w]

(* Find the pitch associated with a twist in R^3 *)
TwistPitch[xi_?VectorQ] := 
  Module[
    {v, w},
    {v, w} = Partition[xi, 3];
    v . w / w.w
  ];
  
WrenchPitch[xi_?VectorQ] := Null;

(* Find the axis of a twist *)
TwistAxis[xi_?VectorQ] := 
  Module[
    {v, w},
    {v, w} = Partition[xi, 3];
    If[(MatchQ[w,{0,0,0}] || MatchQ[w, {{0},{0},{0}}]), 
     {0, v / Sqrt[v.v]}, {AxisToSkew[w] . v / w.w, (w / w.w)}]
  ];

WrenchAxis[xi_?VectorQ] := Null;

(* Find the magnitude of a twist *)
TwistMagnitude[xi_?VectorQ] := 
  Module[
    {v, w},
    {v, w} = Partition[xi, 3];
    If[(MatchQ[w,{0,0,0}] || MatchQ[w, {{0},{0},{0}}]), 
      Sqrt[v.v], Sqrt[w.w]]
  ];
WrenchMagnitude[xi_?VectorQ] := Null;
 

(*
 * Adjoint calculation
 *
 * The adjoint matrix maps twist vectors to twist vectors.
 *
 *)

(* Adjoint matrix calculation *)
RigidAdjoint[g_?MatrixQ] := 
  Module[
    {R = RigidOrientation[g], p = RigidPosition[g]},
    
    BlockMatrix[{{R,                AxisToSkew[p] . R},
                 {ZeroMatrix[3,3],           R       }
                 }
    ]
    
  ];

Ad[g_?MatrixQ] := RigidAdjoint[g];
  
(* Inverse adjoint matrix calculation *)
InverseRigidAdjoint[g_?MatrixQ] := 
  Module[
    {RT = Transpose[RigidOrientation[g]], p = RigidPosition[g]},
    
    BlockMatrix[{{RT,       -RT.AxisToSkew[p]},
                 {ZeroMatrix[3,3],  RT}
                 }
    ]
    
  ];
  
AdInv[g_?MatrixQ] := InverseRigidAdjoint[g];    

TransposeRigidAdjoint[g_?MatrixQ] :=
	Module[
	  {RT = Transpose[RigidOrientation[g]], p = RigidPosition[g]},
    
       BlockMatrix[{
       	         {               RT,        ZeroMatrix[3,3]},
                 {-RT.AxisToSkew[p],        RT             }
                 }
       ]
       
	];

AdTr[g_?MatrixQ] := TransposeRigidAdjoint[g];

InverseTransposeRigidAdjoint[g_?MatrixQ] := 
	Module[
		{R = RigidOrientation[g], p = RigidPosition[g]},
    
         BlockMatrix[{
       	         {               R,        ZeroMatrix[3,3]},
                 { AxisToSkew[p].R,        R              }
                 }
       	 ]
				
	];
	
AdStar[g_?MatrixQ] := InverseTransposeRigidAdjoint[g];	

(* adjoint (small) matrix calculations - aka Lie Bracket *)

ad[xi_?VectorQ] := 
  Module[
    {v, w, vs, ws},
    {v, w} = Partition[xi, 3];
     vs = Hat[v];
     ws = Hat[w];
    
    BlockMatrix[{
       	         {              ws,        vs},
                 { ZeroMatrix[3,3],        ws}
                 }
    ]
    
        ];

adTr[xi_?VectorQ] := 
  Module[
    {v, w, vs, ws},
    {v, w} = Partition[xi, 3];
     vs = Hat[v];
     ws = Hat[w];
    
    BlockMatrix[{
       	         { -ws,        ZeroMatrix[3,3] },
                 { -vs,                    -ws }
                 }
    ]
    
        ];

adStar[xi_?VectorQ] := 
  Module[
    {v, w, vs, ws},
    {v, w} = Partition[xi, 3];
     vs = Hat[v];
     ws = Hat[w];
    
    BlockMatrix[{
       	         {  ws,        ZeroMatrix[3,3] },
                 {  vs,                     ws }
                 }
    ]
    
        ];

LieBracketDull[xia_?VectorQ, xib_?VectorQ] :=
	Module[
		{xias, xibs, xics},
		
		xias = TwistToHomogeneous[xia]; 
		xibs = TwistToHomogeneous[xib];
		
		xics = xias.xibs - xibs.xias;
		
		HomogeneousToTwist[xics]
	];
	
LieBracket[xia_?VectorQ, xib_?VectorQ] :=
	Module[
		{adxia},
		
		adxia = ad[xia];
		
		adxia.xib
	];	

LieBracket[xia_?MatrixQ, xib_?MatrixQ] := LieBracket[HomogeneousToTwist[xia], HomogeneousToTwist[xib]];
  
 (* Calculation of the error twist xi_err and error angle theta_err such that:
    gStart.TwistExp[xi_err, theta_err] = gEnd *)
    
 LocalTranTwist[gStart_?MatrixQ, gEnd_?MatrixQ]:=
 Module[
     {gError, xi, theta, i},
     
     gError = Inverse[gStart].gEnd;
     
     If[gError == IdentityMatrix[4],
       Return[{Table[0,{i,6}],0}];
     ];
     
     {xi,theta} = RigidTwist[gError];
     
     Return[{xi,theta}];
];

(* Calculation of the error twist xi_err and error angle theta_err such that:
    TwistExp[xi_err, theta_err].gStart = gEnd *)
  
 GlobalTranTwist[gStart_?MatrixQ, gEnd_?MatrixQ]:=
 Module[
     {gError, xi, theta, Adg, i},
     
     gError = Inverse[gStart].gEnd;
     
     If[gError == IdentityMatrix[4],
       Return[{Table[0,{i,6}],0}];
     ];
     
     {xi,theta} = RigidTwist[gError];
     
     Adg = RigidAdjoint[gStart];
      xi = Adg.xi;
       
     Return[{xi,theta}];
]; 
 
 (* Systematic Methods to Build Forward Kinematics and Spatial/Body Jacobians *)
 
(* Gives Xi 6 vector given a point on axis and axis unit vector for a Revolute Joint *)
RevoluteTwist[q_, w_]:= Flatten[ { Cross[q, w], w } ];

(* Gives Xi 6 vector given a point on axis and axis unit vector for a Prismatic Joint *)
PrismaticTwist[q_, w_]:= Flatten[ { w, {0,0,0} } ];

(* Gives the homogeneous matrix *)
ForwardKinematics[args__, gst0_]:= 
  Module[
    { g, i,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* Initialize the transformation matrix *)
    g = TwistExp[argList[[1,1]], argList[[1,2]]];   

    (* Build up the transformation matrix joint by joint *)
    For[i = 2, i <= n, i++,
      (* Update the transformation matrix *)
      g = g . TwistExp[ argList[[i,1]], argList[[i,2]] ];
    ];      

    (* Finish by multiplying by the initial tool configuration *)
    g . gst0
  ];

ForwardKinematicsLeft[gst0_, args__]:=
 Module[{g, i, argList={args}, n=Length[{args}]},
g=gst0.TwistExp[argList[[1,1]],argList[[1,2]]];
For[i=2,i<=n,i++,
g=g.TwistExp[argList[[i,1]],argList[[i,2]]];
];
g
];
			
(* Gives the homogeneous matrix via the local POE - R_i = M_i exp(X_i q_i) *)

ForwardKinematicsLocalPOE[args__]:=
	Module[
		{ g, i,
		  argList = {args},
		  n = Length[{args}]
		},
		
		(* Initialize the transformation matrix *)
		g = argList[[ 1, 1 ]] . TwistExp[ argList[[ 1, 2 ]], argList[[ 1, 3]] ];
		
		(* Build up the transformation matrix joint by joint *)
		For[ i = 2, i <= n, i++,
			
			g = g . argList[[ i, 1 ]] . TwistExp[ argList[[ i, 2 ]], argList[[ i, 3 ]] ];
			
		];
		
		(* Return g *)
		
		g
	];

ForwardKinematicsLocalPOERight[args__]:=
Module[{g,i,argList={args},n=Length[{args}]},
g=TwistExp[argList[[1,2]],argList[[1,3]]].argList[[1,1]];
For[i=2,i<=n,i++,
g=g.TwistExp[argList[[i,2]],argList[[i,3]]].argList[[i,1]];
];
g
];

ForwardKinematicsLocalPOETranslated[args__]:=
	Module[
		{ M, Y, g, i,
		  argList = {args},
		  n = Length[{args}]
		},
		
		(* Translation of the quantities performed such that
		   
		   Yi    = Ad[M1 M2 ... Mi] Xi
		   g(0)  = M1 M2 ... Mn
		   
		   g = exp(Y1 q1) exp(Y2 q2) ... exp(Yn qn) g(0)
		
		 *)
		
		(* Initialize the transformation matrix *)
		M = argList[[ 1, 1 ]];
		Y = Ad[M] . argList[[ 1, 2 ]];
		g = TwistExp[ Y , argList[[ 1, 3]] ];
		
		(* Build up the transformation matrix joint by joint *)
		For[ i = 2, i <= n, i++,
			M = M . argList[[ i, 1 ]];
			Y = Ad[M] . argList[[ i, 2 ]];
			g = g . TwistExp[ Y, argList[[ i, 3 ]] ];
			
		];
		
		(* Return g *)
		
		g.M
	];

ForwardKinematicsDH[args__] :=
 Module[{g,i,argList={args},n=Length[{args}]},
(*Initialize the transformation matrix*)
g=RPToHomogeneous[RotZ[argList[[1,1]]],{0,0,argList[[1,2]]}].RPToHomogeneous[RotX[argList[[1,3]]],{argList[[1,4]],0,0}];
For[i=2,i<=n,i++,
g=g.RPToHomogeneous[RotZ[argList[[i,1]]],{0,0,argList[[i,2]]}].RPToHomogeneous[RotX[argList[[i,3]]],{argList[[i,4]],0,0}];
];
(*Return the final transformation*)
g
];

(* Construct the Spatial Jacobian for a robot with any no. of links *)
SpatialJacobianGlobalPOE[args__, gst0_] := 
  Module[
    { i, xi, Js, g,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* First initialize the Jacobian and compute the first column *)
    Js = {argList[[1,1]]};
    g = TwistExp[argList[[1,1]], argList[[1,2]]];   

    (* Build up the Jacobian joint by joint *)
    For[i = 2, i <= n, i++,
      (* Compute this column of the Jacobian and append it to Js *)
      xi = RigidAdjoint[g] . argList[[i,1]];
      Js = Join[Js, {xi}];      

      (* Update the transformation matrix *)
      g = g . TwistExp[argList[[i,1]], argList[[i,2]]];
    ];      

    (* Return the Jacobian *)
    Transpose[Js]
  ];
			
(* Construct the Body Jacobian for a robot with any no. of links *)			
BodyJacobianGlobalPOE[args__, gst0_] := 
  Module[
    { i, xi, Jb, g,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* Initialize the Jacobian and the transformation matrix *)
    Jb = {};
    g = gst0;

    (* Build up the Jacobian joint by joint *)
    For[i = n, i >= 1, i--,
      (* Compute this column of the Jacobian and prepend it to Jb *)
      xi = RigidAdjoint[RigidInverse[g]] . argList[[i,1]];
      Jb = Join[{xi}, Jb];      

      (* Update the transformation matrix *)
      g = TwistExp[argList[[i,1]], argList[[i,2]]] . g;
    ];      

    (* Return the Jacobian *)
    Transpose[Jb]
  ];

BodyJacobianGlobalPOELeft[gst0_, args__] :=
Module[ {i, xi, g, Jb, argList={args}, n=Length[{args}]},
Jb={};
g=IdentityMatrix[4];
For[i = n, i >= 1, i--,
g=TwistExp[argList[[i,1]],argList[[i,2]]].g;
xi=AdInv[g].argList[[i,1]];
Jb=Join[{xi},Jb];
];
Transpose[Jb]
];

SpatialJacobianGlobalPOELeft[gst0_, args__] :=
Module[{i, xi, g, Js, argList={args}, n=Length[{args}]},
g=gst0;
Js={};
For[i=1, i<=n, i++,
g=g.TwistExp[argList[[i,1]],argList[[i,2]]];
xi=Ad[g].argList[[i,1]];
Js=Join[Js, {xi}];
];
Transpose[Js]
];

SpatialJacobianLocalPOERight[args__] :=
Module[{i, xi, g, Js, argList={args}, n=Length[{args}]},
g=IdentityMatrix[4];
Js={};
For[i=1,i<=n,i++,
xi=Ad[g].argList[[i,2]];
Js=Join[Js,{xi}];
g=g.TwistExp[argList[[i,2]],argList[[i,3]]].argList[[i,1]];
];
Transpose[Js]
];

BodyJacobianLocalPOERight[args__] :=
Module[{i, xi, g, Jb, argList={args}, n=Length[{args}]},
g=IdentityMatrix[4];
Jb={};
For[i=n,i>=1,i--,
g=TwistExp[argList[[i,2]],argList[[i,3]]].argList[[i,1]].g;
xi=AdInv[g].argList[[i,2]];
Jb=Join[{xi},Jb];
];
Transpose[Jb]
];

(* Construct the Spatial Jacobian for a robot with any no. of links *)
SpatialJacobianLocalPOE[args__] :=
	Module[
	{ j, xj, Js, g,
      argList = {args},		(* turn arguments into a real list *)
      i = Length[{args}]	(* decide on the number of joints *)
    },
    
    (* First initialize the Jacobian and compute the first column *)
    g  = Eye[4]; 
    Js = {};  

    (* Build up the Jacobian joint by joint *)
    For[j = 1, j <= i, j++,
      
      (* Update the transformation matrix *)
      g = g . argList[[j,1]] . TwistExp[ argList[[j,2]], argList[[j,3]] ];
      
      (* Compute this column of the Jacobian and append it to Js *)
      xj = Ad[g] . argList[[j,2]];
      Js = Join[Js, {xj}];      

    ];      

    (* Return the Jacobian *)
    Transpose[Js]
    
	];

  
(* Construct the Body Jacobian for a robot with any no. of links *)	
BodyJacobianLocalPOE[args__] := 
  Module[
    { j, xj, Jb, g,
      argList = {args},		(* turn arguments into a real list *)
      i = Length[{args}]	(* decide on the number of joints *)
    },

    (* Initialize the Jacobian and the transformation matrix *)
    Jb = {};
     g = Eye[4];

    (* Build up the Jacobian joint by joint *)
    For[j = i, j >= 1, j--,
      (* Compute this column of the Jacobian and prepend it to Jb *)
      xj = AdInv[ g ] . argList[[ j,2 ]];
      Jb = Join[ {xj}, Jb ];      

      (* Update the transformation matrix *)
      g = argList[[j,1]] . TwistExp[ argList[[j,2]], argList[[j,3]] ] . g;
    ];      

    (* Return the Jacobian *)
    Transpose[Jb]
  ];   
 
(* Given a list of twists in the reference global configuration we compute the
   derivative of the twist Vab,t
   {X1, q1, q1p, q1pp}, ..., {Xn, qn, qnp, qnpp}, gst0
    *)
SpatialTwistDerivative[args__, gst0_] :=
	Module[
    { i, k, Xir, Xi, Xkr, Xk, qi, qip, qipp, qk, qkp, gi, gk, Vabp,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* First initialize the quantities and the transformations *)
    Vabp = argList[[1,1]] argList[[1,4]];
    gi   = TwistExp[ argList[[1,1]], argList[[1,2]] ];
       

    (* Build up Vabp joint by joint and by joint dependencies *)
    For[i = 2, i <= n, i++,

	  gk   = Eye[4];
      (* Compute the i-th column of the Jacobian *)
      Xir  = argList[[i,1]];
      Xi   = Ad[gi] . Xir;
      qi   = argList[[i,2]];
      qip  = argList[[i,3]];
      qipp = argList[[i,4]];
      Vabp = Vabp + Xi qipp;
      For[ k = 1, k <= i-1, k++,
               
              Xkr = argList[[k,1]];       	
      	    Xk  = Ad[gk] . Xkr;
      	    qk  = argList[[k,2]];
      	    qkp = argList[[k,3]];
      	      
            Vabp = Vabp +  ad[Xk].Xi qkp qip;
            
           (* Update the transformation matrix in the inner cycle *)
           gk = gk . TwistExp[ Xkr, qk ]; 
 	
          ]; 	

      (* Update the transformation matrix in the outer cycle *)
      gi = gi . TwistExp[ Xir, qi ];
    ];      

    (* Return Vabp *)
       Vabp
	];

SpatialTwistDerivativeLocal[args__] :=
Module[
	{ i, k, Xir, Xi, Xkr, Xk, qi, qip, qipp, qk, qkp, gi, gk, Vabp,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* First initialize the quantities and the transformations *)
    gi   = argList[[1,1]] . TwistExp[ argList[[1,2]], argList[[1,3]] ];
	Vabp = Ad[gi] . argList[[1,2]] argList[[1,5]];
       

    (* Build up Vabp joint by joint and by joint dependencies *)
    For[i = 2, i <= n, i++,

	  gk   = Eye[4];
      (* Compute the i-th column of the Jacobian *)
      Xir  = argList[[i,2]];
      qi   = argList[[i,3]];
      qip  = argList[[i,4]];
      qipp = argList[[i,5]];
	  (* Update the transformation matrix in the outer cycle *)
      gi = gi . argList[[i,1]] . TwistExp[ Xir, qi ];
	  Xi   = Ad[gi] . Xir;
      Vabp = Vabp + Xi qipp;
      For[ k = 1, k <= i-1, k++,
               
              Xkr = argList[[k,2]];      	
      	    qk  = argList[[k,3]];
      	    qkp = argList[[k,4]];
			  (* Update the transformation matrix in the inner cycle *)
              gk = gk . argList[[k,1]] . TwistExp[ Xkr, qk ]; 
			  Xk  = Ad[gk] . Xkr;
      	      
            Vabp = Vabp +  ad[Xk].Xi qkp qip;
            
            	
          ]; 	
      
    ];      

    (* Return Vabp *)
       Vabp
	];

BodyTwistDerivative[args__,gst0_]:=
Module[
	{ i, k, Xir, Xi, Xkr, Xk, qi, qip, qipp, qk, qkp, gi, gk, Vabp,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* First initialize the quantities and the transformations *)
    gi   = gst0;
	Vabp = {0,0,0,0,0,0};
       

    (* Build up Vabp joint by joint and by joint dependencies *)
    For[i = n, i >= 1, i--,

	  gk   = gst0;
      (* Compute the i-th column of the Jacobian *)
      Xir  = argList[[i,1]];
      qi   = argList[[i,2]];
      qip  = argList[[i,3]];
      qipp = argList[[i,4]];
	  Xi   = AdInv[gi] . Xir;
      Vabp = Vabp + Xi qipp;
      For[ k = n, k > i, k--,
               
              Xkr = argList[[k,1]];      	
      	    qk  = argList[[k,2]];
      	    qkp = argList[[k,3]];
			  Xk  = AdInv[gk] . Xkr; 

			  Vabp = Vabp -  ad[Xk].Xi qkp qip;
			  (* Update the transformation matrix in the inner cycle *)
              gk = TwistExp[ Xkr, qk ] . gk; 
          
            	
          ]; 
	  (* Update the transformation matrix in the outer cycle *)
      gi = TwistExp[ Xir, qi ] . gi;	
      
    ];      

    (* Return Vabp *)
       Vabp
	];

BodyTwistDerivativeLocal[args__]:=
Module[
	{ i, k, Xir, Xi, Xkr, Xk, qi, qip, qipp, qk, qkp, gi, gk, Vabp,
      argList = {args},		(* turn arguments into a real list *)
      n = Length[{args}]	(* decide on the number of joints *)
    },

    (* First initialize the quantities and the transformations *)
    gi   = Eye[4];
	Vabp = {0,0,0,0,0,0};

    (* Build up Vabp joint by joint and by joint dependencies *)
    For[i = n, i >= 1, i--,

	  gk   = Eye[4];
      (* Compute the i-th column of the Jacobian *)
      Xir  = argList[[i,2]];
      qi   = argList[[i,3]];
      qip  = argList[[i,4]];
      qipp = argList[[i,5]];
	  Xi   = AdInv[gi] . Xir;
      Vabp = Vabp + Xi qipp;
      For[ k = n, k > i, k--,
               
              Xkr = argList[[k,2]];      	
      	    qk  = argList[[k,3]];
      	    qkp = argList[[k,4]];
			  Xk  = AdInv[gk] . Xkr; 

			  Vabp = Vabp -  ad[Xk].Xi qkp qip;
			  (* Update the transformation matrix in the inner cycle *)
              gk = argList[[k,1]] . TwistExp[ Xkr, qk ] . gk; 
          
            	
          ]; 

    
	  (* Update the transformation matrix in the outer cycle *)
      gi = argList[[i,1]] . TwistExp[ Xir, qi ] . gi;	
      
    ];      

    (* Return Vabp *)
       Vabp
	];

FBSpatialJacobianGlobalPOE[argList_, gst0_, JointChainList_, Nj_]:=
Module[{g, Js, Ji, cond, Yi, qi, i, n},
	n = Length[argList];
	g = RPToHomogeneous[QuatToMat[argList[[1,2]]],argList[[1,1]]];
	Js = Transpose[BlockMatrix[{{Eye[3], Hat[argList[[1,1]]].QuatToSpatialJac[argList[[1,2]]]},{ZeroMatrix[3,3], QuatToSpatialJac[argList[[1,2]]]}}]];
	For[i = 2, i <= n, i++,
		Yi = argList[[i,1]];
		qi = argList[[i,2]];
		cond = Dimensions[Select[JointChainList, # = i&]];
		If[cond[[1]] == 1, g = g . TwistExp[Yi,qi], g = g];
		If[cond[[1]] == 1, Ji = Ad[g].Yi, Ji = First[ZeroMatrix[1,6]]];
		Js = Join[Js, {Ji}];
	];
	For[i = n + 1, i <= Nj, i++,
		Ji = First[ZeroMatrix[1,6]];
		Js = Join[Js, {Ji}];
	];

	(*Return the Spatial Jacobian*)
	Transpose[Js]
];

FBSpatialJacobianLocalPOE[argList_, JointChainList_, Nj_]:=
Module[{g, Js, Ji, cond, Xi, qi, i, n},
	n = Length[argList];
	g = RPToHomogeneous[QuatToMat[argList[[1,2]]],argList[[1,1]]];
	Js = Transpose[BlockMatrix[{{Eye[3], Hat[argList[[1,1]]].QuatToSpatialJac[argList[[1,2]]]},{ZeroMatrix[3,3], QuatToSpatialJac[argList[[1,2]]]}}]];
	For[i = 2, i <= n, i++,
		Xi = argList[[i,2]];
		qi = argList[[i,3]];
		cond = Dimensions[Select[JointChainList, # = i&]];
		If[cond[[1]] == 1, g = g . argList[[i,1]] . TwistExp[Xi,qi], g = g];
		If[cond[[1]] == 1, Ji = Ad[g].Xi, Ji = First[ZeroMatrix[1,6]]];
		Js = Join[Js, {Ji}];
	];
	For[i = n + 1, i <= Nj, i++,
		Ji = First[ZeroMatrix[1,6]];
		Js = Join[Js, {Ji}];
	];

	(*Return the Spatial Jacobian*)
	Transpose[Js]
];

FBBodyJacobianGlobalPOE[argList_, gst0_, JointChainList_, Nj_]:=
Module[{g, Jb, Ji, Yi, qi, cond1, cond2, i, n},
	n = Length[argList];
	g = gst0;
	Jb = {};
	For[i = Nj, i > n, i--,
		Ji = First[ZeroMatrix[1,6]];
		Jb = Join[{Ji}, Jb];
	];
	For[i = n, i >= 2, i--,
		Yi = argList[[i,1]];
		qi = argList[[i,2]];
		cond1 = Dimensions[Select[JointChainList, # = i&]];
		If[cond1[[1]] == 1, Ji = AdInv[g].Yi, Ji = First[ZeroMatrix[1,6]]];
		If[cond1[[1]] == 1, g = TwistExp[Yi,qi] . g, g = g];
		Jb = Join[{Ji}, Jb];
	];
	
	cond2 = Dimensions[Jb];
	If[cond2[[1]] != 0, Jb = Transpose[Jb], Jb = {}];
	Ji = AdInv[g].BlockMatrix[{{Transpose[QuatToMat[argList[[1,2]]]],ZeroMatrix[3,4]},{ZeroMatrix[3,3],QuatToBodyJac[argList[[1,2]]]}}];

	(*Return the Body Jacobian*)
	Jb = Join[{Ji}, Jb]
];

FBBodyJacobianLocalPOE[argList_, JointChainList_, Nj_]:=
Module[{g, Jb, Ji, Xi, qi, cond1, cond2, i, n},
	n = Length[argList];
	g = Eye[4];
	Jb = {};
	For[i = Nj, i > n, i--,
		Ji = First[ZeroMatrix[1,6]];
		Jb = Join[{Ji}, Jb];
	];
	For[i = n, i >= 2, i--,
		Xi = argList[[i,2]];
		qi = argList[[i,3]];
		cond1 = Dimensions[Select[JointChainList, # = i&]];
		If[cond1[[1]] == 1, Ji = AdInv[g].Xi, Ji = First[ZeroMatrix[1,6]]];
		If[cond1[[1]] == 1, g = argList[[i,1]].TwistExp[Xi,qi].g, g = g];
		Jb = Join[{Ji}, Jb];
	];
	
	cond2 = Dimensions[Jb];
	If[cond2[[1]] != 0, Jb = Transpose[Jb], Jb = {}];
	Ji = AdInv[g].BlockMatrix[{{Transpose[QuatToMat[argList[[1,2]]]],ZeroMatrix[3,4]},{ZeroMatrix[3,3],QuatToBodyJac[argList[[1,2]]]}}];

	(*Return the Body Jacobian*)
	Jb = Join[{Ji}, Jb]
];

SpatialJacobianDerivative[Js_,qdotList_]:=
Module[{J, i, Ji, Jdoti, k, Jk, Jdotk, Jdot, n},
	n = Length[qdotList];
	Jdot = ZeroMatrix[1,6];
	J = Transpose[Js];
	For[ i = 2, i <= n, i++,
		Ji = J[[i]];
		Jdotk = ZeroMatrix[6,6];
		For[ k = 1, k < i, k++,
			Jk = J[[k]];
			Jdotk = Jdotk + qdotList[[k]] ad[Jk];
		];
		Jdoti = Jdotk . Ji;
		Jdot = Join[Jdot, {Jdoti}];
	];
	(*Return the derivative*)
	Transpose[Jdot]

];

BodyJacobianDerivative[Jb_,qdotList_]:=
Module[{J, i, Ji, Jdoti, k, Jk, Jdotk, Jdot, n},
	n = Length[qdotList];
	Jdot = ZeroMatrix[1,6];
	J = Transpose[Jb];
	For[ i = n - 1, i >= 1, i--,
		Ji = J[[i]];
		Jdotk = ZeroMatrix[6,6];
		For[ k = i + 1, k <= n, k++,
			Jk = J[[k]];
			Jdotk = Jdotk - qdotList[[k]] ad[Jk];
		];
		Jdoti = Jdotk . Ji;
		Jdot = Join[{Jdoti}, Jdot];
	];
	(*Return the derivative*)
	Transpose[Jdot]

];

FBSpatialJacobianDerivative[FBJs_,qdotList_]:=
Module[{J, i, Ji, Jdoti, k, Jk, Jdotk, Jdot, n},
	n = Length[qdotList];
	Jdot = ZeroMatrix[7,6];
	J = Transpose[FBJs];
	For[ i = 8, i <= n, i++,
		Ji = J[[i]];
		Jdotk = ad[Transpose[Take[J,7]].Take[qdotList,7]];
		For[ k = 8, k < i, k++,
			Jk = J[[k]];
			Jdotk = Jdotk + qdotList[[k]] ad[Jk];
		];
		Jdoti = Jdotk . Ji;
		Jdot = Join[Jdot, {Jdoti}];
	];
	(*Return the derivative*)
	Transpose[Jdot]

];

FBBodyJacobianDerivative[FBJb_,qdotList_]:=
Module[{J, i, Ji, Jdoti, k, Jk, Jdotk, Jdot1, Jdot, n},
	n = Length[qdotList];
	Jdot = ZeroMatrix[1,6];
	J = Transpose[FBJb];
	For[ i = n - 1, i > 7, i--,
		Ji = J[[i]];
		Jdotk = ZeroMatrix[6,6];
		For[ k = i + 1, k <= n, k++,
			Jk = J[[k]];
			Jdotk = Jdotk - qdotList[[k]] ad[Jk];
		];
		Jdoti = Jdotk . Ji;
		Jdot = Join[{Jdoti}, Jdot];
	];
	Jdotk = ZeroMatrix[6,6];
	For[k = 8, k <= n, k++,
		Jk = J[[k]];
		Jdotk = Jdotk - qdotList[[k]] ad[Jk];
	];
	Jdot1 = Transpose[Jdotk . Transpose[Take[J],7]];
	Jdot = Join[{Jdot1}, Jdot];
	(*Return the derivative*)
	Transpose[Jdot]

];

 
End[];
 
EndPackage[];


