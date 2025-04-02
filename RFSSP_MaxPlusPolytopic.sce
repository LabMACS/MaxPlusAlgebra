// Robust Feedback System Synchronization Problem (RFSSP)
// RFSSP for max-plus linear system with polytopic uncertainties

// To be executed with ScicosLab

//Remember: %1 is zero, and %0 is -infinity
clear;
clc;

// To compute the difference of semimodules (ominus)
function Z=ominus(A, B)
    A = full(A);
    B = full(B);
    sA = size(A);
    sB = size(B);
    I = %eye(sA(1), sA(1));
    SX = full(#([I B %zeros(sA(1), sA(2))]));
    DX = full(#([%zeros(sA(1), sA(1) + sB(2)) A]));
    X = mpsolve(SX, DX);
    Z = weakbasis(X(1:sA(1),:));
endfunction


function Z=pullback(A,B)
    A = full(A);
    B = full(B);
    sA = size(A);
    sB = size(B);
    SX = full(#([A %zeros(sB(1), sB(2))]));
    DX = full(#([%zeros(sA(1), sA(2)) B]));
    X = mpsolve(SX, DX);
    Z = weakbasis(X);
endfunction

function Z=invimage(A, V)
     A = full(A);
     V = full(V);
     sA = size(A);
     sV = size(V);
     SX = full(#([A %zeros(sV(1), sV(2))]));
     DX = full(#([%zeros(sA(1), sA(2)) V]));
     X = mpsolve(SX, DX);
     Z = weakbasis(X(1:sA(2),:));
endfunction

function Z=intersection(A,B)
    sA = size(full(A));
    sB = size(full(B));
    SX = full(#([A %zeros(sB(1), sB(2))]));
    DX = full(#([%zeros(sA(1), sA(2)) B]));
    X = mpsolve(SX, DX);
    Z = weakbasis(full(SX * X));
endfunction

//Code for vertical concatenation of matrices given as input (in place of varargin)
function Z=co_i(varargin)
    Z = varargin(1);
    for i = 2:argn(2)
        Z = [Z; varargin(i)];
    end
endfunction

///////////////////////////////////// 


// Plant matrices
disp("Plant Matrices with Polytopic Uncertainties")
disp(ascii(10))


AP1=#([%1 %1; 3 %1])
disp(ascii(10))
AP2=#([2 %1; 4 %1])
disp(ascii(10))

BP1=#([3; %1])  
BP2=#([3; 1])

disp(ascii(10))

CP1=full(#([1 %1]))  
CP2=full(#([1 1]))

disp(ascii(10))


///////////////////////////////////// 

// Now the RFSSP starts.


// Model matrices
disp("Model Matrices")
disp(ascii(10))
AM = #([4])
disp(ascii(10))
BM = #([6])
disp(ascii(10))
CM = #([1])
disp(ascii(10))

[np, mp] = size(full(BP1));
[nm, mm] = size(full(BM));


// Joint dynamics

AE1 = #([AP1 %zeros(np, nm); %zeros(nm, np) AM]);
AE2 = #([AP2 %zeros(np, nm); %zeros(nm, np) AM]);

B1_1 = #([BP1; %zeros(nm, mp)]);
B1_2 = #([BP2; %zeros(nm, mp)]);

B2_1 = #([%zeros(np, mm); BM]);
B2_2 = #([%zeros(np, mm); BM]);


AE_co=co_i(AE1,AE2);
B1_co=co_i(B1_1,B1_2);
B2_co=co_i(B2_1,B2_2);
CP_co=co_i(CP1,CP2);
CM_co=co_i(CM,CM);


// Computations


K = pullback(CP_co, CM_co);

[nK, mK] = size(full(K));

Vk=#([K]);
VkN = #([K %zeros(nK,mK);%zeros(nK,mK) K]);
k = 0;
Vstar_Found = %F;
while Vstar_Found == %F & k < 5
    disp(ascii(10))
    k = k + 1
    disp(ascii(10))
    Vkold = #([Vk]);
    [nVkold, mVkold] = size(full(Vkold));
    VkoldN = #([Vk %zeros(nVkold,mVkold);%zeros(nVkold,mVkold) Vkold]);
    Vk = intersection(Vkold, invimage(AE_co, ominus(VkoldN, B1_co))); 
    Vstar_Found = equalspan(Vkold, Vk) // true -> V* found
        disp("Number of iterations: " + string(k+1))
        disp(ascii(10))
end

Vstar = Vk
[nt, v]=size(Vstar);
VstarN = #([Vstar %zeros(nt,v);%zeros(nt,v) Vstar]);

C1 = includespan(ominus(VstarN, B1_co), full(B2_co)); //Im B2 inside Vstar-ImB1
disp(ascii(10))
Solvable = C1


////// Verifying given matrices F e G /////////////

if equalspan(K, Vstar) == %T then
    disp("Vstar is equivalent to K, since they generate the same space, so:")
    Vstar=K
else
    disp("Vstar is contained in K")
    Vstar=Vstar
  end
  
  
///////////////////////////////////// 

// Inserting matrix F
disp(ascii(10))

F=#([1 1 %1])

disp(ascii(10))
AE1Vstar=AE1*Vstar
disp(ascii(10))
B11FVstar=B1_1*F*Vstar
disp(ascii(10))
AE1Vstar+B11FVstar;

if includespan(Vstar, AE1Vstar+B11FVstar) == %T then
    disp("(AE1+B11F)Vstar is contained in Vstar.")
else
    disp("(AE1+B11F)Vstar is NOT contained in Vstar.")
end

disp(ascii(10))
AE2Vstar=AE2*Vstar
disp(ascii(10))
B11FVstar=B1_1*F*Vstar
disp(ascii(10))
AE2Vstar+B11FVstar;


if includespan(K, AE2Vstar+B11FVstar) == %T then
    disp("(AE2+B11F)Vstar is contained in Vstar.")
else
    disp("(AE2+B11F)Vstar is NOT contained in Vstar.")
end

////////// 
disp(ascii(10))
AE1Vstar=AE1*Vstar
disp(ascii(10))
B12FVstar=B1_2*F*Vstar
disp(ascii(10))
AE1Vstar+B12FVstar;

if includespan(Vstar, AE1Vstar+B11FVstar) == %T then
    disp("(AE1+B12F)Vstar is contained in Vstar.")
else
    disp("(AE1+B12F)Vstar is NOT contained in Vstar.")
end

disp(ascii(10))
AE2Vstar=AE2*Vstar
disp(ascii(10))
B12FVstar=B1_2*F*Vstar
disp(ascii(10))
AE2Vstar+B12FVstar;


if includespan(K, AE2Vstar+B12FVstar) == %T then
    disp("(AE2+B12F)Vstar is contained in Vstar.")
else
    disp("(AE2+B12F)Vstar is NOT contained in Vstar.")
end


///////////////////////////////////// 

// Inserting matrix G
disp(ascii(10))
G=#([3])

Bp1GBM=#([BP1*G; BM]);
disp("Matrix [BP1*G; BM]");
disp(Bp1GBM)


if includespan(K, Bp1GBM) == %T then
    disp("The columns of [BP1*G; BM] are contained in Vstar.")
else
    disp("The columns of [BP1*G; BM] are NOT contained in Vstar")
end
/////
Bp2GBM=#([BP2*G; BM]);
disp("Matrix [BP2*G; BM]");
disp(Bp2GBM)


if includespan(K, Bp2GBM) == %T then
    disp("The columns of [BP2*G; BM] are contained in Vstar.")
else
    disp("The columns of [BP2*G; BM] are NOT contained in Vstar")
end



