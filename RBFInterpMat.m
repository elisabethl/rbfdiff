% -------------------------------------------------------------------------
% RBFInterpMat
% Purpose: Given an RBF, a polynomial degree, a shape parameter and node
%          points, set up an RBF interpolation matrix to be used later by the
%          companion subroutine RBFDiffMat to compute evaluation and
%          differentiation matrices for the RBF approximation. RBFInterpMat
%          provides a unified interface for RBF-Direct and RBF-QR.
% Psi = RBFInterpMat(phi,pdeg,ep,xc,xe)
% Psi = RBFInterpMat(phi,pdeg,ep,xc,C,R)
% Input:  phi    string with RBF name {'gs'|'mq'|'Bmq'|'iq'|'phs'|'rbfqr'} 
%         pdeg	 integer scalar, the polynomial degree to add to the basis
%                use -1 for no polynomial terms.
%         ep     double scalar, the shape parameter of the RBF if applicable
%         xc     double(N,dim) the node points (centers) for the
%                interpolation/approximation
%
% Output: Psi	 structure that contains the parameters of the RBF approximation
%                the node points, and the factorized interpolation matrix.
%                The content of the structure differs between RBF-Direct
%                approximations and RBF-QR approximations.
%
% Example:       Set up the interpolation matrix of RBF-Direct with 21 uniform
%                nodes for a unit patch in 2D using multiquadric RBFs with
%                shape parameter 1.
%                H = 1; R = 1; dim = 2; nodeGen = 'uni'; nLoc = 21;
%                phi = 'mq'; pdeg = -1; ep = 1; 
%                xc = getPtchPts(H,R,dim,nLoc,nodeGen);
%                Psi = RBFInterpMat(phi,pdeg,ep,xc);
%
% Copyright (c) 2024 Elisabeth Larsson <elisabeth.larsson@it.uu.se>
%		       Andreas Michael <andreas.michael@it.uu.se >
% -------------------------------------------------------------------------
function Psi = RBFInterpMat(phi,pdeg,ep,xc,varargin)
narg = length(varargin);
dim = size(xc,2);
ndiff = 0; % We are only going to compute Psi here

if strcmp(phi,'rbfqr')
    callstr1 = '[out1,out2] = RBF_QR_diffmat_';
    callstr2 = 'D(ndiff,xc,ep,';
    if (narg == 2)
        C = varargin{1};
        R = varargin{2};
        callstr3 = 'C,R);';
    elseif (narg==1)
        xe = varargin{1};
        callstr3 = 'xe);';
    end
    callstr = strcat(callstr1,num2str(dim),callstr2,callstr3);
    eval(callstr);

    if narg == 2
        Psi = out1;
    else
        Psi = out2;
    end
    Psi.phi = phi;
else
    %
    % Scale the points and shape parameter based on the given inputs
    %
    if (narg >= 2)
        cc = varargin{1};
        rr = varargin{2};
    elseif (narg==1)
        xe = varargin{1};
        %cc = mean(xe);
        cc = mean(xc); % The centers define the domain. Even for a stencil.
        xe = xe - cc;
        re = sqrt(sum(xe.^2,2));
        rr = max(re);
    end
    xloc = xc - cc; % xloc are the local shifted and scaled centers
    rc = sqrt(sum(xloc.^2,2));
    rr = max(max(rc),rr);
    Psi.rr = rr;
    Psi.cc = cc;

    xloc = xloc./rr;
    ep = ep*rr;
    %
    % RBF-Direct. Compute distance matrix
    %
    rc = xcdist(xloc);
    A = RBFmat(phi,ep,rc,'0');
    if (pdeg>=0)
        %
        % Add polynomial terms
        %
        P = polyMat(xloc,pdeg);
        np = size(P,2);
        A = [A P;P' zeros(np)];
    end
    Psi.A0 = A;
    Psi.ep = ep;
    Psi.phi = phi;
    Psi.pdeg = pdeg;
    Psi.xloc = xloc;
    %
    % Factorize A. A(piv,:) = L*U, piv is a permutation vector. 
    % Cholesky A = L'*L; can be used for pos. def. matrices.
    % But RBFs are not always pos def, e.g., MQ is not
    %
    [Psi.L,Psi.U,Psi.piv] = lu(A,'vector');
end
Psi.xc = xc;

