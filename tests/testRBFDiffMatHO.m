close all
clear all
%
% The purpose here is to test the RBFDiffMat works for the derivatives of order
% 3 and higher 
%
phi = 'mq';
ep = 0.8;
d = 2;
%
% We use the already implemented test of the derivatives by direct call for
% comparison. This only evaluates at one point, but the matrix functionality
% has already been tested in other ways.
%
[phival,xe,xc]=testRBFder(phi,ep,d);
%
% Now, we compute the third order and fourth order derivatives using
% RBFDiffMat
%
pdeg = -1;
Psi = RBFInterpMat(phi,pdeg,ep,xc,xe);
D3 = RBFDiffMat(3,Psi,xe);
D4 = RBFDiffMat(4,Psi,xe);

i3 = dim(2,d)+1:dim(3,d);
i4 = dim(3,d)+1:dim(4,d);

pos=0;
for i=1:d
    for j=i:d
        for k=j:d
            pos = pos+1;
            y1 = phival(:,i3(pos));
            y2 = D3{i,j,k};
            figure
            plot(y1,'-'),hold on
            plot(y2,'--')
            relerr3(i,j,k)=max(abs(y1-y2))/max(abs(y1));
        end
    end
end

pos=0;
for i=1:d
    for j=i:d
        for k=j:d
            for ell=k:d
                pos = pos+1;
                y1 = phival(:,i4(pos));
                y2 = D4{i,j,k,ell};
                relerr4(i,j,k,ell)=max(abs(y1-y2))/max(abs(y1));
            end
        end
    end
end

relerr3
relerr4
%
% Also test the 4th order Laplacian as an example of higher Laplacians
%
r = xcdist(xe,xc);
L2phi = feval(phi,ep,r,'L2',d);

DL2 = RBFDiffMat(3.5,Psi,xe);

relerrL2 = max(abs(L2phi-DL2))/max(abs(L2phi))
