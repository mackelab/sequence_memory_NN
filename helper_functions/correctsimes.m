function cp=correctsimes(pvals,varargin)
% function cp=correctsimons(p,N)
% ncomp = size(p,1), 

[ncomp, ncases]=size(pvals);
cp=zeros(ncomp, ncases);

% simons correction factors
cr=repmat((ncomp./[1:ncomp])',1,ncases);

% sort p values
[sp,ix]=sort(pvals);

% multiply with correction factor
temp_p=sp.*cr;

% bring back in correct order per comparison
for k=1:size(sp, 2)
    cp(ix(:,k),k)=temp_p(:,k);
end