function [unwrap] = cunwrap4D(wrap)

% Source: http://cognitive-eurhythmics.com/code.html
%
% Purpose: Costantini's 2D unwrapping based on Minimum Network Flow
% References:
%   - http://earth.esa.int/workshops/ers97/papers/costantini/
%   - Costantini, M. (1998) A novel phase unwrapping method based on network
%   programming. IEEE Tran. on Geoscience and Remote Sensing, 36, 813-821.
%
% Author: Eric Barnhill
% 
% Dependencies:
% Matlab Optimization Toolbox

if ndims(wrap)~=4
    error('CUNWRAP4D: input must be 4D array.');
end

start = cputime;

[ny, nx, nz, nt] = size(wrap);

% pad

wrap2 = zeros(ny+2, nx+2, nz+2, nt+2);
wrap2(2:end-1, 2:end-1, 2:end-1, 2:end-1) = wrap;
wrap = wrap2;

ny = ny+2;
nx = nx+2;
nz = nz+2;
nt = nt+2;



if nx<2 || ny<2 || nz<2 || nt<2
    error('CUNWRAP4D: size of input must be larger than 2x2x2');
end
                      
% get gradients , use a priori knowledge that most gradients are in [-pi,pi)

gradY = wrap(2:ny, 1:nx, 1:nz, 1:nt) - wrap(1:ny-1, 1:nx, 1:nz, 1:nt);
gradY = mod(gradY+pi, 2*pi)-pi;

gradX = wrap(1:ny, 2:nx, 1:nz, 1:nt) - wrap(1:ny, 1:nx-1, 1:nz, 1:nt);
gradX = mod(gradX+pi,2*pi)-pi;

gradZ = wrap(1:ny, 1:nx, 2:nz, 1:nt) - wrap(1:ny, 1:nx, 1:nz-1, 1:nt);
gradZ = mod(gradZ+pi, 2*pi)-pi;

gradT = wrap(1:ny, 1:nx, 1:nz, 2:nt) - wrap(1:ny, 1:nx, 1:nz, 1:nt-1);
gradT = mod(gradT+pi, 2*pi)-pi;




%% SET UP beq

gradXY = gradY(1:ny-1, 2:nx, 1:nz, 1:nt) - gradY(1:ny-1, 1:nx-1, 1:nz, 1:nt);
gradYX = gradX(2:ny, 1:nx-1, 1:nz, 1:nt) - gradX(1:ny-1, 1:nx-1, 1:nz, 1:nt);

gradZX = gradX(1:ny, 1:nx-1, 2:nz, 1:nt) - gradX(1:ny, 1:nx-1, 1:nz-1, 1:nt);
gradXZ = gradZ(1:ny, 2:nx, 1:nz-1, 1:nt) - gradZ(1:ny, 1:nx-1, 1:nz-1, 1:nt);

gradZY = gradY(1:ny-1, 1:nx, 2:nz, 1:nt) - gradY(1:ny-1, 1:nx, 1:nz-1, 1:nt);
gradYZ = gradZ(2:ny, 1:nx, 1:nz-1, 1:nt) - gradZ(1:ny-1, 1:nx, 1:nz-1, 1:nt);

gradXT = gradT(1:ny, 2:nx, 1:nz, 1:nt-1) - gradT(1:ny, 1:nx-1, 1:nz, 1:nt-1);
gradTX = gradX(1:ny, 1:nx-1, 1:nz, 2:nt) - gradX(1:ny, 1:nx-1, 1:nz, 1:nt-1);

gradYT = gradT(2:ny, 1:nx, 1:nz, 1:nt-1) - gradT(1:ny-1, 1:nx, 1:nz, 1:nt-1);
gradTY = gradY(1:ny-1, 1:nx, 1:nz, 2:nt) - gradY(1:ny-1, 1:nx, 1:nz, 1:nt-1);

gradZT = gradT(1:ny, 1:nx, 2:nz, 1:nt-1) - gradT(1:ny, 1:nx, 1:nz-1, 1:nt-1);
gradTZ = gradZ(1:ny, 1:nx, 1:nz-1, 2:nt) - gradZ(1:ny, 1:nx, 1:nz-1, 1:nt-1);

beq1 = -1/(2*pi)*(gradZY - gradYZ); 
beq2 = -1/(2*pi)*(gradXZ - gradZX);
beq3 = -1/(2*pi)*(gradYX - gradXY);
beq4 = -1/(2*pi)*(gradXT - gradTX);
beq5 = -1/(2*pi)*(gradYT - gradTY);
beq6 = -1/(2*pi)*(gradZT - gradTZ);
beq = [beq1(:); beq2(:); beq3(:); beq4(:); beq5(:); beq6(:)];

masterWeights1 = ones(ny, nx, nz, nt);
masterWeights1(1:ny-1, 1:nx-1, 1:nz-1, 1:nt-1) = masterWeights1(1:ny-1, 1:nx-1, 1:nz-1, 1:nt-1) + abs(beq1(1:ny-1, 1:nx-1, 1:nz-1, 1:nt-1)) + abs(beq2(1:ny-1, 1:nx-1, 1:nz-1, 1:nt-1)) + abs(beq3(1:ny-1, 1:nx-1, 1:nz-1, 1:nt-1));
masterWeights2 = 1./(masterWeights1.^2);
lowpass = ones(3,3,3,3)/81;
masterWeights3 = convn(masterWeights2, lowpass, 'same');



%% SET UP Aeq
% 3D curl of U = (Uzy - Uyz) + (Uxz - Uzx) - (Uyx - Uxy)
% columns - pixels in current gradient matrix // rows - pixels in curl matrix


% ZY
rows = (ny-1)*nx*(nz-1);
columns = (ny-1)*nx*nz;
diagonals = [-ones(rows,1) ones(rows,1)];
ZYpos = spdiags(diagonals, [0 (ny-1)*nx], rows, columns);
ZYpos = kron(eye(nt), ZYpos);

ZYneg = -ZYpos;

% YZ
rows = ny-1;
columns = ny;
diagonals = [-ones(ny-1,1) ones(ny-1,1)];
block = spdiags(diagonals, [0 1], rows, columns);
YZblock1 = -kron(eye(nx), block);
YZpos = kron(eye(nz-1), YZblock1);
YZpos = kron(eye(nt), YZpos);
YZneg = -YZpos;

% XZ
rows = ny*(nx-1);
columns = ny*nx;
diagonals = [-ones(rows,1) ones(rows,1)];
XZblock = spdiags(diagonals, [0 nx-1], rows, columns);
XZpos = kron(eye(nz-1), XZblock);
XZpos = kron(eye(nt), XZpos);
XZneg = -XZpos;

% ZX
rows = ny*(nx-1)*(nz-1);
columns = ny*(nx-1)*nz;
diagonals = [-ones(rows,1) ones(rows,1)];
ZXpos = spdiags(diagonals, [0 ny*(nx-1)], rows, columns);
ZXpos = kron(eye(nt), ZXpos);
ZXneg = -ZXpos;

% YX
rows = ny-1;
columns = ny;
diagonals = [-ones(ny-1,1) ones(ny-1,1)];
block = spdiags(diagonals, [0 1], rows, columns);
YXblock1 = -kron(eye(nx-1), block);
YXpos = kron(eye(nz), YXblock1);
YXpos = kron(eye(nt), YXpos);
YXneg = -YXpos;

% XY
rows = (ny-1)*(nx-1);
columns = (ny-1)*nx;
diagonals = [-ones(rows,1) ones(rows,1)];
XYblock = spdiags(diagonals, [0 nx-1], rows, columns);
XYpos = kron(eye(nz), XYblock);
XYpos = kron(eye(nt), XYpos);
XYneg = -XYpos; 

% XT
rows = ny*(nx-1);
columns = ny*nx;
diagonals = [-ones(rows,1) ones(rows,1)];
XTblock = spdiags(diagonals, [0 nx-1], rows, columns);
XTpos = kron(eye(nz), XTblock);
XTpos = kron(eye(nt-1), XTpos);
XTneg = -XTpos;

% TX
rows = ny*(nx-1)*nz*(nt-1);
columns = ny*(nx-1)*nz*nt;
diagonals = [-ones(rows,1) ones(rows,1)];
TXpos = spdiags(diagonals, [0 ny*(nx-1)*nz*nt], rows, columns);
TXneg = -TXpos;

% YT
rows = ny-1;
columns = ny;
diagonals = [-ones(ny-1,1) ones(ny-1,1)];
block = spdiags(diagonals, [0 1], rows, columns);
YTblock1 = -kron(eye(nx), block);
YTpos = kron(eye(nz), YTblock1);
YTpos = kron(eye(nt-1), YTpos);
YTneg = -YTpos;

% TY
rows = (ny-1)*nx*nz*(nt-1);
columns = (ny-1)*nx*nz*nt;
diagonals = [-ones(rows,1) ones(rows,1)];
TYpos = spdiags(diagonals, [0 (ny-1)*nx*nz*nt], rows, columns);
TYneg = -TYpos;

% ZT
rows = ny*nx*(nz-1);
columns = ny*nx*nz;
diagonals = [-ones(rows,1) ones(rows,1)];
ZTpos = spdiags(diagonals, [0 ny*nx], rows, columns);
ZTpos = kron(eye(nt-1), ZTpos);
ZTneg = -ZTpos;

% TZ
rows = ny*nx*(nz-1)*(nt-1);
columns = ny*nx*(nz-1)*nt;
diagonals = [-ones(rows,1) ones(rows,1)];
TZpos = spdiags(diagonals, [0 ny*nx*(nz-1)*nt], rows, columns);
TZneg = -TZpos;


% Matrix of the LHS of eqt (17)
%Aeq = [ ZYpos ZYneg YZpos YZneg XZpos XZneg ZXpos ZXneg YXpos YXneg XYpos XYneg ];
Aeq = blkdiag([ZYpos ZYneg YZpos YZneg], [XZneg XZpos ZXpos ZXneg], [YXneg YXpos XYneg XYpos], ...
    [XTpos XTneg TXpos TXneg], [YTpos YTneg TYpos TYneg], [ZTneg ZTpos TZpos TZneg] );
nvars = size(Aeq,2);

%% SET UP LP

zySize = (ny-1)*nx*nz*nt;
zyShape = [ny-1 nx nz nt];

yzSize = ny*(nx)*(nz-1)*nt;
yzShape = [ny nx nz-1 nt];

xzSize = (ny)*(nx)*(nz-1)*nt;
xzShape = [ny nx nz-1 nt];

zxSize = ny*(nx-1)*nz*nt;
zxShape = [ny nx-1 nz nt];

yxSize = ny*(nx-1)*(nz)*nt;
yxShape = [ny nx-1 nz nt];

xySize = (ny-1)*nx*(nz)*nt; 
xyShape = [ny-1 nx nz nt];

xtsize = size(XTpos)
xtSize = ny*nx*(nz)*(nt-1)
xtShape = [ny nx nz nt-1];

txSize = ny*(nx-1)*(nz)*nt;
txShape = [ny nx-1 nz nt];

ytSize = ny*nx*(nz)*(nt-1); 
ytShape = [ny nx nz nt-1];

tySize = (ny-1)*nx*(nz)*nt; 
tyShape = [ny-1 nx nz nt];

ztSize = ny*nx*(nz)*(nt-1); 
ztShape = [ny nx nz nt-1];

tzSize = (ny)*(nx)*(nz-1)*nt;
tzShape = [ny nx nz-1 nt];



% Cost function, eqt (16) - Y Z Z X X Y T X T Y T Z

cost1 = masterWeights3(1:ny-1, 1:nx, 1:nz, 1:nt);
cost2 = masterWeights3(1:ny, 1:nx, 1:nz-1, 1:nt);
cost3 = cost2;
cost4 = masterWeights3(1:ny, 1:nx-1, 1:nz, 1:nt);
cost5 = cost4;
cost6 = cost1;
cost7 = masterWeights3(1:ny, 1:nx, 1:nz, 1:nt-1);
cost8 = cost4; 
cost9 = cost7;
cost10 = cost1;
cost11 = cost7;
cost12 = cost2;

cost = [cost1(:); cost1(:); cost2(:); cost2(:); cost3(:); cost3(:); cost4(:); cost4(:); cost5(:); cost5(:); cost6(:); cost6(:); ...
            cost7(:); cost7(:); cost8(:); cost8(:); cost9(:); cost9(:); cost10(:); cost10(:); cost11(:); cost11(:); cost12(:); cost12(:)];

% Lower and upper bound, eqt (18,19)
L = zeros(nvars,1);
U = Inf(size(L)); % No upper bound, U=[];

costsize = size(cost)
aeqsize = size(Aeq)


% Call LP solver

    x = linprog(cost,[],[],Aeq,beq,L,U);


offset = 1;
zyp = reshape(x(offset:zySize), zyShape);
offset = offset + zySize;
zyn = reshape(x(offset:offset+zySize-1), zyShape);
offset = offset + zySize; 
yzp = reshape(x(offset:offset+yzSize-1), yzShape);
offset = offset + yzSize;
yzn = reshape(x(offset:offset+yzSize-1), yzShape);
offset = offset + yzSize;
xzp = reshape(x(offset:offset+xzSize-1), xzShape);
offset = offset + xzSize;
xzn = reshape(x(offset:offset+xzSize-1), xzShape);
offset = offset + xzSize;
zxp = reshape(x(offset:offset+zxSize-1), zxShape);
offset = offset + zxSize;
zxn = reshape(x(offset:offset+zxSize-1), zxShape);
offset = offset + zxSize;
yxp = reshape(x(offset:offset+yxSize-1), yxShape);
offset = offset + yxSize;
yxn = reshape(x(offset:offset+yxSize-1), yxShape);
offset = offset + yxSize;
xyp = reshape(x(offset:offset+xySize-1), xyShape);
offset = offset + xySize;
xyn = reshape(x(offset:offset+xySize-1), xyShape);
offset = offset + xySize;
xtp = reshape(x(offset:offset+xtSize-1), xtShape);
offset = offset + xtSize;
xtn = reshape(x(offset:offset+xtSize-1), xtShape);
offset = offset + xtSize;
txp = reshape(x(offset:offset+txSize-1), txShape);
offset = offset + txSize;
txn = reshape(x(offset:offset+txSize-1), txShape);
offset = offset + txSize;
ytp = reshape(x(offset:offset+ytSize-1), ytShape);
offset = offset + ytSize;
ytn = reshape(x(offset:offset+ytSize-1), ytShape);
offset = offset + ytSize;
typ = reshape(x(offset:offset+tySize-1), tyShape);
offset = offset + tySize;
tyn = reshape(x(offset:offset+tySize-1), tyShape);
offset = offset + tySize;
ztp = reshape(x(offset:offset+ztSize-1), ztShape);
offset = offset + ztSize;
ztn = reshape(x(offset:offset+ztSize-1), ztShape);
offset = offset + ztSize;
tzp = reshape(x(offset:offset+tzSize-1), tzShape);
offset = offset + tzSize;
tzn = reshape(x(offset:offset+tzSize-1), tzShape);

% Compute the derivative jumps, eqt (20,21)
zy = zyp - zyn;
yz = yzp - yzn;
xz = xzp - xzn;
zx = zxp - zxn;
yx = yxp - yxn;
xy = xyp - xyn;
xt = xtp - xtn;
tx = txp - txn;
yt = ytp - ytn;
ty = typ - tyn;
zt = ztp - ztp;
tz = tzp - tzn;

% Sum the jumps with the wrapped partial derivatives, eqt (10,11)

zy = zy(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradY(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
yz = yz(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradZ(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
xz = xz(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradZ(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
zx = zx(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradX(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
yx = yx(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradX(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
xy = xy(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradY(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
xt = xt(1:ny-1,1:nx-1,1:nz-1, 1:nt-1) + gradT(1:ny-1,1:nx-1,1:nz-1, 1:nt-1)/(2*pi);
% add more if necessary

%%
% Integrate the partial derivatives, eqt (6)
% 

innerxy = xy(2:end-1, 2:end, 2:end, 2);
inneryx = yx(2:end-1, 2:end-1, 2:end-1, 2);
inneryz = yz(2:end, 2:end, 2:end-1, 2);
innerxt = xt(2:end, 2:end, 2:end, 2:end-1);
unwrap1 = zeros(ny-2, nx-2, nz-2, 1);
unwrap1(1:end-1,2:end,1:end-1,1) = cumsum(inneryx,2);
unwrap1(2:end, 1:end, 1:end) = innerxy;
unwrap1 = cumsum(unwrap1,1);
unwrap2 = unwrap1;
unwrap2(:,:,2:end, 1) = inneryz;
unwrap3 = cumsum(unwrap2, 3);
unwrap4 = zeros(ny-2, nx-2, nz-2, nt-2);
unwrap4(:,:,:,1) = unwrap3;
unwrap5 = unwrap4;
unwrap5size = size(unwrap5)
xtsize = size(innerxt)

unwrap5(:,:,:,2:end) = innerxt;
unwrap6 = cumsum(unwrap5, 4);

unwrap = 2*pi*unwrap6;

endtime = cputime - start
end % cunwrap



