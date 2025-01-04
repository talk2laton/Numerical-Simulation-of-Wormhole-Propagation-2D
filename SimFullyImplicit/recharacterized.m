function recharacterized(q)
%% Work Title   : Thesis.
% Submited By   : Kareem Lateef Adewale
%..........................................................................
%%
data     = load('Thesiscore_diffussion3D.txt');
L        = data(1);    % core length (inches)
w        = data(2);    % core width(inches)
Ny       = data(3);    % Ny grid number
Nz       = data(4);    % Nz grid number
Pout     = data(5);    % Out let pressure (atm)
miu_acid = data(15);   % viscosity of Acid   (centipoise)
Ks       = data(23);   % Initial Permeability
phio     = data(25);   % Iniital mean porosity
delphi   = data(26);   % heterogeity factor

%%
Lp      = 2.54*L;           % core length (centimeters)
wp      = 2.54*w;           % core width(centimeters)
qp      = 1.66666667e-2*q;  % flow rate (cc/s) converted from (cc/min)

%% Grid size
Nx      = floor(Ny*L/w);
NyNz    = Ny*Nz;
M       = Nx*Ny*Nz;
Mp1     = M+1;

%% Geometry of Core
gridordinate = zeros(M,3);
gridX = zeros(Nx, Ny, Nz);
gridY = zeros(Nx, Ny, Nz);
gridZ = zeros(Nx, Ny, Nz);
gridI = zeros(Nx, Ny, Nz);
for m = 1:M
    nx    = floor((m - 1)/NyNz) + 1;
    Rem   = m - NyNz*(nx - 1);
    ny    = floor((Rem - 1)/Nz) + 1;
    nz    = Rem - Nz*(ny - 1); 
    gridordinate(m,:) = [nx, ny, nz];
    gridX(nx, ny, nz) = nx;
    gridY(nx, ny, nz) = ny;
    gridZ(nx, ny, nz) = nz;
    gridI(nx, ny, nz) = m;
end
Nnx      = gridordinate(:,1);
Nny      = gridordinate(:,2);
Nnz      = gridordinate(:,3);

zmid       = (Nz + 1)/2;
ymid       = (Ny + 1)/2;
radialedge = (Nz - 3)/2;

radialdist  = sqrt((Nnz - zmid).^2 + (Nny - ymid).^2);
indx = (radialdist <= radialedge)&(Nnx < Nx)&(Nnx > 1); 
full_indx = (1:M)';

IDm  = full_indx(indx);   NIDm  = full_indx(~indx);   

Injloc  = Nnx(IDm) == 2;       
Exitloc = Nnx(IDm) == Nx - 1;

%% Core
pore_WH  = 0.4; %porosity of a grid containing a wormhole
Phidist  = phio + delphi*rand(M,1);
aporos   = log(1e4)/(pore_WH - phio);
afrac    = log(1e9)/(1.0000 - pore_WH);
permdist = 1e-6*(exp(aporos*(Phidist - phio)).*(Phidist < pore_WH) + ...
                 exp(afrac*(Phidist - phio)).*(Phidist >= pore_WH));
permdist = permdist.*(radialdist <= radialedge); 
           
permdist(Nnx == Nx & radialdist <= radialedge) = Inf; 
permdist(Nnx == 1 & radialdist <= radialedge) = Inf;
K       = permdist;


%%
dxp      = Lp/Nx;
dyp      = wp/Ny;
dzp      = Hp/Nz;
Axp      = dzp*dyp;
Ayp      = dzp*dxp;
Azp      = dxp*dyp;

mx       = rand(M,1);
my       = rand(M,1);
mz       = rand(M,1);
Kx       = K.*mx;
Ky       = K.*my;
Kz       = K.*mz;

%% Recharacterization

IDmm1    = IDm - 1;      IDmp1    = IDm + 1;
IDmmNz   = IDm - Nz;     IDmpNz   = IDm + Nz;
IDmmNyNz = IDm - NyNz;   IDmpNyNz = IDm + NyNz;
MP1      = repmat(Mp1, sum(Injloc), 1);
MeanK    = @(K, m1, m2) harmmean([K(m1),K(m2)],2);

TZmmm1    = MeanK(Kz, IDm, IDmm1)*Azp/(miu_acid*dzp);    
TZmmp1    = MeanK(Kz, IDm, IDmp1)*Azp/(miu_acid*dzp);    
TYmmmNz   = MeanK(Ky, IDm, IDmmNz)*Ayp/(miu_acid*dyp);    
TYmmpNz   = MeanK(Ky, IDm, IDmpNz)*Ayp/(miu_acid*dyp);    
TXmmmNyNz = MeanK(Kx, IDm, IDmmNyNz)*Axp/(miu_acid*dxp); 
TCOMPXinj = TXmmmNyNz(Injloc);
TXmmpNyNz = MeanK(Kx, IDm, IDmpNyNz)*Axp/(miu_acid*dxp); 
BXmmpNyNz = - TXmmpNyNz(Exitloc)*Pout; 
STmps     = TXmmmNyNz + TYmmmNz + TZmmm1 + TZmmp1 + TYmmpNz + TXmmpNyNz;

AIDnV     = [IDm,          IDm,                    -STmps
             IDm,          IDmm1,                  TZmmm1
             IDm,          IDmp1,                  TZmmp1
             IDm,          IDmmNz,                 TYmmmNz
             IDm,          IDmpNz,                 TYmmpNz
             IDm,          IDmmNyNz,               TXmmmNyNz
             IDm,          IDmpNyNz,               TXmmpNyNz
             IDm(Injloc),  MP1,                    TCOMPXinj
             MP1,          IDm(Injloc),           -TCOMPXinj
             Mp1,          Mp1,                    sum(TCOMPXinj)];
            
BIDnV     = [IDm(Exitloc), ones(sum(Exitloc),1),   BXmmpNyNz
             Mp1,          1,                      qp];

A = sparse(AIDnV(:,1), AIDnV(:,2), AIDnV(:,3), Mp1, Mp1);
b = sparse(BIDnV(:,1), BIDnV(:,2), BIDnV(:,3), Mp1, 1);

A(NIDm,:) = []; A(:,NIDm) = [];
b(NIDm)   = []; b  = full(b);

P = A\b;
Pinj1 = P(end);

perm     = 4*qp*miu_acid*Lp/(pi*wp^2*(Pinj1-Pout));
Correcta = Ks*1e-3/perm;
K        = Correcta*K;

minP = Pout;
maxP = Pinj1;
mnP = minP;
mxP = minP + (maxP-minP)/Correcta;

Pguess = interp1([minP, maxP], [mnP, mxP], P);
qpref  = qp;
    
save Datafile Phidist pore_WH K Correcta aporos afrac mx my mz IDm H Ny ...
     IDmm1 IDmp1 IDmmNz IDmpNz IDmmNyNz IDmpNyNz MP1 Injloc Exitloc L Nx...
     NIDm gridordinate gridX gridY gridZ Pguess gridI qpref Pinj1 Pout
