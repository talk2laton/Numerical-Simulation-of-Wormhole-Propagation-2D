function recharacterized(q)
%% Work Title   : Thesis.
% Submited By   : Kareem Lateef Adewale
%..........................................................................
%%
data     = load('Thesiscore_diffussion_2D.txt');
L        = data(1);    % core length (inches)
w        = data(2);    % core width (inches)
Nx       = data(3);    % Ny grid number
Ny       = data(4);    % Nz grid number
Pout     = data(5);    % Out let pressure (atm)
miu_acid = data(15);   % viscosity of Acid (centipoise)
Ks       = data(22);   % Initial Permeability
phio     = data(24);   % Iniital mean porosity
delphi   = data(25);   % heterogeity factor

%%
Lp      = 2.54*L;           % core length (centimeters)
wp      = 2.54*w;           % core width(centimeters)
qp      = 1.66666667e-2*q;  % flow rate (cc/s) converted from (cc/min)

%% Grid size
M       = Nx*Ny;
Mp1     = M+1;

%% Geometry of Core
gridordinate = zeros(M,2);
gridX = zeros(Nx, Ny);
gridY = zeros(Nx, Ny);
gridI = zeros(Nx, Ny);
for m = 1:M
    nx = floor((m - 1)/Ny) + 1;
    ny = m - Ny*(nx - 1); 
    gridordinate(m,:) = [nx, ny];
    gridX(nx, ny) = nx;
    gridY(nx, ny) = ny;
    gridI(nx, ny) = m;
end
Nnx  = gridordinate(:,1);
Nny  = gridordinate(:,2);

ymid       = (Ny + 1)/2;
radialedge = Ny - ymid;

radialdist = abs(Nny - ymid);
indx       = (radialdist <= radialedge)&(Nnx < Nx)&(Nnx > 1); 
full_indx  = (1:M)';

IDm    = full_indx(indx); 
NIDm   = full_indx(~indx); 

Injloc  = Nnx(IDm) == 2;       
Exitloc = Nnx(IDm) == Nx - 1;

% Index of injection and exit face in the whole grid
INJloc = (radialdist <= radialedge)&(Nnx == 2);     
EXTloc = (radialdist <= radialedge)&(Nnx == Nx - 1);

INJm  = full_indx(INJloc);
EXTm  = full_indx(EXTloc);

INJFace = (radialdist <= radialedge)&(Nnx == 1);     
EXTFace  = (radialdist <= radialedge)&(Nnx == Nx);

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
dzp      = dyp;
Axp      = dzp*dyp;
Ayp      = dzp*dxp;

mx       = rand(M,1);
my       = rand(M,1);
Kx       = K.*mx;
Ky       = K.*my;

%% Recharacterization

IDmm1    = IDm - 1;      IDmp1    = IDm + 1;
IDmmNy   = IDm - Ny;     IDmpNy   = IDm + Ny;
MP1      = repmat(Mp1, sum(Injloc), 1);
MeanK    = @(K, m1, m2) harmmean([K(m1),K(m2)],2);
                           
TYmmm1    = MeanK(Ky, IDm, IDmm1)*Ayp/(miu_acid*dyp);    
TYmmp1    = MeanK(Ky, IDm, IDmp1)*Ayp/(miu_acid*dyp);    
TXmmmNy   = MeanK(Kx, IDm, IDmmNy)*Axp/(miu_acid*dxp); 
TCOMPXinj = TXmmmNy(Injloc);
TXmmpNy   = MeanK(Kx, IDm, IDmpNy)*Axp/(miu_acid*dxp); 
BXmmpNy   = -TXmmpNy(Exitloc)*Pout; 
STmps     = TXmmmNy + TYmmm1 + TYmmp1 + TXmmpNy;

AIDnV     = [IDm,          IDm,                 -STmps
             IDm,          IDmm1,                TYmmm1
             IDm,          IDmp1,                TYmmp1
             IDm,          IDmmNy,               TXmmmNy
             IDm,          IDmpNy,               TXmmpNy
             IDm(Injloc),  MP1,                  TCOMPXinj
             MP1,          IDm(Injloc),         -TCOMPXinj
             Mp1,          Mp1,                  sum(TCOMPXinj)];
            
BIDnV     = [IDm(Exitloc), ones(sum(Exitloc),1), BXmmpNy
             Mp1,          1,                    qp];

A = sparse(AIDnV(:,1), AIDnV(:,2), AIDnV(:,3), Mp1, Mp1);
b = sparse(BIDnV(:,1), BIDnV(:,2), BIDnV(:,3), Mp1, 1);

A(NIDm,:) = []; A(:,NIDm) = [];
b(NIDm)   = []; b  = full(b);

P = A\b;
Pinj1 = P(end);

perm     = qp*miu_acid*Lp/(wp*dzp*(Pinj1-Pout));
Correcta = Ks*1e-3/perm;
K        = Correcta*K;

minP = Pout;
maxP = Pinj1;
mnP = minP;
mxP = minP + (maxP-minP)/Correcta;

Pguess = interp1([minP, maxP], [mnP, mxP], P);
qpref  = qp;
    
save Datafile Phidist pore_WH K Correcta aporos afrac mx my IDm Ny IDmm1...
     IDmp1 IDmmNy IDmpNy MP1 Injloc Exitloc L Nx NIDm gridordinate gridX...
     gridY Pguess gridI qpref Pinj1 Pout
