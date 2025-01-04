function PVfinal = NumCoreSimFullyImplicit_2D(q, PVlimit)
    %% Work Title   : Thesis.
    % Submited By   : Kareem Lateef Adewale
    close all;      % closes all the open figure windows
    clc;            % clear command window
    
    %% Input Data
    data      = load('Thesiscore_diffussion_2D.txt');
    L         = data(1);       % core length (inches)
    w         = data(2);       % core width(inches)
    Nx        = data(3);       % Ny grid number
    Ny        = data(4);       % Nz grid number
    Pout      = data(5);       % Out let pressure (atm)
    Temp      = data(7);       % temperature (K)      
    Sm        = data(10);      % Sepecific area for fast reacting mineral (Calcite-m2/m3)
    DeltaE_R  = data(11);      % constant of arrehnus rate law for (fast-HCl)
    Ef0       = data(12);      % Reaction rate constant Ef0 (fast-HCl) kg-mole mineral/(m2-sec-(kg-mole HCl/m3) acid)
    rhoM      = data(13);      % density of fast reacting mineral (-kg/m3 (density of Calcite))
    Conc_HCl  = data(14);      % concentration of HF (wt%)
    miu_acid  = data(15);      % viscosity of Acid (centipoise)
    rho_acid  = data(16);      % density of acid (kg/m3)
    v1        = data(17);      % volume fraction of fast mineral(Assuming are Calcite)
    MW_HCl    = 1e-3*data(18); % Molecular weight of Hydrochloric acid(Kilograms/mol)
    MW_CaCO3  = 1e-3*data(19); % Molecular weight of Calcite(Kilograms/mol)
    nA        = data(20);      % Stoichiometric Coefficient of HCl in Calcite-HCl reaction
    nM        = data(21);      % Stoichiometric Coefficient of Calcite in Calcite-HCl reaction
    Ks        = data(22);      % Initial Permeability
    zymultply = data(23);      % horizontal permeability factor
    phio      = data(24);      % Iniital mean porosity
    delPhi    = data(25);      % Heterogeity factor
    data2     = load('ReactiveFlow.txt');
    alphax    = data2(1);      % coefficient of diffusion
    Dm        = data2(2);      % Molecular diffusivity (m^2/sec)
    kc        = data2(3);      % Local diffusion constant (sec^-1)
    
    %% Result file names
    if q<1
        ID = ['0_',num2str(10*q),'_'];
    else
        ID = [num2str(q),'_'];
    end
    FileName1 = ['Concentration_',ID];
    FileName2 = ['Porosity_',ID];
    FileName3 = ['VolumeFraction_',ID];
    FileName4 = ['Vectores_',ID];
    
    %% Units
    % Darcy Unit for Pressure Calculations
    Lp   = 2.54*L;          % core length (centimeters)
    wp   = 2.54*w;          % core width(centimeters)
    qp   = 1.66666667e-2*q; % flow rate (cc/s)
    
    % SI Unit for reaction
    Lc   = 0.0254*L;        % core length (meters)
    wc   = 0.0254*w;        % core width(meters)
    qc   = 1.66666667e-8*q; % flow rate (cubic meters/s)
    
    %% Grid size
    M          = Nx*Ny;
    Mp1        = M+1;
    Storagecap = floor(1e8/M);
    
    disp([num2str(M),' grids full geometry'])
    disp(' ');
    %% Geometry of Core
    gridordinate = zeros(M,2);
    for m = 1:M
        nx = floor((m - 1)/Ny) + 1;
        ny = m - Ny*(nx - 1); 
        gridordinate(m,:) = [nx, ny];
    end
    Nnx  = gridordinate(:,1);
    Nny  = gridordinate(:,2);
    
    ymid  = (Ny + 1)/2;
    
    radialedge = Ny - ymid;
    radialdist = abs(Nny - ymid);
    indx       = (radialdist <= radialedge)&(Nnx < Nx)&(Nnx > 1); 
    full_indx  = (1:M)';
    
    IDm     = full_indx(indx); 
    NIDm    = full_indx(~indx); 
    Mnempty = numel(IDm);
    

    % index of injection and exit face in interest region
    Injloc  = Nnx(IDm) == 2;       
    Exitloc = Nnx(IDm) == Nx - 1;

    % Index of injection and exit face in the whole grid
    INJloc = (radialdist <= radialedge)&(Nnx == 2);     
    EXTloc = (radialdist <= radialedge)&(Nnx == Nx - 1);

    INJm  = full_indx(INJloc);
    EXTm  = full_indx(EXTloc);

    INJFace = (radialdist <= radialedge)&(Nnx == 1);     
    EXTFace  = (radialdist <= radialedge)&(Nnx == Nx);
    
    disp([num2str(Mnempty),' grids within region of interest'])
    
    %% load variables
    for i = 1
        Phidist = 0; 
        pore_WH = 0;
        K = 0;
        Kx = 0;
        Ky = 0;
        Kz = 0;
        Correcta = 0;
        aporos = 0;
        afrac = 0;
        mx = 0;
        my = 0; 
        mz = 0;  
        IDm = 0;  
        H = 0;  
        Ny = 0;  
        IDmm1 = 0;  
        IDmp1 = 0;  
        IDmmNz = 0;  
        IDmpNz = 0;  
        IDmmNyNz = 0;  
        IDmpNyNz = 0; 
        MP1 = 0;   
        Injloc = 0;  
        L = 0;  
        Nx = 0;   
        Exitloc = 0;  
        NIDm = 0;  
        gridordinate = 0;  
        gridX = 0; 
        gridY = 0;  
        gridZ = 0; 
        Pguess = 0; 
        gridI = 0;
        qpref = 0;
        Pinj1 = 0;
    end

    %% Generate porosity and Permeability
    if ~exist('Datafile.mat', 'file') 
        recharacterized(q) 
    end
    load ('Datafile.mat') 
    

    minP = Pout;
    maxP = Pinj1;
    mnP = minP;
    mxP = minP + (maxP-minP)*qp/qpref;
    Pguess = interp1([minP, maxP], [mnP, mxP], Pguess);
    
    %% More variables
    C_HCl = Conc_HCl*rho_acid/(100*MW_HCl);
    Ef1   = Ef0*exp(-DeltaE_R/Temp);
    dxp   = Lp/Nx;
    dyp   = wp/Ny;
    dzp   = Hp/Nz;
    Vbp   = dxp*dyp*dzp;
    Axp   = dzp*dyp;
    Ayp   = dzp*dxp;
    Azp   = dxp*dyp;
    
    dxc   = Lc/Nx;
    dyc   = wc/Ny;
    dzc   = Hc/Nz;
    Vbc   = dxc*dyc*dzc;
    Axc   = dzc*dyc;
    Ayc   = dzc*dxc;
    Azc   = dxc*dyc;
    
    alphay  = zymultply*alphax;
    alphaz  = zymultply*alphax;
    
    S1    = Sm*rhoM;
    VB    = pi*wp*Hp*Lp/4;
    disp(['core volume (cc) = ', num2str(VB)])
    PoV   = VB*phio;
    disp(['pore volume (cc) = ', num2str(PoV)])
    dt    = PoV/qp/1000;
    
    %% Result matrix
    Volfrac       = zeros(M, 1);
    Volfrac(IDm)  = v1;
    ConcHCl       = zeros(M, 1);
    Pressure      = Pout*ones(M, 1);
    Pressure(IDm) = Pguess(1:end-1);
    Pinj          = Pguess(end);
    Porosity      = Phidist;
    Permeability  = K;
    
    MeanConc      = 0;
    MaxConcOut    = 0;
    Delta_P       = 0;
    
    Phi           = Phidist;
    K             = Permeability;
    qinj          = 0;   % Volume injected
    Kstim         = Ks;
    nphio         = phio;
    MeanPorosity  = nphio;
    MeanK         = @(K, m1, m2) harmmean([K(m1),K(m2)],2);
    MeanKfactor   = 4e3*qp*miu_acid*Lp/(pi*wp^2);
    
    
    %% Initialization  
    disp(['C_HCl = ', num2str(C_HCl)]);
    Max_Exit_poron  = max(Phidist(EXTloc));
    Max_Exit_poro   = Max_Exit_poron;
    MaxExitPoro     = Max_Exit_poro;
    Max_Exit_PoroL  = 0.55;
    disp(['Maximum exit Porosity = ',num2str(Max_Exit_poron)]);
    disp(['Maximum exit Porosity Limit = ',num2str(Max_Exit_PoroL)]);
    PVinj   = 0;
    disp(' ');
    rT      = 0;
    Time    = 0;
    n       = 1;
    np1     = 2;
    Storage = 1;
    

    %% Indecies for elements of the jacobian
    I_idx = [[IDm; IDm; IDm; IDm; IDm; IDm; IDm; ...    % dRpdp
              IDm; IDm; IDm; IDm; IDm; IDm; IDm; ...    % dRpdc
              IDm;   ...                                % dRpdv
              INJm]; ...                                % dRpdInjp

         M + [IDm; IDm; IDm; IDm; IDm; IDm; IDm; ...    % dRcdp
              IDm; IDm; IDm; IDm; IDm; IDm; IDm; ...    % dRcdc
              IDm;   ...                                % dRcdv
              INJm]; ...                                % dRcdInjp

       2*M + [IDm; IDm; IDm; IDm; IDm; IDm; IDm; ...    % dRvdp
              IDm; IDm; IDm; IDm; IDm; IDm; IDm; ...    % dRvdc
              IDm;   ...                                % dRvdv 
              INJm]; ...                                % dRvdInjp 

       3*M + ones(size(INJm));                          % dRidp
       3*M + 1];                                        % dRidInjp   

              
    J_idx = [[IDmmNyNz; IDmmNz; IDmm1; IDm; IDmp1; IDmpNz; IDmpNyNz; ....         % dRpdp
             M + [IDmmNyNz; IDmmNz; IDmm1; IDm; IDmp1; IDmpNz; IDmpNyNz]; ....    % dRpdc
             2*M + IDm;               ...                                         % dRpdv
             3*M + ones(size(INJm))]; ...                                         % dRpdInjp

             [IDmmNyNz; IDmmNz; IDmm1; IDm; IDmp1; IDmpNz; IDmpNyNz; ....         % dRcdp
             M + [IDmmNyNz; IDmmNz; IDmm1; IDm; IDmp1; IDmpNz; IDmpNyNz]; ....    % dRcdc
             2*M + IDm;               ...                                         % dRcdv
             3*M + ones(size(INJm))]; ...                                         % dRcdInjp

             [IDmmNyNz; IDmmNz; IDmm1; IDm; IDmp1; IDmpNz; IDmpNyNz; ....         % dRvdp
             M + [IDmmNyNz; IDmmNz; IDmm1; IDm; IDmp1; IDmpNz; IDmpNyNz]; ....    % dRvdc
             2*M + IDm;               ...                                         % dRvdv
             3*M + ones(size(INJm))]; ...                                         % dRvdInjp                           

             INJm;    ...                                                         % dRidp
             3*M + 1];                                                            % dRidInjp 

    Jremove = [NIDm; M + NIDm; 2*M + NIDm];
    n_Idx = numel(IDm);
    nIndx = 1:n_Idx;

    % Vectorised Residual function
    function [Rp, Rc, Rv, Ri] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                         CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                         Vm, Pinj)

        PmpNyNz(Exitloc) = Pout;      CmpNyNz(Exitloc) = 0;
        PmmNyNz(Injloc)  = Pinj;      CmmNyNz(Injloc)  = C_HCl;

        % Compute Porosity
        Phi(IDm) = Phidist(IDm) + (1 - Phidist(IDm)).*(v1 - Vm); 

        % Compute Cross Sectional for Local Diffusion Coefficient
        Ls  = (3*Phi(IDm)*Vbc/(4*pi)).^(1/3);    

        % Compute Concentration at the Reaction Surface
        Cs  = (kc./Ls).*Cm./((kc./Ls) + (1 - Phidist(IDm))*Vbc.*Vm*S1*Ef1);  

        % Compute Permeability (Low Perm and High Perm)
        K(IDm)   = 1e-6*Correcta*(exp(aporos*(Phi(IDm) - phio)).*(Phi(IDm) < pore_WH) + ...
                                  exp(afrac*(Phi(IDm) - phio)).*(Phi(IDm) >= pore_WH)); 
        K  = K.*(radialdist <= radialedge); 
        K(INJFace) = Inf;  K(EXTFace) = Inf; 
        
        Kx = K.*mx; Ky = K.*my; Kz = K.*mz;

  %%    %%%%%%%%%%%% Compute Bulk Flow Residual %%%%%%%%%%%%%%    %%
                  
        % Compute voidage rate
        DPhiDT = (1 - Phidist(IDm))*Vbp.*Vm*S1*Ef1.*Cs*(nM/nA)*(MW_CaCO3/rhoM);
                  % + Phidist(IDm)*Vbp*Comp*(P(IDm, n) - Pm)

        % Compute Mean Permeabilities
        KZmmm1    = MeanK(Kz, IDm, IDmm1);    
        KZmmp1    = MeanK(Kz, IDm, IDmp1);    
        KYmmmNz   = MeanK(Ky, IDm, IDmmNz);    
        KYmmpNz   = MeanK(Ky, IDm, IDmpNz);    
        KXmmmNyNz = MeanK(Kx, IDm, IDmmNyNz); 
        KXmmpNyNz = MeanK(Kx, IDm, IDmpNyNz); 
        
        % Compute Transmisibilities
        TZmmm1    = KZmmm1*Azp/(miu_acid*dzp);    
        TZmmp1    = KZmmp1*Azp/(miu_acid*dzp);    
        TYmmmNz   = KYmmmNz*Ayp/(miu_acid*dyp);    
        TYmmpNz   = KYmmpNz*Ayp/(miu_acid*dyp);    
        TXmmmNyNz = KXmmmNyNz*Axp/(miu_acid*dxp); 
        TXmmpNyNz = KXmmpNyNz*Axp/(miu_acid*dxp);  

        % Compute Pressure Difference
        DPmmm1    = Pm - Pmm1;    
        DPmmp1    = Pm - Pmp1;
        DPmmmNz   = Pm - PmmNz;   
        DPmmpNz   = Pm - PmpNz;
        DPmmmNyNz = Pm - PmmNyNz; 
        DPmmpNyNz = Pm - PmpNyNz;

        % Compute Fluxes
        JPmmm1    = TZmmm1.*(Pmm1 - Pm);
        JPmmp1    = TZmmp1.*(Pmp1 - Pm);
        JPmmmNz   = TYmmmNz.*(PmmNz - Pm);
        JPmmpNz   = TYmmpNz.*(PmpNz - Pm);
        JPmmmNyNz = TXmmmNyNz.*(PmmNyNz - Pm);
        JPmmpNyNz = TXmmpNyNz.*(PmpNyNz - Pm);

        % Sum Fluxes
        SumJP = JPmmmNyNz + JPmmmNz + JPmmm1 + JPmmp1 + JPmmpNz + JPmmpNyNz; 

        % Compute Residual
        Rp    = (SumJP - DPhiDT)*dt;

 %%  %%%%%%%%%%%%% Compute Reaction Residual %%%%%%%%%%%%%%  %%

        % Compute Velocities Using Darcy Law (m/s)
        Vzmmm1    = -0.01*KZmmm1.*DPmmm1/(miu_acid*dzp);    
        Vzmmp1    = -0.01*KZmmp1.*DPmmp1/(miu_acid*dzp); 
        VymmmNz   = -0.01*KYmmmNz.*DPmmmNz/(miu_acid*dyp);    
        VymmpNz   = -0.01*KYmmpNz.*DPmmpNz/(miu_acid*dyp); 
        VxmmmNyNz = -0.01*KXmmmNyNz.*DPmmmNyNz/(miu_acid*dxp); 
        VxmmpNyNz = -0.01*KXmmpNyNz.*DPmmpNyNz/(miu_acid*dxp);

        % Applying Upstream Weighting to Determine Flowing Concentration
        Cflowmmm1    = Cm.*(DPmmm1 >= 0) + Cmm1.*(DPmmm1 < 0); 
        Cflowmmp1    = Cm.*(DPmmp1 >= 0) + Cmp1.*(DPmmp1 < 0); 
        CflowmmmNz   = Cm.*(DPmmmNz >= 0) + CmmNz.*(DPmmmNz < 0); 
        CflowmmpNz   = Cm.*(DPmmpNz >= 0) + CmpNz.*(DPmmpNz < 0); 
        CflowmmmNyNz = Cm.*(DPmmmNyNz >= 0) + CmmNyNz.*(DPmmmNyNz < 0); 
        CflowmmpNyNz = Cm.*(DPmmpNyNz >= 0) + CmpNyNz.*(DPmmpNyNz < 0);

        % Compute Fluxes
        JCmmm1    = dt*Azc*((Cmm1 - Cm).*(Dm + alphaz*abs(Vzmmm1))/dzc + Cflowmmm1.*Vzmmm1);
        JCmmp1    = dt*Azc*((Cmp1 - Cm).*(Dm + alphaz*abs(Vzmmp1))/dzc + Cflowmmp1.*Vzmmp1);
        JCmmmNz   = dt*Ayc*((CmmNz - Cm).*(Dm + alphay*abs(VymmmNz))/dyc + CflowmmmNz.*VymmmNz);
        JCmmpNz   = dt*Ayc*((CmpNz - Cm).*(Dm + alphay*abs(VymmpNz))/dyc + CflowmmpNz.*VymmpNz);
        JCmmmNyNz = dt*Axc*((CmmNyNz - Cm).*(Dm + alphax*abs(VxmmmNyNz))/dxc + CflowmmmNyNz.*VxmmmNyNz);
        JCmmpNyNz = dt*Axc*((CmpNyNz - Cm).*(Dm + alphax*abs(VxmmpNyNz))/dxc + CflowmmpNyNz.*VxmmpNyNz);

        % Sum Fluxes
        SumJC = JCmmmNyNz + JCmmmNz + JCmmm1 + JCmmp1 + JCmmpNz + JCmmpNyNz; 

        % Compute Residuals
        Rc = SumJC - (1 - Phidist(IDm)).*Vm.*Cs*Vbc*S1*Ef1*dt - ...
                (Phi(IDm).*Cm - Porosity(IDm, n).*ConcHCl(IDm, n))*Vbc; 

%%  %%%%%%%%%%%%%%%%% Compute Mineral Volume Residuals %%%%%%%%%%%%%%%%%%%%%%%%%% %%
        Rv = (Vm - Volfrac(IDm, n)) + Vm.*Cs*(nM/nA)*(MW_CaCO3/rhoM)*S1*Ef1*dt;

%%  %%%%%%%%%%%%%%%%% Compute Injection Balance Residuals %%%%%%%%%%%%%%%%%%%%%%% %%
        Ri = TXmmmNyNz(Injloc)'*(Pinj - Pm(Injloc)) - qp;
    end

    % Vectorized Jacobian of the Residual
    function J = Jacobian(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                          CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                          Vm, Pinj, Rp, Rc, Rv, Ri)     

        TXmmmNyNz = MeanK(Kx, IDm, IDmmNyNz)*Axp/(miu_acid*dxp);
        %% Derivatives with respect to Pressures
        dl = 2e-8;
        
        % Perturb PmmNyNz
        del = dl*(abs(PmmNyNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz+del, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                       Vm, Pinj);
        dRpdPmmNyNz = (Rp_-Rp)./del;
        dRcdPmmNyNz = (Rc_-Rc)./del;
        dRvdPmmNyNz = (Rv_-Rv)./del;

        % Perturb PmmNz
        del = dl*(abs(PmmNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz+del, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdPmmNz = (Rp_-Rp)./del;
        dRcdPmmNz = (Rc_-Rc)./del;
        dRvdPmmNz = (Rv_-Rv)./del;

        % Perturb Pmm1
        del = dl*(abs(Pmm1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1+del, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdPmm1 = (Rp_-Rp)./del;
        dRcdPmm1 = (Rc_-Rc)./del;
        dRvdPmm1 = (Rv_-Rv)./del;

        % Perturb Pm
        del = dl*(abs(Pm) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm+del, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdPm = (Rp_-Rp)./del;
        dRcdPm = (Rc_-Rc)./del;
        dRvdPm = (Rv_-Rv)./del;

        % Perturb Pmp1
        del = dl*(abs(Pmp1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1+del, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdPmp1 = (Rp_-Rp)./del;
        dRcdPmp1 = (Rc_-Rc)./del;
        dRvdPmp1 = (Rv_-Rv)./del;

        % Perturb PmpNz
        del = dl*(abs(PmpNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz+del, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdPmpNz = (Rp_-Rp)./del;
        dRcdPmpNz = (Rc_-Rc)./del;
        dRvdPmpNz = (Rv_-Rv)./del;

        % Perturb PmpNyNz
        del = dl*(abs(PmpNyNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz+del, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                       Vm, Pinj);
        dRpdPmpNyNz = (Rp_-Rp)./del;
        dRcdPmpNyNz = (Rc_-Rc)./del;
        dRvdPmpNyNz = (Rv_-Rv)./del;


        %% Derivatives with respect to Concentration
        % Perturb CmmNyNz
        del = dl*(abs(CmmNyNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz+del, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                       Vm, Pinj);
        dRpdCmmNyNz = (Rp_-Rp)./del;
        dRcdCmmNyNz = (Rc_-Rc)./del;
        dRvdCmmNyNz = (Rv_-Rv)./del;

        % Perturb CmmNz
        del = dl*(abs(CmmNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz+del, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdCmmNz = (Rp_-Rp)./del;
        dRcdCmmNz = (Rc_-Rc)./del;
        dRvdCmmNz = (Rv_-Rv)./del;

        % Perturb Cmm1
        del = dl*(abs(Cmm1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1+del, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdCmm1 = (Rp_-Rp)./del;
        dRcdCmm1 = (Rc_-Rc)./del;
        dRvdCmm1 = (Rv_-Rv)./del;

        % Perturb Cm
        del = dl*(abs(Cm) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm+del, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdCm = (Rp_-Rp)./del;
        dRcdCm = (Rc_-Rc)./del;
        dRvdCm = (Rv_-Rv)./del;

        % Perturb Cmp1
        del = dl*(abs(Cmp1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1+del, CmpNz, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdCmp1 = (Rp_-Rp)./del;
        dRcdCmp1 = (Rc_-Rc)./del;
        dRvdCmp1 = (Rv_-Rv)./del;

        % Perturb CmpNz
        del = dl*(abs(CmpNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz+del, CmpNyNz, ...
                                      Vm, Pinj);
        dRpdCmpNz = (Rp_-Rp)./del;
        dRcdCmpNz = (Rc_-Rc)./del;
        dRvdCmpNz = (Rv_-Rv)./del;

        % Perturb CmpNyNz
        del = dl*(abs(CmpNyNz) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz+del, ...
                                       Vm, Pinj);
        dRpdCmpNyNz = (Rp_-Rp)./del;
        dRcdCmpNyNz = (Rc_-Rc)./del;
        dRvdCmpNyNz = (Rv_-Rv)./del;

        %% Derivatives with respect to Volume fraction
        % Perturb Vm
        del = dl*(abs(Vm) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                      CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                      Vm+del, Pinj);
        dRpdVm = (Rp_-Rp)./del;
        dRcdVm = (Rc_-Rc)./del;
        dRvdVm = (Rv_-Rv)./del;

        %% Derivatives with respect to injection Pressure
        % Perturb Pinj
        del = dl*(abs(Pinj) + 1);
        [Rp_, Rc_, Rv_, Ri_] = Residual(PmmNyNz, PmmNz, Pmm1, Pm, Pmp1, PmpNz, PmpNyNz, ...
                                        CmmNyNz, CmmNz, Cmm1, Cm, Cmp1, CmpNz, CmpNyNz, ...
                                        Vm, Pinj+del);

        dRpdPinj = (Rp_(Injloc)-Rp(Injloc))./del;
        dRcdPinj = (Rc_(Injloc)-Rc(Injloc))./del;
        dRvdPinj = (Rv_(Injloc)-Rv(Injloc))./del;
        dRidPinj = (Ri_-Ri)./del;

        V_val = [dRpdPmmNyNz; dRpdPmmNz; dRpdPmm1; dRpdPm; dRpdPmp1; dRpdPmpNz; dRpdPmpNyNz; ...
                 dRpdCmmNyNz; dRpdCmmNz; dRpdCmm1; dRpdCm; dRpdCmp1; dRpdCmpNz; dRpdCmpNyNz; ...
                 dRpdVm; dRpdPinj; ...
                 dRcdPmmNyNz; dRcdPmmNz; dRcdPmm1; dRcdPm; dRcdPmp1; dRcdPmpNz; dRcdPmpNyNz; ...
                 dRcdCmmNyNz; dRcdCmmNz; dRcdCmm1; dRcdCm; dRcdCmp1; dRcdCmpNz; dRcdCmpNyNz; ...
                 dRcdVm; dRcdPinj; ...
                 dRvdPmmNyNz; dRvdPmmNz; dRvdPmm1; dRvdPm; dRvdPmp1; dRvdPmpNz; dRvdPmpNyNz; ...
                 dRvdCmmNyNz; dRvdCmmNz; dRvdCmm1; dRvdCm; dRvdCmp1; dRvdCmpNz; dRvdCmpNyNz; ...
                 dRvdVm; dRvdPinj; ...
                 dRpdPinj; dRidPinj];
        if(any(isinf(V_val)))
            lat=1;
        end

        J = sparse(I_idx, J_idx, V_val, 3*M + 1, 3*M + 1, 50*M);
        J(Jremove, :) = []; J(:, Jremove) = [];
    end

    %conditioned = false;
    while (PVinj(end) < PVlimit)
        C  = ConcHCl(:, n); P = Pressure(:, n); V = Volfrac(:, n);

        iter = 0; maxIter = 10; isconverged = 0;
        RelTol = 1e-4; AbsTol = 1e-6;
        while(~isconverged && iter < maxIter)

            % Extract variables
            PmmNyNz_n = P(IDmmNyNz); PmmNz_n = P(IDmmNz); Pmm1_n = P(IDmm1); 
            PmpNyNz_n = P(IDmpNyNz); PmpNz_n = P(IDmpNz); Pmp1_n = P(IDmp1); 

            CmmNyNz_n = C(IDmmNyNz); CmmNz_n = C(IDmmNz); Cmm1_n = C(IDmm1); 
            CmpNyNz_n = C(IDmpNyNz); CmpNz_n = C(IDmpNz); Cmp1_n = C(IDmp1); 

            Pm_n = P(IDm); Cm_n = C(IDm); Vm_n = V(IDm); Pinj_n = Pinj;

           
            % Update the residual
            [Rp, Rc, Rv, RinjP] = Residual(PmmNyNz_n, PmmNz_n, Pmm1_n, Pm_n, Pmp1_n, PmpNz_n, PmpNyNz_n, ...
                                           CmmNyNz_n, CmmNz_n, Cmm1_n, Cm_n, Cmp1_n, CmpNz_n, CmpNyNz_n, ...
                                           Vm_n, Pinj_n);
            
            % Evaluate the Jacobian
            J = Jacobian(PmmNyNz_n, PmmNz_n, Pmm1_n, Pm_n, Pmp1_n, PmpNz_n, PmpNyNz_n, ...
                         CmmNyNz_n, CmmNz_n, Cmm1_n, Cm_n, Cmp1_n, CmpNz_n, CmpNyNz_n, ...
                         Vm_n, Pinj_n, Rp, Rc, Rv, RinjP);

            warning('') % Clear last warning message
            % Compute Newton Step
            r = [Rp; Rc; Rv; RinjP];

            dX = -(J\r);
            % if(~conditioned)
            %     [L,U] = ilu(J);
            %     conditioned = true;
            % end
            % dX = -cgs(J, r, 1e-6, 1000, L, U, 0*r);
            [~, warnId] = lastwarn;
            if(strcmp(warnId, 'MATLAB:singularMatrix') || ...
                    strcmp(warnId, 'MATLAB:nearlySingularMatrix'))
                break;
            end

            % Extract Pressure, Concentration and Volume Step
            dP = dX(0*n_Idx + nIndx);
            dC = dX(1*n_Idx + nIndx);
            dV = dX(2*n_Idx + nIndx);
            dPinj = dX(3*n_Idx + 1);

            % Update Pressure, Concentration and Volume
            P(IDm) = P(IDm) + dP; 
            reldP = relchange(P(IDm), dP, AbsTol);
            if(max(dP)>0.1)
                break;
            end

            C(IDm) = min(C_HCl, max(0, C(IDm) + dC)); 
            reldC = relchange(C(IDm), dC, AbsTol);
           
            V(IDm) = min(v1, max(0, V(IDm) + dV)); 
            reldV = relchange(V(IDm), dV, AbsTol);
            Pinj = Pinj + dPinj; 
            P(NIDm) = Pinj;

            error = norm(r); 
            tolcheck = max([reldP, reldC, reldV]) < RelTol;
            isconverged = error < AbsTol || tolcheck; 
            iter = iter + 1;
        end

        if(~isconverged)
            dt = dt/2
            Clock = datetime;
            disp(Clock)
            continue;
        else
            % log in solutions
            Time(np1)            = Time(n) + dt; 
            Pressure(:, np1)     = P;
            ConcHCl(:, np1)      = C;
            Volfrac(:, np1)      = V;
            Porosity(:, np1)     = Phi;
            Permeability(:, np1) = K;

            MeanConc(np1)        = (C(IDm)'*Phi(IDm))/sum(Phi(IDm));      % Mean concentration in the core
            MaxConcOut(np1)      = max(C(EXTloc));                        % Maximum concentration flowing out at the exit end
            Delta_P(np1)         = Pinj - Pout;                           % Change in pressure accross the core
            Kmean                = MeanKfactor/Delta_P(np1);              % Mean permeability
            MeanPorosity(np1)    = mean(Phi(IDm));                        % Mean porosity
            qinj(np1)            = 80.5196416*Time(np1)*qc;
            Kstim(np1)           = Kmean;
            PVinj(np1)           = Time(np1)*qp/PoV;
            Max_Exit_poro        = max(Phi(EXTloc));
            MaxExitPoro(np1)     = Max_Exit_poro;
            
            disp(char(["Maximum exit Porosity = " + Max_Exit_poro
                       "Mean concentration = " + MeanConc(np1) + " mol/cubic meter"
                       "Maximum outflow concentration = " + MaxConcOut(np1) + " mol/cubic meter"
                       "Pore Volume Injected = " + PVinj(np1)
                       "Injection Pressure = " + Pinj + " atm"
                       "Pressure difference = " + Delta_P(np1) + " atm"
                       "Mean Permebility = " + Kstim(np1) + " darcy"
                       "Time = " + Time(np1)+ "sec"
                       "Flowrate = " + num2str(q) + "cc/min"]));

            if(iter < 4)
                disp('Adjusting time step to optimize to speed up computation')
                dt = 2*dt;
            end
            
            Qdisp = [num2str(q),'cc/min ============='];
            disp('==========================================================');
            disp(['===================== q = ',Qdisp(1:17),'===============']);
            disp('==========================================================');
            
            Clock = datetime;
            disp(Clock)
            disp(' ');
            disp(' ');
            n = np1;
            np1 = n + 1;
            if np1 > Storagecap + 2
                FileName1x  = [FileName1,num2str(Storage)];
                FileName2x  = [FileName2,num2str(Storage)];
                FileName3x  = [FileName3,num2str(Storage)];
                FileName4x  = [FileName4,num2str(Storage)];
                
                ConcHClx = ConcHCl(:,1:Storagecap); 
                ConcHCl(:,1:Storagecap) = [];
                Volfracx = Volfrac(:,1:Storagecap); 
                Volfrac(:,1:Storagecap) = [];
                Porosityx = Porosity(:,1:Storagecap); 
                Porosity(:,1:Storagecap) = [];
                
                
                Kstimx  = Kstim(:,1:Storagecap); 
                Kstim(:,1:Storagecap) = [];
                PVinjx  = PVinj(:,1:Storagecap); 
                PVinj(:,1:Storagecap) = [];
                Delta_Px  = Delta_P(:,1:Storagecap); 
                Delta_P(:,1:Storagecap) = [];
                MeanConcx = MeanConc(:,1:Storagecap); 
                MeanConc(:,1:Storagecap) = [];
                MaxConcOutx = MaxConcOut(:,1:Storagecap); 
                MaxConcOut(:,1:Storagecap) = [];
                Timex  = Time(:,1:Storagecap); 
                Time(:,1:Storagecap) = [];
                
                save (FileName1x, 'ConcHClx')    
                save (FileName2x, 'Porosityx')
                save (FileName3x, 'Volfracx')
                save (FileName4x, 'Kstimx', 'PVinjx', 'Delta_Px', ...
                          'MeanConcx', 'MaxConcOutx', 'Timex')        
                
                clear ConcHClx Porosityx Volfracx Kstimx PVinjx Delta_Px ...
                    MeanConcx MaxConcOutx Timex
                
                np1 = np1 - Storagecap;
                n   = n   - Storagecap;
                
                Storage = Storage + 1;
                save("RestartFile.mat")
            end 
        end
    
        if(mod(np1, 100)==0)
            save("RestartFile.mat")
        end

        if MaxExitPoro(n) > Max_Exit_PoroL
            break;
        end 
    end
    
    PVfinal = PVinj(end);
    Nt = (Storage - 1)*Storagecap + np1;
    FileName0   = ['support_',ID(1:end-1)];
    FileName1x  = [FileName1,num2str(Storage)];
    FileName2x  = [FileName2,num2str(Storage)];
    FileName3x  = [FileName3,num2str(Storage)];
    FileName4x  = [FileName4,num2str(Storage)];

    save (FileName0, 'dxc', 'dyc', 'dzc', 'M','v1', 'Max_Exit_PoroL', ...
                     'C_HCl', 'H', 'w', 'L', 'pore_WH', 'Nx', 'Ny', ...
                     'Nz', 'gridordinate', 'Nt', 'Storage', 'Storagecap')
    
    ConcHClx    = ConcHCl;
    Volfracx    = Volfrac;
    Porosityx   = Porosity;
    Kstimx      = Kstim;
    PVinjx      = PVinj;
    Delta_Px    = Delta_P;
    MeanConcx   = MeanConc;
    MaxConcOutx = MaxConcOut;
    Timex       = Time;
    
    save (FileName1x, 'ConcHClx')    
    save (FileName2x, 'Porosityx')
    save (FileName3x, 'Volfracx')
    save (FileName4x, 'Kstimx', 'PVinjx', 'Delta_Px', 'MeanConcx', ...
                      'MaxConcOutx', 'Timex')     

    delete RestartFile.mat
end 

function rdel = relchange(X, dX, Abstol)
    I = abs(X)>Abstol;
    rdel = max(abs(dX(I)./X(I)));
end
