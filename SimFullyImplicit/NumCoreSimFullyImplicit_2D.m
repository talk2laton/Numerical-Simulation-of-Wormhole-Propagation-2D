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
        IDm = 0;  
        H = 0;  
        Ny = 0;  
        IDmm1 = 0;  
        IDmp1 = 0;  
        IDmmNy = 0;  
        IDmpNy = 0; 
        MP1 = 0;   
        Injloc = 0;  
        L = 0;  
        Nx = 0;   
        Exitloc = 0;  
        NIDm = 0;  
        gridordinate = 0;  
        gridX = 0; 
        gridY = 0;  
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
    dzp   = dyp;
    Vbp   = dxp*dyp*dzp;
    Axp   = dzp*dyp;
    Ayp   = dzp*dxp;
    
    dxc   = Lc/Nx;
    dyc   = wc/Ny;
    dzc   = dyc;
    Vbc   = dxc*dyc*dzc;
    Axc   = dzc*dyc;
    Ayc   = dzc*dxc;
    
    alphay  = zymultply*alphax;
    
    S1    = Sm*rhoM;
    VB    = wp*dzp*Lp;
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
    MeanKfactor   = 1e3*qp*miu_acid*Lp/(wp*dzp);
    
    
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
    I_idx = [[IDm; IDm; IDm; IDm; IDm; ...    % dRpdp
              IDm; IDm; IDm; IDm; IDm; ...    % dRpdc
              IDm;   ...                      % dRpdv
              INJm]; ...                      % dRpdInjp

         M + [IDm; IDm; IDm; IDm; IDm; ...    % dRcdp
              IDm; IDm; IDm; IDm; IDm; ...    % dRcdc
              IDm;   ...                      % dRcdv
              INJm]; ...                      % dRcdInjp

       2*M + [IDm; IDm; IDm; IDm; IDm; ...    % dRvdp
              IDm; IDm; IDm; IDm; IDm; ...    % dRvdc
              IDm;   ...                      % dRvdv 
              INJm]; ...                      % dRvdInjp 

       3*M + ones(size(INJm));                % dRidp
       3*M + 1];                              % dRidInjp   

              
    J_idx = [[IDmmNy; IDmm1; IDm; IDmp1; IDmpNy; ....         % dRpdp
             M + [IDmmNy; IDmm1; IDm; IDmp1; IDmpNy]; ....    % dRpdc
             2*M + IDm;               ...                     % dRpdv
             3*M + ones(size(INJm))]; ...                     % dRpdInjp

             [IDmmNy; IDmm1; IDm; IDmp1; IDmpNy; ....         % dRcdp
             M + [IDmmNy; IDmm1; IDm; IDmp1; IDmpNy]; ....    % dRcdc
             2*M + IDm;               ...                     % dRcdv
             3*M + ones(size(INJm))]; ...                     % dRcdInjp

             [IDmmNy; IDmm1; IDm; IDmp1; IDmpNy; ....         % dRvdp
             M + [IDmmNy; IDmm1; IDm; IDmp1; IDmpNy]; ....    % dRvdc
             2*M + IDm;               ...                     % dRvdv
             3*M + ones(size(INJm))]; ...                     % dRvdInjp                           

             INJm;    ...                                     % dRidp
             3*M + 1];                                        % dRidInjp 

    Jremove = [NIDm; M + NIDm; 2*M + NIDm];
    n_Idx = numel(IDm);
    nIndx = 1:n_Idx;

    % Vectorised Residual function
    function [Rp, Rc, Rv, Ri] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                         CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                         Vm, Pinj)

        PmpNy(Exitloc) = Pout;      CmpNy(Exitloc) = 0;
        PmmNy(Injloc)  = Pinj;      CmmNy(Injloc)  = C_HCl;

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
        
        Kx = K.*mx; Ky = K.*my;

  %%    %%%%%%%%%%%% Compute Bulk Flow Residual %%%%%%%%%%%%%%    %%
                  
        % Compute voidage rate
        DPhiDT = (1 - Phidist(IDm))*Vbp.*Vm*S1*Ef1.*Cs*(nM/nA)*(MW_CaCO3/rhoM);
                  % + Phidist(IDm)*Vbp*Comp*(P(IDm, n) - Pm)

        % Compute Mean Permeabilities  
        KYmmm1  = MeanK(Ky, IDm, IDmm1);    
        KYmmp1  = MeanK(Ky, IDm, IDmp1);    
        KXmmmNy = MeanK(Kx, IDm, IDmmNy); 
        KXmmpNy = MeanK(Kx, IDm, IDmpNy); 
        
        % Compute Transmisibilities  
        TYmmm1  = KYmmm1*Ayp/(miu_acid*dyp);    
        TYmmp1  = KYmmp1*Ayp/(miu_acid*dyp);    
        TXmmmNy = KXmmmNy*Axp/(miu_acid*dxp); 
        TXmmpNy = KXmmpNy*Axp/(miu_acid*dxp);  

        % Compute Pressure Difference
        DPmmm1  = Pm - Pmm1;    
        DPmmp1  = Pm - Pmp1;
        DPmmmNy = Pm - PmmNy; 
        DPmmpNy = Pm - PmpNy;

        % Compute Fluxes
        JPmmm1  = TYmmm1.*(Pmm1 - Pm);
        JPmmp1  = TYmmp1.*(Pmp1 - Pm);
        JPmmmNy = TXmmmNy.*(PmmNy - Pm);
        JPmmpNy = TXmmpNy.*(PmpNy - Pm);

        % Sum Fluxes
        SumJP = JPmmmNy + JPmmm1 + JPmmp1 + JPmmpNy; 

        % Compute Residual
        Rp    = (SumJP - DPhiDT)*dt;

 %%  %%%%%%%%%%%%% Compute Reaction Residual %%%%%%%%%%%%%%  %%

        % Compute Velocities Using Darcy Law (m/s)
        Vymmm1  = -0.01*KYmmm1.*DPmmm1/(miu_acid*dyp);    
        Vymmp1  = -0.01*KYmmp1.*DPmmp1/(miu_acid*dyp); 
        VxmmmNy = -0.01*KXmmmNy.*DPmmmNy/(miu_acid*dxp); 
        VxmmpNy = -0.01*KXmmpNy.*DPmmpNy/(miu_acid*dxp);

        % Applying Upstream Weighting to Determine Flowing Concentration
        Cflowmmm1  = Cm.*(DPmmm1 >= 0) + Cmm1.*(DPmmm1 < 0); 
        Cflowmmp1  = Cm.*(DPmmp1 >= 0) + Cmp1.*(DPmmp1 < 0); 
        CflowmmmNy = Cm.*(DPmmmNy >= 0) + CmmNy.*(DPmmmNy < 0); 
        CflowmmpNy = Cm.*(DPmmpNy >= 0) + CmpNy.*(DPmmpNy < 0); 

        % Compute Fluxes
        JCmmm1  = dt*Ayc*((Cmm1 - Cm).*(Dm + alphay*abs(Vymmm1))/dyc + Cflowmmm1.*Vymmm1);
        JCmmp1  = dt*Ayc*((Cmp1 - Cm).*(Dm + alphay*abs(Vymmp1))/dyc + Cflowmmp1.*Vymmp1);
        JCmmmNy = dt*Axc*((CmmNy - Cm).*(Dm + alphax*abs(VxmmmNy))/dxc + CflowmmmNy.*VxmmmNy);
        JCmmpNy = dt*Axc*((CmpNy - Cm).*(Dm + alphax*abs(VxmmpNy))/dxc + CflowmmpNy.*VxmmpNy);

        % Sum Fluxes
        SumJC = JCmmmNy + JCmmm1 + JCmmp1 + JCmmpNy; 

        % Compute Residuals
        Rc = SumJC - (1 - Phidist(IDm)).*Vm.*Cs*Vbc*S1*Ef1*dt - ...
                (Phi(IDm).*Cm - Porosity(IDm, n).*ConcHCl(IDm, n))*Vbc; 

%%  %%%%%%%%%%%%%%%%% Compute Mineral Volume Residuals %%%%%%%%%%%%%%%%%%%%%%%%%% %%
        Rv = (Vm - Volfrac(IDm, n)) + Vm.*Cs*(nM/nA)*(MW_CaCO3/rhoM)*S1*Ef1*dt;

%%  %%%%%%%%%%%%%%%%% Compute Injection Balance Residuals %%%%%%%%%%%%%%%%%%%%%%% %%
        Ri = TXmmmNy(Injloc)'*(Pinj - Pm(Injloc)) - qp;
    end

    % Vectorized Jacobian of the Residual
    function J = Jacobian(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                          CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                          Vm, Pinj, Rp, Rc, Rv, Ri)     

        %% Derivatives with respect to Pressures
        dl = 2e-8;
        

        % Perturb PmmNy
        del = dl*(abs(PmmNy) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy+del, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdPmmNy = (Rp_-Rp)./del;
        dRcdPmmNy = (Rc_-Rc)./del;
        dRvdPmmNy = (Rv_-Rv)./del;

        % Perturb Pmm1
        del = dl*(abs(Pmm1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1+del, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdPmm1 = (Rp_-Rp)./del;
        dRcdPmm1 = (Rc_-Rc)./del;
        dRvdPmm1 = (Rv_-Rv)./del;

        % Perturb Pm
        del = dl*(abs(Pm) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm+del, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdPm = (Rp_-Rp)./del;
        dRcdPm = (Rc_-Rc)./del;
        dRvdPm = (Rv_-Rv)./del;

        % Perturb Pmp1
        del = dl*(abs(Pmp1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1+del, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdPmp1 = (Rp_-Rp)./del;
        dRcdPmp1 = (Rc_-Rc)./del;
        dRvdPmp1 = (Rv_-Rv)./del;

        % Perturb PmpNy
        del = dl*(abs(PmpNy) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy+del, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdPmpNy = (Rp_-Rp)./del;
        dRcdPmpNy = (Rc_-Rc)./del;
        dRvdPmpNy = (Rv_-Rv)./del;


        %% Derivatives with respect to Concentration
        % Perturb CmmNy
        del = dl*(abs(CmmNy) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy+del, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdCmmNy = (Rp_-Rp)./del;
        dRcdCmmNy = (Rc_-Rc)./del;
        dRvdCmmNy = (Rv_-Rv)./del;

        % Perturb Cmm1
        del = dl*(abs(Cmm1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1+del, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdCmm1 = (Rp_-Rp)./del;
        dRcdCmm1 = (Rc_-Rc)./del;
        dRvdCmm1 = (Rv_-Rv)./del;

        % Perturb Cm
        del = dl*(abs(Cm) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm+del, Cmp1, CmpNy, ...
                                      Vm, Pinj);
        dRpdCm = (Rp_-Rp)./del;
        dRcdCm = (Rc_-Rc)./del;
        dRvdCm = (Rv_-Rv)./del;

        % Perturb Cmp1
        del = dl*(abs(Cmp1) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1+del, CmpNy, ...
                                      Vm, Pinj);
        dRpdCmp1 = (Rp_-Rp)./del;
        dRcdCmp1 = (Rc_-Rc)./del;
        dRvdCmp1 = (Rv_-Rv)./del;

        % Perturb CmpNy
        del = dl*(abs(CmpNy) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy+del, ...
                                      Vm, Pinj);
        dRpdCmpNy = (Rp_-Rp)./del;
        dRcdCmpNy = (Rc_-Rc)./del;
        dRvdCmpNy = (Rv_-Rv)./del;

        %% Derivatives with respect to Volume fraction
        % Perturb Vm
        del = dl*(abs(Vm) + 1);
        [Rp_, Rc_, Rv_, ~] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm+del, Pinj);
        dRpdVm = (Rp_-Rp)./del;
        dRcdVm = (Rc_-Rc)./del;
        dRvdVm = (Rv_-Rv)./del;

        %% Derivatives with respect to injection Pressure
        % Perturb Pinj
        del = dl*(abs(Pinj) + 1);
        [Rp_, Rc_, Rv_, Ri_] = Residual(PmmNy, Pmm1, Pm, Pmp1, PmpNy, ...
                                      CmmNy, Cmm1, Cm, Cmp1, CmpNy, ...
                                      Vm, Pinj+del);

        dRpdPinj = (Rp_(Injloc)-Rp(Injloc))./del;
        dRcdPinj = (Rc_(Injloc)-Rc(Injloc))./del;
        dRvdPinj = (Rv_(Injloc)-Rv(Injloc))./del;
        dRidPinj = (Ri_-Ri)./del;

        V_val = [dRpdPmmNy; dRpdPmm1; dRpdPm; dRpdPmp1; dRpdPmpNy; ...
                 dRpdCmmNy; dRpdCmm1; dRpdCm; dRpdCmp1; dRpdCmpNy; ...
                 dRpdVm; dRpdPinj; ...
                 dRcdPmmNy; dRcdPmm1; dRcdPm; dRcdPmp1; dRcdPmpNy; ...
                 dRcdCmmNy; dRcdCmm1; dRcdCm; dRcdCmp1; dRcdCmpNy; ...
                 dRcdVm; dRcdPinj; ...
                 dRvdPmmNy; dRvdPmm1; dRvdPm; dRvdPmp1; dRvdPmpNy; ...
                 dRvdCmmNy; dRvdCmm1; dRvdCm; dRvdCmp1; dRvdCmpNy; ...
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
            PmmNy_n = P(IDmmNy); Pmm1_n = P(IDmm1); 
            PmpNy_n = P(IDmpNy); Pmp1_n = P(IDmp1); 

            CmmNy_n = C(IDmmNy); Cmm1_n = C(IDmm1); 
            CmpNy_n = C(IDmpNy); Cmp1_n = C(IDmp1); 

            Pm_n = P(IDm); Cm_n = C(IDm); Vm_n = V(IDm); Pinj_n = Pinj;

           
            % Update the residual
            [Rp, Rc, Rv, RinjP] = Residual(PmmNy_n, Pmm1_n, Pm_n, Pmp1_n, PmpNy_n, ...
                                           CmmNy_n, Cmm1_n, Cm_n, Cmp1_n, CmpNy_n, ...
                                           Vm_n, Pinj_n);
            
            % Evaluate the Jacobian
            J = Jacobian(PmmNy_n, Pmm1_n, Pm_n, Pmp1_n, PmpNy_n, ...
                         CmmNy_n, Cmm1_n, Cm_n, Cmp1_n, CmpNy_n, ...
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
                     'gridordinate', 'Nt', 'Storage', 'Storagecap')
    
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
