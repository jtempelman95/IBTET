classdef VI_data_processing < matlab.mixin.Copyable %& VI_LatticeSys % handle% < matlab.mixin.Copyable

properties
    t 
    x
    v
    tc
    xc
    vc
    FVIc
    %spectral
    sys
  %  excitation
    Eband
    Ebandc
    wavelet_params
    decomp_params
    entropydata
    waveletdata
    wbands
    dispersion
    energy
    FFTdata
    FVI
end

methods  
    % ------------------------------------------
    function self = VI_data_processing(t,x,v)
    % ------------------------------------------
    % Intiize the data set
    % ------------------------------------------
        self.t = t;
        self.x = x;
        self.v = v;
        self.wavelet_params = wavelet_parameters();
        self.decomp_params.wavelet_params = self.wavelet_params;
    end

    % ------------------------------------------
    function compress(self,subsamp)
    % ------------------------------------------
    % Create subsampled arrays of data
    % ------------------------------------------
    if nargin <2
        self.tc = self.t(1:10:end);
        self.xc = self.x(1:10:end,:);
        self.vc = self.v(1:10:end,:);
            if ~isempty(self.FVI)
                self.FVIc = self.FVI(1:10:end);
            end
            if ~isempty(self.Eband)
                for k = 1:size(self.Eband,3)
                 self.Ebandc(:,:,k) = self.Eband(1:40:end,:,k);
                end
            end
    else
        self.tc = self.t(1:subsamp:end,:);
        self.xc = self.x(1:subsamp:end,:);
        self.vc = self.v(1:subsamp:end,:);
            if ~isempty(self.FVI)
                self.FVIc = self.FVI(1:subsamp:end);
            end
            if ~isempty(self.Eband)
                for k = 1:size(self.Eband,3)
                 self.Ebandc(:,:,k) = self.Eband(1:subsamp:end,:,k);
                end
            end
    end
    end

    % ------------------------------------------
    function compressEband(self)
    % ------------------------------------------
                 self.Ebandc = self.Eband(1:40:end,:,:);
    end

    % ------------------------------------------
    function clearTS(self)
    % ------------------------------------------
    % Clear time series
    % ------------------------------------------
    self.t = [];
    self.x = [];
    self.v = [];
    self.Eband = [];
    if ~isempty(self.FVI),self.FVI =[];end
    end


    % ------------------------------------------
    function Process(self)
    % ------------------------------------------
    % Perform all data processing
    % ------------------------------------------
        if isempty(self.wavelet_params)
        self.wavelet_params.Fs   =  pi                          ;
        self.wavelet_params.fo   =  2                           ;
        self.wavelet_params.fl   =  0                           ;
        self.wavelet_params.fu   =  pi                          ;
        self.wavelet_params.fn   =  300                         ;
        end
        if isempty(self.decomp_params)
            self.decomp_params.wavelet_params.fo = self.wavelet_params.fo;
            self.decomp_params.wavelet_params.fl = self.wavelet_params.fl;
            self.decomp_params.wavelet_params.fu = self.wavelet_params.fu;
            self.decomp_params.wavelet_params.fn = self.wavelet_params.fn;
            self.decomp_params.wavelet_params.Fs = self.wavelet_params.Fs;
        end
        
            % ------------------------------------------
            % Get Entropy Data
            % ------------------------------------------
            if isempty(self.entropydata)
                samps                   = round(linspace(1,numel(self.t),100));
                [PEntropy,PEinst]       = entropy(self,self.v,samps);
                self.entropydata.PE1    = PEntropy;
                self.entropydata.PE2    = PEinst;
            end
            % ------------------------------------------
            % Get Wavelet Data
            % ------------------------------------------
            if isempty(self.waveletdata)
                tstamps                         = round(linspace(1,numel(self.t),4));
                xinds                           = 1:size(self.x,2);
                [wavelet,frq]                   = self.wavelet(self.v(tstamps,xinds)');
                self.waveletdata.timestamps     = self.t(tstamps);
                self.waveletdata.frq            = frq;
                self.waveletdata.xinds          = xinds;
                self.waveletdata.wdata          = wavelet;
            end
            % ------------------------------------------
            % Get Wavenumber Data
            % ------------------------------------------
            if isempty(self.wbands)
                self.wbands = self.K_evolution(self.v);
            end
            % ------------------------------------------
            % Get Numerical Dispersion Data
            % ------------------------------------------
            if isempty(self.dispersion)
                    [mat,omega,kappa] = self.Fdomain();
                    self.dispersion.mat = mat;
                    self.dispersion.kappa = kappa;
                    self.dispersion.omega = omega;
            end
            % ------------------------------------------
            % Get Linear Eneergy
            % ------------------------------------------
            if isempty(self.energy)
                [E, EPE, EKE, TE] = linenergy(self);
                self.energy.E = E   ;
                self.energy.PE = EPE;
                self.energy.KE = EKE;
                self.energy.TE = TE ;
            end
            % ------------------------------------------
            % Get Truncated Spectrum 
            % ------------------------------------------
            if isempty(self.FFTdata)
                [fftmat,frqs,idx] = FFTs(self);
                self.FFTdata.idx = idx;
                self.FFTdata.frqs = frqs;
                self.FFTdata.fftmat = fftmat;
            end
            % ------------------------------------------
            % Get VI force
            % ------------------------------------------
            if isempty(self.FVI)
                self.FVI = self.FNL(self);
            end
            % ------------------------------------------
            % Get Energy Bands
            % ------------------------------------------
            if isempty(self.Eband)
                self.Eband = self.BandEnergy();
            end
    end

    % ------------------------------------------
    function Eband = BandEnergy(self)
    % ------------------------------------------
    % Compute the spatio-temporal energy on each band
    if ~isempty(self.t)
    tft     = self.t;
    else
    tft     = self.tc;
    end

    if ~isempty(self.v)
    ftdata  = self.v;
    else
    ftdata  = self.vc;
    end

    % Fourier domain properies
    Lt = size(ftdata,1);
    Lx = size(ftdata,1);
    dt = mean(diff(tft));
    dF = 1/(dt*Lt);
    F = (-(Lt/2)*dF):dF:(Lt/2)*dF;F = F(1:Lt);
    k = linspace(-2*pi,2*pi,Lx);
    [DISP,wavevec] = self.sys.dispersion;
    disp = (fft2(ftdata))           ;
    sz = size(disp)      ;
    Euse = [];
    for k = 1:size(DISP,2)+2

    if k <= size(DISP,2)
        optfrqs = [min(min(DISP(:,k))) max(max(DISP(:,k)))  ]/2/pi;
        E1 =[find( F<-optfrqs(1) & F>-optfrqs(2)) find( F>optfrqs(1) & F<optfrqs(2))];
    end

    % Data to pull
    if k == size(DISP,2)+1
        optfrqs = [max(max(DISP(:,size(DISP,2))))  ]/2/pi;
        E1 =[find( F<-optfrqs(1)) find( F>optfrqs(1) )];
    end

    if k == size(DISP,2)+2
       E1= find(~ismember(1:sz(1),Euse)) ;
    end

        dispfilt  = zeros(size(disp))                       ; 
        dispshift = fftshift(disp)                          ; 
        dispfilt(E1,:) = (dispshift(E1,:) )                 ; 
        org = ifft2(ifftshift(dispfilt))                    ; 
        Eband(:,:,k) = real(org);
        Euse = [Euse E1];
    end
    end



    
    % ------------------------------------------
    function [fftmat,frqs,idx] = FFTs(self,LTS,S)
    % ------------------------------------------
   
    % Length of time series to compute fft over
    if nargin<2
        LTS = round(numel(self.t)/2);
    end

    % S is the 'sites' to compute fft over
    if nargin < 3
        S = [31:36];
        asites = self.sys.ndof;
        S = [S self.sys.ndof+[1:6]-3-3*42];
    end

    Fs = 1/mean(diff(self.t));
    L = LTS;                 %(Length of time series)
    F = Fs*(0:(L-1))/L;      % frequency vector
    F = F * 2* pi;
    
    Fc = max(max(self.sys.dispersion))*2;
    for k = 1:length(S)
            if S(k) < self.sys.ndof/2
                D = self.v(1:LTS,S(k)) ;
            else
                D = self.v(end-(LTS-1):end,S(k)) ;
            end
            fdata =  abs(fft(D));
            fftmat(:,k) = fdata(F<Fc);
    end

    frqs = F(F<Fc);
    idx = S;
    end

    % ------------------------------------------
    function [mat,omega,kappa] = Fdomain(self,y,maxF)
    % ------------------------------------------
        % Compute the space-time fft of the data
        % ------------------------------------------
        if ~isempty(self.v)
            if nargin < 3
                maxF = 1.25*max(self.sys.Dispersion.DISP(:));
            end
            if nargin < 2
                y = self.v;
            end
            ndof = self.sys.ndof;
            t    = self.t;
            dt = mean(diff(t));
        else
            if nargin < 3
                maxF = 1.25*max(self.sys.Dispersion.DISP(:));
            end
            if nargin < 2
                y = self.vc;
            end
            ndof = self.sys.ndof;
            t    = self.tc;
            dt = mean(diff(t));
        end
            

        dy = 1;
        if self.sys.config == 3; dy = 3/2; end
        if self.sys.config == 4; dy = 3/2; end

        Lt = size(y,1);
        Ly = size(y,2);
        dF = 1/(dt*Lt);
        dk = 1/(dy*Ly);
        fftmat = fftshift(fft2(y));
        fftmat = fftmat(round(Lt/2:end),round(Ly/2:end));
        F = 0:dF:(Lt/2)*dF;
        k = 0:dk:(Ly/2)*dk;
        k = linspace(0,pi*2,numel(k));
        if self.sys.config ==3; 
        k = linspace(-pi,pi*2,numel(k));end
        if self.sys.config ==4; 
        k = linspace(-pi,pi,numel(k));end
        omega = 2*pi*F;
        kappa = k;
        idx = find(omega<maxF);
        mat = flip(fftmat(idx,:),2);
        omega = omega(idx);
    end


    % ------------------------------------------
    function WnumMat = K_evolution(self,y,samps)
    % ------------------------------------------
    % Compute the evolution of wavenumber
    % ------------------------------------------
    if nargin < 3
        samps = round(linspace(1,numel(self.t),100));
    end
    if nargin < 2
        y = self.v(:,1:2:self.sys.ndof);
    end
        r = 1;
        WnumMat = zeros(size(y,2)/2,numel(samps));
        for k = samps
            fdata   = (fft(y(k,:)));
            WnumMat(:,r) = fdata(1:end/2);
            r = r + 1;
        end
    end
    



    % ------------------------------------------
    function delta = DeltaX(self)
    % ------------------------------------------
    % Compute the VI deflection 
    % ------------------------------------------
        VIlocs = self.sys.VIsites;
        X1 = VIlocs;
        X2 = VIlocs - 1;
        delta = self.x(:,X2) - self.x(:,X1);
    end

    % ------------------------------------------
    function delta = DeltaV(self)
    % ------------------------------------------
    % Compute the VI deflection velocity
    % ------------------------------------------
        VIlocs = self.sys.VIsites;
        X1 = VIlocs+1;
        X2 = VIlocs;
        delta = self.v(:,X2) - self.v(:,X1);
    end


    % ------------------------------------------
    function [mod, freq] = wavelet(self,x)
    % ------------------------------------------
    % Compute the wavelet transform of some data
    % ------------------------------------------


        Fs = self.wavelet_params.Fs                              ;
        fo = self.wavelet_params.fo                              ;
        fl = self.wavelet_params.fl                              ;
        fu = self.wavelet_params.fu                              ;
        fn = self.wavelet_params.fn                              ;
    

        % Single input
        if size(x,2) == 1
        [~,freq,mod] = freq_inst_morlet(x,Fs,fl,fu,fn,fo)        ;
        
        % Mutliple Input
        else
            for k = 1:size(x,2)
                [~,freq,out] = freq_inst_morlet(x(:,k),Fs,fl,fu,fn,fo)        ;
                mod{k} = out;
            end
        end
    end

    % ------------------------------------------
    function [E, EPE, EKE, TE] = linenergy(self,Nt)
    % ------------------------------------------
    % Compute the kinetic, potential, and total energy
    % with respect to time (linear energy only)
    % ------------------------------------------
    if nargin < 2;
        Nt = 100;
    end
    if isempty(self.sys.K)
        [K,M,~] = self.sys.system_mats();
    else
        K = self.sys.K;
        M = self.sys.M;
    end
    if nargin == 1, Nt = 100; end
    % ------------------------
    % Compute energy over time interval
    tsim = self.t;
    R = 1;
    EKE = zeros(Nt,1);EPE = zeros(Nt,1); TE =zeros(Nt,1); E = zeros(Nt,1);
    for ii = round(linspace(1,numel(tsim), Nt))
        EPE(R) =  self.v(ii,:)*M*self.v(ii,:)'/2;
        EKE(R) =  self.x(ii,:)*K*self.x(ii,:)'/2;
        E(R) = EPE(R)+EKE(R);
        TE(R) = tsim(ii);
        R = R + 1;
    end
    % ------------------------
    end

    % ------------------------------------------
    function F = FNL(self,Nt)
    % ------------------------------------------
    % Compute the nonlinear forces experienced d
    % by the VI
    % ------------------------------------------
    if nargin < 2
        Nt = 1000;
    end
    if isempty(self.sys.K)
        [K,M,~] = self.sys.system_mats();
    else
        K = self.sys.K;
        M = self.sys.M;
    end
    kc = self.sys.VIforce.kc;
    d = self.sys.VIforce.d;
    VIsites = self.sys.VIforce.VIsites;
    j =1 ;

    if ~isempty(self.x)
        if numel(kc) == 1
        for SITE = VIsites
            delta = self.x(:,SITE+1)-self.x(:,SITE);
            F(:,j) = -kc*( max(-delta-d(SITE),0) .^(3/2) - max(delta-d(SITE),0) .^(3/2) ) ;
            j =j + 1;
        end
        else
        for SITE = VIsites
            delta = self.x(:,SITE+1)-self.x(:,SITE);
            F(:,j) = -kc(j)*( max(-delta-d(SITE),0) .^(3/2) - max(delta-d(SITE),0) .^(3/2) ) ;
            j =j + 1;
        end
        end
    else
        if numel(kc) == 1
        for SITE = VIsites
            delta = self.xc(:,SITE+1)-self.xc(:,SITE);
            F(:,j) = -kc*( max(-delta-d(SITE),0) .^(3/2) - max(delta-d(SITE),0) .^(3/2) ) ;
            j =j + 1;
        end
        else
        for SITE = VIsites
            delta = self.xc(:,SITE+1)-self.xc(:,SITE);
            F(:,j) = -kc(j)*( max(-delta-d(SITE),0) .^(3/2) - max(delta-d(SITE),0) .^(3/2) ) ;
            j =j + 1;
        end
             end
    end
    % ------------------------
    end


    % ------------------------------------------
    function   [PEntropy,PEinst] =    entropy(self,Y,samps)
    % ------------------------------------------
    % Compute the energy in seperate bands
    % ------------------------------------------
    % Y - time series matrix to compute spectral entropy
    % samps - indicies in time to sample entropy
    
    if nargin < 2
        Y = self.v';
    end
    if nargin < 3
        if ~isempty(self.t)
        samps = round(linspace(1,numel(self.t),100 ));
        else
        samps = round(linspace(1,numel(self.tc),100 ));
        end
   end
    % Get the size of the entropy matrix
    x = Y(1,:)';
    tstsize = numel(pentropy(x,1:length(x)));
    PEntropy = zeros(100,tstsize);
    KK = 1;

    % ----------------
    % Loop through time snapshots
    for k = samps
        x = Y(k,:)';
        PEntropy(KK,:) =  pentropy(x,1:length(x));
        PEinst(KK)     =  pentropy(x,1:length(x),'Instantaneous',false);
        KK = KK + 1;
    end
    end


    % ------------------------------------------
    function   Wband =    Wbands_vel(self,tsamp)
    % ------------------------------------------
                y = self.v(:,:);
                Wband = self.Wbands(y,tsamp);
    end
    % ------------------------------------------
    function   Wband =    Wbands_pos(self,tsamp)
    % ------------------------------------------
                y = self.x(:,:);
                Wband = self.Wbands(y,tsamp);
    end

    % ------------------------------------------
    function   [Wband Xuse ]=    Wbands(self,y,tsamp)
    % ------------------------------------------
    % Compute the energy in seperate bands
    % ------------------------------------------
    
    if nargin<3
        if ~isempty(self.t)
        tsamp = round(linspace(1,numel(self.t),100));
        else
        tsamp = round(linspace(1,numel(self.tc),100));
        end
    end
    if nargin<2
        if ~isempty(self.v)
            y = self.v;
        else
            y= self.vc;
        end
    end


    % --------------------------- ;
    %   Wavelet params            ;
    % --------------------------- ;
    if isempty(self.decomp_params.wavelet_params.Fs)
        self.decomp_params.wavelet_params.Fs = self.wavelet_params.Fs;
    end
    if isempty(self.decomp_params.wavelet_params.fo)
        self.decomp_params.wavelet_params.fo = self.wavelet_params.fo;
    end
    if isempty(self.decomp_params.wavelet_params.fl)
        self.decomp_params.wavelet_params.fl = self.wavelet_params.fl;
    end
    if isempty(self.decomp_params.wavelet_params.fu)
        self.decomp_params.wavelet_params.fu = self.wavelet_params.fu;
    end
    if isempty(self.decomp_params.wavelet_params.fn)
        self.decomp_params.wavelet_params.fn = self.wavelet_params.fn;
    end
    if isempty(self.decomp_params.nbands)
        self.decomp_params.nbands = 12;
    end
    Fo = self.decomp_params.wavelet_params.fo;
    fi = self.decomp_params.wavelet_params.fl;
    fu = self.decomp_params.wavelet_params.fu;
    nF = self.decomp_params.wavelet_params.fn;
    Fs = self.decomp_params.wavelet_params.Fs;
 

    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%
    y = y(tsamp,:);
    Xuse = y;
    % -------------------------------
    % Loop Through Time Snapshots
    for tt = 1:numel(tsamp);%[round(linspace(1,numel(tsim),100))]
        x  = y(tt,:)'              ;
        [~,frequency_WT,MODS] = freq_inst_morlet(x,Fs,fi,fu,nF,Fo);
        % -------------------------------
        % Define patches to take inverse transforms of
        NK =  self.decomp_params.nbands;
        N_modes = NK;
        LX = length(x);
        loc1 = length(x);
        clear POSS
        for kk = 1:NK
            1+(kk-1)*nF/NK;
            POSS{kk}    = [[1;1;loc1;loc1] [1+(kk-1)*nF/NK;1+kk*nF/NK;1+kk*nF/NK;1+(kk-1)*nF/NK]];
        end
        nbands = self.decomp_params.nbands;
        bins = round(linspace(1,nF,nbands));
        idx =1 ;
        lbin = min(diff(bins));
        for r = 1:nbands-1
            Wband(:,:,r,tt) = MODS(:,idx:(idx+lbin));
            idx = idx + lbin;
        end
    end
    % ------------------------------------------
    end


    % ------------------------------------------
    function   TSband =    Wbands_inv(self,Wbands)
    % ------------------------------------------
    % Compute Time series of different spectral partitions
    % ------------------------------------------
    

    % Get the dimensions of the wavelet partition structure
    nsamp = size(Wbands,4);
    npart = size(Wbands,3);
    
    % --------------------------- ;
    %   Wavelet params            ;
    % --------------------------- ;
    y = self.x(:,:);
    Wband = self.Wbands(y,tsamp);
    


            for tt = [round(linspace(1,numel(tsim),100))]
            
            clf
            ndofN = 1;
            XX = ndof+[1:2:ndof]               ;
            x  = y(tt,XX)'              ;
            NDOF = size(x,1)            ;
            l_new = NDOF                ;
            XN = zeros(l_new,1)         ;
            XN(1:size(x,1)) = x         ;
            nfourier = l_new            ;
            npt = nfourier/2            ;
            freq = Fs*([0:nfourier-1])/nfourier       ;                                 
            FREQ = Fs*([0:npt-1 npt:-1:1])/nfourier   ;
            interval_freq = FREQ        ;
            a = Fo./interval_freq       ;
            fft_MW = conj(bsxfun(@times,pi^(1/4)*(2^0.5)*(exp(-0.5*(2*pi*(bsxfun(@times,FREQ',a)-Fo)).^2)-1*exp(-0.5*(2*pi^2*(bsxfun(@times,FREQ',a).^2+Fo.^2)))),sqrt(a)));
            [~,frequency_WT,MODS] = freq_inst_morlet(XN,Fs,fi,fu,nF,Fo);
            
            
            if pltmodes
            subplot(1,2,1)
            hold on
            colormap([ones(3,3); jet.^2]);
            view([0 90])
            light('Position',[-1 0 0],'Style','local')
            ax = gca;
            xlim([1 length(XN)])
            end
            
            NK = 12;
            N_modes = NK;
            LX = length(XN);
            loc1 = length(XN);round(VIsites(1)/2);
            clear POSS
            for kk = 1:NK
                1+(kk-1)*nF/NK;
            POSS{kk}    = [[1;1;loc1;loc1] [1+(kk-1)*nF/NK;1+kk*nF/NK;1+kk*nF/NK;1+(kk-1)*nF/NK]];
            %POSS{kk+NK} = [[loc1; loc1; LX; LX] [1+(kk-1)*nF/NK;1+kk*nF/NK;1+kk*nF/NK;1+(kk-1)*nF/NK]];
            end
            
            
            CMAP = jet(N_modes);
            modes = zeros(length(XN),N_modes,ndofN);
            
            
            for counter = 1:N_modes
                poss = POSS{counter};
                freq_repmat = repmat(freq,length(poss(:,2)),1);
                patch_freq  = frequency_WT(1,ceil(poss(:,2)));
                indices     = sum(freq_repmat<=patch_freq',2);
                Mask1       = poly2mask(POSS{5}(:,1),indices,length(freq)/2,length(XN));
                Mask2       = poly2mask(poss(:,1),indices,length(freq)/2,length(XN));
                
                if counter <5
               % Mask2(Mask2==Mask1) = 0;
                end
                pgon = polyshape(poss);
                p2=plot(pgon);
                p2.FaceAlpha = 0.25;
                p2.FaceColor = CMAP(counter,:);
                p2.LineWidth = 1;
            
                text(poss(1,1)+5,poss(4,2)+5,0, ['$\mathcal{K}_{',num2str(counter),'}$'],'color',CMAP(counter,:)/2,'fontsize',20 )
                for counter_NDOF = 1:ndofN
                    tff = fft(XN(:,counter_NDOF),nfourier);
                    noyau2 = bsxfun(@times,fft_MW,tff);
                    resut2 = ifft(noyau2,nfourier);
                    Patch = [Mask2' fliplr(Mask2')].*resut2;
                    fft_Patch = fft(Patch,nfourier);
                    fft_mode_2D = fft_Patch./fft_MW;
                    fft_mode = diag(fft_mode_2D);
                    fft_mode(isnan(fft_mode)) = 0;
                    modes(:,counter,counter_NDOF) = ifft(fft_mode,nfourier,'symmetric');
                end
            end
         end

    end



end
end


