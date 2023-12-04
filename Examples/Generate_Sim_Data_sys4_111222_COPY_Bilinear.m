clear all
clc


savedata= input('Do you want to save the data? (0) - no, (1) - yes');
% ------------------------------------------------- ;
% Define paramters to loop across                   ;
% ------------------------------------------------- ;
datapath = '..\Data\F22\10-05-22_data\4band_loop_EVI_kVI_COPY\';
if exist(datapath)
    disp('Filepath exists! Press a key to continue')
    pause;
    disp('Are you sure you want to overwrite:')
    disp(datapath)
    pause;
end
addpath('..\Methods\')
addpath('..\Test_Sim_Codes\')
addpath('..\..\format files\')
plothelper.plotformat
clc
%%
Ndof = 1200;
for Ncyc = [15 30]
RVI  = 5e-3;
for knominal = 1e4
for nvii  = [5 1];
datapath =[ '..\Data\F22\11-11-22_data\Fmin1e1_4band_loop_EVI_kVI_Bilinear_knominal_',num2str(log10(knominal)),'_RVI',num2str(RVI),'_Ncyc',num2str(Ncyc),'_Ndof',num2str(Ndof)];
if exist(datapath)
%     disp('Filepath exists! Press a key to continue')
%     pause;
%     disp('Are you sure you want to overwrite:')
%     disp(datapath)
%     pause;
end

eta = 1;
for EVI = flip([1e9])
eta2 = 1;
fampidx = 1;
for Famp = logspace(1,4,50)*knominal/1e4



tic
% ------------------------------------------------- ; 
% Lattice setup                                     ; 
% ------------------------------------------------- ; 
DoF         = 1                                     ; 
Ncell       = Ndof                                  ;     
%if eta >= 1    
m1          = 5e-2                                  ; 
m2          = m1                                    ; 
m3          = m2*eta                                ;
m4          = m1                                    ;
k1          = 4e5                                   ;  
k2          = 4e5                                   ; 
k3          = k2*sqrt(eta)/5                        ; 
k4          = k1                                    ;
Kg          = 2e3                                   ; 
config      = 4                                     ; 
grounded    = 0                                     ;       
ndof        = DoF*Ncell                             ; 
KVI         = round(k2)/1                           ; 
NVI         = nvii                                  ;     
if config <3
VIsites     = round(linspace(Ncell/4,2*Ncell/4,NVI)); 
VIsites     = 2*(rem(VIsites,2)+floor(VIsites/2)) -1; 
elseif config == 3
VIsites     = round(linspace(Ncell/4,2*Ncell/4,NVI)) ; 
VIsites     = 3*(rem(VIsites,3)+floor(VIsites/3))- 1 ; 
NVI = numel(VIsites)                                 ;
elseif config == 4
VIsites     = round(linspace(Ncell/4,2*Ncell/3,NVI)) ; 
if NVI == 1; VIsites = round(Ncell/3);end
VIsites     = 4*(rem(VIsites,4)+floor(VIsites/4))-2  ;

% PARAMETER CASCADE
% -------------------------------------
% %VIsites     = [ 3*(rem(VIsites,3)+floor(VIsites/3))- 2] ;
% bsites = 2:4:ndof;
% %VIsites     = sort([VIsites VIsites+4 VIsites+8 VIsites+12 VIsites+16 VIsites+20 VIsites+24 VIsites+28]);
% VIsites = bsites(40:4:end-40);
% -------------------------------------
NVI = numel(VIsites);
end
VIsites     = sort(VIsites)                         ; 

Cnorm       = 0                                     ; 
m           = [m1 m2 m3 m4]                         ;
k           = [k1 k2 k3 k4] *1e1                    ;

%%% OPIMTAL HEAVY %%%
% lambda = .1;eta = 1.95; %  'optimal'
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1/eta eta 1]*5e-2);

% lambda = 1;eta = 1.95; %  'non optimal'
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1/eta eta 1]*5e-2);
%lambda = .3;eta = .4; % 'non-optimal'

%%% OPIMTAL %%%
lambda = .1;eta = .5;
k = [2 1 lambda 2] * knominal;
m = ([1 1-eta eta 1]*5e-2);
KVI = k(3);



%%% OTEHR OPIMTAL %%%
% lambda = .1;eta = .3;
% k = [1 1 lambda 1] * knominal;
% m = ([1 1-eta eta 1]*5e-2);

%%% OPIMTAL (m2=1) %%%
% lambda = .025;eta = .1050;
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1 eta 1]*5e-2);

% Get the 'optimal set' of lambda and eta per the grid search
% lambda = .01; eta = .25;
% k = [1 1 lambda 1] * 1e4;
% m= ([1 4*eta eta 1]*2e-2);

% ------------------------------------------------- ;
% Define Lattice Object                             ;
% ------------------------------------------------- ;
lattice = VI_LatticeSys(ndof,config,m,k,VIsites,KVI,Kg,grounded);

if fampidx > 0 % Add some damping to the system to stop the constant ringing 
lattice.c = [5 5]*1e-8;
end
% figure(43);clf
% plot(lattice.dispersion)
%yyaxis right;
hold on
[K,M] =lattice.system_mats;
%plot(sqrt(eig(K,M)),'.');pause(.01)
%%
% ------------------------------------------------- ;
% Define location to save data                      ;
% ------------------------------------------------- ;
folder   = ['config_',num2str(config),'_eta_',num2str(eta),'_EVI_',num2str((EVI)),'_NVI_',num2str(nvii)];
fname    = ['latticedata_config',num2str(config),'_famp',num2str(fampidx)];
dataloc  = strcat(datapath,folder);
dataloc_comp  = strcat(datapath,'/compressed/',folder);
dataloc_full  = strcat(datapath,'/full/',folder);
dataloc_minimal  = strcat(datapath,'/minimal/',folder);
if ~isfile([dataloc_minimal,'/',fname,'.mat'])
if ~isdir(dataloc_full)
    mkdir(dataloc_full)
end
if ~isdir(dataloc_comp)
    mkdir(dataloc_comp)
end
if ~isdir(dataloc_minimal)
    mkdir(dataloc_minimal)
end
% ------------------------------------------------- ;
% Compute dispersion relation                       ;
% ------------------------------------------------- ;
dkappa                  = .025                      ;
[DISP,wavevec]          = lattice.dispersion(dkappa);
dispersion.wavevec      = wavevec                   ;  
dispersion.DISP         = DISP                      ;
lattice.Dispersion      = dispersion                ;

% ------------------------------------------ ;
% Profile of excitation                      ;
% ------------------------------------------ ;
band    = 2         ;
BandLoc =  lattice.K_max_group_vel(2)/pi     ; % max Vg
kappa   = wavevec   ;

% ------------------------------------------ ;
%       Select excitation freqs              ;
% ------------------------------------------ ;
bandidx = find(kappa >= pi*BandLoc)                             ;
if band ==  1,   omega =DISP(bandidx(1),1)                      ;     end
if band == 2,    omega = DISP(bandidx(1),2)                     ;     end
if band == 1.5,  omega = ( max(DISP(:,1) ) + min(DISP(:,2)))/2  ;     end
omega = real(omega)                                             ;  

% ------------------------------------------ ;
%       Gauss tone burst options             ;
% ------------------------------------------ ;
periods  = Ncyc                                ;
amp      = Famp                              ;
FLOC     = 1                                 ;
% ------------------------------------------ ;
Fin = excitation(omega,periods,amp,FLOC)     ;
% ------------------------------------------ ;
Fin.band        = band;
Fin.BandLoc     = BandLoc;


% ------------------------------------------     ;
%       Nonlinear Forcing                        ;
% ------------------------------------------     ;
timeVI      = 1e-6                               ;
radius      = RVI                                ;
TimeVI      = 1e-5                               ;
rescoef     = .7                                 ;
E           = EVI                                ;
nu          = 0.3                                ;
dnominal    = 100*ones(ndof,1)                   ;
sided       = 2                                  ;
dVI         =  logspace(-2.65,-2.75,NVI)         ;
% dVI         =  logspace(-3,-3,NVI)         ;

% --------------------------------------------------
% PARAMETER CASCDE
% --------------------------------------------------
% E(1:4:NVI) = 1e6; E(2:4:NVI) = 1e7             ; 
% E(3:4:NVI) = 1e8; E(4:4:NVI) = 1e9             ;
% E           = logspace(6,8,NVI)                           ;
% dVI         =  logspace(-2.65,-2.65,NVI)         ;
% dVI         =  linspace(10^-2,10^-3,NVI)       ;
% --------------------------------------------------
FVI         = VIforce(TimeVI, E, nu, radius ,rescoef , dnominal,sided, VIsites, dVI);
% --------------------------------------------------

% ------------------------------------------ ;
%      Update lattice object                 ;
% ------------------------------------------ ;
lattice.VIforce = FVI                        ;
lattice.excitation = Fin                     ;

% ------------------------------------------------------------------        ;
%      Simulate System                                                      ;
% ------------------------------------------------------------------        ;
S = lattice.system_mats; S(S~=0) = 1; S = sparse(S)                         ;
tstop =  Fin.tforce+lattice.ndof/lattice.group_vel(omega)/config            ;
Simparams = SimParams()                                                     ;
Simparams.tsim  = linspace(0,tstop, 25000)                                  ;                                                        
Simparams.opts = odeset('RelTol',1e-10, 'AbsTol',Fin.amp*1e-9,'JPattern',S) ;
Simparams.IC = zeros(1,lattice.ndof*2)                                      ;
Simparams.Integrator = 'ode78'                                              ;
lattice.SimParams = Simparams                                               ;
% [t,y] = lattice.simulate()                                               ;
[t,y] = lattice.simulate_bilinear()                                               ;


% ------------------------------------------ ;
%      Define Data Class                     ;
% ------------------------------------------ ;
latticedata     = VI_data_processing(t,y(:,1:lattice.ndof), y(:,(1+lattice.ndof):2*lattice.ndof) );
latticedata.sys = lattice;

% ------------------------------------------ ;
%      Define the wavelet parameters         ;
% ------------------------------------------ ;
latticedata.wavelet_params.fl = 0;
latticedata.wavelet_params.fu = pi;
latticedata.wavelet_params.fn = 300;
latticedata.wavelet_params.fo = 2;
latticedata.wavelet_params.Fs = pi;

% ------------------------------------------ ;
%      Compute some basic results            ;
% ------------------------------------------ ;
[EPE, EKE, E, TE]   = latticedata.linenergy();
Y                   = latticedata.v;
samps               = round(linspace(1,numel(latticedata.t),4 ) );
tstamps             =  round(linspace(1,numel(latticedata.t),4 ) );
xinds               = 1:2:ndof;
[wavelet,frq]       = latticedata.wavelet(latticedata.v(tstamps,xinds)');
[PE1,PE2]           = latticedata.entropy(Y,samps);
maxF                = 1.5*max(DISP(:));
[mat,omega,kappa]   = latticedata.Fdomain(latticedata.v,maxF);          % Numerical dispersion
latticedata.entropydata.PE1 = PE1;
latticedata.entropydata.PE2 = PE2;
latticedata.waveletdata.timestamps = tstamps;
latticedata.waveletdata.wdata      = wavelet;
latticedata.waveletdata.frq        = frq;
latticedata.waveletdata.xinds      = xinds;



% ------------------------------------------------------
%  figure(1)
% clf
%subplot(1,2,1)
%imagesc(latticedata.v.^2);%colormap(linspecer.^2)
% set(gca,'colorscale','log')
% title(['(2) FampIDX = ',num2str(fampidx),', Amp = ',num2str(round(Famp))])
% plothelper.tickformat
% yyaxis right
% plot(linspace(1,600,100),latticedata.energy.E)
% figure(1)
% ------------------------------------------------------
disp(['findex ',num2str(fampidx),' Amp = ',num2str(round(Famp))])

% ------------------------------------------------------
if savedata,%savedata = saveobj(latticedata);
save([dataloc_full,'/',fname],'latticedata'),end
% ------------------------------------------------------
latticedata.waveletdata = [];
latticedata.Process;
latticedata.compress;
latticedata.clearTS;
% subplot(1,2,2)
figure(1);clf
imagesc(latticedata.vc.^2)%colormap(linspecer.^2)
set(gca,'colorscale','log')
title(['(2) FampIDX = ',num2str(fampidx),', Amp = ',num2str(round(Famp))])
yyaxis right
plot(linspace(1,600,100),latticedata.energy.E)
pause(0.01)
plothelper.tickformat

% set(gca,'colorscale','log')
% title(['(2) FampIDX = ',num2str(fampidx),', Amp = ',num2str(round(Famp))])
% pause(.01)
% plothelper.tickformat
% ------------------------------------------------------
if savedata,%savedata = saveobj(latticedata);
save([dataloc_comp,'/',fname],'latticedata'),end
% ------------------------------------------------------


latticedata.xc =[latticedata.xc(:,VIsites) latticedata.xc(:,VIsites+1)];
latticedata.vc =[latticedata.vc(:,VIsites) latticedata.vc(:,VIsites+1)];
latticedata.FFTdata         = [];
latticedata.dispersion      = [];
latticeadata.wbands         = [];
latticedata.waveletdata     = [];
latticedata.entropydata     = [];
% ------------------------------------------------------
if savedata,%savedata = saveobj(latticedata);
save([dataloc_minimal,'/',fname],'latticedata'),end
% ------------------------------------------------------
end
toc
fampidx = fampidx +1;
end
end
end
end
end
%%
figure;imagesc(kappa,omega,abs(mat));set(gca,'ydir','normal');set(gca,'colorscale','lin')
hold on
plot(wavevec,DISP,'--','color',[ 1 1 1 .5])

figure;imagesc(abs(y(:,ndof:ndof*2)));set(gca,'colorscale','log')
figure;plot(E)



%%

ld = copy(latticedata)
ld.compress;
ld.clearTS


%%
% sampleload = load('C:\Users\josht\Box\Codes\NL MetaMat\Data\F22\09-18-22_data\PrelimStudy_VaryEta_3band\compressed\config_3_eta_0.2\latticedata_famp25.mat')
% figure;imagesc(sampleload.latticedatacomp.vc.^2);set(gca,'colorscale','log')
% 
