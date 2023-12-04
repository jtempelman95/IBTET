clear all
clc
clear all
close all
%%
figpath = '..\..\..\Josh\Weekly Meetings\10-06-22\figures\masterfig\';
if~isdir(figpath);mkdir(figpath);end
addpath('..\Methods\')
addpath('..\Test_Sim_Codes\')
addpath('..\..\format files\')
%%
% ------------------------------------------------- ;
% Define paramters to loop across                   ;
% ------------------------------------------------- ;
for eta  = 1 ;    [1/10 1/5 1/2 1 2 5 10];
    eta2 = .5;
fampidx  = 1    ;
for Famp = 60;  logspace(.1,1.25,25)    
tic
% ------------------------------------------------- ; 
% Lattice setup                                     ; 
% ------------------------------------------------- ; 
DoF         = 1                                     ; 
Ncell       =  1200                                  ;     
%if eta >= 1
m1          = 5e-2                                  ; 
m2          = m1*eta                                ; 
m3          = m2*eta2                               ;
m4          = m1                                    ;
k1          = 4e3                                   ;
k2          = 4e3*sqrt(eta)                         ; 
k3          = k2*sqrt(eta2)/5                       ; 
k4          = k1                                    ;
Kg          = 2e3                                   ; 
config      = 4                                     ; 
grounded    = 0                                     ;       
ndof        = DoF*Ncell                             ; 
KVI         = round(k2)/1                           ; 
NVI         = 5                                     ;     
if config <3
VIsites     = round(linspace(Ncell/4,2*Ncell/4,NVI)); 
VIsites     = 2*(rem(VIsites,2)+floor(VIsites/2)) -1; 
elseif config == 3
VIsites     = round(linspace(Ncell/4,2*Ncell/4,NVI)) ; 
VIsites     = 3*(rem(VIsites,3)+floor(VIsites/3))- 1 ; 
%VIsites     = [ 3*(rem(VIsites,3)+floor(VIsites/3))- 2] ; 
NVI = numel(VIsites)                                 ;
elseif config == 4
VIsites     = round(linspace(Ncell/4,2*Ncell/3,NVI)) ; 
if NVI == 1; VIsites = Ncell/4;end
if NVI == 1;  bsites = 2:4:ndof; VIsites = bsites(round(end/2));end
VIsites     = 4*(rem(VIsites,4)+floor(VIsites/4))-2  ;

% PARAMETER CASCADE
% -------------------------------------
% VIsites     = [ 3*(rem(VIsites,3)+floor(VIsites/3))- 2] ;
% bsites = 2:4:ndof;
% VIsites     = sort([VIsites VIsites+4 VIsites+8 VIsites+12 VIsites+16 VIsites+20 VIsites+24 VIsites+28]);
% VIsites = bsites(40:end-40);
% -------------------------------------
NVI = numel(VIsites);
end
VIsites     = sort(VIsites)                         ; 
Ccoeff      = [.01 .01]                             ; 
Cnorm       = 0                                     ;  
m           = [m1 m2 m3 m4]                         ;
k           = [k1 k2 k3 k4] *1e1                    ;

%if kk == 1
%%% OPIMTAL %%%
lambda = .1;eta = .5;
k = [2 1 lambda 2] * 1e4;
m = ([1 1-eta eta 1]*5e-2);
KVI = k(3);

% m = ([1 1-eta eta/4 1]*5e-2);
% k = [2 1 lambda/15 2] * 1e4;

% lambda = .5;eta = .2450;
% k = [2 1 lambda 2]*1e4;
% m = ([1 1 eta 1]*5e-2);


% lambda = .025;eta = .1050;
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1 eta 1]*5e-2);

%  lambda = .01;eta = 0.1900;
% k = [1 1 lambda 1] * 1e4;
% m= ([1 5*eta eta 1]*5e-2);

%figpath = '..\..\..\Josh\Weekly Meetings\09-22-22\figures\sys_opt';
%elseif kk==2
%%% OPIMTAL HEAVY %%%
% lambda = .1;eta = 1.95;
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1/eta eta 1]*5e-2);
% figpath = '..\..\..\Josh\Weekly Meetings\09-22-22\figures\sys_opt_heavy';
%elseif kk == 3
% lambda = .4;eta = 1.5; %  'non optimal'
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1/eta eta 1]*5e-2);
% figpath = '..\..\..\Josh\Weekly Meetings\09-22-22\\figures\sys_nom_heavy';
%elseif kk ==4
% lambda = .15;eta = .35; %  'non optimal'
% k = [2 1 lambda 2] * 1e4;
% m = ([1 1/eta eta 1]*5e-2);
% figpath = '..\..\..\Josh\Weekly Meetings\09-22-22\\figures\sys_nom';
%end
%if~isdir(figpath);mkdir(figpath);end

% ------------------------------------------------- ;
% Define Lattice Object                             ;
% ------------------------------------------------- ;
lattice = VI_LatticeSys(ndof,config,m,k,VIsites,KVI,Kg,grounded);

figure;
plot(lattice.dispersion)
%yyaxis right;
hold on
[K,M] =lattice.system_mats;
plot(sqrt(eig(K,M)),'.');pause(.01)
%%
% ------------------------------------------------- ;
% Define location to save data                      ;
% ------------------------------------------------- ;
% datapath = '..\Data\F22\';
% folder   = ['config_',num2str(config),'_eta_',num2str(eta)];
% fname    = ['latticedata_famp',num2str(fampidx)];
% dataloc  = strcat(datapath,folder);
% if ~isdir(dataloc)
%     mkdir(dataloc)
% end

% ------------------------------------------------- ;
% Compute dispersion relation                       ;
% ------------------------------------------------- ;
dkappa                  = .025                      ;
[DISP,wavevec]          = lattice.dispersion(dkappa);
dispersion.wavevec      = wavevec                   ;  
dispersion.DISP         = DISP                      ;
lattice.Dispersion      = dispersion                ;
%%
% ------------------------------------------ ;
% Profile of excitation                      ;
% ------------------------------------------ ;
band    = 2                                  ;
BandLoc =  lattice.K_max_group_vel(2)/pi     ; % max Vg
kappa   = wavevec                            ;
% 
% band = 1;
% BandLoc = 1/10;
%%

% ------------------------------------------ ;
%       Select excitation freqs              ;
% ------------------------------------------ ;
bandidx = find(kappa >= pi*BandLoc)                             ;
if band ==  1,   omega =DISP(bandidx(1),1)                      ;     end
if band == 2,    omega = DISP(bandidx(1),2)                     ;     end
if band == 1.5,  omega = ( max(DISP(:,1) ) + min(DISP(:,2)))/2  ;     end
if band >2,    omega = DISP(bandidx(1),band)                    ;     end
omega = real(omega)                                             ;  

% ------------------------------------------ ;
%       Gauss tone burst options             ;
% ------------------------------------------ ;
periods  = 15                                 ;
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
radius      = .001                               ;
TimeVI      = 1e-10                              ;
rescoef     = .7                                 ;
E           = 1e16                               ;
nu          = 0.3                                ;
dnominal    = 100*ones(ndof,1)                   ;
sided       = 2                                  ;
dVI         =  logspace(-2.65,-2.75,NVI)         ;

% --------------------------------------------------
% PARAMETER CASCDE
% --------------------------------------------------
% E(1:4:NVI) = 1e6; E(2:4:NVI) = 1e7               ; 
% E(3:4:NVI) = 1e8; E(4:4:NVI) = 1e9               ;
% E = logspace(6,8,NVI)                            ;
% dVI         =  logspace(-2.65,-2.65,NVI)         ;
% dVI         =  linspace(10^-2.5,10^-3,NVI)       ;
% --------------------------------------------------
FVI         = VIforce(TimeVI, E, nu, radius ,rescoef , dnominal,sided, VIsites, dVI);
% --------------------------------------------------

% ------------------------------------------ ;
%      Update lattice object                 ;
% ------------------------------------------ ;
lattice.VIforce = FVI                        ;
lattice.excitation = Fin                     ;

% ------------------------------------------------------------------ ;
%      Simulate System                                               ;
% ------------------------------------------------------------------ ;
S = lattice.system_mats; S(S~=0) = 1; S = sparse(S)                         ;
tstop =  Fin.tforce+lattice.ndof/lattice.group_vel(omega)/config            ;
Simparams = SimParams()                                                     ;
Simparams.tsim  = linspace(0,tstop, 30000)                                  ;                                                        
Simparams.opts = odeset('RelTol',1e-8, 'AbsTol',Fin.amp*1e-7,'JPattern',S);
Simparams.IC = zeros(1,lattice.ndof*2)                                      ;
Simparams.Integrator = 'ode45'                                              ;
lattice.SimParams = Simparams                                               ;
% [t,y] = lattice.simulate_bilinear()                                                  ;
[t,y] = lattice.simulate()                                                  ;

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
samps               = round(linspace(1,numel(latticedata.t),100 ) );
tstamps             = round(linspace(1,numel(latticedata.t),10 ) )
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
savedata = saveobj(latticedata);
% save([dataloc,'/',fname],'latticedata')


% intermediate plotting
if (mod(fampidx,5)-1)==0
figure(1)
clf
imagesc(latticedata.v.^2);colormap(linspecer.^2)
title(['FampIDX=',num2str(fampidx)])
pause(.01)
end
disp(['findex ',num2str(fampidx)])
toc
fampidx = fampidx +1;
end
end


%%
Ebands = latticedata.BandEnergy;
latticedata.Eband = Ebands;
latticedata.compress

%% Paper Figure
addpath('..\')
subplot = @(m,n,p) subtightplot (m, n, p, [0.14 0.065]*1, [0.1 0.1]*1, [0.05 0.05]*1);

figfontsize = 11;
clc
latticedata.decomp_params.nbands = 12;
asites = 1:4:ndof;
y = latticedata.v(:,asites);
latticedata.decomp_params.wavelet_params.Fs = pi*8;
latticedata.decomp_params.wavelet_params.fu = pi*4;
[wb xb] = latticedata.Wbands(y);
sz = size(wb);
CM = [1-pink(170).^3];


figure(44)
clf

% -------------------------------- ;
% Subplot: dispersion input        ;
% -------------------------------- ;
subplot(2,4,5)
set(gca,'fontsize', figfontsize);pause(.01);
 set(gca,'fontsize', figfontsize)

L1 = size(latticedata.v,1);
L2 = size(latticedata.v,2);
S10 = round(1:L1/3.8);
S20 = round(1:L2/3.8);

S10 = round(1e3:L1/3.5);
S20 = round(8:L2/3-8);
[mat,omega,kappa]   = latticedata.Fdomain(latticedata.v(S10,S20),maxF);          % Numerical dispersion
plothelper.plotformat
kappa = kappa*2;
R = 5;
pltd = abs(mat)/max(abs(mat(:)));
imagesc(linspace(kappa(1),kappa(end),numel(kappa)*R),...
linspace(omega(1),omega(end),numel(omega)*R),imresize(pltd,R));set(gca,'ydir','normal');
 set(gca,'fontsize', figfontsize)
set(gca,'colorscale','log')
ax = gca; pos = ax.Position;
cb=colorbar;
cb.Position = [pos(1)+pos(3)/2 pos(2)+pos(4)*.85 .05 .01];
cb.Orientation = 'horizontal';
cb.AxisLocation = 'in';
cb.Color = 'k';
caxis([.005 1])
hold on
plot(wavevec,DISP,'--','color',[ [1 1 1]/1e3 .5])
plothelper.dispersionplot
title('Input Segments')
plothelper.tickformat
ax = gca;
text(ax.XLim(1),ax.YLim(2)*1.05,'(e)','fontsize',12)
    set(gca,'fontsize', figfontsize)


    
% -------------------------------- ;
% Subplot: dispersion output       ;
% -------------------------------- ;
subplot(2,4,6)
set(gca,'fontsize', figfontsize)

L1 = size(latticedata.v,1);
L2 = size(latticedata.v,2);
S1 = round(L1*4/7:L1);
S2 = round(1*L2/4:L2);
S1 = round(L1*4/7:L1-1e3);
S2 = round(12:L2-12);
[mat,omega,kappa]   = latticedata.Fdomain(latticedata.v(S1,S2),maxF);          % Numerical dispersion
plothelper.plotformat
kappa = kappa*2;
R = 5;
pltd = abs(mat)/max(abs(mat(:)));
imagesc(linspace(kappa(1),kappa(end),numel(kappa)*R),...
linspace(omega(1),omega(end),numel(omega)*R),imresize(pltd,R));set(gca,'ydir','normal');
set(gca,'colorscale','log')
ax = gca; pos = ax.Position;
cb=colorbar;
cb.Position = [pos(1)+pos(3)/2 pos(2)+pos(4)*.85 .05 .01];
cb.Orientation = 'horizontal';
cb.AxisLocation = 'in';
cb.Color = 'k';
caxis([.005 1])
hold on
set(gca,'fontsize', figfontsize)
plot(wavevec,DISP,'--','color',[ [1 1 1]/1e3 .5])
plothelper.dispersionplot
xlim([-2 2]*pi/2)
title('Output Segment')
ax = gca;
plothelper.tickformat

set(gca,'fontsize', figfontsize)
set(gca,'fontsize', figfontsize)
ax = gca;
text(ax.XLim(1),ax.YLim(2)*1.05,'(f)','fontsize',12)

t = latticedata.t;
tc = latticedata.tc;

% ax = gca;
% text(ax.XLim(1),ax.YLim(2)*1.22,'(a)','fontsize',12)


% -------------------------------- ;
% Subplot: Energy                  ;
% -------------------------------- ;
subplot(2,4,1)
set(gca,'fontsize', figfontsize)

pltd = movmean(movmean(abs(latticedata.v).^1,60),10,2)';
pltd = pltd/max(pltd(:));
axx =max(pltd(:));
imagesc(tc,1:ndof,pltd)
caxis([axx/100 axx])
%caxis([-axx axx]/4)
colormap(CM(1:170,:))
    ax = gca; pos = ax.Position;
    cb=colorbar;
    cb.Position = [pos(1) pos(2)+pos(4)*1.05 pos(3) .01];
    cb.Orientation = 'horizontal';
    cb.AxisLocation = 'out';
    cb.Color = 'k';
    cb.TickLength = .01;
    cb.Label.String = 'Energy';
    cb.Label.Interpreter='Latex';
    cb.Label.FontSize = figfontsize;
%cb.Ticks = [1e-2 1e-1 1e0];
set(gca,'ydir','normal')
xlabel('Time')
ylabel('Oscillator Index')
set(gca,'colorscale','log')
set(gca,'TickLength',[.01 .1])
grid on

set(gca,'fontsize', figfontsize)
ax = gca;
text(ax.XLim(1),ax.YLim(2)*1.25,'(a)','fontsize',12)

r = rectangle('Position',[t(S1(1)) S2(1) t(S1(end))-t(S1(1)) S2(end)-S2(1)]);% r.FaceColor = [0 0 0 .1]; r.EdgeColor = 'none';
text((max(t(S1))/1.85+min(t(S1))/2)*.8,max(S2)-35,'Output','color','r')
r = rectangle('Position',[t(S10(1)) S20(1) t(S10(end))-t(S10(1)) S20(end)-S20(1)]);%r.FaceColor = [0 0 0 .1];r.EdgeColor = 'none';
text((max(t(S10))/4+min(t(S10))/2),max(S20)-35,'Input','color','r')
plothelper.tickformat

% -------------------------------- ;
% Subplot: ENergy in bands           ;
% -------------------------------- ;
tc = latticedata.tc;
EB = latticedata.Ebandc;
ii = round(linspace(10,80,2));
r = [3:6];
TXT = {'(b)','(c)','(d)','(e)'}
TXT = {'(b)','(b)','(c)','(d)'}
LEG = {'Acoustic Band (A)',' 1st Optical Band (O1)', '2nd Optical Band (O2)', '3rd Optical Band (O3)'};
for k =2:4
    subplot(2,4,r(k)-2)
    set(gca,'fontsize', figfontsize);pause(.01);
    set(gca,'fontsize', figfontsize)

    pltd = abs(EB(:,:,k));
    pltd = movmean(movmean(pltd,10),10,2)';
    imagesc(tc,1:ndof,pltd)
     set(gca,'fontsize', figfontsize)
    %yticks([pi/4 pi/2 3*pi/4 pi]);yticklabels({'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});set(gca,'ydir','normal')
    set(gca,'colorscale','lin')
    colormap(CM(1:170,:))
    caxis([max(pltd(:))/1e2 max(pltd(:)) ])
    ax = gca; pos = ax.Position;
    cb=colorbar;
    cb.Position = [pos(1) pos(2)+pos(4)*1.05 pos(3) .01];
    cb.Orientation = 'horizontal';
    cb.AxisLocation = 'out';
    cb.Color = 'k';
    cb.TickLength = .01;
    cb.Ticks =[.01 .1 1];
    %cb.Label.String = ['Band ',num2str(k)];   
    cb.Label.String = [LEG{k}];
    cb.Label.Interpreter='Latex';
    cb.Label.FontSize = figfontsize;
    set(gca,'YDir','normal')
    ylabel('Oscillator Index')
    xlabel('Time')
    caxis([max(max(max(EB(:,:,:))))/200 max(max(max(EB(:,:,:))))/2])
    set(gca,'colorscale','log')
    % grid on
   % % cb.LineWidth = 1;
    grid on
    %text(tc(round(end/10)),ndof/1.1,['band ',num2str(k)])
set(gca,'fontsize', figfontsize)
    plothelper.tickformat

set(gca,'fontsize', figfontsize)
ax = gca;
text(ax.XLim(1),ax.YLim(2)*1.25,TXT{k},'fontsize',12)
end



%figure(44);clf
% -------------------------------- ;
% Subplot: FFTs and time series    ; 
% -------------------------------- ;
YM = 0;
TXT2 = {'(g)','(h)'}
for SB= 1:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,6+SB)
%subplot(1,2,SB)

set(gca,'fontsize', figfontsize)

n = round(ndof*5/10);
ts = 15000;
ts2 =  size(latticedata.v,1);
plothelper.plotformat
d11=latticedata.v(1:ts,5);          d11(12500:end) = 0;
d12=latticedata.v(1:ts,6);          d12(12500:end) = 0;
d13=latticedata.v(1:ts,7);          d13(12500:end) = 0;
d14 =latticedata.v(1:ts,8);         d14(12500:end) = 0;
asites = 1:4:Ncell;
bsites = 2:4:Ncell;
csites = 3:4:Ncell;;
NN = 2*ndof/4;
d21 = latticedata.v(1+end-ts2:end,33+NN);
d22 = latticedata.v(1+end-ts2:end,33+1+NN);
d23 = latticedata.v(1+end-ts2:end,33+2+NN);
d24 = latticedata.v(1+end-ts2:end,33+3+NN);
Fs = 1/mean(diff(t));
L = length(d11);
F = Fs*(0:(L-1))/L;
F = F * 2* pi;
L2 = ts2;
F2 = Fs*(0:(L2-1))/L2;
F2 = F2 * 2* pi;
hold on
C = linspecer(4);

if SB == 2
plot(F2,abs(fft(d21)),'-','color',[C(1,:) 1],'LineWidth',1)
plot(F2,abs(fft(d22)),'-','color',[C(2,:) 1],'LineWidth',1)
plot(F2,abs(fft(d23)),'-','color',[C(3,:) 1],'LineWidth',1)
plot(F2,abs(fft(d24)),'-','color',[C(4,:) 1],'LineWidth',1)
title('Output Spectra')
%legend('Excitation','Host Mass','Resonator Mass')
else
plot(F,abs(fft(d11)),'color',C(1,:),'LineWidth',1)
plot(F,abs(fft(d12)),'color',C(2,:),'LineWidth',1)
plot(F,abs(fft(d13)),'color',C(3,:),'LineWidth',1)
plot(F,abs(fft(d14)),'color',C(4,:),'LineWidth',1)
title('Input Spectra')
end
ym = max([abs(fft(d11));abs(fft(d12));abs(fft(d13))]);
xlim([10 2400]);
ylim([1e-2 ym])

if SB ==1
    Leg = legend('$\dot{u}_1^{5}$', '$\dot{u}_2^{5}$', '$\dot{u}_3^{5}$','$\dot{u}_4^{5}$');
else
Leg = legend('$\dot{u}_1^{150}$', '$\dot{u}_2^{150}$', '$\dot{u}_3^{150}$','$\dot{u}_4^{150}$');
end
Leg.NumColumns = 1;
ax = gca;POS =ax.Position;
Leg.Position = [POS(1)+POS(3)/1.9 POS(2)+POS(4)/50 .1 .1];
Leg.Color ='none';
Leg.EdgeColor = 'none';
% Input force
Fi = Fin.Ffunction(t(1:ts));
scale = abs(fft(d11));
scale = max(scale);
P=abs(fft(Fi));
%plot(F,envelope(P/max(P)*scale,3,'peaks'),'k-')
set(gca,'yscale','lin')
L01 = min(DISP(:,1));
L02 = max(DISP(:,1));
L1 = min(DISP(:,2));
L2 = max(DISP(:,2));
L3 = min(DISP(:,3));
L4 = max(DISP(:,3));
L5 = min(DISP(:,4));
L6 = max(DISP(:,4));
ax = gca;
Y = ax.YLim;
r = rectangle('Position',[L01 Y(1)*.1  L02-L01  Y(2)*2-Y(1)]);
r.EdgeColor = 'none';r.FaceColor = [CM(end,:) .1]; r.HandleVisibility='off';
r = rectangle('Position',[L1 Y(1)*.1  L2-L1  Y(2)*2-Y(1)]);
r.EdgeColor = 'none';r.FaceColor = [CM(end-20,:) .1]; r.HandleVisibility='off';
r = rectangle('Position',[L3 Y(1)*.1  L4-L3  Y(2)*2-Y(1)]);
r.EdgeColor = 'none';r.FaceColor = [CM(end-40,:) .1]; r.HandleVisibility='off';
r = rectangle('Position',[L5 Y(1)*.1  L6-L5  Y(2)*2-Y(1)]);
r.EdgeColor = 'none';r.FaceColor = [CM(end-60,:) .1]; r.HandleVisibility='off';
xlabel('$\omega$')
ylabel('$|F(\omega)|$')
%set(gca,'children',flipud(get(gca,'children')))
ax = gca;
ax.GridColor= [1 1 1]/1.25e0;ax.GridAlpha = 1;
ax.MinorGridColor=ax.GridColor;ax.MinorGridAlpha = ax.GridAlpha ;
plothelper.tickformat;
ylim(Y)
ax = gca;
text(ax.XLim(1),ax.YLim(2)*1.05,TXT2{SB},'fontsize',12)


pos = ax.Position;
ax3 = axes;
ax3.Position = [pos(1)+pos(3)/2 pos(2)+pos(4)/2 pos(3)/2 pos(4)/2];
if SB == 2
hold on
plot(d21,'-','color',[C(1,:) .5],'LineWidth',1);hold on
plot(d22,'-','color',[C(2,:) .5],'LineWidth',1)
plot(d23,'-','color',[C(3,:) .5],'LineWidth',1)
plot(d24,'-','color',[C(4,:) .5],'LineWidth',1)
% plot(envelope(d21),'color',[C(1,:) 1],'LineWidth',1);hold on;
% plot(envelope(d22),'color',[C(2,:) 1],'LineWidth',1)
% plot(envelope(d23),'color',[C(3,:) 1],'LineWidth',1)
% plot(envelope(d24),'color',[C(4,:) 1],'LineWidth',1)
% plot(-envelope(d21),'color',[C(1,:) 1],'LineWidth',1);hold on;
% plot(-envelope(d22),'color',[C(2,:) 1],'LineWidth',1)
% plot(-envelope(d23),'color',[C(3,:) 1],'LineWidth',1)
% plot(-envelope(d24),'color',[C(4,:) 1],'LineWidth',1)
%title('Time Histories')
xticks([]);%yticks([])
xlabel('Time')
ylabel('Velocity')
else,
hold on
plot(d11,'-','color',[C(1,:) .5],'LineWidth',1);hold on
plot(d12,'-','color',[C(2,:) .5],'LineWidth',1)
plot(d13,'-','color',[C(3,:) .5],'LineWidth',1)
plot(d14,'-','color',[C(4,:) .5],'LineWidth',1)
% plot(envelope(d11),'color',[C(1,:) .3],'LineWidth',1);hold on;
% plot(envelope(d12),'color',[C(2,:) .3],'LineWidth',1)
% plot(envelope(d13),'color',[C(3,:) .3],'LineWidth',1)
% plot(envelope(d14),'color',[C(4,:) .3],'LineWidth',1)
% plot(-envelope(d11),'color',[C(1,:) .3],'LineWidth',1);hold on;
% plot(-envelope(d12),'color',[C(2,:) .3],'LineWidth',1)
% plot(-envelope(d13),'color',[C(3,:) .3],'LineWidth',1)
% plot(-envelope(d14),'color',[C(4,:) .3],'LineWidth',1)
%title('Time Histories')
%title('Time Histories')
xticks([]);%yticks([])xticks([]);%yticks([])
xlabel('Time')
ylabel('Velocity')
end
ym = max([abs((d11));abs((d12));abs((d13))]);
%xlim([10 1400]);
ylim([-ym ym])
ax3.Color = [1 1 1 .15];
plothelper.tickformat

set(gca,'fontsize', figfontsize)
ax = gca;
end

colormap((CM(1:170,:)))

colormap(1-pink.^2)
% set(ax2,'colormap',movmean( [flip(CM(20:end,:)); ones(3,3);flip(1-CM(20:end,:)) ].^1,2))


%set(ax2,'colormap',movmean( [ones(15,3)*0+[0 0 1]; ones(3,3); ones(15,3)*0+[1 0 0] ].^1,2))
ax = gcf;
ax.Position = [4  41  1200  600];
% save2pdf('..\..\..\Josh\Weekly Meetings\11-3-22\figures\paper\VI_plt_all_hertzian')

save2pdf('C:\Users\josht\Box\Josh\Papers\VI Resonator Chain\APS\First Revsion\VI_plt_all_hertzian_R2')

%% Wave reconstruction

% ---------------------- ;
% % Intialize the data     ;
% % ---------------------- ;
% ftdata  = latticedata.v;
% tft     = latticedata.t;
% Lt = size(ftdata,1);
% Lx = size(ftdata,1);
% dt = mean(diff(tft));
% dF = 1/(dt*Lt);
% F = (-(Lt/2)*dF):dF:(Lt/2)*dF;F = F(1:Lt);
% k = linspace(-2*pi,2*pi,Lx);
% [DISP,wavevec] = lattice.dispersion;
% disp = (fft2(ftdata))           ;
% sz = size(disp)                         ;
% 
% 
% figure(3)
% clf
% imagesc(fftshift(abs(disp)));
% ylim([0 400])
% figure(4)
% clf
% for k = 1:6
% 
% if k < 5
% optfrqs = [min(min(DISP(:,k))) max(max(DISP(:,k)))  ]/2/pi;
% E1 =[find( F<-optfrqs(1) & F>-optfrqs(2)) find( F>optfrqs(1) & F<optfrqs(2))];
% end
% % Data to pull
% if k == 5
% optfrqs = [ max(max(DISP(:,2))) min(min(DISP(:,3))) ]/2/pi;
% E1 =[find( F<-optfrqs(1) & F>-optfrqs(2)) find( F>optfrqs(1) & F<optfrqs(2))];
% end
% if k == 6
%     optfrqs = [max(max(DISP(:,4)))  ]/2/pi;
%     E1 =[find( F<-optfrqs(1)) find( F>optfrqs(1) )];
% end
% 
% 
% dispfilt  = zeros(size(disp))                       ; 
% dispshift = fftshift(disp)                          ; 
% dispfilt(E1,:) = (dispshift(E1,:) )                 ; 
% org = ifft2(ifftshift(dispfilt))                    ; 
% 
% % subplot(4,2,k)
% % imagesc(k,2*pi*F,(abs(dispfilt)))
% % colorbar
% % set(gca,'colorscale','log')
% % hold on
% % plot(wavevec/2,DISP,'r')
% % ylim([-max(DISP(:))*0 max(DISP(:))])
% % set(gca,'YDir','normal')
% 
% subplot(2,3,k)
% pltd = real(org);
% pltd = movmean(movmean(pltd,20),5,2);
% imagesc(pltd)
% colorbar
% ax = gca;
% %caxis([max(pltd(:))/100 max(pltd(:))/1 ])
% set(gca,'colorscale','lin')
% %set(gca,'colorscale','log')
% %caxis([max(abs(ftdata(:)))/1e3 max(abs(ftdata(:)))])
% %caxis([-max(abs(ftdata(:))) max(abs(ftdata(:)))]/30)
% caxis([-max(pltd(:)) max(pltd(:))])
% colormap(movmean( [flip(CM(20:end,:)); ones(3,3);flip(1-CM(20:end,:)) ].^1,2))
% 
% %colormap(flip(CM(1:end,:),2))
% %set(ax2,'colormap',movmean( [flip(CM(20:end,:)); ones(3,3);flip(1-CM(20:end,:)) ].^1,2))
% %ax.Position = ax.Position;
% %plothelper.tickformat
% end
% %%
% figure(3)
% clf
% imagesc(ifftshift(abs(dispshift)))
% set(gca,'colorscale','log')
% ylim([0 400])
% 
% %%
% clc
% Eband = latticedata.BandEnergy;
% 
% figure(3422)
% clf
% pltd =Eband(:,:,2);
% pltd = movmean(movmean(pltd,100),10,2);
% 
% imagesc(pltd)
% colormap(movmean( [flip(CM(20:end,:)); ones(3,3);flip(1-CM(20:end,:)) ].^1,2))
% caxis([-max(pltd(:)) max(pltd(:))])
% colorbar
% 
% figure(344)
% clf
% subplot(1,3,1)
% pltd =sum(Eband,3);
% imagesc(pltd)
% caxis([-max(pltd(:)) max(pltd(:))])
% colorbar
% 
% subplot(1,3,2)
% pltd = latticedata.v;
% imagesc(pltd)
% caxis([-max(pltd(:)) max(pltd(:))])
% 
% colorbar
% 
% 
% subplot(1,3,3)
% pltd = latticedata.v - sum(Eband,3);
% imagesc(pltd)
% caxis([-max(pltd(:)) max(pltd(:))])
% colorbar
% 
% colormap(movmean( [flip(CM(20:end,:)); ones(3,3);flip(1-CM(20:end,:)) ].^1,2))
% 
% latticedata.Eband = Eband;
% %%
% 
% lb.Eband = [];
% lb.Ebandc = []
% lb.Process
% 
% 
% %%
% 
% rec = sum(Eband,3);
% figure;imagesc(abs(fft2(rec)));set(gca,'colorscale','log')
%% Getting the inverse of the trnasforms

%   NDOF = size(y,2);
% 
%  Fo = latticedata.decomp_params.wavelet_params.fo;
%  nsamp = size(wb,4);
%  npart = size(wb,3) ;
%  nfourier = NDOF                                ;
%    npt = nfourier/2;
% 
%  freq = Fs*([0:nfourier-1])/nfourier            ;
%  FREQ = Fs*([0:npt-1 npt:-1:1])/nfourier        ;
% 
% 
% 
%   interval_freq = FREQ;
%   a = Fo./interval_freq;
% 
% 
%  fft_MW = conj(bsxfun(@times,pi^(1/4)*(2^0.5)*(exp(-0.5*(2*pi*(bsxfun(@times,FREQ',a)-Fo)).^2)-1*exp(-0.5*(2*pi^2*(bsxfun(@times,FREQ',a).^2+Fo.^2)))),sqrt(a))) ;
% 
% 
%  Patch = wb(:,:,5,20);
% 
%    tff = fft(x_new(:,counter_NDOF),nfourier);
%    noyau2 = bsxfun(@times,fft_MW,tff);
%    resut2 = ifft(noyau2,nfourier);
% 
%    fft_Patch = fft(Patch,nfourier);
%    fft_Patch = fft(Patch,nfourier);
%    fft_mode_2D = fft_Patch./fft_MW;
% 
% 
%    for counter_NDOF = 1:ndofN
%         tff = fft(XN(:,counter_NDOF),nfourier);
%         noyau2 = bsxfun(@times,fft_MW,tff);
%         resut2 = ifft(noyau2,nfourier);
%         Patch = [Mask2' fliplr(Mask2')].*resut2;
%         fft_Patch = fft(Patch,nfourier);
%         fft_mode_2D = fft_Patch./fft_MW;
%         fft_mode = diag(fft_mode_2D);
%         fft_mode(isnan(fft_mode)) = 0;
%         modes(:,counter,counter_NDOF) = ifft(fft_mode,nfourier,'symmetric');
%    end


