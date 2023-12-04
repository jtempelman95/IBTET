clear all
close all
data_3bandImpact = [];


addpath('..\Methods\')
addpath('..\Test_Sim_Codes\')
addpath('..\..\format files\')
addpath('..\')
plothelper.plotformat



% ------------------------------------------------- ;
% Define location to load data                      ;
% ------------------------------------------------- ;
etaidx              = 1;
% config              = 2;datapath        = '..\Data\F22\PrelimStudy_VaryEta_sep7\';figpath = '..\..\..\Josh\Weekly Meetings\09-08-22\figures\HighF';
% config              = 3;datapath        = '..\Data\F22\PrelimStudy_VaryEta_3band\'; figpath = '..\..\..\Josh\Weekly Meetings\09-08-22\figures\3band';
% config              = 2;datapath        = '..\Data\F22\PrelimStudy_VaryEta\'; figpath = '..\..\..\Josh\Weekly Meetings\09-08-22\figures';
config              = 3;datapath        = '..\Data\F22\09-18-22_data\PrelimStudy_VaryEta_3band\'; figpath = '..\..\..\Josh\Weekly Meetings\09-18-22\figures\3band_091322';


% -------------------------------------------------  ;
%  Set directory to save the figure and format plots ;
% -------------------------------------------------  ;
if ~isdir(figpath);mkdir(figpath);end
plothelper.plotformat



% ------------------------------------------------- ;
%            Import all data in parfor              ;
% ------------------------------------------------- ;
clear data amp disp 
for eta             = [.1 .2 .5 1 2 5 10];
    tic
    for fampidx      = 1:25
      
        % ------------------------------------------------- ;
        %            Load Data               ;
        % ------------------------------------------------- ;
        folder          = ['config_',num2str(config),'_eta_',num2str(eta)];
        fname           = ['latticedata_famp',num2str(fampidx)];
        dataloc         = strcat(datapath,folder);
        dataloc_comp    =strcat(datapath,'/compressed/',folder);
        load([dataloc,'/',fname,'.mat']);
     
        % ------------------------------------------------- ;
        %       Proecess and   Compress Data                ;
        % ------------------------------------------------- ;
        latticedatacomp = copy(latticedata);
        latticedatacomp.waveletdata=[];
        latticedatacomp.Process;
        latticedatacomp.compress;
        latticedatacomp.clearTS;
        if~isdir([dataloc_comp]);
        mkdir([dataloc_comp])
        end
        save([dataloc_comp,'/',fname,'.mat'],'latticedatacomp')
        disp(fampidx)
    end
end
