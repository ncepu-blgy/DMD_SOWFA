%Economic Model Predictive Control for Wake Steering: an Extended Dynamic Mode Decomposition Approach
%Master Thesis Dissertation
%Author: Nassir Rodrigues Cassamo
%Supervisors: Professor Jan Willem Van Wingerden and Professor João Sousa

%% RELEVAT INFORMATION
%olha aqui um grande testinho
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamica Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data
% (4) Analyses the models, mainly the intrinsic dynamics  modes,
% damping, frequencies, energy) 
% (5) Reconstructs flow based on models an computes deviations from real 
% (6) Designs an Economic Model Predictive Control 

%This script requires the following functions to be added to MATLAB's path
% (1) Functions (LTI toolbox from DCSC, post processing tools from NRL
% (2) cbrewer: color scaling
% (3) altmany-export_fig: to export figures directly to specified directories

%This script requires data in the following fashion
% (1) 2 folders of name (2.1) steps_yaw and (2.2) steps_yaw_val, that
% contain the data non processed directly from CFD simulation from SOWFA
% (2) Data vectors with post processed information: the post processing is
% performed in the cluster and these vectors contain the resampled
% flowfield and resampled grid points. The grid points were resampled every
% Decimate (variable should be in these data vectors)
    % (2.1) U_data_complete_vec
    % (2.2) U_data_complete_vec_val

%% (0) INITIALISE
    %define type of simulation
part=0; subpart=1; [f]= MPC_progress(part,subpart,{},{},{});

    clc
    close all
    addpath('Functions')
    
    maindir='/Volumes/NASSIR/MATLAB/'; %DEFINE MAIN DIRECTORY IN USERS COMPUTER TO STORE ALL RESULTS
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    pitchmode=0; %0 for yaw control and 1 for pith control
    
    if pitchmode==0
        dirName={'steps_yaw'}; %directory for identification data: full results exported from SOWFA
        dirName_val={'steps_yaw_val'}; %directory for identification data: full resulsts exported from SOWFA
    elseif pitchmode==1
        dirName={'steps_theta_YT';}; %directory for identification data: full results exported from SOWFA
        dirName_val={'steps_theta_YT_val'}; %directory for identification data: full resulsts exported from SOWFA
    end
    
    detrendingstates=1; %1 to take mean flow 
    method=3; %0: DMD ; 1:DMDc; 2:ioDMD; 3:extioDMD
    koopman=1;
    videos=0;
    snapshots=0;
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{});
    % Turbine and flow characteristics to be used 
    rho=1.225; %air density in [kg m^-3]
    D=178; %Rotor Diameter used in simulations= 178 [m]
    % Simulation characterisitc (resampling)
    dt=2; %time sampling

%% (1) ASSESS DATA
part=1; subpart=1; [f]= MPC_progress(part,subpart,f,{},{});

    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    if pitchmode==0
        %analysis='Yaw/';
        analysis='Yaw/';
        filename='U_data_complete_vec_yaw.mat';
        filenamevalid='U_data_complete_vec_yaw_val.mat';
        load(filename) ;
        QQ_u=double(QQ_u);
        QQ_v=double(QQ_v);
    elseif pitchmode==1
        analysis='Helix/';
        filename='U_data_complete_vec_pitch.mat';
        load(filename) ;
        QQ_u=full(QQ);
        QQ_u=double(QQ_u);
    end
    
    %define main directory
    maindir=strcat(maindir,analysis);
    
    % IDENTIFICATION DATA
    visualisefirstresults(dirName,rho,0,maindir) %0: skip this; %1: see graphs
    % VALIDATION DATA
    visualisefirstresults(dirName_val,rho,0,maindir) %0: skip this; %1: see graphs
 
    % Make movie from results (identification data, steady state (1500-3000)
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    if videos==1
        %second input argument: dirctory to save movies
        [dirpathwake,dirpathveldef]=makeframes(D,strcat(maindir,'MovieWake'),filename);
        [dirpathfinal]=makefinalframes(D,strcat(maindir,'Movie_combined'),filename);
        [dirpathpowerinsights]=makeframespowerinsights(D,rho,strcat(maindir,'Power_insights'),filename);
        [dirpathvelfield]=plotvecfield(D,strcat(maindir,'power_velfield'),filename);
        [dirpathcuthubheightvec]=cuthubheightvec(D,strcat(maindir,'cuthubheight'),filename);
        
        %Movie making. Please specify Frames Per Second (FPS) inside
        makemoviewake('wake_deflection_yawcontrol_sowfa',dirpathwake);
        makemoviewake('wake_deflection_velocitydefifcitslice_yawcontrol_sowfa',dirpathveldef);
        makemoviewake('wake_deflection_combined',dirpathfinal);
        makemoviewake('wake_deflection_powerdelta',dirpathpowerinsights);
        makemoviewake('wake_deflection_velfield',dirpathvelfield);
        makemoviewake('wake_deflection_hubcut',dirpathcuthubheightvec);
        
        %Converts avi video to mp4. Also specify FPS inside
        movietomp4() 
     else
     end
     
    % Make figures to assess wake evolution (snapshots in sequential time
    % instants)
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    if snapshots==1
        wake_vorticity_deflection
        wakevorticity_secondyaw
        hubheightcut
        cuthubtheightsecond
        vefield5D1
        velfield5D2
    else 
    end
    
    %finding wake center at hub height
    %wakecenter(x,y,z,Decimate, QQ_u,D)
 
%% (2) DYNAMIC MODE DECOMPOSITION 
part=2; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 

    begin=750;
    beg=750;
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Read and process identification data
    [rotSpeed, nacelleYaw, time1,rotorAzimuth,pitch]=readdmdinformation(dirName); %read information from simulation
    [Inputs, Outputs, Deterministic,scalingfactors]=preprocessdmdid(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw, pitchmode,pitch ); %preprocess information (resample and only relevant data)
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Read and process validation data
    [rotSpeed_val, nacelleYaw_val, time1_val,rotorAzimuth_val,pitch_val]=readdmdinformation(dirName_val); %read information from simulation
    [Inputs_val, Outputs_val, Deterministic_val]=preprocessdmdval(beg, rotSpeed_val,time1_val,rotorAzimuth_val,nacelleYaw_val,pitchmode,pitch_val,scalingfactors); %preprocess information (resample and only relevant data)

    % Define states to be used for DMD
    states=QQ_u(:,(begin-beg)+1:end); % define states: first hypothesis 
    n=size(states,1);
    
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    if detrendingstates
        [states,meansteadystate,scalingfactor]=preprocessstates(states);
    else
    end
    
    %include non linear observables - Koopman extensions to better recover
    %non linear dynamics
    if koopman
        [nonlobs]=koopmanstateextension(QQ_u, QQ_v, QQ_w,rho);
        states=[states;nonlobs*ones(1,751)*rho];
    else
    end
    
    subpart=5; [f]= MPC_progress(part,subpart,f,{},{}); 
    r=200; %define truncation level for Singular Value Decomposition 
    
    [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd,xstates]=dynamicmodedecomposition(states,Inputs, Outputs, Deterministic,method,r,maindir,f); 

%% (3) DATA VALIDATION 
    part=3; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Validate Models from validation data set
    [FITje_val,dirdmd_val,xstatesvalid]=validatemodels(sys_red, Inputs_val,Outputs_val,r,strcat(dirdmd, '/val'),f);
    save(strcat(dirdmd,'/FIT.mat'),'FITje_val','FITje');
     
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Identification and Validation taks overview
    [modelVAF_val]=idvaloverview(FITje,FITje_val,dirdmd);
    
    if detrendingstates
        save(strcat(dirdmd,'/RESULTS.mat'),'sys_red',...
            'X','X_p','Xd','Inputs','Outputs','Deterministic','Inputs_val',...
            'Outputs_val','Deterministic_val','U','S','V','dt','r',...
            'dirdmd','method','meansteadystate','scalingfactor','x','y','z',...
            'D','xstates','xstatesvalid','n','Decimate');
    else
        save(strcat(dirdmd,'/RESULTS.mat'),'sys_red',...
            'X','X_p','Xd','Inputs','Outputs','Deterministic','Inputs_val',...
            'Outputs_val','Deterministic_val','U','S','V','dt','r',...
            'dirdmd','method','x','y','z','D','xstates','xstatesvalid','n','Decimate');
    end
    
%% (4) DYNAMICAL ANALYSIS
    part=4; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    [freq,LambdaDiag, P, phi,damping,b]=dynamicalanalysis(sys_red, U, S,V, dt,X_p,X, method,length(sys_red),0,D,9,Deterministic,r,dirdmd,n,Xd);
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    visualisepodmodes(phi,freq, P,x,y,z,Decimate,D,LambdaDiag,damping,method,Xd,dirdmd) 
    
    %podmodes3dim(x,y,z,Xd, P,phi,freq,damping,Decimate,D,dirdmd,f)
    %modeanimation([4],x,y,z,P,phi,freq,damping,LambdaDiag, b, Decimate,D,dirdmd,f,1,scalingfactor,meansteadystate,Xd,X)
%% (5) REBUILD FLOW FIELD AND ASSESS DEVIATIONS
    part=5;subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    [maxval,modeltouse]=max(FITje_val(2,:));
    
    if method==0
        [statesrebuild]=rebuild(phi,b,LambdaDiag,r,X,Xd); %rebuild states with highest order model when no external forcing is done
    else
        model=length(sys_red);
        if isempty(Xd)
            statesrebuild=U(1:n,1:model)*xstates{model}';
            statesrebuildvalid=U(1:n,1:modeltouse)*xstatesvalid{modeltouse}';
        else
            statesrebuild=U(1:n,1:model)*(xstates{model}(:,size(Xd,1)+1:end))';
            statesrebuildvalid=U(1:n,1:modeltouse)*(xstatesvalid{modeltouse}(:,size(Xd,1)+1:end))';
        end
    end
    
    if detrendingstates
        for i=1:size(statesrebuild,2)
            statesrebuild(:,i)=statesrebuild(:,i)*scalingfactor+meansteadystate;
            statesrebuildvalid(:,i)=statesrebuildvalid(:,i)*scalingfactor+meansteadystate;
        end
    else
    end
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    X=QQ_u;
    comparereconstruction(X, statesrebuild,D,dirdmd,x,y,z,Decimate,dirName)
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    valid=load(filenamevalid);
    evauatemodelerror(X, statesrebuild,D,dirdmd,filename,dirName,x,y,z,Decimate)
    [nrmse, nrmsevalid]=evaluatetimevaryingerror(X,statesrebuild,valid.QQ_u, statesrebuildvalid,dirdmd);
    save(strcat(dirdmd,'/RECONSTRUCTION.mat'),'nrmse','nrmsevalid');
    
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    
    close all
    close(f)
%% (6) ECONOMIC MODEL PREDICTIVE CONTROL DESIGN 
 %   
 %   [yf1, yf2]=evaluatepredictionpower(sys_red, modeltouse,Inputs_val, Outputs_val,dt, [1 350],400,200);