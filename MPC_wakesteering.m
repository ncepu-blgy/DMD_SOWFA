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
    
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    pitchmode=1;
    
    if pitchmode==0
        %dirName={'steps_yaw'}; %directory for identification data: full results exported from SOWFA
        dirName={'/Volumes/NASSIR/Data_Sets/Steps_yaw_low_freq'}; %low frequency simulation
        dirName_val={'steps_yaw_val'}; %directory for identification data: full resulsts exported from SOWFA
    elseif pitchmode==1
        dirName={'steps_theta_YT';}; %directory for identification data: full results exported from SOWFA
        dirName_val={'steps_theta_YT_val'}; %directory for identification data: full resulsts exported from SOWFA
    end
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
        analysis='YawLowFreq/';
        %filename='U_data_complete_vec_yaw.mat';
        filename='/Volumes/NASSIR/Data_Sets/Steps_yaw_low_freq/U_data_complete_vec.mat'; %low freq vector
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
    maindir='/Volumes/NASSIR/MATLAB/'; %define directoty in user's computer to store all results
    maindir=strcat(maindir,analysis);
    
    % IDENTIFICATION DATA
    visualisefirstresults(dirName,rho,1,maindir) %0: skip this; %1: see graphs
    % VALIDATION DATA
    visualisefirstresults(dirName_val,rho,1,maindir) %0: skip this; %1: see graphs
 
    % Make movie from results (identification data, steady state (1500-3000)
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    if videos==1
        %second input argument: dirctory to save movies
        [dirpathwake,dirpathveldef]=makeframes(D,strcat(maindir,'MovieWake'),filename);
        [dirpathfinal]=makefinalframes(D,'/Volumes/NASSIR/MATLAB/Movie_combined');
        [dirpathpowerinsights]=makeframespowerinsights(D,rho,'/Volumes/NASSIR/MATLAB/power_insights');
        [dirpathvelfield]=plotvecfield(D,'/Volumes/NASSIR/MATLAB/power_velfield');
        [dirpathcuthubheightvec]=cuthubheightvec(D,'/Volumes/NASSIR/MATLAB/cuthubheight');
        
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

    detrendingstates=1;
    begin=750;
    beg=750;
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Read and process identification data
    [rotSpeed, nacelleYaw, time1,rotorAzimuth,pitch]=readdmdinformation(dirName); %read information from simulation
    [Inputs, Outputs, Deterministic]=preprocessdmd(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw, pitchmode,pitch ); %preprocess information (resample and only relevant data)
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Read and process validation data
    [rotSpeed_val, nacelleYaw_val, time1_val,rotorAzimuth_val,pitch_val]=readdmdinformation(dirName_val); %read information from simulation
    [Inputs_val, Outputs_val, Deterministic_val]=preprocessdmd(beg, rotSpeed_val,time1_val,rotorAzimuth_val,nacelleYaw_val,pitchmode,pitch_val); %preprocess information (resample and only relevant data)

    % Define states to be used for DMD
    states=QQ_u(:,(begin-beg)+1:end); % define states: first hypothesis 
    %states=[QQ_u(:,(begin-beg)+1:end);QQ_v(:,(begin-beg)+1:end);QQ_w(:,(begin-beg)+1:end)];
    
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    if detrendingstates
        [states,meansteadystate,scalingfactor]=preprocessstates(states);
    else
    end
    
    subpart=5; [f]= MPC_progress(part,subpart,f,{},{}); 
    r=200; %define truncation level for Singular Value Decomposition 
    
    [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd]=dynamicmodedecomposition(states,Inputs, Outputs, Deterministic,3,r,maindir,f); 

%% (3) DATA VALIDATION 
    part=3; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Validate Models from validation data set
    [FITje_val,dirdmd_val]=validatemodels(sys_red, Inputs_val,Outputs_val,r,strcat(dirdmd, '/val'),f);

    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Identification and Validation taks overview
    [modelVAF_val]=idvaloverview(FITje,FITje_val,dirdmd);

%% (4) DYNAMICAL ANALYSIS
    part=4; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    [freq,LambdaDiag, P, phi,damping,b]=dynamicalanalysis(sys_red, U, S,V, dt,X_p,X, method,length(sys_red),0,D,9,Deterministic,r,dirdmd);
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    visualisepodmodes(phi,freq, P,x,y,z,Decimate,D,LambdaDiag,damping,method,Xd,dirdmd) 
    
    podmodes3dim(x,y,z,Xd, P,phi,freq,damping,Decimate,D,dirdmd,f)
    modeanimation([1],x,y,z,P,phi,freq,damping,LambdaDiag, b, Decimate,D,dirdmd,f,1,scalingfactor,meansteadystate,Xd,X)
%% (5) REBUILD FLOW FIELD AND ASSESS DEVIATIONS

    part=5;subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    [statesrebuild]=rebuild(phi,b,LambdaDiag,r,X,Xd); %rebuild states with highest order model
  
    if detrendingstates
        for i=1:size(statesrebuild,2)
            statesrebuild(:,i)=statesrebuild(:,i)*scalingfactor+meansteadystate;
        end
    else
    end
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    X=QQ_u(:,(begin-beg)+1:end-1);
    comparereconstruction(X, statesrebuild,D,dirdmd,x,y,z,Decimate,dirName)
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    evauatemodelerror(X, statesrebuild,D,dirdmd,filename,dirName,x,y,z,Decimate)
    [statesfit]=evaluatetimevaryingerror(X,statesrebuild,dirdmd);
   
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    save(strcat(dirdmd,'/RESULTS.mat'),'sys_red','FITje','FITje_val','X','statesrebuild','statesfit','Inputs','Outputs','Deterministic','Inputs_val','Outputs_val','Deterministic_val');
    close all
    close(f)
%% (6) ECONOMIC MODEL PREDICTIVE CONTROL DESIGN 

    [maxval,modeltouse]=max(FITje_val(2,:));
    [yf1, yf2]=evaluatepredictionpower(sys_red, modeltouse,Inputs_val, Outputs_val,dt, [1 350],400,200);


