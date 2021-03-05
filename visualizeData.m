% Code for visualization purpose
% Supplementary code to the paper:
% "Reference in-vitro dataset for inertial-sensor-to-bone alignment 
% applied to the tibiofemoral joint"

% Dependencies: 
% -------------
% Matlab Aerospace toolbox:
%   https://www.mathworks.com/products/aerospace-toolbox.html
% STL Handling functions:
%   Eric Johnson (2021). 
%   https://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader 
%   Add the following functions to path: stlread (differs from Matlab function), stlbinary

% Data:
% -----


% Note that alignment of sensor axes will yield high errors for the functional
% trials where a segments is fixated and gyroscope measure are close to
% zero. We recomend to use the averaged constant alignments that are provided.

% Copyright (C) 2021 by Ive Weygers.
%% Variables
clear; clc; close all;
frame = 1;                              % The frame that will be plotted
landm_sz = 12;                          % Size of the anatomical landmarks (visual purpose)
bone_color = [0.89 0.855 0.788];        % Color of the bones               (visual purpose)
alignRangeStart = 700;                  % Start range data for alignement of inertial sensor axes (Movement is necessary)
alignRangeStop = 1500;                  % Stop range data for alignement of inertial sensor axes  (Movement is necessary) 

%% 1. Load Data
load('V_30_s_110.mat');                 %  Data-file
load('ct.mat');                         %  Identifed anatomical landmarks
load('align.mat');                      %  Average inertial sensor alignment rotations
% Uncomment to visualize STL (requires toolbox)
% fv_femur = stlread('femur_red.stl');    %  STL 3D surface model femur (reduced vertex count)
% fv_tibia = stlread('tibia_red.stl');    %  STL 3D surface model femur (reduced vertex count)

%% 2. Define 'Common coordinate system' on the base of O1,O2,O3 in CT-scan (M)
[R_MFm, H_FmM] = cartCoSys(CT.femur.o.O1(1:3),CT.femur.o.O2(1:3),CT.femur.o.O3(1:3));
[R_MTm, H_TmM] = cartCoSys(CT.tibia.o.O1(1:3),CT.tibia.o.O2(1:3),CT.tibia.o.O3(1:3));

%% 3. CT-scan (M) -> 'Common coordinate system' on the base of O1,O2,O3 in CT-scan
% Uncomment to visualize STL (requires toolbox)
% [fv_T_m] = vert2int(fv_tibia,R_MTm,H_TmM);  
% [fv_F_m] = vert2int(fv_femur,R_MFm,H_FmM);
FMCC_Fm = R_MFm'*(CT.femur.a.FMCC(1:3)' + H_FmM); FLCC_Fm = R_MFm'*(CT.femur.a.FLCC(1:3)' + H_FmM);
FHC_Fm = R_MFm'*(CT.femur.a.FHC(1:3)' + H_FmM);   FKC_Fm = R_MFm'*(CT.femur.a.FKC' + H_FmM);
FME_Fm = R_MFm'*(CT.femur.a.FME' + H_FmM);        FLE_Fm = R_MFm'*(CT.femur.a.FLE' + H_FmM);
TAC_Tm = R_MTm'*(CT.tibia.a.TAC(1:3)' + H_TmM);   TKC_Tm = R_MTm'*(CT.tibia.a.TKC' + H_TmM);
TLCC_Tm = R_MTm'*(CT.tibia.a.TLCC(1:3)' + H_TmM); TMCC_Tm = R_MTm'*(CT.tibia.a.TMCC(1:3)' + H_TmM);

%% 4.'Common coordinate system' on the base of O1,O2,O3 in CT-scan (M) -> OMC (G)
[R_VFm, H_VFm] = cartCoSys(traj.femur.o.O1(frame,:),traj.femur.o.O2(frame,:),traj.femur.o.O3(frame,:));
[R_VTm, H_VTm] = cartCoSys(traj.tibia.o.O1(frame,:),traj.tibia.o.O2(frame,:),traj.tibia.o.O3(frame,:));
FMCC_V = (R_VFm*FMCC_Fm) - H_VFm;   FLCC_V = (R_VFm*FLCC_Fm) - H_VFm;
FHC_V = (R_VFm*FHC_Fm) - H_VFm;     FKC_V = (R_VFm*FKC_Fm) - H_VFm;
FME_V = (R_VFm*FME_Fm) - H_VFm;     FLE_V = (R_VFm*FLE_Fm) - H_VFm;
TAC_V = (R_VTm*TAC_Tm) - H_VTm;     TKC_V = (R_VTm*TKC_Tm) - H_VTm;
TLCC_V = (R_VTm*TLCC_Tm) - H_VTm;   TMCC_V = (R_VTm*TMCC_Tm) - H_VTm;
% Uncomment to visualize STL (requires toolbox)
% [fv_F_V] = vert2opt(fv_F_m,R_VFm,H_VFm); [fv_T_V] = vert2opt(fv_T_m,R_VTm,H_VTm);

%% 5. Plot data
figure(22); subplot(5,2,[1, 3, 5, 7,9]);  
scat_mF1 = scatter3(traj.femur.o.O1(frame,1),traj.femur.o.O1(frame,2),traj.femur.o.O1(frame,3),landm_sz,'b','filled');hold on;
daspect([1 1 1]);               % SCALING
scat_mF2 = scatter3(traj.femur.o.O2(frame,1),traj.femur.o.O2(frame,2),traj.femur.o.O2(frame,3),landm_sz,'b','filled');hold on;
scat_mF3 = scatter3(traj.femur.o.O3(frame,1),traj.femur.o.O3(frame,2),traj.femur.o.O3(frame,3),landm_sz,'b','filled');hold on;
scat_mF4 = scatter3(traj.femur.o.O4(frame,1),traj.femur.o.O4(frame,2),traj.femur.o.O4(frame,3),landm_sz,'b','filled');hold on;
scat_mT1 = scatter3(traj.tibia.o.O1(frame,1),traj.tibia.o.O1(frame,2),traj.tibia.o.O1(frame,3),landm_sz,'r','filled');hold on;
scat_mT2 = scatter3(traj.tibia.o.O2(frame,1),traj.tibia.o.O2(frame,2),traj.tibia.o.O2(frame,3),landm_sz,'r','filled');hold on;
scat_mT3 = scatter3(traj.tibia.o.O3(frame,1),traj.tibia.o.O3(frame,2),traj.tibia.o.O3(frame,3),landm_sz,'r','filled');hold on;
scat_mT4 = scatter3(traj.tibia.o.O4(frame,1),traj.tibia.o.O4(frame,2),traj.tibia.o.O4(frame,3),landm_sz,'r','filled');hold on; 

% Transformed STL 
% Uncomment to visualize STL (requires toolbox)
% F_patch = patch(fv_F_V,'FaceColor',      bone_color, ...
% 'EdgeColor',       'none',        ...
% 'FaceLighting',    'gouraud',     ...
% 'FaceAlpha',        .4, ...
% 'AmbientStrength', 0.15); hold on;
% T_patch =patch(fv_T_V,'FaceColor',       bone_color, ...
% 'EdgeColor',       'none',        ...
% 'FaceLighting',    'gouraud',     ...
% 'FaceAlpha',        .4, ...
% 'AmbientStrength', 0.15);
% camlight('headlight');
% material('dull');   

% Anatomical landmarks
scat_FMCC = scatter3(FMCC_V(1),FMCC_V(2),FMCC_V(3),landm_sz,'b','filled');hold on;
scat_FLCC = scatter3(FLCC_V(1),FLCC_V(2),FLCC_V(3),landm_sz,'b','filled');hold on;
scat_FKC = scatter3(FKC_V(1),FKC_V(2),FKC_V(3),landm_sz,'b','filled');hold on;
scat_FME = scatter3(FME_V(1),FME_V(2),FME_V(3),landm_sz,'b','filled');hold on;
scat_FLE = scatter3(FLE_V(1),FLE_V(2),FLE_V(3),landm_sz,'b','filled');hold on;
scat_FHC = scatter3(FHC_V(1),FHC_V(2),FHC_V(3),landm_sz,'b','filled');hold on;
scat_TAC =scatter3(TAC_V(1),TAC_V(2),TAC_V(3),landm_sz,'r','filled');hold on;
scat_TKC =scatter3(TKC_V(1),TKC_V(2),TKC_V(3),landm_sz,'r','filled');hold on;
scat_TLCC =scatter3(TLCC_V(1),TLCC_V(2),TLCC_V(3),landm_sz,'r','filled');hold on;
scat_TMCC =scatter3(TMCC_V(1),TMCC_V(2),TMCC_V(3),landm_sz,'r','filled');hold on;
% Base vectors for the femoral and tibia Cartesian coordinate systems
quiv_i = quiver3(TKC_V(1),TKC_V(2),TKC_V(3),kin.tibia.i(frame,1),kin.tibia.i(frame,2),kin.tibia.i(frame,3),70,'b','linewidth',1); hold on;
quiv_j = quiver3(TKC_V(1),TKC_V(2),TKC_V(3),kin.tibia.j(frame,1),kin.tibia.j(frame,2),kin.tibia.j(frame,3),70,'g','linewidth',1); hold on; 
quiv_k = quiver3(TKC_V(1),TKC_V(2),TKC_V(3),kin.tibia.k(frame,1),kin.tibia.k(frame,2),kin.tibia.k(frame,3),70,'r','linewidth',1);hold on;
quiv_I = quiver3(FKC_V(1),FKC_V(2),FKC_V(3),kin.femur.I(frame,1),kin.femur.I(frame,2),kin.femur.I(frame,3),70,'b','linewidth',1); hold on;
quiv_J = quiver3(FKC_V(1),FKC_V(2),FKC_V(3),kin.femur.J(frame,1),kin.femur.J(frame,2),kin.femur.J(frame,3),70,'g','linewidth',1); hold on;
quiv_K = quiver3(FKC_V(1),FKC_V(2),FKC_V(3),kin.femur.K(frame,1),kin.femur.K(frame,2),kin.femur.K(frame,3),70,'r','linewidth',1); hold on;

% ** Sensor alignment femur-attached sensor**
% A) Sensor alignment on inertial data of the current trial
[q_al_f] = findMisalignment(traj.femur.o.O1(alignRangeStart:alignRangeStop,:)',traj.femur.o.O2(alignRangeStart:alignRangeStop,:)',traj.femur.o.O3(alignRangeStart:alignRangeStop,:)', ...
    [imu.femur.Gyr_X(alignRangeStart:alignRangeStop,:) imu.femur.Gyr_Y(alignRangeStart:alignRangeStop,:) imu.femur.Gyr_Z(alignRangeStart:alignRangeStop,:)]');
R = R_VFm*quat2dcm(q_al_f')'; x = R(:,1);y = R(:,2);z = R(:,3);
% B) Use averaged constant misalignment align.femur
% R = R_VFm*quat2dcm(align.femur)'; x = R(:,1);y = R(:,2);z = R(:,3);
quiv_IMU_Femur_b_x = quiver3(traj.femur.o.O1(frame,1),traj.femur.o.O1(frame,2),traj.femur.o.O1(frame,3),x(1),x(2),x(3),50,'r');hold on;
quiv_IMU_Femur_b_y = quiver3(traj.femur.o.O1(frame,1),traj.femur.o.O1(frame,2),traj.femur.o.O1(frame,3),y(1),y(2),y(3),50,'g');hold on;
quiv_IMU_Femur_b_z = quiver3(traj.femur.o.O1(frame,1),traj.femur.o.O1(frame,2),traj.femur.o.O1(frame,3),z(1),z(2),z(3),50,'b');hold on;
% ** Sensor alignment tibia-attached sensor**
% A) Sensor alignment on inertial data of the current trial
[q_al_t] = findMisalignment(traj.tibia.o.O1(alignRangeStart:alignRangeStop,:)',traj.tibia.o.O2(alignRangeStart:alignRangeStop,:)',traj.tibia.o.O3(alignRangeStart:alignRangeStop,:)', ...
    [imu.tibia.Gyr_X(alignRangeStart:alignRangeStop,:) imu.tibia.Gyr_Y(alignRangeStart:alignRangeStop,:) imu.tibia.Gyr_Z(alignRangeStart:alignRangeStop,:)]');
R = R_VTm*quat2dcm(q_al_t')'; x = R(:,1);y = R(:,2);z = R(:,3);
% B) Use averaged constant misalignment align.tibia
% R = R_VTm*quat2dcm(align.tibia)'; x = R(:,1);y = R(:,2);z = R(:,3);   
quiv_IMU_Tibia_b_x = quiver3(traj.tibia.o.O1(frame,1),traj.tibia.o.O1(frame,2),traj.tibia.o.O1(frame,3),x(1),x(2),x(3),50,'r');hold on;
quiv_IMU_Tibia_b_y = quiver3(traj.tibia.o.O1(frame,1),traj.tibia.o.O1(frame,2),traj.tibia.o.O1(frame,3),y(1),y(2),y(3),50,'g');hold on;
quiv_IMU_Tibia_b_z = quiver3(traj.tibia.o.O1(frame,1),traj.tibia.o.O1(frame,2),traj.tibia.o.O1(frame,3),z(1),z(2),z(3),50,'b');hold on;
xlabel('x-pos (mm)'); ylabel('y-pos (mm)'); zlabel('z-pos (mm)');
title('Visualization');

% Reference kinematics
subplot(5,2,2);
endFrame=size(kin.flexion_dh,1);
time = frame/imu.femur.fs:1/imu.femur.fs:endFrame/imu.femur.fs;
plot_flex = plot(time',kin.flexion_dh','k');hold on; 
plot_rot = plot(time',kin.rotation','k-.');hold on;
plot_add = plot(time',kin.abduction','k:');                  
ylim([-50 130]); xlim([frame time(end)]);
xlabel('Time (s)'); ylabel('Angle (deg)');     
title('Tibiofemoral reference kinematics'); 
legend('Flexion','External rotation','abduction');
subplot(5,2,4);
endFrame=size(imu.femur.Acc_X,1);
time = frame/imu.femur.fs:1/imu.femur.fs:endFrame/imu.femur.fs;
plot_Acc_X= plot(time',imu.femur.Acc_X','r');hold on;
plot_Acc_Y= plot(time',imu.femur.Acc_Y','g');hold on;
plot_Acc_Z= plot(time',imu.femur.Acc_Z','b');hold on;
xlim([frame time(end)]); xlabel('Time (s)'); ylabel('(m/s^2)');  
title('Accelerometer measurements S_F'); legend('Acc_X','Acc_Y','Acc_Z');
subplot(5,2,6);
endFrame=size(imu.femur.Acc_X,1);
time = frame/imu.femur.fs:1/imu.femur.fs:endFrame/imu.femur.fs;
plot_Acc_X= plot(time',imu.femur.Gyr_X','r');hold on;
plot_Acc_Y= plot(time',imu.femur.Gyr_Y','g');hold on;
plot_Acc_Z= plot(time',imu.femur.Gyr_Z','b');hold on;
xlim([frame time(end)]); xlabel('Time (s)'); ylabel('(rad/s)');
title('Gyroscope measurements S_F'); legend('Gyr_X','Gyr_Y','Gyr_Z');
subplot(5,2,8);
endFrame=size(traj.femur.a.FHC,1);
time = frame/imu.femur.fs:1/imu.femur.fs:endFrame/imu.femur.fs;
plot_Acc = plot(time',traj.femur.a.FHC');
xlim([frame time(end)]); xlabel('Time (s)'); ylabel('(mm)');     
title('Trajectory FHC'); legend(plot_Acc, 'x-pos','y-pos','z-pos');

%% Functions
function [R, H] = cartCoSys(A,B,C)
% Creates a right handed rotation matrix R and translation vector H 
% from three three-dimentional non-collinear points A,B,C
% Copyright (C) 2021 by Ive Weygers 

    X = (B-A)/norm(B-A);
    Z = cross((C-A),X)/norm(cross((C-A),X));
    Y = cross(X,Z);
        %  Z := cross(X, Y)
        if (norm(cross(X,Y) - Z)> 1e-3)
            Z = -Z;             
        end 
    % Rotation 
    R = [X; Y;Z]'; 
    % Translation
    H = -A';
end
function [fv_trnsf] = vert2int(fv,R,H)
% Transform coordinates from STL-vertices fv into a common coordiante system  
% Input fv: vertices, R: DCM, H: Translation vector
% Output: fv_trnsf: transformed vertices  
% Copyright (C) 2021 by Ive Weygers 

for i=1:size(fv.vertices,1)
    fv_trnsf.vertices(i,:) = R'*(fv.vertices(i,:)' + H);
end
fv_trnsf.faces = fv.faces;
end
function [fv_trnsf] = vert2opt(fv,R,H)
% Transform coordinates from STL-vertices fv into a the optical reference
% coordinate frame G
% Input fv: vertices, R: DCM, H: Translation vector
% Output: fv_trnsf: transformed vertices  
% Copyright (C) 2021 by Ive Weygers 

for i=1:size(fv.vertices,1)       
    fv_trnsf.vertices(i,:) = (R*fv.vertices(i,:)') - H; 
end
fv_trnsf.faces = fv.faces;
end
function [q_JDH] = findMisalignment(I1,I2,I3, gyr)
% Calculates a constant misalignment quaternion from approximated 
% (marker trajectories I1,I2,I3) and measured angular velocities (gyr)
% Copyright (C) 2021 by Ive Weygers 

[R_VM, H_VM] = cartCoSys(I1(:,1),I2(:,1),I3(:,1));
[qVM] = pos2quat(I1',I2',I3');
dt = 0.01;                    % For a sample frequency of 100Hz
wI = gyr';   
% Approximate angular velocities from optical marker trajectories I1,I2,I3
for i = 2:size(wI,1)
    qdM(i,:) = (qVM(:,i)' - qVM(:,i-1)')/dt;                     
    wM_ = 2*quatmultiply(quatconj(qVM(:,i)'),qdM(i,:));  
    wM(i,:) = wM_(2:4);
end
% Search for misalignment
q_JDH = JDH_4_1(wM',wI');
end
function [q] = pos2quat(m1,m2,m3)
% Creates a time-depending quaternion qVM for three given time-depending 
% position vectors m1,m2,m3 
% Copyright (C) 2021 by Ive Weygers

q = zeros(4,size(m1,2));
for a=1:size(m1,1)  
    [R, H] = cartCoSys(m1(a,:),m2(a,:),m3(a,:));
    q(:,a) = dcm2quat(R');                   
end
end
function [q_UV] = JDH_4_1(y_u,y_v)
% Calculates a constant misalignment quaternion q_uv from two vectors y_u, y_v 
% From, Jeroen Hol, Sensor fusion and calibration of inertial sensors, vision, ultra-wideband and GPS 
% Theorem 4.1 (Relative orientation from vector measurements)
% https://www.diva-portal.org/smash/get/diva2:417835/FULLTEXT03
% Copyright (C) 2021 by Ive Weygers

A = zeros(4,4);
for i=1:size(y_u,2)
    A = A + (quatL([0;y_u(:,i)])*quatR([0;y_v(:,i)]));
end

[V,D] = eig(-A);
[~,pos]=max(abs(diag(D)));  
q_UV= V(:,pos);
end
function [qL] = quatL(q)
% Calculates left multiplication matrix of a quaternion 
% https://arxiv.org/abs/1704.06053
% Copyright (C) 2021 by Ive Weygers

if size(q,2)==3
    q = [0 q];
end
    qv = q(2:4);
    q0 = q(1);
    qL = zeros(4,4);

    qL(1,1) = q(1);
    qL(1,2:4) = -qv;
    qL(2:4,1) = qv;
    qL(2:4,2:4) = q0*eye(3) + crossM(q);
end
function [qR] = quatR(q)
% Calculates left multiplication matrix of a quaternion 
% https://arxiv.org/abs/1704.06053
% Copyright (C) 2021 by Ive Weygers

if size(q,2)==3
    q = [0 q];
end
    qv = q(2:4);
    q0 = q(1);
    qR = zeros(4,4);
    qR(1,1) = q(1);
    qR(1,2:4) = -qv;
    qR(2:4,1) = qv;
    qR(2:4,2:4) = q0*eye(3) - crossM(q);
end
function [qx] = crossM(q)
% Computes  matrix cross product of vector or vector part of a quaternion
% input q=[w v]
% Copyright (C) 2021 by Ive Weygers

if(size(q,2)==3)
    q = [0 q];
end
if(size(q,1)==3)
    q = [0 q'];
end

qx = zeros(3,3);
qx(1,2) = -q(4);
qx(1,3) = q(3);
qx(2,1) = q(4);
qx(2,3) = -q(2);
qx(3,1) = -q(3);
qx(3,2) = q(2);
end

