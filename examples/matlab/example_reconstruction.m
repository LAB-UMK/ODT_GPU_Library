%% Example: Reconstruction using ODT_GPU DLL library
% This script demonstrates how to use the GPU-accelerated DLL library
% for optical diffraction tomography (ODT) in MATLAB.
%
% Steps:
% 1. Load sample dataset and parameters
% 2. Load the DLL and header file
% 3. Preprocess hologram and generate K-space
% 4. Run reconstruction (Direct Inverse or Gerchberg–Papoulis)
% 5. Retrieve and display the result
% 6. Unload the library

clear all; clc;

%% 1. Load workspace with data and parameters
load dataAndParams.mat   % provides: data, data_ref, Nx, Ny, nproj, NA, lambda, etc.
%% set number of iterations
nGPi = 100;

%% 2. Load the DLL and header
currentDir = pwd;
currentDir  = fileparts(currentDir);        % go one level up
parentDir  = fileparts(currentDir);        % go next level up
headerPath = fullfile(parentDir, 'include', 'ODT_GPU.h');
dllPath    = fullfile(parentDir, 'bin', 'v-1_0', 'cuda-12_6', 'ODT_GPU.dll');
loadlibrary(dllPath, headerPath);

%% 3. Load and preprocess reference hologram
% Convert MATLAB matrix to C-compatible format (int16)
hologramInt16_ref = int16(convertIntMatrixForC(data_ref));

% Send reference data to GPU
err = calllib('ODT_GPU','HL_addReference', hologramInt16_ref, Nx, Ny, nproj, ...
              NA, lambda, cam_pix, M, n_imm, do_NNC, fftWindowScale);
checkError(err);

%% 4. Preprocess measurement hologram and generate K-space
hologramInt16 = int16(convertIntMatrixForC(data));
Kxy_dim_p = libpointer('int32Ptr', 1);

tic;
err = calllib('ODT_GPU','HL00to02_FromPreprocToGenKO', hologramInt16, Nx, Ny, nproj, ...
              NA, lambda, cam_pix, M, n_imm, do_NNC, ...
              Kxy_dim_p, Kspace_oversampling_z, cos_factor, fftWindowScale, ...
              strcmp(Approx,'Born'));
checkError(err);

% Get K-space dimensions
N_Kspace_xy_padded = Kxy_dim_p.Value;
N_Kspace_z_padded  = round(Kspace_oversampling_z * N_Kspace_xy_padded / 2) * 2;

%% 5. Run reconstruction
% Allocate host memory for the reconstructed volume
pNrec = libpointer('singlePtr', ...
    zeros(2 * N_Kspace_xy_padded * N_Kspace_xy_padded * N_Kspace_z_padded, 1, 'single'));

% Call reconstruction (GP or DI algorithm depending on iteration count)
err = calllib('ODT_GPU','HL03_setParamsAndStartDIandGP', nGPi, do_TC, do_NNC, ...
              relaxM_cuda, betaM, -1, objshift, convertFloatMatrixForC(object_support));
checkError(err);

% Retrieve reconstruction
err = calllib('ODT_GPU','HL04_takeReconstructionAndFreeMemory', pNrec);
checkError(err);
toc;

%% 6. Convert reconstruction back to MATLAB format
n_rec = convertMatrixFromCToMatlab(pNrec.Value, ...
          N_Kspace_xy_padded, N_Kspace_xy_padded, N_Kspace_z_padded);

%% 7. Unload the library
unloadlibrary('ODT_GPU');

%% 8. Display central slice of reconstruction
figure;
imagesc(abs(n_rec(:,:,N_Kspace_xy_padded/2)));
axis image; colormap hot; colorbar;
title('Central slice of reconstructed volume');
