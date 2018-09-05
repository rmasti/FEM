%Plotting Script
% Written By: Robert Masti
% 09/01/2018

clc, clear, close all

format long;


%% Geometry Load
xc_g = load('xc_g-0.txt'); % good
yc_g = load('yc_g-0.txt'); % good 
xc = load('xc-0.txt'); % good 
yc = load('yc-0.txt'); % good


%% Areas normal vectors and Volume
Aj = load('Aj-0.txt'); % good
Ai = load('Ai-0.txt'); % good


njx = load('njx-0.txt');
njy = load('njy-0.txt');
nix = load('nix-0.txt'); %check
niy = load('niy-0.txt'); %checked 
Volume = load('volume-0.txt'); %checked

%% Initialized vars

URho = load('V_rho-0.txt');
Uu = load('V_u-0.txt');
Uv = load('V_v-0.txt');
Up = load('V_p-0.txt');
Ubx = load('V_bx-0.txt');
Uby = load('V_by-0.txt');

