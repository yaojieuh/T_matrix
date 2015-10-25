% rickerwavelet
% 
% Call:
%    [w,t]=rickerwavelet(f0,dt,Nt,t0)
%
% f0: is the central frequency
% dt: is the temporal sample interval
% Nt: is the number of smaples in the wavelet
% t0: The time center of the ricker wavelet
% 
% Example : 
%   f0=100*10^6;
%   dt=4.2456*10^-11;
%   time=0:dt:50*10^-9; % 
%   [w,t]=rickerwavelet(f0,dt);
%   plot(t,w);xlabel('time(s)');ylabel('wavelet amplitude')
%

% 27-01-2010 : TMH : time output, auto calc of Nt and t0 if not set 
%
% Update history
% 27-12-2011 : TyZhu : add zero-phase Ricker wavelet 
% whichRicker =! 0 to chose zero-phase ricker
% 23-05-2012 : TyZhu : add Gaussian wavelet 
%

function [w,time,t0,it0]=rickerwavelet(f0,dt,Nt,t0,whichRicker)

if nargin<3;
    Nt=(3/f0)./dt;
end

if nargin<4
    it0=ceil(Nt/2);
    t0=it0*dt;
else
    it0=floor(t0/dt);
end
if nargin<5;
    whichRicker=1;
    
end
%time=dt:dt:Nt*dt;
%w=(1-2*pi^2*f0^2*(time-t0).^2).*exp(-pi^2*f0^2*(time-t0).^2);
time=(dt:dt:Nt*dt)-t0;w0=2*pi*f0;
if whichRicker==0
    w0=w0/pi; % cutoff frequency(Ty,Jan 2012)
    w=cos(pi*w0*time).*exp(-0.5*w0^2*time.^2);
elseif whichRicker==1
    % w0: center frequency
    w=(1-2*pi^2*f0^2*(time).^2).*exp(-pi^2*f0^2*(time).^2);
elseif whichRicker==2
    % Gaussian
    w=exp(-2*pi^2*f0^2*(time).^2);    
end