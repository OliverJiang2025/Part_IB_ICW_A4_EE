% A matlab program to plot the theoretical response of the building in the
% 1A vibration lab.
% Based on code written by Penny Cox, now maintained by Aidan Reilly.
%  Tidied up by Jim Woodhouse October 2012
%
close all
clear all

%  USER MUST SET TWO VALUES TO DETERMINE PLOT OPTIONS

% Set value for plottype as follows:
%       1 -> use 'ezplot' to plot a graph on normal axis (in rad/s)
%       2 -> produce a semilog plot (in Hz)
%       3 -> produce a graph on linear axis using 'plot' (in Hz)
plottype=1;

% Set value to be positive to produce plots of the modeshapes
modeshape_visualisation = 0;

m = 1.83; % mass of one floor
L = 0.2; % length
N = 3; % number of degrees of freedom
b = 0.08; % width
E = 210E9; % Young's Modulus
d = 0.001; % thickness
I = b*d*d*d/12; % second moment of area
k = (24*E*I)/(L*L*L); % static stiffness for each floor
M = m*eye(N); % create the mass matrix
K = k*[2 -1 0;-1 2 -1;0 -1 1]; % create the stiffness matrix

% To include vibration absorbers, you will need to modify
%   the mass and stiffness matrices (above)

num_absorbers = 1;  % Set the number of absorbers to zero
absorber_frequency_shift = 0.5; % Small shift in stiffness to slightly change frequency
absorber_mass = 1/num_absorbers;  % Assigning a small mass to each absorber

% Loop over each absorber and extend the mass and stiffness matrices
for i = 1:num_absorbers
    % Create absorber stiffness (tuned slightly for each absorber)
    absorber_stiffness = k * (1 + i * absorber_frequency_shift); % Slightly different stiffness for each absorber
    
    % Extend the mass and stiffness matrices
    M = blkdiag(M, absorber_mass);  % Add a new row and column for each absorber
    K(3, 3) = K(3, 3) + absorber_stiffness;
    new_column = [0; 0; -absorber_stiffness; zeros(N+i-4, 1)];
    K = [K,new_column];  % Add corresponding stiffness
    % K
    new_row = [0,0,-absorber_stiffness,zeros(1, N+i-4),absorber_stiffness];
    new_row
    K = [K; new_row];
    K
end

[V,D] = eig(K,M);
syms w;

for imode=1:(N+num_absorbers)
  freqs(imode) = sqrt(D(imode,imode));
end
freqs

%  Print natural frequencies and mode vectors in command window
%  natural frequencies  3.3933    9.5078   13.7391

hertz = freqs/(2*pi)
modeshapes = V

B = K - ((w*w)*M); 
F = [1; zeros(N + num_absorbers - 1, 1)];
% harmonic solution for unit force at floor 1
disp = B\F;

%start of ezplot section
if (plottype == 1)
  hold on
  
  ifloor=1;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line'),'Color','k')
  ifloor=2;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','g')
  ifloor=3;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','r')

  set(findobj('Type','line','Color','k'),'Color','b')

  set(findobj('Type','line'),'LineStyle','-')
end

% Calculate frequency response functions
all_disp = [];
for w = 1:130
    % w
  B = K - ((w*w)*M); 
  % B
  F = [1; zeros(N + num_absorbers - 1, 1)];
  % harmonic solution for unit force at floor 1
  disp = B\F;
  all_disp = [all_disp disp];
end

w = 1:130;

% Log plot
if (plottype == 2)
  semilogy((w./(2*pi)),abs(all_disp),'-');

% Linear plot
elseif (plottype == 3)
  plot((w./(2*pi)),(all_disp),'-');
end

% Plot modeshapes

if (modeshape_visualisation > 0 )
  V = [0 0 0; V];
  V_ = V + 0.25;
  V = V - 0.25;
  for imode=1:3
    figure
    axis([-5 5 0 3.5])
    title1 = ['Mode ' int2str(imode)];
    title(title1)
    hold on
    plot((V(:,imode)),([0 1 2 3]))
    plot([0 0 0 0],[0 1 2 3],'k')
    plot((V_(:,imode)),([0 1 2 3]))
    for jmode=1:3
      plot([V(jmode+1,imode) V_(jmode+1,imode)],[jmode jmode])
    end
  end
end