%% MACRO-DEFINITIONS:

GENERATE_OWN_PATTERN = 0;  % Generates a gain pattern based on the user's location w.r.t the antenna
GENERATE_DOLPH_CHEBYSHEV_PATTERN = 1; % Generates a Dolph-Chebyshev pattern for linear arrays



%% Initial Variable definition: 

c = 3e8;
f0 = 77e9;
BW = 4e9;

% Define the given parameters-locations
lambda = 3e8/f0;
s = lambda/2;          %% Distance between two antennas
Nant = 8;

AP = (0:1:Nant-1)*s;

theta_vals = (-pi/2:0.01:pi/2);

if(GENERATE_OWN_PATTERN)

    r1 = 5;
    theta1 = deg2rad(-30);    %% The Mainlobe direction.

    Ant_Index = 0:1:Nant-1;
    
    d1 = r1;

    H1 = [];

    user1 = r1*exp(1i*theta1);
    
    frequency_axis = f0;

    % An user is assumed to be present at (r1, theta1). We compute the channel
    % matrix between the elements of the Antenna and the user. The Antenna gain
    % pattern is precoded based on theta1 and the magnitude based on channel
    % vector

    for NAntennas = 2:Nant
        d1 = [d1 sqrt((AP(NAntennas) - real(user1))^2 + ((AP(NAntennas) - imag(user1))^2))];
    end

    for NAntennas = 1:Nant
        H1 = [H1 ((3e8./(frequency_axis.*d1(NAntennas))).*(exp(-1i*2*pi*d1(NAntennas)*frequency_axis/c))).'];
    end

    % Estimating the individual antenna gain patterns
    gain_pattern_1 = compute_antenna_gain_pattern(H1.',theta_vals,s,lambda, Nant);

end
if(GENERATE_DOLPH_CHEBYSHEV_PATTERN == 1) % Compute gain pattern from the Dolph-Chebyshev weights
   
    M = Nant;
    sldb = 30;
    coeffs = chebarray(M, sldb).';
    gain_pattern_1 = compute_antenna_gain_pattern(coeffs, theta_vals, s, lambda, Nant);
    
    writematrix(gain_pattern_1, 'gain_pattern.csv');

end
    
%% Plotting the gain pattern of the Antenna
    figure(2), plot(rad2deg(theta_vals),db(gain_pattern_1)),...
        title('Dolph-Chebyshev gain pattern'), xlabel('Orientation (\circ)'),...
        ylabel('Gain (db)');
    axis([-90 90 -50 50]);
    grid on

%% Computing Antenna Gain Pattern here:
function [gain_pattern] = compute_antenna_gain_pattern(H,theta_vals,s,lambda, Nant)
    % IFFT is performed here.
    gain_pattern = exp(1i*2*(pi/lambda)*s*sin(theta_vals).*(0:1:Nant-1)').' * H;
end
