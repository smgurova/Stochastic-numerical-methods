% estimate the probability to obtain 8 or
% more heads, if a coin is tossed 10 times
% the real value is approx. 0.0547

clear, clc, close all
% initialization
c = 0; % number of head of the coin
N = 20e3; % number steps
h = waitbar(0, 'Please wait...');
% run the simulation
for n = 1:N
    
    % generate a set of random numbers x(0, 1)
    x = rand(1, 10);
        % check if the event occurs
    if length(x(x>0.5)) >= 8
                % find the number of occurrences
        c = c + 1;
        end 
       % waitbar
   waitbar(n/N, h, ['Number of the test: ' num2str(n)])
    end
% close the waitbar
close(h)
% estimate the probability to obtain 8 or
% more heads, if a coin is tossed 10 times
% the real value is approx. 0.0547
Pest = c/N;
err = (Pest - 0.0547)/0.0547*100;
disp(['The estimated value is ' num2str(Pest)])
disp(['The error of the estimation is ' num2str(err) ' %'])
commandwindow