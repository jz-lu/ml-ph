clear all
close all
f_size = 22;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

mat = 'MoS2_3.1';
fpath = ['/Users/jonathanlu/Documents/MoS2/coeff' mat '/']; 
fname1 = 'stretch1.txt';
fname2 = 'stretch2.txt';
fname3 = 'shear.txt';
A = importdata([fpath fname1]);
B = importdata([fpath fname2]);
C = importdata([fpath fname3]);

if mat(1) == 'm'
    a0 = sqrt(3)*1.42; 
elseif mat(1) == 'W' 
    a0 = 3.28; 
elseif mat(1) == 'M'
    a0 = 3.122;
elseif mat(1) == 'J' && mat(end) == '1'
    a0 = 3.1;
end 
A0 = sqrt(3)/2*a0^2;

E1 = A.data(:, end);
E2 = B.data(:, end);
E3 = C.data(:, end);
eta = linspace(-0.04, 0.04, 11);

figure;
box on
hold all;
plot(eta, E1, 's');
plot(eta, E2, 's');
plot(eta, E3, 's');

p1 = polyfit(eta, E1', 2);
p2 = polyfit(eta, E2', 2);
p3 = polyfit(eta, E3', 2);

plot(eta, p1(1).*eta.^2 + p1(2).*eta + p1(3));
plot(eta, p2(1).*eta.^2 + p2(2).*eta + p2(3));
plot(eta, p3(1).*eta.^2 + p3(2).*eta + p3(3));

gamma11 = p1(1)/A0*2;
gamma12 = p2(1)/A0-gamma11;
gamma66 = p3(1)/A0/2;

% in eV per cell 
G = gamma66*A0;
K = (gamma11+gamma12)/2*A0;
disp(['G = ' num2str(G) ' eV/cell'])
disp(['K = ' num2str(K) ' eV/cell'])