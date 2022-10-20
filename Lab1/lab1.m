clear

%VARIABLES REGION ---------------------------------------------------------

%load the values from the .txt file
y = load('f14.txt');

%Setting values, that are given in the 
%task's description
T = 5;
dt = 0.01;
df = 1/T;
t = 0:dt:T;
N = length(y);
n = length(t);
f = 0:df:round(n/2) * df;

%END VARIABLES REGION------------------------------------------------------

%FUNCTIONS REGION ---------------------------------------------------------

%"manual" fourier's function calculation
yf = zeros(1, N);
for m = 1:N
  for j = 1:N
    yf(m) = yf(m) + 1/N*y(j)*exp(1)^(-1i*2*pi/N*m*j);
  end
end

%"manual" search for extremums
yf = abs(yf);
extremums = zeros(2,1);
counter = 0;
for j = 3:round(N/2)-1
  if (yf(j) > yf(j+1) && yf(j) > yf(j-1) && abs(yf(j)-yf(j+1)) > 1)
    counter = counter + 1;
    extremums(counter) = j*df;
  end
end

%solving the system of equations
f_sin = sin(2*pi*extremums(1)*t);
A = [sum(t.^6), sum(t.^5), sum(t.^4), sum(f_sin.*t.^3), sum(t.^3);
     sum(t.^5), sum(t.^4), sum(t.^3), sum(f_sin.*t.^2), sum(t.^2);
     sum(t.^4), sum(t.^3), sum(t.^2), sum(f_sin.*t),    sum(t);
     sum(f_sin.*t.^3), sum(f_sin.*t.^2), sum(f_sin.*t), sum(f_sin.*f_sin), sum(N*f_sin);
     sum(t.^3), sum(t.^2), sum(t), sum(N*f_sin), N];
     
c = [sum(y.*t.^3), sum(y.*t.^2), sum(y.*t), sum(y.*f_sin),  sum(y)];

a = inv(A)*c';

disp(a)

%END FUNCTIONS REGION -----------------------------------------------------

%SHOW FIGURES REGION-------------------------------------------------------

%show value's distribution
figure('Name', 'Initial values plot')
plot(t, y), grid

%show Fourier's transform
figure('Name', 'Fourier transform')
plot(abs(yf)), grid

%show Fourier's transform extremums
figure('Name', 'Extremums')
plot(f, abs(yf(1:round(n/2)+1))), grid
