function [xvec,yvec] = natCubicSpline(X,Y)
% Different from normal spline by using cubic instead of not-a-knot. That's right, we better than matlab.
%% Initial Formatting
n = length(X);
A_size = 4*n-4;

intervals = n-1;
knots  = n - 2;
A=zeros(A_size,A_size);k=0;j=0;
Y_out = zeros(A_size,1);
j = 2;
Y_out(1) = Y(1);
Y_out(A_size/2) = Y(n);

%% Populate B Output Matrix
for i = 2:1:n-1
        Y_out(j) = Y(i);
        j = j+1;
        Y_out(j) = Y(i);
        j = j+1;
end
j=0;

%% Populate A Matrix Interval Calculations
for i = 1:1:intervals
    A(1+k,1+j) = X(i)^3; A(1+k,2+j) = X(i)^2; A(1+k,3+j) = X(i);A(1+k,4+j)=1;
    A(2+k,1+j) = X(i+1)^3; A(2+k,2+j) = X(i+1)^2; A(2+k,3+j) = X(i+1);A(2+k,4+j)=1;
    k=k+2;j=j+4;
end
j=0;

%% Populate A Matrix Knot Calculations
for i=2:1:knots+1
    A(1+k,1+j) = 3*X(i)^2; A(1+k,2+j) = 2*X(i); A(1+k,3+j) = 1;A(1+k,5+j)= -3*X(i)^2;A(1+k,6+j)= -2*X(i);A(1+k,7+j)= -1;
    A(2+k,1+j) = 6*X(i); A(2+k,2+j) = 2; A(2+k,5+j)= -6*X(i);A(2+k,6+j)= -2;
    k=k+2;j=j+4;
end
A(A_size-1,1) = 6*X(1); A(A_size-1,2) = 2;
A(A_size,A_size-3) = 6*X(n); A(A_size,A_size-2) = 2;

g = size(Y_out);

%% Solve for spline coefficients
if g(1) == 1;
    Coefs = A\transpose(Y_out);
else 
    Coefs = A\Y_out;
end

%% Reformat to ascending order for spline equation
Coef_Matrix = zeros(n-1,4);
k=1;
for i = 1:1:n-1
    Coef_Matrix(i,1) = Coefs(k);
    k=k+1;
    Coef_Matrix(i,2) = Coefs(k);
    k=k+1;
    Coef_Matrix(i,3) = Coefs(k);
    k=k+1;
    Coef_Matrix(i,4) = Coefs(k);
    k=k+1;
end

Coefs_Matrix = Coef_Matrix;

%% Plot Prepping
% This plots all the splines for you on a single graph
x_low = X(1);
x_up = X(n);
y_low = min(Y)-1;
y_up = max(Y)+1;
%axis([x_low x_up y_low y_up]) 
k=1;
p = zeros(1,4);

xvec = [];
yvec = [];
for i = 1:1:intervals
    p(1) = Coefs(k);
    k=k+1;
    p(2) = Coefs(k);
    k=k+1;
    p(3) = Coefs(k);
    k=k+1;
    p(4) = Coefs(k);
    k=k+1;
    x1 = linspace(X(i),X(i+1));
    y1 = polyval(p,x1);
    if i ==1
        y_min = min(y1);
        y_max = max(y1);
    else
        if min(y1)<y_min
            y_min=min(y1);
        elseif max(y1) >y_max
            y_max = max(y1);
        end
    end
    xvec = [xvec x1];
    yvec = [yvec y1];
end
end