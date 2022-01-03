function [X,XM,XP,f]=tfour(x,fm,nf)
% calcula la FFT de x (vector columna)
% nf es el número de puntos de la transformada
% figure
% hold all;
dt=1/fm;
N=length(x);
w=hamming(N);           % ventana Hamming
xw=x.*w;                % aplica ventana
X=fft(xw,nf)*dt;        % multiplica por dt para tener la TF 'continua'
X=X(1:nf/2);            % primera mitad
XM=abs(X);XP=angle(X);  % módulo y fase
f=fm*(0:nf/2-1)/nf;     % frecuencias
f=f';                   % no incluye la frecuencia de Nyquist, pero si 0
                        % el primer valor X(1) debe ser N*mean(xw)*dt por
                        % la normalización que usa Matlab
