function [f,X]=trfour(x,N,dt)
% devuelve la transformada compleja de Fourier X de la serie real x,
% muestrada con dt intervalo entre muestras tomando N puntos
% si N<length(x), lo trunca y si N>length(x) añade ceros
% también devuelve el vector de frecuencias desde fmin hasta la de Nyquist
% excluida f=0. Por tanto X no incluye la componente DC
if(N/2>floor(N/2));     % si nf es impar
    N=N+1;              % le añade 1 para hacerlo par 
end
X=fft(x,N);
df=1/(dt*(N-1)); % df = 1/ventana
f=(1:N/2-1)*df;
X=X(2:N/2)*dt;   % multiplica por dt para aproximar la transformada continua