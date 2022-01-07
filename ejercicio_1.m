% ===========================================================
% Ejercicio de Análisis y Tratamiento de Datos en Geofísica y 
% Meteorología. Curso 2021-22.
% ----------------------------
% Autor: Calzada Chávez, Alberto.
% Email: alcalzada95@correo.ugr.es
% ===========================================================

% Ejercicio 1: Análisis de señales sísmicas
% ==========================================
%% Ejercicio 1.1: Análisis de la señal de un vibrador sísmico
% ===========================================================
clear, clc

% Definimos el parámetro 'Position' de las figuras
pos = [200 200 900 300];

% Cargamos los datos
load vibroseis.mat;

%% 1.1.1 Representar xw y yr en función del tiempo:
% -------------------------------------------------
% Creamos vectores tiempo
tx = (1:length(xw))*dt;
ty = (1:length(yr))*dt;

% Representamos xw
figure('Position', pos)
plot(tx,xw,'b')
ylim([-1.4 1.4])
xlabel('Time [s]')
ylabel('xw')
title('xw signal')
xlim([0 10])

% Representamos yr
figure('Position', pos)
plot(ty,yr,'Color',[1 0.745 0.039])
xlabel('Time [s]')
ylabel('yr')
title('yr signal')

%% 1.1.2 Trans. Fourier de xw. Representar módulo en función de la
%  frecuencia:
% ----------------------------------------------------------------
% Calculamos Transformada de Fourier con trfour
[f,XW] = trfour(xw,length(xw),dt);

% Calculamos su módulo
mXW = abs(XW);

% Representamos el módulo de XW en función de la frecuencia
figure('Position',[200 200 600 300])
plot(f,mXW,'b')
xlim([0 100])
xlabel('Frequency [Hz]')
ylabel('abs(XW)')
title('Module of the F.T. of xw')

%% 1.1.3 Representación del espectrograma de xw:
% ----------------------------------------------
% Representamos el espectrograma
figure('Position',[200 200 800 400]);
spectrogram(xw,1024,900,1024,1/dt,'yaxis');  % $\Delta f=\sfrac{1}{1024dt}=0.98Hz$
cb = colorbar;
ylabel(cb,'Power/Frequency [dB/Hz]')
title('xw spectrogram')
colormap(parula)

%% 1.1.4 Cálculo de la correlación entre yr y xw:
% -----------------------------------------------
% Calculamos la correlación cruzada
[cc,lags]=xcorr(yr,xw);

% Representamos
figure('Position',[200 200 800 400])
plot(lags,cc,'b')
xlabel('Time [s]')
ylabel('Cross-correlation')
title('Cross-correlation between yr and xw')
xlim([0 4e3])
ylim([-5500 10000])

% Tiempos de viaje de las reflexiones:
% t1 = 1000s
% t2 = 1500s
% t3 = 2300s
% t4 = 3000s

%% Cerramos todas las figuras
close all

%% Ejercicio 1.2: Análisis sismograma vertical del terremoto de Hokaido
% =====================================================================
clear, clc

% Definimos el parámetro 'Position' de las figuras
pos = [200 200 900 300];

% Cargamos los datos
load japon.mat;

% Convertimos el número de cuentas a velocidad del suelo
f_conv = 1500;
wx = wx/f_conv;   % $\mu m/s$

% Convertimos minutos a segundos
t = t*60;
dt = t(2)-t(1);

% Representamos wx
figure('Position', pos)
plot(t,wx)
xlabel('Time [s]')
ylabel('wx [\mum/s]')
title('Vertical ground velocity - Hokaido')

%% 1.2.1 Estimación el rango de periodos del tren de ondas superficiales
%  y diseño de un filtro FIR paso-baja
% ----------------------------------------------------------------------
% Transformada de Fourier
[f,WX] = trfour(wx,length(wx),dt);

% Calculamos el módulo de la TF de wx
mWX = abs(WX);

% Represtamos el módulo de la TF frente a la frecuencia
figure('Position', pos)
plot(f,mWX)
xlabel('Frequency [Hz]')
ylabel('abs(WX)')
title('Module of F.T. of the vertical ground velocity - Hokaido')
xlim([0,0.0012])

% Estimación periodos (y frecuencias asociadas) de ondas superficiales
% a partir de la figura de wx:
% Periodo 1: $t_2-t_1=1.534*10^5 - 1.510*10^5$ -> $T_1 = 2.4*10^3\;s$ -> $f_1 = 4.1667*10^{-4}\;Hz$
% Periodo 2: $t_2-t_1=1.818*10^5 - 1.806*10^5$ -> $T_2 = 1.2*10^3\;s$ -> $f_2 = 8.3333*10^{-4}\;Hz$

% Frecuencia de corte estimada a partir del periodo estimado
fc = 0.0005;     % frecuencia de corte
fny = 1/dt/2;    % frecuencia de Nyquist (=fs/2)

% Un periodo de la frecuencia de corte corresponde a un periodo
% de 1/0.0005 s = 2e4 s, lo que corresponde a 2e4/dt = 1.6667e3 puntos.
% Tenemos que coger un filtro de mayor longitud, por ejemplo de 3e3 puntos.
% ----
% Un periodo de la frecuencia de corte corresponde a un periodo
% de 1/0.0009 s = 1.1111e3 s, lo que corresponde a 1.1111e3/dt = 0.9259e3 puntos.
% Tenemos que coger un filtro de mayor longitud, por ejemplo de 2e3 puntos.
l_filtro = 3e3;

% Creamos el filtro
b = fir1(l_filtro,fc/fny);

% Aplicamos el filtro
wx_f = filter(b,1,wx);

% Representamos la señal y la señal filtrada
figure('Position',pos)
plot(t,wx,'b',t,wx_f,'r')
legend('wx','wx_f')
xlabel('Time [s]')
ylabel('wx | wx_f [\mum/s]')
title(sprintf('wx and wx filtered with fc=%.4f Hz', fc))

% 1.2.2 Representación de la señal filtrada junto son su envolvente
% ------------------------------------------------------------------
% Calculamos la señal analítica
wx_a = hilbert(wx_f);

% Calculamos su envolvente
mwx_a = abs(wx_a);

% Representamos la señal filtrada junto con su envolvente
figure('Position', pos)
plot(t,wx_f,'b')
hold on
plot(t,-1*mwx_a,'r',t,mwx_a,'r','Linewidth',1.5)
hold off
legend('wx_f','\pm mwx_f (envelope)')
xlabel('Time [s]')
ylabel('wx_f | mwx_f [\mum/s]')
title(sprintf('wx_f and its envelope with fc=%.4f Hz', fc))

%% 1.2.3 Diseño filtro Butterworth paso-alta para simular la respuesta
%  sismómetro y representación de la señal filtrada
% ----------------------------------------------------------------------
% Redefinimos la frecuencia de corte
fc = 1/60;  % 1 ciclo/min * 1 min/60s = 1/60 Hz

% Definimos el filtro de Butterworth
[b,a] = butter(4,fc/fny,'high');

% Aplicamos el filtro
wx_fb = filter(b,a,wx);

% Representamos la señal filtrada
figure('Position', pos)
plot(t,wx_fb)
xlabel('Time [s]')
ylabel('wx_{fb} [\mum/s]')
title(sprintf('wx high pass filtered with fc=%.4f Hz', fc))
xlim([0 2.24e5])

%% Cerramos todas las figuras
close all