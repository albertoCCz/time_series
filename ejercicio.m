% ===========================================================
% Ejercicio de An�lisis y Tratamiento de Datos en Geof�sica y 
% Meteorolog�a. Curso 2021-22.
% ----------------------------
% Autor: Calzada Ch�vez, Alberto.
% Email: albertocalzada95@correo.ugr.es
% ===========================================================

%% Ejercicio 1: An�lisis de se�ales s�smicas
% ==========================================
%% Ejercicio 1.1: An�lisis de la se�al de un vibrador s�smico
% ===========================================================
clear, clc

% Definimos el par�metro 'Position' de las figuras
pos = [200 200 900 300];

% Cargamos los datos
load vibroseis.mat;

% 1.1.1 Representar xw y yr en funci�n del tiempo:
% ------------------------------------------------
% Creamos vectores tiempo
tx = (1:length(xw))*dt;
ty = (1:length(yr))*dt;

% Representamos xw
figure('Position', pos)
plot(tx,xw,'b')
ylim([-1.4 1.4])
xlabel('time [s]')
ylabel('xw')
title('xw signal')

% Representamos yr
figure('Position', pos)
plot(ty,yr,'r')
xlabel('time [s]')
ylabel('yr')
title('yr signal')

% 1.1.2 Trans. Fourier de xw. Representar m�dulo en funci�n de la
% frecuencia:
% ---------------------------------------------------------------
% Calculamos Transformada de Fourier con trfour
[f,XW] = trfour(xw,length(xw),dt);

% Calculamos su m�dulo
mXW = abs(XW);

% Representamos el m�dulo de XW en funci�n de la frecuencia
figure
plot(f,mXW,'b')
xlim([0 100])
xlabel('frequency [Hz]')
ylabel('abs(XW)')
title('Module of the F.T. of wx')

% 1.1.3 Representaci�n del espectrograma de xw:
% ---------------------------------------------
% Representamos el espectrograma
figure('Position',[200 200 800 400]);
spectrogram(xw,1024,900,1024,1/dt,'yaxis');  % delta(f)=1/dt/1024=0.98Hz
cb = colorbar;
ylabel(cb,'Power/Frequency [dB/Hz]')
title('xw spectrogram')
colormap(parula)

% 1.1.4 C�lculo de la correlaci�n entre yr y xw:
% ----------------------------------------------
% Calculamos la correlaci�n cruzada
[cc,lags]=xcorr(yr,xw);

% Representamos
figure('Position',[200 200 800 400])
plot(lags,cc,'b')
xlabel('time offset [s]')
ylabel('Cross-correlation')
title('Cross-correlation between yr and xw')
xlim([0 4e3])
ylim([-5500 10000])

% Tiempos de viaje de las reflexiones:
% t1 = 1000s
% t2 = 1500s
% t3 = 2300s
% t4 = 3000s

close all

%% Ejercicio 1.2: An�lisis sismograma vertical del terremoto de Hokaido
% ====================================================================
clear, clc

% Definimos el par�metro 'Position' de las figuras
pos = [200 200 900 300];

% Cargamos los datos
load japon.mat;

% Convertimos el n�mero de cuentas a velocidad del suelo
f_conv = 1500;
wx = wx/f_conv;   % 10^(-6)m/s

% Convertimos minutos a segundos
t = t*60;
dt = t(2)-t(1);

% Representamos wx
figure('Position', pos)
plot(t,wx)
xlabel('time [s]','Interpreter','latex')
ylabel('ground velocity [$\mu m/s$]','Interpreter','latex')
title('Vertical ground velocity - Hokaido')

% 1.2.1 Estimaci�n el rango de periodos del tren de ondas superficiales
% y dise�o de un filtro FIR paso-baja
% ---------------------------------------------------------------------
% Transformada de Fourier
[f,WX] = trfour(wx,length(wx),dt);

% Calculamos el m�dulo de la TF de wx
mWX = abs(WX);

% Represtamos el m�dulo de la TF frente a la frecuencia
figure('Position', pos)
plot(f,mWX)
xlabel('Frequency [Hz]')
ylabel('abs(WX)')
title('Module of F.T. of the vertical ground velocity - Hokaido')
xlim([0,0.0012])

% Estimaci�n periodos (y frecuencias asociadas) de ondas superficiales
% a partir de la figura de wx:
% Periodo 1: 1.534e5 - 1.510e5 -> T = 2.4e3 [s] -> f = 4.1667-4 [Hz]
% Periodo 2: 1.818e5 - 1.806e5 -> T = 1.2e3 [s] -> f = 8.3333-4 [Hz]

% Frecuencia de corte estimada a partir del periodo estimado
fc = 0.0005;     % frecuencia de corte
fny = 1/dt/2;    % frecuencia de Nyquist (=fs/2)

% Un periodo de la frecuencia de corte corresponde a un periodo
% de 1/0.0005 s = 2e4 s, lo que corresponde a 2e4/dt = 1.6667e3 puntos.
% Tenemos que coger un filtro de mayor longitud, por ejemplo de 3e3 puntos.
% ----
% Un periodo de la frecuencia de corte corresponde a un periodo
% de 1/0.0009 s = 1.1111e3 s, lo que corresponde a 2e4/dt = 0.9259e3 puntos.
% Tenemos que coger un filtro de mayor longitud, por ejemplo de 2e3 puntos.
l_filtro = 3e3;

% Creamos el filtro
b = fir1(l_filtro,fc/fny);

% Aplicamos el filtro
wx_f = filter(b,1,wx);

% Representamos la se�al y la se�al filtrada
figure('Position',pos)
plot(t,wx,'b',t,wx_f,'r')
legend('wx','wx_f')
xlabel('time [s]','Interpreter','latex')
ylabel('signals [$\mu m/s$]','Interpreter','latex')
title(sprintf('wx and wx filtered with fc=%.4f Hz', fc))

% 1.2.2 Representaci�n de la se�al filtrada junto son su envolvente
% -----------------------------------------------------------------
% Calculamos la se�al anal�tica
wx_a = hilbert(wx_f);

% Calculamos su envolvente
mwx_a = abs(wx_a);

% Representamos la se�al filtrada junto con su envolvente
figure('Position', pos)
plot(t,wx_f,'b')
hold on
plot(t,-1*mwx_a,'r',t,mwx_a,'r','Linewidth',1.5)
hold off
legend('wx_f','\pm mwx_f (envelope)')
xlabel('time [s]','Interpreter','latex')
ylabel('signals [$\mu m/s$]','Interpreter','latex')
title(sprintf('wx_f and its envelope with fc=%.4f Hz', fc))

% 1.2.3 Dise�o filtro Butterworth paso-alta para simulacar la respuesta
% sism�metro y representaci�n de la se�al filtrada
% ---------------------------------------------------------------------
% Redefinimos la frecuencia de corte
fc = 1/60;  % 1 ciclo/min * 1 min/60s = 1/60 Hz

% Definimos el filtro de Butterworth
[b,a] = butter(4,fc/fny,'high');

% Aplicamos el filtro
wx_fb = filter(b,a,wx);

% Representamos la se�al filtrada
figure('Position', pos)
plot(t,wx_fb)
xlabel('time [s]','Interpreter','latex')
ylabel('$wx_{fb}$ [$\mu m/s$]','Interpreter','latex')
title(sprintf('wx high pass filtered with fc=%.4f Hz', fc))

close all

%% Ejercicio 2: An�lisis de la se�al meteorol�gica
% ------------------------------------------------
clear, clc

% Definimos el par�metro 'Position' de las figuras
pos = [200 200 900 300];

% Cargamos los datos
load meteo1.mat;

% Definimos vector de tiempo tipo string
ts = datestr(ti+td); % en d�as

% Calculamos el intervalo de muestreo
dt = td(2)-td(1);

% 2.1 Representar la presi�n y la temperatura
% -------------------------------------------
% Representamos la presi�n (izqda.) y la temperatura (dcha.)
figure('Position',pos)
plot(td,p)
xlabel('Time [day]')
ylabel('Pressure [mbar]')
title('Pressure in Granada')

% Representamos la temperatura
figure('Position',pos)
plot(td,te)
xlabel('Time [day]')
ylabel('Temperature [\circ C]')
title('Temperature in Granada')

% 2.2 Diezmar se�al a un periodo de muestreo de 1 hora
% ----------------------------------------------------
% (Entiendo que se refiere a la temperatura)
% Diezmamos la temperatura a intervalos de 1 hora. Puesto que la tenemos en
% intervalos de 15min -> submuestreamos a 1/4
ted = decimate(te,4);    % Submuestreo temperatura
tdd = downsample(td,4);  % Submuestreo tiempo

% Calculamos el nuevo intervalo de muestreo
dtd = tdd(2)-tdd(1);

figure('Position',pos)
plot(tdd,ted)
xlabel('Time [day]')
ylabel('Temperature [\circ C]')
title('Temperature in Granada')

% 2.3 Filtrado paso-alta de la temperatura para eliminar tendencia
% estacional
% ----------------------------------------------------------------
% Calculamos y restamos la media
mted = mean(ted);
ted = ted-mted;

% Definimos un filtro tipo Butterworth para eliminar tendencia estacional
% Queremos filtrar, por ejemplo, frecuencias por debajo de 1 ciclo/mes.
% 1ciclo/mes * 1mes/30d�as = 1/30 ciclos/d�a
fc = 1/30;
fny = 1/dtd/2;
[b,a] = butter(6,fc/fny,'high');

% Representamos la respuesta en frecuencia del filtro
figure
freqz(b,a,500,1/dtd)

% Aplicamos el filtro y le sumamos la media
ted_f = filter(b,a,ted)+mted;

% Representamos la se�al filtrada
figure('Position',pos)
plot(tdd,ted_f)
ylim([-10 40])
xlabel('Time [day]')
ylabel('Temperature [\circ C]')
title('Temperature in Granada after seasonal detrending')

% 2.4 Filtrado paso baja de la temperatura para eliminar la fluctuaci�n
% diaria
% ---------------------------------------------------------------------
% Calculamos y restamos la media
mted_f = mean(ted_f);
ted_ff = ted_f-mted_f;

% Definimos un filtro tipo Butterworth para eliminar fluctuaci�n diaria
% Queremos filtrar, por ejemplo, frecuencias por encima de 1 ciclo/semana.
% 1ciclo/semana * 1semana/7d�as = 1/7 ciclos/d�a
fc = 1/7;
fny = 1/dtd/2;
[b,a] = butter(6,fc/fny,'low');

% Representamos la respuesta en frecuencia del filtro
figure
freqz(b,a,500,1/dtd)

% Aplicamos el filtro y le sumamos la media
ted_ff = filter(b,a,ted_ff)+mted_f;

% Representamos la se�al filtrada
figure('Position',pos)
plot(tdd,ted_ff)
ylim([-10 40])
xlabel('Time [day]')
ylabel('Temperature [\circ C]')
title('Temperature in Granada after seasonal and daily detreding')

% 2.5 Envolvente de la fluctuaci�n diaria de la temperatura
% ---------------------------------------------------------
% Creamos un filtro opuesto al del apartado anterior
fc = 1/7;
fny = 1/dtd/2;
[b,a] = butter(8,fc/fny,'high');

% Representamos la respuesta en frecuencia del filtro
figure
freqz(b,a,500,1/dtd)

% Aplicamos el filtro y sumamos la media
ted_f = filter(b,a,ted);

% Calculamos la se�al anal�tica
ted_fa = hilbert(ted_f);   % Volvemos a utilizar la se�al no diezmada

% Calculamos su envolvente
mted_fa = abs(ted_fa);

% Sumamos la media a la se�al filtrada
ted_f = ted_f+mted;

% Suavizamos la envolvente con filtro butter (no altera la potencia de
% la se�al en la banda de paso)
fc = 1/2;   % filtramos frecuencias mayores de 1ciclo/2d�as=0.5ciclos/d�a
fny = 1/dtd/2;
[b,a] = butter(4,fc/fny,'low');

% Aplicamos el filtro butter
mted_fa = filter(b,a,mted_fa);

% Estudiamos la correlaci�n entre la fluctuaci�n diaria y su envolvente
% para corregir el desfase introducido por el filtro
[cc,lags] = xcorr(ted_f,mted_fa);

% Print lag de la correlaci�n m�xima
[mccv,mcci] = max(cc);
lag = lags(mcci);
sprintf('lag of max cross-correlation: %d puntos',lag)

% Representamos la correlaci�n
figure('Position',pos)
plot(lags,cc)
xlabel('Time lag [points]')
ylabel('Cross-correlation')
title('Cross-correlation of ted_f and mted_{fa}')

% Corregimos el lag de la envolvente
mted_fa = mted_fa(abs(lag)+1:end);

% Representamos la se�al filtrada junto con su envolvente (sumando la
% media de ted, es decir, sumando mted) corrigiendo el lag
figure('Position', pos)
plot(tdd,ted_f,'b')
hold on
plot(tdd(1:end-abs(lag)),-mted_fa+mted,'r','Linewidth',1,...
     tdd(1:end-abs(lag)), mted_fa+mted,'r','Linewidth',1)
legend('ted','\pmmted_{fa} (envelope)')
xlabel('Time [day]')
ylabel('Temperature [\circC]')
title(sprintf('Daily temperature fluctuation and its envelope with fc=%.4f cicles/day = 1/7 cicles/day', fc))

% 2.6 Transformada de Fourier de la fluctuaci�n diaria de la temperatura
% ----------------------------------------------------------------------
% Calculamos la TF
[f,TED_F] = trfour(ted_f,length(ted_f),dtd);

% Representamos la transformada frente a la frecuencia
figure('Position',pos)
plot(f,TED_F,'b')
xlabel('Frequency [cicles/day]')
ylabel('TED_F')
title('Fourier Transform of the temperature fluctuations in Granada')
xlim([0 5.2])

% 2.7 Predicci�n lineal de la temperatura con un modelo autorregresivo y
% representar el error de predicci�n
% ----------------------------------------------------------------------
% Calculamos los coeficientes del modelo autorregresivo
a = lpc(te,20);  % probamos predictor con 20 coeficientes

% Realizamos la predicci�n
te_est = filter([0 -a(2:end)],1,te);

% Calculamos el error de predicci�n
te_est_err = te-te_est;  % Error=diferencia entre se�al real y estimada

% Representamos la temperatura, la temperatura estimada y el error de
% predicci�n en la misma figura
figure('Position',[100 100 900 500])
subplot(4,1,1:3)
plot(td,te,'b')
hold on
plot(td,te_est,'Color',[0.8500 0.3250 0.0980])
hold off
legend('te','te_{est}')
ylabel('te | te_{est} [\circC]')
title('Temperature, its predicted value and the prediction error')
xlim([60 70])

subplot(4,1,4)
plot(td,te_est_err,'g')
xlabel('Time [day]')
ylabel('te_{est_ err} [\circC]')
xlim([60 70])

%%
close all
