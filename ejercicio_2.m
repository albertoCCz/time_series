% ===========================================================
% Ejercicio de An�lisis y Tratamiento de Datos en Geof�sica y 
% Meteorolog�a. Curso 2021-22.
% ----------------------------
% Autor: Calzada Ch�vez, Alberto.
% Email: alcalzada95@correo.ugr.es
% ===========================================================

%% Ejercicio 2: An�lisis de la se�al meteorol�gica
% ------------------------------------------------
clear, clc

% Definimos el par�metro 'Position' de las figuras
pos = [200 200 900 300];

% Cargamos los datos
load meteo1.mat;

% Definimos vector de tiempo tipo string
td_str = datestr(td+ti);

% Calculamos el intervalo de muestreo
dt = td(2)-td(1);

%% 2.1 Representar la presi�n y la temperatura
% --------------------------------------------
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
ylabel('Temperature [\circC]')
title('Temperature in Granada')

%% 2.2 Diezmar se�al a un periodo de muestreo de 1 hora
% -----------------------------------------------------
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
ylabel('Temperature [\circC]')
title('Temperature in Granada')

%% 2.3 Filtrado paso-alta de la temperatura para eliminar tendencia
%  estacional
% -----------------------------------------------------------------
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
%figure
%freqz(b,a,500,1/dtd)

% Aplicamos el filtro y le sumamos la media a la se�al filtrada y a la
% original
ted_f = filter(b,a,ted)+mted;
ted = ted+mted;

% Representamos la se�al filtrada
figure('Position',pos)
plot(tdd,ted_f)
ylim([-10 40])
xlabel('Time [day]')
ylabel('Temperature [\circC]')
title('Temperature in Granada after seasonal detrending')

%% 2.4 Filtrado paso baja de la temperatura para eliminar la fluctuaci�n
%  diaria
% ----------------------------------------------------------------------
% Calculamos y restamos la media
mted = mean(ted);
ted = ted-mted;

% Definimos un filtro tipo Butterworth para eliminar fluctuaci�n diaria
% Queremos filtrar, por ejemplo, frecuencias por encima de 1 ciclo/semana.
% 1ciclo/semana * 1semana/7d�as = 1/7 ciclos/d�a
fc = 1/7;
fny = 1/dtd/2;
[b,a] = butter(6,fc/fny,'low');

% Representamos la respuesta en frecuencia del filtro
%figure
%freqz(b,a,500,1/dtd)

% Aplicamos el filtro y le sumamos la media a la se�al filtrada y a la
% original
ted_ff = filter(b,a,ted)+mted;
ted = ted+mted;

% Representamos la se�al filtrada
figure('Position',pos)
plot(tdd,ted_ff)
%ylim([-10 40])
xlabel('Time [day]')
ylabel('Temperature [\circC]')
title('Seasonal temperature trend in Granada')

%% 2.5 Envolvente de la fluctuaci�n diaria de la temperatura
% ----------------------------------------------------------
% Calculamos y restamos la media
mted = mean(ted);
ted = ted-mted;

% Creamos un filtro opuesto al del apartado anterior
fc = 1/7;
fny = 1/dtd/2;
[b,a] = butter(6,fc/fny,'high');

% Representamos la respuesta en frecuencia del filtro
%figure
%freqz(b,a,500,1/dtd)

% Aplicamos el filtro y sumamos la media
ted_f = filter(b,a,ted);   % ted ya tiene su media restada

% Calculamos la se�al anal�tica
ted_fa = hilbert(ted_f);

% Calculamos su envolvente
mted_fa = abs(ted_fa);

% Sumamos la media a la se�al filtrada
ted_f = ted_f+mted;

% SUAVIZADO OPCIONAL ------
% =========================
quiere_suavizar = 0; % 1 si queremos suavizar la envolvente; 0 si no queremos.  
if quiere_suavizar
    % Suavizamos la envolvente con filtro butter (no altera la potencia de
    % la se�al en la banda de paso). Filtramos frecuencias mayores de
    % 1ciclo/1.5d�as = 0.667ciclos/d�a
    fc_s = 1/1.5;
    fny = 1/dtd/2;
    [b,a] = butter(4,fc_s/fny,'low');

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
    xlabel('Lag [points ~ 15min]')
    ylabel('Cross-correlation')
    title('Cross-correlation of ted_f and mted_{fa}')

    % Corregimos el lag de la envolvente
    if lag>0 % si la envolvente suavizada se retrasa
        mted_fa = mted_fa(1+abs(lag):end);
    else     % si la envolvente suavizada se adelanta
        mted_fa = mted_fa(1:end-abs(lag));
    end
    
    % Definimos el vector de tiempo para la se�al filtrada como
    tdd_f = tdd(1:end-abs(lag));
else
    % Definimos el vector de tiempo para la se�al filtrada como
    tdd_f = tdd;
end
% =========================

% Representamos la se�al filtrada junto con su envolvente (sumando la
% media de ted, es decir, sumando mted)
figure('Position', pos)
plot(tdd,ted_f,'b')
hold on
plot(tdd_f,-mted_fa+mted,'r','Linewidth',1)
plot(tdd_f, mted_fa+mted,'r','Linewidth',1)
hold off
ylim([-10 40])
legend('ted','\pmmted_{fa} (envelope)')
xlabel('Time [day]')
ylabel('Temperature [\circC]')
title(sprintf('Daily temperature fluctuation and its envelope with fc=%.4f cicles/day = 1/7 cicles/day',fc))

%% 2.6 Transformada de Fourier de la fluctuaci�n diaria de la temperatura
% -----------------------------------------------------------------------
% Calculamos la TF
[f,TED_F] = trfour(ted_f,length(ted_f),dtd);

% Representamos la transformada frente a la frecuencia
figure('Position',pos)
plot(f,TED_F,'b')
xlabel('Frequency [cicles/day]')
ylabel('TED_F')
title('Fourier Transform of the daily temperature fluctuations in Granada')
xlim([0 5.2])

%% 2.7 Predicci�n lineal de la temperatura con un modelo autorregresivo y
%  representar el error de predicci�n
% -----------------------------------------------------------------------
% Calculamos los coeficientes del modelo autorregresivo
[a,g] = lpc(te,20);  % probamos predictor con 20 coeficientes

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
legend_str_1 = sprintf('Mean error = %.4f', mean(te_est_err));
legend_str_2 = sprintf(' Error variance = %.4f', g);
legend(strcat(legend_str, ' [\circC] | ', legend_str_2, ' [\circC]'))
xlim([60 70])
ylim([-10 15])

%% Cerramos todas las figuras
close all