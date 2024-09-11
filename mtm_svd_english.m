
clear all

% Directorio con rutinas mtm-svd y datos a utilizar
% Directory with routines for mtm-svd and data to use
path1 = '/Users/sck117/Box Sync/MM_NSF2018/MTM-SVD/';
path(path,[path1,'Code']);  % incorpora el codigo al path/the path to the code
path2 = [path1,'Data/SLA/'];% Camino de archivos con datos de nivel del mar/file path with sea %level data
path3 = [path1,'Data/WIN/'];% camino de archivos con datos de viento/file path with wind data
path4 = [path1,'Output/'];% camino para guardar archivos de resultados/file path for Output

%--------------------------------------------------------------------------
% Lectura de los datos/Reading the data
%
% los datos deben ser interpolados antes de comenzar el analisis
% the data must be interpolated before starting the analysis
%--------------------------------------------------------------------------

% lectura de las matrices de fecha(fec) Latitud(lat) y Longitud(lon)
% reading the matrices of date(fec), Latitude(lat), and Longitude(lon)
%
% Todos los archivos tienen la misma dimension y resolucion espacial
% fec corresponde a la fecha, en formato matlab, de cada uno de los campos
% satelitales (lat,lon) de nivel del mar, obtenidos semanalmente durante
% el periodo 1992-2011
% All files have the same spatial resolution and dimension
% fec corresponds to the date, in matlab format, of each of the fields
% satellite (lat, lon) of sea level, obtained weekly during
% the period 1992-2011
%
load ([path2,'fecha_lat_lon.mat']); 

% Chequeo de lectura de un archivo de anomalias de nivel del mar(sla)
% Reading check of a file of sea level anomalies (sla)
%load ([path2,'sla_19921014.mat']);
%
% Visualizacion de los datos de un dia
% plot of the one day of data
%figure; pcolor(lon, lat, sla); shading flat; 
%title(['SLA ' datestr(fec(1))])

% Construccion de la matriz espacio-tiempo/Construction of the space-time Matrix
%-----------------------------------------

% lectura de un archivo de datos de sla para obtener la posicion de los
% datos validos (datos sobre el oceano)
% Reading one sla data file to obtain the positions of the valid data (data about the ocean)
load ([path2,'sla_19921014.mat']);

% Indices de las posiciones con datos validos/indices of the valid data point)
id = find(~isnan(sla));

% Pre-definir el tamaño final de la matriz 2D [tiempo-espacio]
% (esto permite ahorrar tiempo de computo)
% Pre-define the final size of the 2D matrix [time-space]
% (this saves computing time)
H = nan(length(fec), length(id));

% lectura secuencial de los datos sin brechas
% Sequential reading of the data without gaps
  for i = 1:length(fec); 
      sla = load ([path2,'sla_',datestr(fec(i), 'yyyymmdd'),'.mat']);
      sla=sla.sla;
      H(i, :) = sla(id)';
  end
  
% Verifica que los datos han sido leidos/Verify that the data has been read
% Grafia de las primeras 100 posiciones en la matriz tiempo-espacio/plot the first 100 points
% figure: pcolor(1:100, fec, H(:,1:100)); shading flat 
% datetick('y')
% ylabel('Tiempo')
% xlabel('Posiciones')

% Guarda la Matriz 2D /save the 2D matrix
save ([path4,'SLA_2D'], 'H','id','lon','lat','fec')


%% ------------------------------------------------------------------------
% Construccion de la Matriz 2D ([tiempo, espacio])
% del esfuerzo del viento
%los datos deben ser interpolados antes de comenzar el analisis
% Construction of the 2D Matrix ([time,space]) of wind speed
% the data must be interpolated prior to analysis
%--------------------------------------------------------------------------

% lectura de archivos con fecha (fec) Latitud(lat) y Longitud(lon)
% Reading the file with date (fec), latitude (lat), and longitude (lon)
load ([path2,'fecha_lat_lon.mat']); 

% lectura del primer archivo de vientos/reading the first wind file
load ([path3,'win_19921014.mat']);

% Chequeo de lectura de un archivo del esfuerzo del viento
% Reading check of a wind file
% es=sqrt(eu.^2+ev.^2); % magnitud del esfuerzo del viento/wind speed
% figure; pcolor(lon, lat, es); shading flat; colorbar
% hold on; 
% quiver(lon, lat, eu, ev,1.5, 'k'); %plots arrows of wind

% Indices de posiciones con datos validos./Indices of valid data
% Como los campos de SLA y vientos son del mismo tamaño y poseen igual
% resolucion espacial (1/4°), podemos utilizar el indice de datos valios de
% SLA para construir la matriz 2d de vientos, y asi asegurarnos de utilizar
% las mismas posiciones en los analisis
% As the fields of SLA and winds are the same size and have the same
% spatial resolution (1/4 °), we can use the valuation data index of
% SLA to build the 2d matrix of winds, and so make sure to use
% the same positions in the analyzes

load ([path4,'SLA_2D'], 'id');      % id = find(~isnan(es));

% Pre-definir el tamaño final de la matriz 2D [tiempo-espacio]
% (esto permite ahorrar tiempo de computo)
% Pre-define the final size of the 2D Matrix [time-space] (this speeds up computing time)
EU = nan(length(fec), length(id));
EV = EU;

% lectura secuenacial de los datos semanales sin brechas
% Reading weekly data without gaps
for i = 1:length(fec)
      win = load ([path3,'win_',datestr(fec(i), 'yyyymmdd'),'.mat']);
      eu=win.eu;
      EU(i, :) = eu(id)';
      ev=win.ev;
      EV(i, :) = ev(id)';
  end
  
% Verifica que los datos han sido leidos/Verify the data have been read
% Grafia de las primeras 100 posiciones en la matriz tiempo-espacio
% plot the first 100 points of the Matrix
%  figure 
%  subplot(1,2,1)
%  pcolor(1:100, fec, EU(:,1:100)); shading flat ; colorbar
%  datetick('y'); 
%  title('Esfuerzo Viento Zonal') 
%  ylabel('Tiempo')
%  xlabel('Posiciones')
% 
%  subplot(1,2,2)
%  pcolor(1:100, fec, EV(:,1:100)); shading flat ; colorbar
%  datetick('y'); 
%  title('Esfuerzo Viento Meridional') 
%  ylabel('Tiempo')
%  xlabel('Posiciones')

% Guarda la Matriz 2D/ Save the 2D Matrix
  save ([path4, 'WIN_2D'],'EU','EV','id','lon','lat','fec')

%% ------------------------------------------------------------------------
% Reconstruccion Campos y Caracterisiticas medias de los datos
% ejemplos para manipular las matrices generadas anteriormente
% Reconstruction Fields and Average data characteristics
% examples to manipulate the matrices generated previously
%--------------------------------------------------------------------------

clear all

%---------------------
% Seteo de directorios/ Set the directories
%---------------------

% Directorio con rutinas mtm-svd y datos a utilizar
% Directory with mtm-svd routines and data
path1 = '/Users/sck117/Box Sync/MM_NSF2018/MTM-SVD/';
path(path,[path1,'Code']);  % incorpora el codigo al path/path to code
path4 = [path1,'Output/'];% camino para archivos de datos /path for data files

 
% lectura de los datos de sla en formato 2D [tempop-espacio]
% ademas de la fecha(fec), Latitud(lat) y Longitud(lon)
% reading the 2D SLA data [time-space], fields are date (fec), latitude (lat), and Longitude (lon)
load ([path4,'SLA_2D.mat']);

% construye marices regulares de coordenadas lon, lat
% de igual tamaño que el campo original
% construct regular lat/lon matrices of equal size as the original field
[x,y]=meshgrid(lon, lat); 


%-----------------------------------------------------------------
% encontrando el indice del dia mas cercano al 20 de julio de 2010
% Finding the index of the closest day to 20 July 2010
%-----------------------------------------------------------------
dia=datenum(2010,7,20)
idia=near(fec,dia)

% dia mas cercano al 20 de julio de 2010/closest day to 20 July 2010
datestr(fec(idia))

% reconstruyendo campo nivel del mar mas cercano al 20 de julio de 2010
% rebuilding the sea level field closest to July 20, 2010
sla=x*nan;          % generando la matriz vacia
sla(id)=H(idia,:);  % asignando los valores validos del campo a la matriz

% Graficnado campo mas cercano al 20 de julio de 2010/plotting the data closest to 20 July 2010
figure;
pcolor(lon, lat, sla); shading flat; colorbar
title(['SLA ' datestr(fec(idia))])

%-----------------------------------------------------------------
% Promedio y desviacion estandar de las anomalias del nivel del mar
% Mean and stdev of SLA
%-----------------------------------------------------------------
msla = mean(H);         % promedio temporal
MSLA = x*nan;           % generando la matriz vacia
MSLA(id) = msla;        % asignando los valores validos del campo a la matriz

dsla = std(H);          % desviacion estandar temporal
DSLA = x*nan;           % generando la matriz vacia
DSLA(id) = dsla;        % asignando los valores validos del campo a la matriz

% Graficando promedio y desviacion/Plotting the mean and stdev
figure
subplot(1,2,1)
pcolor(lon, lat, MSLA); shading flat; colorbar
title('Promedio Anomalias Nivel del Mar/Mean Sea Level Anomalies')
subplot(1,2,2)
pcolor(lon, lat, DSLA); shading flat; colorbar
title('Desviacion Estandar Anomalias Nivel del Mar/Standard Deviation SLA')

%-----------------------------------------------------------------
% Promedio y desviacion estandar del esfuerzo del viento/ mean and stdev of wind
%-----------------------------------------------------------------
load ([path4,'WIN_2D.mat']);
% promedios temporales/temporary means
meu = mean(EU);  
mev = mean(EV);     

% generando la matriz vacia y asignando los valores promedios
% generating the empty matrix and assigning mean values
MEU = x*nan; MEU(id) = meu; 
MEV = x*nan; MEV(id) = mev; 

% desviacion estandar temporal/temporary stdev
deu = std(EU);
dev = std(EV);

% generando la matriz vacia y asignando los valores de desviacion
% generating the empty matrix and assigning stdev values
DEU = x*nan; DEU(id) = deu; 
DEV = x*nan; DEV(id) = dev; 


% Graficando promedio y desviacion/plotting mean and stdev
figure
subplot(2,2,1)
pcolor(lon, lat, MEU); shading flat; colorbar
title('Promedio/Mean')
ylabel('Esfuerzo viento zonal/Zonal Wind')
subplot(2,2,2)
pcolor(lon, lat, DEU); shading flat; colorbar
title('Desviacion Estandar/Standard Deviation')

subplot(2,2,3)
pcolor(lon, lat, MEV); shading flat; colorbar
title('Promedio/Mean')
ylabel('Esfuerzo viento meridional/Meridional Wind')
subplot(2,2,4)
pcolor(lon, lat, DEV); shading flat; colorbar
title('Desviacion Estandar/Standard Deviation')


%% -----------------------------------------------------------------------
% Estimacion de los Espectros MTM-SVD
% Estimation of the MTM-SVD spectrum
%-------------------------------------------------------------------------

clear all

% Directorio con rutinas mtm-svd y datos a utilizar
% Directory with MTM-SVD routines and data
path1 = '/Users/sck117/Box Sync/MM_NSF2018/MTM-SVD/';
path(path,[path1,'Code']);  % incorpora el codigo al path/path to code
path4 = [path1,'Output/'];  % camino para archivos de datos/path to data

% Lectura de datos/Reading the data
%-----------------
% lectura de archivos de datos de nivel del mar formato 2D [tiempo-espacio]
% nivel del mar(sla), fecha(fec), Latitud(lat), Longitud(lon), indice(id)
% Reading the 2D SLA file [time-space], sea level anomaly (sla), date (fec),
% latitude (lat), longitude (lon), index (id)
load ([path4,'SLA_2D.mat']);

% paremetros para el calculo de los espectros/parameters for the spectrum calculation
%--------------------------------------------

% Padding con ceros/padding with zeros
% Permite definir la resolucion en la frecuencia. Se utiliza una o varias
% potencias de dos mayores a la longitud temporal de la serie
% Define the frequency resolution. One or several is used
% powers of two greater than the temporary length of the series.

lst = length(fec);            % longitud de la serie/length of the series
npad = 2^(nextpow2(lst)+2);   % segunda potencia mayor a lst/second power greater than lst

% Determina la frecuencia de muestreo/Determine the sampling frequency
anos_st = (fec(end)-fec(1))/365.25; % longitud de la serie en años /series length in years
T = lst/anos_st; % periodo de muestreo (en años)/samping period (in years)
fca = anos_st/lst; % frecuencia muestreo (ciclos por año)/sampling frequency (in years)

% Determina las frecuencias de la estimacion espectral/Determine the spectral estimation frequency
% nf = npad/2; % periodo fundamental para la estimacion espectral/fundamental period for the % spectral estimation
% ddf = 1./(npad*fca); % frecuencia de muestreo en npad intervalos/Frequency in npad intervals
% fr = ([1:nf]-1')*ddf; % vector de frecuencias/frequency vector

% Parametros para el espectro MTM/Parameters for the MTM spectrum
% Para minimizar la perdida de potencia espectral se recomienda que
% el numero de ventanas ortogonales sea igual a 2*bw-1
% To minimize the loss of spectral power it is recommended that
% the number of orthogonal windows equals 2 * bw-1
bw = 2; % ancho de banda bw=2 es el recomendado en la literatura
        % bandwidth=2 recommended in literature
k = 3;  % Numero de ventanas ortogonales/number of orthogonal windows

% Ponderacion por el area representada por cada serie (propocional a la
% latitud de cada posicion). Esto es opcional, y se utiliza para compensar
% la diferencia de los tamaños de los pixeles al considerar grandes areas,
% o cuando se quiere dar mas peso a una serie que otra (e.g. solo una serie
% en una region y multiples series en otras).
% Weighting by the area represented by each series (proportional to the
% latitude of each point). This is optional, and is used to compensate
% the difference in pixel sizes when considering large areas,
% or when you want to give more weight to one series than another (e.g. only one series
% in one region and multiple series in others).
        
[ x,y] = meshgrid(lon,lat);
w = cosd(y(id));

% Espectro nivel del mar/Spectrum of Sea Level anomaly
%---------------------------------------

% calculo del espectro de varianza fraccional local (LFV)
% calculation of the local fractional variance (LFV) spectrum
[LFVH] = mtm_svd_lfv(H,bw,k,fca,w,npad); 
LFVH.name = 'Espectro del nivel del mar/Sea Level Spectrum'; % cambia el nombre a la variable/%Change the variable name
LFVH % Despliega el nombre de las variables/Display the variable names

% niveles de significancia/significance levels
ql=[.99 .95 .90 .80 .50]; % niveles a calcular/levels to calculate
iter=100; % numero de iteraciones. Aconsejable minimo 1000/Number of iterations. minimum 1000
[fn, SLn]=mtm_svd_conf(H,bw,k,fca,iter,12,ql,w,0);

% Incorpora los resultados del limite de conmfianza en el archivo
% estructurado LFV, de tal forma de dejar organizado los resultados en
% una sola matriz
% Incorporates the results of the confidence limit in the
% structure LFV, organizes the results in a single matrix
LFVH.time=fec;
LFVH.quantil=ql;
LFVH.conflevel=SLn;
LFVH.confleveldomain=fn;

% grafico del espectro LFV/plot the LFV spectrum
%---------------------------
figure, %subplot(3,1,1)
semilogx(LFVH.specdomain,LFVH.spectrum, 'b')
xl = ([1/anos_st 12]); %frecuencia minima y maxima
xt = 1./[ 18 15 10 7 5 3 2 1 .5 .3 .25 .17 .08 ]; % frecuencias para labels
set(gca, 'xlim',(xl),'xtick',xt,'xticklabel', 1./xt,'ylim',([0.3 1]))
title('Espectro SLA/SLA Spectrum'), ylabel('LFV'), xlabel('Periodo [años]/Period [years]')
grid on, hold on
plot(LFVH.confleveldomain, LFVH.conflevel, 'r'); % niveles de significancia
text([xl(2) xl(2) xl(2) xl(2) xl(2)], LFVH.conflevel(:,end), {'99%'; '95%'; '90%'; '80%'; '50%'}) 

% Guarda las matrices de los espectros LFV del nivel del mar
save ([path4, 'espectros_LFVH'],'LFVH')

%% -----------------------------------------------------------------
% Espectro de varianza LFV Esfuerzo del Viento/Spectrum of the wind LFV
%-----------------------------------------------------------------

% lectura de datos del esfierzo del viento formato 2D [tiempo-espacio]
% Read the 2D wind data [time-space]
% esfuerzo del viento(EU,EV), fecha(fec), Latitud(lat), Longitud(lon), indice(id)
% wind componenets (EU,EV), date(fec), latitude(lat), longitude(lon), index (id)
        
load ([path4,'WIN_2D.mat']);

% Calculaola matriz 2D de la magnitud del esfuerzo del viento
% calculate the 2D wind speed matrix
W = sqrt(EU.^2 + EV.^2);

% calcula el espectro LFV/Calculate the LFV spectrum
[LFVW] = mtm_svd_lfv(W,bw,k,fca,w,npad); 
LFVW.name='Espectro del esfuerzo del viento/Wind Spectrum';% cambia el nombre de la variable
            %change the variable name
LFVW % Despliega el nombre de las variables/Display the variable names

% niveles de significancia/significance levels
ql=[.99 .95 .90 .80 .50]; % niveles a calcular/levels to Calculate
        iter=100; % numero de iteraciones. Aconsejable 1000/number of iterations.  minimum 1000
[fn, SLw]=mtm_svd_conf(W,bw,k,fca,iter,12,0,w,0);

% organizando en una sola matriz/organize in one matrix structure
LFVW.time=fec;
LFVW.quantil=ql;
LFVW.conflevel=SLw;
LFVW.confleveldomain=fn;

% Despliega el grafico del espectro LFV/display the spectrum plot
figure, % subplot(3,1,2)
semilogx(LFVW.specdomain,LFVW.spectrum, 'k')
xl = ([1/anos_st 12]); %frecuencia minima y maxima
xt = 1./[ 18 15 10 7 5 3 2 1 .5 .3 .25 .17 .08 ]; % frecuencias para labels
set(gca, 'xlim',(xl),'xtick',xt,'xticklabel', 1./xt,'ylim',([0.4 1]))
title('Espectro Esfuerzo del Viento/wind spectrum'), ylabel('LFV'), xlabel('Periodo [años]/Period [years]')
grid on, hold on
plot(LFVW.confleveldomain, LFVW.conflevel, 'r')% nivel de significancia
text([xl(2) xl(2) xl(2) xl(2) xl(2)], LFVW.conflevel(:,end), {'99%'; '95%'; '90%'; '80%'; '50%'}) 

% Guarda las matrices de espectro LFV del viento
save ([path4, 'espectros_LFVW'],'LFVW')

%-----------------------------------------------------------------
% Espectro Conjunto nivel del mar - Esfuerzo del Viento
% combined SLA-wind spectrum
%-----------------------------------------------------------------

% agrupa el campo del nivel del mar con las componentes del viento en
% una sola matriz
% group the SLA and wind fields in a single matrix
HW = [H W];

% calcula el espectro LFV comjunto/calculate the combined LFV spectrum
[LFVHW] = mtm_svd_lfv(HW,bw,k,fca,[w;w],npad); 
LFVHW.name='Espectro conjunto nivel del mar - esfuerzo del viento';
LFVHW

% niveles de significancia/significance levels
ql=[.99 .95 .90 .80 .50]; % niveles a calcular/levels to calculate
iter=100; % numero de iteraciones. Aconsejable 1000/number of iterations.  minimum 1000.
[fn, SLn]=mtm_svd_conf(HW,bw,k,fca,iter,12,0,[w;w],0);
LFVHW.time=fec;
LFVHW.quantil=ql;
LFVHW.conflevel=SLn;
LFVHW.confleveldomain=fn;

% grafico del espectro LFV conjunto/plot the combined LFV spectrum
figure, %subplot(3,1,3)
semilogx(LFVHW.specdomain,LFVHW.spectrum, 'k')
set(gca, 'xlim',(xl),'xtick',xt,'xticklabel', 1./xt,'ylim',([0.3 1]))
title('Espectro Conjunto SLA - Esfuerzo del Viento/Combined spectrum'),ylabel('LFV'),xlabel('Periodo [años]/Period [years]')
grid on, hold on
plot(LFVHW.confleveldomain, LFVHW.conflevel, 'r') % nivel de significancia/significance level
text([xl(2) xl(2) xl(2) xl(2) xl(2)],LFVHW.conflevel(:,end), {'99%'; '95%'; '90%'; '80%'; '50%'}) 

% Guarda las matrices de los espectros LFV conjunto/save the combined LFV spectrum matrix
save ([path4, 'espectros_LFVHW'],'LFVHW')

%% -------------------------------------------------------------------------
% Recontruccion del patron espacial para la frecuencia anual
% del nivel del mar
% Reconstruction of the spatial pattern for the annual frequency
% of sea level
%--------------------------------------------------------------------------

clear all

% Directorio con rutinas mtm-svd y datos a utilizar
% Directory with the MTM-SVD code and data
path1 = '/Users/sck117/Box Sync/MM_NSF2018/MTM-SVD/';
path(path,[path1,'Code']);  % incorpora el codigo al path/path to code
path4 = [path1,'Output/'];  % camino para archivos de datos/path to data

% Lectura de datos/Reading the data
%-----------------
% lectura de archivos de datos de nivel del mar formato 2D [tiempo-espacio]
% nivel del mar(sla), fecha(fec), Latitud(lat), Longitud(lon), indice(id)
% reading the 2D SLA files [time-space]
% Sea level anomaly (sla), date (fec), latitude (lat), longitude(lon), index (id)
  load ([path4,'SLA_2D.mat']);
  
% paremetros para las reconstrucciones/reconstruction parameters
%--------------------------------------------

% Padding con ceros/padding with zeros
% Permite definir la resolucion en la frecuencia. Se utiliza una o varias
% potencias de dos mayores a la longitud temporal de la serie
% Define the resolution in the frequency.
% One or several powers of two greater than the temporary length of the series
lst = length(fec);            % longitud de la serie/length of the series
npad = 2^(nextpow2(lst)+2);   % segunda potencia mayor a lst/power of 2 of length

% Determina la frecuencia de muestreo/determine the sampling frequency
anos_st = (fec(end)-fec(1))/365.25; % longitud de la serie en años/length in years
T = lst/anos_st; % periodo de muestreo (en años)/sampling period (years)
fca = anos_st/lst; % frecuencia muestreo (ciclos por año)/sampling frequency (years)

% Parametros para el espectro MTM/MTM spectrum parameters
% Para minimizar la perdida de potencia espectral se recomienda que
% el numero de ventanas ortogonales sea igual a 2*bw-1
% To minimize the loss of spectral power it is recommended that
% the number of orthogonal windows equals 2 * bw-1

bw = 2; % ancho de banda bw=2 es el recomendado en la literatura
        %bandwidth bw=2 is recommended in literature
k = 3;  % Numero de ventanas ortogonales/number of orthogonal windows

% Ponderacion por el area representada por cada serie (propocional a la
% latitud de cada posicion). Esto es opcional, y se utiliza para compensar
% la diferencia de los tamaños de los pixeles al considerar grandes areas,
% o cuando se quiere dar mas peso a una serie que otra (e.g. solo una serie
% en una region y multiples series en otras).  
% Weighting by the area represented by each series (proportional to the
% latitude of each point). This is optional, and is used to compensate
% the difference in pixel sizes when considering large areas,
% or when you want to give more weight to one series than another (e.g. only one series
% in one region and multiple series in others).

[ x,y] = meshgrid(lon,lat); % Campos de coordenadas/coordinate fields
w = cosd(y(id));

% frecuencias a reconstruir/reconstruction frequencies
% fo puede ser un vector con una o varias frecuencias que se deseen
% reconstruir
% fo is a vector with one or more reconstruction frequencies
        
fo = [1]; %(fo == 1 frecuancia anual)/(fo == 1 annual frequency)

% coordenadas de las series de tiempo en la matriz 2D
% time series coordinates in the 2D matrix
x2d = x(id); % latitud de c/u de las series/latitudes
y2d = y(id); % longitud de c/u de las series/longitudes

% reconstruccion del ciclo anual/reconstruction of the annual cycle
%-----------------------------------

% R  serie reconstruida para la(s) frecuencia(s) fo
% RC ciclo canonico de la(s) freciencia(s) fo
% RP matriz de estructura que contiene varias variables (ver help)
% U,S,V matrices de la descomposicion ortogonal
% R series reconstructed for the frequency(ies) fo
% RC canonical cycle of the frequency(ies) fo
% RP structure matrix that contains several variables (see help)
% U, S, V matrices of orthogonal decomposition

[R,RC,RP,U,S,V]=mtm_svd_bandrecon(H,bw,k,fca,fo,1,npad,0,w);

% Varianza total explicada por la(s) frecuencia(s) fo en la region (%)
% Total explained variance for the frequency(ies) fo in the region (%)
sprintf('Porcentaje de varianza explicada =  %4.2f',RP.totvarexplained)
sprintf('Percentage of variance explained =  %4.2f',RP.totvarexplained)


% Figuras asociadas a la recontruccion/figures for the reconstruction
%-------------------------------------

% Figura de la varianza explicada para la(s) frecuencia(s) fo en cada una
% de las series de toda la region de estudio
% Figure of the explained variance for the frequency(ies) fo in each of the series
% of the entire study region
%-----------------------------------------------------------------
RV=x*nan; RV(id)=RP.varexplained; % campo de varainza explicada/explained variance
figure, colormap jet % escala de colores lineal (empieza en azul y termina en rojo)
        %lineal color scale (blue to red)
pcolor(x, y, RV); shading flat; 
hold on; [c,h]=contour(x, y, RV,  [10:10:100], 'k')
clabel(c,h),title(['Varainza explicada - frecuencia' num2str(fo) ' año'])

% Figura de la Fase (en grados), asociada a la(s) frecuencia(s) fo en cada
% una de las series de toda la region de estudio
%  La fase es relativa a la primera posicion de la matriz tiempo-espacio
% Figure of the Phase (in degrees), associated to the frequency(ies) fo in each
% of the series of the entire study region
% The phase is relative to the first position of the time-space matrix
%-------------------------------------------------------------------------
RF=x*nan; RF(id)=RP.phase;   % campo de varainza explicada/explained variance
figure;
colormap hsv;
% escala de colores circular (empieza en rojo y termina en rojo) para
% representar mejor el cambio de fase entre 0 y 360
% circular color scale (red to red) for representing the entire phase from
% 0 to 360
pcolor(x, y, RF); shading flat; 
hold on; [c,h]=contour(x, y, RF,  [30:30:360], 'k');
clabel(c,h); title(['Fase [°]'])

% Figura de la diferencia de fase (en dias) respecto a la primera
% posicion de la matriz 2D tiempo-espacio
% Figure of the phase difference (in days) with respect to the first
% point of the 2D time-space matrix
%-----------------------------------------------------------------
RL = x*nan; RL(id)=RP.lag;   % campo de la diferencia de fase (años)/phase diff (years)
RL = RL.*365.25;             % campo de la diferencia de fase (dias)/phase diff (days)
figure; colormap jet; % escala de colores lineal (de azul a rojo)
        $lineal colormap (blue to red)
pcolor(x, y, RL); shading flat; 
hold on; [c,h]=contour(x, y, RL,  [0:30:365], 'k');
clabel(c,h); title(['Desfase [dias]'])

% Reconstruccion de la serie de tiempo en una posicion determinada para
% la(s) frecuencia(s) fo
% Reconstruction of the time series in a certain location for the frequency(ies) fo
%----------------------------------------------------------------------

lons=-76; lats=-34; % Posisciond e la serie de tiempo/position of the time series
ist=find(x2d==lons & y2d==lats);% indice de la serie en la matriz 2D
        %index of the series in the 2D matrix
st=H(:,ist); % obtiene la serie de tiempo original /get the original time series
sr=R(ist, :)+RP.temporalmean(ist); % serie reconstruida mas su valor medio
        % reconstructed series plus its mean value

% Figura con serie original y serie reconstruida con la(s) frecuancia(s) fo
%Figure with the original and reconstructed series for the frequency(ies) fo
figure;
plot(fec, st, 'k'); datetick('x')
title(['SLA (latitud=' num2str(lats) ' longitud=' num2str(lons) ')'] )
ylabel('[cm]')
hold on; 
plot(fec, sr, 'r');

% Reconstruccion del patron espacial asociado a la(s) frecuencia(s) fo
% para diferentes fases (en grado) del ciclo
%reconstruction of the spatial pattern associated with frequency(ies) fo
% for different phases of the cycle (in degrees)
%----------------------------------------------------------------------

pha = [0:45:360]; % fase (en grados) del ciclo que se desean graficar
        %cycle phase (in degrees) to plot
figure;
for i=1:8 % ejemplo para ocho fases diferentes en una pagina
    %example for 8 different phases in one page
    tic
    subplot(2,4,i);
    % Reconstruccion del campo espacial en la fase pha(i)+1 grados
    % Reconstruction of the spatial field in phase pha(i)+1 degrees
    A=x*nan; A(id) = RC(:, pha(i)+1);
    pcolor(x, y, (A)); shading flat;
    caxis([-5 5]);
    hold on;
    [c,h] = contour(x, y, A, [ -5:1:5], 'k');
    clabel(c,h);
    title(['Fase ' num2str(pha(i))]);
    toc
end




