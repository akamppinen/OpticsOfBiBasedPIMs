% Optical modeling of Bi based PIMs in solar cells %
% Aleksi Kamppinen
% Solar Energy Materials and Systems, University of Turku

% Clear workspace
clear all;

% Active material
% activeMat = 'CABI';
activeMat = 'CBI';

% Data storage folder
dateToday = datetime("today");
dataFolder = ['.\',num2str(year(dateToday)),'_',num2str(month(dateToday)),'_',num2str(day(dateToday)),'_',activeMat,'_Revision'];

[~,msg,msgID] = mkdir(dataFolder);
if ~isempty(msg)
    error('Check the proposed folder.')
end

%% Define settings 

% Natural constants
consts = ([]);
consts.('h') = 6.626*10^-34;    % Js
consts.('h_eV') = 4.136*10^-15; % eVs
consts.('c') = 2.998*10^8;      % m/s
consts.('kB') = 1.381*10^-23;   % J/K
consts.('q') = 1.602*10^-19;    % C
consts.('SB') = 5.670*10^-8;    % W/(m^2K^4);

% Visibility factor of the active layer that defines partial coherence
ActiveLayerVis = 0;

% Sensitivity analysis of extinction (/absorption) coefficient
kScale = 0.5:0.05:1.5;
% kScale = 1;

% Modelled wavelengths
dLambda = 1; % nm
lambda = transpose(300:dLambda:1300);
lambdaDataLims = transpose(360:dLambda:830); % Data limits for indoor spectrum

% Solar irradiation data
AM15 = load('AM15.mat').AM15;

% Irradiance for defined wavelength range
wl = transpose(300:dLambda:4000);
Bsolar = spline(AM15{:,1},AM15{:,2},wl);
Fsolar = (wl*1e-9)/(consts.h*consts.c).*Bsolar; % Photon flux (s^-1m^-2nm^-1)

% Indoor irradiance
filepathSpectra = '.\';
filenameSpectra = {'solar',...
    '2700K_1000lux_wLED.txt','4000K_1000lux_wLED_Spectroradiometer.txt',...
    '5000K_1010lux_wLED','6500K_1010lux_wLED'};

% Spectral luminous efficiency for photopic vision, V(lambda) 
% doi.org/10.25039/CIE.DS.dktna2s3
effPhotopicVision = table2array(readtable([filepathSpectra,'\','CIE_sle_photopic']));

% Illuminance 
desiredIlluminance = [100,500,1000]; % lux
refIlluminance = desiredIlluminance(end);

% Vary and scale spectra
IntensityFactors = cell(length(filenameSpectra),1);
B = cell(length(filenameSpectra),1);
Btot = cell(length(filenameSpectra),1);
for idxSpectrum = 1:length(filenameSpectra) 
    switch filenameSpectra{idxSpectrum}
        case 'solar'
            IntensityFactors{1} = 1;
            Btot{1} = sum(dLambda*Bsolar);
            B{1} = {Bsolar(ismembertol(wl,lambda,1e-5))}; % Cut solar spectrum to the modeled wl range
        otherwise
            dataIndoorSpectrum = readtable([filepathSpectra,'\',...
                filenameSpectra{idxSpectrum}]);
            B_indoor = zeros(length(lambda),1);
            switch filenameSpectra{idxSpectrum}
                case 'B4_LED_spectrum_at_1000lux'
                    Bdata = cell2mat(cellfun(@(x) str2double(x),...
                        table2cell(dataIndoorSpectrum),'UniformOutput',false));
                    wlsData = Bdata(:,1);                                       % Wavelength (nm)
                    B_indoor(ismembertol(lambda,wlsData,1e-5)) = Bdata(:,2);    % Radiance (Wm^-2nm^-1)
                otherwise
                    wlsData = cell2mat(cellfun(@(x) str2double(x),...
                        dataIndoorSpectrum{185:655,1},'UniformOutput',false));  % Wavelength (nm)
                    B_indoor(ismembertol(lambda,wlsData,1e-5)) = ...
                        dataIndoorSpectrum{185:655,2};                          % Radiance (Wm^-2nm^-1)
            end
            B{idxSpectrum} = cell(length(desiredIlluminance),1);
            Btot{idxSpectrum} = zeros(length(desiredIlluminance),1);
            IntensityFactors{idxSpectrum} = zeros(length(desiredIlluminance),1);
            for idxIlluminance = 1:length(desiredIlluminance)
                Illuminance = Illuminance([lambda,B_indoor],effPhotopicVision); % lux
                IntensityFactors{idxSpectrum}(idxIlluminance) = ...
                    desiredIlluminance(idxIlluminance)/Illuminance;
                B{idxSpectrum}{idxIlluminance} = B_indoor*...
                    IntensityFactors{idxSpectrum}(idxIlluminance);
                Btot{idxSpectrum}(idxIlluminance) = sum(dLambda*...
                    B{idxSpectrum}{idxIlluminance});
            end
    end
end

%%% Define settings
% Materials 
materials = {'SiO2','FTO','TiO2compact',activeMat,'SpiroOMeTAD','Au'};
indActiveLayer = find(strcmp(activeMat,materials));
% Semiconductor layers
semicondLayers = {'TiO2compact',activeMat,'SpiroOMeTAD'};
semicondLayerInds = [];
for idxMat = 1:length(semicondLayers)
    semicondLayerInds = [semicondLayerInds,...
        find(strcmp(semicondLayers(idxMat),materials))];
end
% Materials considered in optical model
lastOpticsLayer = 'Au';
indLastOpticsLayer = find(strcmp(lastOpticsLayer,materials));
% Coarse optics layers for sparse spatial discretization
coarseOpticsLayers = {'PVglass','EVA','SiO2'};

% Layer thickness variation: initial values and settings
tInitial = [2*10^6,400,50,400,50,100]; % (nm)
varyLayers = {activeMat};
tVarLayerInds = find(ismember(materials,varyLayers));
tStepSizes = 0.025*tInitial(tVarLayerInds);
NtStepsLow = 39;
NtStepsHigh = 40;

ActiveMATthicknesses = (tInitial(indActiveLayer)-NtStepsLow*tStepSizes):...
    tStepSizes:(tInitial(indActiveLayer)+NtStepsHigh*tStepSizes);

% Reference band gap
switch activeMat
    case 'CABI'
        EgReal = 1.98; % eV
        % EgActiveMATref = EgReal;
        EgActiveMATref = sort([1.9:0.05:2.1,EgReal]); % sensitivity analysis
    case 'CBI'
        EgReal = 2.33; % eV
        % EgActiveMATref = EgReal; 
        EgActiveMATref = sort([2.2:0.05:2.45,EgReal]); % sensitivity analysis
end

% Temperature cycling 
Tstc = 298.15; % K
Tlattice = Tstc; % 

%%% Material properties
% Optical data
opticalData = readtable('materialOpticalConstants_CABI_CBI_deviceModels.xlsx');

% Reference index of refraction data
wlActiveRef = rmmissing(opticalData.(strcat(activeMat,'_wl')));
nRealActiveRef = rmmissing(opticalData.(strcat(activeMat,'_n')));
nImagActiveRef = rmmissing(opticalData.(strcat(activeMat,'_k')));

%%% Set cell properties
[~,~,nReal,nImag,N_EDphases,W_EDphases,Nlayers] = SetUpCell(materials,...
    tInitial,indActiveLayer,ActiveLayerVis,indLastOpticsLayer,lambda,AM15,opticalData,consts);
clear('opticalData')

%% Data storage
% Error log
folderPathErrors = strcat(dataFolder,'\Errors');
mkdir(folderPathErrors);

% Varied k spectra (imaginary part of refractive index)
kSpectra = zeros(size(nImag,1),length(kScale));

% Photogenerations from modelled optics
StoreTotGphPred = zeros(length(desiredIlluminance),length(filenameSpectra),length(ActiveMATthicknesses),length(kScale),length(EgActiveMATref));
StoreGphPred = zeros(length(filenameSpectra),length(ActiveMATthicknesses),length(lambda),length(kScale));
StoreGphVsPosition = cell(length(filenameSpectra),length(ActiveMATthicknesses),length(kScale),length(EgActiveMATref));

% Normalized (device) optics
StoreReflection = zeros(length(ActiveMATthicknesses),length(kScale),length(lambda));
StoreAbsorption = zeros(length(ActiveMATthicknesses),length(kScale),length(lambda),length(materials));
StoreTransmission = zeros(length(ActiveMATthicknesses),length(kScale),length(lambda));

%% Simulation

% Loop over different absorption spectra
kInitial = nImag(:,indActiveLayer+1);
for idxkScale = 1:length(kScale)
    % Progress report
    disp(['kScale: ',num2str(kScale(idxkScale)),' (', ...
        num2str(idxkScale),'/',num2str(length(kScale)),')'])
    tStartkScale = tic;

    kSpectra(:,idxkScale) = kScale(idxkScale)*kInitial;
    nImag(:,indActiveLayer+1) = kSpectra(:,idxkScale);

    %%% Loop over active layer thicknesses
    idxActiveMATt = 1;
    errorCount = 0;
    while idxActiveMATt <= length(ActiveMATthicknesses)
        try
            % % Progress report
            % disp(['ActiveMAT thickness: ',num2str(ActiveMATthicknesses(idxActiveMATt)),' nm (', ...
            %     num2str(idxActiveMATt),'/',num2str(length(ActiveMATthicknesses)),')'])
            % tStartActiveMATt = tic;

            % Update thickness
            t = tInitial;
            t(indActiveLayer) = ActiveMATthicknesses(idxActiveMATt);

            % Spatial discretization
            N_opticalSublayers = round(t);
            N_opticalSublayers(ismember(materials,coarseOpticsLayers)) = 100;
            p = zeros(1,length(materials));
            [xPos,xMat,~,tCumsum] = xDiscretization(t,N_opticalSublayers,p);

            % Optics
            [E2,~,R,A,T,~] = CalcSolarCellOptics(...
                nReal(:,1:(indLastOpticsLayer+2)),nImag(:,1:(indLastOpticsLayer+2)),...
                t(1:indLastOpticsLayer),xPos(ismember(xMat,(1:indLastOpticsLayer))),...
                xMat(ismember(xMat,(1:indLastOpticsLayer))),lambda,...
                N_EDphases(1:indLastOpticsLayer),W_EDphases);

            StoreReflection(idxActiveMATt,idxkScale,:) = R;
            StoreAbsorption(idxActiveMATt,idxkScale,:,:) = A;
            StoreTransmission(idxActiveMATt,idxkScale,:) = T;

            % Loop over spectral ratio factor
            for idxSpectrum = 1:length(filenameSpectra)
                % Progress report
                % tStartSpectrum = tic;
                % disp(['Indoor spectrum ratio: ',num2str(SRfactor(idxSRfactor)),' (', ...
                %     num2str(idxSpectrum),'/',num2str(length(filenameSpectra)),')'])

                % Loop over intensity factors
                for idxIfactor = 1:length(IntensityFactors{idxSpectrum})

                    % Vary band gap for sensitivity analysis
                    for idxEg = 1:length(EgActiveMATref)
                        EgWl = consts.h_eV*consts.c*10^9/EgActiveMATref(idxEg); % nm

                        % Photogeneration
                        [Gx,Glx,~,~,~,~,~,~] = ...
                            absorption(E2,xPos(ismember(xMat,(1:indLastOpticsLayer))),...
                            xMat(ismember(xMat,(1:indLastOpticsLayer))),...
                            nReal(:,1:(indLastOpticsLayer+2)),nImag(:,1:(indLastOpticsLayer+2)),...
                            EgActiveMATref(idxEg),EgWl,indActiveLayer,...
                            B{idxSpectrum}{idxIfactor},lambda,...
                            indLastOpticsLayer,consts,Tlattice);

                        % Total photogeneration
                        StoreTotGphPred(idxIfactor,idxSpectrum,idxActiveMATt,idxkScale,idxEg) = ...
                            numInt(tCumsum(indActiveLayer)*1e-9,tCumsum(indActiveLayer+1)*1e-9,...
                            xPos(xMat==indActiveLayer)*1e-9,Gx);

                        if ismembertol(idxIfactor,length(IntensityFactors{idxSpectrum}),1e-5)
                            % Photogeneration profiles vs position (only
                            % for the last intensity value, photogeneration
                            % can be scaled to different intensity conditions
                            % afterwards if necessary)
                            StoreGphVsPosition{idxSpectrum,idxActiveMATt,idxkScale,idxEg} = ...
                                [transpose(xPos(xMat==indActiveLayer))-tCumsum(indActiveLayer),transpose(Gx)];
                            if ismembertol(EgActiveMATref(idxEg),EgReal,1e-4)
                                % Photogeneration vs wavelength for the
                                % cumulative Jph (only for the determined
                                % band gap)
                                idxLambda = 1;
                                while lambda(idxLambda)<=EgWl && lambda(idxLambda)<length(lambda)
                                    StoreGphPred(idxSpectrum,idxActiveMATt,idxLambda,idxkScale) = ...
                                        numInt(tCumsum(indActiveLayer)*1e-9,tCumsum(indActiveLayer+1)*1e-9,...
                                        xPos(xMat==indActiveLayer)*1e-9,Glx(idxLambda,:));
                                    idxLambda = idxLambda+1;
                                end
                            end
                        end
                    end
                end
            end
            % % Progress report
            % % tElapsed = toc(tStart);
            % tElapsedActiveMATt = toc(tStartActiveMATt);
            % % disp(['Runtime since beginning: ', num2str(tElapsed/(60)),' min']);
            % disp(['Runtime of last ActiveMATt: ', num2str(tElapsedActiveMATt),' s']);

            % Increase index
            idxActiveMATt = idxActiveMATt+1;
            errorCount = 0;

        catch ME
            errorCount = errorCount+1;
            % Write the error to a file
            fid = fopen([folderPathErrors,...
                '\errorlog_ActiveMATt',num2str(round(ActiveMATthicknesses(idxActiveMATt))),...
                'nm_errorCount',num2str(errorCount),'.txt'],'w');
            fprintf(fid, '%s', ME.getReport('extended', 'hyperlinks','off'));
            fclose(fid);

            pause(5*errorCount) % pause
            if errorCount>5
                pause
            end
        end
    end
    % Progress report
    % tElapsed = toc(tStart);
    tElapsedkScale = toc(tStartkScale);
    % disp(['Total runtime so far: ', num2str(tElapsed/(60)),' min']);
    disp(['Runtime of last kScale: ', num2str(tElapsedkScale/60),' min']);
end

% Save workspace
save([dataFolder,'\workspace.mat']);
