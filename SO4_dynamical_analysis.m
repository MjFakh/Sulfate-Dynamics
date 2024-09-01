%%   Dynamical analysis of SO4
 
clear all

%------------------------------------------
% -------------  SET PARAMETERS -----------
%------------------------------------------

n_size = 3000;

for i=1:n_size

SO4_1(:,1:10) = [0.1 1 2 3 4 5 6 7 8 9];
SO4_1(:,11:3010) = 10:15:45000;

SO4_modern = 28000; %ppm
Ca_modern  = 10000;

power_gyp_MC = [0.1 0.2 0.3 0.4 0.5];
power_gyp    = (power_gyp_MC(randi(numel(power_gyp_MC))));  % submerged area power coefficient 

S_gyp_modern_MC = [150 200 250];
S_gyp_modern = (S_gyp_modern_MC(randi(numel(S_gyp_modern_MC)))).*1E18;  % Modern crustal reservoir size of gypsum S mol

S_pyr_modern_MC = [150 200 250];
S_pyr_modern = (S_pyr_modern_MC(randi(numel(S_pyr_modern_MC)))).*1E18;  % Modern crustal reservoir size of pyrite S mol

S_total = S_gyp_modern + S_pyr_modern;

F_gypsum_MC   = [1500 1700 2000];    % Modern gypsum flux              
F_gypsum  = (F_gypsum_MC(randi(numel(F_gypsum_MC)))).* 1E4; 

F_wpyrite_bio_MC   = [200 250 300];    % Modern weathering pyrite biotic flux              
F_wpyrite_bio  = (F_wpyrite_bio_MC(randi(numel(F_wpyrite_bio_MC)))).* 1E4; 

F_wpyrite_abio_MC   = [300 350 400];    % Modern weathering pyrite abiotic flux              
F_wpyrite_abio  = (F_wpyrite_abio_MC(randi(numel(F_wpyrite_abio_MC)))).* 1E4;

Fvolc_MC   = [500 600 800 1000];     % Total volcanic flux           
Fvolc  = (Fvolc_MC(randi(numel(Fvolc_MC)))).* 1E4; 

Vmax1_MC   = [0.05 0.07 0.1 0.2 0.3 0.4 0.5];     % Max SR in WC 1E-8 mmol/cm2/yr
Vmax1  = (Vmax1_MC(randi(numel(Vmax1_MC)))).*1E-3.*1E4.*1E-9.*1E4; 

AREAtot = 3.6*1E14;

km_MC   = [5 10 20 50 70 100 200 500];     % half saturation for sulfate      
km  = (km_MC(randi(numel(km_MC)))); 

km_Fe_MC   = [50 70 100];     % half saturation for iron         
km_Fe  = (km_Fe_MC(randi(numel(km_Fe_MC)))); 

Ca_MC   = [10 15 20 25 30];     % Calcium concentration          
Ca_1  = (Ca_MC(randi(numel(Ca_MC)))).*1000; 

% ----------------------------------------------------------
% ------------ parameters for different runs ---------------
% ----------------------------------------------------------

% --------------------- Modern -----------------------------

% gypsum weathering

f_eros_modern_MC   = [0.5 0.7 1];     % erosion efficiency         
f_eros_modern      = (f_eros_modern_MC(randi(numel(f_eros_modern_MC))));

F_gyp_weath_modern_MC   = [1500 1700 2000];     % Modern gypsum weathering               
F_gyp_weath_modern      = (F_gyp_weath_modern_MC(randi(numel(F_gyp_weath_modern_MC)))).* 1E4; 

F_gyp_weath_1   = f_eros_modern.*F_gyp_weath_modern;

% pyrite weathering

f_O2_abio_pyrite1  = 1;

F_pyr_weath_1   = f_eros_modern.*f_O2_abio_pyrite1.*F_wpyrite_abio + f_eros_modern.*F_wpyrite_bio;

% total weathering

F_weath_1 = F_gyp_weath_1 + F_pyr_weath_1;

% pyrite burial 

Fe_modern_MC   = [1 5 10];     % Fe concentration        
Fe_modern      = (Fe_modern_MC(randi(numel(Fe_modern_MC)))); 

f_area_anoxia_modern_MC   = [0.0001 0.0005 0.001 0.003 0.005];     % areal fraction of anoxia         
f_area_anoxia_modern      = (f_area_anoxia_modern_MC(randi(numel(f_area_anoxia_modern_MC)))); 

fprod_modern    = 1;            % fraction of modern primary production for SR

% gypsum burial

Area_submerg_modern = 5; 
Area_submerg1_MC    = [5 7 10 12 15 17 20 23 25 27 30 32 35 37 40];     % percent submerged continental area      
Area_submerg1       = (Area_submerg1_MC(randi(numel(Area_submerg1_MC)))); 

F_gypsum_modern = F_gypsum.*(SO4_1./SO4_modern).*(Ca_1./Ca_modern).*((Area_submerg1./Area_submerg_modern).^power_gyp);

 
% --------------------- Precambrian ------------------------

% gypsum weathering

S_gyp_precam_MC = [25 50 75];
S_gyp_precam = (S_gyp_precam_MC(randi(numel(S_gyp_precam_MC)))).*1E18;  % Modern crustal reservoir size of gypsum S mol

f_eros_pre_MC   = [0.1 0.3 0.5 0.7 1];     % erosion efficiency         
f_eros_pre  = (f_eros_pre_MC(randi(numel(f_eros_pre_MC)))); 

F_gyp_weath_modern_MC   = [1500 1700 2000];     % Modern pyrite weathering               
F_gyp_weath_modern      = (F_gyp_weath_modern_MC(randi(numel(F_gyp_weath_modern_MC)))).* 1E4; 

F_gyp_weath_2   = f_eros_pre.*F_gyp_weath_modern.*(S_gyp_precam./S_gyp_modern);

% pyrite weathering

S_pyr_precam = S_total - S_gyp_precam;  % Modern crustal reservoir size of pyrite S mol

f_O2_abio_pyrite2_MC = [0.5 0.7 1];    % Modern weathering pyrite abio flux              
f_O2_abio_pyrite2  = (f_O2_abio_pyrite2_MC(randi(numel(f_O2_abio_pyrite2_MC))));

F_pyr_weath_2   = f_eros_pre.*f_O2_abio_pyrite2.*F_wpyrite_abio.*(S_pyr_precam./S_pyr_modern) + ...
                  f_eros_pre.*F_wpyrite_bio.*(S_pyr_precam./S_pyr_modern);


% total weathering

F_weath_2 = F_gyp_weath_2 + F_pyr_weath_2;


% pyrite burial 

Fe_precambrian_MC   = [10 50 100 200];     % Fe concentration        
Fe_precambrian      = (Fe_precambrian_MC(randi(numel(Fe_precambrian_MC)))); 

fprod_precambrian_MC   = [0.01 0.03 0.05 0.07 0.1];     % NPP factor       
fprod_precambrian  = (fprod_precambrian_MC(randi(numel(fprod_precambrian_MC)))); 

f_area_anoxia_precambrian_MC   = [0.1 0.2 0.3 0.4 0.5 0.7 0.8];     % areal fraction of anoxia         
f_area_anoxia_precambrian  = (f_area_anoxia_precambrian_MC(randi(numel(f_area_anoxia_precambrian_MC)))); 

% gypsum burial

Area_submerg2_MC   = [10 15 20 25 30 35 40 50];     % percent submerged continental area      
Area_submerg2      = (Area_submerg2_MC(randi(numel(Area_submerg2_MC)))); 

F_gypsum_precambrian = F_gypsum.*(SO4_1./SO4_modern).*(Ca_1./Ca_modern).*((Area_submerg2./Area_submerg_modern).^power_gyp);

% --------------------- OAE --------------------------------

% gypsum weathering

f_eros_OAE_MC   = [0.5 0.7 1 1.3 1.5];     % erosion efficiency         
f_eros_OAE      = (f_eros_OAE_MC(randi(numel(f_eros_OAE_MC))));

F_gyp_weath_3   = f_eros_OAE.*F_gyp_weath_modern;

% pyrite weathering

f_O2_abio_pyrite1  = 1;

F_pyr_weath_3   = f_eros_OAE.*f_O2_abio_pyrite1.*F_wpyrite_abio + f_eros_OAE.*F_wpyrite_bio;

% total weathering

F_weath_3 = F_gyp_weath_3 + F_pyr_weath_3;

% pyrite burial

Fe_OAE_1_MC   = [10 50 100];     % Fe concentration        
Fe_OAE_1  = (Fe_precambrian_MC(randi(numel(Fe_precambrian_MC)))); 

f_area_anoxia_OAE_1_MC   = [0.01 0.02 0.05 0.1];     % areal fraction of anoxia         
f_area_anoxia_OAE_1  = (f_area_anoxia_OAE_1_MC(randi(numel(f_area_anoxia_OAE_1_MC)))); 

fprod_OAE_1_MC   = [0.5 1 1.5 2];     % NPP factor       
fprod_OAE_1  = (fprod_OAE_1_MC(randi(numel(fprod_OAE_1_MC)))); 


% ----------------------------------------------------------
% --------- Sulfate reduction RESPIRATION RATES ------------
% ----------------------------------------------------------


% --------------------- Modern -----------------------------

F_SRR_1   = AREAtot.* fprod_modern.* f_area_anoxia_modern.*Vmax1.*(SO4_1./(SO4_1 + km)).*(Fe_modern./(Fe_modern + km_Fe));

% --------------------- Precambrian ------------------------

F_SRR_2   = AREAtot.* fprod_precambrian.* f_area_anoxia_precambrian.*Vmax1.*(SO4_1./(SO4_1 + km)).*...
           (Fe_precambrian./(Fe_precambrian + km_Fe));

% --------------------- OAE --------------------------------

F_SRR_3   = AREAtot.* fprod_OAE_1.* f_area_anoxia_OAE_1.*Vmax1.*(SO4_1./(SO4_1 + km)).*(Fe_OAE_1./(Fe_OAE_1 + km_Fe));

% ----------------------------------------------------------
% ----------------------- SO4 ODE --------------------------
% ----------------------------------------------------------

% --------------------- Modern -----------------------------

dSO4dt_1  = F_weath_1 + Fvolc - F_SRR_1 - F_gypsum_modern;

% --------------------- Precambrian ------------------------

dSO4dt_2  = F_weath_2 + Fvolc - F_SRR_2 - F_gypsum_precambrian;

% --------------------- OAE --------------------------------

dSO4dt_3  = F_weath_3 + Fvolc - F_SRR_3 - F_gypsum_modern;

% ------------------- Storing Data -----------------

dSO4dt_1_MC(i,:) = dSO4dt_1;
dSO4dt_2_MC(i,:) = dSO4dt_2;
dSO4dt_3_MC(i,:) = dSO4dt_3;

end

% -------------- PLOT RESULTS  -------------


hold on
linewidth = 2;
Fontsize = 12;

% Modern

subplot(3,3,2)
fanChart(SO4_1, dSO4dt_1_MC'.*1E-4.*1E-3,'median', 5:10:93, ...
    'alpha', .5, 'colormap', {'shadesOfColor', [0.04,0.30,0.48]});

hold on

plot(SO4_1,zeros(1,size(SO4_1,2)),'--k','lineWidth',linewidth);

set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
set(gca,'xscale','log')
ylabel ('d([SO4])/dt (Tmol/year)');
xlim([1 45000]);
ylim([-10 10]);

xticks([1 10 100 1000 10000])
xticklabels({'','','','','',})

box on
grid on
grid minor
ax.LineWidth = 2;


% Precambrian

subplot(3,3,5)

fanChart(SO4_1, dSO4dt_2_MC'.*1E-4.*1E-3,'median', 5:10:82, ...
    'alpha', .5, 'colormap', {'shadesOfColor', [0.04,0.30,0.48]});

hold on

plot(SO4_1,zeros(1,size(SO4_1,2)),'--k','lineWidth',linewidth);

set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
set(gca,'xscale','log')
ylabel ('d([SO4])/dt (Tmol/year)');
xlim([1 45000]);
ylim([-10 10]);

xticks([1 10 100 1000 10000])
xticklabels({'','','','','',})

box on
grid on
grid minor
ax.LineWidth = 2;


% OAE1

subplot(3,3,8)

fanChart(SO4_1, dSO4dt_3_MC'.*1E-4.*1E-3,'median', 5:10:82, ...
    'alpha', .5, 'colormap', {'shadesOfColor', [0.04,0.30,0.48]});

hold on

plot(SO4_1,zeros(1,size(SO4_1,2)),'--k','lineWidth',linewidth);

set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
set(gca,'xscale','log')
xlabel ('[SO_4] (\muM)' );
ylabel ('d([SO4])/dt (Tmol/year)');
xlim([1 45000]);
ylim([-10 10]);

xticks([1 10 100 1000 10000])
xticklabels({'1','10','100','1000','10000'})

box on
grid on
grid minor
ax.LineWidth = 2;
