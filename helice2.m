%% Projeto 2 - Curvas de desempenho de hélice
% Propulsão 

% Marina Matsudo Miyasato
% Talita Vitoria De Freitas Martin
% Gabriel Guimarães Pereira 
 
clear variables
close all 

%% Variáveis de entrada 

% Aerofólio: Clark-Y (ao longo de toda a pá)

% Parâmetros atmosféricos
H = 0;

[~, Vsom, ~, rho, NU] = atmosisa(H,"extended",true); %tool box = atmosfera padrão (o matlab tem implementado)

% T - Temperatura (K)               % rho - Densidade (kg/m^3)
% Vsom - Velocidade do som (m/s)    % NU - Viscosidade cinemática (m^2/s)
% P - Pressão (Pa)                  % H - Altitude (m)

D_tip = 1.93;       % Diâmetro da hélice - metros 
D_rub = 0.30;       % Diâmetro do cubo (hub) - metros 
R_tip = D_tip/2;    % Raio da hélice - metros 
R_hub = D_rub/2;    % Raio do cubo (hub) - metros

% Determinação da posição radial
N_elementos = 50;   % Divisão ao longo do aerofólio
dr = (R_tip - R_hub)/N_elementos;     % largura do elemento da pá
for i = 1:N_elementos
    r(i) = R_hub + (i - 1/2).*dr;     % Posição radial - distância arbritrária de dentro do cubo ao elemento de pá
end

C = 0.1524 - 0.0508 .* (r./R_tip);    % Corda (em metros) em função da posição radial r; R é o raio da hélice

beta_tip2 = 15;     % Ângulo de passo na ponta – hélice 2
beta_hub2 = 60;     % Ângulo de passo no cubo – hélice 2

RPM = 2500;         % Rotação por minuto
n = RPM/60;         % Rotação da hélice em RPS
Npas = 3;           % Número de pás

%% Funções

beta_hub = beta_hub2;
beta_tip = beta_tip2;
beta2 = calcula_beta(beta_hub, beta_tip, R_hub, R_tip,r); %para a segunda hélice

% Passo (pitch distance) 
P_D = tan(beta2)*2*pi*(0.75*R_tip);

% Blade Element Theory
w = 0;              % velocidade induzida no elemento
alpha_i = 0;        % ângulo de ataque induzido devido à velocidade induzida

J = 0:0.01:1.5;        % Razão de Avanço

% Cálculo de V0 a partir da razão de avanço
for i = 1:length(J) 
V0 = (J(i)*n*D_tip);                    % Velocidade de voo da aeronave na direção de eixo da hélice

Omega = 2*pi*(n);                       % velocidade angular da hélice (rad/s)

phi = (atand(V0./(Omega.*r)));          % ângulo de espiral (helix angle)

V_E = sqrt((w+ V0).^2 + (Omega.*r).^2); % Velocidade resultante efetiva

Re = rho*(V_E.*C)./NU;                  % Número de Reynolds

Ma = V_E./Vsom;                         % Número de Mach

 %% Equações aerodinâmicas 

dCl_dalpha = 0.00078.*(Re./10^6) + 0.0704;
Cl_max = 0.0645.*(Re./10^6) + 1.1699;

% Ângulo de ataque de sustentação nula
alpha_zl = -1.358*10^(-4).*(Re./10^6).^5 + 3.180*10^(-3).*(Re./10^6).^4 - 0.0242.*(Re./10^6).^3 + 0.0366.*(Re./10^6).^2 + 0.3421.*(Re./10^6) - 6.586;

C_D0 = 1.129*10^(-6).*(Re./10^6).^4 - 2.682*10^(-5).*(Re./10^6).^3 + 2.299*10^(-4).*(Re./10^6).^2 - 9.154*10^(-4).*(Re./10^6) + 0.01034;
C_f = 1.7*10^(-7).*(Re./10^6).^5 + 4.73*10^(-6).*(Re./10^6).^4 - 5.17*10^(-5).*(Re./10^6).^3 + 2.85*10^(-4).*(Re./10^6).^2 - 8.43*10^(-4).*(Re./10^6) + 0.00539;

% Correções de compressibilidade 

Cl_alpha_compressivel = (dCl_dalpha)./(sqrt(1 - Ma.^2));
Cf_compressivel = C_f.*(0.000162.*(Ma.^5) - 0.00383.*(Ma.^4) + 0.0332.*(Ma.^3) - 0.118.*(Ma.^2) + 0.0204.*Ma + 0.996);

Cd_forma = C_D0 - C_f;

C_D0_corrigido = Cf_compressivel + Cd_forma;

alpha = beta2 - phi + alpha_zl;         % ângulo de ataque no elemento

%%  Forças aerodinâmicas no elemento

Cl = Cl_alpha_compressivel.*(beta2 - phi + alpha_zl - alpha_i); % coeficiente de sustentação do perfil

Cd = C_D0_corrigido + 0.06262.*(Cl.^2); % coeficiente de arrasto do perfil

dL = 1/2.*rho.*(V_E.^2).*C.*Cl.*dr;     % Força de sustentação no elemento

dD = 1/2.*rho.*(V_E.^2).*C.*Cd.*dr;     % Força de arrasto no elemento

%% Correções de ponta e de cubo de Prandtl 

% Parâmetro de correção de ponta
P_tip = (Npas.*(R_tip-r))./(2.*r.*sind(phi));

% Fator de correção de ponta 
F_tip = (2/pi)*acos(exp(-P_tip));

% Parâmetro de correção de cubo
P_hub = (Npas.*(r - R_hub))./(2.*r.*sind(phi));

% Fator de correção de cubo 
F_hub = (2/pi)*acos(exp(-P_hub));

% Fator de correção de Prandtl
F_P = F_hub.*F_tip;

%% Empuxo, torque e potência total da hélice
% Considerando as conrreções de Prandtl

% Empuxo 
T(i) = Npas * (sum(F_P.*dL.*cosd(phi + alpha_i))) - Npas * (sum(F_P.*dD.*sind(phi + alpha_i)));

% Torque
Q(i) = Npas * (sum(F_P.*r.*dL.*sind(phi + alpha_i))) + Npas * (sum(F_P.*r.*dD.*cosd(phi + alpha_i)));

% Potência
P(i) = Npas * (sum(F_P.*Omega.*r.*dL.*sind(phi + alpha_i))) + Npas * (sum(F_P.*Omega.*r.*dD.*cosd(phi + alpha_i)));

end 

%% Cálculo dos coeficientes adimensionais

% Coeficiente de empuxo 
C_T = T./(rho.*(n^2)*(D_tip^4));

% Coeficiente de torque 
C_Q = Q./(rho.*(n^2)*(D_tip^5));

% Coeficiente de potência
C_P = P./(rho.*(n^3)*(D_tip^5));

% Eficiência da hélice 
Eta = (C_T.*J)./C_P;

% Figura de mérito 
FOM = sqrt(2/pi).*((C_T.^1.5)./C_P);

% Fator de atividade 
AF = (100000/16)*(1/D_tip)*0.027921;

TAF = Npas .* AF;

% Solidez 
sigma = Npas.*(C./2*pi*0.7*R_tip);

%% Plots

% Coeficiente de empuxo x Razão de Avanço
figure(1);
plot(J, C_T,'LineWidth',1); 
xlabel('{V}/{n\cdotD}');
ylabel('C_{T}');
title('Razão de Avanço x Coeficiente de empuxo')
grid on;

% Coeficiente de empuxo x Razão de Avanço
figure(2);
plot(J, C_Q,'LineWidth',1); 
xlabel('{V}/{n\cdotD}');
ylabel('C_{Q}');
title('Razão de Avanço x Coeficiente de torque')
grid on;

% Coeficiente de potência x Razão de Avanço
figure(3);
plot(J, C_P,'LineWidth',1); 
xlabel('{V}/{n\cdotD}');
ylabel('C_{P}');
title('Razão de Avanço x Coeficiente de potência')
grid on;

% Eficiência x Razão de Avanço
figure(4);
plot(J, Eta,'LineWidth',1); 
xlabel('{V}/{n\cdotD}');
ylabel('\eta');
title('Razão de Avanço x Eficiência da hélice')
grid on;
