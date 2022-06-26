clc; clear; close all
%% Space Propulsion
%%% ASSUMPTIONS %%%
% air is a ideal gaz -> adiabatic
% variation of bottle's pressure in an adiabatic process
% water is incompressible
% no losses at the throat (-> Bernouilli)
% no thrust given by the air at the end of burn
% negligeable mass of air 
% the rocket is going straight following the z axis 
%% Constants
% mass
number_booster=2; %[-]
mass_parachute=0.03;
dry_mass_booster=0.075; % [kg]
dry_mass_core=0.2+mass_parachute;
mass_water_booster=0.5;
mass_water_core=0.7;
dry_mass_total=number_booster*dry_mass_booster+dry_mass_core+mass_parachute; % with 2 boosters
wet_mass_total=dry_mass_total+number_booster*mass_water_booster+mass_water_core; % with 2 boosters

% volume
rho_water=1000; %[kg/m^3]
booster_volume=0.0015; %[m^3]
core_volume=0.002;
vol_air_booster_init=booster_volume-mass_water_booster/rho_water;
vol_air_core_init=core_volume-mass_water_core/rho_water;

% proportion water/air initial
Prop_w_a_core=mass_water_core/rho_water/vol_air_core_init;
Prop_w_a_booster=mass_water_booster/rho_water/vol_air_booster_init;

% pressure
init_press_booster=7*10^5; % [Pa]
init_press_core=7*10^5; 
P_atm=10^5;

% throat area 
booster_throat=(0.005)^2*pi; %[m^2]
core_throat=(0.005)^2*pi;

% drag force and gravity constants 
A_1=pi*(0.05)^2+number_booster*pi*(0.05)^2; %[m^2] (core 10cm of diameter, booster 10cm diam)
A_2=pi*(0.05)^2; %[m^2] (core 8cm of large)
Cd=0.4; %(inbetween drop and rectangle)
rho_air=1.23; %[kg/m^3] @ 20 degress
g=9.81; %[m/s^2]

% parachute descent
S_eff=pi*0.55^2/4;
Cd_para=1.5;

% others
gamma=1.4; % coeff adiabatic air
lambda=0.9; % coeff efficiency nozzle

%% Calculs
% initialisation 
V_b(1)=vol_air_booster_init;
P_b(1)=init_press_booster;
V_c(1)=vol_air_core_init;
P_c(1)=init_press_core;
md_b(1)=0;
md_c(1)=0;
M(1)=wet_mass_total;
alt(1)=0;
vel(1)=0;
event=1; % event=1 -> burn of booster and core, event=2 -> ejection boosters and core burn, event=3 -> core balistic flight 
i=0;
step=0.001;
for dt=0:step:25
    i=i+1;
    if event==1 % burn of booster and core
        % pression bottle booster
        V_b(i+1)=V_b(i)+md_b(i)/rho_water*step;
        P_b(i+1)=P_b(i)*(V_b(i)/V_b(i+1))^gamma;
        
        % pression bottle center core
        V_c(i+1)=V_c(i)+md_c(i)/rho_water*step;
        P_c(i+1)=P_c(i)*(V_c(i)/V_c(i+1))^gamma;
    
        % debit mass booster
        md_b(i+1)=sqrt((P_b(i)-P_atm)/(0.5*rho_water))*booster_throat*rho_water; % in sqrt=velocity booster water
        
        % debit mass core
        md_c(i+1)=sqrt((P_c(i)-P_atm)/(0.5*rho_water))*core_throat*rho_water; % in sqrt=velocity center core water
        
        % mass rocket at anytime
        mass_loss=(number_booster*md_b(i)+md_c(i))*step;
        M(i+1)=M(i)-mass_loss;
        mass_loss=0;
        
        % drag
        Drag(i)=0.5*Cd*rho_air*vel(i)^2*A_1; % vel=velocity of rocket 
        % gravity 
        Gravity(i)=M(i+1)*g; % M=mass of the rocket at anytime 
        % thrust
        Thrust(i)=lambda*(number_booster*2*(P_b(i+1)-P_atm)*booster_throat+2*(P_c(i+1)-P_atm)*core_throat);
        
        % Situation of the rocket when booster are stil burning
        acc(i)=(Thrust(i)-Drag(i)-Gravity(i))/M(i);
        vel(i+1)=vel(i)+acc(i)*step;
        alt(i+1)=alt(i)+0.5*acc(i)*step^2+vel(i)*step;
    end

    if event==2 % ejection booster and core burn        
        % pression bottle center core
        V_c(i+1)=V_c(i)+md_c(i)/rho_water*step;
        P_c(i+1)=P_c(i)*(V_c(i)/V_c(i+1))^gamma;
        
        % debit mass core
        md_c(i+1)=sqrt((P_c(i)-P_atm)/(0.5*rho_water))*core_throat*rho_water; % in sqrt=velocity center core water
        
        % mass rocket at anytime
        mass_loss=md_c(i)*step;
        M(i+1)=M(i)-mass_loss;
        mass_loss=0;
        
        % drag
        Drag(i)=0.5*Cd*rho_air*vel(i)^2*A_2; % vel=velocity of rocket 
        % gravity 
        Gravity(i)=M(i+1)*g; % M=mass of the rocket at anytime 
        % thrust
        Thrust(i)=lambda*(2*(P_c(i+1)-P_atm)*core_throat);
    
        % Situation of the rocket when booster are stil burning
        acc(i)=(Thrust(i)-Drag(i)-Gravity(i))/M(i);
        vel(i+1)=vel(i)+acc(i)*step;
        alt(i+1)=alt(i)+0.5*acc(i)*step^2+vel(i)*step;
    end

    if event==3 % core balistic flight until apogee
        % mass rocket at anytime
        M(i+1)=M(i);

        % drag
        Drag(i)=0.5*Cd*rho_air*vel(i)^2*A_2; % vel=velocity of rocket 
        % gravity 
        Gravity(i)=M(i)*g; % M=mass of the rocket at anytime 
        % thrust
        Thrust(i)=0;
    
        % Situation of the rocket when booster are stil burning
        acc(i)=(Thrust(i)-Drag(i)-Gravity(i))/M(i);
        vel(i+1)=vel(i)+acc(i)*step;
        alt(i+1)=alt(i)+0.5*acc(i)*step^2+vel(i)*step;
    end

    if event==4 % core under parachute at constant speed (-> no external force)
        % mass rocket at anytime
        M(i+1)=M(i);

        % drag
        Drag(i)=0; 
        % gravity 
        Gravity(i)=0;
        % thrust
        Thrust(i)=0;
    
        % Situation of the rocket when booster are stil burning
        acc(i)=0;
        vel(i+1)=-sqrt(2*M(i)*g/(rho_air*S_eff*Cd_para));
        alt(i+1)=alt(i)+vel(i)*step;
    end

%%%%%%%%% Event'conditions%%%%%%%%%
    if (sum(md_b*step)>=mass_water_booster) && (event==1) 
        event=2;
        M(i+1)=M(i+1)-number_booster*dry_mass_booster;
        burn_time_booster=i*step;
    end

    if (sum(md_c*step)>=mass_water_core) && (event==2) 
        event=3;
        M(i+1)=dry_mass_core;
        burn_time_core=i*step;
    end

    if (alt(i+1)<alt(i)) && (event==3) 
        event=4;
        M(i+1)=dry_mass_core;
        apogee_time=i*step;
    end

    if alt(i)<0
        break
    end
end 

%% Plots
t=0:step:i*step;
t2=0:step:(i-1)*step;
t3=0:step:(24.016-step);
t1000=1:1000;
ta=1:5000;

figure
plot(t,alt)
title('Altitude')
xlabel('time')
ylabel('m')

figure
plot(ta*step,vel(1:5000))
title('Velocity')
xlabel('s')
ylabel('m/s')

figure
plot(ta*step,acc(1:5000)/g)
title('Acceleration')
xlabel('s')
ylabel('g')

figure
plot(t1000*step,Thrust(1:1000))
title('Thrust')
xlabel('s')
ylabel('N')

figure
plot(t1000*step,M(1:1000)*1000)
title('Mass')
xlabel('s')
ylabel('g')

