clear all;close all;clc
number_of_Parameters=1e5;
number_of_iterations =1e2;
Font_size = 12;

%% Heavy Shadowing
m = 1; b = 0.063; omega = 0.0007;
    %% Ns = 4
    normalized_SNR_dB=linspace(2,6,20);
    Ns=4; N= 8; M=8 ; 
    E=2.6309; E2 =10.3823;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle=":",Color=[0 0.4470 0.7410],LineWidth=3,DisplayName="$N_S=4$")
    ylim([1e-5 1])
    grid on
    hold on
    
    Pout_Sim=zeros(length(normalized_SNR),number_of_iterations);
    for iter=1:number_of_iterations
        Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
        A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
        Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
        Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
        Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
        for kk=1:length(normalized_SNR)
            Pout_Sim(kk,iter) = Pout_Sim(kk,iter) + sum(normalized_SNR(kk)<1*Phi_Rv)/number_of_Parameters;
        end
    end
    scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0 0.4470 0.7410],"filled","diamond")
    pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Ns = 8
    normalized_SNR_dB=linspace(0,1.6,20);
    Ns=8; N=8;M=8; 
    E=1.1275; E2 =1.4832;   
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle="-.",Color=[0.8500 0.3250 0.0980],LineWidth=3,DisplayName="$N_S=8$")
    ylim([1e-5 1])
    grid on
    hold on
    
    Pout_Sim=zeros(length(normalized_SNR),number_of_iterations);
    for iter=1:number_of_iterations
        Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
        A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
        Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
        Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
        Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
        for kk=1:length(normalized_SNR)
            Pout_Sim(kk,iter) = Pout_Sim(kk,iter) + sum(normalized_SNR(kk)<1*Phi_Rv)/number_of_Parameters;
        end
    end
    scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0.8500 0.3250 0.0980],"filled","diamond") 
    pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Ns = 16
    normalized_SNR_dB=linspace(-3,-2,20);
    Ns=16; N=8;M=8; 
    E=0.5262; E2 =0.2966;
      
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M); 
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle="--",Color=[0.4660 0.6740 0.1880],LineWidth=3,DisplayName="$N_S=16$")
    ylim([1e-5 1])
    grid on
    hold on
    
    Pout_Sim=zeros(length(normalized_SNR),number_of_iterations);
    for iter=1:number_of_iterations
        Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
        A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
        Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
        Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
        Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
        for kk=1:length(normalized_SNR)
            Pout_Sim(kk,iter) = Pout_Sim(kk,iter) + sum(normalized_SNR(kk)<1*Phi_Rv)/number_of_Parameters;
        end
    end
    scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0.4660 0.6740 0.1880],"filled","diamond")
    pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

%% Karasawa Shadowing
m=2; b=0.0158; omega=0.1230;
alpha =(2*b*m/(2*b*m+omega))^m/(2*b);
Beta = 1/(2*b);
c = omega/(2*b*(2*b*m+omega));
ksi = @(k) (-1)^k * pochhammer(1-m,k) * c^k / factorial(k)^2 ;
    %% Ns = 4
    normalized_SNR_dB=linspace(1,5,20);
    Ns=4; N=8;M=8; 
    E=1.9657; E2 =5.0167;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle=":",Color=[0 0.4470 0.7410],LineWidth=3)  
    ylim([1e-5 1])
    grid on
    hold on
    
    Pout_Sim=zeros(length(normalized_SNR),number_of_iterations);
    for iter=1:number_of_iterations
        Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
        A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
        Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
        Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
        Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
        for kk=1:length(normalized_SNR)
            Pout_Sim(kk,iter) = Pout_Sim(kk,iter) + sum(normalized_SNR(kk)<1*Phi_Rv)/number_of_Parameters;
        end
    end
    scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0 0.4470 0.7410],"filled","diamond")
    pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Ns = 8
    normalized_SNR_dB=linspace(-1.5,0.5,20);
    Ns=8; N=8;M=8; 
    E=0.8854; E2 =0.8682;  
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle="-.",Color=[0.8500 0.3250 0.0980],LineWidth=3)
    ylim([1e-5 1])
    grid on
    hold on
    
    Pout_Sim=zeros(length(normalized_SNR),number_of_iterations);
    for iter=1:number_of_iterations
        Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
        A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
        Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
        Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
        Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
        for kk=1:length(normalized_SNR)
            Pout_Sim(kk,iter) = Pout_Sim(kk,iter) + sum(normalized_SNR(kk)<1*Phi_Rv)/number_of_Parameters;
        end
    end
    scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0.8500 0.3250 0.0980],"filled","diamond")   
    pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Ns = 16
    normalized_SNR_dB=linspace(-4.5,-3,20);
    Ns=16; N=8;M=8; 
    E=0.4225; E2 =0.1869;
      
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M); 
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle="--",Color=[0.4660 0.6740 0.1880],LineWidth=3)
    ylim([1e-5 1])
    grid on
    hold on
    
    Pout_Sim=zeros(length(normalized_SNR),number_of_iterations);
    for iter=1:number_of_iterations
        Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
        A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
        Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
        Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
        Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
        for kk=1:length(normalized_SNR)
            Pout_Sim(kk,iter) = Pout_Sim(kk,iter) + sum(normalized_SNR(kk)<1*Phi_Rv)/number_of_Parameters;
        end
    end
    scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0.4660 0.6740 0.1880],"filled","diamond")
    pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

%% Plotting fixis
dummyh = scatter(nan, nan,[],[0 0 0], "filled","diamond","DisplayName","Simulation");
lgnd = legend(["$N_S=4$","","$N_S=8$","","$N_S=16$","","","","","","","","Simulations"]);
set(lgnd,'Interpreter','latex','FontSize',Font_size)
set(gca,"FontSize",Font_size)





