clear ;close all;clc
number_of_Parameters=1e5;
number_of_iterations =1e2;
Font_size = 12;

%% Heavy Shadowing
m = 1; b = 0.063; omega = 0.0007;
    %% N = M = 4
    normalized_SNR_dB=linspace(-4,-1.5,40);
    Ns=16; N= 4; M=4 ; 
    E=0.5262; E2 =0.2966;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle=":",Color=[0 0.4470 0.7410],LineWidth=3,DisplayName="$N=M=4$")
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

    %% N = M = 8
    normalized_SNR_dB=linspace(-4,-2,40);
    Ns=16; N=8;M=8; 
    E=0.5262; E2 =0.2966;   
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle="-.",Color=[0.8500 0.3250 0.0980],LineWidth=3,DisplayName="$N=M=8$")
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

    %% N = M = 16
    normalized_SNR_dB=linspace(-4,-2.4,40);
    Ns=16; N=16;M=16; 
    E=0.5262; E2 =0.2966;
      
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M); 
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout,LineStyle="--",Color=[0.4660 0.6740 0.1880],LineWidth=3,DisplayName="$N=M=16$")
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
    %% N = M = 4
    normalized_SNR_dB=linspace(-4.5,-2.8,40);
    Ns=16; N=4;M=4; 
    E=0.4225; E2 =0.1869;

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

    %% N = M = 8
    normalized_SNR_dB=linspace(-4.5,-3.2,40);
    Ns=16; N=8;M=8; 
    E=0.4225; E2 =0.1869 ;  
    
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

    %% N = M = 16
    normalized_SNR_dB=linspace(-4.5,-3.4,40);
    Ns=16; N=16;M=16; 
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
lgnd = legend(["$N=M=4$","","$N=M=8$","","$N=M=16$","","","","","","","","Simulations"]);
set(lgnd,'Interpreter','latex','FontSize',Font_size)
set(gca,"FontSize",Font_size)





