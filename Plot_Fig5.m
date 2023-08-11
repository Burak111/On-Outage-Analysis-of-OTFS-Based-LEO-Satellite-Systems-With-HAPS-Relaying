clear all;close all;clc
%% Ns = 8, N = M = 8, NKGM-m = 8
number_of_iterations =1e5;
number_of_Parameters =1e2;
SNR_Resolution = 20;
    %% Channel 1 
    normalized_SNR_dB=linspace(0,2,SNR_Resolution);
    Ns=8; N= 8; M=8 ; 
    m = 1; b = 0.063; omega = 0.0007;
    E=1.1275; E2 =1.4832;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
 
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout_Teorik = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
%     semilogy(normalized_SNR_dB,Pout,LineStyle=":",Color=[0 0.4470 0.7410],LineWidth=3,DisplayName="$N=M=4$")

    
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
    Pout_Sim = mean(Pout_Sim,2).';
%     scatter(normalized_SNR_dB,mean(Pout_Sim,2),[],[0 0.4470 0.7410],"filled","diamond")
%     pause(1)
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv
    %% Channel 2
    m = 8;
    Gamma_Rate = N*M*(m-1)*(m-2)/m;
    Gamma_Scale = 1./Gamma_Rate;
    Gamma_Shape = N*M*(m-2);

    for kk=1:length(normalized_SNR)
        P_Nout_Teorik(kk) = gammainc(Gamma_Rate*normalized_SNR(kk),Gamma_Shape);
        for iter=1:number_of_iterations
            nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
            phi_sim = mean(abs(nkgm).^-2,2);
            P_Out_Sim(kk,iter) = sum(normalized_SNR(kk)<phi_sim)./number_of_Parameters;
        end
    end
    P_Out_Sim = mean(P_Out_Sim,2).';
    P_out_Teorik = 1-P_Nout_Teorik;
    %% Plotting Outage
    Pout_E2E_Teorik =  P_out_Teorik + Pout_Teorik - P_out_Teorik .* Pout_Teorik;
    Pout_E2E_Sim = Pout_Sim + P_Out_Sim - Pout_Sim .* P_Out_Sim;
%     figure;
    semilogy(normalized_SNR_dB,Pout_E2E_Teorik,LineStyle="-",Color=[0 0.4470 0.7410],LineWidth=3)
    hold on
    ylim([1e-5 1])
    grid on
    scatter(normalized_SNR_dB,Pout_E2E_Sim,[],[0 0 0],"filled","diamond")

%% Ns = 8, N = M = 16, NKGM-m = 8
clear all;clc
number_of_iterations =1e5;
number_of_Parameters =1e2;
SNR_Resolution = 20;
    %% Channel 1 
    normalized_SNR_dB=linspace(0,1.5,SNR_Resolution);
    Ns = 8; N = 16; M = 16; 
    m = 1; b = 0.063; omega = 0.0007;
    E=1.1275; E2 =1.4832;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout_Teorik = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));

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
    Pout_Sim = mean(Pout_Sim,2).';
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Channel 2
    m = 8;
    Gamma_Rate = N*M*(m-1)*(m-2)/m;
    Gamma_Scale = 1./Gamma_Rate;
    Gamma_Shape = N*M*(m-2);
    
    for kk=1:length(normalized_SNR)
        P_Nout_Teorik(kk) = gammainc(Gamma_Rate*normalized_SNR(kk),Gamma_Shape);
        for iter=1:number_of_iterations
            nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
            phi_sim = mean(abs(nkgm).^-2,2);
            P_Out_Sim(kk,iter) = sum(normalized_SNR(kk)<phi_sim)./number_of_Parameters;
        end
    end
    P_Out_Sim = mean(P_Out_Sim,2).';
    P_out_Teorik = 1-P_Nout_Teorik;
    %% Plotting Outage
    Pout_E2E_Teorik =  P_out_Teorik + Pout_Teorik - P_out_Teorik .* Pout_Teorik;
    Pout_E2E_Sim = Pout_Sim + P_Out_Sim - Pout_Sim .* P_Out_Sim;
%     figure;
    semilogy(normalized_SNR_dB,Pout_E2E_Teorik,LineStyle="--",Color=[0.9290 0.6940 0.1250]	,LineWidth=3)
    hold on
    scatter(normalized_SNR_dB,Pout_E2E_Sim,[],[0 0 0]	,"filled","diamond")
    ylim([1e-5 1])
    grid on

%% Ns = 16, N = M = 8, NKGM-m = 8
clear all;clc
number_of_iterations =1e5;
number_of_Parameters =1e2;
SNR_Resolution = 20;
    %% Channel 1 
    normalized_SNR_dB=linspace(0,4,SNR_Resolution);
    Ns = 16; N = 8; M = 8; 
    m = 1; b = 0.063; omega = 0.0007;
    E=0.5262; E2 =0.2966;

    m = 1; b = 0.063; omega = 0.0007;
    alpha =(2*b*m/(2*b*m+omega))^m/(2*b);
    Beta = 1/(2*b);
    c = omega/(2*b*(2*b*m+omega));
    ksi = @(k) (-1)^k * pochhammer(1-m,k) * c^k / factorial(k)^2 ;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout_Teorik = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));

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
    Pout_Sim = mean(Pout_Sim,2).';
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Channel 2
    m = 8;
    Gamma_Rate = N*M*(m-1)*(m-2)/m;
    Gamma_Scale = 1./Gamma_Rate;
    Gamma_Shape = N*M*(m-2);
    
    for kk=1:length(normalized_SNR)
        P_Nout_Teorik(kk) = gammainc(Gamma_Rate*normalized_SNR(kk),Gamma_Shape);
        for iter=1:number_of_iterations
            nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
            phi_sim = mean(abs(nkgm).^-2,2);
            P_Out_Sim(kk,iter) = sum(normalized_SNR(kk)<phi_sim)./number_of_Parameters;
        end
    end
    P_Out_Sim = mean(P_Out_Sim,2).';
    P_out_Teorik = 1-P_Nout_Teorik;
    %% Plotting Outage
    Pout_E2E_Teorik =  P_out_Teorik + Pout_Teorik - P_out_Teorik .* Pout_Teorik;
    Pout_E2E_Sim = Pout_Sim + P_Out_Sim - Pout_Sim .* P_Out_Sim;
%     figure;
    semilogy(normalized_SNR_dB,Pout_E2E_Teorik,LineStyle=":",Color=[0.4940 0.1840 0.5560]	,LineWidth=3)
    hold on
    scatter(normalized_SNR_dB,Pout_E2E_Sim,[],[0 0 0]	,"filled","diamond")
    ylim([1e-5 1])
    grid on

%% Ns = 16, N = M = 16, NKGM-m = 8
clear all;clc
number_of_iterations =1e5;
number_of_Parameters =1e2;
SNR_Resolution = 20;
    %% Channel 1 
    normalized_SNR_dB=linspace(0,2.5,SNR_Resolution);
    Ns = 16; N = 16; M = 16; 
    m = 1; b = 0.063; omega = 0.0007;
    E=0.5262; E2 =0.2966;

    m = 1; b = 0.063; omega = 0.0007;
    alpha =(2*b*m/(2*b*m+omega))^m/(2*b);
    Beta = 1/(2*b);
    c = omega/(2*b*(2*b*m+omega));
    ksi = @(k) (-1)^k * pochhammer(1-m,k) * c^k / factorial(k)^2 ;

    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout_Teorik = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));

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
    Pout_Sim = mean(Pout_Sim,2).';
    clear Uniform_RV A Z Shadowed_Rician_Variable Phi_Rv

    %% Channel 2
    m = 8;
    Gamma_Rate = N*M*(m-1)*(m-2)/m;
    Gamma_Scale = 1./Gamma_Rate;
    Gamma_Shape = N*M*(m-2);

    for kk=1:length(normalized_SNR)
        P_Nout_Teorik(kk) = gammainc(Gamma_Rate*normalized_SNR(kk),Gamma_Shape);
        for iter=1:number_of_iterations
            nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
            phi_sim = mean(abs(nkgm).^-2,2);
            P_Out_Sim(kk,iter) = sum(normalized_SNR(kk)<phi_sim)./number_of_Parameters;
        end
    end
    P_Out_Sim = mean(P_Out_Sim,2).';
    P_out_Teorik = 1-P_Nout_Teorik;
    %% Plotting Outage
    Pout_E2E_Teorik =  P_out_Teorik + Pout_Teorik - P_out_Teorik .* Pout_Teorik;
    Pout_E2E_Sim = Pout_Sim + P_Out_Sim - Pout_Sim .* P_Out_Sim;
%     figure;
    semilogy(normalized_SNR_dB,Pout_E2E_Teorik,LineStyle="-.",Color=[0.6350 0.0780 0.1840]	,LineWidth=3)
    hold on
    scatter(normalized_SNR_dB,Pout_E2E_Sim,[],[0 0 0]	,"filled","diamond")
    ylim([1e-5 1])
    grid on
%% Plotting fixis
Font_size = 12;
dummyh = scatter(nan, nan,[],[0 0 0], "filled","square");
lgnd = legend(["$N_S=8, N=M=8, m=8$","","$N_S=8, N=M=16, m=8$","","$N_S=16, N=M=8, m=8$","","$N_S=16, N=M=16, m=8$","","Simulations"],Interpreter="latex");
set(gca,"FontSize",Font_size)





