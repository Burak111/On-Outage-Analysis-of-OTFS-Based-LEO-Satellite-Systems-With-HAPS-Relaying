clear ;close all;clc
number_of_Parameters=1e5;
number_of_iterations =1e2;
N=8;M=8;
normalized_SNR_dB=-4:0.2:10;

%% Heavy Shadowing
m = 1; b = 0.063; omega = 0.0007;
alpha =(2*b*m/(2*b*m+omega))^m/(2*b);
Beta = 1/(2*b);
c = omega/(2*b*(2*b*m+omega));
ksi = @(k) (-1)^k * pochhammer(1-m,k) * c^k / factorial(k)^2 ;
    %% Ns = 4
    Ns=4;
    E=0;E2=0;
    for k1=0:m-1
        for k2=0:m-1
            for k3=0:m-1
                for k4=0:m-1
                    E = E + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*alpha^4*beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*...
                        gamma(k1+k2+k3+k4+4-1)/(Beta-c)^(k1+k2+k3+k4+4-1);
    
                    E2 = E2 + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*alpha^4*beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*...
                        gamma(k1+k2+k3+k4+4-2)/(Beta-c)^(k1+k2+k3+k4+4-2);
                end
            end
        end
    end
    
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout)
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
    scatter(normalized_SNR_dB,mean(Pout_Sim,2))
    %% Ns = 8
    Ns=8;
    E=0;E2=0;
    for k1=0:m-1
        for k2=0:m-1
            for k3=0:m-1
                for k4=0:m-1
                    for k5=0:m-1
                        for k6=0:m-1
                            for k7=0:m-1
                                for k8=0:m-1
    E = E + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*alpha^8*...
        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*...
        gamma(k1+k2+k3+k4+k5+k6+k7+k8+8-1)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+8-1);
    
    E2 = E2 + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*alpha^8*...
        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*...
        gamma(k1+k2+k3+k4+k5+k6+k7+k8+8-2)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+8-2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout)
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
    scatter(normalized_SNR_dB,mean(Pout_Sim,2))
    %% Ns = 16
    Ns=16;
    E=0;E2=0;
    for k1=0:m-1
        for k2=0:m-1
            for k3=0:m-1
                for k4=0:m-1
                    for k5=0:m-1
                        for k6=0:m-1
                            for k7=0:m-1
                                for k8=0:m-1
                                    for k9=0:m-1
                                        for k10=0:m-1
                                            for k11=0:m-1
                                                for k12=0:m-1
                                                    for k13=0:m-1
                                                        for k14=0:m-1
                                                            for k15=0:m-1
                                                                for k16=0:m-1
                                                                    E = E + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*ksi(k9)*ksi(k10)*ksi(k11)*ksi(k12)*ksi(k13)*ksi(k14)*ksi(k15)*ksi(k16)*alpha^16*...
                                                                        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+8,k9+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+9,k10+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+10,k11+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+11,k12+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+12,k13+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+13,k14+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+14,k15+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+15,k16+1)*...
                                                                        gamma(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-1)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-1);
    
                                                                    E2 = E2 + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*ksi(k9)*ksi(k10)*ksi(k11)*ksi(k12)*ksi(k13)*ksi(k14)*ksi(k15)*ksi(k16)*alpha^16*...
                                                                        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+8,k9+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+9,k10+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+10,k11+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+11,k12+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+12,k13+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+13,k14+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+14,k15+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+15,k16+1)*...
                                                                        gamma(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-2)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-2);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout)
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
    scatter(normalized_SNR_dB,mean(Pout_Sim,2))


%%  Karasawa
m = 2; b = 0.0158; omega = 0.123;
alpha =(2*b*m/(2*b*m+omega))^m/(2*b);
Beta = 1/(2*b);
c = omega/(2*b*(2*b*m+omega));
ksi = @(k) (-1)^k * pochhammer(1-m,k) * c^k / factorial(k)^2 ;
normalized_SNR_dB=-4:0.2:10;
    %% Ns = 4
    Ns=4;
    E=0;E2=0;
    for k1=0:m-1
        for k2=0:m-1
            for k3=0:m-1
                for k4=0:m-1
                    E = E + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*alpha^4*beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*...
                        gamma(k1+k2+k3+k4+4-1)/(Beta-c)^(k1+k2+k3+k4+4-1);
    
                    E2 = E2 + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*alpha^4*beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*...
                        gamma(k1+k2+k3+k4+4-2)/(Beta-c)^(k1+k2+k3+k4+4-2);
                end
            end
        end
    end
    
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout)
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
    scatter(normalized_SNR_dB,mean(Pout_Sim,2))
    
    %% Ns = 8
    Ns=8;
    E=0;E2=0;
    for k1=0:m-1
        for k2=0:m-1
            for k3=0:m-1
                for k4=0:m-1
                    for k5=0:m-1
                        for k6=0:m-1
                            for k7=0:m-1
                                for k8=0:m-1
    E = E + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*alpha^8*...
        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*...
        gamma(k1+k2+k3+k4+k5+k6+k7+k8+8-1)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+8-1);
    
    E2 = E2 + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*alpha^8*...
        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*...
        gamma(k1+k2+k3+k4+k5+k6+k7+k8+8-2)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+8-2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout)
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
    scatter(normalized_SNR_dB,mean(Pout_Sim,2))
    %% Ns = 16
    Ns=16;
    E=0;E2=0;
    for k1=0:m-1
        for k2=0:m-1
            for k3=0:m-1
                for k4=0:m-1
                    for k5=0:m-1
                        for k6=0:m-1
                            for k7=0:m-1
                                for k8=0:m-1
                                    for k9=0:m-1
                                        for k10=0:m-1
                                            for k11=0:m-1
                                                for k12=0:m-1
                                                    for k13=0:m-1
                                                        for k14=0:m-1
                                                            for k15=0:m-1
                                                                for k16=0:m-1
                                                                    E = E + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*ksi(k9)*ksi(k10)*ksi(k11)*ksi(k12)*ksi(k13)*ksi(k14)*ksi(k15)*ksi(k16)*alpha^16*...
                                                                        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+8,k9+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+9,k10+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+10,k11+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+11,k12+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+12,k13+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+13,k14+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+14,k15+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+15,k16+1)*...
                                                                        gamma(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-1)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-1);
    
                                                                    E2 = E2 + ksi(k1)*ksi(k2)*ksi(k3)*ksi(k4)*ksi(k5)*ksi(k6)*ksi(k7)*ksi(k8)*ksi(k9)*ksi(k10)*ksi(k11)*ksi(k12)*ksi(k13)*ksi(k14)*ksi(k15)*ksi(k16)*alpha^16*...
                                                                        beta(k1+1,k2+1)*beta(k1+k2+2,k3+1)*beta(k1+k2+k3+3,k4+1)*beta(k1+k2+k3+k4+4,k5+1)*beta(k1+k2+k3+k4+k5+5,k6+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+6,k7+1)*beta(k1+k2+k3+k4+k5+k6+k7+7,k8+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+8,k9+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+9,k10+1)*beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+10,k11+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+11,k12+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+12,k13+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+13,k14+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+14,k15+1)*...
                                                                        beta(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+15,k16+1)*...
                                                                        gamma(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-2)/(Beta-c)^(k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11+k12+k13+k14+k15+k16+16-2);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    
    
    normalized_SNR = 10.^(normalized_SNR_dB./10);
    Pout = qfunc((normalized_SNR-Mean_Phi_Theoretical)/sqrt(Var_Phi_Theoretical));
    semilogy(normalized_SNR_dB,Pout)
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
    scatter(normalized_SNR_dB,mean(Pout_Sim,2))


