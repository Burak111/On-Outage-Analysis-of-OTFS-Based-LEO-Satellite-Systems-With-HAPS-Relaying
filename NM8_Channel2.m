clear ;close all;clc
m = 8; %m>2
number_of_Parameters=1e5;
number_of_iterations =1e2;
N=8;M=8;
normalized_SNR_dB=-0:0.1:4;
normalized_SNR=10.^(normalized_SNR_dB./10);
Gamma_Rate = N*M*(m-1)*(m-2)/m;
Gamma_Scale = 1./Gamma_Rate;
Gamma_Shape = N*M*(m-2);

nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
nkgm2 = abs(nkgm).^2;
nkgm_2 = nkgm2.^-1;
phi_sim = mean(nkgm_2,2);

P_out_Teorik = zeros(length(normalized_SNR),1);
P_Out_Sim = zeros(length(normalized_SNR),1);
for kk=1:length(normalized_SNR)
    P_out_Teorik(kk) = 1 - gammainc(Gamma_Rate*normalized_SNR(kk),Gamma_Shape);

    P_Out_Sim(kk) = sum(normalized_SNR(kk)<phi_sim)./number_of_Parameters;
end

semilogy(normalized_SNR_dB,P_out_Teorik)
hold on
grid on
scatter(normalized_SNR_dB,P_Out_Sim)
ylim([1e-5,1])



