clear all;close all;clc
%% Not Fitting Shadowed Rician
clear all
N=4;M=4;    Ns=4;
m = 1; b = 0.063; omega = 0.0007;
number_of_Parameters=1e5;
Line_Width = 2;

    % Simulation
    Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
    A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
    Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
    Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
    clear A Z Uniform_RV
    Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
    clear Shadowed_Rician_Variable

    % Theoretical
    E=2.6309; E2 =10.3823;
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    PDF = @(plot_range) 1./(sqrt(Var_Phi_Theoretical*2*pi)) * exp(-1/2*((plot_range-Mean_Phi_Theoretical)./sqrt(Var_Phi_Theoretical)).^2);



    % Plotting
    figure;
    h = histogram(Phi_Rv,"Normalization","pdf");
    xlim([1 5.5])
    hold on
    plot_range=linspace(0,6,1000);
    plot(plot_range, PDF(plot_range), "LineWidth", Line_Width)
    xlabel("$\phi$","Interpreter","latex")
    ylabel("Probability Density Function")
    legend(["Simulation", "Gaussian"])
    

    % Error Calculation
    Range1 = linspace(plot_range(1),h.BinEdges(1)-h.BinWidth/2, 100 );
    Range2 = h.BinEdges(1)+h.BinWidth/2 : h.BinWidth : h.BinEdges(end);
    Range3 = linspace(h.BinEdges(end)+h.BinWidth/2, plot_range(end), 100 );
    Range = [Range1, Range2, Range3];
    PDF_Values = PDF(Range);
    Normalized_PDF_Values = PDF_Values./ sum(PDF_Values);
    Hist_Values = [zeros(1,length(Range1)), h.Values, zeros(1,length(Range3))]; 
    Normalized_Hist_Values = Hist_Values ./ sum(Hist_Values);

    Err =  sum(  (Normalized_PDF_Values-Normalized_Hist_Values).^2  ) / ( sum( (Normalized_PDF_Values-mean(Normalized_PDF_Values)).^2 ) ) ;
    
    % KL Divergence
    KL_Div = nansum( Normalized_Hist_Values.*log(Normalized_Hist_Values./Normalized_PDF_Values) );
%     title(["NMSE = ", num2str(Err),"KLDiv = ",num2str(KL_Div)]);



    
%% Not Fitting Nakagami
clear all
N=4;M=4;    m=3;
number_of_Parameters=1e5;
Line_Width = 2;

    % Simulation 
    nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
    phi_sim = mean(abs(nkgm).^-2,2);

    % Theoretical
    Gamma_Rate = N*M*(m-1)*(m-2)/m;
    Gamma_Scale = 1./Gamma_Rate;
    Gamma_Shape = N*M*(m-2);
    PDF = @(x) 10^(Gamma_Shape*log10(Gamma_Rate)-sum(log10(1:1:Gamma_Shape-1))).* x.^(Gamma_Shape-1) .* exp(-Gamma_Rate*x);

    % Plotting
    figure;
    h = histogram(phi_sim,"Normalization","pdf");
    xlim([0.4 4])
    hold on
    plot_range=linspace(0,4,1000);
    plot(plot_range,PDF(plot_range),"LineWidth",Line_Width)
    xlabel("$\phi$","Interpreter","latex")
    ylabel("Probability Density Function")
    legend(["Simulation", "Gamma"])

    % Error Calculation
    Range1 = linspace(plot_range(1),h.BinEdges(1)-h.BinWidth/2, 100 );
    Range2 = h.BinEdges(1)+h.BinWidth/2 : h.BinWidth : h.BinEdges(end);
    Range3 = linspace(h.BinEdges(end)+h.BinWidth/2, plot_range(end), 100 );
    Range = [Range1, Range2, Range3];
    PDF_Values = PDF(Range);
    Normalized_PDF_Values = PDF_Values./ sum(PDF_Values);
    Hist_Values = [zeros(1,length(Range1)), h.Values, zeros(1,length(Range3))]; 
    Normalized_Hist_Values = Hist_Values ./ sum(Hist_Values);

    Err =  sum(  (Normalized_PDF_Values-Normalized_Hist_Values).^2  ) / ( sum( (Normalized_PDF_Values-mean(Normalized_PDF_Values)).^2 ) ) ;
    

    %KLDiv
    KLDiv = nansum( Normalized_Hist_Values.*log(Normalized_Hist_Values./Normalized_PDF_Values) );
%     title(["NMSE = ", num2str(Err),"KLDiv = ",num2str(KLDiv)]);


%% Fitting Shadowed Rician
clear all
N=8;M=8;    Ns=16;
m = 1; b = 0.063; omega = 0.0007;
number_of_Parameters=1e5;
Line_Width = 2;

    % Simulation
    Uniform_RV = rand(number_of_Parameters,Ns,N*M)*2*pi;
    A = randraw('rayleigh',sqrt(b),[number_of_Parameters,Ns,N*M]);
    Z = randraw('nakagami',[m,omega],[number_of_Parameters,Ns,N*M]);
    Shadowed_Rician_Variable = A.*exp(1i*Uniform_RV) + Z*exp(1i*pi);
    clear A Z Uniform_RV
    Phi_Rv = mean((1./sum(abs(Shadowed_Rician_Variable).^2,2)),3);
    clear Shadowed_Rician_Variable

    % Theoretical
    E=0.5262; E2 =0.2966;
    Mean_Phi_Theoretical = E;
    Var_Phi_Theoretical=(E2-E^2)/(N*M);
    PDF = @(plot_range) 1./(sqrt(Var_Phi_Theoretical*2*pi)) * exp(-1/2*((plot_range-Mean_Phi_Theoretical)./sqrt(Var_Phi_Theoretical)).^2);


    % Plotting
    figure;
    h = histogram(Phi_Rv,"Normalization","pdf");
    xlim([0.45 0.6])
    hold on
    plot_range=linspace(0.4,0.65,1000);
    plot(plot_range,PDF(plot_range),"LineWidth",Line_Width)
    xlabel("$\phi$","Interpreter","latex")
    ylabel("Probability Density Function")
    legend(["Simulation", "Gaussian"])

    % Error Calculation
    Range1 = linspace(plot_range(1),h.BinEdges(1)-h.BinWidth/2, 100 );
    Range2 = h.BinEdges(1)+h.BinWidth/2 : h.BinWidth : h.BinEdges(end);
    Range3 = linspace(h.BinEdges(end)+h.BinWidth/2, plot_range(end), 100 );
    Range = [Range1, Range2, Range3];
    PDF_Values = PDF(Range);
    Normalized_PDF_Values = PDF_Values./ sum(PDF_Values);
    Hist_Values = [zeros(1,length(Range1)), h.Values, zeros(1,length(Range3))]; 
    Normalized_Hist_Values = Hist_Values ./ sum(Hist_Values);

    Err =  sum(  (Normalized_PDF_Values-Normalized_Hist_Values).^2  ) / ( sum( (Normalized_PDF_Values-mean(Normalized_PDF_Values)).^2 ) ) ;

    % KL Divergence
    KL_Div = nansum( Normalized_Hist_Values.*log(Normalized_Hist_Values./Normalized_PDF_Values) );

%     title(["NMSE = ", num2str(Err),"KLDiv = ",num2str(KL_Div)]);

%% Fitting Nakagami
clear all
N=8;M=8;    m=8;
number_of_Parameters=1e5;
Line_Width = 2;

    % Simulation 
    nkgm = randraw('nakagami',[m,1],[number_of_Parameters,N*M]);
    phi_sim = mean(abs(nkgm).^-2,2);

    % Theoretical
    Gamma_Rate = N*M*(m-1)*(m-2)/m;
    Gamma_Scale = 1./Gamma_Rate;
    Gamma_Shape = N*M*(m-2);
%     PDF = @(x) Gamma_Rate^Gamma_Shape/gamma(Gamma_Shape) .* x.^(Gamma_Shape-1) .* exp(-Gamma_Rate*x);
    PDF = @(x) 10^(Gamma_Shape*log10(Gamma_Rate)-sum(log10(1:1:Gamma_Shape-1))).* x.^(Gamma_Shape-1) .* exp(-Gamma_Rate*x);

    % Plotting
    figure;
    h = histogram(phi_sim,"Normalization","pdf");
    xlim([0.9 1.4])
    hold on
    plot_range=linspace(0.9,1.5,1000);
    plot(plot_range,PDF(plot_range),"LineWidth",Line_Width)
    xlabel("$\phi$","Interpreter","latex")
    ylabel("Probability Density Function")
    legend(["Simulation", "Gamma"])

    % Error Calculation
    Range1 = linspace(plot_range(1),h.BinEdges(1)-h.BinWidth/2, 100 );
    Range2 = h.BinEdges(1)+h.BinWidth/2 : h.BinWidth : h.BinEdges(end);
    Range3 = linspace(h.BinEdges(end)+h.BinWidth/2, plot_range(end), 100 );
    Range = [Range1, Range2, Range3];
    PDF_Values = PDF(Range);
    Normalized_PDF_Values = PDF_Values./ sum(PDF_Values);
    Hist_Values = [zeros(1,length(Range1)), h.Values, zeros(1,length(Range3))]; 
    Normalized_Hist_Values = Hist_Values ./ sum(Hist_Values);

    Err =  sum(  (Normalized_PDF_Values-Normalized_Hist_Values).^2  ) / ( sum( (Normalized_PDF_Values-mean(Normalized_PDF_Values)).^2 ) ) ;
    

    %KLDiv
    KLDiv = nansum( Normalized_Hist_Values.*log(Normalized_Hist_Values./Normalized_PDF_Values) );
%     title(["NMSE = ", num2str(Err),"KLDiv = ",num2str(KLDiv)]);
%%










    