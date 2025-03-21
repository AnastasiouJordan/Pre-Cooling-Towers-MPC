function [p,u] = ARModels(t, u, p)

% ARIMA Models

% ARIMA MODEL FOR F_inPT
% Load the data
Raw.F_in  = u.F_outSD(t) + u.F_AntiSurge(t);
Size_F_in = size(Raw.F_in,1);
t1 = t(1:51086,1);

% Create Model Template
%Mdl_F_in = arima(1,0,0); % First-order Auto-Regressive Model plus constant
dep = 0.97; % Dependence of the next L measurement on the previous one
Mdl_F_in = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degreee
% D-value: MA - the moving average degree of integration
% q-value: ARMA - the moving average polynomial degree

% Here, we have:
% 1 nonseasonal AR polynomial lag
% 0 degree nonseasonal integration polynomial
% 0 nonseasonal MA polynomial lags


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_in = 1:Mdl_F_in.P;               % Contains 2 observations
estsample_F_in = (Mdl_F_in.P + 1):Size_F_in; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_in = estimate(Mdl_F_in, Raw.F_in(estsample_F_in), 'Y0', Raw.F_in(presample_F_in));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_F_in = cell2mat(EstMdl_F_in.AR(1)); % AR coefficient, alpha
p.w_F_in = EstMdl_F_in.Variance;        % Variance which is rooted to obtain error
p.c_F_in = EstMdl_F_in.Constant;        % Constant or intercept


% ARIMA MODEL FOR F_outPT

% Load the data
Raw.F_out = u.F_outPT(t);
Size_F_out = size(Raw.F_out,1);

% Create Model Template
Mdl_F_out = arima(1,0,0); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degreee
% D-value: MA - the moving average degree of integration
% q-value: ARMA - the moving average polynomial degree

% Here, we have:
% 1 nonseasonal AR polynomial lag
% 0 degree nonseasonal integration polynomial
% 0 nonseasonal MA polynomial lags


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_out = 1:Mdl_F_out.P;                % Contains 2 observations
estsample_F_out = (Mdl_F_out.P + 1):Size_F_out; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_out = estimate(Mdl_F_out, Raw.F_out(estsample_F_out), 'Y0', Raw.F_out(presample_F_out));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_F_out = cell2mat(EstMdl_F_out.AR(1)); % AR coefficient, alpha
p.w_F_out = EstMdl_F_out.Variance;        % Variance which is rooted to obtain error
p.c_F_out = EstMdl_F_out.Constant;        % Constant or intercept


% ARIMA MODEL FOR L_PT

% Load the data
Raw_L = (u.L_PTA(t) + u.L_PTB(t))./2; % Raw level measurements as percentages
Raw.L = Raw_L/100*p.height_PT;        % Raw level measurements as heighs in meters
Size_L = size(Raw.L,1);

% Create Model Template
dep = 0.8; % Dependence of the next L measurement on the previous one
%Mdl_L_A = arima(1,0,0); % This gives an alpha value of 0.9999 and poor fit
Mdl_L = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degreee
% D-value: MA - the moving average degree of integration
% q-value: ARMA - the moving average polynomial degree

% Here, we have:
% 1 nonseasonal AR polynomial lag
% 0 degree nonseasonal integration polynomial
% 0 nonseasonal MA polynomial lags
% We also specify the alpha value as 0.99


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_L = 1:Mdl_L.P;            % Contains 2 observations
estsample_L = (Mdl_L.P + 1):Size_L; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_L = estimate(Mdl_L, Raw.L(estsample_L), 'Y0', Raw.L(presample_L));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_L    = cell2mat(EstMdl_L.AR(1)); % AR coefficient, alpha
p.w_L    = EstMdl_L.Variance;        % Variance which is rooted to obtain error
p.c_L    = EstMdl_L.Constant;        % Constant or intercept
deltaT = t(2) - t(1);                % s, Delta time in seconds
p.C_L    = (deltaT/p.area_PT)*0.001;   % m, Coefficient to the flowrates.
                                     % This calculation takes delta time and
                                     % divides it by the area of the Chill Dam
                                     % in m2. This is then multiplied by 0.001
                                     % to adjust for the flowrate which is in
                                     % L/s. This gives a constant in meters,
                                     % equivalent to the units in which level
                                     % is measured.


% AR Model Predictions
k.F_in(1)  = Raw.F_in(1);  % Specify the inital value of the inlet flowrate
k.F_out(1) = Raw.F_out(1); % Specify the inital value of the outlet flowrate
k.L(1)     = Raw.L(1);     % Specify the inital value of the level
for i = 1:51086            % From here onwards, make predictions
    e.F_in(i)    = randn*sqrt(p.w_F_in);
    k.F_in(i+1)  = p.c_F_in + p.a_F_in*k.F_in(i) + e.F_in(i);
    e.F_out(i)   = randn*sqrt(p.w_F_out);
    k.F_out(i+1) = p.c_F_out + p.a_F_out*k.F_out(i) + e.F_out(i);
    e.L(i)       = randn*sqrt(p.w_L);
    k.L(i+1)     = p.a_L*k.L(i) + p.C_L*(k.F_in(i) - k.F_out(i)) + e.L(i) + p.c_L;
end

g.F_in(1)  = Raw.F_in(1);  % Specify the inital value of the inlet flowrate
g.F_out(1) = Raw.F_out(1); % Specify the inital value of the outlet flowrate
g.L(1)     = Raw.L(1);     % Specify the inital value of the level
for i = 1:51086           
    e.F_in(i)    = randn*sqrt(p.w_F_in);
    g.F_in(i+1)  = p.c_F_in + p.a_F_in*g.F_in(i) + e.F_in(i);    % Generated data used to test the model in simulation
    e.F_out(i)   = randn*sqrt(p.w_F_out);
    g.F_out(i+1) = p.c_F_out + p.a_F_out*g.F_out(i) + e.F_out(i); % Generated data used to test the model in simulation
    e.L(i)       = randn*sqrt(p.w_L);
    g.L(i+1)     = p.a_L*g.L(i) + p.C_L*(g.F_in(i) - g.F_out(i)) + e.L(i) + p.c_L; % Generated data used to test the model in simulation
end
% Where in the above loop:
% e(i) is the error obtained from the variance calculated by the model for each variable.
% c is the constant or intercept calculated by the model for each variable.
% a is the 'alpha' co-efficient for the dependence of the next value on the previous one for each variable.
% k(i) is the prediction for each variable.
% C is the constant used only in the level predictions, derived from the relationship between level and flowrates.


% g_F_in = g.F_in';
% Size_g_F_in = size(g_F_in,1);
% ipt_g = findchangepts(g_F_in,'Statistic','std','MinThreshold',150);
% intvls_g = [1 ipt_g' length(g_F_in)];
% for k1_g = 1:length(intvls_g)-2
%     interval_means_g(k1_g) = mean(g_F_in(intvls_g(k1_g):intvls_g(k1_g+1)-1));
% end
% interval_means_g(k1_g+1) = mean(g_F_in(intvls_g(k1_g):intvls_g(end)));
% xval_g = [];
% yval_g = [];
% for k1_g = 1:length(intvls_g)-1
%     xval_g = intvls_g(k1_g):intvls_g(k1_g+1);
% ynext_g = interval_means_g(k1_g)*ones(1,(length(xval_g)-1));
% yval_g = [yval_g ynext_g];
% end
% yval_g(end+1) = 0;
% g_F_in = yval_g;

u.F_in_generated = griddedInterpolant(t, g.F_in', "previous");
u.F_out_generated = griddedInterpolant(t, g.F_out');
% % INLET FLOWRATE
% % Plot Actual Data vs AR Model
figure(1)
title('Inlet Flowrate AR Model');
plot(t, Raw.F_in, 'r', t, u.F_in_generated(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_i_n_P_T (L/s)')
% 
% % Plot Observations and Fitted Values
% subplot(3,1,2)
% resid_F_in = infer(EstMdl_F_in, Raw.F_in(estsample_F_in), ...
%              'Y0', Raw.F_in(presample_F_in));      % Infer residuals from
%                                                    % the estimated model
% yhat_F_in = Raw.F_in(estsample_F_in) - resid_F_in; % Compute the fitted values
% plot(t1, Raw.F_in(estsample_F_in),'r', t1, yhat_F_in,'b--','LineWidth',1)
% legend('Observations', 'Fitted Values')
% 
% % Plot Residuals vs Fitted Values
% subplot(3,1,3)
% plot(yhat_F_in,resid_F_in,'b.')
% ylabel('Residuals')
% xlabel('Fitted Values')
% 
% % OUTLET FLOWRATE
% % Plot Actual Data vs AR Model
% figure(2)
% title('Outlet Flowrate AR Model');
% subplot(3,1,1)
% plot(t, Raw.F_out, 'r', t, k.F_out', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (s)')
% ylabel('F_o_u_t_P_T (L/s)')
% 
% % Plot Observations and Fitted Values
% subplot(3,1,2)
% resid_F_out = infer(EstMdl_F_out, Raw.F_out(estsample_F_out), ...
%              'Y0', Raw.F_out(presample_F_out));        % Infer residuals from
%                                                        % the estimated model
% yhat_F_out = Raw.F_out(estsample_F_out) - resid_F_out; % Compute the fitted values
% plot(t1, Raw.F_out(estsample_F_out),'r', t1, yhat_F_out,'b--','LineWidth',1)
% legend('Observations', 'Fitted Values')
% 
% % Plot Residuals vs Fitted Values
% subplot(3,1,3)
% plot(yhat_F_out,resid_F_out,'b.')
% ylabel('Residuals')
% xlabel('Fitted Values')
% 
% % LEVEL
% % Plot Actual Data vs AR Model
% figure(3)
% title('Level AR Model');
% subplot(3,1,1)
% plot(t, Raw.L, 'r', t, k.L', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (s)')
% ylabel('L_P_T (m)')
% 
% % Plot Observations and Fitted Values
% subplot(3,1,2)
% resid_L = infer(EstMdl_L, Raw.L(estsample_L), ...
%              'Y0', Raw.L(presample_L)); % Infer residuals from
%                                         % the estimated model
% yhat_L = Raw.L(estsample_L) - resid_L;  % Compute the fitted values
% plot(t1, Raw.L(estsample_L),'r', t1, yhat_L,'b--','LineWidth',1)
% legend('Observations', 'Fitted Values')
% 
% % Plot Residuals vs Fitted Values
% subplot(3,1,3)
% plot(yhat_L, resid_L,'b.')
% ylabel('Residuals')
% xlabel('Fitted Values')


% Convert Level back to %
% figure(4)
% title('Level AR Model in %');
% k_L_percent     = k.L*100/p.height_PT;
% L_Data_percent  = Raw.L*100/p.height_PT;
% plot(t, L_Data_percent, 'r', t, k_L_percent', 'b');
% legend('Raw Data', 'AR Model');
% xlabel('Time (s)');
% ylabel('L_P_T (%)');



end
