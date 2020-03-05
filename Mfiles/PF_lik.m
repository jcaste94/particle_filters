function [lik, all_s_up, Neff] = PF_lik(A, B, H, Phi, R, S2, N, yt, x0, P0, naive)
% -----------------------------------------------------------------------
% particle filtering for DSGE model with measurement error
%
% Author: Minsu Chang        minsuc@sas.upenn.edu
% Last modified: 2/24/2016
% -----------------------------------------------------------------------
%                        Model( Linear Model)
%               y(t) = A + B*s(t) + u(t), u(t) ~ N(0,H)
%           s(t) = Phi*s(t-1) + R*e(t), e(t) ~ iid N(0,S2)
% -----------------------------------------------------------------------
% Input
%   A, B, H    : parameters for measurement equation
%   Phi, R, S2 : parameters for transition equation
%   N          : number of particles
%   yt         : data
%   x0, P0     : initial prior for the filter
%   naive      : ==1 (bootstrap PF), !=1 (conditionally optimal PF)
% -----------------------------------------------------------------------
% Output
%   lik       : likelihood (T by 1)
%   all_s_up  : particles resampled (updated) (T by N)
%   Neff      : Effective sample size
% -----------------------------------------------------------------------
 rng(111122) % For replicability of results
    
% housekeeping
ne        = size(S2,1);
[~, ns] = size(B);
T         = size(yt,1);
sqrtS2    = R*chol(S2)';

% matrix for store
all_s_up  = zeros(T, ns, N);   % resampled 
lik       = zeros(T,1);
Neff      = zeros(T,1);

% initialization
temp_s = x0;
temp_P = P0;
s_up   = repmat(temp_s, 1, N) + chol(temp_P)'*randn(ns, N);


% Rest of Steps
for tt=1:1:T
    
    yy = yt(tt,:);
    
    if naive == 1 
    % Propagate
    s_fore   = Phi*s_up + sqrtS2*randn(ne, N);
    
    % Un-normalized weights (naive boostrap particle filter)
    perror  = repmat(yy'-A, 1, N) - B*s_fore;
    density = mvnpdf(perror', zeros(1,size(yy,2)), H);
    
    % Sum of density (normalized weights)
    weights = density/mean(density);
    
    
    else
        
    % Propagate  (p.186 of the textbook)  
        
    s_pred           = Phi*s_up;
    P_pred           = R*S2*R';
    P_pred           = nearestSPD(P_pred);
    
    v  = repmat(yy'-A, 1, N) - B*s_pred;
    F  = B*P_pred*B' + H;
    s_upd = s_pred + P_pred * B' * (F\v);
    P_upd = P_pred - P_pred * B' * (F\B) * P_pred;
    P_upd = nearestSPD(P_upd);
    
    s_fore       = s_upd + chol(P_upd)' * randn(ns,N); 
    
    % Weights Calculation 
    perror  = repmat(yy'-A, 1, N) - B*s_fore;
    adjustment        = mvnpdf(s_fore',s_pred',P_pred) ./ mvnpdf(s_fore',s_upd',P_upd);
    density        = mvnpdf(perror',zeros(1,size(yy,2)), H) .* adjustment; 
    
    % Sum of density (normalized weights)
    weights = density/mean(density);   
        
    end
    
   % Effective sample size
   Neff(tt,1) = (N^2)/sum(weights.^2);
   
    
   if Neff(tt,1)>=N/2
       s_up = s_fore; 
   else
       s_up = s_fore(:,randsample(N,N,true,weights)');  % Resampling if ESS falls below a threshold
   end    
    
    % Store results
    lik(tt,1)        = log(mean(density));
    all_s_up(tt,:,:) = s_up;
    
  
    
end


    
  