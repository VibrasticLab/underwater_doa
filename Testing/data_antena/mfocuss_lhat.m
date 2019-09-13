function [X, gamma_ind, gamma_est, count, lambda_history, dmu_history2] = mfocuss_lhat(A, Y, varargin) 
% MFOCUSS algorithm for the MMV model 
%    Mainbody from http://dsp.ucsd.edu/~zhilin/MFOCUSS.m 
% ============================== INPUTS ==============================  
%   A             : N X M dictionary matrix % 
%   Y             : N X L measurement matrix, i.e. Y = Phi * X + V. % 
%  'p'            : p-norm. p lies in [0,1].  Default value: p = 0.8 % 
%  'PRUNE_GAMMA'  : Threshold for prunning small gamma_i. 
%                   In noisy cases, you can set PRUNE_GAMMA = 1e-3 or 1e-4. 
%                   In strongly noisy cases (SNR<=5 dB), suggest to set PRUNE_GAMMA = 0.01; 
%                   Default value: MIN_GAMMA = 1e-4.  %                                    
%  'MAX_ITERS'    : Maximum number of iterations. 
%                   Default value: MAX_ITERS = 800 % 
%  'EPSILON'      : Threshold to stop the whole algorithm.  
%                   Default value: EPSILON = 1e-8 % 
%  'PRINT'        : Display flag. If PRINT = 1: show output; If PRINT = 0: supress output 
%                   Default value: PRINT = 0 % 
% ==============================  OUTPUTS ==============================  
%   X            : Estimated solution matrix(size: M X L) 
%   gamma_ind    : Indexes of nonzero gamma_i 
%   gamma_est    : Final values of all the gamma_i (including zeros). An M X 1 vector 
%   count        : number of iterations used % 
% ============== Examples of Commands =============== 
% [Example 1] 
%   lambda = 1e-3; 
%   [X,gamma_ind,gamma_est,count] ... 
%       = MFOCUSS(Phi,Y,lambda,'p',0.8,'prune_gamma',1e4,'max_iters',500,'epsilon',1e-8,'print',0); % 
% [Example 2] 
%      lambda = 1e-5;   
%      [X,gamma_ind,gamma_est,count] = MFOCUSS(Phi,Y, lambda); % 
% ==============================  Reference ============================= 
%   [1] Cotter, S.F.;   Rao, B.D.;   Kjersti Engan;   Kreutz-Delgado, K.; 
%       Sparse solutions to linear inverse problems with multiple measurement vectors 
%   [2] Zdunek, R.; Cichocki A.;  
%       Improved M-FOCUSS Algorithm with overlapping blocks for locally 
%       smooth sparse signals. 
% ==============================  Author ==============================  
%   Zhilin Zhang (z4zhang@ucsd.edu) 
%   Mainbody was written by David Wipf 
%   Renaming parameters cf ref [1] and adding the Generalized Cross 
%   Valitation Technique to determine lambda cf ref [2] by Michiel van Tent 
%   Beking 
% % ==============================  Version ==============================  
%   1.0 (04/03/2016) %   

%INITIALIZATION FOR THE UPDATE LAMBDA PART 
%Mel=16; %define number of working elements for the golden section search 
fgcv = @GCV %define funciton handle for GS search 
lambda_min=10^-10; %lower boundry value for lambda 
lambda_max=0.99;   %upper boundry value for lambda   

% Dimension of the Problem 
[N M] = size(A);  
[N L] = size(Y);     

% Default Control Parameters  
PRUNE_GAMMA = 1e-4;        % threshold for prunning small gamma_i 
p           = 0.8;         % p-norm 
EPSILON     = 1e-8;        % threshold for stopping iteration.  
MAX_ITERS   = 800;         % maximum iterations 
PRINT       = 0;           % not show progress information   

% get input argument values 
if(mod(length(varargin),2)==1)     
    error('Optional parameters should always go by pairs\n'); 
else
    for i=1:2:(length(varargin)-1)         
        switch lower(varargin{i})             
            case 'p'                 
                p = varargin{i+1};              
            case 'prune_gamma'                 
                PRUNE_GAMMA = varargin{i+1};              
            case 'epsilon'                    
                EPSILON = varargin{i+1};              
            case 'print'                     
                PRINT = varargin{i+1};              
            case 'max_iters'                 
                MAX_ITERS = varargin{i+1};               
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);        
        end
    end
end

if (PRINT) fprintf('\nRunning M-FOCUSS for the MMV Problem...\n'); 
end

% Initializations  
Wk = ones(M,1);         % initialization of gamma_i 
keep_list = [1:M]';      % record the index of nonzero gamma_i 
m = length(keep_list);     % number of nonzero gamma_i 
Xk = zeros(M,L);           % initialization of the solution matrix 
count = 0;                 % record iterations     

% Learning loop  
while (1)       
    % =========== Prune weights as their hyperparameters go to zero ===========
 
    if (min(Wk) < PRUNE_GAMMA )         
        index = find(Wk > PRUNE_GAMMA); %Find the columns in gamma which are <PRUNE_GAMMA         
        Wk = Wk(index);           %only columns with values > PRUNE_GAMMMA remain        
        A = A(:,index);            % corresponding columns in Phi         
        keep_list = keep_list(index);   %store the columns that remain after removing columns with values <PRUNE_GAMMA         
        m = length(Wk);           
        
        if (m == 0)   break; 
        end;    
    end;        
    
    % ====== Compute new weights ======     
    Wk = repmat(sqrt(Wk)',N,1);   %=Wk+1 combined with taking the sqrt of abs(sum(Xk[i])=Ck[i]^2)  
    %Wk = diag(sqrt(Wk));   %=Wk+1 combined with taking the sqrt of abs(sum(Xk[i])=Ck[i]^2)     
    Ak = A.*Wk;  %=AWk+1     
    [U,S,V] = svd(Ak,'econ'); %needed to determine AWk+1 dagger = PhiG dagger       [d1,d2] = size(S);     if (d1 > 1)     diag_S = diag(S);     else            diag_S = S(1);      end;       %NEW PART TO UPDATE LAMBDA PER ITERATION     %find new lambda using Golden section search and Generalized CrossValidation      %see section between equation (21) and (25) of "Improved M-FOCUSS     %algorithm with overlapping blocks for locally smooth sparse signals"     %compute GCV function with V(lambda) eq(22) is done in the (self written) function     %GCV.m Then the minimum lambda_hat is obtained using the Golden Section     %Search function (fminbnd) of Matlab:     lambda_hat = fminbnd(@(lambda_hat) fgcv(lambda_hat,U,Y,diag_S),lambda_min,lambda_max);     lambda_hat=real(lambda_hat); %use real values only and discard complex values          %now continue to calculate X using lambda_hat, according to equation     %(18) sqrt(lambda_hat) + 1e-16     %thikonov regularization with solving 2-norm: PhiGdagger=G'*V*(sigma/(sigma^2+1))*U     U_scaled = U(:,1:min(N,m)).*repmat((diag_S./(diag_S.^2 + sqrt(lambda_hat) + 1e-16))',N,1); %the part (sigma/(sigma^2+lambda))*U = U_scaled     Qk_noy_withWk = Wk'.*(V*U_scaled');  %The part Wk'*V*U_scaled makes the total Wk'*V*(sigma/(sigma^2+1))*U (=Adagger)         %      %=Akdagger     Xk_old = Xk; %store previous solution before writing over the new one     Xk = Qk_noy_withWk*Y;  %=Xk+1 total solution to y=Ax using the tikhnov SVD ==> x=Adagger*y=Xk+1 (see step 3 algortihm.)       % *** Update hyperparameters ***     Wk_old = Wk;    %store the old weighing matrix     ck_nosqrt = sum(abs(Xk).^2,2); %=Ck[i] without taking the sqrt (=Ck[i]^2) 

     [d1,d2] = size(S);     
     if (d1 > 1)     diag_S = diag(S);     
     else            diag_S = S(1);     
     end;       
     
     %NEW PART TO UPDATE LAMBDA PER ITERATION     
     %find new lambda using Golden section search and Generalized CrossValidation      
     %see section between equation (21) and (25) of "Improved M-FOCUSS     
     %algorithm with overlapping blocks for locally smooth sparse signals"     
     %compute GCV function with V(lambda) eq(22) is done in the (self written) function     
     %GCV.m Then the minimum lambda_hat is obtained using the Golden Section     
     %Search function (fminbnd) of Matlab:     
     lambda_hat = fminbnd(@(lambda_hat)fgcv(lambda_hat,U,Y,diag_S),lambda_min,lambda_max);     
     lambda_hat=real(lambda_hat); %use real values only and discard complex values          
     
     %now continue to calculate X using lambda_hat, according to equation     
     %(18) sqrt(lambda_hat) + 1e-16     
     %thikonov regularization with solving 2-norm: PhiGdagger=G'*V*(sigma/(sigma^2+1))*U     
     U_scaled = U(:,1:min(N,m)).*repmat((diag_S./(diag_S.^2 + sqrt(lambda_hat) + 1e-16))',N,1); %the part (sigma/(sigma^2+lambda))*U = U_scaled     
     Qk_noy_withWk = Wk'.*(V*U_scaled');  %The part Wk'*V*U_scaled makes the total Wk'*V*(sigma/(sigma^2+1))*U (=Adagger)         %      
     %=Akdagger    
     Xk_old = Xk; %store previous solution before writing over the new one     
     Xk = Qk_noy_withWk*Y;  %=Xk+1 total solution to y=Ax using the tikhnov SVD ==> x=Adagger*y=Xk+1 (see step 3 algortihm.)       
     
     % *** Update hyperparameters ***     
     Wk_old = Wk;    %store the old weighing matrix     
     ck_nosqrt = sum(abs(Xk).^2,2); %=Ck[i] without taking the sqrt (=Ck[i]^2) 
    
     Wk = (ck_nosqrt/L).^(1-p/2); %=Ck[1[^1-p/2] %create the new weighing vector: SQRT MUST STILL BE TAKEN this hapens in the first step of the new iteration.        
     % ========= Check stopping conditions, etc. =========      
     count = count + 1;     
     lambda_history(count,1)=lambda_hat; %store the lambdas used     
     if (PRINT) disp(['iters: ',num2str(count),'   num coeffs: ',num2str(m),'gamma change: ',num2str(max(abs(Wk - Ck_old)))]); 
     end;     
     
     if (count >= MAX_ITERS) break; 
     end;       
     
     if (size(Xk) == size(Xk_old))        
         dmu = max(max(abs(Xk_old - Xk)));        
         dmu_history2(count,1)=dmu; %store the deltas between the solutions 
         if (dmu < EPSILON)  break;  
         end;    
     end;   
end;     

gamma_ind = sort(keep_list); 
gamma_est = zeros(M,1); 
gamma_est(keep_list,1) = Wk;     

% expand the final solution 
X = zeros(M,L); 
X(keep_list,:) = Xk;      
if (PRINT) fprintf('\nM-FOCUSS finished !\n'); 
end
return;
 
 