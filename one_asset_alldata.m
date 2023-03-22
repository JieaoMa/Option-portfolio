clear
%--------Loading data------------
load('ETF50.mat')



%---------Parameter settings-----------
rf = 0.02/12;      
W0 = 1;        
rho = 0.05;  
L = 5;
weight = [];
return_option = [];
sigma = 0.0733;
mon = [10;12;12;12;12;12;10];
mon_conut = [0;10;12;12;12;12;12];
final_value2 = [];
count = 1;

for year = 1:7
    fix_year = year+14;
    dataname   = strcat('data_50_20',num2str(fix_year));
    datafile   = strcat(dataname,'.mat');
    load(datafile);
    eval(['load FS_predict_50_20',num2str(fix_year),'.mat'])
    final_value = [];
    

    for i = 1:mon(year,1)
        S0_50 = CP(i,1);

        eval(['K50=','K50_',num2str(i)]);
        eval(['OF_c50=','OF_c50_',num2str(i)]);
        eval(['OF_p50=','OF_p50_',num2str(i)]);
    
        n_50 = size(K50,1);      
        v_OF_c = [];v_OF_p = [];
        C1 = []; P1 = []; OF_c = []; OF_p = []; r_call = []; r_call_50 = []; r_put = []; r_put_50 = [];
        C = []; OF = []; rr = [];C50=[];P50=[];CF_50=[];PF_50=[];rc50=[];rp50=[];  S_L= S0_50*0.07;
        
    %---------Calculating the expected value of a contract using forecast prices-----------
    
        S50 = FS_predictive_50(i,1);
        for i_2 = 1:n_50
            if S50 > K50(i_2)
                C50(i_2,1) = S50-K50(i_2,1); 
                P50(i_2,1) = 0; 
            else
                C50(i_2,1) = 0; 
                P50(i_2,1) = K50(i_2,1)-S50; 
            end
        end
    
    
    %---------MC-----------
        scenarios = 100;          
        mu_50 = 0.000288219104703902;
        sig_50 = 0.0148524693516507;
        [S_mc50] = MC(S0_50,mu_50,sig_50);    
        beta = 0.05; 
        alpha = 1.96*sig_50; 
        for j = 1:scenarios
            S_mc1 = S_mc50(j);
            for i_1 = 1:n_50
                if S_mc1 > K50(i_1)
                    r_call_50(j,i_1) = S_mc1-K50(i_1,1); 
                    r_put_50(j,i_1) = 0; 
                else
                    r_call_50(j,i_1) = 0; 
                    r_put_50(j,i_1) = K50(i_1,1)-S_mc1; 
                end
            end 
        end
        p = ones(scenarios,1)/scenarios;
    
    
    %---------solve-----------
        v_OF_c = repmat(OF_c50',100,1);v_OF_p = repmat(OF_p50',100,1);
        C1 = C50; P1 = P50; OF_c = OF_c50; OF_p = OF_p50; r_call = r_call_50-v_OF_c; r_put = r_put_50-v_OF_p;
        C = [C1;P1]; OF = [OF_c;OF_p];rr = [r_call,r_put];
    
        % proposed model
        cvx_begin
            variables xf(1) x(n_50*2) d(scenarios) alp(1) z(n_50*2)
            maximize ((1+rf)*xf+C'*x)
            xf >= 0;
            alp + (1/(1-beta))*p'*d <= rho;
            xf + OF'*x == W0; 
            d >= -rr*x-ones(scenarios,1)*alp;
            d >= 0;
            ones(1,n_50*2)*(2*z-x) <= L*S_L;
            z >= x; 
            z >= 0;
        cvx_end
 
        
        
    %----------out-of-sample returns-----------
        FS = FS_50(i,1);
    
        for ii = 1:n_50
                if FS > K50(ii)
                    CF_50(ii,1) = FS-K50(ii,1); 
                    PF_50(ii,1) = 0; 
                else
                    CF_50(ii,1) = 0;   
                    PF_50(ii,1) = K50(ii,1)-FS; 
                end 
        end
        
        CF = [CF_50;PF_50];
        
        
        eval(['x_',num2str(i),'=x']);
        xx = (abs(x)+x)/2;

        return_option = [return_option;CF-OF];
        
        x_c = x(1:n_50);
        x_p = x(n_50+1:n_50*2);
        final_value(i,1) = (rf+1)*xf+CF_50'*x_c+PF_50'*x_p; 
        eval(['xf',num2str(i),'=xf']);
        eval(['x_c',num2str(i),'=x_c']);
        eval(['x_p',num2str(i),'=x_p']);
    
        count = count+1;
    
    end

    final_value2 = [final_value2;final_value];  %Portfolio value per period
    
end