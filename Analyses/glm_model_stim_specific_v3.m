function res=glm_model_stim_specific_v3(x)
%Need to do this so the code works properly
global f_spikes f_spikes_evenodd t c train_test_flag with_T lb ub trial_length ut st

%Transforms data
[y] = vec2struct(lb,x);

%This if/ifelse is set up to determine which equation the spikes are being
%Fit with
if with_T==1 % Fits Exponentially Modified Gaussian
    T=y.peak*((1/2)*exp((y.tau/2).*(2*y.mu + y.tau*y.sig^2-2*t)).*erfc((y.mu + y.tau*y.sig^2 - t)/(sqrt(2)*y.sig)));
    
elseif with_T==2 % Fits Guassian
    T=y.peak*exp((-(t - y.mu).^2)/(2*y.sig^2));
end

%compute log likelihood or produce time series of the fit
if train_test_flag==1 %if fitting compute log likelihood
    if with_T==1 ||  with_T==2 ||  with_T==3 || with_T==4 || with_T==5 || with_T==6
        res=sum(f_spikes.*(-log(y.o+T))+...
            (1-f_spikes).*(-log(1-(y.o+T))));
    elseif with_T==0
        res=sum(f_spikes.*(-log(y.o))+...
            (1-f_spikes).*(-log(1-y.o)));
    end
elseif train_test_flag>=2 %compute LL per individual trial
    tl=trial_length; %trial length in samples
    no_trials=length(f_spikes_evenodd)/tl;
    n_split=train_test_flag; %in how many blocks to split the data
    tln=tl*floor(no_trials/n_split); %lenght in ms of each block
    if with_T==1
        for trial=1:n_split%length(P)/tln
            res(trial)=sum(f_spikes_evenodd(tln*(trial-1)+1:tln*(trial-1)+tln).*(-log((y.o+T(tln*(trial-1)+1:tln*(trial-1)+tln))))+...
                (1-f_spikes_evenodd(tln*(trial-1)+1:tln*(trial-1)+tln)).*(-log(1-(y.o+T(tln*(trial-1)+1:tln*(trial-1)+tln)))));
        end
    elseif with_T==0
        for trial=1:n_split%length(P)/tln
            res(trial)=sum(f_spikes_evenodd(tln*(trial-1)+1:tln*(trial-1)+tln).*(-log((y.o)))+...
                (1-f_spikes_evenodd(tln*(trial-1)+1:tln*(trial-1)+tln)).*(-log(1-(y.o))));
        end
    end
elseif train_test_flag==0 %if done fitting produce time series of the fit
    if with_T==1 || with_T==2 || with_T==3 || with_T==4 || with_T==5 || with_T==6
        res=y(1).o+T;
    elseif with_T==0
        res=y(1).o+zeros(length(f_spikes),1);
    end
end

