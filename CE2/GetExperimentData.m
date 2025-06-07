function [y,u,Ts] = GetExperimentData(path_)

    run(path_)

    ind = find(Y3,1,'first')
    Y30 = Y3(ind:end);
    [p,nr] = seqperiod(Y30)


    periods = 3;
    u = Y30(1,1:p*periods)';
    y = Y2(1,ind:ind+p*periods-1)';

    %Ts =  T(ind+1)-T(ind);
    Ts = 0.01;
end
