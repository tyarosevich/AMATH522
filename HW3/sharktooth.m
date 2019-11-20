function dy = sharktooth(t,y,p)
    % Condition on the threshhold equation that immediately starts 
    % production of Z1.
    if y(1) > 0
        z1=1;
    else
        z1=0;
    end
    % Condition on threshhold equation to start production of Y1 
    % a little later
    if y(1) > 1.5
        z2 = 1;
    else
        z2=0;
    end
    % If condition to start X2 as soon as X1 and Y1 are being made
    if y(3) > 0
        z3 = 1;
    else 
        z3=0;
    end
    % Threshhold condition for Y2
    if y(4) > 1.5
        z4 = 1;
    else
        z4=0;
    end
    dy = zeros(7,1);
    % Threshhold equation to represent X1 activating 
    % Z1 and then afterward Y1
    dy(1) = 1 - .3*y(1);
    % Corresponds to Z1
    dy(2) = p(1) / (.0001 + y(3)^2) - y(2);
    % Corresponds to Y1
    dy(3) = -z2*p(2) * (y(3) - y(2)) ;
    % Threshhold equation to represent X2
    dy(4) = z3 - .3*y(4);
    % Corresponds to Z2
    dy(5) = z3*p(1) / (.0001 + y(6)^2) - y(5);
    % Corresponds to Y2
    dy(6) = -z4*p(2) * (y(6) - y(5));
    % Corresponds to Z3
    dy(7) = 8.8*z4 - y(7);


%Maybe try having the threshhold turn z1 off as it turns z2 on? I dunno.
%Curernt form gives the tooth.