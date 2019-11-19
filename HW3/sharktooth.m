function dy = sharktooth(t,y,p)
    if y(1) > 0
        z1=1;
    else
        z1=0;
    end
    if y(1) > 1.5
        z2 = 1;
    else
        z2=0;
    end
    dy = zeros(3,1);
    dy(1) = 1 - .3*y(1);
    dy(2) = p(1) / (.01 + y(3)^.2) - y(2); 
    dy(3) = -z2*p(2) * (y(3) - y(2)) ;
    
%Maybe try having the threshhold turn z1 off as it turns z2 on? I dunno.
%Curernt form gives the tooth.