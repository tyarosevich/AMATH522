function dy = part_c(t,y,p)
    if (0<t) && (t<2)
        z = 1;
    else
        z = 0;
    end    

    dy = zeros(2,1);
    dy(1) = -p(1) / (1+y(2)) + (p(2)*z^p(4))/(1 + z^p(4)) - p(3)*y(1); 
    dy(2) = (p(2)*z^p(4))/(1 + z^p(4));



