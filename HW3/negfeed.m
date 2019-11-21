function dy = negfeed(t,y,p)
 
    if (2<t) && (t<4)
        z = 1;
    else
        z = 0;
    end
    dy = zeros(2,1);
    dy(1) = ( p(1) * (z + y(1))^p(4) ) /(1 + (y(1) + z)^p(4))*(p(1) / (1+y(2)^2)) - p(3)*y(1); 
    dy(2) = - (y(1) + z)* (y(2))*(p(1) / (1+(z +y(1))^2));



