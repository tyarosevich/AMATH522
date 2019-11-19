function dy = posfeed(t,y,p)
 
   dy = zeros(3,1);
   if t~=10
    dy(1) = 1;
   else:
    dy(10) = 0;
   end
   dy(2) = p(1)*y(1) + p(2)*y(3);
   dy(3) = p(1)*y(1) + p(3)*y(2);



