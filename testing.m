numsteps = 100; %number of timesteps simulated

A = [0.90, 0.07 ; 
     0.10, 0.93] ;
 
%list of states on this realization. xlist(k)=1 means in state 1 at
%timestep k, etc
states=zeros(1,numsteps);  

%initial state
states(1)=1;

for k=1:numsteps-1

    %uniformly distributed random number - will use for transitions from
    %timestep k to current timestep k+1
    rd=rand ;
    
    if rd < A(1,states(k))  %for transition FROM states(k) to state 1
        states(k+1)=1;
    else
        states(k+1)=2;
    end

end;

%----
figure
set(gca,'FontSize',18)
plot(1:numsteps,states,'.','MarkerSize',20)
xlabel('timestep','FontSize',16)
ylabel('state','FontSize',16)