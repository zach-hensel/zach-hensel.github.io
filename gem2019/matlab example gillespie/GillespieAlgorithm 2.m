%Gillespie Algorithm:

%In this computational exploration we will solve the simple case of a
%promoter expressing mRNA at a constant rate where the mRNA can in turn
%decay.

%Let's start by defining some parameters:
%For a typical mRNA production rate we take the lac promoter in its induced
%state. Here, it has been shown that a new mRNA is produced approximately
%every 3.3 seconds. As a result, our probability of mRNA production per unit
%time (in min^-1) is
r=20;
%In bacteria, mRNA life times are of just a few minutes. Here we take a
%life time of 1.5 minutes measured for the mRNA of a fusion of the yellow
%fluorescent protein to the membrane protein Tsr. This probability of mRNA
%decay per mRNA and unit time (in min^-1) is given by
gamma=1/1.5;

%First, we want to run some simulations where we start with no mRNA
%molecules. As a result, we set the initial condition:
m(1)=0;

%Next, we need to specify how many iterations of the algorithm we want to
%make. For each iteration we will flip a coin.
MaxT=500;

for i=1:MaxT
    %Calculate k1, k2 and k0, the sum of all probabilities per unit time.
    %Notice how k2 actually depends on time and, as a result, k0 does as
    %well.
    k1=r;
    k2=m(i)*gamma;
    k0=k1+k2;
    
    %Next, we want to determine the time to the next reaction, which we
    %will store in the vector tau.
    %In order to determine the time to the next time point we need to draw
    %it from the distribution shown in Figure 19.38(A). As shown in the
    %book one can relate such a distribution to a uniform distribution
    %between 0 and 1 through
    tau(i)=1/k0*log(1/rand);
    %where the function "rand" gives a random number with uniform
    %probability between 0 and 1.
    
    
    %Now, flip a coin to determine which one of the reactions will take
    %place. First, generate a random number between 0 and 1 using the
    %function "rand"
    CoinFlip=rand;
    
    %If the coin landed in the interval [0,k1/k0] then go for reaction 1,
    %corresponding to an mRNA production. If it landed in the interval
    %[k1/k0,1] then do reaction 2, corresponding to an mRNA decay.
    if CoinFlip<=k1/k0
        m(i+1)=m(i)+1;
    else
        m(i+1)=m(i)-1;
    end
end


%Now, in order to plot the mRNA as a function of time,  let's first 
%calculate the time axis. Remember
%that the tau's are just the time increments. Here, we sum over all tau's
%up to each time point.
T(1)=0;
for i=1:length(tau)
    T(i+1)=sum(tau(1:i));
end

%Finally, plot the results of the simulation together with the
%deterministic solution for the mean. Notice that if you run this code
%several times you'll get different traces!
figure(1)
plot(T,m)
hold on
plot(T,r/gamma*(1-exp(-T*gamma)),'-k')
hold off
xlabel('time (min)')
ylabel('number of mRNA molecules')


%Finally, we want to calculate the distribution in steady state. In order
%to do this we will start the simulation with the initial condition
%corresponding to the mean steady state number of mRNA molecules.

%First, clear all variables
clear all

%Set up the parameters once again:
r=20;
gamma=1/1.5;

%Notice the different initial condition
m(1)=round(r/gamma);

%Also, we're going to need more data points in order to calculate the
%moments of the distribution.
MaxT=10000;


for i=1:MaxT
    %Calculate k1, k2 and k0, the sum of all probabilities per unit time.
    %Notice how k2 actually depends on time and, as a result, k0 does as
    %well.
    k1=r;
    k2=m(i)*gamma;
    k0=k1+k2;
    
    %Instead of flipping a coin to calculate the next time point at which
    %there is a transition we just take the average. tau is the time for the 
    %next reaction. If we were solving the full problem we'd flip a coin
    %in order to determine the value of tau and calculate it
    %from: tau(i)=1/k0*log(1/CoinFlip);
    tau(i)=1/k0;
    
    %Now, flip a coin to determine which one of the reactions can take
    %place. First, generate a random number between 0 and 1 using the
    %function "rand"
    CoinFlip=rand;
    
    %If the coin landed in the interval [0,k1/k0] then go for reaction 1,
    %corresponding to an mRNA production. If it landed in the interval
    %[k1/k0,1] then do reaction 2, corresponding to an mRNA decay.
    
    if CoinFlip<=k1/k0
        m(i+1)=m(i)+1;
    else
        m(i+1)=m(i)-1;
    end
end


%One quick way to check our simulation is to see if it follows a Poisson
%distribution. Of course, we can just plot our distribution together with
%the expectation (as we will do later). Still, we can look for one of the
%main features of the Poisson distribution: the variance should be equal to
%the mean.

%In order to calculate the mean value we need to remember that each
%iteration is associated with a different time interval. As a result, the number
%of mRNA molecules stored in m will have existed for different time
%periods, which are specified by tau. This means that when calculating an
%average we need to weight each mRNA occurence in m using the time over
%which it existed.

MeanM=sum(m(1:end-1).*tau)/sum(tau)
SecondMoment=sum((m(1:end-1).^2).*tau)/sum(tau);

VarianceM=SecondMoment-MeanM^2

%We see that indeed the variance and mean are comparable. The more
%iterations we do in our algorithm, the closer they should get to each
%other.

%We end this exploration by computing the probability distribution. For
%each number of mRNA found throughout the simulation we want to ask how
%often we found it. In other words, the probability of finding a certain
%number of mRNA molecules is related to the period of time that number of
%molecules existed in our simulation.

%Let's start by initializing a vector that will keep track of how long each
%number of mRNA molecules was found throughout the simulation/
MaxmRNA = max(m); %This is just the maximum number of mRNA molecules observed.

%p will be the vector were we store the information. It is a vector with
%length MaxmRNA such that each index along the vector corresponds to an
%observed number of mRNA molecules. We start by setting every element in p
%to zero.
p = zeros(1,MaxmRNA);

%Now, let's go through each time step. For each time step we know the
%number of mRNA molecules m(i) and the time that number existed until the
%next reaction, tau(i). We store that information in p.
for i=1:MaxT
   p(m(i)) = p(m(i)) + tau(i);
end
%Finally, we normalize the distribution
p = p/sum(p);


%Now, let's plot a histogram of p together with the expected Poisson
%distribution
figure(1)
bar([1:MaxmRNA],p);
hold on
X = 0:MaxmRNA;
Y = poisspdf(X,r/gamma);
plot(X,Y,'-r');
hold off
xlabel('number of mRNA molecules')
ylabel('probability')
legend('Simulation','Poisson distribution')


%The agreement might not be perfect. At least, it is probably not as good
%as the agreement in the figure associated with this computational
%exploration in the chapter. The main reason for the disagreement is the
%fact that we did not perform a coin flip to determine the values of tau
%(unlike what we actually did to generate the figure). Instead we only took
%its average for each iteration. A second cause of
%disagreement is the number of iterations employed in our simulation. Both
%of these limitations are easily solvable and are, in fact, taken into
%account when this algorithm is used in the context of research.


