import numpy as np
import matplotlib.pyplot as plt
import math


def dice(sides): #Generates distribution of single die

    dice=np.array([0.0]*(sides+1));
    for s in range(sides):
        dice[s+1]=1/sides;
    return dice;

def shift(d,F): #Shifts function F units, effectively adding constant damage

    s=np.array([0.0]*F);
    shifted=np.append(s,d);
    return shifted;

def matfold(d1,d2): #Adds distributions
    folded=np.array([0.0]*(d1.size+d2.size))
    for s in range(d1.size):
        b=d1[s]*d2;
        b=shift(b,s);
        b0=np.array([0.0]*(d1.size+d2.size-b.size));
        b=np.append(b,b0);
        folded=folded+b;
    return folded;
    
def activate(attack, AC, d, flat): #Modifies dmg distribution to account for hits, critical hits and misses
    C=1/20; #Critical hit probability (cnst)
    M=(AC-attack-1)/20; #Miss probability
    H=19/20-M; #Ordinary hit probability
    if attack>=AC-1: #Minimal amount of misses
        M=1/20;
        H=18/20;
    if AC-attack-1>=19: #Maximal amount of misses
        M=19/20;
        H=0;

    d0=shift(d,flat); #Damage distribution if hit
    d=matfold(d,d);
    d=shift(d,flat); #Damage distribution if critical hit
    d0=np.append(d0,np.array([0.0]*(d.size-d0.size)));
    dist=H*d0+C*d;
    dist[0]=dist[0]+M; #Missing does 0 damage
    return dist; #Final distribution

def Ndice(N,sides): #Generates distribution for N dice
    d0=dice(sides);
    d=d0;
    if N>=2:
        for s in range(N-1):
            d=matfold(d,d0);
    return d;

def sneak(distribution,sneakdice): #Appends sneak attack to an attack
    b=sneakdice;
    oh=np.array([0.0]);
    dis=np.array([0.0]*sneakdice.size)
    dis[0]=distribution[0]*np.sum(sneakdice);
    for s in range(distribution.size-1):
        b=np.append(oh,b);
        b2=sneakdice[s+1]*b;
        dis=np.append(dis,oh);
        dis=dis+b2;
    dis=dis/np.sum(dis);
    return dis;

def averages(distribution): #Returns the average damage
    tau=np.array([0.0]*distribution.size);
    total=np.sum(distribution);
    for s in range(distribution.size):
        tau[s]=s*distribution[s];
    average=np.sum(tau)/total;
    return average;

def spread(distribution,average): #Returns the standard deviation
    tau=np.array([0.0]*distribution.size);
    total=np.sum(distribution);
    for s in range(distribution.size):
        tau[s]=(s-average)*(s-average)*distribution[s];
    spr=(np.sum(tau))/total;
    spr=math.sqrt(spr)
    return spr;

def survprops(distribution,Edistribution,rounds,HP,EHP):
    #Returns probability for victory, assuming constant dmg distribution
    d=distribution;
    Ed=Edistribution;
    surv=1.0;
    Esurv=1.0;

    if d.size>=HP:
        surv=0.0;
        for i in range(HP):
            surv=surv+d[i];#probability that enemy survives
    if Ed.size>=EHP:
        Esurv=0.0;
        for i in range(EHP):
            Esurv=Esurv+Ed[i];#probability that I survive enemy atk
    M=np.matrix([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
    H=np.matrix([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
    M[0,1]=1-surv;
    M[1,1]=surv;
    H[2,1]=1-Esurv;
    H[1,1]=Esurv;
    P0=np.matmul(H,M,out=None);
    for s in range(rounds-1):
        surv=1.0;
        Esurv=1.0;
        d=matfold(d,distribution);
        Ed=matfold(Ed,Edistribution);
        if d.size>=HP:
            surv=0.0;
            for i in range(HP-1):
                surv=surv+d[i];#probability that enemy survives
        if Ed.size>=EHP:
            Esurv=0.0;
            for i in range(EHP-1):
                Esurv=Esurv+Ed[i];#probability that I survive enemy atk
        M=np.matrix([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        H=np.matrix([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        M[0,1]=1-surv;
        M[1,1]=surv;
        H[2,1]=1-Esurv;
        H[1,1]=Esurv;
        P1=np.matmul(H,M,out=None);
        P0=np.matmul(P1,P0,out=None);
    
    return P0;




        


#atk=#input("Input attack score: ");
#atk=int(atk);
#AC=input("Input enemy AC: ");
#AC=int(AC);
#sides=input("Input number of sides on dmg dice: ");
#sides=int(sides);
#N=input("Input number of this dmg dice: ");
#N=int(N);
#flat=input("Input flat damage (no dice): ");
#flat=int(flat);
#natk=input("How many attacks do you do with these? (at least 1)");
#natk=int(natk);
#atk2=input("Do you have a second type of dmg dice? (1 for yes, 0 for no)");
#atk2=int(atk2);


#User's parameters (this guy goes first)
atk=10; #attack modifier
AC=21; #Armor class
sides=8; #Sides of dice
N=1; #Number of dice/attack
flat=7; #number added to each attack if it hits
natk=4; #Number of identical attacks
HP=84; #HP of the enemy (what he has to kill)

atk2=0; #Whether or not he has a second type of attack

sides2=0; #sides on second attack's dice
N2=0; #Number of dice on second attack


#Similarly for the opponent
Eatk=10;
EAC=21;
Esides=8;
EN=1;
Eflat=7;
Enatk=4;
EHP=84;


#------------Ignore this, this is for other cases than the current one---     
if atk2==1:
    sides2=input("Input number of sides on dmg dice: ");
    sides2=int(sides2);
    N2=input("Input number of this dmg dice: ");
    N2=int(N2);
     
#sneak=input("Do you have sneak attack, or similar? (1 for yes, 0 for no)");
#sneak=int(sneak);
sneak=0;
if sneak==1:
    sidesneak=input("How many sides on extra dice? ");
    sidesneak=int(sidesneak);
    Nsneak=input("How many extra dice? ");
    Nsneak=int(Nsneak);
#-----------------------------------------------------------------------


Ed=Ndice(EN,Esides);
d=Ndice(N,sides);


#----------Ignore this too---------------------------------------------
if atk2==1:
    if sides2>=1:
        if N2>=1:
            d2=Ndice(N2,sides);
            d=matfold(d,d2);
#--------------------------------------------------------------------            
Edistribution=activate(Eatk,EAC,Ed,Eflat);
distribution=activate(atk,AC,d,flat);


if natk>=2:
    dist2=distribution;
    for s in range(natk-1):
        distribution=matfold(dist2,distribution);
if Enatk>=2:
    Edist2=Edistribution;
    for s in range(Enatk-1):
        Edistribution=matfold(Edist2,Edistribution);

#----------------This isn't used at the moment---------------------
if sneak==1:
    sneakdice=Ndice(Nsneak,sidesneak);
    distribution=sneak(distribution,sneakdice);





damage=np.array([0.0]*distribution.size);
for s in range(distribution.size):
    damage[s]=s;

#-----------------------------------------------------------------

z=np.array([0.0]*10); #Turn-counter
w=np.array([0.0]*10); #prepares victory probability for guy 1
L=np.array([0.0]*10)  #prepares loss probability for guy 1
c=np.array([0.0]*10)  #Prepares probability of continued battle
for rounds in range(10):
    P0=survprops(distribution,Edistribution,rounds+2,HP,EHP);
    z[rounds]=rounds+2 #Sets turn counter
    w[rounds]=P0[0,1] #writes win probability
    L[rounds]=P0[2,1] #Writes loss probability
    c[rounds]=P0[1,1] #Writes continued battle probability

#Printing final results:
print(w[9])
print(c[9])
print(L[9])


#Plotting as a function of turns:    
plt.plot(z,w,label="victory")
plt.plot(z,L,label="Failure")
plt.plot(z,c,label="battle goes on")
plt.xlabel('Turn')
plt.ylabel('Probability')
plt.title('monk mirror match')
plt.legend()
plt.show()

'''
average=averages(distribution);
spr=spread(distribution,average);
str1="average = "+str(average)+" (+/-) "+str(spr);
plt.plot(damage,distribution)
plt.ylabel('Probability')
plt.xlabel('Damage')
plt.title(str1)
plt.show()
'''

