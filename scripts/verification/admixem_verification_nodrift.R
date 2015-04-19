# calculate change in frequency of a single locus under selection or drift# define the selection coefficient (0-1)
s<-0.1#starting allele frequency for the simulation
q<-0.5
#dominance (0-1)h<-0#set number of generations to simulate 
num_gen<-100

#initialize frequency arrayfreq<-{}

#run selection for num_gen
for (x in 1:num_gen){q = (q*(1-q)*(1+h*s) + q^2*(1+s))/((1-q)^2 + 2*q*(1-q)*(1+h*s) + q^2*(1+s))freq<-c(freq,q)} 