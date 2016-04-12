simple_positive_feedback <- function(k_1,k_2,k_3,maxtime,timestep,numberofrealisations,initial_condition=0){

output_GFP = matrix(0,nrow=(maxtime/timestep)+1,ncol=numberofrealisations);

for (realisation in seq(1,numberofrealisations,by=1)){

#print(realisation)
GFP= initial_condition;

a = rep(0,3);
a[1] = k_1;
a[2] = k_2*GFP;
a[3] = k_3*GFP;
asum = sum(a);

time=0.0;

count = 1;

for (ts in seq(0,maxtime,by=timestep)) {

while (time < ts){

tau = log(1/runif(1))/asum;

time = time + tau;

#DETERMINE WHICH REACTION WILL OCCUR
psm = 0;
drxn = asum*runif(1);
j = 0;

while (psm < drxn){
j = j + 1;
psm = psm+ a[j];
}

#UPDATE SYSTEM CONFIGURATION AND FURTHER PARAMETERS
if (j == 1){
GFP = GFP + 1;
a[2] = a[2] + k_2;
a[3] = a[3] + k_3;
asum = asum + k_2 + k_3;}
else if (j == 2){
GFP = GFP - 1;

#HGF diffusion
a[2] = a[2] - k_2;
a[3] = a[3] - k_3;
asum = asum - k_2 - k_3;
}
else if (j == 3){
GFP = GFP + 1;

a[2] = a[2] + k_2;
a[3] = a[3] + k_3;
asum = asum + k_2 + k_3;
}
} # End of integration loop
output_GFP[count,realisation] = GFP;  
count = count + 1;
}
}
return(output_GFP)
}