
G_PROB = function(pheno_prev,df_SNPs, df_ORs, N){
	#pheno_prev data.frame with on each row of column1 the  disease name and the prevalence of the diseases in column 2 
	#df_SNPs = data.frame wit column 1: PT_id, columns 2 till k: headers are SNP_names, each field contains the number of minor allele a patient has (0,1 or 2) for that particular SNP. 
	#df_ORs column 1: SNP_names, column 2 till k: each columns representing the ORs for the minor allele of SNPs for one disease 


#############CHECK FILES########################

#check df_ORs
if(ncol(df_ORs)<3){
stop("df_ORs should contain at least an ID column and two phenotype column ")	
}

#check df_SNPs
if(ncol(df_SNPs)<2 | nrow(df_SNPs)<2){
stop("df_SNP should contain at least two IDs and at least two genetic variants")}

#check pheno_prev
if(ncol(pheno_prev)!=2 | nrow(pheno_prev)<2){
stop("pheno_prev should contain one columns with at least two phenotypes and one columns with the corresponding population prevalences.")}

#disease names should match in  pheno_prev and df_ORs
phenos = names(df_ORs[,2:ncol(df_ORs)])
pheno_count  = ncol(df_ORs)-1


y = phenos %in% pheno_prev[,1]
if( FALSE  %in% (phenos == pheno_prev[,1])) {
stop("Phenotypes differ between files")} else {
print(paste("Phenotypes match between the files"))}

#remove rows that contain incorrect number of minor alleles
df_SNPs_up = df_SNPs
for (i in 2:(ncol(df_SNPs_up))){
df_SNPs_up= df_SNPs_up [!(df_SNPs_up[,i] >2 | df_SNPs_up[,i] <0),] 
}

ID_included = df_SNPs_up[,1] %in% df_SNPs[,1]
ID_removed = df_SNPs[!ID_included,1]
print("The following ID was removed because the number of minor alleles were less than 0 or more than 2:")
print(ID_removed)

#remove ids with missing data
df_SNPs_noMiss = df_SNPs_up[complete.cases(df_SNPs_up),]
list = df_SNPs_up[,1] %in% df_SNPs_noMiss[,1]
incomplete = df_SNPs_up[!list,1]
if (!is.na(incomplete) ){
print("The following IDs were removed because of incomplete data:")
print(incomplete)}

print(paste("Input contains",pheno_count,"diseases"))

#############MATCH FILES########################
#extract the relevant variants
x = colnames(df_SNPs_noMiss) %in% df_ORs[,1]
df_SNPs_sel = df_SNPs_noMiss[,x]
z = df_ORs[,1] %in% colnames(df_SNPs_noMiss[,-1])
df_ORs_sel = df_ORs[z,]
missing_SNPs = as.character(df_ORs[!z,1])

print("The genotype file does not contain information on the following genetic variants:") 
print(missing_SNPs)
print("These variants were ignored in the calculations.")#**make this an if statement describing wether the missing snps or no missing SNPs


#sort genotype data in similar order as ORs dataframe 
target = df_ORs_sel[,1]
df_SNPs_sel_sort = as.matrix(df_SNPs_sel[,target]) 

#############components of GRS########################

#components of GRS
OR_list =as.matrix(df_ORs_sel[,-1])
prev_list = pheno_prev[,2]

#make the genotype matrices 
geno = as.matrix(df_SNPs_sel_sort)

#############CALCULATIONS########################

#define theta
cal_theta = function(prev, OR_list, geno, theta,N){
 newton <- function(f, tol=1E-12,x0=theta,N=N) {
     h <- 0.00000001 
     i <- 1; x1 <- x0
     p <- numeric(N)
     while (i<=N) {
      print(paste(i,x0,x1))
         df.dx <- (f(x0+h)-f(x0-h))/(2*h)
         x1 <- (x0 - (f(x0)/df.dx))
         p[i] <- x1
         i <- i + 1
         if (abs(x1-x0) < tol) break
         x0 <- x1
     }
     return(p[1:(i-1)])
 }
        f = function(theta){ (mean( 1/ (1 + exp(- (theta + geno  %*% log(OR_list))))) - prev ) ^ 2 }
        q <- newton(f(theta), x0=theta, N=N)
        v=q[length(q)]
        return(v)
}

theta_a = 0.000001 #**** 

theta = rep(theta_a, length(phenos))
print("The theta input values are:") #**remove in final function
print(theta) #**remove in final function
i = 0
theta_run = as.data.frame(matrix(NA,ncol(OR_list),2))
theta_run[,2]= pheno_prev[,1]
 

#calcalute the theta 
while (i < (ncol(OR_list))){
  i = i + 1
  print(i) #**remove in final function
  theta_run[i,1] = cal_theta(as.numeric(as.character(prev_list[i])), OR_list[,i], geno, theta[i],N)
} #note : gives a problem when I use a higher prevalence than 0.001....**

# GRS  resulting in population probability
GRS_pop = data.frame(sapply(1:pheno_count, function(x) 1/(1+exp(-(theta_run[x,1] + geno %*% log(OR_list[,x]))))))
index_pop_probs = paste0("Population_prob_", phenos)
names(GRS_pop) = index_pop_probs
#within cases G-PROB
gg =  sapply(1:pheno_count, function(x,o)  o[,x],  o = GRS_pop[, index_pop_probs])
#row wise summation of gg to obtain he total of each individual
sum_GG =  apply(gg,1, function(x) sum(x))

index_cases_probs = paste0("WithinCases_prob_", phenos)


# gg[,x]/sum_GG for each field
GRS_pop[, index_cases_probs] =  sapply(1:length(index_pop_probs), function(x) gg[,x] / sum_GG)

G_PROB_final = cbind(df_SNPs_noMiss[,1], GRS_pop)
names(G_PROB_final)[1]= names(df_SNPs_noMiss)[1]

print("You've calculated G-PROB for the following disease :")
print(phenos)
print("The final list of genetic variants used to calculate the GRS were:")
print(as.character(target))


return(G_PROB_final)


}

