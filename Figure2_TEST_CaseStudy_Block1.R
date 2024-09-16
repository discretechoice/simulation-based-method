################################################################################
#
#
#  This script implements the TEST and reproduce FIGURE 2 presented in the paper by 
#  Mariel and Artabe titled 
#
# "A simulation-based method for the identification of lexicographic, non-trading, 
#  or non-compensatory behaviour in stated choice studies" 
#
#  Date: 2024-09-09
#
#  TESTED PATTERNS IN BLOCK 1:
#                  (1,1,1,1,1)
#                  (2,2,2,2,2)
#                  (3,3,3,3,3)
#                  (4,4,4,4,4)
#                  (2,3,1,1,1)
#                  (2,3,1,1,3)
#                  (2,3,2,1,1)
#                  (2,3,2,1,3)
#
################################################################################

# Preliminaries
rm(list = ls(all = TRUE))             # Clear workspace
Sys.setenv(LANG = "en")               # Set language
setwd(".")                            # Set working directory             
set.seed(12345)

# Loading R packages
library(apollo)
library(evd)                          # package for functions for Extreme Value 
                                      # distributions (needed for "rgumbel")
library(ggplot2)
library(patchwork)
library(plotrix)

# Reading the design matrix
design.5.rows <- read.table("data.block1.txt", header=TRUE) 


# Setting parameters
n.individuals         <- 242                  # Number of individuals
n.choices             <- 5                    # Number of choice occasions of one individual
n.blocks		          <- 1 						        # Number of clocks in each design
n.filas 		          <- n.choices*n.blocks		# Number of rows in each design 
n.designs		          <- 1                    # Number of designs in the data
n.ite.R               <- 700                  # Number of iterations to simulate R first stage draws for Beta
n.ite.B               <- 700                  # Number of iterations to simulate B second stage draws for the error term
n.iteraciones         <- n.ite.R * n.ite.B    # Number of total iterations
n.parametros.estimados<- 7                    # Number of parameters of the model

#######################
##
##   OPENS THE DESIGN AND CREATES THE DATA FRAME WITH "n.individuals*n.choices" ROWS
##
#######################

k                     <-1
filamin               <- n.filas*(k-1)+1
filamax               <- n.filas*k 
design.matrix         <- design.5.rows[filamin:filamax,]
n.rows.design.matrix  <- length(design.matrix[,1])             # Number of rows in design matrix

# Replicates the design matrix to create a data frame with "n.individuals * n.choices" rows
data.RPL      <- do.call("rbind", replicate((n.individuals*n.choices)/n.rows.design.matrix, design.matrix, simplify = FALSE))
n.rows        <- length(data.RPL[,1])        # Number of rows 
cat("Design matrix has been replicated","\n",
    "Number of individuals: ", n.individuals       ,"\n",
    "Number of rows:        ", length(data.RPL[,1]),"\n")

# Creates new variables
data.RPL$choice.card   <- rep(seq(1,n.choices)    , n.individuals)     # ID of choice card 
data.RPL$id.individual <- rep(seq(1,n.individuals),each=n.choices)     # ID of individual


################
##
#  STEP 1: COMPUTES OBSERVED RELATIVE FREQUENCY OF EACH TESTED CHOICE PATTERN
##
################

OBSERV.vector.pattern.11111 <- 20/ n.individuals   # Tested pattern (1,1,1,1,1)  
OBSERV.vector.pattern.22222 <- 11/ n.individuals   # Tested pattern (2,2,2,2,2)
OBSERV.vector.pattern.33333 <- 10/ n.individuals   # Tested pattern (3,3,3,3,3)
OBSERV.vector.pattern.44444 <- 25/ n.individuals   # Tested pattern (4,4,4,4,4)
OBSERV.vector.pattern.23111 <- 6 / n.individuals   # Tested pattern (2,3,1,1,1)
OBSERV.vector.pattern.23113 <- 1 / n.individuals   # Tested pattern (2,3,1,1,3)
OBSERV.vector.pattern.23211 <- 2 / n.individuals   # Tested pattern (2,3,2,1,1)
OBSERV.vector.pattern.23213 <- 0 / n.individuals   # Tested pattern (2,3,2,1,3)

################################################
##
##  STEP 2 : THE ESTIMATED PARAMETERS THROUGH THE USE OF CORRELATED-RPL
##
################################################

mu_asc1                       <-  7.584976
mu_asc2                       <-  6.856404
mu_asc3                       <-  7.045594
mu_infmortality               <- -0.320598
mu_reducedvis                 <- -0.053152
mu_morbidity                  <- -0.338522
mu_cost                       <- -0.067447 

chol_asc1                     <- -8.301113
chol_asc2                     <-  2.061707
chol_asc3                     <-  1.533441
chol_asc1_asc2                <- -8.389787
chol_asc1_asc3                <- -8.137979
chol_asc1_infmortality        <- -0.098207
chol_asc2_asc3                <-  1.368488
chol_asc2_infmortality        <-  0.144170
chol_asc3_infmortality        <- -0.008104
chol_infmortality             <-  0.370757
chol_reducedvis               <-  0.008898
chol_morbidity                <- -0.015961
chol_asc1_reducedvis          <- -0.037465
chol_asc1_morbidity           <-  0.042747
chol_asc1_cost                <- -0.054779
chol_asc2_reducedvis          <-  0.022215
chol_asc2_morbidity           <-  0.214573
chol_asc2_cost                <-  0.016017
chol_asc3_reducedvis          <-  0.006635
chol_asc3_morbidity           <- -0.135766
chol_asc3_cost                <-  0.005929
chol_infmortality_reducedvis  <-  0.071227
chol_infmortality_morbidity   <-  0.428432
chol_infmortality_cost        <-  0.016922
chol_reducedvis_morbidity     <-  0.072000
chol_reducedvis_cost          <-  0.024227
chol_morbidity_cost           <- -0.007731
chol_cost                     <- -0.024953  

choleski.cov <- t(matrix (c(chol_asc1             , 0                      , 0                     , 0                           , 0                        ,0                  ,0      ,
                            chol_asc1_asc2        , chol_asc2              , 0                     , 0                           , 0                        ,0                  ,0      ,
                            chol_asc1_asc3        , chol_asc2_asc3         , chol_asc3             , 0                           , 0                        ,0                  ,0      ,
                            chol_asc1_infmortality, chol_asc2_infmortality , chol_asc3_infmortality, chol_infmortality           , 0                        ,0                  ,0      ,
                            chol_asc1_reducedvis  , chol_asc2_reducedvis   , chol_asc3_reducedvis  , chol_infmortality_reducedvis, chol_reducedvis          ,0                  ,0      ,
                            chol_asc1_morbidity   , chol_asc2_morbidity    , chol_asc3_morbidity   , chol_infmortality_morbidity , chol_reducedvis_morbidity,chol_morbidity     ,0      ,
                            chol_asc1_cost        , chol_asc2_cost         , chol_asc3_cost        , chol_infmortality_cost      , chol_reducedvis_cost     ,chol_morbidity_cost,chol_cost )
                          , 7))


######################################
##
#   STEP 3: SIMULATES HYPOTHETICAL CHOICES REPEATEDLY USING THE ATTRIBUTES VALUES, 
#           NORMALLY DISTRIBUTED DRAWS FOR THE BETAS (based on the estimated parameters: means, standard deviations and correlations) 
#           AND GUMBEL DISTRIBUTED RANDOM DRAWS FOR THE EPSILON 
##
######################################

# Creates vectors to store the count of how many times each tested choice pattern appears in each iteration
vector.pattern.11111 <- rep(0,n.iteraciones)  # Tested pattern (1,1,1,1,1)
vector.pattern.22222 <- rep(0,n.iteraciones)  # Tested pattern (2,2,2,2,2)
vector.pattern.33333 <- rep(0,n.iteraciones)  # Tested pattern (3,3,3,3,3)
vector.pattern.44444 <- rep(0,n.iteraciones)  # Tested pattern (4,4,4,4,4)
vector.pattern.23111 <- rep(0,n.iteraciones)  # Tested pattern (2,3,1,1,1)
vector.pattern.23113 <- rep(0,n.iteraciones)  # Tested pattern (2,3,1,1,3)
vector.pattern.23211 <- rep(0,n.iteraciones)  # Tested pattern (2,3,2,1,1)
vector.pattern.23213 <- rep(0,n.iteraciones)  # Tested pattern (2,3,2,1,3)


simu.normal.0.1.ite <- cbind(rnorm(n.individuals*n.ite.R),
                             rnorm(n.individuals*n.ite.R),
                             rnorm(n.individuals*n.ite.R),
                             rnorm(n.individuals*n.ite.R),
                             rnorm(n.individuals*n.ite.R),
                             rnorm(n.individuals*n.ite.R),
                             rnorm(n.individuals*n.ite.R))

simu.beta.attr    <- cbind(rep(mu_asc1        ,n.individuals),
                           rep(mu_asc2        ,n.individuals),
                           rep(mu_asc3        ,n.individuals),
                           rep(mu_infmortality,n.individuals),
                           rep(mu_reducedvis  ,n.individuals),
                           rep(mu_morbidity   ,n.individuals),
                           rep(mu_cost        ,n.individuals))

###
# OBTAINs THE EMPIRICAL DISTRIBUTION OF OBSERVING EACH OF THE TESTED PATTERNS
###

##
# OPENS A LOOP TO REPEAT "n.ite.R" TIMES  
###

for(ind.R in 1:n.ite.R){
  
  simu.normal.0.1   <- simu.normal.0.1.ite[(((n.ite.R-1)*n.individuals)+1):(n.ite.R*n.individuals),]
 
  # Generates random MVN (beta.attr, cov.matrix)
  simu.random.attr         <- simu.normal.0.1 %*% t(choleski.cov) + simu.beta.attr
  simu.random.asc1         <- simu.random.attr[,1]
  simu.random.asc2         <- simu.random.attr[,2]
  simu.random.asc3         <- simu.random.attr[,3]
  simu.random.infmortality <- simu.random.attr[,4]
  simu.random.reducedvis   <- simu.random.attr[,5]
  simu.random.morbidity    <- simu.random.attr[,6]
  simu.random.cost         <- simu.random.attr[,7]
  
  # each individual chooses in "n.choices" choice cards, so the parameter of each individual must be repeated "n.choices" times
  data.RPL$rand.asc1         <- rep(simu.random.asc1        , each = n.choices)
  data.RPL$rand.asc2         <- rep(simu.random.asc2        , each = n.choices)
  data.RPL$rand.asc3         <- rep(simu.random.asc3        , each = n.choices)
  data.RPL$rand.infmortality <- rep(simu.random.infmortality, each = n.choices)
  data.RPL$rand.reducedvis   <- rep(simu.random.reducedvis  , each = n.choices)
  data.RPL$rand.morbidity    <- rep(simu.random.morbidity   , each = n.choices)
  data.RPL$rand.cost         <- rep(simu.random.cost        , each = n.choices)
  
  for(ind.B in 1:n.ite.B){
  
      # Random errors generation Gumbel(0,1)
      data.RPL$error1 <- rgumbel(n.individuals * n.choices ,loc=0, scale=1)
      data.RPL$error2 <- rgumbel(n.individuals * n.choices ,loc=0, scale=1)
      data.RPL$error3 <- rgumbel(n.individuals * n.choices ,loc=0, scale=1)
      data.RPL$error4 <- rgumbel(n.individuals * n.choices ,loc=0, scale=1)
      
      
      ## Utility and choice generation
      data.RPL$utility1.Step3 <- (                                       # Creating utility 1 
                                     data.RPL$rand.asc1
                                   + data.RPL$rand.infmortality * data.RPL$infmortality.1 
                                   + data.RPL$rand.reducedvis   * data.RPL$reducedvis.1
                                   + data.RPL$rand.morbidity    * data.RPL$morbidity.1
                                   + data.RPL$rand.cost         * data.RPL$cost.1
                                   + data.RPL$error1
                                 )
      
      data.RPL$utility2.Step3 <- (                                       # Creating utility 2 
                                     data.RPL$rand.asc2
                                   + data.RPL$rand.infmortality * data.RPL$infmortality.2 
                                   + data.RPL$rand.reducedvis   * data.RPL$reducedvis.2
                                   + data.RPL$rand.morbidity    * data.RPL$morbidity.2
                                   + data.RPL$rand.cost         * data.RPL$cost.2
                                   + data.RPL$error2
                                 )
      
      data.RPL$utility3.Step3 <- (                                       # Creating utility 3 
                                     data.RPL$rand.asc3
                                   + data.RPL$rand.infmortality * data.RPL$infmortality.3
                                   + data.RPL$rand.reducedvis   * data.RPL$reducedvis.3
                                   + data.RPL$rand.morbidity    * data.RPL$morbidity.3
                                   + data.RPL$rand.cost         * data.RPL$cost.3
                                   + data.RPL$error3
                                 )
      
      data.RPL$utility4.Step3 <- (                                       # Creating utility 4    
                                     data.RPL$rand.infmortality * data.RPL$infmortality.4 
                                   + data.RPL$rand.reducedvis   * data.RPL$reducedvis.4
                                   + data.RPL$rand.morbidity    * data.RPL$morbidity.4
                                   + data.RPL$rand.cost         * data.RPL$cost.4
                                   + data.RPL$error4
                                 )
      
      data.RPL$choice.Step3   <- apply(cbind(data.RPL$utility1.Step3, data.RPL$utility2.Step3, data.RPL$utility3.Step3, data.RPL$utility4.Step3),1,which.max)    # Generating choice
      
      data.RPL$patternStep3 <- 0
      
      ## Identifies the choice pattern of each individual  
      for (ind in (1:n.individuals)){
        rows.individuo <- which(data.RPL$id.individual == ind)  
        if(all(data.RPL$choice.Step3[rows.individuo] == c(1,1,1,1,1))){data.RPL$patternStep3[rows.individuo] <- 11111}
        if(all(data.RPL$choice.Step3[rows.individuo] == c(2,2,2,2,2))){data.RPL$patternStep3[rows.individuo] <- 22222}
        if(all(data.RPL$choice.Step3[rows.individuo] == c(3,3,3,3,3))){data.RPL$patternStep3[rows.individuo] <- 33333}
        if(all(data.RPL$choice.Step3[rows.individuo] == c(4,4,4,4,4))){data.RPL$patternStep3[rows.individuo] <- 44444}
        
        if(all(data.RPL$choice.Step3[rows.individuo] == c(2,3,1,1,1))){data.RPL$patternStep3[rows.individuo] <- 23111}
        if(all(data.RPL$choice.Step3[rows.individuo] == c(2,3,1,1,3))){data.RPL$patternStep3[rows.individuo] <- 23113}
        if(all(data.RPL$choice.Step3[rows.individuo] == c(2,3,2,1,1))){data.RPL$patternStep3[rows.individuo] <- 23211}
        if(all(data.RPL$choice.Step3[rows.individuo] == c(2,3,2,1,3))){data.RPL$patternStep3[rows.individuo] <- 23213}
      }                                   
      
      # Counts how many times each tested choice pattern appears and stores the result in the corresponding vector. 
      vector.pattern.11111[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==11111))/n.choices
      vector.pattern.22222[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==22222))/n.choices
      vector.pattern.33333[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==33333))/n.choices
      vector.pattern.44444[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==44444))/n.choices
      vector.pattern.23111[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==23111))/n.choices
      vector.pattern.23113[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==23113))/n.choices
      vector.pattern.23211[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==23211))/n.choices
      vector.pattern.23213[(((ind.R-1)*n.ite.B)+ind.B)] <-  length (which(data.RPL$patternStep3==23213))/n.choices
      
      
      #####
      # CLOSES A LOOP TO REPEAT "n.ite.B" TIMES FOR THE r-th BETA  
      #####
  }
  
  #####
  # CLOSES A LOOP TO REPEAT "n.ite.R" TIMES 
  #
  # OBTAINs THE EMPIRICAL DISTRIBUTION OF OBSERVING EACH OF THE TESTED PATTERNS
  #####
  
}

###
##   CALCULATES THE EMPIRICAL DISTRIBUTION OF THE PROBABILITY OF OBSERVING EACH OF THE TESTED PATTERNS
###
PROB.vector.pattern.11111  <- vector.pattern.11111/ n.individuals 
PROB.vector.pattern.22222  <- vector.pattern.22222/ n.individuals
PROB.vector.pattern.33333  <- vector.pattern.33333/ n.individuals
PROB.vector.pattern.44444  <- vector.pattern.44444/ n.individuals
PROB.vector.pattern.23111  <- vector.pattern.23111/ n.individuals
PROB.vector.pattern.23113  <- vector.pattern.23113/ n.individuals
PROB.vector.pattern.23211  <- vector.pattern.23211/ n.individuals
PROB.vector.pattern.23213  <- vector.pattern.23213/ n.individuals


################################################
##
#   STEP 4 : TESTS WHETHER THE OBSERVED FREQUENCY OF THE PATTERN IS LESS THAN OR EQUAL TO THE SIMULATED FREQUENCY? 
##
################################################

alpha <- 0.05      # significance level
alpha1 <- 1-alpha

# Pattern (1,1,1,1,1)
if (quantile(PROB.vector.pattern.11111,0)<=OBSERV.vector.pattern.11111 & quantile(PROB.vector.pattern.11111,alpha1)>=OBSERV.vector.pattern.11111){
  TEST.PROB.pattern11111 <- "YES"
}else{
  TEST.PROB.pattern11111 <- "NO"
}

#Pattern (2,2,2,2,2)
if (quantile(PROB.vector.pattern.22222,0)<=OBSERV.vector.pattern.22222 & quantile(PROB.vector.pattern.22222,alpha1)>=OBSERV.vector.pattern.22222){
  TEST.PROB.pattern22222 <- "YES"
}else{
  TEST.PROB.pattern22222 <- "NO"
}

#Pattern (3,3,3,3,3)
if (quantile(PROB.vector.pattern.33333,0)<=OBSERV.vector.pattern.33333 & quantile(PROB.vector.pattern.33333,alpha1)>=OBSERV.vector.pattern.33333){
  TEST.PROB.pattern33333 <- "YES"
}else{
  TEST.PROB.pattern33333 <- "NO"
}

#Pattern (4,4,4,4,4)
if (quantile(PROB.vector.pattern.44444,0)<=OBSERV.vector.pattern.44444 & quantile(PROB.vector.pattern.44444,alpha1)>=OBSERV.vector.pattern.44444){
  TEST.PROB.pattern44444 <- "YES"
}else{
  TEST.PROB.pattern44444 <- "NO"
}

#Pattern (2,3,1,1,1)
if (quantile(PROB.vector.pattern.23111,0)<=OBSERV.vector.pattern.23111 & quantile(PROB.vector.pattern.23111,alpha1)>=OBSERV.vector.pattern.23111){
  TEST.PROB.pattern23111 <- "YES"
}else{
  TEST.PROB.pattern23111 <- "NO"
}

#Pattern (2,3,1,1,3)
if (quantile(PROB.vector.pattern.23113,0)<=OBSERV.vector.pattern.23113 & quantile(PROB.vector.pattern.23113,alpha1)>=OBSERV.vector.pattern.23113){
  TEST.PROB.pattern23113 <- "YES"
}else{
  TEST.PROB.pattern23113 <- "NO"
}

#Pattern (2,3,2,1,1)
if (quantile(PROB.vector.pattern.23211,0)<=OBSERV.vector.pattern.23211 & quantile(PROB.vector.pattern.23211,alpha1)>=OBSERV.vector.pattern.23211){
  TEST.PROB.pattern23211 <- "YES"
}else{
  TEST.PROB.pattern23211 <- "NO"
}

#Pattern (2,3,2,1,3)
if (quantile(PROB.vector.pattern.23213,0)<=OBSERV.vector.pattern.23213 & quantile(PROB.vector.pattern.23213,alpha1)>=OBSERV.vector.pattern.23213){
  TEST.PROB.pattern23213 <- "YES"
}else{
  TEST.PROB.pattern23213 <- "NO"
}

# Saves the results of the test for each of the patter
RESULTS <- matrix(NA,8,3 )
rownames(RESULTS) <- c(
  "pattern.11111",
  "pattern.22222",
  "pattern.33333",
  "pattern.44444",
  "pattern.23111",
  "pattern.23113",
  "pattern.23211",
  "pattern.23213"
)
colnames(RESULTS) <- c("observed frequency", "95% quantile of simulated the frequency", "fail to reject Ho?")
RESULTS["pattern.11111",] <- c(OBSERV.vector.pattern.11111,quantile(PROB.vector.pattern.11111,0.95),TEST.PROB.pattern11111)
RESULTS["pattern.22222",] <- c(OBSERV.vector.pattern.22222,quantile(PROB.vector.pattern.22222,0.95),TEST.PROB.pattern22222)
RESULTS["pattern.33333",] <- c(OBSERV.vector.pattern.33333,quantile(PROB.vector.pattern.33333,0.95),TEST.PROB.pattern33333)
RESULTS["pattern.44444",] <- c(OBSERV.vector.pattern.44444,quantile(PROB.vector.pattern.44444,0.95),TEST.PROB.pattern44444)
RESULTS["pattern.23111",] <- c(OBSERV.vector.pattern.23111,quantile(PROB.vector.pattern.23111,0.95),TEST.PROB.pattern23111)
RESULTS["pattern.23113",] <- c(OBSERV.vector.pattern.23113,quantile(PROB.vector.pattern.23113,0.95),TEST.PROB.pattern23113)
RESULTS["pattern.23211",] <- c(OBSERV.vector.pattern.23211,quantile(PROB.vector.pattern.23211,0.95),TEST.PROB.pattern23211)
RESULTS["pattern.23213",] <- c(OBSERV.vector.pattern.23213,quantile(PROB.vector.pattern.23213,0.95),TEST.PROB.pattern23213)

# DO WE FAIL TO REJECT THE NULL HYPOTHESIS, Ho?
RESULTS


#################
##
#   PLOTS IN FIGURE 2:
#      - The empirical distribution of the probability of observing each of the tested choice patterns
#      - The dashed vertical line represents the 95th percentile of the simulated empirical distribution
#      - the dotted vertical line the observed relative frequency
##
#################

weight <- rep(1/n.iteraciones, length(PROB.vector.pattern.11111))
resultados_prob <- data.frame(PROB.vector.pattern.11111,PROB.vector.pattern.22222,PROB.vector.pattern.33333,PROB.vector.pattern.44444,PROB.vector.pattern.23111,PROB.vector.pattern.23113,PROB.vector.pattern.23211,PROB.vector.pattern.23213, weight)
breaks1 <- c(seq(from=0, to=0.246, by=0.0041))

# Pattern (1,1,1,1,1)
p1 <- ggplot(resultados_prob, aes(PROB.vector.pattern.11111, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 11111", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.11111, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.11111,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        )
p1

# Pattern (2,2,2,2,2)
p2 <- ggplot(resultados_prob, aes(PROB.vector.pattern.22222, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 22222", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.22222, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.22222,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p2

# Pattern (3,3,3,3,3)
p3 <- ggplot(resultados_prob, aes(PROB.vector.pattern.33333, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 33333", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.33333, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.33333,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p3

# Pattern (4,4,4,4,4)
p4 <- ggplot(resultados_prob, aes(PROB.vector.pattern.44444, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 44444", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.44444, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.44444,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p4

# Pattern (2,3,1,1,1,4)
p5 <- ggplot(resultados_prob, aes(PROB.vector.pattern.23111, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 23111", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.23111, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.23111,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p5

# Pattern (2,3,1,1,3)
p6 <- ggplot(resultados_prob, aes(PROB.vector.pattern.23113, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 23113", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.23113, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.23113,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p6

# Pattern (2,3,2,1,1)
p7 <- ggplot(resultados_prob, aes(PROB.vector.pattern.23211, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 23211", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.23211, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.23211,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p7

# Pattern (2,3,2,1,3)
p8 <- ggplot(resultados_prob, aes(PROB.vector.pattern.23213, weight = weight)) + geom_histogram(color="black", fill="white", breaks=breaks1)+  
  labs(x="Pattern 23213", y="Frequency") + xlim(c(0,0.25)) + ylim(c(0,1))+
  geom_vline(xintercept = OBSERV.vector.pattern.23213, linetype="dotted", color = "black", size=1.5)+
  geom_vline(xintercept = quantile(PROB.vector.pattern.23213,0.95), linetype="longdash",color = "black", size=1.5)+ 
  theme_classic()+
  theme(text = element_text(size=25),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  )
p8









