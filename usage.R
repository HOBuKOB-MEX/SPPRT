source("SPPRT.R")


# Cost function. Returns m by default
cost <- function(m){
  1 + 0.01 * m
}

# Action set for each step of evaluation
step_fn<-function(step){
  print(paste("step",step))
}

# counts groups
count_groups_fn<-function(m)
  1

# counts number of observations
count_observations_fn <- function(m)
  m


# input parameters can be set here

th0 = 0.52  # theta0
th1 = 0.48  # theta1
gamma = 0.5
gsizes = seq(10, 600, by=10)  # eligible group sizes
H = 15  # maximum number of steps
h = 0.1  # grid size
l0 = 44  # lambda0
l1 = 44  # lambda1


print("Designing a test")

# designing a test with given input parameters
test = design_test(l0, l1, th0, th1, H, gsizes, gamma, h, step_fn=step_fn)


print("Calculating error probabilities")

# type I error calculation
print(paste("alpha=",1 - operating_characteristic(test, th0, step_fn)))

# type II error calculation
print(paste("beta=", operating_characteristic(test, th1, step_fn)))


print("Calculating average sampling cost (ASC)")
print(paste("ASC under H0 =", sampling_report(test, th0, step_fn=step_fn)))
print(paste("ASC under H1 =", sampling_report(test, th1)))
print(paste("ASC under theta=0.5 =", sampling_report(test, 0.5)))


print("Calculating average number of groups (ANG)")
print(paste("ANG under H0 =", sampling_report(test, th0, step_fn=step_fn, accounting_fn=count_groups_fn)))
print(paste("ANG under H1 =", sampling_report(test, th1, accounting_fn=count_groups_fn)))


print("Calculating average number of observations (ANO)")
print(paste("ANO under H0 =", sampling_report(test, th0, accounting_fn=count_observations_fn)))
print(paste("ANO under H1 =", sampling_report(test, th1, step_fn=step_fn, accounting_fn=count_observations_fn)))
