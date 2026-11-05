###############################################################################  
# Our sample size formulas for comparing first-stage and second-stage options 
# rely on a two-sample z-test. We should expect generally higher values than 
# those from a two-sample t-test formula, such as the code below.  
#  
# Kelley Kidwell and team’s app 
# https://umich-biostat-smartsize.share.connect.posit.cloud/
# appear to use formulas based on a two-sample t-test. 
# Their results for comparing first-stage and second-stage options are 
# typically a bit higher than those from our app.  
###############################################################################


# Calculate sample size for a two-sample t-test
result <- power.t.test(
  n = NULL,            # Leave NULL to solve for sample size
  delta = 0.5,         # True difference in means (mu1 - mu2)
  sd = 1.0,            # Common standard deviation
  sig.level = 0.05,    # Significance level (Type I error rate, alpha)
  power = 0.80,        # Statistical power (1 - Type II error rate, beta)
  type = "two.sample", # Independent two-sample t-test
  alternative = "two.sided"
)

# Print the results
print(result)
# Access the required sample size per group
result$n 
