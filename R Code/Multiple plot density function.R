carrots <- data.frame(length = rnorm(100000, 6, 2))
cukes <- data.frame(length = rnorm(50000, 7, 2.5))
tomatos <- data.frame(length = rnorm(30000,3,8))
# Now, combine your two dataframes into one.  
# First make a new column in each that will be 
# a variable to identify where they came from later.
carrots$veg <- 'carrot'
cukes$veg <- 'cuke'
tomatos$veg <- 'tomato'

# and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes,tomatos)
ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)