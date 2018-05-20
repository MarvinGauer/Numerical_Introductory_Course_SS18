# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This is the R-Code with respect to the Topic "Numericla Integration"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clean the vurrent workspace
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

rm(list = ls(all = TRUE))
graphics.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Install and load packages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

libraries = c("ggplot2")
lapply(libraries, function(x) 
  if (!(x %in% installed.packages())) {
    install.packages(x)
  }
)
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Input Values
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

u    = 4     # Upper Boundary
l    = -4    # Lower Boundary

# Function to numerically integrate
pol  = function(x){
  y  = exp(x) * x^2 + 3*x + 4
  return(y)
}

n    = 100 # Max Number of Interations

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

MonteCarloIntegration = function(x = NULL, FUN = dnorm, n = 100000, graphic = TRUE){
  
  # x contains the x values and FUN is a function
  # n is an integer and represents the number of iterations
  
  # create DataFrame
  df = data.frame(1.5*x, y = FUN(1.5*x))
  
  # Generate Random Points
  mcp_x = runif(n, min = min(x), max = max(x))
  mcp_y = runif(n, min = 0, max = max(FUN(x)))
  
  # Calculate the area
  area = (max(x) - min(x))*(max(FUN(x)))
  
  # Approx. Area
  dfp = data.frame(mcp_x, mcp_y)
  dfp["Within"] = ifelse(dfp$mcp_y <= FUN(dfp$mcp_x), TRUE, FALSE)
  
  # Percentage of points <= function
  perc = length(dfp$mcp_x[dfp["Within"] == TRUE])/n
  
  if (graphic == TRUE){
    # Visualization of the graph
    p = ggplot(aes(df[,1], df[,2]), data=df) +
      geom_line(size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      geom_point(aes(x=mcp_x, y=mcp_y, colour = Within), data = dfp, size = 1, show.legend=F) + 
      ylab("f(y)") + xlab("x")
    
    print(p)
  }
  
  sol = list(as.integer(n),area * perc)
  names(sol) = c("Iterations","Area")
  
  return(sol)
}

MonteCarloIteration = function(x = NULL, FUN = dnorm, n = 100000, graphic = TRUE){
  
  # x contains the x values and FUN is a function
  # n is an integer and represents the number of iterations (min is 10)
  
  # Convergence Data.Frame
  dfg = data.frame("Iteration"     = as.integer(), 
                   "Approx. Value" = as.numeric(), 
                   "Real Value"    = as.numeric(),
                   "Difference"    = as.numeric())
  
  for(i in seq(10,n,10)){
    dfg[nrow(dfg) + 1,] = c(i,
                            MonteCarloIntegration(x = x, FUN = FUN, n = i, graphic = F)$Area,
                            integrate(pol,min(x),max(x))$value,
                            integrate(pol,min(x),max(x))$value - MonteCarloIntegration(seq(l,u,0.01), FUN = pol, n = i, graphic = F)$Area)
  }
  
  if (graphic == TRUE){
    
    # Visualization of the Function
    MonteCarloIntegration(x = x, FUN = pol, n = n, graphic = T)
    
    # Visualization of the convergence
    g = ggplot(aes(x = dfg[,1], y = dfg[,2]), data=dfg) +
      geom_line() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      geom_line(aes(x = dfg[,1], y = dfg[,3], colour = 'red'), data = dfg, show.legend=F) +
      xlab("Iterations") + ylab("Area")
    
    print(g)
  }
  
  return(dfg)
}


MidpointIntegration = function(l = NULL, u = NULL, n = 10, FUN = dnorm, graphic = T){
  
  # FUN is the function of interest
  
  x_mid = head(filter(seq(l,u,(u-l)/(n-1)),c(0.5,0.5)),-1)
  y_mid = FUN(x_mid)
  
  # create DataFrame
  df   = data.frame(seq(l,u,(u-l)/(n-1)), y = FUN(seq(l,u,(u-l)/(n-1))))
  rect = data.frame(xl = seq(l,u,(u-l)/(n-1))[1:length(seq(l,u,(u-l)/(n-1)))-1], xr = seq(l,u,(u-l)/(n-1))[2:length(seq(l,u,(u-l)/(n-1)))], yu = rep(0,length(seq(l,u,(u-l)/(n-1)))-1), yo = y_mid )
  if (graphic == TRUE){
    # Visualization of the approximation
    g = ggplot() +
      geom_line(aes(x = df[,1], y = df[,2]), data=df) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      geom_rect(data=rect, aes(xmin=rect[,1], 
                               xmax=rect[,2], 
                               ymin=rect[,3], 
                               ymax=rect[,4]),
                fill = 'red', alpha = 0.2, col = 'blue') +
      xlab("x") + ylab("f(x)")
    
    print(g)
    
  }
  
  sol   = sum(y_mid %*% diff(seq(l,u,(u-l)/(n-1))))
  
  return(sol)
}

MidpointIteration = function(l = NULL, u = NULL, n = 10, FUN = dnorm, graphic = T){
  
  # FUN is the function of interest
  
  steps = length(head(filter(seq(l,u,(u-l)/(n-1)),c(0.5,0.5)),-1))
  
  # Convergence Data.Frame
  dfg = data.frame("Iteration"     = as.integer(), 
                   "Approx. Value" = as.numeric(), 
                   "Real Value"    = as.numeric(),
                   "Difference"    = as.numeric())
  
  for(i in seq(2,steps,1)){
    dfg[nrow(dfg) + 1,] = c(i,
                            MidpointIntegration(l = l, u = u, n = i, FUN = FUN, graphic = F),
                            integrate(pol,l,u)$value,
                            integrate(pol,l,u)$value - MidpointIntegration(l = l, u = u, n = n, FUN = FUN, graphic = F))
  }
  if (graphic == TRUE){
    # Visualization of the Function
    MidpointIntegration(l = l, u = u, n = n, FUN = FUN, graphic = T)
    
    # Visualization of the convergence
    g = ggplot(aes(x = dfg[,1], y = dfg[,2]), data=dfg) +
      geom_line() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      geom_line(aes(x = dfg[,1], y = dfg[,3], colour = 'red'), data = dfg, show.legend=F) +
      xlab("Iterations") + ylab("Area")
    
    print(g)
  }
  
  return(dfg)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculations & Results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# Midpoint Approach
MidpointIntegration(l = l, u = u, n = n, FUN = pol, graphic = F)
MidpointIteration(l = l, u = u, n = n, FUN = pol, graphic = T)

# MonteCarlo
MonteCarloIntegration(seq(l,u,0.01), FUN = pol, n = n*100, graphic = F)$Area
MonteCarloIteration(seq(l,u,0.01), FUN = pol, n = n*100, graphic = T)



