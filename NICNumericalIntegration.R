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
  y  = x^2 + 3*x + 4
  return(y)
}

n    = 10 # Max Number of Iterations
m    = 10 # Number of Bins

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Crude_MonteCarloIntegration = function(l = NULL, u = NULL, FUN = dnorm, n = 100, m = 10, graphic = TRUE){
  
  # l,u see other functions
  # n is an integer and represents the number of iterations per bin
  # and m is the number of bins
  
  x = seq(l,u,0.01)
  
  # create DataFrame
  df = data.frame(1.5*x, y = FUN(1.5*x))
  
  # create equidistant breaks
  l      = max(x)-min(x)
  step   = l/m
  breaks = seq(min(x),max(x),step)
  
  # Generate Random Points
  mcp_x = runif(n*m, min = min(x), max = max(x))
  
  # Calc corresponding y's
  mcp_y = FUN(mcp_x)
  
  # Calc y's mean within bins
  dfp   = data.frame(mcp_x, mcp_y, "Mean" = vector(length = length(mcp_y)))
  means = vector(length = length(breaks)-1)
  
  for(i in 1:length(breaks)-1){
    means[i] = mean(dfp$mcp_y[dfp$mcp_x >= breaks[i] & dfp$mcp_x <= breaks[i+1]])
    dfp$Mean[dfp$mcp_x >= breaks[i] & dfp$mcp_x <= breaks[i+1]] = means[i]
  }
  
  if (graphic == TRUE){
    
    rect = data.frame(xr = breaks[2:length(breaks)],
                      xl = breaks[1:length(breaks)-1],
                      yu = rep(0,length(breaks)-1),
                      yo = means)
    
    # Visualization of the graph
    p = ggplot(aes(df[,1], df[,2]), data=df) +
      geom_line(size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      ylab("f(y)") + xlab("x") +
      geom_point(data = dfp, aes(x = mcp_x, y = mcp_y, col = 'blue'),size = 3,inherit.aes = F, shape = 20, show.legend = F) + 
      geom_point(data = dfp, aes(x = mcp_x, y = rep(0,length(mcp_x)), col = 'green'),size = 3,inherit.aes = F, shape = 20, show.legend = F) + 
      geom_rect(data=rect, inherit.aes = F,aes(xmin=rect$xl, 
                               xmax=rect$xr, 
                               ymin=rect$yu, 
                               ymax=rect$yo),
                fill = 'red', alpha = 0.2, col = 'red')
    
    print(p)
  }
  
  sol = rep(step,length(means)) %*% means
  
  return(sol)
}

Crude_MonteCarloIteration = function(l = NULL, u = NULL, FUN = dnorm, n = 100, m = 10, graphic = TRUE){
  
  # l,u see other functions
  # n is an integer and represents the number of iterations (min is 10) per bin
  # m is the number of bins and is fixed
  
  x = seq(l,u,0.01)
  
  # Convergence Data.Frame
  dfg = data.frame("Iteration"     = as.integer(), 
                   "Approx. Value" = as.numeric(), 
                   "Real Value"    = as.numeric(),
                   "Difference"    = as.numeric())
  
  for(i in seq(10,n,10)){
    dfg[nrow(dfg) + 1,] = c(i,
                            Crude_MonteCarloIntegration(l = l, u = u, FUN = FUN, n = i, m = m, graphic = F),
                            integrate(pol,min(x),max(x))$value,
                            integrate(pol,min(x),max(x))$value - Crude_MonteCarloIntegration(l = l, u = u, FUN = FUN, n = i, m = m, graphic = F))
  }
  
  if (graphic == TRUE){
    
    # Visualization of the Function
    Crude_MonteCarloIntegration(l = l, u = u, FUN = FUN, n = n, m = m, graphic = T)
    
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

Hit_Miss_MonteCarloIntegration = function(l = NULL, u = NULL, FUN = dnorm, n = 100000, graphic = TRUE){
  
  # l,u see other functions
  # n is an integer and represents the number of iterations
  
  x = seq(l,u,0.01)
  
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

Hit_Miss_MonteCarloIteration = function(l = NULL, u = NULL, FUN = dnorm, n = 100000, graphic = TRUE){
  
  # l,u see other functions
  # x contains the x values and FUN is a function
  # n is an integer and represents the number of iterations (min is 10)
  
  x = seq(l,u,0.01)
  
  # Convergence Data.Frame
  dfg = data.frame("Iteration"     = as.integer(), 
                   "Approx. Value" = as.numeric(), 
                   "Real Value"    = as.numeric(),
                   "Difference"    = as.numeric())
  
  for(i in seq(10,n,10)){
    dfg[nrow(dfg) + 1,] = c(i,
                            Hit_Miss_MonteCarloIntegration(l = l, u = u, FUN = FUN, n = i, graphic = F)$Area,
                            integrate(pol,min(x),max(x))$value,
                            integrate(pol,min(x),max(x))$value - Hit_Miss_MonteCarloIntegration(l = l, u = u, FUN = pol, n = i, graphic = F)$Area)
  }
  
  if (graphic == TRUE){
    
    # Visualization of the Function
    Hit_Miss_MonteCarloIntegration(l = l, u = u, FUN = pol, n = n, graphic = T)
    
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

SimpsonIntegration = function(l = NULL, u = NULL, n = 10, FUN = dnorm, graphic = T){
  
  # FUN is the function of interest
  
  x_mid = head(filter(seq(l,u,(u-l)/(n-1)),c(0.5,0.5)),-1)
  y_mid = FUN(x_mid)
  
  # create DataFrame
  df   = data.frame("interval" = 1:length(x_mid), "x_i" = seq(l,u,(u-l)/(n-1))[1:(n-1)],
                    "x_mid" = x_mid, "x_i+1" = seq(l,u,(u-l)/(n-1))[2:(n)],
                    "y_i" = FUN(seq(l,u,(u-l)/(n-1))[1:(n-1)]), "y_mid" = FUN(x_mid),
                    "y_i+1" = FUN(seq(l,u,(u-l)/(n-1))[2:(n)]))
  
  # Function to fit: y = a*x^2 + b*x + c 
  
  para  = list()
  fct   = list()
  fct_y = list()
  fct_x = list()
  area  = list()
  
  for(i in 1:length(x_mid)){
    A = matrix(c(df$x_i[i]^2, df$x_i[i], 1,
                 df$x_mid[i]^2, df$x_mid[i], 1, 
                 df$x_i.1[i]^2, df$x_i.1[i],1), ncol = 3, byrow = T)
    y = c(df$y_i[i],df$y_mid[i],df$y_i.1[i])
    
    z = solve(A) %*% y
    
    para[[i]] = z
    fct[[i]]  = function(x, c = z){
      v = c[1] * x^2 + c[2] * x + c[3]
      return(v)
    }
    fct_y[[i]] = sapply(seq(df$x_i[i],df$x_i.1[i],0.01), fct[[i]], c = z)
    fct_x[[i]] = seq(df$x_i[i],df$x_i.1[i],0.01)
    
    #fct_y[[i]] = sapply(c(df$x_i[i],df$x_mid[i], df$x_i.1[i]), fct[[i]], c = z)
    #fct_x[[i]] = c(df$x_i[i],df$x_mid[i], df$x_i.1[i])
    area[i]   = integrate(fct[[i]],lower = df$x_i[i],upper = df$x_i.1[i])$value
  }
  
  df1    = data.frame(seq(l,u,0.01), y = FUN(seq(l,u,0.01)))
  df_fct = data.frame(matrix(c(unlist(fct_x),unlist(fct_y)), ncol = 2, byrow=F))
  
  
  if (graphic == TRUE){
    # Visualization of the approximation
    g = ggplot() +
      geom_line(aes(x = df1[,1], y = df1[,2]), data=df1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      xlab("x") + ylab("f(x)") + 
      geom_line(aes(x = df_fct[,1], y = df_fct[,2], col = 'red'), data=df_fct, show.legend = F)
      df2 = data.frame(x = seq(l,u,(u-l)/(n-1))[1:n])
      for(t in 1:length(df2$x)){
        g = g + geom_vline(xintercept = df2$x[t], linetype='dashed', alpha = 0.4)
      }
    print(g)
    
  }
  
  sol = sum(unlist(area))
  return(sol)
}
#AdaptiveIntegration

IntervalShifter = function(FUN = NULL, b = as.vector(length = 2)){
  
  if(b[1] == -Inf & b[2] == Inf){
    y = function(t){
      z = FUN(t/(1-t^2)) * (1+t^2)/((1-t^2)^2)
      return(z)
    }
  } else if(as.numeric(b[1]) == -Inf & b[2] != Inf){
    y = function(t){
      z = FUN(as.numeric(b[1]) - ((1-t)/(t))) * (1)/((t^2))
      return(z)
    }
  } else if(b[1] != -Inf & b[2] == Inf){
    y = function(t){
      z = FUN(as.numeric(b[1]) + (t/(1-t))) * (1)/((1-t)^2)
      return(z)
    }
  }
  
  return(y)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculations & Results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Midpoint Approach
MidpointIntegration(l = l, u = u, n = n/5, FUN = pol, graphic = F)
MidpointIteration(l = l, u = u, n = n/5, FUN = pol, graphic = T)

# MonteCarlo
Hit_Miss_MonteCarloIntegration(l,u, FUN = pol, n = n*100, graphic = F)$Area
Hit_Miss_MonteCarloIteration(l,u, FUN = pol, n = n*100, graphic = T)

Crude_MonteCarloIntegration(l,u,pol,n = n*10, m = 15)
Crude_MonteCarloIteration(l,u,pol,n = n*10,m = 15,T) # m is fixed here

Crude_MonteCarloIntegration(l,u,pol,n = n*10, m = 30)
Crude_MonteCarloIteration(l,u,pol,n = n*10,m = 30,T) # m is fixed here

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Application
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

x      = seq(-3,3,1)
bins   = 15*(1:7)
cutoff = -6
y      = vector(length=length(x)) 

for(i in 1:length(x)){
  y[i] = MidpointIntegration(l = cutoff, u = x[i], n = bins[i], FUN = dnorm, graphic = T)
}

fct    = lm(y ~ poly(x, 3, raw=TRUE))
x_new  = runif(6,-4,4)
d      = data.frame(x = x_new, r_F = pnorm(x_new), e_F = predict(fct,data.frame(x =x_new)), Diff = pnorm(x_new)-predict(fct,data.frame(x =x_new)))
d_plot = data.frame(x = seq(-3,3,0.01), r_F = pnorm(seq(-3,3,0.01)), e_F = predict(fct,data.frame(x =seq(-3,3,0.01))), Diff = pnorm(seq(-3,3,0.01))-predict(fct,data.frame(x =seq(-3,3,0.01))))

g = ggplot() +
  geom_line(aes(x = x, y = r_F, col = 'black'), data=d_plot, show.legend = F) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
  geom_line(aes(x = x, y = e_F, colour = 'red'), data = d_plot, show.legend = F) +
  xlab("x") + ylab("F(x)") +
  geom_point(aes(x = x, y = y, colour = 'black'), data = data.frame(x = x, y = y), shape = 8,show.legend = F)


print(g)