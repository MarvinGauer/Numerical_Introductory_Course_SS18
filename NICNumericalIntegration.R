# This is the R-Code with respect to the Topic "Numericla Integration"

# Install and load packages
libraries = c("ggplot2")
lapply(libraries, function(x) 
  if (!(x %in% installed.packages())) {
    install.packages(x)
  }
)
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Define Functions
MonteCarloIntegration = function(x = NULL, FUN = dnorm, n = 100000, graphic = TRUE){
  
  # x contains the x values and FUN is a function
  # n is an integer and represents the number of iterations
  
  # create DataFrame
  df = data.frame(1.5*x, y = FUN(1.5*x))
  
  # Generate Random Points
  mcp_x = runif(n, min = min(x), max = max(x))
  mcp_y = runif(n, min = min(FUN(x)), max = max(FUN(x)))
  
  # Calculate the area
  area = (max(x) - min(x))*(max(FUN(x)) - min(FUN(x)))
  
  # Approx. Area
  dfp = data.frame(mcp_x, mcp_y)
  dfp["Within"] = ifelse(dfp$mcp_y <= FUN(dfp$mcp_x), TRUE, FALSE)
  
  # Percentage of points <= function
  perc = length(dfp$mcp_x[dfp["Within"] == TRUE])/n
  
  if (graphic == TRUE){
    # Visualization
    p = ggplot(aes(x, y), data=df) +
      geom_line(size = 1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black", arrow = arrow(length = unit(0.25, "cm")))) +
      ylab("f(y)") + ylab("x") + 
      geom_point(aes(x=mcp_x*(2/3), y=mcp_y, colour = Within), data = dfp, size = 1, show.legend=F)
      
    print(p)
  }
  
  return(area * perc)
}

pol = function(x){
  y = exp(x) + 3
  return(y)
}

MonteCarloIntegration(seq(-4,4,0.01), FUN = pol, n = 100, graphic = T)


