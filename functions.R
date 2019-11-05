library(ggplot2)

#
# t-test
#
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) {
    # welch-satterthwaite equation
    se = sqrt( (s1^2/n1) + (s2^2/n2) )
    df = ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else {
    # pooled standard deviation, scaled by the sample sizes
    se = sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df = n1+n2-2
  }      
  
  t = (m1-m2-m0)/se 
  dat = data.frame("Difference of means"=m1-m2, "Std Error"=se, "t"=t, "p_value"=2*pt(-abs(t),df))
  return(dat) 
}

myTheme <- list(theme_minimal(base_size = 15), 
                theme(panel.background = element_rect(color = "grey90"), 
                      panel.grid.minor = element_line(colour = "white"), 
                      panel.grid.major = element_line(colour = "grey95")),
                geom_vline(xintercept = 0, color = "grey50") , 
                geom_hline(yintercept = 0, color = "grey50"),
                theme(axis.title.y = element_text(size = rel(0.8), angle = 90)),
                theme(axis.title.x = element_text(size = rel(0.8), angle = 00)),
                theme(legend.title = element_text(size = rel(0.8))),
                theme(legend.text = element_text(size = rel(0.8))))



log.na = function(x) log10(ifelse(x>0, x, NA))
log.zero = function(x) log10(ifelse(x>0, x, 1))
glog2 = function(x) ((asinh(x)-log(2))/log(2))


strset = function(x, f=T) {
  x = sort(unique(na.omit(x[f])))
  if(length(x)==0) return(NA_character_)
  paste0(x, collapse="|")
}

chunk = function(x,n) {
  split(x, rep(1:ceiling(length(x)/n), each=n)[1:length(x)])
}