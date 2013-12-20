# nominal associations
# (C) Antti Arppe 2007-2011
# E-mail: antti.arppe@helsinki.fi

# This is an R function that calculates a range of measures of
# association for nominal variables.  In addition to eachstatistic
# value, the variance, ASE and significance (P-level) is computed The
# formulas are based on the following sources:
 
# Cohen, Jacob. 1988. Statistical power analysis for the behavioral
# sciences, (2nd edition). Lawrence Erlbaum Associates, Hillsdale, New
# Jersey, United States.

# Garson, G. David. 2007. Statnotes: Topics in Multivariate
# Analysis. URL:
# http://www2.chass.ncsu.edu/garson/pa765/statnote.htm. Visited Spring
# 2006 -- Summer 2007.

# Goodman, Leo A. and William H. Kruskal. 1954. Measures of Association
# for Cross-Classifications. Journal of the American Statistical
# Association, Vol. 49, No. 268 (December 1954), pp. 732-764.

# Liebetrau, Albert M. 1983. Measures of Association. Sage University
# Paper series on Quantitative Applications in the Social Sciences,
# 07-032. Sage Publications, Beverly Hills and London, United
# States/England.

# Margolin, Barry H. and Richard J. Light. 1974. An Analysis of Variance
# for Categorical Data II: Small Sample Comparisons with Chi Square and
# Other Competitors. Journal of the American Statistical Association,
# Vol. 69, No. 347 (September 1974), pp. 755-764.

# Reynolds, H. T. 1977. Analysis of Nominal Data. Sage University Paper
# series on Quantitative Applications in the Social Sciences, 08-007,
# Sage Publications, Beverly Hills/London, California/UK.

# SAS Institute. 2007. Measures of Association
# http://support.sas.com/onlinedoc/913/getDoc/en/statug.hlp/freq_sect20.htm
# Visited January 2007.

# Theil, Henri. 1970. On the Estimation of Relationships Involving
# Qualitative Variables.  The American Journal of Sociology, Vol. 76,
# No. 1 (July 1970), pp. 103-154.

# Furthermore, I have had access to a similar function script
# <measures.R> by Marc Schwartz, from whom I have also received valuable
# assistance in finding sources for the computation of the variances,
# and thus the other statistics based on them.

# N.B. One should use the values for the significance of the
# Goodman-Kruskal lambda and Theil's UC with reservation, as these
# have been modeled to mimic the the behavior of the same statistics
# in SPSS.

associations <- function(ctable, alpha=0.05, p.zero.correction=1/sum(ctable)^2)
{ N <- sum(ctable);
  ctable.rows <- nrow(ctable);
  ctable.cols <- ncol(ctable);
  max.col <- sum.col <- L.col <- matrix(,ctable.cols);
  max.row <- sum.row <- L.row <- matrix(,ctable.rows);
  for(i in 1:ctable.cols)
     { sum.col[i] <- sum(ctable[,i]); max.col[i] <- max(ctable[,i]); }
  for(i in 1:ctable.rows)
     { sum.row[i] <- sum(ctable[i,]); max.row[i] <- max(ctable[i,]); }
  max.row.margin <- max(apply(ctable,1,sum));
  max.col.margin <- max(apply(ctable,2,sum));

# slightly nudge zero values so that their logarithm can be calculated (cf. Theil 1970: x->0 => xlogx->0)
  for(i in 1:ctable.rows)
     for(j in 1:ctable.cols)
        if(ctable[i,j]==0)
          ctable[i,j]=p.zero.correction;

# Pearson chi-squared (X2) test of independence (homogeneity)
  chisq.results <- suppressWarnings(chisq.test(ctable));
  alpha.X2 <- chisq.results$p.value;

# Log-likelihood chi-squared (G2) test of independence (homogeneity)
  likelihood.ratio <- 2*sum(chisq.results$observed*log(chisq.results$observed/chisq.results$expected))
  alpha.G2 <- pchisq(likelihood.ratio,df=(ctable.rows-1)*(ctable.cols-1),lower.tail=FALSE);

# Cramér's V
  cramers.v <- sqrt(chisq.results$statistic/(N*(min(ctable.rows,ctable.cols)-1)));
  names(cramers.v) <- NULL;

# Goodman-Kruskal lambda (1954)
  lambda.RC <- (sum(max.col)-max.row.margin)/(N-max.row.margin);
  lambda.CR <- (sum(max.row)-max.col.margin)/(N-max.col.margin);
  L.col.max <- min(which(sum.col==max.col.margin));
  for(i in 1:ctable.rows)
     { if(length(which(ctable[i,intersect(which(ctable[i,]==max.col.margin),which(ctable[i,]==max.row.margin))]==N))>0)
         L.col[i] <- min(which(ctable[i,intersect(which(ctable[i,]==max.col.margin),which(ctable[i,]==max.row.margin))]==N))
       else
         if(ctable[i,L.col.max]==max.col.margin)
           L.col[i]=L.col.max
         else
           L.col[i]=min(which(ctable[i,]==max.row[i]));
     }
  var.lambda.CR <- (N-sum(max.row))*(sum(max.row)+max.col.margin-2*(sum(max.row[which(L.col==L.col.max)])))/(N-max.col.margin)^3;
  ASE.lambda.CR <- sqrt(var.lambda.CR);
  L.row.max <- min(which(sum.row==max.row.margin));
  for(i in 1:ctable.cols)
     { if(length(which(ctable[intersect(which(ctable[,i]==max.row.margin),which(ctable[,i]==max.col.margin)),i]==N))>0)
         L.row[i] <- min(which(ctable[i,intersect(which(ctable[i,]==max.col.margin),which(ctable[i,]==max.row.margin))]==N))
       else
         if(ctable[L.row.max,i]==max.row.margin)
           L.row[i]=L.row.max
         else
           L.row[i]=min(which(ctable[,i]==max.col[i]));
     }
  var.lambda.RC <- (N-sum(max.col))*(sum(max.col)+max.row.margin-2*(sum(max.col[which(L.row==L.row.max)])))/(N-max.row.margin)^3;
  ASE.lambda.RC <- sqrt(var.lambda.RC);
  z.lambda.RC <- lambda.RC/ASE.lambda.RC;
  z.lambda.CR <- lambda.CR/ASE.lambda.CR;
  p.lambda.RC <- (pnorm(z.lambda.RC,lower.tail=FALSE))*2;
  p.lambda.CR <- (pnorm(z.lambda.CR,lower.tail=FALSE))*2;
#  p.lambda.RC <- 1-pnorm(z.lambda.RC);
#  p.lambda.CR <- 1-pnorm(z.lambda.CR);

# Goodman-Kruskal tau (Liebetrau 1983)

# Tau Column|Row
  n.err.unconditional <- N^2;
  for(i in 1:ctable.rows)
     n.err.unconditional <- n.err.unconditional-N*sum(ctable[i,]^2/sum.row[i]);
  n.err.conditional <- N^2-sum(sum.col^2);
  tau.CR <- 1-(n.err.unconditional/n.err.conditional);
  v <- n.err.unconditional/(N^2);
  d <- n.err.conditional/(N^2);
  f <- d*(v+1)-2*v;
  var.tau.CR <- 0;
  for(i in 1:ctable.rows)
     for(j in 1:ctable.cols)
        var.tau.CR <- var.tau.CR + ctable[i,j]*(-2*v*(sum.col[j]/N)+d*((2*ctable[i,j]/sum.row[i])-sum((ctable[i,]/sum.row[i])^2))-f)^2/(N^2*d^4);
  ASE.tau.CR <- sqrt(var.tau.CR);
  z.tau.CR <- tau.CR/ASE.tau.CR;
#  p.tau.CR <- 1-pnorm(z.tau.CR);
  U.tau.CR <- (N-1)*(ctable.cols-1)*tau.CR; # Chi-squared approximation for H0 according to Margolin & Light (1974), see also Liebetrau (1983)
  p.tau.CR <- pchisq(U.tau.CR,df=(ctable.rows-1)*(ctable.cols-1),lower.tail=FALSE);

# Tau Row|Column
  n.err.unconditional <- N^2;
  for(j in 1:ctable.cols)
     n.err.unconditional <- n.err.unconditional-N*sum(ctable[,j]^2/sum.col[j]);
  n.err.conditional <- N^2-sum(sum.row^2);
  tau.RC <- 1-(n.err.unconditional/n.err.conditional);
  v <- n.err.unconditional/(N^2);
  d <- n.err.conditional/(N^2);
  f <- d*(v+1)-2*v;
  var.tau.RC <- 0;
  for(i in 1:ctable.rows)
     for(j in 1:ctable.cols)
        var.tau.RC <- var.tau.RC + ctable[i,j]*(-2*v*(sum.row[i]/N)+d*((2*ctable[i,j]/sum.col[j])-sum((ctable[,j]/sum.col[j])^2))-f)^2/(N^2*d^4);
  ASE.tau.RC <- sqrt(var.tau.RC);
  z.tau.RC <- tau.CR/ASE.tau.RC;
#  p.tau.RC <- 1-pnorm(z.tau.RC);
  U.tau.RC <- (N-1)*(ctable.rows-1)*tau.RC; # Chi-squared approximation for H0 according to Margolin & Light 1974, see also Liebetrau 1983
  p.tau.RC <- pchisq(U.tau.RC,df=(ctable.rows-1)*(ctable.cols-1),lower.tail=FALSE);

# Theil's UC (1970)
  hx= -sum((apply(ctable,1,sum)*log(apply(ctable,1,sum)/N))/N);
  hy= -sum((apply(ctable,2,sum)*log(apply(ctable,2,sum)/N))/N);
  hxy= -sum(apply(ctable,c(1,2),sum)*log(apply(ctable,c(1,2),sum)/N)/N);
  uc.RC <- (hx+hy-hxy)/hx;
  uc.CR <- (hx+hy-hxy)/hy;
  uc.sym <- 2*(hx+hy-hxy)/(hx+hy);
  var.uc.RC <- var.uc.CR <- 0;
  for(i in 1:ctable.rows)
     for(j in 1:ctable.cols)
        { var.uc.RC <- var.uc.RC + ctable[i,j]*(hx*log(ctable[i,j]/sum.col[j])+((hy-hxy)*log(sum.row[i]/N)))^2/(N^2*hx^4);
          var.uc.CR <- var.uc.CR + ctable[i,j]*(hy*log(ctable[i,j]/sum.row[i])+((hx-hxy)*log(sum.col[j]/N)))^2/(N^2*hy^4);
        }
  ASE.uc.RC <- sqrt(var.uc.RC);
  ASE.uc.CR <- sqrt(var.uc.CR);
  z.uc.RC <- uc.RC/ASE.uc.RC;
  z.uc.CR <- uc.CR/ASE.uc.CR;
  p.uc.RC <- pnorm(z.uc.RC,lower.tail=FALSE)*2; # N.B. SPSS appears to use the log-likelihood chi-squared test (G2) to assess the significance of UC.
  p.uc.CR <- pnorm(z.uc.CR,lower.tail=FALSE)*2; # Thus, the significance value is the same in both directions.
#  p.uc.RC <- 1-pnorm(z.uc.RC);
#  p.uc.CR <- 1-pnorm(z.uc.CR);

# Cohen's Effect Size (1988)
  e0 <- matrix(,ctable.rows,ctable.cols)
  for(i in 1:ctable.rows)
     for(j in 1:ctable.cols)
        e0[i,j] <- sum.row[i]*sum.col[j]/N;
  p0 <- e0/N;
  p1 <- ctable/N;
  effect.size <- sqrt(sum(((p1-p0)^2)/p0));
  noncentrality <- N*(effect.size^2);
  d.f=(ctable.rows-1)*(ctable.cols-1);
  beta <- pchisq(qchisq(alpha,df=d.f,lower.tail=FALSE),df=d.f,ncp=noncentrality);
  power <- 1-beta;

  results <- list(alpha.X2 = alpha.X2, alpha.G2 = alpha.G2, beta = beta, power = power, effect.size = effect.size, likelihood.ratio = likelihood.ratio, cramers.v = cramers.v, lambda.RC = lambda.RC, lambda.CR = lambda.CR, tau.RC = tau.RC, tau.CR = tau.CR, uc.RC = uc.RC, uc.CR = uc.CR, uc.sym = uc.sym, p.lambda.RC = p.lambda.RC, p.lambda.CR = p.lambda.CR, p.tau.RC = p.tau.RC, p.tau.CR = p.tau.CR, p.uc.RC = p.uc.RC, p.uc.CR = p.uc.CR, var.lambda.RC = var.lambda.RC, var.lambda.CR = var.lambda.CR, var.tau.RC = var.tau.RC, var.tau.CR = var.tau.CR, var.uc.RC = var.uc.RC, var.uc.CR = var.uc.CR, ASE.lambda.RC = ASE.lambda.RC, ASE.lambda.CR = ASE.lambda.CR, ASE.tau.RC = ASE.tau.RC, ASE.tau.CR = ASE.tau.CR, ASE.uc.RC = ASE.uc.RC, ASE.uc.CR = ASE.uc.CR, noncentrality = noncentrality);

  return(results);
}
