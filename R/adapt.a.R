#' Adaptive Alpha
#' @description Compute an alpha value adjusted for sample size. The adjusted value is based on
#' Perez and Pericchi's (2014) formula (equation 11, see below) using a reference sample, which can be
#' defined a priori or estimated using the sample size calculation from power.
#' 
#' \deqn{\frac{\alpha * \sqrt{n_0 \times (log(n_0) + \chi^{2}_{\alpha}(1))}}{\sqrt{n^* \times (log(n^*) + \chi^{2}_{\alpha}(1))}}}{\alpha * \sqrt(n0 times (log(n0) + \chi^2_\alpha(1))) / \sqrt(n* times (log(n*) + \chi^2_\alpha(1)))}
#' 
#' 
#' @param test Type of statistical test being used.
#' Can be any of the tests listed
#' 
#' @param ref.n \emph{n0} in the above equation.
#' Reference sample size.
#' If sample size was determined a priori, then the reference
#' number of participants can be set. This removes the calculation of sample
#' size based on power
#' 
#' @param n \emph{n*} in the above equation.
#' Number of participants in the experiment sample (or per group)
#' 
#' @param alpha \eqn{\alpha} in the above equation.
#' Alpha value to adjust.
#' Defaults to \code{.05}
#' 
#' @param power Power (\eqn{1 - \beta}) value.
#' Used to estimate the reference sample size (n0).
#' Defaults to \code{.80}
#' 
#' @param efxize Effect size to be used to estimate the reference sample size.
#' Effect sizes are based on Cohen (1992).
#' Numeric values can be used.
#' Defaults to \code{"medium"}
#' 
#' @param groups Number of groups (only for \code{test = "anova"})
#' 
#' @param df Number of degrees of freedom (only for \code{test = "chisq"})
#' 
#' @return A list containing the following objects:
#' 
#' \item{adapt.a}{The adapted alpha value}
#' 
#' \item{crit.value}{The critical value associated with the adapted alpha value}
#' 
#' \item{orig.a}{The original alpha value}
#' 
#' \item{ref.n}{The reference sample size based on alpha, power, effect size, and test}
#' 
#' \item{exp.n}{The sample size of the experimental sample}
#' 
#' \item{power}{The power used to determine the reference sample size}
#' 
#' \item{test}{The type of statistical test used}
#' 
#' @examples
#' #ANOVA
#' adapt.anova <- adapt.a(test = "anova", n = 200, alpha = .05, power = .80, groups = 3)
#' 
#' #Chi-square
#' adapt.chisq <- adapt.a(test = "chisq", n = 200, alpha = .05, power = .80, df = 3)
#' 
#' #Correlation
#' adapt.cor <- adapt.a(test = "cor", n = 200, alpha = .05, power = .80)
#' 
#' #One-sample t-test
#' adapt.one <- adapt.a(test = "one.sample", n = 200, alpha = .05, power = .80)
#' 
#' #Two-sample t-test
#' adapt.two <- adapt.a(test = "two.sample", n = 200, alpha = .05, power = .80)
#' 
#' #Paired sample t-test
#' adapt.paired <- adapt.a(test = "paired", n = 200, alpha = .05, power = .80, efxize = "medium")
#' 
#' @references
#' Cohen, J. (1992).
#' A power primer.
#' \emph{Psychological Bulletin}, \emph{112}, 155-159.
#' 
#' Perez, M. E., & Pericchi, L. R. (2014).
#' Changing statistical significance with the amount of information: The adaptive \emph{a} significance level.
#' \emph{Statistics & Probability Letters}, \emph{85}, 20-24.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats qf
#' 
#' @export
#Adaptive Alpha----
adapt.a <- function (test = c("anova","chisq","cor","one.sample","two.sample","paired"),
                     ref.n = NULL, n = NULL, alpha = .05, power = .80,
                     efxize = c("small","medium","large"), groups = NULL, df = NULL)
{
    if(missing(test))
    {stop("test must be selected")
    }else{test <- match.arg(test)}
    
    if(missing(efxize))
    {
        efxize <- "medium"
        message("No effect size selected. Medium effect size computed.")
    }else{efxize <- efxize}
    
    if(test=="anova")
    {
        if(is.null(groups))
        {stop("ANOVA is selected. Number of groups must be set")}
        
        if(efxize=="small")
        {efxize <- .10
        }else if(efxize=="medium")
        {efxize <- .25
        }else if(efxize=="large")
        {efxize <- .40}
        
        if(!is.numeric(efxize))
        {stop("Effect size must be numeric")}
        
        if(is.null(ref.n))
        {
            ref.n <- pwr::pwr.anova.test(f=efxize,power=power,sig.level=alpha,k=groups)$n
            message("ref.n is observations per group")
        }
        
        num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    }else if(test=="chisq")
    {
        if(is.null(df))
        {stop("Chi-square is selected. Degrees of freedom must be set")}
        
        if(efxize=="small")
        {efxize <- .10
        }else if(efxize=="medium")
        {efxize <- .30
        }else if(efxize=="large")
        {efxize <- .50}
        
        if(!is.numeric(efxize))
        {stop("Effect size must be numeric")}
        
        if(is.null(ref.n))
        {ref.n <- pwr::pwr.chisq.test(w=efxize,df=df,power=power,sig.level=alpha)$N}
        
        num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    }else if(test=="cor")
    {
        if(efxize=="small")
        {efxize <- .10
        }else if(efxize=="medium")
        {efxize <- .30
        }else if(efxize=="large")
        {efxize <- .50}
        
        if(!is.numeric(efxize))
        {stop("Effect size must be numeric")}
        
        if(is.null(ref.n))
        {ref.n <- pwr::pwr.r.test(r=efxize,power=power,sig.level=alpha)$n}
        
        num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    }else if(any(c("one.sample","two.sample","paired") %in% test))
    {
        if(efxize=="small")
        {efxize <- .20
        }else if(efxize=="medium")
        {efxize <- .50
        }else if(efxize=="large")
        {efxize <- .80}
        
        if(!is.numeric(efxize))
        {stop("Effect size must be numeric")}
        
        if(is.null(ref.n))
        {ref.n <- pwr::pwr.t.test(d=efxize,power=power,sig.level=alpha,type=test)$n}
        
        num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    }else{stop("test does not exist")}
 
    #denominator
    denom <- (sqrt(n*(log(n)+qchisq((1-alpha),1))))
    #adjusted alpha calculation
    adj.a <- alpha*num/denom
    
    #critical values
    if(test=="anova")
    {
        critical.f <- function (groups, n, a)
        {
            df1 <- groups - 1
            df2 <- n - groups
            cvf <- qf(a, df1, df2, lower.tail = FALSE)
            return(cvf)
        }
        
        cv <- critical.f(groups, n, adj.a)
    }else if(test=="chisq")
    {
        critical.chi <- function (df, a)
        {
            cvchi <- qchisq(a, df, lower.tail = FALSE)
            return(cvchi)
        }
        
        cv <- critical.chi(df, adj.a)
    }else if(test=="cor")
    {
        critical.r <- function (n, a)
            {
                df <- n - 2
                critical.t <- qt( a/2, df, lower.tail = FALSE )
                cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
                return(cvr)
            }
        
        cv <- critical.r(n, adj.a)
    }else if(any(c("one.sample","two.sample","paired") %in% test))
    {
        critical.t <- function (n, a)
            {
                df <- n - 2
                cvt <- qt( a/2, df, lower.tail = FALSE )
                return(cvt)
            }
        
        cv <- critical.t(n, adj.a)
    }
    
    #output
    output <- list()
    output$adapt.a <- adj.a
    output$crit.value <- cv
    output$orig.a <- alpha
    output$ref.n <- ref.n
    output$exp.n <- n
    if(test=="anova")
    {
        output$groups <- groups
        output$df <- c((groups - 1), (n - groups))
    }
    if(test=="chisq")
    {output$df <- df}
    output$power <- power
    output$efxize <- efxize
    output$test <- test
    
    return(output)
}
#----
