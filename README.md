# levels: Sufficient Representations of Categorical Variables

This package implements the methods for providing sufficient representations of categorical variables mentioned in Johannemann et al. (2019).


To install this package in R, run the following commands:
```R
library(devtools)
install_github(gsbDBI/sufrep)
```

Example usage:
```
library(sufrep)

n = 5000
p = 20

X <- matrix(rnorm(n*p),n,p)
G <- sample(5,size=n,replace=TRUE)
df <- data.frame(cbind(X,G))

enc <- encoder(X=df,G="G",method = "one_hot")

train.df <- enc(df)

head(train.df)

```

#### References
Jonathan Johannemann, Vitor Hadad, Susan Athey, and Stefan Wager. Sufficient Representations of Categorical Variables. 2019.
