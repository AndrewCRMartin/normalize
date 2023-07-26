normalize
=========

(c) 2009 UCL, Andrew C.R. Martin
--------------------------------

`normalize` takes a set of datapoints, a target mean and a target
standard deviation and normalizes the datapoints to follow that
distribution.

Algorithm
---------

Given a set of datapoints *X*, a target mean, *mu*, and target
standard deviation, *s*:

```
foreach datapoint *Xi*
{
   Calculate absolute *|z| = |(X_i - u) / s|*
   Calculate probability *(p)* of a value *(x > |z|)*
   Generate a random number *0 <= r <= 1*
   if *(p >= r)* then place *Xi* in the output set
}
```


*p* is calculated as follows:

```
*D(x) = 0.5[1 + erf((x-u}/(s sqrt(2)))]*

 D(x) = 0.5[1 + erf((z)/(sqrt(2)))]
```

So:
```
p(|x-u| > |z|) = [1 + erf((-z)/(sqrt(2)))]
```

This is taken from
http://mathworld.wolfram.com/NormalDistribution.html and checked
http://www.stat.wvu.edu/SRS/Modules/Normal/males.html
