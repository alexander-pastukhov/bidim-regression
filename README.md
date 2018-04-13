# BiDimRegression

Package to calculate the bidimensional regression between two 2D configurations following the approach by Tobler (1965).

The package is described in detail in: Carbon, C. C. (2013). BiDimRegression: Bidimensional Regression Modeling Using R. Journal of Statistical Software,\ Code Snippets, 52(1), 1-11, doi: [10.18637/jss.v052.c01](http://dx.doi.org/10.18637/jss.v052.c01).

## Installation from Github

```
library("devtools"); install_github("alexander-pastukhov/bidim-regression",dependencies=TRUE)
```

## What's new in version 2
Added `lm2` S3-class to implement an interface consistent with that of R in general and of `lm` in particular.

## Examples
* Using legacy `BiDimRegression` function:
```
print(BiDimRegression(NakayaData))
```

* Using `lm2` S3-class:
```
lm2euc <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'euclidean')
lm2aff <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'affine')
lm2prj <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'projective')
anova(lm2euc, lm2aff, lm2prj
predict(lm2euc)
summary(lm2euc)
```

## Thanks

The author (CCC) is grateful to Waldo R. Tobler, now Professor Emeritus at the Department of Geography, University of California, Santa Barbara, for providing his original publications and his helpful correspondence and to Dan Montello for calling the author's attention to Tobler's work many years ago. I would like to thank Alinda Friedman, Gregory Francis, Jan de Leeuw, Achim Zeileis, two anonymous reviewers and Arne Terkowski for valuable comments on an earlier version of this paper, and Andrea Lyman and Vera M. Hesslinger for proofreading the manuscript. Last but not least, I am very indepted to Tomoki Nakaya, who has developed the original inference statistics of the overall models and the referring parameters and who helped me with reanalyzing these statistics as well as with ensuring the reliability of the used methods. Thank you!

## References
* Tobler, W. R. (1965). Computation of the corresponding of geographical patterns. Papers of the Regional Science Association, 15, 131-139.
* Tobler, W. R. (1966). Medieval distortions: Projections of ancient maps. Annals of the Association of American Geographers, 56(2), 351-360.
* Tobler, W. R. (1994). Bidimensional regression. Geographical Analysis, 26(3), 187-212.
* Friedman, A., & Kohler, B. (2003). Bidimensional regression: Assessing the configural similarity and accuracy of cognitive maps and other two-dimensional data sets. Psychological Methods, 8(4), 468-491.
* Nakaya, T. (1997). Statistical inferences in bidimensional regression models. Geographical Analysis, 29(2), 169-186.
* Waterman, S., & Gordon, D. (1984). A quantitative-comparative approach to analysis of distortion in mental maps. Professional Geographer, 36(3), 326-337.

## License
All code is licensed under the [GPL 3.0](https://opensource.org/licenses/GPL-3.0) license.
