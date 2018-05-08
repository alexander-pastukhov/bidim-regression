## Test environments
* local Windows 10 install, R 3.5.0
* ubuntu 16.04 (on VirtualBox), R 3.4.4

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

1)  * checking R code for possible problems ... NOTE
anova.lm2: no visible binding for global variable 'is'

However, this depedency was not included as it is not required and generates an error otherwise: "Namespace dependency not required: 'methods'".
