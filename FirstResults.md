First Estimates
================
Lucas Ramalho Anderson
21/10/2020

``` r
library ( dplyr )
library ( ggplot2 )

tempDir = "/scratch/genevol/users/lucas/"
saveDir = "/raid/genevol/users/lucas/heritability/plots/"
```

Introduction
------------

Step 1
------

After the removal of monomorphisms, filtration of the desired samples and after obtaining the list of non-correlated snp's per chromosome (considering correlation value of \(\sqrt{0.1}\)), it is now desired to calculate the GRM matrix.

``` r
# Read file with all chromosomes
# allChrFile = SeqArray::seqOpen ( paste0 ( tempDir , "allChr.gds" ) )
# List of all genes of interest (after pruning)
listSnps = readRDS ( paste0 ( tempDir , "fullPrunedList.rds" ) )

# GRM - calculated as defined in CGTA
# grm_obj = SNPRelate::snpgdsGRM( allChrFile , snp.id = listGenes , method = "GCTA")

# Estimating through "gaston" package
altReadSnps = gaston::read.vcf( paste0 ( tempDir , "allChr.vcf.gz" ) )
```

    ## ped stats and snps stats have been set. 
    ## 'p' has been set. 
    ## 'mu' and 'sigma' have been set.

``` r
# setting "p" parameter - correction with mean "p" and std sqrt(2p(1-2p))
gaston::standardize( altReadSnps ) <- "p"
grm_matrix = gaston::as.matrix ( altReadSnps )
# grm_scaled = scale( grm_matrix , center = T , scale = T )
# grm_scaled = readRDS (paste0(tempDir , "scaledMatrixBk.rds"))


# manual_GRM = ( 1 / nrow ( grm_scaled ) ) * grm_scaled %*% t ( grm_scaled )
# GRM matrix calculation (GCTA)
grm_alt_p = gaston::GRM ( altReadSnps , which.snps = listSnps )
```

    ## Warning in which.snps & is.autosome(x@snps$chr): longer object length is not a
    ## multiple of shorter object length

``` r
# transform matrix into dataframe (3 columns - col1 = samples each row, col2 = samples each column ,  col3 = values for each pair)
dfGrm = reshape2::melt(grm_alt_p)

# indexing with numeric values each sample (columns and rows)
# dfGrm$sampLines = rep ( seq ( 1 , nrow ( grm_alt_p ) ) , nrow ( grm_alt_p ) ) 
# dfGrm$sampCols = sort ( rep ( seq ( 1 , nrow ( grm_alt_p ) ) , nrow ( grm_alt_p ) ) )


# To calculate the correlation between individuals, the calculation A_ij/sqrt(A_ii)sqrt(A_jj) will be done
# dataframe with only diag. values
dfGrmDiag = dfGrm[ dfGrm$Var1 == dfGrm$Var2,]
# sqrt of those values
dfGrmDiag = dfGrmDiag %>% mutate ( sqrtVal = sqrt ( value ) , sqrtVal2 = sqrt ( value ) )

# merging each A_ii for each row and col
dfGrmM = merge ( dfGrm , dfGrmDiag[ ,c ( "sqrtVal" , "Var1" ) ] , on = c ( "Var1" ) )
dfGrmM2 = merge ( dfGrmM , dfGrmDiag[ ,c ( "sqrtVal2" , "Var2" ) ] , on = c ( "Var2" ) )

# Calculating A_ij/(sqrt(A_ii)sqrt(A_jj))
dfGrmFinal = dfGrmM2 %>% mutate ( corrIndividuals = value / ( sqrtVal * sqrtVal2 ) ) %>% arrange ( Var1 , Var2 )

# plot heatmap - correlation between individuals
dfGrmFinal %>% ggplot( aes ( x = Var1 , y = Var2 , fill = corrIndividuals ) ) + 
geom_tile() +
theme( axis.text.x = element_text(angle = 90, hjust = 1) , text = element_text (size = 5) ) +
labs ( x = "Sample ID" , y = "Sample ID" )
```

![](FirstResults_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# It seems there are blocks with higher correlation between individuals between the samples

# Filter of all correlation values between individuals
dfUniqueCorr = dfGrmFinal %>% filter ( corrIndividuals < .9999 ) %>% distinct ( corrIndividuals , .keep_all = TRUE)

# Histogram and density of correlation values
dfUniqueCorr %>% ggplot ( aes ( x = corrIndividuals ) ) +
geom_histogram ( aes(y=..density..) , bins = 100 ) +
geom_density ( ) +
labs ( x = "Correlation between individuals" , y = "Density" , title = "Histogram of correlation between distinct individuals" )
```

![](FirstResults_files/figure-markdown_github/unnamed-chunk-2-2.png)

``` r
# The correlation blocks are bolder in this plot 

# Readind file with HLA expressions and ancestry information
hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   subject_id = col_character(),
    ##   continental_pop = col_character(),
    ##   population = col_character(),
    ##   sex = col_character(),
    ##   gene_name = col_character(),
    ##   NumReads = col_double(),
    ##   TPM = col_double()
    ## )

``` r
# Ancestry of all samples
ancestry = unique ( hlaExp[ , c ( "subject_id" , "continental_pop" )] )

# Merging ancestry info with correlation dataframe
check = merge ( dfUniqueCorr , ancestry , by.x = c ( "Var1" ) , by.y = c ( "subject_id" ) )
check2 = merge ( check , ancestry , by.x = c ( "Var2" ) , by.y = c ( "subject_id" ) )

tableAncestry = unique ( check[,c("continental_pop" , "Var1")] ) %>% select ( continental_pop ) %>% table() %>% as.data.frame ( ) %>%  mutate ( relFreq = paste0 ( 100 * round ( Freq / sum ( Freq ) , 4 ) , "%") ) %>% rename ( "Ancestry" = "." )

knitr::kable( tableAncestry )
```

| Ancestry |  Freq| relFreq |
|:---------|-----:|:--------|
| AFR      |    87| 19.59%  |
| EUR      |   357| 80.41%  |

``` r
# Approximately 20% of the 444 individuals are African, while the other 80% are European



# Checking the amount of comparisons between individuals with same ancestry and different ones
checkFin = check2 %>% mutate ( ancestries = ifelse ( continental_pop.x == continental_pop.y , continental_pop.x , "Diff" ) )


tableComparisons = table ( checkFin$ancestries ) %>% as.data.frame() %>% mutate ( freqRel = Freq/ sum ( Freq ) ) %>% rename ( "Ancestry" = "Var1" , "NumComparisons" = "Freq" )

knitr::kable ( tableComparisons )
```

| Ancestry |  NumComparisons|    freqRel|
|:---------|---------------:|----------:|
| AFR      |            3741|  0.0378682|
| Diff     |           31146|  0.3152748|
| EUR      |           63903|  0.6468570|

``` r
checkFin %>% ggplot ( aes ( x = corrIndividuals , fill =  ancestries ) ) +
geom_histogram ( aes(y=..density..) , bins = 10 ) +
geom_density ( ) +
facet_wrap ( ~ancestries ) +
theme(panel.spacing = unit (2, "lines") ) +
labs ( x = "Correlation between individuals" , y = "Density" , title = "Histogram of correlation between distinct individuals" )
```

![](FirstResults_files/figure-markdown_github/unnamed-chunk-2-3.png)

``` r
# display individuals with correlation greater than 10% in the sample
listGreatCorr = checkFin[ ( checkFin$corrIndividuals > .1 ) & ( checkFin$corrIndividuals < .999 ) , ] %>% distinct( corrIndividuals , .keep_all = TRUE)

listGreatCorr
```

    ##      Var2    Var1      value   sqrtVal  sqrtVal2 corrIndividuals
    ## 1 HG00120 HG00116 0.09008931 0.7819242 0.7784499       0.1480055
    ## 2 HG00240 HG00238 0.07855307 0.7790604 0.7979312       0.1263649
    ##   continental_pop.x continental_pop.y ancestries
    ## 1               EUR               EUR        EUR
    ## 2               EUR               EUR        EUR

``` r
grm = grm_alt_p
# rownames ( grm ) = altReadSnps
# colnames ( grm ) = altReadSnps$sample.id

eigenValuesGrm = eigen ( grm )
dfEigen = eigenValuesGrm$values %>% 
as.data.frame ( ) %>% 
mutate ( order = 1:n() ) %>%  
rename ( "Value" = '.' ) %>%  
mutate ( neg = ifelse ( Value < 0 , "Negative" , "Positive" ) )


dfEigen %>% ggplot ( aes ( x = order , y = Value , colour = neg )  ) + 
geom_point ( ) +
labs ( x = "Order" , y = "Eigen Value" , title = "Eigen values plot" , colour = "Sign" )
```

![](FirstResults_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
matrixCorrection = eigenValuesGrm$vectors %*% diag( eigenValuesGrm$values + abs ( min ( eigenValuesGrm$values ) ) ) %*% t ( eigenValuesGrm$vectors )
# rownames ( matrixCorrection ) = grm_obj$sample.id
# colnames ( matrixCorrection ) = grm_obj$sample.id

rownames ( matrixCorrection ) = rownames ( grm )
colnames ( matrixCorrection ) = colnames ( grm )


eigenCorr = eigen ( matrixCorrection )

dfEigenCorr = eigenCorr$values %>% 
  as.data.frame ( ) %>% 
  mutate ( order = 1:n() ) %>%  
  rename ( "Value" = '.' ) %>%  
  mutate ( neg = ifelse ( Value < 0 , "Negative" , "Positive" ) )


dfEigenCorr %>% ggplot ( aes ( x = order , y = Value , colour = neg )  ) + 
  geom_point ( ) +
  labs ( x = "Order" , y = "Eigen Value" , title = "Eigen values plot" , colour = "Sign" )
```

![](FirstResults_files/figure-markdown_github/unnamed-chunk-3-2.png)

``` r
expressionInterest = hlaExp %>% filter ( subject_id %in% colnames ( grm ) )

mainInfo = expressionInterest %>% distinct( subject_id , continental_pop ,population )
numEigen = 2
print ( paste0 ( "Total variation explained by the first ", numEigen , " eigen values: " , 100*round ( sum ( eigenCorr$values[1:numEigen] )/ sum ( eigenCorr$values ) , 4 ) , "%" ) )
```

    ## [1] "Total variation explained by the first 2 eigen values: 5.61%"

``` r
vectors_ = eigenCorr$vectors[,1:numEigen]
calcScores = matrixCorrection %*% vectors_ %>% 
  as.data.frame() %>% 
  rename ( "PC1" = "V1" , "PC2" = "V2" ) %>% 
  mutate ( subject_id = rownames ( matrixCorrection ) )
pcaPlot = merge ( mainInfo , calcScores )

pcaPlot %>% ggplot ( aes ( x = PC1 , y = PC2  , colour = continental_pop ) ) +
  geom_point ( ) +
  labs ( title = "PCA plot (first 2 dimensions)" , colour = "Origin" )
```

![](FirstResults_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
  # geom_text ( )
```

``` r
simpleModels = function ( exp_ , df ){
  
  dfFilter = df  %>% filter ( gene_name == exp_ )
  
  fixed0 = lm ( TPM ~ 1 , data = dfFilter )
  sum0 = summary ( fixed0 )
  fixedEffectSigma = sum0$sigma^2

  mixedModel = coxme::lmekin( dfFilter$TPM ~ 1 + (1|dfFilter$subject_id) , data=dfFilter, varlist=list(matrixCorrection), vinit=2)
  
  mixedEffectSigma = mixedModel$sigma^2
  sigmaA = as.numeric(mixedModel$vcoef)
  
  # comparison = mixedEffectSigma/fixedEffectSigma
  
  # h = sigmaA / ( sigmaA + mixedEffectSigma)
  
  
  modelExpanded = coxme::lmekin( dfFilter$TPM ~ 1 + dfFilter$PC1 + dfFilter$PC2 + (1|dfFilter$subject_id), data=dfFilter, varlist=list(matrixCorrection), vinit=2)

  mixedEffectSigmaExp <- modelExpanded$sigma^2
  # comparisonExp = modelExpanded/fixedEffectSigma
  sigmaAExp = as.numeric(modelExpanded$vcoef)
  
  # hExp = sigmaAExp / (sigmaAExp + mixedEffectSigmaExp )
  
  return ( c ( exp_ , fixedEffectSigma , mixedEffectSigma , sigmaA , mixedEffectSigmaExp , sigmaAExp ) )
  
}

listNames = unique ( expressionInterest$gene_name )
modelDf = merge ( expressionInterest , calcScores )

requiredInfo = NULL
for ( name_ in listNames ){
  
  requiredInfo = rbind ( requiredInfo , simpleModels ( exp_ = name_ ,df = modelDf ) )
  
}

finalDf = requiredInfo %>% as.data.frame ( ) %>% rename ("Gene" = "V1" , 
                                                          "fixedSigma" = "V2" ,
                                                          "residualMixedSigma" = "V3" , 
                                                          "randomEffectSigma" = "V4",
                                                          "residualMixedSigmaExp" = "V5" ,
                                                          "randomEffectSigmaExp" = "V6") %>%
mutate ( fixedSigma = as.numeric ( as.character ( fixedSigma ) ) ,
    residualMixedSigma = as.numeric ( as.character (residualMixedSigma)) ,
    randomEffectSigma =  as.numeric ( as.character (randomEffectSigma)) ,
    residualMixedSigmaExp = as.numeric ( as.character (residualMixedSigmaExp)),
    randomEffectSigmaExp = as.numeric ( as.character (randomEffectSigmaExp))
    ) %>%
  mutate ( comparisonNull = residualMixedSigma/fixedSigma , 
           comparisonNullExp = residualMixedSigmaExp/fixedSigma ,
           hSimple = randomEffectSigma / ( randomEffectSigma + residualMixedSigma ) ,
           hExpanded = randomEffectSigmaExp / ( randomEffectSigmaExp + residualMixedSigmaExp ))

finalDf %>% knitr::kable()
```

<table>
<colgroup>
<col width="6%" />
<col width="7%" />
<col width="12%" />
<col width="11%" />
<col width="14%" />
<col width="13%" />
<col width="9%" />
<col width="11%" />
<col width="6%" />
<col width="6%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Gene</th>
<th align="right">fixedSigma</th>
<th align="right">residualMixedSigma</th>
<th align="right">randomEffectSigma</th>
<th align="right">residualMixedSigmaExp</th>
<th align="right">randomEffectSigmaExp</th>
<th align="right">comparisonNull</th>
<th align="right">comparisonNullExp</th>
<th align="right">hSimple</th>
<th align="right">hExpanded</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">HLA-A</td>
<td align="right">180900.15</td>
<td align="right">126067.66</td>
<td align="right">51701.260</td>
<td align="right">132658.36</td>
<td align="right">43406.933974</td>
<td align="right">0.6968908</td>
<td align="right">0.7333237</td>
<td align="right">0.2908341</td>
<td align="right">0.2465388</td>
</tr>
<tr class="even">
<td align="left">HLA-B</td>
<td align="right">526383.05</td>
<td align="right">471645.07</td>
<td align="right">48733.998</td>
<td align="right">511004.80</td>
<td align="right">362.711598</td>
<td align="right">0.8960111</td>
<td align="right">0.9707850</td>
<td align="right">0.0936510</td>
<td align="right">0.0007093</td>
</tr>
<tr class="odd">
<td align="left">HLA-C</td>
<td align="right">127438.40</td>
<td align="right">85467.99</td>
<td align="right">40356.415</td>
<td align="right">88511.16</td>
<td align="right">36480.976408</td>
<td align="right">0.6706612</td>
<td align="right">0.6945408</td>
<td align="right">0.3207360</td>
<td align="right">0.2918662</td>
</tr>
<tr class="even">
<td align="left">HLA-DPA1</td>
<td align="right">31587.29</td>
<td align="right">29133.91</td>
<td align="right">2064.005</td>
<td align="right">30542.58</td>
<td align="right">5.088393</td>
<td align="right">0.9223301</td>
<td align="right">0.9669265</td>
<td align="right">0.0661584</td>
<td align="right">0.0001666</td>
</tr>
<tr class="odd">
<td align="left">HLA-DPB1</td>
<td align="right">37625.10</td>
<td align="right">25389.65</td>
<td align="right">9592.266</td>
<td align="right">33403.72</td>
<td align="right">9.749311</td>
<td align="right">0.6748061</td>
<td align="right">0.8878041</td>
<td align="right">0.2742065</td>
<td align="right">0.0002918</td>
</tr>
<tr class="even">
<td align="left">HLA-DQA1</td>
<td align="right">51654.46</td>
<td align="right">41390.59</td>
<td align="right">8644.455</td>
<td align="right">48501.41</td>
<td align="right">28.894018</td>
<td align="right">0.8012975</td>
<td align="right">0.9389589</td>
<td align="right">0.1727680</td>
<td align="right">0.0005954</td>
</tr>
<tr class="odd">
<td align="left">HLA-DQB1</td>
<td align="right">38524.93</td>
<td align="right">36032.89</td>
<td align="right">2319.584</td>
<td align="right">38005.48</td>
<td align="right">27.989865</td>
<td align="right">0.9353135</td>
<td align="right">0.9865164</td>
<td align="right">0.0604807</td>
<td align="right">0.0007359</td>
</tr>
<tr class="even">
<td align="left">HLA-DRA</td>
<td align="right">390199.52</td>
<td align="right">371872.42</td>
<td align="right">15914.581</td>
<td align="right">381119.41</td>
<td align="right">32.216643</td>
<td align="right">0.9530315</td>
<td align="right">0.9767296</td>
<td align="right">0.0410395</td>
<td align="right">0.0000845</td>
</tr>
<tr class="odd">
<td align="left">HLA-DRB1</td>
<td align="right">148507.21</td>
<td align="right">65927.22</td>
<td align="right">69784.460</td>
<td align="right">83398.53</td>
<td align="right">46952.682502</td>
<td align="right">0.4439328</td>
<td align="right">0.5615790</td>
<td align="right">0.5142112</td>
<td align="right">0.3602014</td>
</tr>
</tbody>
</table>
