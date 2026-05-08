---
title: "Moreton Bay Seagrass Resilience Modelling"
subtitle: "DRFA 21/22 Biodiversity Conservation Marine Final Report Project"
author: "Michael Dumelle and Erin E. Peterson"
format:
  pdf:
    toc: true
    cite-method: natbib
    biblio-style: plainnat
    number-sections: true
    keep-tex: true
    keep-md: true
  html:
    toc: true
    cite-method: natbib
    biblio-style: plainnat
    number-sections: true
    keep-tex: true
    keep-md: true
execute: 
  echo: false
  message: false
  warning: false
  error: true
  cache: false
embed-resources: true
editor_options: 
  chunk_output_type: console
bibliography: references.bib
---


::: {.cell}

:::


# Introduction

Over the past few decades, various scientific studies have informed our current understanding of trends, behaviors, and flood responses of seagrass cover in Moreton Bay [@abal1996seagrass; @udy1997growth; @samper2016organic; @ovsyanikova2026spatial]  as well as Australia more widely [@kilminster2015unravelling; @gilby2018seagrass; @o2018seagrass]. Remote sensing has offered new capabilities for seagrass cover mapping [@roelfsema2009integrated], while innovative statistical approaches have offered novel methods for studying resilience at site-specific scales [@gilby2023drivers]. Our work adds to this body of literature by applying spatially explicit statistical models [@cressie1993statistics] to study total percent seagrass cover (i.e., seagrass cover) drivers, resilience, and recovery in Moreton Bay and Wide Bay (Dumelle and Peterson 2026b). Specifically, we provide and contextualize answers to four specific research questions, applicable to Moreton Bay:

1. Can we identify water quality covariates that act as triggers for seagrass cover management responses? Are these triggers different based on whether they result from a specific event (e.g., flood) or ambient conditions?
2. Did seagrass cover improve after the flood event(s)? If so, how and where?
3. Can we identify site-specific resilience indicators, which describe where recovery is likely to be enhanced or impeded? If so, what factors contribute to seagrass resilience?


# Seagrass Cover Model
## Methods

### Data


::: {.cell}

:::


The seagrass cover dataset was provided by Science Under Sail Australia (SUSA; https://www.scienceundersail.com.au/). It included 7657 observations across three financial years (2020-2021, 2022-2023, and 2023-2024), where water quality and geographic covariates were also available. Of the 7657 observations, 1995 were collected in 2020-2021 (pre-flood), 3329 in 2022-2023 (1-year post-flood), and 2333 in 2023-2024 (2-years post-flood) (@tbl-n_fy_mb). 


::: {#tbl-n_fy_mb .cell tbl-cap='The sample size and percent of total sample size for total percent seagrass cover observations in Moreton Bay, by financial year.'}
::: {.cell-output-display}


|Financial Year | Sample Size|Percent of Total |
|:--------------|-----------:|:----------------|
|2020-2021      |       1,995|26.1%            |
|2022-2023      |       3,329|43.5%            |
|2023-2024      |       2,333|30.5%            |
|Total          |       7,657|100.0%           |


:::
:::


Seagrass cover is a percentage that ranges from 0 (no seagrass) to 100 (complete seagrass cover).  shows total cover by financial year throughout the Bay. Total cover tends to be highest near the eastern shoreline of the Bay and lowest in the middle of the Bay (@fig-cover_fy_mb). Note that 2022-2023 (post-flood) had the the largest number of seagrass cover samples, as well as the only samples collected in the middle of the Bay.


::: {.cell}
::: {.cell-output-display}
![Total percent seagrass cover in Moreton Bay by financial year.](MBreport_ep_files/figure-pdf/fig-cover_fy_mb-1.pdf){#fig-cover_fy_mb}
:::
:::


EcoFutures Consulting used eReefs data [@steven2019ereefs] to generate a broad suite of ambient and flood-related turbidity, salinity, and temperature variables at biologically and ecologically relevent temporal windows (@tbl-wq_vars_mb). Ambient water quality variables were based on the mean, minimum, and maximum values for three periods: 7, 30, and 90 days prior to sampling. Ambient variables were also calculated representing the frequency of days water quality thresholds were exceeded (above or below) and the maximum spell (i.e. duration) of days those thresholds were continuously exceeded during the time periods. Similar water quality variables representing flood-event conditions were also calculated for a 30-day flood period. 



::: {#tbl-wq_vars_mb .cell tbl-cap='Moreton Bay water quality variables considered in the exploratory analyses, including the variable name and measurement unit (in parentheses), the type of water-quality variable (Ambient, Flood), the time period in days, statistic calculated, and the threshold value used (if any).'}
::: {.cell-output-display}


|Variable            |Type    |Period    |Statistic |Threshold |
|:-------------------|:-------|:---------|:---------|:---------|
|Turbidity (NTU)     |Ambient |30, 7, 90 |Mean      |          |
|Turbidity (NTU)     |Ambient |30, 7, 90 |Frequency |> 3       |
|Turbidity (NTU)     |Ambient |30, 7, 90 |Frequency |> 5       |
|Turbidity (NTU)     |Ambient |30, 7, 90 |Frequency |> 6       |
|Salinity (PSU)      |Ambient |30, 7, 90 |Frequency |< 20      |
|Salinity (PSU)      |Ambient |30, 7, 90 |Frequency |< 25      |
|Salinity (PSU)      |Ambient |30, 7, 90 |Frequency |< 30      |
|Temperature (deg C) |Ambient |30, 7, 90 |Mean      |          |
|Turbidity (NTU)     |Flood   |30        |Mean      |          |
|Turbidity (NTU)     |Flood   |30        |Max       |          |
|Turbidity (NTU)     |Flood   |30        |Frequency |> 3       |
|Turbidity (NTU)     |Flood   |30        |Frequency |> 5       |
|Turbidity (NTU)     |Flood   |30        |Frequency |> 6       |
|Turbidity (NTU)     |Flood   |30        |Duration  |> 3       |
|Turbidity (NTU)     |Flood   |30        |Duration  |> 5       |
|Turbidity (NTU)     |Flood   |30        |Duration  |> 6       |
|Salinity (PSU)      |Flood   |30        |Min       |          |
|Salinity (PSU)      |Flood   |30        |Mean      |          |
|Salinity (PSU)      |Flood   |30        |Frequency |< 20      |
|Salinity (PSU)      |Flood   |30        |Frequency |< 25      |
|Salinity (PSU)      |Flood   |30        |Frequency |< 30      |
|Salinity (PSU)      |Flood   |30        |Duration  |< 20      |
|Salinity (PSU)      |Flood   |30        |Duration  |< 25      |
|Salinity (PSU)      |Flood   |30        |Duration  |< 30      |


:::
:::


Additional covariates were also provided by EcoFutures Consulting, which are used to represent georaphic features, proximity to rivers, and depth, as well as the years since the flood and the seagrass growing season. These variables have the following structure:

* Sediment Type includes six categorical levels: Sand, Sandy Mud, Muddy Sand, Mud, Rock, and Other.
* Tidal Range has three categorical levels: Shallow Subtidal, Intertidal, and Deep Subtidal.
* Years Since Flood is numeric ranging from zero (pre flood) to 1.93 years post-flood (February, 2024). 
* Nearest River includes six categorical levels: Brisbane River, Caboolture River, Logan River, Ningi Creek, Pimpama Coomera and Nerang Rivers, and Pine River.
* Distance to (Nearest River) Outlet is numeric and ranges from 0.24 to 34.97 kilometers.
* Growing Season is a categorical variable with two levels: Growing (September - December) and Not Growing (January - August). 
* EHMP subzone has eleven categorical levels: Bramble Bay, Broadwater, Central Bay, Deception Bay, Eastern Banks, Moreton Bay Open Coastal Waters, North-Eastern Moreton Bay, Northern Banks Open Coastal Waters, South-Eastern Moreton Bay, Southern Bay, and Waterloo Bay (@fig-ehmp_subzone). 


::: {.cell}

:::




::: {.cell}
::: {.cell-output-display}
![Observations and their EHMP Subzones in Moreton Bay.](MBreport_ep_files/figure-pdf/fig-ehmp_subzone-1.pdf){#fig-ehmp_subzone}
:::
:::




### Exploratory Analyses
Extensive exploratory analyses were undertaken to 1) understand the data structure, 2) identify outliers/errors, 3) examine correlations between potential covariates and the response (i.e. seagrass cover), as well as multicollinearity in the covariates, 4) evaluate model assumptions, and 5) visualise spatial and temporal patterns in the the data. 

Results revealed substantial collinearity among ambient and flood-related water-quality variables, which was expected given their nested temporal windows and threshold values. We retained only one variable from each correlated group. We also removed variables showing no meaningful association with seagrass cover or where the relationship lacked biological plausibility. Geographic variables exhibited less collinearity than water quality variables; however, we removed the distance to the nearest river outlet and the nearest river variables due to weak associations with seagrass cover. Growing season was excluded because a simple monthly proxy inadequately captures seasonal seagrass growth dynamics.

The final set of potential covariates considered in the modelling included mean turbidity in the 90 days prior to sampling, mean turbidity during the flood event, mean temperature in the 90 days prior to sampling, sediment type, tidal range, years since flood, and EHMP zone. 

### Statistical Modeling {#sec-mb_model}

We considered multiple modeling approaches in this study, which each accounted for the data distribution, nonlinear relationships between seagrass cover and covariates, and spatial dependence inherent in biological data collected through space and time to different degrees (Appendix 1). We used the full set of potential covariates to fit a nonspatial linear model, spatial linear model (i.e., spatial model), spatial beta regression model, spatial gamma regression model, and a random forest model. 

There were several reasons for fitting these initial models. The nonspatial linear model provided a useful baseline to compare the spatial linear model against, elucidating the impact of spatial dependence on model fit. The spatial beta regression model is a generalized linear model used to characterize proportion data ranging between zero and one. To adhere to this restriction, total seagrass cover was 1) scaled to a proportion and 2) if exactly zero or one, set to 0.01 or 0.99, respectively. The spatial gamma regression model is a generalized linear model used to characterize skewed positive data. To adhere to this restriction, a value of 0.01 was added to seagrass cover observations that were equal to zero. Finally, the random forest model provided a useful assessment of a machine learning approach, which requires fewer distributional assumptions than the aforementioned statistical approaches and more naturally describes nonlinear relationships between variables (Appendix 1). We also characterized random forest effects using partial dependence plots [@greenwell2017pdp], which isolate the marginal contribution of a variable while holding other variables at their mean values. Leave-one-out cross-validation predictions were then used to calculate the mean bias (MBias) and root mean square prediction error (RMSPE) for each model. 


::: {.cell}

:::


The results of this initial model comparison showed that the spatial linear model had the lowest RMSPE, followed by the random forest, generalized linear models, and the nonspatial linear model (@tbl-mb_cv). The nonspatial and spatial linear models also had much lower mean bias values than the random forest model, while the two generalized linear models produced unacceptably biased predictions. Nevertheless, all five approaches identified similar relationships between seagrass cover and the covariates, in the directions expected. Random forest partial dependence plots [@greenwell2017pdp] suggested that linear effects of the model variables was reasonable. Given these results, we chose to use a spatial linear modeling framework for the Moreton Bay seagrass cover model.


::: {#tbl-mb_cv .cell tbl-cap='Moreton Bay out of sample prediction metrics for the various approaches. Mean bias (MBias) tended to be small relative to root-mean-squared-prediction error (RMSPE).'}
::: {.cell-output-display}


|Approach    | MBias| RMSPE|
|:-----------|-----:|-----:|
|SPLM        | -0.02| 18.41|
|RF          | -0.11| 18.84|
|SPGLM-Beta  | -2.76| 19.60|
|SPGLM-Gamma |  5.42| 20.07|
|LM          |  0.00| 22.54|


:::
:::



::: {.cell}
::: {.cell-output-display}
![Random forest partial dependence for seagrass cover as a function of mean turbidity in the 90 days prior to sampling.](MBreport_ep_files/figure-pdf/fig-rf_turb-1.pdf){#fig-rf_turb width=50%}
:::
:::


A suite of spatial linear models (Appendix 1, @eq-splm1), with exponential spatial covariance functions, were fit using the `splm()` function in `spmodel`. Parameters were estimated using restricted maximum likelihood. Models were compared using XXX.

The final total seagrass cover model for Moreton Bay included the following covariates: tidal range, sediment type, 90-day mean temperature (prior to sampling), pre- versus post-flood, years since the flood, 90 day turbidity (prior to sampling), and mean turbidity during the flood period scaled by time since the flood squared (@tbl-mb_tc_splm_vars). Notice that we scaled mean turbidity during the flood period so that the effect of turbidity during the flood event was allowed to decay with time since the flood. Intuitively, this means that as the time since the flood increases, the impact of turbidity during the flood event decreases. We also considered scaling by time since flood (not squared), which produced similar results. A random effect for EHMP subzone was also included in the final model. 


::: {#tbl-mb_tc_splm_vars .cell tbl-cap='Covariates in the final Moreton Bay total cover model, their type (numeric or categorical), and if categorical, the number of unique levels.'}
::: {.cell-output-display}


|Variable                                 |Type        |No. Levels |
|:----------------------------------------|:-----------|:----------|
|90 Day Temperature                       |Numeric     |NA         |
|90 Day Turbidity                         |Numeric     |NA         |
|Mean Flood Turbidity/Years Since Flood^2 |Numeric     |NA         |
|Sediment Type                            |Categorical |6          |
|Tidal Range                              |Categorical |3          |


:::
:::



::: {.cell}

:::



::: {#tbl-mb-covariates .cell tbl-cap='ANOVA table of covariates in the Moreton Bay water quality and geographic model.'}
::: {.cell-output-display}


|Effects                                  | df|  Chi-sq| p.value|
|:----------------------------------------|--:|-------:|-------:|
|Tidal Range                              |  2| 608.900|   0.000|
|Sediment Type                            |  5| 113.720|   0.000|
|90 Day Temperature                       |  1|  16.967|   0.000|
|Mean Flood Turbidity/Years Since Flood^2 |  1|  11.966|   0.001|
|90 Day Turbidity                         |  1|   7.854|   0.005|


:::
:::


@tbl-mb-covariates shows the results from an analysis of variance applied to the fitted model. All covariates had large test statistics and small $p$-values ($p$-value $<$ 0.01), suggesting that they had a significant relationship with total seagrass cover; here are some relevant takeaways:

* Tidal Range: The estimated (marginal) average total cover is largest in the intertidal (26.89%), followed by shallow subtidal (16.51%) and deep subtidal (4.88%).
* Sediment Type: The estimated (marginal) average total cover is largest in the other category (25.6%), followed by muddy sand (19.1%), sandy mud (17.8%), mud (14.3%), sand (14.2%), and rock (5.62%). 
* 90 Day Temperature: A one-degree (C) increase in 90-day water temperature prior to sampling is associated with a decrease in average seagrass cover by 0.88%.
* Mean Flood Turbidity/(Years Since Flood)$^2$: A one-unit increase in mean turbidity during the flood is associated with a decrease in average seagrass cover of 2.05%/(Years Since Flood)$^2$ (@fig-mb_turb_90d_fld).
* 90 Day Turbidity: A one-unit increase in 90-day turbidity (prior to sampling) is associated with a decrease in average seagrass cover by 7.18% (@fig-mb_turb_90d_fld).

Note that the estimated marginal means for tidal range and sediment type were calculated using the `emmeans` **R** package [@lenth2025emmeans]. 


::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![Average total seagrass percent cover as a function of 90-day mean turbidity (left) and flood turbidity (right). Estimates (and confidence intervals) were evaluated at one year post flood, at the means of all other numeric variables, and for the intertidal depth range and sand sediment type (the same trends occur for both variables among the remaining species and subzones).](MBreport_ep_files/figure-pdf/fig-mb_turb_90d_fld-1.pdf){#fig-mb_turb_90d_fld}
:::
:::



We examined model diagnostics to better understand the overall model fit. First, we partitioned the model into different components of variability (@tbl-mb_varcomp) and found that 11.2% is attributable to the covariates (this quantity is sometimes called the Pseudo R-squared), 31.6% to the spatial random error, 7.3% to the EHMP subzone, and the remaining variability is independent random error. Second, we examined the decay behaviour of spatial autocorrelation. The estimated exponential spatial range ($\phi$) is 745.5 meters, which suggests two total seagrass cover observations are approximately (spatially) uncorrelated at a distance of $3\phi \approx 2$km, or two kilometers. Third, we analyzed the standardized model residuals, which have been decorrelated, to assess model assumptions and identify influential observations that had a sizable impact on model fit. Generally, the standardized residuals are centered around zero, though there is a slight right skew suggesting some observed total cover values are larger than expected (according to the model). Observations are referred to as "influential" if they highly influence model fit. Influence is often assessed a Cook's distance threshold value of one [@cook1982residuals]. No observations exceed this threshold.


::: {#tbl-mb_varcomp .cell tbl-cap='Variance components in the Moreton Bay seagrass model.'}
::: {.cell-output-display}


|Variance Component |Percent (of Total Variability) | Parameter Value|
|:------------------|:------------------------------|---------------:|
|Covariates         |11.2%                          |              NA|
|Spatial Effect     |31.6%                          |           182.1|
|Independent Error  |49.9%                          |           287.5|
|EHMP Subzone       |7.3%                           |           745.5|


:::
:::


### Block Kriging

A tremendous benefit of using the spatial linear model is that block kriging [@ver2002sampling; @cressie2006block] can be used to predict the mean total seagrass cover bay-wide for any combination of covariates and geographic regions, even if these variables were not directly used in the model (see Appendix 1 for additional details). Block kriging was undertaken, by financial year, for EHMP subzone, Marine Park zone, tidal range and for the whole of Moreton Bay.  


@fig-bk_mb_ehmp and @fig-bk_mb_zone show there are several EHMP subzones and Marine Park zones that have recovered to pre-flood levels (e.g., Deception Bay, Marine National Park Zone), while others have not (e.g., Southern Bay, General Use Zone). When we examine recovery by tidal range (@fig-bk_mb_tidal), we see that seagrass cover is recovering less effectively in the intertidal range than the shallow or deep subtidal ranges. When we examine recovery for the whole of Moreton Bay (@fig-bk_mb_all), the model suggests that seagrass cover has improved since the 2022 flood. 


::: {.cell}
::: {.cell-output-display}
![Average total cover predictions by EHMP subzone. Red dots indicate the mean from the raw data. Black dots indicate the predicted mean in the region alongside a prediction interval.](MBreport_ep_files/figure-pdf/fig-bk_mb_ehmp-1.pdf){#fig-bk_mb_ehmp}
:::
:::



::: {.cell}
::: {.cell-output-display}
![Average total cover predictions by conservation zone. Red dots indicate the mean from the raw data. Black dots indicate the predicted mean in the region alongside a prediction interval.](MBreport_ep_files/figure-pdf/fig-bk_mb_zone-1.pdf){#fig-bk_mb_zone}
:::
:::



::: {.cell}
::: {.cell-output-display}
![Average total cover predictions by tidal range. Red dots indicate the mean from the raw data. Black dots indicate the predicted mean in the region alongside a prediction interval.](MBreport_ep_files/figure-pdf/fig-bk_mb_tidal-1.pdf){#fig-bk_mb_tidal}
:::
:::



::: {.cell}
::: {.cell-output-display}
![Average total cover predictions in Moreton Bay. Red dots indicate the mean from the raw data. Black dots indicate the predicted mean in the region alongside a prediction interval.](MBreport_ep_files/figure-pdf/fig-bk_mb_all-1.pdf){#fig-bk_mb_all}
:::
:::



# Resilience Model {#sec-mb_resilience}

The seagrass cover model quantifies bay-wide effects of influential drivers on seagrass cover and enables assessment of recovery trajectories at management-relevant scales. However, it does not quantify site-specific resilience. @gilby2023drivers propose identifying hotspots (i.e. brightspots, sites performing better than expected) and coldspots (sites performing worse than expected) based on model residuals exceeding confidence interval bounds. They argue that hotspots and coldspots offer significant utility for understanding drivers of site-specific resilience and informing management decisions.

Building on @gilby2023drivers, we use the raw seagrass model residuals (@fig-resid_fy_mb) as a response variable in a second spatial linear model used to study site-specific resilience in Moreton Bay. Model covariates included variables describing seagrass diversity and dominant species, and bathymetry, as well as proximity to and area of coral, seagrass, mangrove, and gastropod habitat (@tbl-mb-res-covariates). EHMP subzone was again included as a random effect. We used B-spline basis functions [@deboor1978practical; @eilers1996flexible] to describe nonlinear relationships between the response and potential covariates. The degree of these functions was informed by random forest partial dependence plots derived from a random forest model fitted to the seagrass model residuals.  



::: {.cell}
::: {.cell-output-display}
![Seagrass model residuals for resilience analysis in Moreton Bay, by financial year.](MBreport_ep_files/figure-pdf/fig-resid_fy_mb-1.pdf){#fig-resid_fy_mb}
:::
:::






::: {#tbl-res_vars_mb .cell tbl-cap='Moreton Bay resilience variables considered for modeling.'}
::: {.cell-output-display}


|Variable                      |Effect |Included | Nonlinearity Degree|
|:-----------------------------|:------|:--------|-------------------:|
|Coral Habitat Total Area      |Fixed  |Yes      |                   2|
|Seagrass Habitat Total Area   |Fixed  |Yes      |                   2|
|Mangrove Habitat Total Area   |Fixed  |Yes      |                   2|
|Gastropod Habitat Total Area  |Fixed  |Yes      |                   2|
|Distance to Coral Habitat     |Fixed  |Yes      |                   3|
|Distance to Seagrass Habitat  |Fixed  |Yes      |                   2|
|Distance to Mangrove Habitat  |Fixed  |Yes      |                   2|
|Distance to Gastropod Habitat |Fixed  |Yes      |                   2|
|Species Richness              |Fixed  |Yes      |                   1|
|Habitat Count                 |Fixed  |Yes      |                   2|
|Dominant Species              |Fixed  |Yes      |                  NA|
|Bathymetry Mean Depth         |Fixed  |Yes      |                   2|
|Bathymetry Mean Curvature     |Fixed  |Yes      |                   2|
|Bathymetry Mean Slope         |Fixed  |Yes      |                   2|
|EHMP Subzone                  |Random |Yes      |                  NA|


:::
:::



@gilby2023drivers primarily focussed on specific hotspots, linking environmental variables to those sites. Here, our aim is to characterize the most important drivers of resilience, on average. Because we were interested in hypotheses testing rather than parsimony, we included all potential covariates as fixed effects in the model.  (@tbl-res_vars_mb), which have the following structure:


* The "Coral Habitat Total Area" is numeric from 0 msq to 616,180 msq.
* The "Seagrass Habitat Total Area" is numeric from 0 msq to 773,161 msq.
* The "Mangrove Habitat Total Area" is numeric from 0 msq to 561,317 msq.
* The "Gastropod Habitat Total Area" is numeric from 0 msq to 114,795 msq.
* The "Distance to Coral Habitat" is numeric from 0 m to 35,387 m.
* The "Distance to Seagrass Habitat" is numeric from 0 m to 11,382 m.
* The "Distance to Mangrove Habitat" is numeric from 0 m to 27,617 m.
* The "Distance to Gastropod Habitat" is numeric from 0 m to 32,034 m.
* The "Species Richness" is numeric (integer) from 0 unique species to 4 unique species.
* The "Habitat Count" is numeric (integer) from 0 unique habitats to 4 unique habitats.
* The "Dominant Species" has seven levels: "Absent", "CS", "HD", "HO", "HS", "SI", and "ZMHU".
* The "Bathymetry Mean Depth" is numeric from -35.95 m to 0.79 m.
* The "Bathymetry Mean Curvature" is numeric from -0.002 degrees to 0.0001 degrees.
* The "Bathymetry Mean Slope" is numeric from 0.009 degrees to 4.27 degrees.
* The "EHMP subzone" variable has eleven levels: Bramble Bay, Broadwater, Central Bay, Deception Bay, Eastern Banks, Moreton Bay Open Coastal Waters, North-Eastern Moreton Bay, Northern Banks Open Coastal Waters, South-Eastern Moreton Bay, Southern Bay, and Waterloo Bay (@fig-ehmp_subzone). 

## Model Fit


::: {#tbl-mb-res-covariates .cell tbl-cap='ANOVA table of covariates in the Moreton Bay resilience model.'}
::: {.cell-output-display}


|Effects                       | df|   Chi-sq| p.value|
|:-----------------------------|--:|--------:|-------:|
|Dominant Species              |  6| 1637.726|   0.000|
|Species Richness              |  1|   22.188|   0.000|
|Bathymetry Mean Depth         |  2|   21.876|   0.000|
|Bathymetry Mean Curvature     |  2|    8.999|   0.011|
|Gastropod Habitat Total Area  |  2|    5.155|   0.076|
|Bathymetry Mean Slope         |  2|    4.768|   0.092|
|Coral Habitat Total Area      |  2|    2.845|   0.241|
|Seagrass Habitat Total Area   |  2|    1.971|   0.373|
|Distance to Mangrove Habitat  |  2|    1.871|   0.392|
|Distance to Seagrass Habitat  |  2|    1.836|   0.399|
|Habitat Count                 |  2|    0.722|   0.697|
|Distance to Gastropod Habitat |  2|    0.613|   0.736|
|Distance to Mangrove Habitat  |  2|    0.536|   0.765|
|Distance to Coral Habitat     |  3|    0.718|   0.869|


:::
:::


From @tbl-mb-res-covariates, we find strong evidence ($p$-value < 0.001) that dominant species, species richness, and bathymetry mean depth are associated with resilience: 

* Dominant Species: The estimated (marginal) average resilience score is largest for CS (38.4), followed by ZMHU (18.8), HS (10.8), SI (10.7), HO (3.7), and HD (-3.1).
* Species Richness: An increase of species richness by one species was associated with an average increase in resilience by approximately two points (@fig-mb_res_rich_bath).
* Bathymetry Mean Depth: Deeper sites tended to be more resilient than shallower sites, on average (@fig-mb_res_rich_bath). 

We found moderate evidence ($p$-value < 0.1) that bathymetry mean curvature, gastropod habitat total area, and bathymetry mean slope are associated with resilience: 

* Bathymetry Mean Curvature: Sites with more curvature tended to be less resilient than sites with less curvature.
* Gastropod Habitat Total Area: Sites with more gastropod habitat area tended to be more resilient than sites with less gastropod habitat area.
* Bathymetry Mean Slope: Sites with less extreme slopes (between zero and one) tended to have increasing resilience with slope, while sites with more extreme slopes (greater than one) tended to have decreasing resilience with slope.

We found little evidence ($p$-value > 0.1) the remaining variables were associated with resilience, after controlling for other modeled variables.



::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![Average species richness (left) and bathymetry mean depth (right) effects on resilience. Predictions (and confidence intervals) were evaluated at the means of all other numeric variables and for the ZMHU species in the Southern Bay (similar trends occur for both variables among the remaining species and subzones).](MBreport_ep_files/figure-pdf/fig-mb_res_rich_bath-1.pdf){#fig-mb_res_rich_bath}
:::
:::


## Alternative Approaches

We also studied resilience using a random forest fit to the same residuals. Largely, we found similar results as the spatial linear model. Notably, for the random forest:

* Dominant species was most important, with a similar ranking of resilient species as for the spatial linear model approach (e.g., CS, ZMHU most resilient and Absent least resilient).
* Species richness was second-most important, with a general increasing trend in resilience with richness.
* Distance to gastropod habitat was third-most important, with a decreasing trend at small distances (less than 15,000 meters) and an increasing trend at large distances (greater than 15,000 meters).
* Bathymetry mean depth was fourth-most important, with a general decreasing trend with decreasing depth (@fig-rf_bathy), similar to the spatial linear model.
* Distance to gastropod habitat, distance to mangrove habitat, distance to coral habitat, and seagrass area had similar importance as the remaining bathymetry variables.
* The remaining variables (distance to seagrass habitat, mangrove area, EHMP subzone, habitat count, coral area, and gastropod area) were least important.


::: {.cell}
::: {.cell-output-display}
![Random forest partial dependence for resilience as a function of bathymetry mean depth.](MBreport_ep_files/figure-pdf/fig-rf_bathy-1.pdf){#fig-rf_bathy width=50%}
:::
:::



# Limitations and Future Suggestions

While the spatial linear model was useful for characterizing broad patterns related to tidal range, sediment type, temperature, and turbidity and site-specific patterns capturing resilience, all while accounting for spatial dependence, there are some limitations. First, the eReefs water quality covariates are based on modeled data that are not directly observed in the field, potentially propagating modeling error. Second, there are many correlated water quality variables, making it hard (for any model) to isolate the effects of each component of water quality pre- and post-flood. Third, the lack of baseline mean salinity in the 90 days prior to sampling makes it challenging to directly study the effect of mean salinity during the flood period. Fourth, there is a lack of seagrass cover data immediately after the flood -- the first post-flood sampling effort occurred four months after the flood. Fifth, there is a lack of additional data sources that could help quantify the impact of local-scale drivers of total cover pre- and post-flood, potentially varying across both space and time (e.g., light attenuation, biological interactions, etc.). Sixth, the locations at which data were collected were not randomly sampled, which helps to guard against preferential sampling [@diggle2010preferential] and provides more robust and representative data sources [@dumelle2022comparison].

# Revisiting the Research Questions

**Research Question 1**: Can we identify water quality covariates that act as triggers for seagrass cover management responses? Are these triggers different based on whether they result from a specific event (e.g., flood) or are ambient?

We identified ambient 90-day turbidity, turbidity during the flood event, and ambient 90-day water temperature as significant water quality triggers for seagrass cover responses. Our results suggest that ambient and flood-related turbidity act through distinct mechanisms and at different magnitudes; a one-unit increase in 90-day mean turbidity was associated with a 7.18 percentage point decrease in total cover, while the effect of flood turbidity decayed with time since the flood. Across tidal ranges, intertidal sites supported the highest average total cover (26.89%), followed by shallow subtidal (16.51%) and deep subtidal (4.88%), suggesting that management thresholds may need to be calibrated differently depending on depth zone, as deeper sites with inherently lower cover may be disproportionately sensitive to turbidity increases. Similarly, sediment type affected baseline cover levels, as sites with muddy sand and sandy mud substrates supported higher average cover than sand or rock sites. Together, these results suggest that management triggers should differ depending on whether the driver is chronic ambient turbidity or an acute flood event, and should further account for the tidal and sediment context of the sites in question (e.g., perhaps turbidity-driven management responses are most important in areas of densely populated seagrass).


**Research Question 2**: Did seagrass cover improve after the flood event(s)? If so, how and where?

Total seagrass cover improved bay-wide between 2022–2023 and 2023–2024, suggesting the system has been recovering since the flood. However, recovery was spatially heterogeneous — Deception Bay and the Marine National Park Zone returned to pre-flood cover levels, while Southern Bay and the General Use Zone lagged considerably behind. Recovery was also relatively slow compared to the magnitude of the impact in the intertidal tidal range compared to shallow and deep subtidal ranges, suggesting that intertidal sites may warrant prioritized management attention immediately following floods.

**Research Question 3**: Can we identify site-specific resilience indicators, which describe where recovery is likely to be enhanced or impeded?  What factors contribute to seagrass resilience?

We identified several site-specific indicators of seagrass resilience. Dominant species, species richness, and bathymetry mean depth were the strongest predictors of resilience, with sites dominated by CS or ZMHU species, higher species richness, and greater water depth tending to outperform model expectations and hence, infer greater resilience. Gastropod habitat area, bathymetric curvature, and slope provided moderate additional explanatory power. These findings suggest that management actions aimed at protecting structurally diverse, species-rich meadows at deeper sites are likely to be most effective at improving resilience, while sites with low species richness dominated by other species (or no seagrass) may require more active intervention to support recovery.

# Discussion

In this report we apply a rigorous, spatially explicit statistical framework to characterize the drivers, recovery, and resilience of seagrass cover in Moreton Bay across a three-year monitoring period that included the 2022 Moreton Bay flood event. The analytical approach taken here is well-suited to the structure and scale of the data, and the findings have direct implications for adaptive management of the bay's seagrass ecosystem.

The use of spatial linear models as the primary analytical framework is well-justified from both statistical and ecological perspectives. Total percent seagrass cover observations collected across the bay are inherently spatially structured, as nearby sites share similar environmental conditions, water quality and disturbance histories, as well as ecological communities. Models that ignore this dependence may well provide misleading results and inference. By explicitly modeling spatial dependence via a spatial covariance function, the analyses presented here produce more reliable estimates of covariate effects and more honest characterizations of uncertainty than a conventional non-spatial model would provide. These models are typically also more interpretable and characterize uncertainty more effectively than machine learning models. The inclusion of EHMP subzone as a random effect further accounts for structured variation among management-relevant regions of the bay, ensuring that covariate estimates reflect within-subzone patterns rather than being confounded by between-subzone differences in baseline cover. The block kriging helps understand bay-wide trends and can immediately characterize bay-wide predictions on the original total seagrass cover scale (0-100), rather than a transformed scale; providing a flexible framework for scenario modeling under future flood conditions and climate scenarios.

The use of a two-stage modeling framework -- first modeling total cover as a function of water quality and geographic drivers, then modeling the residuals from that model as a function of site-level resilience variables -- is a structured approach to separating bay-wide environmental patterns from site-specific ecological resilience capacity. This builds directly on the framework proposed by @gilby2023drivers to help understand drivers of average resilience and provide useful insight needed to inform broader conservation strategies. The robustness of the findings is further supported by the consistency of results across multiple modeling approaches (e.g., spatial beta regression model, random forest, etc.), which provides evidence that the relationships identified are genuine rather than resulting solely from specific modeling response frameworks. 

Overall, the results described in this report provide a quantitative, spatially structured foundation for adaptive decision making. The understanding of important water quality and geographic variables, average predictions via block Kriging, and the study of site-level resilience indicators provides evidence that can inform  temporally targeted, context-specific interventions throughout the bay.

\clearpage

# References {.unnumbered}

::: {#refs}
:::

# Acknowledgements
The eReefs datasets (model simulations, satellite or in-situ observations) and software were produced as part of the eReefs project (ereefs.org.au), which is a collaboration between Australia’s national science agency CSIRO, the Australian Institute of Marine Science (AIMS) and the Queensland Government, with observations obtained through the Integrated Marine Observing System (IMOS) and the Great Barrier Reef Marine Park Authority (GBRMPA). eReefs is funded by the Australian Government’s Reef Trust.

#Appendix 1. Spatial Linear Models and Random Forest Models {.unnumbered}
The spatially explicit statistical models (i.e., spatial models) we used to answer these research questions are related to commonly used linear and generalized linear models but incorporate _spatial dependence_. Spatial dependence measures how similar two sites are based on their distance. Formally, we can build spatial dependence directly into spatial models by leveraging Tobler's First Law of Geography [@tobler1970computer], which states succinctly that "nearby sites tend to be more similar than distant sites". Not only is this intuitive from a "first principles" perspective, but spatial models also provide many advantages for modeling spatial data compared to models that assume independence among observations (i.e., nonspatial models). First, nonspatial models tend to underestimate uncertainty, yielding covariate standard errors and p-values for covariates that are too small and corresponding confidence intervals that are too narrow. Second, spatial models improve prediction at unsampled locations by borrowing strength from nearby observations, yielding more precise predictions and corresponding standard errors that vary by location. Third, spatial models decompose model error into separate spatial and nonspatial components, informing residual and resilience structure in ways that nonspatial models cannot. Fourth, spatial models can provide ecologically meaningful insight into the spatial dependence structure of the process, informing relationships among measured and unmeasured variables as well as future sampling plans. @zimmerman2024spatial provide a thorough review of spatial models for ecological and environmental data, providing more detail and justification for the claims made above.

Spatial models also offer many benefits over commonly used machine learning algorithms like random forests [@breiman2001random; @cutler2007random; @james2017introduction; @kuhn2022tidy]. First, traditional machine learning algorithms do not incorporate spatial dependence and, by consequence, often perform worse than spatial models at characterizing important covariates and predicting at unobserved locations [@fox2020comparing]. Second, machine learning algorithms struggle to quantify uncertainty in a principled manner, something spatial models excel at. Third, covariate inference is more straightforward using spatial models, as they naturally provide p-values that enable structured testing of hypotheses. Fourth, spatial models provide a rich set of diagnostics like leverage and Cook's distances [@montgomery2021introduction], quantities not generally defined for machine learning algorithms. This is not to say that machine learning algorithms cannot be very useful in ecological data analyses, but rather that spatial models are sometimes a more effective tool for modeling spatial data. 

Though there are many benefits of spatial models over machine learning algorithms, the reverse is also true. It is much simpler to capture complex nonlinearities using "out-of-the-box" machine learning algorithms, bypassing the tuning often required for spatial models that involve nonlinear components like splines. This can be especially useful for ecological data, which often exhibit nonlinear relationships [@fox2017assessing]. Machine learning algorithms also typically rely on few, if any, distributional assumptions and can be more robust against unusual observations like outliers [@hastie2009elements]. Machine learning algorithms also tend to be much more computational efficient than spatial models. Importantly, machine learning algorithms and spatial models can be complementary tools that help shed insight into complicated ecological processes (something we reiterate in the upcoming analyses). An exciting and developing area of statistical research involves incorporating spatial modeling strategies directly into machine learning algorithms [@saha2023random; @heaton2025scalable].

The first spatially explicit model we define is the spatial linear model. The spatial linear model is a _mixed-effects model_ [@pinheiro2000mixed; @bolker2015linear] written as

$$
y_i = \beta_0 + \beta_1x_{1, i} + \beta_2x_{2, i} + ... + \beta_p x_{p, i} + \tau_i + \epsilon_i,
$$ {#eq-splm1}
where $i = 1, 2, ..., n$ indexes the $n$ observations. The quantity $y_i$ is the response variable (for observation $i$), the $\beta$ parameters are slope parameters controlling the covariates, $x_{j, i}$ (for $j = 1, 2, ..., p$). The random error $\tau_i$ is a spatial random error, and the random error $\epsilon_i$ is the nonspatial (i.e., independent) random error. The spatial random error is spatially structured, capturing spatial residual effects like environmental gradients we do not have data to fully represent. The nonspatial random error is not spatially structured and captures nonspatial residual effects like measurement error and microscale variation. The difference between a traditional (i.e., nonspatial) linear model and a spatial linear model is the inclusion of $\tau_i$ (in the spatial linear model). The spatial random error implied by $\tau_i$ is governed by a _spatial covariance function_. The spatial covariance function allows the model to incorporate spatial dependence, yielding the improved model performance previously described. @fig-spcov_fns shows three different spatial covariance functions that describe different ways the spatial covariance decays with distance. Model parameters and standard errors are estimated from data using likelihood-based [@patterson1971recovery; @harville1977maximum; @wolfinger1994computing] or semivariogram-based [@cressie1985fitting; @curriero1999composite] estimation methods. Often, the spatial linear model in @eq-splm1 is written more compactly using _matrix notation_:
$$
\mathbf{y} = \mathbf{X} \boldsymbol{\beta} + \boldsymbol{\tau} + \boldsymbol{\epsilon},
$$ {#eq-splm2}
where $\mathbf{X}$ is a matrix that holds the $j$ covariates for each of $n$ observations, and $\mathbf{y}$, $\boldsymbol{\beta}$, $\boldsymbol{\tau}$, and $\boldsymbol{\epsilon}$ are vectors that contain each of the $n$ or $j$ quantities from @eq-splm1. The spatial linear model can be adjusted to account for nonspatial random effects just as one would add random effects to a nonspatial linear mixed model:
$$
\mathbf{y} = \mathbf{X} \boldsymbol{\beta} + \mathbf{Z}\mathbf{u} + \boldsymbol{\tau} + \boldsymbol{\epsilon},
$$ {#eq-splm3}
where $\mathbf{Z}$ is a random effect design matrix (dimension $n \times k$) that indexes each of the $k$ random effects in $\mathbf{u}$. For more information about specific forms of various spatial covariance functions, see @schabenberger2017statistical and @zimmerman2024spatial.


::: {.cell}
::: {.cell-output-display}
![Three spatial covariance functions that represent different spatial dependence behavior.](MBreport_ep_files/figure-pdf/fig-spcov_fns-1.pdf){#fig-spcov_fns width=75%}
:::
:::


Like the traditional linear model, the spatial linear model can capture complex nonlinearities in the data through the use of interaction effects [@faraway2002practical; @lenth2025emmeans], splines [@chambers1992statistical; @ruppert2003semiparametric], additive models [@wood2017generalized; @pedersen2019hierarchical], and penalized regression models like the ridge and lasso [@tibshirani1996regression]. Generally, this is accomplished by leveraging _basis functions_ that break down nonlinear behavior into a set of additive linear components that may include penalty terms which guard against overfitting. Furthermore, there is no assumption that $y_i$ itself is normally distributed; rather, distributional assumptions are placed on the random errors after incorporating the slope parameters. This clarification implies the spatial (and nonspatial) linear models can be quite effective at modeling a range of data that are not normally distributed. For highly nonnormal data, tools like generalized linear models, which we discuss later, can be helpful. Finally, under general conditions, the slope parameters of the spatial (and nonspatial) linear models are approximately normally distributed with a large sample size (see @casella2024statistical for a general justification via the Central Limit Theorem and @zhang2005towards for nuances regarding spatial data). This is beneficial because it justifies valid hypothesis testing on the slope coefficients (often of primary ecological interest), no matter the distribution of the data (given certain assumptions are satisfied).

One of the most useful applications of the spatial linear model is predicting the response variable at unobserved locations. A very useful spatial prediction method called Kriging has a rich theoretical justification [@cressie1990origins; @stein1999interpolation]. Kriging is the spatial analogue to best linear unbiased prediction (BLUP) for mixed models [@henderson1975best]. Via Kriging, the spatial linear model produces BLUPs at each unobserved location of interest, alongside prediction intervals to characterize uncertainty. While useful, sometimes the goal is to create a BLUP for an entire region of interest. This technique, which builds upon the theoretical foundations of Kriging, is called block Kriging [@ver2002sampling; @cressie2006block]. The intuition is to make point predictions at a very dense grid in the study region and then average these predictions, appropriately characterizing uncertainty. Together, Kriging and block Kriging provide a complete set of tools for spatial prediction, given the fit of a spatial linear model.

Spatial generalized linear models are an extension of spatial linear models to highly nonnormal data [@christensen2002bayesian; @zhang2002estimation; @diggle2007model; @bonat2016practical]. @verhoef2024marginal proposed a novel application of the Laplace approximation to enable restricted maximum likelihood estimation of spatial generalized linear models. Their model formulation builds upon @eq-splm1 and is given by
$$
g(\mu_i) = \beta_0 + \beta_1x_{1, i} + \beta_2x_{2, i} + ... + \beta_p x_{p, i} + \tau_i + \epsilon_i,
$$ {#eq-spglm1}
where $g(\mu_i)$ is a "link function" that links the distribution of $y_i$ to a linear function of its mean, $\mu_i$. Common link functions include the log link for count and skewed data and the logit link for binomial and proportion data (@tbl-glm-families). Using the same argument as for spatial linear models, nonspatial random effects can be added to spatial generalized linear models.


::: {#tbl-glm-families .cell tbl-cap='Common generalized linear model families, link functions, and appropriate data types.'}
::: {.cell-output-display}


|Family            |Link Function                    |Link Name |Data Type            |
|:-----------------|:--------------------------------|:---------|:--------------------|
|Binomial          |$f(\mu) = \log(\mu / (1 - \mu))$ |Logit     |Binary; Binary Count |
|Poisson           |$f(\mu) = \log(\mu)$             |Log       |Count                |
|Negative Binomial |$f(\mu) = \log(\mu)$             |Log       |Count                |
|Beta              |$f(\mu) = \log(\mu / (1 - \mu))$ |Logit     |Proportion           |
|Gamma             |$f(\mu) = \log(\mu)$             |Log       |Skewed               |
|Inverse Gaussian  |$f(\mu) = \log(\mu)$             |Log       |Skewed               |


:::
:::


Spatial models are complicated, both theoretically and computationally. This has contributed to their relative lack of use in the ecological community. However, recent advances in open-source software have made spatial models more accessible to ecologists who have strong subject matter expertise but may lack sufficient statistical background or computational skills to implement them from scratch. One such advancement is the `spmodel` **R** package [@dumelle2023spmodel], which emulates base-**R** functions like `lm()` and `glm()` to account for spatial dependence via `splm()` and `spglm()`, respectively. Familiar **R** functions like `summary()` and `predict()` can be used directly on models fit via `spmodel`. Importantly, `spmodel` also has tools for spatial random forest residual modeling [@fox2020comparing], applications to large data of several thousand observations via spatial indexing [@ver2023indexing], and cross validation [@stone1974cross; @roberts2017cross].

#Appendix 2 {.unnumbered}

This is an appendix

::: {.cell}

:::



::: {.cell}

:::



::: {.cell}

:::




::: {.cell}

:::

