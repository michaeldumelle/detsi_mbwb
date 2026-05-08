
# MB Questions
* I'm a little confused about the modelling workflow. Did you do any model selection after you removed collinear covariates? Or did you simply keep all the variables in the model? I've done a lot of rearranging and editing and I'm not sure that the statistical methods workflow captures what you actually did. I did my best...

* You say that: "Random forest partial dependence plots [@greenwell2017pdp] suggested that linear effects of the model variables was reasonable." But then the partial dependence plot between total cover and 90 day mean turbidity does not look linear?

* Were the numeric covariates centered and scaled? If so, let's say that. Will be important if we include the fixed effects table. 

* Appendix 2. there are some figures here that aren't being printed or referenced. Mentioning it in case they should be!


# MB Sections for Mike to edit

* Exploratory analysis: Please proof read and edit this whole section. 

* Statistical modeling: Describe model selection - how did you narrow it down to the final set of covariates? I've put some place holder text in this section. 

* I think it's important to show the fixed effects table since one of the main goals is to interpret the relationships between the covariates and response. I suggest either 1) replacing Table 4 (tbl-mb_tc_splm_vars) with a table showing the fixed effects. The info in the current Table 4 is a repeat of what's described in the previous sections. OR - I would remove table 4 and put the FE table in Appendix 2 and then just refer readers there when you discuss the results. I would still keep the ANOVA table in the main document. If you include the FE table for the seagrass model in Appendix 2, I think we should also add the FE table for the resilience model. 

* Re-read the RF model results dot points. Distance to gastropod habitat is listed in dot points 3 and 5.

# Minor comments

* Add affiliations and date to title page

* Do you know the spatial scale over which the area was calculated for the resilience covariates? If not, we should ask EF on Wednesday and add this to the definitions in the Resilience modelling section. Also what is habitat count? And is distance to, the distance to the nearest coral/mangrove/etc. habitat patch? 

# MB Things to check in final pdf

* Read the comments I sent previously in the word document in case any are still relevant.

* Check all the figures and tables to ensure they are fully visible and are not split over pages. 

* I added two references to references.bib for the two technical reports (MB and WB). I haven't used them yet and they're a bit incomplete I think, but they are there if we need them. 

* Figures 5-8 are getting rendered in the Resilience model section instead of the end of the Block kriging section. Can we force them to print before this section heading?