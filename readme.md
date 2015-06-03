Readme
========================================================

Date: czerwiec 03 2015

### List of files

<!-- html table generated in R 3.2.0 by xtable 1.7-4 package -->
<!-- Wed Jun  3 15:32:26 2015 -->
<table border=1>
<tr> <th> File name </th> <th> Description </th> <th> Sourced files </th>  </tr>
  <tr> <td> cs_analysis.R </td> <td> Analysis of cleavage sites amino acid composition. </td> <td> start.R </td> </tr>
  <tr> <td> cs_exploratory_analysis.R </td> <td> Easy plots of amino acids versus their group in the cleavage sites. Use this code to create better groups for cleavage site prediction. </td> <td> start.R </td> </tr>
  <tr> <td> cs_ngrams.R </td> <td> Important n-gram in cleavage sites. </td> <td> start.R </td> </tr>
  <tr> <td> cs_ngrams.RData </td> <td> Data created by cs_ngrams.R </td> <td>  </td> </tr>
  <tr> <td> kmer_cs_prediction.R </td> <td> Adding n-gram to hsmm to predict cleavage site more accurately. </td> <td> start.R </td> </tr>
  <tr> <td> os_AUC50.csv </td> <td> AUC for different species (only species with 50 and more sequences </td> <td>  </td> </tr>
  <tr> <td> os_AUC.csv </td> <td> AUC for different species. Created by taxonomy_prediction.R </td> <td>  </td> </tr>
  <tr> <td> plasmodium_k_mers_signalHsmm.html </td> <td> Parsed plasmodium_k_mers_signalHsmm.Rmd </td> <td>  </td> </tr>
  <tr> <td> plasmodium_k_mers_signalHsmm.Rmd </td> <td> Report: incorporating n-grams to obtain more precise predictions of the cleavage site for Plasmodiidae. </td> <td>  </td> </tr>
  <tr> <td> readme.md </td> <td> Parsed readme.Rmd </td> <td>  </td> </tr>
  <tr> <td> readme.Rmd </td> <td> A template for automatically generated readme. Date and list of file included. </td> <td>  </td> </tr>
  <tr> <td> signalHsmm_with_k_mers.html </td> <td> Parsed signalHsmm_with_k_mers.Rmd </td> <td>  </td> </tr>
  <tr> <td> signalHsmm_with_k_mers.Rmd </td> <td> Report: incorporating n-grams to obtain more precise predictions of the cleavage site. </td> <td>  </td> </tr>
  <tr> <td> signalp_benchmark_fun.R </td> <td> Set of functions to benchmark signalHsmm versus signalP. </td> <td>  </td> </tr>
  <tr> <td> signalp_benchmark.R </td> <td> Benchmark of signalHsmm and signalP on all SPs from UniProt (Data downloaded in June 2015). </td> <td> signalp_benchmark_fun.R </td> </tr>
  <tr> <td> start.R </td> <td> All packages, variables and so on universally needed for the research. Should be sourced every time. </td> <td>  </td> </tr>
  <tr> <td> taxonomy_plots.R </td> <td> Plot which species are recognised by signalHsmm. </td> <td>  </td> </tr>
  <tr> <td> taxonymy_prediction.R </td> <td> signalP and signalHsmm predictions for different species </td> <td>  </td> </tr>
  <tr> <td> whole_taxonomy_prediction.R </td> <td> signalP and signalHsmm predictions for all eukaryotic taxonomic groups (not species). </td> <td>  </td> </tr>
   </table>

### Repository rules
1. Each file must have a desctiption in readme using the R code in *readme.Rmd*.
2. Use comments frequently.
