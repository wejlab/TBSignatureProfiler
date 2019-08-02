#' TB Signatures
#'
#' A set of 34 Tuberculosis gene signatures from various publications. This set
#' of signatures uses gene symbols. Attempts have been made to use updated gene
#' symbols and remove symbols that did not match the most recent annotation.
#' Additional sets for Entrez IDs and Ensembl IDs are forthcoming. The assigned
#' signature name ends with the number of genes in the original signature.
#' Note that in some cases signatures will be positive identifiers of TB
#' whereas others are negative identifiers; this should be taken into account
#' when creating ROC curves and computing any AUC estimates.
#'
#' @name TBsignatures
#' @docType data
#' @format list
#' @source
#' \itemize{
#'  \item{\strong{Anderson_42}}{: Anderson, Suzanne T., Myrsini Kaforou, Andrew J. Brent, Victoria J. Wright, Claire M. Banwell, George Chagaluka, Amelia C. Crampin, et al. 2014. "Diagnosis of Childhood Tuberculosis and Host RNA Expression in Africa." The New England Journal of Medicine 370 (18): 1712–23. \href{http://dx.doi.org/10.1056/NEJMoa1303657}{10.1056/NEJMoa1303657}}
#'  \item{\strong{Anderson_OD_51}}{: Anderson, Suzanne T., Myrsini Kaforou, Andrew J. Brent, Victoria J. Wright, Claire M. Banwell, George Chagaluka, Amelia C. Crampin, et al. 2014. "Diagnosis of Childhood Tuberculosis and Host RNA Expression in Africa." The New England Journal of Medicine 370 (18): 1712–23. \href{http://dx.doi.org/10.1056/NEJMoa1303657}{10.1056/NEJMoa1303657}}
#'  \item{\strong{Berry_393}}{: Berry, Matthew P. R., Christine M. Graham, Finlay W. McNab, Zhaohui Xu, Susannah A. A. Bloch, Tolu Oni, Katalin A. Wilkinson, et al. 2010. "An Interferon-Inducible Neutrophil-Driven Blood Transcriptional Signature in Human Tuberculosis." Nature 466 (7309): 973–77. \href{http://dx.doi.org/10.1038/nature09247}{10.1038/nature09247}}
#'  \item{\strong{Berry_OD_86}}{: Berry, Matthew P. R., Christine M. Graham, Finlay W. McNab, Zhaohui Xu, Susannah A. A. Bloch, Tolu Oni, Katalin A. Wilkinson, et al. 2010. "An Interferon-Inducible Neutrophil-Driven Blood Transcriptional Signature in Human Tuberculosis." Nature 466 (7309): 973–77. \href{http://dx.doi.org/10.1038/nature09247}{10.1038/nature09247}}
#'  \item{\strong{Blankley_5}}{: Blankley, Simon, Christine M. Graham, Joe Levin, Jacob Turner, Matthew P. R. Berry, Chloe I. Bloom, Zhaohui Xu, et al. 2016. "A 380-Gene Meta-Signature of Active Tuberculosis Compared with Healthy Controls." The European Respiratory Journal: Official Journal of the European Society for Clinical Respiratory Physiology 47 (6): 1873–76. \href{http://dx.doi.org/10.1183/13993003.02121-2015}{10.1183/13993003.02121-2015}}
#'  \item{\strong{Blankley_380}}{: Blankley, Simon, Christine M. Graham, Joe Levin, Jacob Turner, Matthew P. R. Berry, Chloe I. Bloom, Zhaohui Xu, et al. 2016. "A 380-Gene Meta-Signature of Active Tuberculosis Compared with Healthy Controls." The European Respiratory Journal: Official Journal of the European Society for Clinical Respiratory Physiology 47 (6): 1873–76. \href{http://dx.doi.org/10.1183/13993003.02121-2015}{10.1183/13993003.02121-2015}}
#'  \item{\strong{Bloom_OD_140}}{: Bloom, Chloe I., Christine M. Graham, Matthew P. R. Berry, Fotini Rozakeas, Paul S. Redford, Yuanyuan Wang, Zhaohui Xu, et al. 2013. "Transcriptional Blood Signatures Distinguish Pulmonary Tuberculosis, Pulmonary Sarcoidosis, Pneumonias and Lung Cancers." PloS One 8 (8): e70630. \href{http://dx.doi.org/10.1371/journal.pone.0070630}{10.1371/journal.pone.0070630}}
#'  \item{\strong{Esmail_82}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964–73. \href{http://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Esmail_203}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964–73. \href{http://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Esmail_893}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964–73. \href{http://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Gliddon_OD_3}}{: Gliddon, Harriet D., Kaforou, Myrsini, Alikian, Mary, Habgood-Coote, Dominic, Zhou, Chenxi, Oni, Tolu, Anderson, Suzanne T., Brent, Andrew J., Crampin, Amelia C., Eley, Brian, Kern, Florian, Langford, Paul R., Ottenhoff, Tom H. M., Hibberd, Martin L., French, Neil, Wright, Victoria J., Dockrell, Hazel M., Coin, Lachlan J., Wilkinson, Robert J., Levin, Michael. 2019 "Identification of reduced host transcriptomic signatures for tuberculosis and digital PCR-based validation and quantification" biorxiv.org: . \href{http://dx.doi.org/10.1101/583674}{10.1101/583674}}
#'  \item{\strong{Gliddon_OD_4}}{: Gliddon, Harriet D., Kaforou, Myrsini, Alikian, Mary, Habgood-Coote, Dominic, Zhou, Chenxi, Oni, Tolu, Anderson, Suzanne T., Brent, Andrew J., Crampin, Amelia C., Eley, Brian, Kern, Florian, Langford, Paul R., Ottenhoff, Tom H. M., Hibberd, Martin L., French, Neil, Wright, Victoria J., Dockrell, Hazel M., Coin, Lachlan J., Wilkinson, Robert J., Levin, Michael. 2019 "Identification of reduced host transcriptomic signatures for tuberculosis and digital PCR-based validation and quantification" biorxiv.org: . \href{http://dx.doi.org/10.1101/583674}{10.1101/583674}}
#'  \item{\strong{Jacobsen_3}}{: Jacobsen, Marc, Dirk Repsilber, Andrea Gutschmidt, Albert Neher, Knut Feldmann, Hans J. Mollenkopf, Andreas Ziegler, and Stefan H. E. Kaufmann. 2007. "Candidate Biomarkers for Discrimination between Infection and Disease Caused by Mycobacterium Tuberculosis." Journal of Molecular Medicine  85 (6): 613–21. \href{http://dx.doi.org/10.1007/s00109-007-0157-6}{10.1007/s00109-007-0157-6}}
#'  \item{\strong{Kaforou_27}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{http://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kaforou_OD_44}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{http://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kaforou_OD_53}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{http://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Lee_4}}{: Lee, Shih-Wei, Lawrence Shih-Hsin Wu, Guan-Mau Huang, Kai-Yao Huang, Tzong-Yi Lee, and Julia Tzu-Ya Weng. 2016. "Gene Expression Profiling Identifies Candidate Biomarkers for Active and Latent Tuberculosis." BMC Bioinformatics 17 Suppl 1 (January): 3. \href{http://dx.doi.org/10.1186/s12859-015-0848-x}{10.1186/s12859-015-0848-x}}
#'  \item{\strong{Maertzdorf_4}}{: Maertzdorf, Jeroen, Gayle McEwen, January Weiner 3rd, Song Tian, Eric Lader, Ulrich Schriek, Harriet Mayanja-Kizza, Martin Ota, John Kenneth, and Stefan He Kaufmann. 2016. "Concise Gene Signature for Point-of-Care Classification of Tuberculosis." EMBO Molecular Medicine 8 (2): 86–95. \href{http://dx.doi.org/10.15252/emmm.201505790}{10.15252/emmm.201505790}}
#'  \item{\strong{Maertzdorf_OD_100}}{: Maertzdorf, Jeroen, January Weiner 3rd, Hans-Joachim Mollenkopf, TBornot TB Network, Torsten Bauer, Antje Prasse, Joachim Müller-Quernheim, and Stefan H. E. Kaufmann. 2012. "Common Patterns and Disease-Related Signatures in Tuberculosis and Sarcoidosis." Proceedings of the National Academy of Sciences of the United States of America 109 (20): 7853–58. \href{http://dx.doi.org/10.1073/pnas.1121072109}{10.1073/pnas.1121072109}}
#'  \item{\strong{Rajan_HIV_5}}{: Rajan, Jayant V., Semitala, Fred C., Kamya, Moses R., Yoon, Christina., Mehta, Tejas., Cattamanchi, Adithya., Seielstad, Mark., Montalvo, Lani., Andama, Alfred., Katende, Jane., Asege, Lucy., Nakaye, Martha., Mwebe, Sandra. 2018 "A Novel, 5-Transcript, Whole-blood Gene-expression Signature for Tuberculosis Screening Among People Living With Human Immunodeficiency Virus" Clinical Infectious Diseases XX(XX): 1-7. \href{https://doi.org/10.1093/cid/ciy835}{10.1093/cid/ciy835}}
#'  \item{\strong{Roe_COR_16}}{: Roe, Jennifer, Venturini, Cristina, Gupta, Rishi K., Gurry, Celine, Chain, Benjamin M., Sun, Yuxin, Southern, Jo, Jackson, Charlotte, Lipman, Marc, C., Miller, Robert F., Martineau, Adrian R., Abubakar, Ibrahim, Noursadeghi, Mahdad. 2019 "T1 Blood transcriptomic stratification of short-term risk in contacts of tuberculosis" XX(XX): . \href{https://doi.org/10.1093/cid/ciz252}{10.1093/cid/ciz252}}
#'  \item{\strong{Roe_OD_4}}{: Roe, Jennifer K., Niclas Thomas, Eliza Gil, Katharine Best, Evdokia Tsaliki, Stephen Morris-Jones, Sian Stafford, et al. 2016. "Blood Transcriptomic Diagnosis of Pulmonary and Extrapulmonary Tuberculosis." JCI Insight 1 (16): e87238. \href{http://dx.doi.org/10.1172/jci.insight.87238}{10.1172/jci.insight.87238}}
#'  \item{\strong{Sambarey_HIV_10}}{: Sambarey, Awanti, Abhinandan Devaprasad, Abhilash Mohan, Asma Ahmed, Soumya Nayak, Soumya Swaminathan, George D'Souza, et al. 2017. "Unbiased Identification of Blood-Based Biomarkers for Pulmonary Tuberculosis by Modeling and Mining Molecular Interaction Networks." EBioMedicine 15 (February): 112–26. \href{http://dx.doi.org/10.1016/j.ebiom.2016.12.009}{10.1016/j.ebiom.2016.12.009}}
#'  \item{\strong{Singhania_OD_20}}{: Singhania, Akul, Raman Verma, Christine M. Graham, Jo Lee, Trang Tran, Matthew Richardson, Patrick Lecine, et al. 2018. "A Modular Transcriptional Signature Identifies Phenotypic Heterogeneity of Human Tuberculosis Infection." Nature Communications 9 (1): 2308. \href{http://dx.doi.org/10.1038/s41467-018-04579-w}{10.1038/s41467-018-04579-w}}
#'  \item{\strong{Sloot_HIV_2}}{: Sloot, Rosa, Maarten F. Schim van der Loeff, Erik W. van Zwet, Mariëlle C. Haks, Sytze T. Keizer, Maarten Scholing, Tom H. M. Ottenhoff, Martien W. Borgdorff, and Simone A. Joosten. 2015. "Biomarkers Can Identify Pulmonary Tuberculosis in HIV-Infected Drug Users Months Prior to Clinical Diagnosis." EBioMedicine 2 (2): 172–79. \href{http://dx.doi.org/10.1016/j.ebiom.2014.12.001}{10.1016/j.ebiom.2014.12.001}}
#'  \item{\strong{Suliman_RISK_4}}{: Suliman, Sara, Ethan Thompson, Jayne Sutherland, January Weiner Rd, Martin O. C. Ota, Smitha Shankar, Adam Penn-Nicholson, et al. 2018. "Four-Gene Pan-African Blood Signature Predicts Progression to Tuberculosis." American Journal of Respiratory and Critical Care Medicine, April. https://doi.org/10.1164/rccm.201711-2340OC. \href{http://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{Sweeney_OD_3}}{: Sweeney, Timothy E., Lindsay Braviak, Cristina M. Tato, and Purvesh Khatri. 2016. "Genome-Wide Expression for Diagnosis of Pulmonary Tuberculosis: A Multicohort Analysis." The Lancet. Respiratory Medicine 4 (3): 213–24. \href{http://dx.doi.org/10.1016/S2213-2600(16)00048-5}{10.1016/S2213-2600(16)00048-5}}
#'  \item{\strong{Thompson_9}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48–58. \href{http://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Thompson_FAIL_13}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48–58. \href{http://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Thompson_RES_5}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48–58. \href{http://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Walter_51}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274–82. \href{http://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Walter_PNA_47}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274–82. \href{http://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Walter_PNA_119}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274–82. \href{http://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Zak_RISK_16}}{: Zak, Daniel E., Adam Penn-Nicholson, Thomas J. Scriba, Ethan Thompson, Sara Suliman, Lynn M. Amon, Hassan Mahomed, et al. 2016. "A Blood RNA Signature for Tuberculosis Disease Risk: A Prospective Cohort Study." The Lancet 387 (10035): 2312–22. \href{http://dx.doi.org/10.1016/S0140-6736(15)01316-1}{10.1016/S0140-6736(15)01316-1}}
#'  }
#'
#' @keywords datasets
#' @examples
#' data("TBsignatures")
"TBsignatures"


#' TB Example dataset with Indian population data
#'
#' An example dataset containing the gene expression and metadata in a
#' SummarizedExperiment object for an Indian population. Active TB contamination
#' of the 44 subjects is denoted for each as a "1"(active) or "0"
#' (latent/not present), and can be accessed via \code{TB_indian$label}. The
#' SummarizedExperiment object contains 2 assays (counts and log(counts)),
#' and the column names give the unique subject identification number along
#' with the subject's gender.
#'
#' This dataset was published as part of a study to assess performance of
#' published TB signatures in a South Indian population (Leong et. al 2018).
#' RNA sequencing was performed on whole blood PAX gene samples collected
#' from 28 TB patients and 16 latent TB infected (LTBI) subjects enrolled
#' as part of an ongoing household contact study.
#'
#' @name TB_indian
#' @docType data
#' @format SummarizedExperiment
#' @keywords datasets
#' @references
#' Leong S., Zhao Y., et. al. (2018). Existing blood transcriptional classifiers
#' accurately discriminate active tuberculosis from latent infection in
#' individuals from south India. \emph{Tuberculosis} \strong{109}, 41-51.
#' doi: \href{https://doi.org/10.1016/j.tube.2018.01.002}{10.1016/j.tube.2018.01.002}.
#'
#' @examples
#' data("TB_indian")
"TB_indian"


#' TB Example dataset with TB/HIV data
#'
#' An example dataset containing the gene expression and metadata in a
#' SummarizedExperiment object for 31 subjects with HIV and/or Tuberculosis
#' diseases. Information on subject infection status can be accessed with
#' \code{TB_hiv$Disease}. Samples with both TB and HIV contamination are
#' marked as \code{tb_hiv}, while samples with HIV and no TB are marked
#' as \code{hiv_only}.
#'
#' This dataset was published as part of a study to assess assess if gene expression
#' signatures and cytokine levels would distinguish active TB in advanced HIV
#' in a cohort residing in Sub-Saharan Africa (Verma et. al 2018).
#' Participants were severely immunosuppressed TB-HIV patients who had
#' not yet received TB treatment or anti-retroviral therapy (ART). The dataset included
#' in this package has been lightly edited from the originally published dataset
#' due to the removal of one participant who was HIV positive, on ART and developed
#' TB during follow-up.
#'
#' @name TB_hiv
#' @docType data
#' @format SummarizedExperiment
#' @keywords datasets
#' @references
#' Verma S., Du P., et. al. (2018). Tuberculosis in advanced HIV infection is
#' associated with increased expression of IFN and its downstream targets.
#' \emph{BMC Infectious Diseases} \strong{18:220}.
#' doi: \href{https://doi.org/10.1186/s12879-018-3127-4}{10.1186/s12879-018-3127-4}.
#'
#' @examples
#' data("TB_hiv")
"TB_hiv"

