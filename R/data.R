globalVariables(c("BS_AUC", "FPR", "LowerTPR", "Signatures",
                  "TBsignatures", "TPR", "UpperTPR", "sigAnnotData"))

#' A list of published TB signatures.
#'
#' A set of Tuberculosis gene signatures from various publications. This set
#' of signatures uses gene symbols. Attempts have been made to use updated gene
#' symbols and remove symbols that did not match the most recent annotation.
#' Additional sets for Entrez IDs and Ensembl IDs are forthcoming.
#'
#' Signature names are composed of the last name of the primary author,
#' followed by
#' a possible context for the signature, and ending with either the number of
#' gene transcripts or genes in the signature, with respect to however
#' it was described in the signature in the original publication.
#'
#' Possible signature contexts:
#' \itemize{
#' \item{OD: Other diseases}
#' \item{HIV: Human Immunodeficiency Virus}
#' \item{PNA: Pneumonia}
#' \item{RISK: Risk of developing active TB}
#' \item{RES: Response to TB treatment}
#' \item{FAIL: Failure of TB treatment}
#' }
#'
#' Note that in some cases signatures will be positive identifiers of TB
#' whereas others are negative identifiers; this should be taken into account
#' when creating ROC curves and computing any AUC estimates.
#'
#' @name TBsignatures
#' @docType data
#' @format list
#' @source
#' \itemize{
#'  \item{\strong{Anderson_42}}{: Anderson, Suzanne T., Myrsini Kaforou, Andrew J. Brent, Victoria J. Wright, Claire M. Banwell, George Chagaluka, Amelia C. Crampin, et al. 2014. "Diagnosis of Childhood Tuberculosis and Host RNA Expression in Africa." The New England Journal of Medicine 370 (18): 1712-23. \href{https://dx.doi.org/10.1056/NEJMoa1303657}{10.1056/NEJMoa1303657}}
#'  \item{\strong{Anderson_OD_51}}{: Anderson, Suzanne T., Myrsini Kaforou, Andrew J. Brent, Victoria J. Wright, Claire M. Banwell, George Chagaluka, Amelia C. Crampin, et al. 2014. "Diagnosis of Childhood Tuberculosis and Host RNA Expression in Africa." The New England Journal of Medicine 370 (18): 1712-23. \href{https://dx.doi.org/10.1056/NEJMoa1303657}{10.1056/NEJMoa1303657}}
#'  \item{\strong{Berry_393}}{: Berry, Matthew P. R., Christine M. Graham, Finlay W. McNab, Zhaohui Xu, Susannah A. A. Bloch, Tolu Oni, Katalin A. Wilkinson, et al. 2010. "An Interferon-Inducible Neutrophil-Driven Blood Transcriptional Signature in Human Tuberculosis." Nature 466 (7309): 973-77. \href{https://dx.doi.org/10.1038/nature09247}{10.1038/nature09247}}
#'  \item{\strong{Berry_OD_86}}{: Berry, Matthew P. R., Christine M. Graham, Finlay W. McNab, Zhaohui Xu, Susannah A. A. Bloch, Tolu Oni, Katalin A. Wilkinson, et al. 2010. "An Interferon-Inducible Neutrophil-Driven Blood Transcriptional Signature in Human Tuberculosis." Nature 466 (7309): 973-77. \href{https://dx.doi.org/10.1038/nature09247}{10.1038/nature09247}}
#'  \item{\strong{Blankley_380}}{: Blankley, Simon, Christine M. Graham, Joe Levin, Jacob Turner, Matthew P. R. Berry, Chloe I. Bloom, Zhaohui Xu, et al. 2016. "A 380-Gene Meta-Signature of Active Tuberculosis Compared with Healthy Controls." The European Respiratory Journal: Official Journal of the European Society for Clinical Respiratory Physiology 47 (6): 1873-76. \href{https://dx.doi.org/10.1183/13993003.02121-2015}{10.1183/13993003.02121-2015}}
#'  \item{\strong{Blankley_5}}{: Blankley, Simon, Christine M. Graham, Joe Levin, Jacob Turner, Matthew P. R. Berry, Chloe I. Bloom, Zhaohui Xu, et al. 2016. "A 380-Gene Meta-Signature of Active Tuberculosis Compared with Healthy Controls." The European Respiratory Journal: Official Journal of the European Society for Clinical Respiratory Physiology 47 (6): 1873-76. \href{https://dx.doi.org/10.1183/13993003.02121-2015}{10.1183/13993003.02121-2015}}
#'  \item{\strong{Bloom_OD_144}}{: Bloom, Chloe I., Christine M. Graham, Matthew P. R. Berry, Fotini Rozakeas, Paul S. Redford, Yuanyuan Wang, Zhaohui Xu, et al. 2013. "Transcriptional Blood Signatures Distinguish Pulmonary Tuberculosis, Pulmonary Sarcoidosis, Pneumonias and Lung Cancers." PloS One 8 (8): e70630. \href{https://dx.doi.org/10.1371/journal.pone.0070630}{10.1371/journal.pone.0070630}}
#'  \item{\strong{Bloom_RES_268}}{: Bloom CI, Graham CM, Berry MP, et al. Detectable changes in the blood transcriptome are present after two weeks of antituberculosis therapy. PLoS One. 2012;7(10):e46191. \href{https://dx.doi.org/10.3389/fmicb.2021.650567}{10.3389/fmicb.2021.650567}}
#'  \item{\strong{Bloom_RES_558}}{: Bloom CI, Graham CM, Berry MP, et al. Detectable changes in the blood transcriptome are present after two weeks of antituberculosis therapy. PLoS One. 2012;7(10):e46191. \href{https://dx.doi.org/10.3389/fmicb.2021.650567}{10.3389/fmicb.2021.650567}}
#'  \item{\strong{Chen_HIV_4}}{: Chen Y, Wang Q, Lin S, et al. Meta-Analysis of Peripheral Blood Transcriptome Datasets Reveals a Biomarker Panel for Tuberculosis in Patients Infected With HIV. Front Cell Infect Microbiol. 2021;11:585919. Published 2021 Mar 19. \href{https://dx.doi.org/10.3389/fcimb.2021.585919}{10.3389/fcimb.2021.585919}}
#'  \item{\strong{Chendi_HIV_2}}{: Chendi BH, Tveiten H, Snyders CI, et al. CCL1 and IL-2Ra differentiate Tuberculosis disease from latent infection Irrespective of HIV infection in low TB burden countries [published online ahead of print, 2021 Jul 29]. J Infect. 2021;S0163-4453(21)00379-0. \href{https://dx.doi.org/10.1016/j.jinf.2021.07.036}{10.1016/j.jinf.2021.07.036}}
#'  \item{\strong{Darboe_RISK_11}}{: Darboe, F. et al. Diagnostic performance of an optimized transcriptomic signature of risk of tuberculosis in cryopreserved peripheral blood mononuclear cells. Tuberculosis 108, 124-126 (2018). \href{https://dx.doi.org/ 10.1016/j.tube.2017.11.001}{ 10.1016/j.tube.2017.11.001}}
#'  \item{\strong{Dawany_HIV_251}}{: Dawany, N. et al. Identification of a 251 gene expression signature that can accurately detect M. tuberculosis in patients with and without HIV co-infection. PLoS One 9, (2014). \href{https://dx.doi.org/10.1371/journal.pone.0089925}{10.1371/journal.pone.0089925}}
#'  \item{\strong{Duffy_23}}{: Duffy FJ, Olson GS, Gold ES, Jahn A, Aderem A, Aitchison J, Rothchild AC, Diercks AH, Nemeth J. A contained Mycobacterium tuberculosis mouse infection model predicts active disease and containment in humans. The Journal of Infectious Diseases. 2021 Mar 10. \href{https://dx.doi.org/10.1093/infdis/jiab130}{10.1093/infdis/jiab130}}
#'  \item{\strong{Esmail_203}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964-73. \href{https://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Esmail_82}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964-73. \href{https://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Esmail_893}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964-73. \href{https://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Estevez_133}}{: Estévez O, Anibarro L, Garet E, et al. An RNA-seq Based Machine Learning Approach Identifies Latent Tuberculosis Patients With an Active Tuberculosis Profile. Front Immunol. 2020;11:1470. Published 2020 Jul 14. \href{https://dx.doi.org/10.3389/fimmu.2020.01470}{10.3389/fimmu.2020.01470}}
#'  \item{\strong{Estevez_259}}{: Estévez O, Anibarro L, Garet E, et al. An RNA-seq Based Machine Learning Approach Identifies Latent Tuberculosis Patients With an Active Tuberculosis Profile. Front Immunol. 2020;11:1470. Published 2020 Jul 14. \href{https://dx.doi.org/10.3389/fimmu.2020.01470}{10.3389/fimmu.2020.01470}}
#'  \item{\strong{Gjoen_10}}{: Gjøen, J.E., Jenum, S., Sivakumaran, D. et al. 'Novel transcriptional signatures for sputum-independent diagnostics of tuberculosis in children.' Sci Rep 7, 5839 (2017). \href{https://doi.org/10.1038/s41598-017-05057-x}{10.1038/s41598-017-05057-x}}
#'  \item{\strong{Gjoen_7}}{: Gjøen, J.E., Jenum, S., Sivakumaran, D. et al. 'Novel transcriptional signatures for sputum-independent diagnostics of tuberculosis in children.' Sci Rep 7, 5839 (2017). \href{https://doi.org/10.1038/s41598-017-05057-x}{10.1038/s41598-017-05057-x}}
#'  \item{\strong{Gliddon_2_OD_4}}{: Gliddon HD, Kaforou M, Alikian M, et al. Identification of Reduced Host Transcriptomic Signatures for Tuberculosis Disease and Digital PCR-Based Validation and Quantification. Front Immunol. 2021;12:637164. Published 2021 Mar 2. \href{https://dx.doi.org/10.3389/fimmu.2021.637164}{10.3389/fimmu.2021.637164}}
#'  \item{\strong{Gliddon_HIV_3}}{: Gliddon HD, Kaforou M, Alikian M, et al. Identification of Reduced Host Transcriptomic Signatures for Tuberculosis Disease and Digital PCR-Based Validation and Quantification. Front Immunol. 2021;12:637164. Published 2021 Mar 2. \href{https://dx.doi.org/10.3389/fimmu.2021.637164}{10.3389/fimmu.2021.637164}}
#'  \item{\strong{Gliddon_OD_3}}{: Gliddon, Harriet D., Kaforou, Myrsini, Alikian, Mary, Habgood-Coote, Dominic, Zhou, Chenxi, Oni, Tolu, Anderson, Suzanne T., Brent, Andrew J., Crampin, Amelia C., Eley, Brian, Kern, Florian, Langford, Paul R., Ottenhoff, Tom H. M., Hibberd, Martin L., French, Neil, Wright, Victoria J., Dockrell, Hazel M., Coin, Lachlan J., Wilkinson, Robert J., Levin, Michael. 2019 "Identification of reduced host transcriptomic signatures for tuberculosis and digital PCR-based validation and quantification" biorxiv.org:\href{https://dx.doi.org/10.1101/583674}{10.1101/583674}}
#'  \item{\strong{Gliddon_OD_4}}{: Gliddon, Harriet D., Kaforou, Myrsini, Alikian, Mary, Habgood-Coote, Dominic, Zhou, Chenxi, Oni, Tolu, Anderson, Suzanne T., Brent, Andrew J., Crampin, Amelia C., Eley, Brian, Kern, Florian, Langford, Paul R., Ottenhoff, Tom H. M., Hibberd, Martin L., French, Neil, Wright, Victoria J., Dockrell, Hazel M., Coin, Lachlan J., Wilkinson, Robert J., Levin, Michael. 2019 "Identification of reduced host transcriptomic signatures for tuberculosis and digital PCR-based validation and quantification" biorxiv.org:\href{https://dx.doi.org/10.1101/583674}{10.1101/583674}}
#'  \item{\strong{Gong_OD_4}}{: Gong Z, Gu Y, Xiong K, Niu J, Zheng R, Su B, Fan L and Xie J (2021) The Evaluation and Validation of Blood-Derived Novel Biomarkers for Precise and Rapid Diagnosis of Tuberculosis in Areas With High-TB Burden. Front. Microbiol. 12:650567. \href{https://dx.doi.org/10.3389/fmicb.2021.650567}{10.3389/fmicb.2021.650567}}
#'  \item{\strong{Heycken_FAIL_22}}{: Heyckendorf J, Marwitz S, Reimann M, et al. Prediction of anti-tuberculosis treatment duration based on a 22-gene transcriptomic model [published online ahead of print, 2021 Feb 11]. Eur Respir J. 2021;2003492. \href{https://dx.doi.org/10.1183/13993003.03492-2020}{10.1183/13993003.03492-2020}}
#'  \item{\strong{Hoang_OD_13}}{: Hoang, Long & Jain, Pooja & Pillay, Timesh & Tolosa-Wright, Mica & Niazi, Umar & Takwoingi, Yemisi & Halliday, Alice & Berrocal-Almanza, Luis & Deeks, Jonathan & Beverley, Peter & Kon, Onn & Lalvani, Ajit. (2021). Transcriptomic signatures for diagnosing tuberculosis in clinical practice: a prospective, multicentre cohort study. The Lancet Infectious Diseases. \href{https://dx.doi.org/10.1016/S1473-3099(20)30928-2}{10.1016/S1473-3099(20)30928-2}}
#'  \item{\strong{Hoang_OD_20}}{: Hoang, Long & Jain, Pooja & Pillay, Timesh & Tolosa-Wright, Mica & Niazi, Umar & Takwoingi, Yemisi & Halliday, Alice & Berrocal-Almanza, Luis & Deeks, Jonathan & Beverley, Peter & Kon, Onn & Lalvani, Ajit. (2021). Transcriptomic signatures for diagnosing tuberculosis in clinical practice: a prospective, multicentre cohort study. The Lancet Infectious Diseases. \href{https://dx.doi.org/10.1016/S1473-3099(20)30928-2}{10.1016/S1473-3099(20)30928-2}}
#'  \item{\strong{Hoang_OD_3}}{: Hoang, Long & Jain, Pooja & Pillay, Timesh & Tolosa-Wright, Mica & Niazi, Umar & Takwoingi, Yemisi & Halliday, Alice & Berrocal-Almanza, Luis & Deeks, Jonathan & Beverley, Peter & Kon, Onn & Lalvani, Ajit. (2021). Transcriptomic signatures for diagnosing tuberculosis in clinical practice: a prospective, multicentre cohort study. The Lancet Infectious Diseases. \href{https://dx.doi.org/10.1016/S1473-3099(20)30928-2}{10.1016/S1473-3099(20)30928-2}}
#'  \item{\strong{Huang_13}}{: Huang, Hai-Hui et al. 'Identification of 13 Blood-based Gene Expression Signatures to Accurately Distinguish Tuberculosis from Other Pulmonary Diseases and Healthy Controls'. 1 Jan. 2015 : S1837 - S1843.\href{https://doi.org/10.3233/BME-151486}{10.3233/BME-151486}}
#'  \item{\strong{Jacobsen_3}}{: Jacobsen, Marc, Dirk Repsilber, Andrea Gutschmidt, Albert Neher, Knut Feldmann, Hans J. Mollenkopf, Andreas Ziegler, and Stefan H. E. Kaufmann. 2007. "Candidate Biomarkers for Discrimination between Infection and Disease Caused by Mycobacterium Tuberculosis." Journal of Molecular Medicine  85 (6): 613-21. \href{https://dx.doi.org/10.1007/s00109-007-0157-6}{10.1007/s00109-007-0157-6}}
#'  \item{\strong{Jenum_8}}{: Jenum, S., Dhanasekaran, S., Lodha, R. et al. Approaching a diagnostic point-of-care test for pediatric tuberculosis through evaluation of immune biomarkers across the clinical disease spectrum. Sci Rep 6, 18520 (2016). \href{https://doi.org/10.1038/srep18520}{10.1038/srep18520}}
#'  \item{\strong{Kaforou_27}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{https://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kaforou_OD_44}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{https://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kaforou_OD_53}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{https://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kulkarni_HIV_2}}{: Kulkarni V, Queiroz ATL, Sangle S, et al. A Two-Gene Signature for Tuberculosis Diagnosis in Persons With Advanced HIV. Front Immunol. 2021;12:631165. Published 2021 Feb 22. \href{https://dx.doi.org/10.3389/fimmu.2021.631165}{10.3389/fimmu.2021.631165}}
#'  \item{\strong{LauxdaCosta_OD_3}}{: Laux da Costa L, Delcroix M, Dalla Costa ER, et al. A real-time PCR signature to discriminate between tuberculosis and other pulmonary diseases. Tuberculosis (Edinb). 2015;95(4):421-425. \href{https://dx.doi.org/10.1016/j.tube.2015.04.008}{10.1016/j.tube.2015.04.008}}
#'  \item{\strong{Lee_4}}{: Lee, Shih-Wei, Lawrence Shih-Hsin Wu, Guan-Mau Huang, Kai-Yao Huang, Tzong-Yi Lee, and Julia Tzu-Ya Weng. 2016. "Gene Expression Profiling Identifies Candidate Biomarkers for Active and Latent Tuberculosis." BMC Bioinformatics 17 Suppl 1 (January): 3. \href{https://dx.doi.org/10.1186/s12859-015-0848-x}{10.1186/s12859-015-0848-x}}
#'  \item{\strong{Leong_24}}{: Leong, Samantha, Yue Zhao, Noyal M. Joseph, Natasha S. Hochberg, Sonali Sarkar, Jane Pleskunas, David Hom, et al. 2018. "Existing blood transcriptional classifiers accurately discriminate active tuberculosis from latent infection in individuals from south India." Tuberculosis (109): 41-51. \href{https://doi.org/10.1016/j.tube.2018.01.002}{10.1016/j.tube.2018.01.002}}
#'  \item{\strong{Leong_RISK_29}}{: Leong, S., Zhao, Y., Ribeiro-Rodrigues, R., Jones-López, E. C., Acuña-Villaorduña, C., Rodrigues, P. M., Palaci, M., Alland, D., Dietze, R., Ellner, J. J., Johnson, W. E., Salgame, P., Cross-validation of existing signatures and derivation of a novel 29-gene transcriptomic signature predictive of progression to TB in a Brazilian cohort of household contacts of pulmonary TB. Tuberculosis (Edinb). 2020 Jan;120:101898. \href{https://dx.doi.org/10.1016/j.tube.2020.101898}{10.1016/j.tube.2020.101898}}
#'  \item{\strong{Maertzdorf_15}}{: Maertzdorf J, McEwen G, Weiner J 3rd, et al. Concise gene signature for point-of-care classification of tuberculosis. EMBO Mol Med. 2016;8(2):86-95. \href{https://dx.doi.org/10.15252/emmm.201505790}{10.15252/emmm.201505790}}
#'  \item{\strong{Maertzdorf_4}}{: Maertzdorf, Jeroen, Gayle McEwen, January Weiner 3rd, Song Tian, Eric Lader, Ulrich Schriek, Harriet Mayanja-Kizza, Martin Ota, John Kenneth, and Stefan He Kaufmann. 2016. "Concise Gene Signature for Point-of-Care Classification of Tuberculosis." EMBO Molecular Medicine 8 (2): 86-95. \href{https://dx.doi.org/10.15252/emmm.201505790}{10.15252/emmm.201505790}}
#'  \item{\strong{Maertzdorf_OD_100}}{: Maertzdorf, Jeroen, January Weiner 3rd, Hans-Joachim Mollenkopf, TBornot TB Network, Torsten Bauer, Antje Prasse, Joachim Müller-Quernheim, and Stefan H. E. Kaufmann. 2012. "Common Patterns and Disease-Related Signatures in Tuberculosis and Sarcoidosis." Proceedings of the National Academy of Sciences of the United States of America 109 (20): 7853-58. \href{https://dx.doi.org/10.1073/pnas.1121072109}{10.1073/pnas.1121072109}}
#'  \item{\strong{PennNich_RISK_6}}{: Penn-Nicholson, A. et al. RISK6, a 6-gene transcriptomic signature of TB disease risk, diagnosis and treatment response. Sci. Rep. 10, (2020). \href{https://dx.doi.org/10.1038/s41598-020-65043-8}{10.1038/s41598-020-65043-8}}
#'  \item{\strong{Qian_OD_17}} {: Qian, Zhongqing et al. "Expression of nuclear factor, erythroid 2-like 2-mediated genes differentiates tuberculosis." Tuberculosis (Edinburgh, Scotland) vol. 99 (2016): 56-62. \href{https://doi.org/10.1016/j.tube.2016.04.008}{10.1016/j.tube.2016.04.008}}
#'  \item{\strong{Rajan_HIV_5}}{: Rajan, Jayant V., Semitala, Fred C., Kamya, Moses R., Yoon, Christina., Mehta, Tejas., Cattamanchi, Adithya., Seielstad, Mark., Montalvo, Lani., Andama, Alfred., Katende, Jane., Asege, Lucy., Nakaye, Martha., Mwebe, Sandra. 2018 "A Novel, 5-Transcript, Whole-blood Gene-expression Signature for Tuberculosis Screening Among People Living With Human Immunodeficiency Virus" Clinical Infectious Diseases: 1-7. \href{https://doi.org/10.1093/cid/ciy835}{10.1093/cid/ciy835}}
#'  \item{\strong{Roe_3}}{: Roe, Jennifer, Venturini, Cristina, Gupta, Rishi K., Gurry, Celine, Chain, Benjamin M., Sun, Yuxin, Southern, Jo, Jackson, Charlotte, Lipman, Marc, C., Miller, Robert F., Martineau, Adrian R., Abubakar, Ibrahim, Noursadeghi, Mahdad. 2019 "T1 Blood transcriptomic stratification of short-term risk in contacts of tuberculosis": . \href{https://doi.org/10.1093/cid/ciz252}{10.1093/cid/ciz252}}
#'  \item{\strong{Roe_OD_4}}{: Roe, Jennifer K., Niclas Thomas, Eliza Gil, Katharine Best, Evdokia Tsaliki, Stephen Morris-Jones, Sian Stafford, et al. 2016. "Blood Transcriptomic Diagnosis of Pulmonary and Extrapulmonary Tuberculosis." JCI Insight 1 (16): e87238. \href{https://dx.doi.org/10.1172/jci.insight.87238}{10.1172/jci.insight.87238}}
#'  \item{\strong{Sambarey_HIV_10}}{: Sambarey, Awanti, Abhinandan Devaprasad, Abhilash Mohan, Asma Ahmed, Soumya Nayak, Soumya Swaminathan, George D'Souza, et al. 2017. "Unbiased Identification of Blood-Based Biomarkers for Pulmonary Tuberculosis by Modeling and Mining Molecular Interaction Networks." EBioMedicine 15 (February): 112-26. \href{https://dx.doi.org/10.1016/j.ebiom.2016.12.009}{10.1016/j.ebiom.2016.12.009}}
#'  \item{\strong{Singhania_OD_20}}{: Singhania, Akul, Raman Verma, Christine M. Graham, Jo Lee, Trang Tran, Matthew Richardson, Patrick Lecine, et al. 2018. "A Modular Transcriptional Signature Identifies Phenotypic Heterogeneity of Human Tuberculosis Infection." Nature Communications 9 (1): 2308. \href{https://dx.doi.org/10.1038/s41467-018-04579-w}{10.1038/s41467-018-04579-w}}
#'  \item{\strong{Sivakumaran_11}}{: Sivakumaran D, Ritz C, Gjøen JE, et al. Host Blood RNA Transcript and Protein Signatures for Sputum-Independent Diagnostics of Tuberculosis in Adults. Front Immunol. 2021;11:626049. Published 2021 Feb 4. \href{https://dx.doi.org/10.3389/fimmu.2020.626049}{10.3389/fimmu.2020.626049}}
#'  \item{\strong{Sloot_HIV_2}}{: Sloot, Rosa, Maarten F. Schim van der Loeff, Erik W. van Zwet, Mariëlle C. Haks, Sytze T. Keizer, Maarten Scholing, Tom H. M. Ottenhoff, Martien W. Borgdorff, and Simone A. Joosten. 2015. "Biomarkers Can Identify Pulmonary Tuberculosis in HIV-Infected Drug Users Months Prior to Clinical Diagnosis." EBioMedicine 2 (2): 172-79. \href{https://dx.doi.org/10.1016/j.ebiom.2014.12.001}{10.1016/j.ebiom.2014.12.001}}
#'  \item{\strong{Suliman_4}}{: Suliman, Sara, Ethan Thompson, Jayne Sutherland, January Weiner Rd, Martin O. C. Ota, Smitha Shankar, Adam Penn-Nicholson, et al. 2018. "Four-Gene Pan-African Blood Signature Predicts Progression to Tuberculosis." American Journal of Respiratory and Critical Care Medicine, April. https://doi.org/10.1164/rccm.201711-2340OC. \href{https://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{Suliman_RISK_2}}{: Suliman, S. et al. Four-gene pan-African blood signature predicts progression to tuberculosis. Am. J. Respir. Crit. Care Med. 197, 1198-1208 (2018). \href{https://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{Suliman_RISK_4}}{: Suliman, Sara, Ethan Thompson, Jayne Sutherland, January Weiner Rd, Martin O. C. Ota, Smitha Shankar, Adam Penn-Nicholson, et al. 2018. "Four-Gene Pan-African Blood Signature Predicts Progression to Tuberculosis." American Journal of Respiratory and Critical Care Medicine, April. https://doi.org/10.1164/rccm.201711-2340OC. \href{https://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{Sweeney_OD_3}}{: Sweeney, Timothy E., Lindsay Braviak, Cristina M. Tato, and Purvesh Khatri. 2016. "Genome-Wide Expression for Diagnosis of Pulmonary Tuberculosis: A Multicohort Analysis." The Lancet. Respiratory Medicine 4 (3): 213-24. \href{https://dx.doi.org/10.1016/S2213-2600(16)00048-5}{10.1016/S2213-2600(16)00048-5}}
#'  \item{\strong{Thompson_9}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48-58. \href{https://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Thompson_FAIL_13}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48-58. \href{https://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Thompson_RES_5}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48-58. \href{https://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Tornheim_71}}{: Tornheim, Jeffrey A., Anil K. Madugundu, Mandar Paradkar, Kiyoshi F. Fukutani, Artur TL Queiroz, Nikhil Gupte, Akshay N. Gupte et al. 2020. "Transcriptomic Profiles of Confirmed Pediatric Tuberculosis Patients and Household Contacts Identifies Active Tuberculosis, Infection, and Treatment Response Among Indian Children." The Journal of Infectious Diseases 221(10): 1647-1658. \href{https://doi.org/10.1093/infdis/jiz639}{10.1093/infdis/jiz639}}
#'  \item{\strong{Tornheim_RES_25}}{: Tornheim, Jeffrey A., Anil K. Madugundu, Mandar Paradkar, Kiyoshi F. Fukutani, Artur TL Queiroz, Nikhil Gupte, Akshay N. Gupte et al. 2020. "Transcriptomic Profiles of Confirmed Pediatric Tuberculosis Patients and Household Contacts Identifies Active Tuberculosis, Infection, and Treatment Response Among Indian Children." The Journal of Infectious Diseases 221(10): 1647-1658. \href{https://doi.org/10.1093/infdis/jiz639}{10.1093/infdis/jiz639}}
#'  \item{\strong{Verhagen_10}} {: Verhagen, L.M., Zomer, A., Maes, M. et al. A predictive signature gene set for discriminating active from latent tuberculosis in Warao Amerindian children. BMC Genomics 14, 74 (2013). \href{https://doi.org/10.1186/1471-2164-14-74}{10.1186/1471-2164-14-74}}
#'  \item{\strong{Walter_51}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274-82. \href{https://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Walter_PNA_119}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274-82. \href{https://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Walter_PNA_47}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274-82. \href{https://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Zak_RISK_16}}{: Zak, Daniel E., Adam Penn-Nicholson, Thomas J. Scriba, Ethan Thompson, Sara Suliman, Lynn M. Amon, Hassan Mahomed, et al. 2016. "A Blood RNA Signature for Tuberculosis Disease Risk: A Prospective Cohort Study." The Lancet 387 (10035): 231222. \href{https://dx.doi.org/10.1016/S0140-6736(15)01316-1}{10.1016/S0140-6736(15)01316-1}}
#'  \item{\strong{Zhao_NANO_6}}{: To be available when published \href{https://dx.doi.org/To be available when published}{To be available when published}}
#'  \item{\strong{Zimmer_RES_3}}{: Zimmer, A.J., Schumacher, S.G., Södersten, E. et al. A novel blood-based assay for treatment monitoring of tuberculosis. BMC Res Notes 14, 247 (2021). \href{https://dx.doi.org/10.1186/s13104-021-05663-z}{10.1186/s13104-021-05663-z}}
#'  }
#'
#' @keywords datasets
#' @examples
#' data("TBsignatures")
"TBsignatures"

#' A list of published TB signatures, using author-given names.
#'
#' A set of Tuberculosis gene signatures from various publications. This set
#' of signatures uses gene symbols. Attempts have been made to use updated gene
#' symbols and remove symbols that did not match the most recent annotation.
#' Additional sets for Entrez IDs and Ensembl IDs are forthcoming.
#'
#' This list differs from \code{TBsignatures} in that signatures with names
#' specified in their originating publication (or that of a peer)
#' are given that common name rather than using the \code{TBSignatureProfiler}
#' naming system. Otherwise, signature names are composed of the last name
#' of the primary author, followed by
#' a possible context for the signature, and ending with either the number of
#' gene transcripts or genes in the signature with respect to however
#' it was described in the original publication.
#'
#' Possible signature contexts:
#' \itemize{
#' \item{OD: Other diseases}
#' \item{HIV: Human Immunodeficiency Virus}
#' \item{PNA: Pneumonia}
#' \item{RISK: Risk of developing active TB}
#' \item{RES: Response to TB treatment}
#' \item{FAIL: Failure of TB treatment}
#' }
#'
#' Note that in some cases signatures will be positive identifiers of TB
#' whereas others are negative identifiers; this should be taken into account
#' when creating ROC curves and computing any AUC estimates.
#'
#' @name TBcommon
#' @docType data
#' @format list
#' @source
#' \itemize{
#'  \item{\strong{Anderson_42}}{: Anderson, Suzanne T., Myrsini Kaforou, Andrew J. Brent, Victoria J. Wright, Claire M. Banwell, George Chagaluka, Amelia C. Crampin, et al. 2014. "Diagnosis of Childhood Tuberculosis and Host RNA Expression in Africa." The New England Journal of Medicine 370 (18): 1712-23. \href{https://dx.doi.org/10.1056/NEJMoa1303657}{10.1056/NEJMoa1303657}}
#'  \item{\strong{Anderson_OD_51}}{: Anderson, Suzanne T., Myrsini Kaforou, Andrew J. Brent, Victoria J. Wright, Claire M. Banwell, George Chagaluka, Amelia C. Crampin, et al. 2014. "Diagnosis of Childhood Tuberculosis and Host RNA Expression in Africa." The New England Journal of Medicine 370 (18): 1712-23. \href{https://dx.doi.org/10.1056/NEJMoa1303657}{10.1056/NEJMoa1303657}}
#'  \item{\strong{Berry_393}}{: Berry, Matthew P. R., Christine M. Graham, Finlay W. McNab, Zhaohui Xu, Susannah A. A. Bloch, Tolu Oni, Katalin A. Wilkinson, et al. 2010. "An Interferon-Inducible Neutrophil-Driven Blood Transcriptional Signature in Human Tuberculosis." Nature 466 (7309): 973-77. \href{https://dx.doi.org/10.1038/nature09247}{10.1038/nature09247}}
#'  \item{\strong{Berry_OD_86}}{: Berry, Matthew P. R., Christine M. Graham, Finlay W. McNab, Zhaohui Xu, Susannah A. A. Bloch, Tolu Oni, Katalin A. Wilkinson, et al. 2010. "An Interferon-Inducible Neutrophil-Driven Blood Transcriptional Signature in Human Tuberculosis." Nature 466 (7309): 973-77. \href{https://dx.doi.org/10.1038/nature09247}{10.1038/nature09247}}
#'  \item{\strong{Blankley_380}}{: Blankley, Simon, Christine M. Graham, Joe Levin, Jacob Turner, Matthew P. R. Berry, Chloe I. Bloom, Zhaohui Xu, et al. 2016. "A 380-Gene Meta-Signature of Active Tuberculosis Compared with Healthy Controls." The European Respiratory Journal: Official Journal of the European Society for Clinical Respiratory Physiology 47 (6): 1873-76. \href{https://dx.doi.org/10.1183/13993003.02121-2015}{10.1183/13993003.02121-2015}}
#'  \item{\strong{Blankley_5}}{: Blankley, Simon, Christine M. Graham, Joe Levin, Jacob Turner, Matthew P. R. Berry, Chloe I. Bloom, Zhaohui Xu, et al. 2016. "A 380-Gene Meta-Signature of Active Tuberculosis Compared with Healthy Controls." The European Respiratory Journal: Official Journal of the European Society for Clinical Respiratory Physiology 47 (6): 1873-76. \href{https://dx.doi.org/10.1183/13993003.02121-2015}{10.1183/13993003.02121-2015}}
#'  \item{\strong{Bloom_OD_144}}{: Bloom, Chloe I., Christine M. Graham, Matthew P. R. Berry, Fotini Rozakeas, Paul S. Redford, Yuanyuan Wang, Zhaohui Xu, et al. 2013. "Transcriptional Blood Signatures Distinguish Pulmonary Tuberculosis, Pulmonary Sarcoidosis, Pneumonias and Lung Cancers." PloS One 8 (8): e70630. \href{https://dx.doi.org/10.1371/journal.pone.0070630}{10.1371/journal.pone.0070630}}
#'  \item{\strong{Bloom_RES_268}}{: Bloom CI, Graham CM, Berry MP, et al. Detectable changes in the blood transcriptome are present after two weeks of antituberculosis therapy. PLoS One. 2012;7(10):e46191. \href{https://dx.doi.org/10.3389/fmicb.2021.650567}{10.3389/fmicb.2021.650567}}
#'  \item{\strong{Bloom_RES_558}}{: Bloom CI, Graham CM, Berry MP, et al. Detectable changes in the blood transcriptome are present after two weeks of antituberculosis therapy. PLoS One. 2012;7(10):e46191. \href{https://dx.doi.org/10.3389/fmicb.2021.650567}{10.3389/fmicb.2021.650567}}
#'  \item{\strong{Chen_HIV_4}}{: Chen Y, Wang Q, Lin S, et al. Meta-Analysis of Peripheral Blood Transcriptome Datasets Reveals a Biomarker Panel for Tuberculosis in Patients Infected With HIV. Front Cell Infect Microbiol. 2021;11:585919. Published 2021 Mar 19. \href{https://dx.doi.org/10.3389/fcimb.2021.585919}{10.3389/fcimb.2021.585919}}
#'  \item{\strong{Chendi_HIV_2}}{: Chendi BH, Tveiten H, Snyders CI, et al. CCL1 and IL-2Ra differentiate Tuberculosis disease from latent infection Irrespective of HIV infection in low TB burden countries [published online ahead of print, 2021 Jul 29]. J Infect. 2021;S0163-4453(21)00379-0. \href{https://dx.doi.org/10.1016/j.jinf.2021.07.036}{10.1016/j.jinf.2021.07.036}}
#'  \item{\strong{RISK11}}{: Darboe, F. et al. Diagnostic performance of an optimized transcriptomic signature of risk of tuberculosis in cryopreserved peripheral blood mononuclear cells. Tuberculosis 108, 124-126 (2018). \href{https://dx.doi.org/ 10.1016/j.tube.2017.11.001}{ 10.1016/j.tube.2017.11.001}}
#'  \item{\strong{Dawany_HIV_251}}{: Dawany, N. et al. Identification of a 251 gene expression signature that can accurately detect M. tuberculosis in patients with and without HIV co-infection. PLoS One 9, (2014). \href{https://dx.doi.org/10.1371/journal.pone.0089925}{10.1371/journal.pone.0089925}}
#'  \item{\strong{CMTB_CT}}{: Duffy FJ, Olson GS, Gold ES, Jahn A, Aderem A, Aitchison J, Rothchild AC, Diercks AH, Nemeth J. A contained Mycobacterium tuberculosis mouse infection model predicts active disease and containment in humans. The Journal of Infectious Diseases. 2021 Mar 10. \href{https://dx.doi.org/10.1093/infdis/jiab130}{10.1093/infdis/jiab130}}
#'  \item{\strong{Esmail_203}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964-73. \href{https://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Esmail_82}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964-73. \href{https://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Esmail_OD_893}}{: Esmail, Hanif, Rachel P. Lai, Maia Lesosky, Katalin A. Wilkinson, Christine M. Graham, Stuart Horswell, Anna K. Coussens, Clifton E. Barry 3rd, Anne O'Garra, and Robert J. Wilkinson. 2018. "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis." Proceedings of the National Academy of Sciences of the United States of America 115 (5): E964-73. \href{https://dx.doi.org/10.1073/pnas.1711853115}{10.1073/pnas.1711853115}}
#'  \item{\strong{Estevez_133}}{: Estévez O, Anibarro L, Garet E, et al. An RNA-seq Based Machine Learning Approach Identifies Latent Tuberculosis Patients With an Active Tuberculosis Profile. Front Immunol. 2020;11:1470. Published 2020 Jul 14. \href{https://dx.doi.org/10.3389/fimmu.2020.01470}{10.3389/fimmu.2020.01470}}
#'  \item{\strong{Estevez_259}}{: Estévez O, Anibarro L, Garet E, et al. An RNA-seq Based Machine Learning Approach Identifies Latent Tuberculosis Patients With an Active Tuberculosis Profile. Front Immunol. 2020;11:1470. Published 2020 Jul 14. \href{https://dx.doi.org/10.3389/fimmu.2020.01470}{10.3389/fimmu.2020.01470}}
#'  \item{\strong{Gjoen_10}}{: Gjøen, J.E., Jenum, S., Sivakumaran, D. et al. 'Novel transcriptional signatures for sputum-independent diagnostics of tuberculosis in children.' Sci Rep 7, 5839 (2017). \href{https://doi.org/10.1038/s41598-017-05057-x}{10.1038/s41598-017-05057-x}}
#'  \item{\strong{Gjoen_7}}{: Gjøen, J.E., Jenum, S., Sivakumaran, D. et al. 'Novel transcriptional signatures for sputum-independent diagnostics of tuberculosis in children.' Sci Rep 7, 5839 (2017). \href{https://doi.org/10.1038/s41598-017-05057-x}{10.1038/s41598-017-05057-x}}
#'  \item{\strong{Gliddon_2_OD_4}}{: Gliddon HD, Kaforou M, Alikian M, et al. Identification of Reduced Host Transcriptomic Signatures for Tuberculosis Disease and Digital PCR-Based Validation and Quantification. Front Immunol. 2021;12:637164. Published 2021 Mar 2. \href{https://dx.doi.org/10.3389/fimmu.2021.637164}{10.3389/fimmu.2021.637164}}
#'  \item{\strong{Gliddon_HIV_3}}{: Gliddon HD, Kaforou M, Alikian M, et al. Identification of Reduced Host Transcriptomic Signatures for Tuberculosis Disease and Digital PCR-Based Validation and Quantification. Front Immunol. 2021;12:637164. Published 2021 Mar 2. \href{https://dx.doi.org/10.3389/fimmu.2021.637164}{10.3389/fimmu.2021.637164}}
#'  \item{\strong{Gliddon_OD_3}}{: Gliddon, Harriet D., Kaforou, Myrsini, Alikian, Mary, Habgood-Coote, Dominic, Zhou, Chenxi, Oni, Tolu, Anderson, Suzanne T., Brent, Andrew J., Crampin, Amelia C., Eley, Brian, Kern, Florian, Langford, Paul R., Ottenhoff, Tom H. M., Hibberd, Martin L., French, Neil, Wright, Victoria J., Dockrell, Hazel M., Coin, Lachlan J., Wilkinson, Robert J., Levin, Michael. 2019 "Identification of reduced host transcriptomic signatures for tuberculosis and digital PCR-based validation and quantification" biorxiv.org: . \href{https://dx.doi.org/10.1101/583674}{10.1101/583674}}
#'  \item{\strong{Gliddon_OD_4}}{: Gliddon, Harriet D., Kaforou, Myrsini, Alikian, Mary, Habgood-Coote, Dominic, Zhou, Chenxi, Oni, Tolu, Anderson, Suzanne T., Brent, Andrew J., Crampin, Amelia C., Eley, Brian, Kern, Florian, Langford, Paul R., Ottenhoff, Tom H. M., Hibberd, Martin L., French, Neil, Wright, Victoria J., Dockrell, Hazel M., Coin, Lachlan J., Wilkinson, Robert J., Levin, Michael. 2019 "Identification of reduced host transcriptomic signatures for tuberculosis and digital PCR-based validation and quantification" biorxiv.org: . \href{https://dx.doi.org/10.1101/583674}{10.1101/583674}}
#'  \item{\strong{Gong_OD_4}}{: Gong Z, Gu Y, Xiong K, Niu J, Zheng R, Su B, Fan L and Xie J (2021) The Evaluation and Validation of Blood-Derived Novel Biomarkers for Precise and Rapid Diagnosis of Tuberculosis in Areas With High-TB Burden. Front. Microbiol. 12:650567. \href{https://dx.doi.org/10.3389/fmicb.2021.650567}{10.3389/fmicb.2021.650567}}
#'  \item{\strong{Heycken_FAIL_22}}{: Heyckendorf J, Marwitz S, Reimann M, et al. Prediction of anti-tuberculosis treatment duration based on a 22-gene transcriptomic model [published online ahead of print, 2021 Feb 11]. Eur Respir J. 2021;2003492. \href{https://dx.doi.org/10.1183/13993003.03492-2020}{10.1183/13993003.03492-2020}}
#'  \item{\strong{Hoang_OD_13}}{: Hoang, Long & Jain, Pooja & Pillay, Timesh & Tolosa-Wright, Mica & Niazi, Umar & Takwoingi, Yemisi & Halliday, Alice & Berrocal-Almanza, Luis & Deeks, Jonathan & Beverley, Peter & Kon, Onn & Lalvani, Ajit. (2021). Transcriptomic signatures for diagnosing tuberculosis in clinical practice: a prospective, multicentre cohort study. The Lancet Infectious Diseases. \href{https://dx.doi.org/10.1016/S1473-3099(20)30928-2}{10.1016/S1473-3099(20)30928-2}}
#'  \item{\strong{Hoang_OD_20}}{: Hoang, Long & Jain, Pooja & Pillay, Timesh & Tolosa-Wright, Mica & Niazi, Umar & Takwoingi, Yemisi & Halliday, Alice & Berrocal-Almanza, Luis & Deeks, Jonathan & Beverley, Peter & Kon, Onn & Lalvani, Ajit. (2021). Transcriptomic signatures for diagnosing tuberculosis in clinical practice: a prospective, multicentre cohort study. The Lancet Infectious Diseases. \href{https://dx.doi.org/10.1016/S1473-3099(20)30928-2}{10.1016/S1473-3099(20)30928-2}}
#'  \item{\strong{Hoang_OD_3}}{: Hoang, Long & Jain, Pooja & Pillay, Timesh & Tolosa-Wright, Mica & Niazi, Umar & Takwoingi, Yemisi & Halliday, Alice & Berrocal-Almanza, Luis & Deeks, Jonathan & Beverley, Peter & Kon, Onn & Lalvani, Ajit. (2021). Transcriptomic signatures for diagnosing tuberculosis in clinical practice: a prospective, multicentre cohort study. The Lancet Infectious Diseases. \href{https://dx.doi.org/10.1016/S1473-3099(20)30928-2}{10.1016/S1473-3099(20)30928-2}}
#'  \item{\strong{Huang_13}}{: Huang, Hai-Hui et al. 'Identification of 13 Blood-based Gene Expression Signatures to Accurately Distinguish Tuberculosis from Other Pulmonary Diseases and Healthy Controls'. 1 Jan. 2015 : S1837 - S1843.\href{https://doi.org/10.3233/BME-151486}{10.3233/BME-151486}}
#'  \item{\strong{Jacobsen_3}}{: Jacobsen, Marc, Dirk Repsilber, Andrea Gutschmidt, Albert Neher, Knut Feldmann, Hans J. Mollenkopf, Andreas Ziegler, and Stefan H. E. Kaufmann. 2007. "Candidate Biomarkers for Discrimination between Infection and Disease Caused by Mycobacterium Tuberculosis." Journal of Molecular Medicine  85 (6): 613-21. \href{https://dx.doi.org/10.1007/s00109-007-0157-6}{10.1007/s00109-007-0157-6}}
#'  \item{\strong{Jenum_8}}{: Jenum, S., Dhanasekaran, S., Lodha, R. et al. Approaching a diagnostic point-of-care test for pediatric tuberculosis through evaluation of immune biomarkers across the clinical disease spectrum. Sci Rep 6, 18520 (2016). \href{https://doi.org/10.1038/srep18520}{10.1038/srep18520}}
#'  \item{\strong{Kaforou_27}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{https://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kaforou_OD_44}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{https://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kaforou_OD_53}}{: Kaforou, Myrsini, Victoria J. Wright, Tolu Oni, Neil French, Suzanne T. Anderson, Nonzwakazi Bangani, Claire M. Banwell, et al. 2013. "Detection of Tuberculosis in HIV-Infected and -Uninfected African Adults Using Whole Blood RNA Expression Signatures: A Case-Control Study." PLoS Medicine 10 (10): e1001538. \href{https://dx.doi.org/10.1371/journal.pmed.1001538}{10.1371/journal.pmed.1001538}}
#'  \item{\strong{Kulkarni_HIV_2}}{: Kulkarni V, Queiroz ATL, Sangle S, et al. A Two-Gene Signature for Tuberculosis Diagnosis in Persons With Advanced HIV. Front Immunol. 2021;12:631165. Published 2021 Feb 22. \href{https://dx.doi.org/10.3389/fimmu.2021.631165}{10.3389/fimmu.2021.631165}}
#'  \item{\strong{LauxdaCosta_OD_3}}{: Laux da Costa L, Delcroix M, Dalla Costa ER, et al. A real-time PCR signature to discriminate between tuberculosis and other pulmonary diseases. Tuberculosis (Edinb). 2015;95(4):421-425. \href{https://dx.doi.org/10.1016/j.tube.2015.04.008}{10.1016/j.tube.2015.04.008}}
#'  \item{\strong{Lee_4}}{: Lee, Shih-Wei, Lawrence Shih-Hsin Wu, Guan-Mau Huang, Kai-Yao Huang, Tzong-Yi Lee, and Julia Tzu-Ya Weng. 2016. "Gene Expression Profiling Identifies Candidate Biomarkers for Active and Latent Tuberculosis." BMC Bioinformatics 17 Suppl 1 (January): 3. \href{https://dx.doi.org/10.1186/s12859-015-0848-x}{10.1186/s12859-015-0848-x}}
#'  \item{\strong{Leong_24}}{: Leong, Samantha, Yue Zhao, Noyal M. Joseph, Natasha S. Hochberg, Sonali Sarkar, Jane Pleskunas, David Hom, et al. 2018. "Existing blood transcriptional classifiers accurately discriminate active tuberculosis from latent infection in individuals from south India." Tuberculosis (109): 41-51. \href{https://doi.org/10.1016/j.tube.2018.01.002}{10.1016/j.tube.2018.01.002}}
#'  \item{\strong{PREDICT29}}{: Leong, S., Zhao, Y., Ribeiro-Rodrigues, R., Jones-López, E. C., Acuña-Villaorduña, C., Rodrigues, P. M., Palaci, M., Alland, D., Dietze, R., Ellner, J. J., Johnson, W. E., Salgame, P., Cross-validation of existing signatures and derivation of a novel 29-gene transcriptomic signature predictive of progression to TB in a Brazilian cohort of household contacts of pulmonary TB. Tuberculosis (Edinb). 2020 Jan;120:101898. \href{https://dx.doi.org/10.1016/j.tube.2020.101898}{10.1016/j.tube.2020.101898}}
#'  \item{\strong{Maertzdorf_15}}{: Maertzdorf J, McEwen G, Weiner J 3rd, et al. Concise gene signature for point-of-care classification of tuberculosis. EMBO Mol Med. 2016;8(2):86-95. \href{https://dx.doi.org/10.15252/emmm.201505790}{10.15252/emmm.201505790}}
#'  \item{\strong{DIAG4}}{: Maertzdorf, Jeroen, Gayle McEwen, January Weiner 3rd, Song Tian, Eric Lader, Ulrich Schriek, Harriet Mayanja-Kizza, Martin Ota, John Kenneth, and Stefan He Kaufmann. 2016. "Concise Gene Signature for Point-of-Care Classification of Tuberculosis." EMBO Molecular Medicine 8 (2): 86-95. \href{https://dx.doi.org/10.15252/emmm.201505790}{10.15252/emmm.201505790}}
#'  \item{\strong{Maertzdorf_OD_100}}{: Maertzdorf, Jeroen, January Weiner 3rd, Hans-Joachim Mollenkopf, TBornot TB Network, Torsten Bauer, Antje Prasse, Joachim Müller-Quernheim, and Stefan H. E. Kaufmann. 2012. "Common Patterns and Disease-Related Signatures in Tuberculosis and Sarcoidosis." Proceedings of the National Academy of Sciences of the United States of America 109 (20): 7853-58. \href{https://dx.doi.org/10.1073/pnas.1121072109}{10.1073/pnas.1121072109}}
#'  \item{\strong{RISK6}}{: Penn-Nicholson, A. et al. RISK6, a 6-gene transcriptomic signature of TB disease risk, diagnosis and treatment response. Sci. Rep. 10, (2020). \href{https://dx.doi.org/10.1038/s41598-020-65043-8}{10.1038/s41598-020-65043-8}}
#'  \item{\strong{Qian_OD_17}} {: Qian, Zhongqing et al. "Expression of nuclear factor, erythroid 2-like 2-mediated genes differentiates tuberculosis." Tuberculosis (Edinburgh, Scotland) vol. 99 (2016): 56-62. \href{https://doi.org/10.1016/j.tube.2016.04.008}{10.1016/j.tube.2016.04.008}}
#'  \item{\strong{Rajan_HIV_5}}{: Rajan, Jayant V., Semitala, Fred C., Kamya, Moses R., Yoon, Christina., Mehta, Tejas., Cattamanchi, Adithya., Seielstad, Mark., Montalvo, Lani., Andama, Alfred., Katende, Jane., Asege, Lucy., Nakaye, Martha., Mwebe, Sandra. 2018 "A Novel, 5-Transcript, Whole-blood Gene-expression Signature for Tuberculosis Screening Among People Living With Human Immunodeficiency Virus" Clinical Infectious Diseases: 1-7. \href{https://doi.org/10.1093/cid/ciy835}{10.1093/cid/ciy835}}
#'  \item{\strong{Roe_3}}{: Roe, Jennifer, Venturini, Cristina, Gupta, Rishi K., Gurry, Celine, Chain, Benjamin M., Sun, Yuxin, Southern, Jo, Jackson, Charlotte, Lipman, Marc, C., Miller, Robert F., Martineau, Adrian R., Abubakar, Ibrahim, Noursadeghi, Mahdad. 2019 "T1 Blood transcriptomic stratification of short-term risk in contacts of tuberculosis": . \href{https://doi.org/10.1093/cid/ciz252}{10.1093/cid/ciz252}}
#'  \item{\strong{Roe_OD_4}}{: Roe, Jennifer K., Niclas Thomas, Eliza Gil, Katharine Best, Evdokia Tsaliki, Stephen Morris-Jones, Sian Stafford, et al. 2016. "Blood Transcriptomic Diagnosis of Pulmonary and Extrapulmonary Tuberculosis." JCI Insight 1 (16): e87238. \href{https://dx.doi.org/10.1172/jci.insight.87238}{10.1172/jci.insight.87238}}
#'  \item{\strong{Sambarey_HIV_10}}{: Sambarey, Awanti, Abhinandan Devaprasad, Abhilash Mohan, Asma Ahmed, Soumya Nayak, Soumya Swaminathan, George D'Souza, et al. 2017. "Unbiased Identification of Blood-Based Biomarkers for Pulmonary Tuberculosis by Modeling and Mining Molecular Interaction Networks." EBioMedicine 15 (February): 112-26. \href{https://dx.doi.org/10.1016/j.ebiom.2016.12.009}{10.1016/j.ebiom.2016.12.009}}
#'  \item{\strong{Singhania_OD_20}}{: Singhania, Akul, Raman Verma, Christine M. Graham, Jo Lee, Trang Tran, Matthew Richardson, Patrick Lecine, et al. 2018. "A Modular Transcriptional Signature Identifies Phenotypic Heterogeneity of Human Tuberculosis Infection." Nature Communications 9 (1): 2308. \href{https://dx.doi.org/10.1038/s41467-018-04579-w}{10.1038/s41467-018-04579-w}}
#'  \item{\strong{Sivakumaran_11}}{: Sivakumaran D, Ritz C, Gjøen JE, et al. Host Blood RNA Transcript and Protein Signatures for Sputum-Independent Diagnostics of Tuberculosis in Adults. Front Immunol. 2021;11:626049. Published 2021 Feb 4. \href{https://dx.doi.org/10.3389/fimmu.2020.626049}{10.3389/fimmu.2020.626049}}
#'  \item{\strong{Sloot_HIV_2}}{: Sloot, Rosa, Maarten F. Schim van der Loeff, Erik W. van Zwet, Mariëlle C. Haks, Sytze T. Keizer, Maarten Scholing, Tom H. M. Ottenhoff, Martien W. Borgdorff, and Simone A. Joosten. 2015. "Biomarkers Can Identify Pulmonary Tuberculosis in HIV-Infected Drug Users Months Prior to Clinical Diagnosis." EBioMedicine 2 (2): 172-79. \href{https://dx.doi.org/10.1016/j.ebiom.2014.12.001}{10.1016/j.ebiom.2014.12.001}}
#'  \item{\strong{Suliman_4}}{: Suliman, Sara, Ethan Thompson, Jayne Sutherland, January Weiner Rd, Martin O. C. Ota, Smitha Shankar, Adam Penn-Nicholson, et al. 2018. "Four-Gene Pan-African Blood Signature Predicts Progression to Tuberculosis." American Journal of Respiratory and Critical Care Medicine, April. https://doi.org/10.1164/rccm.201711-2340OC. \href{https://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{Suliman_RISK_2}}{: Suliman, S. et al. Four-gene pan-African blood signature predicts progression to tuberculosis. Am. J. Respir. Crit. Care Med. 197, 1198-1208 (2018). \href{https://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{RISK4}}{: Suliman, Sara, Ethan Thompson, Jayne Sutherland, January Weiner Rd, Martin O. C. Ota, Smitha Shankar, Adam Penn-Nicholson, et al. 2018. "Four-Gene Pan-African Blood Signature Predicts Progression to Tuberculosis." American Journal of Respiratory and Critical Care Medicine, April. https://doi.org/10.1164/rccm.201711-2340OC. \href{https://dx.doi.org/10.1164/rccm.201711-2340OC}{10.1164/rccm.201711-2340OC}}
#'  \item{\strong{DIAG3}}{: Sweeney, Timothy E., Lindsay Braviak, Cristina M. Tato, and Purvesh Khatri. 2016. "Genome-Wide Expression for Diagnosis of Pulmonary Tuberculosis: A Multicohort Analysis." The Lancet. Respiratory Medicine 4 (3): 213-24. \href{https://dx.doi.org/10.1016/S2213-2600(16)00048-5}{10.1016/S2213-2600(16)00048-5}}
#'  \item{\strong{DISEASE}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48-58. \href{https://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{FAILURE}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48-58. \href{https://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{RESPONSE5}}{: Thompson, Ethan G., Ying Du, Stephanus T. Malherbe, Smitha Shankar, Jackie Braun, Joe Valvo, Katharina Ronacher, et al. 2017. "Host Blood RNA Signatures Predict the Outcome of Tuberculosis Treatment." Tuberculosis  107 (December): 48-58. \href{https://dx.doi.org/10.1016/j.tube.2017.08.004}{10.1016/j.tube.2017.08.004}}
#'  \item{\strong{Tornheim_71}}{: Tornheim, Jeffrey A., Anil K. Madugundu, Mandar Paradkar, Kiyoshi F. Fukutani, Artur TL Queiroz, Nikhil Gupte, Akshay N. Gupte et al. 2020. "Transcriptomic Profiles of Confirmed Pediatric Tuberculosis Patients and Household Contacts Identifies Active Tuberculosis, Infection, and Treatment Response Among Indian Children." The Journal of Infectious Diseases 221(10): 1647-1658. \href{https://doi.org/10.1093/infdis/jiz639}{10.1093/infdis/jiz639}}
#'  \item{\strong{Tornheim_RES_25}}{: Tornheim, Jeffrey A., Anil K. Madugundu, Mandar Paradkar, Kiyoshi F. Fukutani, Artur TL Queiroz, Nikhil Gupte, Akshay N. Gupte et al. 2020. "Transcriptomic Profiles of Confirmed Pediatric Tuberculosis Patients and Household Contacts Identifies Active Tuberculosis, Infection, and Treatment Response Among Indian Children." The Journal of Infectious Diseases 221(10): 1647-1658. \href{https://doi.org/10.1093/infdis/jiz639}{10.1093/infdis/jiz639}}
#'  \item{\strong{Verhagen_10}}{: Verhagen, L.M., Zomer, A., Maes, M. et al. A predictive signature gene set for discriminating active from latent tuberculosis in Warao Amerindian children. BMC Genomics 14, 74 (2013). \href{https://doi.org/10.1186/1471-2164-14-74}{10.1186/1471-2164-14-74}}
#'  \item{\strong{Walter_51}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274-82. \href{https://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Walter_PNA_119}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274-82. \href{https://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{Walter_PNA_47}}{: Walter, Nicholas D., Mikaela A. Miller, Joshua Vasquez, Marc Weiner, Adam Chapman, Melissa Engle, Michael Higgins, et al. 2016. "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: A Case-Control Study with Systematic Cross-Classifier Evaluation." Journal of Clinical Microbiology 54 (2): 274-82. \href{https://dx.doi.org/10.1128/JCM.01990-15}{10.1128/JCM.01990-15}}
#'  \item{\strong{ACS_COR}}{: Zak, Daniel E., Adam Penn-Nicholson, Thomas J. Scriba, Ethan Thompson, Sara Suliman, Lynn M. Amon, Hassan Mahomed, et al. 2016. "A Blood RNA Signature for Tuberculosis Disease Risk: A Prospective Cohort Study." The Lancet 387 (10035): 231222. \href{https://dx.doi.org/10.1016/S0140-6736(15)01316-1}{10.1016/S0140-6736(15)01316-1}}
#'  \item{\strong{Zimmer_RES_3}}{: Zimmer, A.J., Schumacher, S.G., Södersten, E. et al. A novel blood-based assay for treatment monitoring of tuberculosis. BMC Res Notes 14, 247 (2021). \href{https://dx.doi.org/10.1186/s13104-021-05663-z}{10.1186/s13104-021-05663-z}}
#'  }
#'
#' @keywords datasets
#' @examples
#' data("TBcommon")
"TBcommon"


#' An example TB dataset with Indian population data.
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
#' as part of an ongoing household contact study. Whole blood RNA-Seq analysis
#' was performed on all 44 participants.
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


#' An example TB dataset with TB/HIV data.
#'
#' An example dataset containing the gene expression and metadata in a
#' SummarizedExperiment object for 31 subjects with HIV and/or Tuberculosis
#' diseases. Information on subject infection status can be accessed with
#' \code{TB_hiv$Disease}. Samples with both TB and HIV contamination are
#' marked as \code{tb_hiv}, while samples with HIV and no TB are marked
#' as \code{hiv_only}.
#'
#' This dataset was published as part of a study to assess whether gene expression
#' signatures and cytokine levels would distinguish active TB in advanced HIV
#' in a cohort residing in Sub-Saharan Africa (Verma et. al 2018).
#' Participants were severely immunosuppressed TB-HIV patients who had
#' not yet received TB treatment or anti-retroviral therapy (ART). The dataset included
#' in this package has been lightly edited from the originally published dataset
#' due to the removal of one participant who was HIV positive, on ART and developed
#' TB during follow-up. Whole blood RNA-Seq analysis was performed on all
#' 31 participants.
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

#' Annotation information for published TB signatures.
#'
#' A data.frame of annotation information for published tuberculosis signatures.
#' Currently, this table includes two variables, \code{disease} and
#' \code{tissue type}.
#'
#' The \code{disease} variable indicates whether the signature was developed
#' to distinguish TB from LTBI ("Disease"), TB from some combination of other
#' diseases and possibly LTBI ("OD"), TB from Human Immunodeficiency Virus ("HIV"),
#' TB from pneumonia ("PNA"), or identify risk of progression to TB ("RISK"),
#' risk of TB treatment failure ("FAIL"), or classify treatment responses
#' (i.e., failures from cures, "RES").
#'
#' The \code{tissue type} variable denotes whether the signature was developed
#' using samples of either whole blood/paxgene or peripheral blood mononuclear
#' cells (PBMCs). Due to the manipulation of cells inherently required to obtain
#' PBMCs, many scientists prefer to use only whole blood samples for analysis.
#'
#' @name sigAnnotData
#' @docType data
#' @format data.frame
#' @keywords datasets
#' @source
#' See \code{?TBsignatures} for reference information.
#'
#' @examples
#' data("sigAnnotData")
"sigAnnotData"

#' Annotation information for published TB signatures.
#'
#' A \code{data.frame} of annotation information for published tuberculosis
#' signatures. This table differs from that of \code{sigAnnotData} as it
#' refers to signatures via the name given in scientific publications, and
#' via a consistent naming system otherwise.
#' Currently, this table includes two variables, \code{disease} and
#' \code{tissue type}.
#'
#' The \code{disease} variable indicates whether the signature was developed
#' to distinguish TB from LTBI ("Disease"), TB from some combination of other
#' diseases and possibly LTBI ("OD"), TB from Human Immunodeficiency Virus ("HIV"),
#' TB from pneumonia ("PNA"), or identify risk of progression to TB ("RISK"),
#' risk of TB treatment failure ("FAIL"), or classify treatment responses
#' (i.e., failures from cures, "RES").
#'
#' The \code{tissue type} variable denotes whether the signature was developed
#' using samples of either whole blood/paxgene or peripheral blood mononuclear
#' cells (PBMCs). Due to the manipulation of cells inherently required to obtain
#' PBMCs, many scientists prefer to use only whole blood samples for analysis.
#'
#' @name common_sigAnnotData
#' @docType data
#' @format data.frame
#' @keywords datasets
#' @source
#' See \code{?TBcommon} for reference information.
#'
#' @keywords datasets
#' @examples
#' data("common_sigAnnotData")
"common_sigAnnotData"

#' A list of published/pre-print COVID-19 signatures.
#'
#' A set of 47 COVID-19 gene signatures from various single-cell and bulk
#' RNA-seq publications and preprint manuscripts from early- to mid-2020.
#' This set of signatures uses gene symbols.
#'
#' Signature names are composed of the last name of the primary author, followed
#' by the type of sequencing data from which the signature was derived, the
#' tissue from which the signature was derived, and ending with a reference to
#' the cell type, infection status, or disease to which the signature belongs,
#' as defined in the original publication.
#'
#' Note that in some cases signatures will be positive identifiers of COVID-19
#' as positive markers of immune cell clusters are often provided for
#' single-cell RNA-seq data; this should be taken into account when creating
#' ROC curves and computing any AUC or disease risk estimates.
#'
#' @name COVIDsignatures
#' @docType data
#' @format list
#' @source
#' \itemize{
#'  \item{\strong{Wilk_sc_PBMC_monocytes_up}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076. {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Wilk_sc_PBMC_monocytes_up}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076.  {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Wilk_sc_PBMC_monocytes_down}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076.  {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Wilk_sc_PBMC_NK_cells_up}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076.  {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Wilk_sc_PBMC_NK_cells_down}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076.  {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Wilk_sc_PBMCs_ISG_signature}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076.  {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Wilk_sc_PBMC_activated_granulocytes}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "A single-cell atlas of the peripheral immune response in patients with severe COVID-19." Nature Medicine 26 (7): 1070-1076.  {https://doi.org/10.1038/s41591-020-0944-y}}
#'  \item{\strong{Huang_sc_PBMC_IFN_signature}}{: Wilk, A.J., Rustagi, A., Zhao, N.Q. et al. 2020. "Blood single cell immune profiling reveals the interferon-MAPK pathway mediated adaptive immune response for COVID-19." medRxiv.org: {https://doi.org/10.1101/2020.03.15.20033472}}
#'  \item{\strong{Wen_sc_PBMC_monocytes}}{: Wen, W., Su, W., Tang, H. et al. 2020. "Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing." Cell Discovery 6 (31).  {https://doi.org/10.1038/s41421-020-0168-9}}
#'  \item{\strong{Wen_sc_PBMC_NK_cells}}{: Wen, W., Su, W., Tang, H. et al. 2020. "Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing." Cell Discovery 6 (31).  {https://doi.org/10.1038/s41421-020-0168-9}}
#'  \item{\strong{Wen_sc_PBMC_CD4_T_cells}}{: Wen, W., Su, W., Tang, H. et al. 2020. "Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing." Cell Discovery 6 (31).  {https://doi.org/10.1038/s41421-020-0168-9}}
#'  \item{\strong{Wen_sc_PBMC_CD8_T_cells}}{: Wen, W., Su, W., Tang, H. et al. 2020. "Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing." Cell Discovery 6 (31).  {https://doi.org/10.1038/s41421-020-0168-9}}
#'  \item{\strong{Wen_sc_PBMC_B_cells}}{: Wen, W., Su, W., Tang, H. et al. 2020. "Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing." Cell Discovery 6 (31).  {https://doi.org/10.1038/s41421-020-0168-9}}
#'  \item{\strong{Xiong_bulk_PBMC_gene_signature_up}}{: Xiong Y, Liu Y, Cao L, et al. 2020. "Transcriptomic characteristics of bronchoalveolar lavage fluid and peripheral blood mononuclear cells in COVID-19 patients." Emerging Microbes & Infections 9 (1):761-770.  {https://doi/org/10.1080/22221751.2020.1747363}}
#'  \item{\strong{Xiong_bulk_PBMC_gene_signature_down}}{: Xiong Y, Liu Y, Cao L, et al. 2020. "Transcriptomic characteristics of bronchoalveolar lavage fluid and peripheral blood mononuclear cells in COVID-19 patients." Emerging Microbes & Infections 9 (1):761-770.  {https://doi/org/10.1080/22221751.2020.1747363}}
#'  \item{\strong{Xiong_sc_PBMC_cytokines_up}}{: Xiong Y, Liu Y, Cao L, et al. 2020. "Transcriptomic characteristics of bronchoalveolar lavage fluid and peripheral blood mononuclear cells in COVID-19 patients." Emerging Microbes & Infections 9 (1):761-770.  {https://doi/org/10.1080/22221751.2020.1747363}}
#'  \item{\strong{Xiong_sc_PBMC_cytokines_down}}{: Xiong Y, Liu Y, Cao L, et al. 2020. "Transcriptomic characteristics of bronchoalveolar lavage fluid and peripheral blood mononuclear cells in COVID-19 patients." Emerging Microbes & Infections 9 (1):761-770.  {https://doi/org/10.1080/22221751.2020.1747363}}
#'  \item{\strong{Liao_sc_BALF_G1_macrophages}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Liao_sc_BALF_G1_2_macrophages}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Liao_sc_BALF_G2_macrophages}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Liao_sc_BALF_G3_macrophages}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Liao_sc_BALF_G4_macrophages}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Liao_sc_BALF_CD8_T_cells_up}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Liao_sc_BALF_CD8_T_cells_down}}{: Liao, M., Liu, Y., Yuan, J. et al. 2020. "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19." Nature Medicine 26 (6): 842-844.  {https://doi.org/10.1038/s41591-020-0901-9}}
#'  \item{\strong{Hadjadj_nanostring_WB_gene_signature_up}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_WB_gene_signature_down}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_WB_ISG_signature}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_WB_mild_moderate_up}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_WB_mild_moderate_down}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_WB_severe_up}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_WB_severe_down}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_critical_up}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Hadjadj_nanostring_critical_down}}{: Hadjadj J, Yatim N, Barnabei L, et al. 2020. "Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients." Science 369 (6504): 718-724.  {https://doi.org/10.1126/science.abc6027}}
#'  \item{\strong{Wei_sc_PBMC_inactivated_monocytes}}{: Wei et al. 2020. "Viral Invasion and Type I Interferon Response Characterize the Immunophenotypes during COVID-19 Infection." SSRN: {https://dx.doi.org/10.2139/ssrn.3555695}}
#'  \item{\strong{Wei_sc_PBMC_classical_monocytes}}{: Wei et al. 2020. "Viral Invasion and Type I Interferon Response Characterize the Immunophenotypes during COVID-19 Infection." SSRN: {https://dx.doi.org/10.2139/ssrn.3555695}}
#'  \item{\strong{Wei_sc_PBMCs_T_cells}}{: Wei et al. 2020. "Viral Invasion and Type I Interferon Response Characterize the Immunophenotypes during COVID-19 Infection." SSRN: {https://dx.doi.org/10.2139/ssrn.3555695}}
#'  \item{\strong{Wei_sc_PBMC_B_cells}}{: Wei et al. 2020. "Viral Invasion and Type I Interferon Response Characterize the Immunophenotypes during COVID-19 Infection." SSRN: {https://dx.doi.org/10.2139/ssrn.3555695}}
#'  \item{\strong{Silvin_sc_WB_combined_signature}}{: Silvin A, Chapuis N, Dunsmore G, et al. 2020. "Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19." Cell 182 (6): 1401-1418.E18.  {https://doi.org/10.1016/j.cell.2020.08.002}}
#'  \item{\strong{Silvin_sc_WB_monocytes_up}}{: Silvin A, Chapuis N, Dunsmore G, et al. 2020. "Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19." Cell 182 (6): 1401-1418.E18.  {https://doi.org/10.1016/j.cell.2020.08.002}}
#'  \item{\strong{Silvin_sc_WB_monocytes_down}}{: Silvin A, Chapuis N, Dunsmore G, et al. 2020. "Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19." Cell 182 (6): 1401-1418.E18.  {https://doi.org/10.1016/j.cell.2020.08.002}}
#'  \item{\strong{Silvin_sc_WB_neutrophils_up}}{: Silvin A, Chapuis N, Dunsmore G, et al. "Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19." Cell 182 (6): 1401-1418.E18.  {https://doi.org/10.1016/j.cell.2020.08.002}}
#'  \item{\strong{Silvin_sc_WB_neutrophils_down}}{: Silvin A, Chapuis N, Dunsmore G, et al. "Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19." Cell 182 (6): 1401-1418.E18.  {https://doi.org/10.1016/j.cell.2020.08.002}}
#'  \item{\strong{Arunachalam_bulk_PBMC_blood_modules}}{: Arunachalam PS, Wimmers F, Mok CKP, et al. 2020. "Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans." Science 369 (6508): 1210-1220.  {https://doi.org/10.1126/science.abc6261}}
#'  \item{\strong{Arunachalam_bulk_PBMC_covid_combined}}{: Arunachalam PS, Wimmers F, Mok CKP, et al. 2020. "Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans." Science 369 (6508): 1210-1220.  {https://doi.org/10.1126/science.abc6261}}
#'  \item{\strong{Arunachalam_bulk_PBMC_moderate}}{: Arunachalam PS, Wimmers F, Mok CKP, et al. 2020. "Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans." Science 369 (6508): 1210-1220.  {https://doi.org/10.1126/science.abc6261}}
#'  \item{\strong{Arunachalam_bulk_PBMC_severe}}{: Arunachalam PS, Wimmers F, Mok CKP, et al. 2020. "Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans." Science 369 (6508): 1210-1220.  {https://doi.org/10.1126/science.abc6261}}
#'  \item{\strong{Arunachalam_bulk_PBMC_intensive_care}}{: Arunachalam PS, Wimmers F, Mok CKP, et al. 2020. "Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans." Science 369 (6508): 1210-1220.  {https://doi.org/10.1126/science.abc6261}}
#'  \item{\strong{Dunning_bulk_WB_flu}}{: Dunning J, Blankley S, Hoang LT, et al. 2018. "Progression of whole-blood transcriptional signatures from interferon-induced to neutrophil-associated patterns in severe influenza." Nature Immunology 19 (6): 625-635.  {https://doi.org/10.1038/s41590-018-0111-5}}
#'  }
#'
#' @keywords datasets
#' @examples
#' data("COVIDsignatures")

"COVIDsignatures"

#' Up/Down-regulated genes information for selected TB signatures.
#' @name TBsignaturesSplit
#' @docType data
#' @format list
#' @source
#' See \code{?TBsignatures} for reference information.
#' @keywords datasets
#' @examples
#' data("TBsignaturesSplit")
"TBsignaturesSplit"

#' Discovery datasets for corresponding gene signatures.
#' @name OriginalTrainingData
#' @docType data
#' @format list
#' @source
#' See \code{?TBsignatures} for reference information.
#' @keywords datasets
#' @examples
#' data("OriginalTrainingData")
"OriginalTrainingData"

