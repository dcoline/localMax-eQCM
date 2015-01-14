# localMax-eQCM



THESIS


Presented in Partial Fulfillment of the Requirements for the Degree Master of Science in the Graduate School of The Ohio State University

By
Daniel C. Morgan
The Ohio State University
2014

Master's Examination Committee:
Kun Huang, Ph.D. Thesis Advisor
Albert Lai, Ph.D. Committee Member
James Chen, M.D. Committee Member
Peter Embi, M.D. Committee Member 
















Copyright by
Daniel C. Morgan
2014



 


Abstract
Biomarkers are the actionable factors differentiating any condition states, here disease response and non-response to treatment. They have been used as screening measures to indicate surrogate or subclinical manifestations, and are also prime targets of drug repositioning approaches. We aim to identify gene sub-networks involved in response to drug by implementing a data-mining algorithm that identifies genes with highly correlated expression across the population, concomitantly significantly variable between response and non-response cohorts to allow differentiation. Our greedy Quasi-Clique Merger implementation, localMax-eQCM, is tuned to identifying unique sub-network modules by restricting initial gene pairs to only those of highest local co-expression, a rare event, thereafter keeping with the original QCM, that each and every gene within the sub-network be co-expressed to a high degree (Figure 2.3). We are then able to differentiate differentially expressed modules between response subpopulations by t-test among respective eigengene. This method of biomarker evaluation could ultimately be employed to assign patients to receive chemotherapy, predicted to productively response based upon biomarker screening. We implement our algorithm across CCLE gene expression from a NSCLC large and squamous cell subpopulation, of which cell lines were immortalized and subsequent Taxol response profile (IC50) measured. We identified three sub-network modules highly correlated to Taxol response in this way.











Dedication

For Melvin












Vita

2006	St. Ignatius HS
2011	B.S. Microbiology, Miami University



Fields of Study

Major Field:  Public Health

 


Table of Contents
Abstract	ii
Dedication	iii
Vita	iv
Fields of Study	iv
Table of Contents	v
List of Tables	vii
List of Equations	vii
List of Figures	viii
Chapter 1:  Background & Significance	1
1.1 Individual Variation in Drug Response	1
1.2 Cancer Cell Line Encyclopedia and Non-Small Cell Lung Cancer	2
1.3 Co-expression Analysis & Data Mining	4
1.4 Overview of  QCM Workflow	5
Chapter 2: Workflow for Gene Co-expression Network Analysis on CCLE data	8
2.1 Weighted Gene Co-expression Network Analysis	8
2.2 WGCN Mining vs Clustering	10
2.3 Quasi-Clique Merger Algorithm	13
2.4 Review before Differentiation	19
Chapter 3: Differentially Expressed Quasi-Cliques	21
3.1 Quasi-Clique Formation & Merging	21
3.2 Quasi-Clique Differential Expression	22
Quasi-Clique 12	25
Quasi-Clique 34	25
Quasi-Clique 126	25
Chapter 4: Discussion and Future Work	26
 


List of Tables

Table 1.1 Phase II Trials of Single-Agent Taxol8.	7
Table 2.1 Comparison of Correlation metrics.	9
Table 3.1. Significant sub-networks acting in response to Taxol	24
Table 4.1 Continued Study	27

List of Equations
Equation 2.1 The PCC For any given pair of genes (X,Y)	8
Equation 2.2 Density of Subnetwork	11
Equation 2.3 Criteria for addition	15
Equation 2.4 Expanded criteria for addition	15
Equation 2.5 an is bound by K.	15
Equation 2.6 Density must remain above threshold parameter β	18
Figure 2.6 Merging sub-network matrixes	18
 

List of Figures
Figure 2.1 Adjacency and Distance matrices, resultant weighted graph.	10
Figure 2.2 A set of vertices S is called a clique	11
Figure 2.3 Evolution of sub-network creation	13
Figure 2.4 Evolution of QCM. The previous implementation of eQCM	14
Figure 2.6 Merging sub-network matrixes	18
 
 

Chapter 1:  Background & Significance
Drug repositioning matches underutilized treatments to applications for which they were not designed. We propose a method to screen for cancer patients responsive to chemotherapy through differential expression analysis among subpopulations of variable response to a treatment. The presence or absence of certain pathways dictates the measured response variability, and we aim to identify the gene sub-networks involved using the novel starting point of maximum neighboring co-expression, from which unique sub-networks fitting a number of criteria are assembled.
1.1 Individual Variation in Drug Response 
	Bioinformatics approaches have lead to a transition from the traditional study of individual disease-gene associations to a more physiologically relevant interest in the interaction among groups of genes. We are concerned with the comparison of the expression of gene groups within patients whose immortalized cells have been measured to be responsive to drug versus those characterized as nonresponsive. We hypothesize that such gene groups, with relationships quantified and thusly represented as sub-networks, will illuminate the underlying mechanisms responsible for the variable response. 
	It is undeniable that just as patients vary in response to drug, they vary in adverse reaction to drug therapy1. We identify gene networks calculated to be significantly associated with response to a drug, with the intention that continuing research might lead to new clinical decision support systems (CDSS) that aid clinicians in choosing productive therapies. The foresight of gene networks predicted to be essential for positive response to the any given therapeutic compound could be utilized in determining most likely productive therapy to apply to patient presenting any given set of gene modules. In an extension of drug repositioning, we apply a data mining approach to matching the drug back to the patient’s disease response subpopulation. This would allow clinicians of patients presenting requisite modules to go ahead with otherwise ‘suboptimal at population level’ treatments, as well as to seek other courses of treatment for patients lacking the respective biomarker(s), mitigating side effects of an unproductive therapy. The focus of this thesis is to determine if samples responsive to chemotherapy present such discriminating gene sub-network(s) absent in non-responsive samples through an approach whose novelty lies in the strict calculation of unique, highly co-expressed starting gene pairs from which biologically similar sub-networks are assembled. 

1.2 Dataset: Cancer Cell Line Encyclopedia: Large & Squamous Cell Non-Small Cell Lung Cancer
Non-small-cellular lung cancers (NSCLC) are a non-uniform group of tumors regarding the histological structure and clinical course. Majority constituents’ squamous cell carcinoma and adenocarcinoma are characterized by ill-defined prognosis, which pose problems in the selection of effective post-operative therapy. Descriptive biological factors that could help in choosing an effective method of treatment after primary tumor removal have been searched for for several years2. We aim to identify multiple such gene sub-network modules that is targetable by chemotherapeutic regimen in a pre-determined responsive subpopulation.
We utilized 66 NSCLC immortalized cell lines’ gene expression dataset contained within the Cancer Cell Line Encyclopedia (CCLE) to calculate biologically similar gene modules from their correlated expression patterns, a subset of which were differentially expressed among response subpopulations, indicative of functionality in response to chemotherapy. We utilized the Pearson correlation coefficient (PCC) rather than the non-linear, non-Gaussian Spearman metric of co-expression to sub-network genes for its proven track record in microarray analysis (table 2.1). All 947 human cancer cell lines have been characterized by gene expression using Affymetrix U133 plus 2.0 array, chromosomal copy number using high-density Affymetrix SNP 6.0 array, and mutational status using massively parallel sequencing data and mass spectrometric genotyping. To ensure the validity of these measures, each of these metrics were previously compared to those of cell lineages in other databases, which showed strong positive correlation overall3. 
The same cannot be said of the pharmacological profiling4, measures of dose-response gathered in vitro for 24 anticancer drugs against a subset of 479 cell lines, recorded as IC50 values. However, the pharmacological profiles of the NSCLC subset encompassing squamous and large cell carcinoma enabled us to stratify patient populations based on calculated IC50. We calculated our quasi-cliques and the eigengenes of their representative expression based on this 66 large and squamous cell NSCLC cell line dataset, but only used the top and bottom quartile for differentiating populations in the unpaired two-sample t-test. This study aims to return predictions of which modules’ expression constitute a positive response to chemotherapy.

1.3 Co-expression Analysis & Data Mining
Physiology is rarely dictated by singular components, whether of ecological interaction or biochemical pathway. More approximate relationships can best be represented using statistics, here correlation, stored in matrixes, transformed and mined by clinical questioning, and visualized using graphs. Here, data mining of the correlation matrix allows for the identification of sub-network modules of some specific, shared characteristic, which can then be differentiated using drug profile data; here sub-network density as a function of correlation is used to mine similar genes, which modules’ expression are then differentiated among subpopulations. Implementations of clustering have been used in social science, networking, transportation, image processing and bioinformatics, offering insight to highly interconnected relationships, among other things. Biologically relevant examples include modeling metabolic networks, transcriptional regulatory networks, and protein-protein interaction networks via two-hybrid screening. The two-hybrid system works to identify virtually all protein-protein interaction by measuring reporter gene transcription, a product of interaction between eukaryotic transcription factors’ binding and activation domains. Beyond physical interaction, detection techniques based on genomic information can elucidate such important relationships as co-expression, co-localization, pathway prediction, and shared protein domain. However, these methods have the drawbacks of exhaustive inclusion into clusters and non-gene overlap between clusters. Mining correlated RNA expression is based upon the widely supported hypothesis that similar gene expression is due to functional association5, and was used here under such an assumption. The large amount of publically available microarray data makes this approach stand out because it allows functional associations and relationships (conceptual networks6) to be formed based on patterns otherwise hidden within the data. These patterns are mined by mathematical weighting of similar expression. This allows for unique modules formation of only highly similar elements, and does not restrict the number of sub-networks any gene can find itself comprising, reflective of a nature physiology where any gene might have many function.
	Correlated expression begins by representing any gene pair as a graph G=(V,E) formed by a set of vertices, V={v1,v2,…,,vN} connected by a unique Pearson correlation value for each pair of genes as a set of edges, E= {e1, e2,…,eM}7. Initial pre-processing is often required to format the data from the raw expression output. We pulled GC-RMA normalized read values of 18,898 unique genes characterized within the CCLE (figure 2.4, step 1). We selected genes of high variability, discarding the lowest 15%, returning a pool of 16140 genes. We then performed calculations comparing average expression of genes, after transformation by adjacency function (AF=1-|PCC|), used here as the measure of similarity between the genes, the weight of the edges connecting the vertices, w(e), transformed to length. This process creates a diagonally reflective distance matrix in which every gene’s is represented by a node/vertex, and expression is correlated between every gene-gene pair’s AF-transformed co-expression. It is from this point that dense sub-networks can be mined, differentiated upon, and ultimately association to drug response strata can be characterized. 
1.4 Overview of QCM Workflow
	We improve upon machine learning based clustering efforts exhaustively incorporating every gene into a cluster while not allowing genes to hold position in multiple clusters (and thus hold multiple functions). We also improve upon previous implementations of the Quasi-Clique Merger (QCM) algorithm, with initial focuses in node strength and edge length, in order to identify gene networks associated with response to Taxol in large and squamous cell NSCLC cell lines contained within the CCLE. We begin by ranking gene expression in order to capture the highest quality reads and discard those of lowest 15% variability to ensure a Pearson correlation coefficient (PCC) indicative of high co-expression, not simply low variance (Equation 2.1). We subsequently implement our localMax-eQCM, with the added initial step of restricting potential maximum edges to those of highest correlation only. For there we follow the eQCM, with a γ threshold step to evaluate neighboring edge weight so to begin propagation from a certain level of co-expressed gene-pairs, thereby guaranteeing unique sub-network formation. Two other threshold steps invoking λ and t (ie K), and β ensure genes added to this initial pair result in dense sub-networks, and that merging small sub-network of sufficient density result in larger dense sub-networks, respectively. We used the CCLE pharmacologic profiles to identify cell lines responsive and non-responsive to Taxol (Paclitaxel, Docetaxel,)8, a drug routinely used in NSCLC large and squamous cell combinatorial-chemotherapy for its promotion of microtubule assembly while inhibiting disassembly9. Taxol is no longer used as single-agent chemotherapy, in part due to an overall response rate of 22% (Figure 1.1). However, we aim to identify the subpopulation that does productively respond, for both their benefit of treatment as well as preventing an unresponsive population’s unproductive treatment. We calculated the eigengene of the resultant sub-network module from within which we selected responsive cell lines to be the lowest quartile of calculated IC50, and likewise the non-responsive strata as the highest IC50 quartile from 66 large and squamous cell NSCLC based upon CCLE calculated pharmacologic profiles. We compared expression of these eigengenes, here the trice-transformed representation of the module’s gene expression, resulting in three sub-networks of significance level alpha= 0.1. Further study could show these to be biomarkers of positive, productive response to Taxol chemotherapy.
 
Table 1.1 Phase II Trials of Single-Agent Taxol10. Taken from Socinski 1999. Low response rate indicates opportunity for screening measure to mitigate side effects of low response rate administrations. 
Chapter 2: Workflow for Gene Co-expression Network Analysis on CCLE data

2.1 Weighted Gene Co-expression Network Analysis (WGCNA)
	Co-expression network analysis upon microarray datasets has enabled the identification of strongly connected sub-regions of interaction, subgraphs of biological relevance due to their hyper connected nature which are novel due to the evaluation of initial node pair linkage among other hyper connected graphs. Notably Horvath et al. have used machine-learning methods of clustering to draw conclusion from these WGCN. Here we propose a more organic method of mining the data matrix. We refer to the mined gene groups as sub-network modules for means of visualization ease, as their relationships are multiple and easily manipulated in such a medium. The genes of any biological pathway are represented by nodes of an undirected graph, and connected by edge lengths, absolute correlation coefficients transformed through adjacency function11,12.
〖PCC〗_i=(∑_(j=1)^(n_i)▒〖(x_ij- μ ̂_X1)(y_ij-μ ̂_Y1)〗)/(√(∑_(j=1)^(n_i)▒〖(x_ij-μ ̂_X1)〗^2 )*√(∑_(j=1)^(n_i)▒〖(y_ij-μ ̂_X1)〗^2 ))
Equation 2.1 The PCC For any given pair of genes (X,Y) in the group is equal to the covariance of the gene set over the standard deviation of the composing set, a standard metric when measuring gene co-expression.
The CCLE provides 18,898 GC-RMA normalized genes for which expression is measured in patient non-small cell lung cancer samples grown up to culture. GC-RMA normalization takes into account an estimate of background cross-hybridization for each probe individually to counteract the fact that GC-rich probes seem to have higher non-specific signal (mismatch) among the 25 base Affymetrix oligomers13. This normalization provided for us, our first step in the analysis was to ensure that each correlation calculation was representing what we intended – similar expression levels among representative cell lines. For this we needed to control for low variance genes that would produce large PCCs but not because of high correlation (Equation 2.1). For this we selected the genes of highest average variance, discarding the bottom 15%. From here we were able to calculate the ~260,000,000 PCC for the correlation of 16,140 genes to themselves. The Pearson correlation coefficient metric (Equation 2.1) for proximity was decided upon over the Spearman rank because of the prior’s ability to identify highly similar expression profiles. We are aware of the Pearson correlation coefficient’s shortcomings (Table 2.1), including being too stringent leading to outlier bias, and would have adapted to use the Spearman metric if too many outliers were detected. The Pearson remained optimal, however, due to its assumptions of a Gaussian distribution as well as its history of use in microarray data analysis for resultant tight linear relationships 14, rather than the polynomial nature of RNA-seq data analyzed largely by the Spearman rank correlation coefficient.





Pearson’s Correlation Coefficient	Spearman’s Rank Correlation Coefficient
Gaussian	Non-Gaussian
Scale invariant re: linear transformation
	Tight relationships
	Highly similar	Non-linear
Proven in microarray analysis	RNA-Seq
Is sensitive to outlier observations	Is a Pearson correlation of their ranks

Table 2.1 Comparison of Correlation metrics.
The Pearson Correlation Coefficient measures the similarity of gene expression variation (variance of the pair / SD of each gene) as you progress through samples of a population. Our AF takes the absolute value of the raw PCC (co-expression similarity values) forming an adjacency matrix, which is subsequently subtracted from 1 (ie 1-|PCC|), creating an edge distance matrix (Figure 2.1).
 



Adjacency matrix
■( &■(  1  &2  &■(3&  ■(4& 5)))@■(1@2@■(3@■(4@5)))&[■(■(0&.6@.6&0)&■(.4&.3&.6@.7&.3&.1)@■(.4&.7@.3&.3@.6&.1)&■(0&.5&.4@.5& 0&.6@.4&.6& 0))] )





      

	
 


	

Figure 2.1 Adjacency and Distance matrices, resultant weighted graph. Every gene pair’s correlation is stored in a matrix, which is transformed to a matrix of connection lengths; the adjacency function (AF) is used to convert “a matrix values function that maps an n x n dimensional adjacency matrix onto a new n x n dimensional network adjacency”15.
2.2 WGCN Clustering vs Mining
	It is understood that a group of genes that exhibits a similar expression profile, expression pattern across cell lines, is likely to have similar function12. Grouped co-expression genes contained in matrixes can be represented as graphs; each node corresponds to a gene pair, connected by an edge of non-arbitrary length determined by the corresponding AF-transformed gene-pair correlation value. Groups of these hyper-connected co-expressed genes are termed cliques16 (Figure 2.2A), where every gene pair is joined. A clique is defined as the subgraph H of the un-weighted graph G=(V,E,W), where every pair of vertices is joined by one edge. We are interested in a less stringent quasi-clique (Figure 2.2B) of density only slightly deviated from maximum (Equation 2.2), which deviates from and relaxes this strict criterion of absolute hyper-connectivity by three threshold parameters. Briefly, they are γ, to evaluate neighboring edge weight so to select only the highest co-expressed gene-pairs, λ and t (ie K), to ensure added genes result in dense sub-networks, and β, to restrict merging to that resulting in larger dense sub-networks. A quasi-clique is then the expanded graph C where the connectivity is not absolute but relaxed and allows for varying degrees of (AF transformed) correlation, so capturing more gene-pair relationships; while every node is connected to multiple others, it is not necessary to connect to all other nodes, as is required by the definition of clique.  Unlike clustering, this method allows for gene overlap among resultant sub-network modules, as well limiting gene representation in modules to only those of highest co-expression.
d(C)=(∑_(i=1)^(N-1)▒∑_(j=i+1)^N▒w_ij )/(N(N-1)/2)
where vectorizeMatrix(w)=〖(w〗_12 〖,w〗_13,…,w_(n-1,n))
and n=(n_1,…,n_n )denotes the vector of node connectivities
Equation 2.2 Density of Subnetwork. Used in second threshold step comparing contribution of new gene to cluster’s density to initial density.
  
Figure 2.2 A set of vertices S is called a clique if there is an edge between any two vertices in the subgraph G(S). Our interest in dense subgraphs allows us to generalize this definition and look at less rigid quasi-clique Cg - a subset of vertices V with density function ~ sum(edge weight/ no.nodes) > threshold

Steve Horvath’s group has developed several hierarchical clustering based approaches for identifying such highly correlated gene sub-networks within WGCN (Figure 2.3). Although it is an improvement upon un-weighted graphs, the hierarchical clustering method often identifies large clusters that cannot directly dictate the intra-cluster connectivity17, or real biological systems for that matter. Instead there needs to be a limit of relatedness which dictates when clustering ends; our approach offers such a threshold. Other drawbacks of a clustering approach include (i) not allowing intracluster connectivity-direct control resulting in high average within-cluster vertices correlation; (ii) the inability to share genes between two sub-networks, a phenomena witnessed in natural systems; (iii) clusters identified often exceed 100 genes, overlooking smaller gene networks containing subtle functional information18. “We remove the WGCN hierarchical construction that does not contribute to our dense-sub-network finding.”18 We identify relatively co-expressed, ie highly connected gene networks (C), termed quasi-cliques, through a modified Quasi-Clique Merger algorithm18,19, the localMax-eQCM (section 2.3). In short, a set of comparative density functions’ threshold parameters mine gene pairs based on the matrix ranking step, the addition of genes to maximally correlated pairs and the subsequent merger of overlapping sub-networks, growing the sub-network modules in an organic manner, gene by gene, to a limit, merging sub-networks limited by their overlap.
 
Figure 2.3 Evolution of sub-network creation, increasing by physiology representation; The Horvath et al. approach do not incorporate a threshold in the formation of their dendrograms, so highly dissimilar objects merge the higher the graph progresses. This is remedied in QCM, and enables overlapping sub-networks to be merged as long as density criterion are met.
 

2.3 Quasi-Clique Merger Algorithm
	We hypothesize that datasets of NSCLC gene expression contain sub-networks of genes are in fact likely pathways20 whose expression patterns are highly correlative, which we can investigate through a data mining approach. We implemented a modified version of Ou and Zhang’s Quasi-Clique Merger algorithm more in keeping with the eQCM18 (Figure 2.4).  We chose this non-exact method to solve the NP-complete problem2119,21–23 of evaluating ~260,000,000 gene relationships because it guarantees mathematically good enough gene sub-network modules19. The mathematic proof of Ou and Zhang is enacted by the threshold parameters, which are necessary at three points along the localMax-eQCM implementation. The algorithm starts by exploring each gene’s 16139 PCC values, only stopping when the highest value matches with another gene’s highest correlation value, and only if it is above the γ threshold parameter. In this way, the addition of any other highly correlated gene will decrease the overall correlation, referred to here as sub-network density, by definition that the first gene pair is of highest correlation value. The second threshold parameter then determines subsequent node/vertices addition, insisting that the new density or overall correlation is high enough relative to the initial density relative to the an threshold parameter (Equation 2.3, 2.4 and 2.5). Furthermore, the quasi-cliques are unique because the merger step reduces overlapping gene sub-network above a third density threshold, β (Equation 2.6).
 
Figure 2.4 Evolution of QCM. The previous implementation of eQCM (grey) would form quasi-cliques based upon edges of lesser-correlated nodes. Our localMax-eQCM (red) is strict in that in only selects the highest correlation pair among neighboring nodes, reducing overlap among sub-networks to be merged.

The QCM method necessitates lowered densities bounds for nodes added after the initial pair because the first step selects the highest correlated pair available among the two genes correlation signatures (figure 2.5). This step theoretically begets the formation of more biologically meaningful sub-networks compared to hierarchical clustering due to the allowance of genes to be present multiple modules, which can be verified during enrichment. Following from the definition of the quasi-clique (section 2.3), each gene within the sub-network is correlated to every other gene to a high degree, dictated by the various user-defined thresholds through the steps of the algorithm. Our localMax-eQCM is able to maintain all weighted edges in order to return high-density sub-network modules, more reminiscent of biological networks than what a classical clustering algorithm could produce.
Every QCM derivation decides whether to include a gene to a defined gene sub-network based on the degree to which its expression correlates to other genes, rooted in the hypothesis co-expressed genes function together. To be consistent with the original QCM algorithm, the contribution of any vertex to the cluster (v(C)) is defined as the ratio of G(C) edge weight increases by adding the vertex v,18 over the size of C, and ultimately, a vertex v is added to C if their density remains above a user specified threshold, ie
c(v,C)  ≥a_n d(C)
Equation 2.3 Criteria for addition
(∑_(u∈V(C))▒〖w(uv)〗)/n  ≥[1-1/2λ(n+t) ]*(∑_(i=1)^(N-1)▒∑_(j=i+1)^N▒w_ij )/(N(N-1)/2)
Equation 2.4 Expanded criteria for addition
K=(f(n+1))/(f(μ))=(2(n+2))/(3(n+1))
Equation 2.5 an is bound by K.
where n=|V(C_p)|. The constant K is invoked to bound an to a function of its component nodes (Equation 2.5). This is mathematically proven to produce clusters of density sufficient to guarantee biological similarity within the group (see Ou and Zhang 2007).

After calculating the matrix of Pearson correlation coefficients for each gene to every other gene, we proceed with the heart of the localMax-eQCM, as follows:

 

Algorithm
	Select probes of top 85% variance, calculate PCC for every gene pair
	SORT the set of |PCC| edges
	If the highest correlation of one gene (x) is same cell as another (y), and
	If the edge is greater than the threshold parameter γ, and
	If the vertices are not already in a sub-network C, form new sub-network, where density is a function of their edge weight divided by their constituent vertices
	ADD another vertex if the contribution of the vertex to the sub-network only decreases the density by less than the alphaN factor α_n=1-  1/2λ(n+t)  
	Else the sub-network is at minimum density and therefore addition of any other vertex would return sub-network of insufficient density, so therefore the sub-network is complete; move to new starting edge (2). 
	SORT all sub-networks by number of vertices in a descending order
	If the sub-network has at least 3 vertices, and
	If the absolute density of the intersection of two sub-networks decreases the absolute density of each separately by a small enough margin (beta)
	MERGE the two sub-networks and remove the constituent sub-networks from the sorted list
	Run two-sample, unpaired t-test between the eigengene of each module related to response and non-response strata populations.
	Enrich significantly differentially expressed sub-networks.
The small sub-networks can be merged (Figure 2.6) if their shared genes are of such density above the third and final threshold, β in Equation 2.5, after which point the component sub-networks are removed from the module list and the new sub-network, a quasi-clique, is added. If there is no possible merger, sub-networks on this list are themselves quasi-cliques. Depending on the correlation between eigengene amongst response sub-populations, these modules can be selected as significant sub-networks. Once similarity is assessed in this manner, we set to answer questions pertinent to our hypothesis -- do NSCLC patients responsive to chemotherapy present some gene sub-network absent in non-responsive patients? Could these sub-networks be used as a screening biomarker treatment determination a priori?
              |C_j∩C_h |>βmin⁡(|C_j |,|C_h |)
Equation 2.6 Density must remain above threshold parameter β when forming quasi-cliques from smaller sub-networks of significant gene overlap. This utilizes the weight of the edges, creates similar, highly dense subnetworks. The result of seeing them together increases biological relevance.  
Figure 2.5 Merging sub-network matrixes
2.4 Review and Differentiation
Our localMax-eQCM algorithm, based on the data-mining pipeline developed by Ou and Zhang, enabled us to investigate genes of highly correlated expression across the NSCLC cell line population. These resultant groups are more biologically similar than methods used elsewhere by process of the initial exploratory step to find the maximal local gene correlation. A gene-pair’s distance metric, one minus the absolute value of the Pearson coefficient, is taken as the input for the algorithm, are filtered according to the first threshold parameters γ. This threshold defines local edge weight taken into consideration for local maximum as well as stepwise additions to said maximum pair. Once a maximum pair is located, a gene with high correlation to both is added so long as the density of the resultant sub-network remained above the an threshold, bound by K where t and λ here = 1. This proceeds until the addition trips the threshold, and on again until initial gene pairs exhausted. Finally, if sub-networks share genes at high enough density, the last step merges sub-networks so the new sub-network results in some overall density only some threshold β less dense than the contributing sub-network modules. Horvath et al. suggest that “Eigengene make more sense for network module”24 representation, and so were single value decomposition (SVD) was carried out per module, the third transformation of which produces an eigengene per cell line (66 x 66 NSCLC cell lines). After eigengene calculation representative of expression among quasi-clique modules among 66 large and squamous patients’ immortalized cell lines NSCLC, pharmacologic Taxol profiles were used to stratify the 66 cell lines, based upon top and bottom response quartile. We differentiated among the sub-networks’ response population eigengene by two-sample, unpaired t-test. We found three sub-networks (Table 3.1) to be significantly differentially expressed at an alpha level of 0.05 between the responsive strata (lowest 11 IC50 LC cell lines) and non-responsive strata (highest 11 IC50 LC cell lines).
 
Chapter 3: Differentially Expressed Quasi-Cliques
3.1 Quasi-Clique Formation & Merging
	We utilized the complete gene expression and pharmacologic profiles of the Cancer Cell Line Encyclopedia’s non-small cell lung cancer dataset for our differential expression approach to identify pathway variance among conditions responsible for differential response to chemotherapeutics. Translation to clinical care could involve screening for the biomarker gene sets found to be responsible for response to treatment. This could also aid in prediction of drug affectivity in patient subpopulation, matching individual to response subtype a priori.
	We present a quantitative experimental research design that stratifies condition groups by pharmacologic profile after quasi-clique formation based upon overall CCLE large and squamous cell NSCLC gene expression data (Figure 2.4). The analysis pipeline encompasses ranking the normalized dataset, establishing co-expression network via PCC calculation and transformation, assembling sub-network modules of highly correlative constitution via localMax-eQCM, and finally differentiating among Taxol response sub-population. We use a two-sample, unpaired t-test to identify quasi-cliques assembled from 66 large and squamous cell NSCLC expression, and looked for differentiated expression among these modules by selecting a subset of their representative eigengene; eigengene values were stratified among top and bottom pharmacologic response strata quartiles (11 each) after each module was calculated over the entire 66 cell lines. 
We can identify key pathways of drug response by comparison to non-response, isolating groups of gene created among the whole but investigated in a subset of extremes. We utilized CCLE-assembled large and squamous NSCLC gene signatures to identify cell line response to Taxol (Paclitaxel). Among the 24 drug profiles available, Taxol had the most complete testing among the NSCLC cell lines, as well as a wide enough variance (100 fold) among these cell lines to illicit statistically meaningful differentiation by which to define response/ non-response strata. Taxol response was characterized by fitting its multi-dosages kill curve across the cell lines, measured as the concentration at which the drug response reached an absolute inhibition of 50% of the culture (IC50)3. We selected responsive cell lines from the lowest calculated IC50 quartile, and likewise the non-responsive strata of the highest IC50 quartile. We used these strata to perform a two-sample, unpaired t-test for statistically significant differentiation among the quasi-cliques, represented as eigenegenes, resulting in sub-network modules found to be differentially expressed by large and squamous cell NSCLC lines responsive to Taxol treatment.

3.2 Quasi-Clique Differential Expression
	A patient population of adequate sample size to power a statistically significant result was found in the NSCLC patient cell lines assembled in the CCLE. To summarize, briefly, the first step in identifying gene co-expression networks responsible for response to NSCLC therapy was identifying the correct patient population; contained within the CCLE, 66 patient cell lines identified as large and squamous non-small cell lung carcinoma, 14 large cell carcinoma, 28 squamous cell carcinoma, and 24 general (Figure 1.3). We began by reading in and formatting the normalized gene expression data into a matrix within MATLAB. We indexed the matrix by gene ID and probe ID, and then ranked them by average expression and expression variation, restricting to genes of top 85% variance. Matrix manipulations are based upon this ranked expression matrix of 66 NSCLC cell lines by some 16,140 probes. A correlation matrix (cMatrix) of concordance (network distance) among the cell lines’ expression profiles was calculated using the covariance of any two genes’ expression divided by the product of their expression’s standard deviation (correlation). Our co-expression measure utilized the Pearson correlation coefficient that has become standard practice in microarray WCGNA25,26. We treated all 66 large and squamous cell NSCLC cell lines as a single population and calculated Pearson-correlation coefficients among gene pairs, which we arranged into a matrix of 16,140 x 16,140 gene-pair correlation values, some 260,499,600 correlation calculations. 
	The localMax-eQCM algorithm is tasked with assessing the averaged gene expression correlation values (section 2.3), to find the highest correlation per gene and determine if additional genes result in sub-networks of sufficient density. We next attempt to merge the resultant small sub-networks of sufficient overlap, returning biologically similar, dense gene sub-networks characteristic of NSCLC gene expression when administered Taxol. We begin here by checking the edge coverage by finding the local maximal edge and only iteratively evaluating those with correlation values larger than the weight threshold parameter γ (γ =0.7, t=1, λ=1, β=0.4) initially for 66 cell lines. Edges are added or discarded, until the entire correlation matrix is assessed, at which point sorting of the criteria-meeting edges begins and is iterated for each vertex assembled from within the matrix18. We represent the modules’ gene expression as eigengenes and performed a two-sample, unpaired t-test of genes composing the quasi-cliques we have formed amongst response strata of dissimilar pharmacologic profiles, ie low IC50 (responsive) and high IC50 (non-responsive). We concluded with an enrichment analysis among several online databases for known associations, including ToppGene27, GeneMania and FunCoup. 

1.86E-08		3.25E-17	1.92E-10	6.02E-10	2.88E-07	6.78E-05
module12		module34	module126	mod158	mod172	mod254
NQO1	ZNF323	PGD	ZNF738	CXCL1	ZNF738	CSF3
TDP2	LOC344887	GCLM	ZNF708	IL8	ZNF708	IL1B
SRXN1	OSGIN1	AKR1C3	ZNF430	CXCL2	ZNF430	IL1A
PGD	NR0B1	AKR1C2	LOC100132815	CXCL3	ZNF100	
ME1	EYA4	AKR1C1	ZNF100			
AKR1C3	CYP4F11	EID3				
EPHX1	KIAA0319	AKR1B10				
AKR1C2	EGF	LOC344887				
AKR1C1	SLC6A6	OSGIN1				
SPP1	JAKMIP3	NR0B1				
AKR1B10	AR	CYP4F11				

Table 3.1. Significant sub-networks acting in response to Taxol
 
Figure 3.1. Differential eigengenes expression of significantly differentially expressed module 12. Along the x-axis, non-response strata are cell lines 1-11 and response are cell lines 12-22.

Enrichment for modules 12 and 34 return ToppGene co-expression to “human lung_not_cancer_Harvery06”, as well as PubMed reference to smokers at risk for lung cancer; GeneMania returns aldo-keto/ oxidoreductase activity. Module 126 is composed of zinc finger genes and thusly enriches for zinc and transition metal ion binding, as well as a generic KRAB transcription pathway identified upstream from the zinc finger region.

 

Chapter 4: Discussion and Future Work

	The quasi-cliques of response found to be significantly differentially expressed between the strata gives proof to the adage of each patient being a microcosm of variability, that no cancer responds exactly the same and that therefore no repurposed drug can be expected to cover 100% of any population’s variability. The gene sub-networks involved in the response to Taxol in NSCLC patients present a unique opportunity to study specific exploits of a proven chemotherapy and thus expand the knowledge of its role in NSCLC therapy. 
A similar clinical decision support workflow could make its way into clinical implementation through the current practice of capturing and classifying any given patient to the most detailed condition whereby they are screened for the specific response gene sub-network(s), circulating or otherwise, to determine likelihood of response. The cost-benefit dilemma surely exists inside an oncologist mind when deciding a course of treatment for any given patient based upon population statistics of the drug’s effectiveness.
One could continue the work in both a computation mode as well as a verification mode in the wet lab (Table 4.1). Computationally, it would be interesting to process and form quasi-cliques from another data set, perhaps breast cancer vs triple negative breast cancer, posing such questions as what scale can this work on?, or more tailored clinical questions. 
Computational	Verification
ANOTHER DATA SET
Breast cancer vs triple negative:
	What scale can this work on?
	Tailored clinical questions	NSCLC
	NSCLC non-adenocarcinoma patient treated with Taxol data to run against
	Invstigate targeted glycolysis enzyme and inhibitor
	Test if cluster 3 is biomarker for response in patients
Table 4.1 Continued Study
 
 
References
1. 	Wang L, McLeod HL, Weinshilboum RM. Genomics and drug response. N Engl J Med. 2011;364:1144–1153. doi:10.1056/NEJMra1010600.
2. 	Dziegiel P, Jeleń M, Muszczyńska B, et al. Role of metallothionein expression in non-small cell lung carcinomas. Rocz Akad Med Bialymst. 2004;49 Suppl 1:43–45.
3. 	Barretina J, Caponigro G, Stransky N, et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature. 2012;483:603–307. doi:10.1038/nature11003.
4. 	Haibe-Kains B, El-Hachem N, Birkbak NJ, et al. Inconsistency in large pharmacogenomic studies. Nature. 2013;504(7480):389–93. doi:10.1038/nature12831.
5. 	Colizza V, Flammini A, Maritan A, Vespignani A. Characterization and modeling of protein-protein interaction networks. Phys A Stat Mech its Appl. 2005;352:1–27. doi:10.1016/j.physa.2004.12.030.
6. 	Hu H, Yan X, Huang Y, Han J, Zhou XJ. Mining coherent dense subgraphs across massive biological networks for functional discovery. Bioinformatics. 2005;21 Suppl 1:i213–21. doi:10.1093/bioinformatics/bti1049.
7. 	Costa L da F, Rodrigues FA, Cristino AS. Complex networks: the key to systems biology. Genet Mol Biol. 2008;31:591–601. doi:10.1590/S1415-47572008000400001.
8. 	Socinski MA. Single-agent paclitaxel in the treatment of advanced non-small cell lung cancer. Oncologist. 1999;4:408–416.
9. 	Pazdur R, Kudelka AP, Kavanagh JJ, Cohen PR, Raber MN. The taxoids: Paclitaxel (Taxol) and docetaxel (Taxotere). Cancer Treat Rev. 1993;19:351–386. doi:10.1016/0305-7372(93)90010-O.
10. 	Socinski MA. Single-agent paclitaxel in the treatment of advanced non-small cell lung cancer. Oncologist. 1999;4:408–416.
11. 	Butte AJ, Tamayo P, Slonim D, Golub TR, Kohane IS. Discovering functional relationships between RNA expression and chemotherapeutic susceptibility using relevance networks. Proc Natl Acad Sci U S A. 2000;97:12182–12186. doi:10.1073/pnas.220392197.
12. 	D’haeseleer P, Liang S, Somogyi R. Genetic network inference: from co-expression clustering to reverse engineering. Bioinformatics. 2000;16:707–726. doi:10.1093/bioinformatics/16.8.707.
13. 	Irizarry RA, Wu Z, Jaffee HA. Comparison of Affymetrix GeneChip expression measures. Bioinformatics. 2006;22:789–794. doi:10.1093/bioinformatics/btk046.
14. 	Hauke J, Kossowski T. Comparison of Values of Pearson’s and Spearman's Correlation Coefficients on the Same Sets of Data. Quaest Geogr. 2011;30(2). doi:10.2478/v10117-011-0021-1.
15. 	Horvath S. Weighted Network Analysis: Applications in Genomics and Systems Biology. May 2011. Springer; 2011:77–78. Available at: http://medcontent.metapress.com/index/A65RM03P4874243N.pdf. Accessed May 21, 2014.
16. 	Lee HK, Hsu AK, Sajdak J, Qin J, Pavlidis P. Coexpression analysis of human genes across many microarray data sets. Genome Res. 2004;14:1085–1094. doi:10.1101/gr.1910904.
17. 	Mahanta P, Ahmed HA, Bhattacharyya DK, Kalita JK. An effective method for network module extraction from microarray data. BMC Bioinformatics. 2012;13:S4. doi:10.1186/1471-2105-13-S13-S4.
18. 	Xiang Y, Zhang C-Q, Huang K. Predicting glioblastoma prognosis networks using weighted gene co-expression network analysis on TCGA data. BMC Bioinformatics. 2012;13 Suppl 2:S12. doi:10.1186/1471-2105-13-S2-S12.
19. 	Industrial JOF, Optimization M. A NEW MULTIMEMBERSHIP CLUSTERING METHOD Yongbin Ou and Cun-Quan Zhang. 2007;3(4):619–624.
20. 	Jiang D, Pei J. Mining frequent cross-graph quasi-cliques. ACM Trans Knowl Discov Data. 2009;2:1–42. doi:10.1145/1460797.1460799.
21. 	Pattillo J, Veremyev A, Butenko S, Boginski V. On the maximum quasi-clique problem. Discret Appl Math. 2013;161(1-2):244–257. doi:10.1016/j.dam.2012.07.019.
22. 	Tsourakakis C, Bonchi F, Gionis A. Denser than the densest subgraph: extracting optimal quasi-cliques with quality guarantees. Proc 19th …. 2013.
23. 	Horvath S. Weighted Network Analysis: Applications in Genomics and Systems Biology. May 2011. Springer; 2011:77–78.
24. 	Langfelder P, Horvath S. Eigengene networks for studying the relationships between co-expression modules. BMC Syst Biol. 2007;1:54. doi:10.1186/1752-0509-1-54.
25. 	Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. Stat Appl Genet Mol Biol. 2005;4:Article17. doi:10.2202/1544-6115.1128.
26. 	Butte AJ, Tamayo P, Slonim D, Golub TR, Kohane IS. Discovering functional relationships between RNA expression and chemotherapeutic susceptibility using relevance networks. Proc Natl Acad Sci U S A. 2000;97:12182–12186. doi:10.1073/pnas.220392197.
27. 	Chen J, Bardes EE, Aronow BJ, Jegga AG. ToppGene Suite for gene list enrichment analysis and candidate gene prioritization. Nucleic Acids Res. 2009;37:W305–W311. doi:10.1093/nar/gkp427. 

