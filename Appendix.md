# Appendix: Analysis Methods of Geometric Morphometrics

---

This Appendix describes the process of exploring appropriate metrics for analyzing developmental shape changes and their variation across populations. The exploration consists of the following steps:

**A1**.  Analysis methods for testing overall hypotheses regarding 'Treatment' and 'Population' factors based on generalized Procrustes coordinate shape data.

**A2**.  Potential limitations of analysis due to the high dimensionality of the generalized Procrustes coordinate shape space.

**A3-A5**. Exploration and determination of appropriate metrics.

**A6**.  Visualization methods for shape variation using appropriate metrics.



## A1. Shape analysis using Procrustes ANOVA with RRPP

RRPP (Residual Randomization in Permutation Procedures) offers substantial advantages for dealing with high-dimensional geometric morphometric data (Adams and Collyer, 2016; Collyer and Adamus, 2024; Collyer et al., 2015). We first analyzed generalized Procrustes coordinate shape data using RRPP for an overall hypothesis model incorporating 'Treatment type' as a fixed effect and 'Population (Location)' as a random effect. The results of the analysis using this code are presented in Table 1 of the main text. The R code for the analysis is shown below. 

[source data on GipHub]('./Symm.BE.GPA.AllInd.txt')  (データファイルをGitHubにアップして、そのリンクにするとアップしたGitHubが開くようになるらしい)

```
# Data loading
setwd('/Users/.../...')																# Path to the location of the data file
data=read.delim("Symm.BE.GPA.AllInd.txt", head=T)			# data file

# Data preparation
shape_vars <- grep("Symm", names(data), value = TRUE)
shape_data <- as.matrix(data[, shape_vars])
```

```
# Required libraries
library(geomorph)
library(RRPP)
library(ggplot2)

# Model specifying a random effect in procD.lm from 'geomorph'
fit_proc <- procD.lm(shape_data ~ Treatments + Location + Treatments:Location,
                     data = data,
                     iter = 999,
                     RRPP = TRUE,
                     print.progress = TRUE,
                     effect.type = "F",
                     SS.type = "I",
                     random.effect = "Location")
# Basic ANOVA results
anova_result <- anova(fit_proc)                                     # ANOVA table

# The ANOVA table from anova(fit_proc) does not reflect the degrees of freedom lost due to GPA (a loss of df = 4) in the Residuals Df.
# The following code corrects the ANOVA results to account for this.

################################################################ Create the corrected ANOVA table
# Start
# Degrees of freedom lost due to GPA (4 for 2D data)
gpa_df_loss <- 4

# Get the original ANOVA table and convert it to a data frame
anova_df <- as.data.frame(anova_result$table)

# Correct the degrees of freedom
# Correct the Residuals Df
anova_df$Df[4] <- anova_df$Df[4] - gpa_df_loss  # 317 -> 313

# Also correct the Total Df
anova_df$Df[5] <- anova_df$Df[5] - gpa_df_loss  # 331 -> 327

# Recalculate Mean Squares (MS = SS/Df)
anova_df$MS[4] <- anova_df$SS[4] / anova_df$Df[4]

# Recalculate F-values (as needed)
# Divide the MS of each effect by the new Residual MS
for (i in 1:3) {
  anova_df$F[i] <- anova_df$MS[i] / anova_df$MS[4]
}

# Recalculate p-values (assuming Z-scores are not affected)
# Since p-values are based on permutations (iter), they cannot be simply recalculated from the F-distribution.
# Note: When the Df correction is small, the impact on p-values is expected to be minimal.

# Display the corrected ANOVA table
cat("Corrected ANOVA table (accounting for Df loss from GPA):\n")
print(anova_df)

# For comparison, also display the original ANOVA table
cat("\nOriginal ANOVA table (not accounting for Df loss from GPA):\n")
print(anova_result$table)                                         # <--- ANOVA Table for GPA data on a hypersphere
```

Although some researchers have advocated for analytical methods such as pairwise and trajectory analysis implemented in RRPP using the GEOMORPH functions "pairwise" and "trajectory", which compare distances and angles among multiple objects in high-dimensional generalized Procrustes coordinate space following an overall hypothesis test (Collyer and Adams 2013; Collyer et al. 2015), we recognized that this approach was not suitable for our intended analysis of specific hypotheses and involved potential limitations. Consequently, we explored methods to extract information appropriate for such specific hypothesis analyses from the high-dimensional generalized Procrustes coordinate space.




> **References:**
>
> Adams, D.C. & Collyer, M.L. (2016). On the comparison of the strength of morphological integration across morphometric datasets. *Evolution*, 70, 2623-2631.
>
> Collyer, M. L., and D. C. Adams. (2013). Phenotypic trajectory analysis: comparison of shape change patterns in evolution and ecology. *Hystrix-Italian Journal of Mammalogy* 24:75-83.
>
> Collyer, M.L., Sekora, D.J. & Adams, D.C. (2015). A method for the analysis of phenotypic change for phenotypes described by high-dimensional data. *Heredity*, 115, 357-365.
>
> Collyer, M. L., and D. C. Adams. (2024) RRPP: Linear Model Evaluation with Randomized Residuals in a Permutation Procedure. R package version 2.1.0. 

---



## A2. General difficulties in analyzing different trajectories in high dimensions

### Landmark data on the Shape space

The generalized Procrustes coordinates shape data, recorded in the number of $K$ landmarks in $M$-dimensional coordinates, is defined on a unit hypersphere with $MK-(M+1)-M(M-1)/2$ degrees of freedom in $MK$-dimensional space. The shapes of S5 (pre-feeding Start 5-day), W14 (Worm-feeding 14-day), and T14 (Tad-feeding 14-day) are mapped onto the unit hypersphere (shape space).  The developing shape change caused by each of the two different diets, worm or tadpole, is characterized in this shape space by the directions of the vectors S5 → W14 and S5 → T14 and, in particular, their cross angle, and by the lengths of the vectors S5 → W14 and S5 → T14. The location of the innate shape, S5, also characterizes the shape development reaction norm.

#### *- Induced developmental shapes on the shape space: Vector sets representation*

In high-dimensional shape spaces, comparing and evaluating the characteristics of multiple vector sets involves substantial geometric complexity. We will explain how this causes difficulties in comparing different sets of induced developmental shape trajectories. For example, assessing whether the relationship between vector pair S5_1→W14_1 and S5_1→T14_1 shares the same intrinsic properties as the relationship between vector pair S5_2→W14_2 and S5_2→T14_2 is difficult to determine. This difficulty arises because the subspaces (planes) in which each vector set is embedded do not necessarily coincide with one another. 

We can understand that, from the simplified example in a three-dimensional space, which is actually smaller than the minimal model case for the landmark method (Figure A2).



<div style="display: flex; justify-content: center;">
<img src="./images/dazzle.png" alt="dazzle" style="zoom:40%;" />
</div>
<div style="margin-left: 180px; margin-right: 180px; font-size: 14px;">
<b>Fig. A2.</b> Schematic diagram showing shape development in an imaginary quasi-shape space due to different diets in two populations (red and blue). The shape changes due to different diets are represented as vectors originating from vertex S5 of each triangle to W14 and T14, respectively. Population 1 consists of S5_1 = (2, 3, 0.8), W14_1 = (-0.01, -0.01, 1), and T14_1 = (4, 0, 1.21), while population 2 consists of S5_2 = (6.5, 4, 2), W14_2 = (6, 1.01, 0.01), and T14_2 = (6.19, 1, 4.01) in the shape space coordinates. Geometrically, the two triangles △(S5_1, W14_1, T14_1) and △(S5_2, W14_2, T14_2) are congruent. However, the visual representation of these geometrically identical triangles changes substantially depending on the viewing perspective. Panel (A) shows these triangles in nearly parallel-translated positions, appearing congruent and similarly oriented. In contrast, panel (B) from a different viewing angle reveals that the vectors W14 and T14 emanating from S5 appear mirror-reversed, and the triangles no longer appear congruent. This apparent discrepancy arises because the two triangles are embedded in different planes (subspaces) within the 3D space. The dihedral angle between these two planes is 83.55°, which explains why the triangular configurations appear fundamentally different depending on the viewing perspective.
</div>



Even in a schematic quasi-shape space (Figure A2), where the differences or similarities between two vector sets are examined, substantial geometric complexity persists. As this example illustrates, in actual high-dimensional shape space, when the subspaces embedding different sets of vectors are ambiguous in their similarity or distinctness, determining whether these vector sets have geometrically similar or dissimilar properties, and interpreting their biological significance, becomes problematic.

This geometric complexity reveals a fundamental challenge: assessing morphological congruence or divergence in high-dimensional generalized Procrustes space requires integration of complex geometric information. Within individual subspaces, the apparent similarity or difference between patterns depends critically on the geometric structure of the embedding subspace. Consequently, researchers are unable to derive clear biological conclusions from multidimensional shape data unless they can explicitly consider the structural relationship between the morphology itself and the Procrustes space used to map it.

---



## A3. Exploring biologically informative subspace

To overcome the difficulty in identifying geometric differences or similarities between two or more vector sets using information in all spatial directions of the high-dimensional shape space, it is required to discover a shared subspace of all vector sets that contains summary geometric information, facilitating a comparative analysis.

### **Shape-space dimensionality reduction approaches: PCA vs. bgPCA**

We will examine two methods, PCA and bgPCA, to discover an appropriate shared subspace, given that we are analyzing the shape development trajectories S5 → W14 and S5 → T14 from multiple observations at S5, W14, and T14 across populations. 

**P**rincipal **C**omponent **A**nalysis (**PCA**) extracts subspaces by maximizing the total variance across all observed samples, sequentially identifying orthogonal axes that capture the largest remaining variance. By construction, the resulting principal components are orthogonal and uncorrelated. However, this variance-maximization approach does not prioritize biological or experimental distinctions. The extracted subspace reflects the dominant patterns of morphological variation present in the dataset, regardless of group structure; these patterns may or may not correspond to biologically meaningful group differences or 'Treatment' effects. 

In contrast, **b**etween-**g**roup **PCA** (**bgPCA**) constructs subspaces specifically from group mean coordinates, focusing exclusively on inter-group morphological differences. By prioritizing mean-level separation between groups, bgPCA enhances the detection of group-specific morphological divergence and patterns that may be biologically meaningful. However, a critical limitation emerges when individual sample points are projected onto bgPC axes: the resulting bgPC scores retain non-zero covariance structures. Consequently, features captured by different bgPC axes cannot be independently interpreted, as they exhibit statistical interdependence. Moreover, bgPC axes do not necessarily maximize the total variance structure of the original shape data, potentially missing important within-group variation or overall geometric structure of morphological change.

For studies investigating 'Treatment' and 'Population' effects on organismal form, (1) bgPCA provides a more focused extraction of the relevant morphological subspace, (2) while traditional PCA offers complementary insights into overall shape variation structure. The mathematical properties reflected in these features can be evaluated as (1) the disadvantage of non-orthogonality of the data distribution in bgPC coordinates and (2) the advantage in terms of the degree of separation between groups in bgPCA compared to PCA. We diagnose these points. 



### PCA vs. bgPCA comprehensive comparison

#### *- Assessment of Orthogonality in bgPCA* 

In bgPCA, the bgPC scores for all samples are not necessarily distributed orthogonally to the bgPC axis, and the off-diagonal elements of the covariance matrix are not degenerate to zero. We assess the degree of covariance generation and variance reduction due to deviations from the orthogonal distribution along the bgPC axis for our data. 

Lets ${\bold C} =(c_{ij})$ be the covariance matrix for the bgPC scores of all samples.  Here,  for $j=i, c_{ii}={\mathrm{Var}}({\rm{bgPC}}_i)$, $c_{ij}={\mathrm{Cov}}({\rm{bgPC}}_i,{\rm{bgPC}}_j)$. Let the contribution vector $\mathbf{w} = (w_1, w_2, \ldots)^T $  be 
$$
w_i = \frac{c_{ii}}{\text{tr}(\mathbf{C})} = \frac{c_{ii}}{\sum_{j=1}^{} c_{jj}}.
$$
We evaluate non-orthogonality by two measures, the **Variance retention ratio** ($R$) and ***W*eighted *O*rthogonality *R*etention** ($WOR$).

● **Variance retention ratio** ($R$) is defined as, 
$$
R = \frac{\sum_{i=1}^{} c_{ii}}{\sigma^2_{\text{total}}}.
$$
The numerator represents the sum of variances of all bgPC components, which equals the trace of the covariance matrix, $\text{tr}(\mathbf{C})$. The denominator represents the total variance in the original dataset before any transformation. $R$ means the proportion of variance in the bgPC scores relative to the total variance in the original data. Values approaching 1.0 indicate minimal information loss, while values significantly below 0.95 typically suggest substantial loss of variance that may compromise the validity of subsequent analyses.  In the context of bgPCA evaluation, this high retention ratio demonstrates that despite the non-orthogonal nature of bgPCA, the method successfully preserves nearly all variance information from the original data, supporting its statistical viability as an alternative to traditional PCA. In our dataset, $R=0.9868$.

● ***W*eighted *O*rthogonality *R*etention** ($WOR$) is defined as, 
$$
WOR=1-\frac{\sum_{i=1}^{} \sum_{j \neq i} |c_{ij}| \sqrt{w_i w_j}}{\sum_{i=1}^{} c_{ii} w_i}.
$$
It measures how much orthogonality is preserved when component importance is weighted by variance contributions. While the variance retention ratio compares total variances (bgPC variance sum vs. original data variance), $WOR$ evaluates structural relationships between components. The numerator quantifies variance-weighted off-diagonal interactions, emphasizing relationships involving high-contribution components. The denominator provides variance-weighted normalization. The ratio itself represents weighted non-orthogonality intensity. Subtracting from 1 converts this into an intuitive retention measure: higher values indicate better orthogonality preservation. 

This contrasts sharply with the variance retention ratio, which focuses on total information preservation. $WOR$ specifically assesses whether components maintain independence when weighted by their analytical importance, making it sensitive to violations involving dominant components while de-emphasizing interactions among minor components. The metric answers "how orthogonal do the important components remain?" rather than "how much variance is preserved? In our case $WOR = 0.876$. The level represents a borderline acceptable loss of orthogonality. 



#### *- Separation Performance Comparison: PCA vs. bgPCA*

We compared the group separation performance of the 'Treatment' and 'Population' factors between bgPCA and traditional PCA using our dataset. 

The separation performance evaluation employs a one-way ANOVA framework to quantify how effectively each principal component axis discriminates between experimental groups. For each component axis and grouping factor ('Treatment' or 'Population'), the analysis calculates F-ratios by dividing the between-group mean square by the within-group mean square, where the between-group variance reflects how much group means differ from the overall mean, and the within-group variance captures variability within each group. Higher *F*-ratios indicate better separation, as they represent cases where groups are well-separated relative to their internal scatter. The comparison between PCA and bgPCA involves computing F-ratios for corresponding axes (PC1 vs bgPC1, PC2 vs bgPC2, etc.) and forming improvement ratios by dividing bgPCA F-ratios by PCA F-ratios. Individual axis improvements are weighted by each component's variance contribution, ensuring that separation gains in dominant components receive proportionally greater influence than those in minor components. This weighting scheme acknowledges that, for example, enhanced discrimination in a component explaining 60% of variance has greater practical relevance than improvements in a component explaining only 3% of variance, thereby evaluating bgPCA's performance relative to traditional PCA based on variance-weighted considerations.

<div style="margin-left: 10px; margin-right: 20px; font-size: 14px;">
<b>Table A3. </b> <i>F</i>-ratios represent separation effectiveness calculated via one-way ANOVA (between-group variance divided by within-group variance) for each principal component axis. Improvement ratios indicate bgPCA performance relative to PCA, where values >1.0 favor bgPCA. Variance weights show each component's contribution to total explained variance. <i>F</i>-statistic degrees of freedom are not reported because this analysis uses <i>F</i>-ratios as descriptive effect size measures rather than conducting formal hypothesis tests. The improvement ratios (bgPCA <i>F</i>-ratio / PCA <i>F</i>-ratio) lack established statistical distributions for parametric testing. No significance testing is performed on these ratios, only descriptive comparison of separation magnitudes across methods.
</div>


| Axis | PCA Treatment | PCA Population | bgPCA Treatment | bgPCA Population | Improvement Treatment | Improvement Population | bgPC Variance Weight |
|------|----------------|--------------|------------------|----------------|----------------------|---------------------|-----------------|
| 1    | 40.067         | 57.977       | 43.980           | 58.415         | 1.098                | 1.008               | 64.7%           |
| 2    | 57.911         | 17.248       | 53.100           | 27.930         | 0.917                | 1.619               | 10.3%           |
| 3    | 33.229         | 10.024       | 32.259           | 10.995         | 0.971                | 1.097               | 7.7%            |
| 4    | 17.676         | 8.785        | 10.812           | 3.327          | 0.612                | 0.379               | 6.8%            |
| 5    | 8.337          | 15.188       | 5.991            | 18.767         | 0.719                | 1.236               | 2.3%            |
| 6    | 4.615          | 6.809        | 11.490           | 9.953          | 2.489                | 1.462               | 2.0%           |



The analysis reveals inconsistent performance across individual axes, with bgPCA performing worse than PCA on several components (Axes 2, 3, 4, and 5 for Treatment; Axis 4 for Population). While Axis 4 shows substantial degradation for both factors (0.612× for Treatment, 0.379× for Population), its relatively modest variance contribution (6.8%) limits the practical impact of this decline.

To assess overall performance, individual axis improvements were weighted by their variance contributions and summed to create a composite measure. This variance-weighted assessment revealed that bgPCA provides only marginal advantages over PCA: 5.3% improvement for Treatment discrimination and 5.2% for Population discrimination. These modest gains suggest that, despite bgPCA's theoretical focus on between-group separation, its practical benefit over traditional PCA is limited for this dataset when considering the full component spectrum.　



#### *-Final evaluation on PCA vs. bgPCA*

***Orthogonality***:  Principal Component Analysis (PCA) maintains orthogonal component axes by mathematical construction. In contrast, between-group PCA (bgPCA) exhibits variable orthogonality depending on the metric used. When evaluated using variance-weighted orthogonality retention, bgPC scores show degraded orthogonality (87.6% retention), indicating that components are not fully independent. This degradation arises from bgPCA's focus on between-group separation rather than variance maximization, resulting in non-orthogonal geometric relationships among bgPC axes.

***Separation***: The separation effectiveness of bgPCA relative to PCA varied inconsistently across components. While bgPCA showed improvements on Axis 1 (1.098× for Treatment, 1.008× for Population), it exhibited substantial degradation on Axes 4 and 5, particularly for Treatment discrimination (0.612× and 0.719×, respectively). When all axes were considered together using variance-weighted averaging, bgPCA provided only marginal improvements: 5.3% for Treatment and 5.2% for Population. These minimal gains do not justify the adoption of bgPCA over traditional PCA.

**Conclusion:** Based on a comprehensive analysis of both orthogonality and separation performance, traditional PCA is more appropriate for this dataset than bgPCA. Although bgPCA theoretically targets between-group separation, its practical benefits are negligible for the current data. Moreover, bgPCA's reduced orthogonality compromises the interpretability of individual components, as their independence cannot be assumed. Therefore, we preferred PCA over bgPCA, as it provides mathematically rigorous, orthogonal component axes while achieving comparable group separation with greater simplicity and interpretability. 

---



## A4. Developmental shapes on the PC-subspace planes

We are exploring appropriate subspaces within the generalized Procrustes shape space to conduct comparative analyses of induced shape development—specifically, S5 to W14 and S5 to T14—across populations. To identify these subspaces and quantify shape changes, we employed Principal Component Analysis (PCA) as an initial exploratory approach. PCA was performed on the Generalized Procrustes shape dataset, and for each PC score, the mean scores by 'Population' and 'Treatment' were obtained.  

As can be inferred from the schematic example in Figure A2, it is difficult to understand the developmental shape changes induced by two different diets across different populations by examining the respective vector lengths and directional angles of multiple vectors, and the included angles between the S5→W14 and S5→T14 vectors within each population in the subspace derived from orthogonal projection of generalized Procrustes shape data onto principal component axes. In fact, for each of the projected planes PC1-PC2, PC1-PC3, and PC2-PC3, the lengths, azimuth angles, and intercept angles of each S5→W14 and S5→T14 vector, and the included angles between S5→W14 and S5→T14 vectors in each population were diverse (Figure A4).



<div style="display: flex; justify-content: center;">
  <img src="./images/image-20250819132653738.png" alt="image-20250819132653738" style="zoom:100%;" />
</div>
<div style="margin-left: 180px; margin-right: 120px; font-size: 14px;">
<b>Fig. A4.</b>  Scatter plots for each pair of PC1 - PC3, derived by performing a principal component analysis (PCA) for the Generalized Procrustes shape dataset. In each principal component pair, individual data points are represented by small dots and group means by larger dots. Arrows indicate developmental shape changes.
</div>



Figure A4 shows the magnitude and directional angles of developmental shape changes across multiple subspaces formed by the dominant principal components to identify characteristic patterns within each subspace. However, this analysis revealed that consistent, unified information could not be extracted across the multiple subspace representations. This demonstrates that the structure of shape variation is inherently multidimensional, and no single low-dimensional projection can capture the unified structure of shape variation.

---



## A5. Variance Component Analysis of Principal Components

To partition the morphological variation captured by each principal component axis into distinct factors, we conducted a variance component analysis using linear mixed models. For each principal component score, we fitted a mixed model incorporating 'Treatment' as a fixed effect, 'population (Location)' as a random effect, and their interaction as an additional random term. This model structure allows us to decompose the total variance into four orthogonal components: (1) Treatment variance, representing systematic morphological differences attributable to experimental feeding conditions; (2) Population variance, capturing population-level differences; (3) Interaction variance, reflecting population-specific responses to treatment; and (4) Residual variance, representing unexplained within-group variation.

Variance components were extracted from the mixed model using the estimated variance-covariance matrix. For each PC axis, we calculated the proportion of variance explained by each component relative to the total variance within that axis, yielding contribution ratios for 'Treatment', 'Population', 'Interaction', and 'Residual' effects. These ratios were then scaled by the variance explained by each principal component to estimate the contribution of each component to the total morphological variation across the entire dataset.

Two weighting schemes were applied to assess the relative importance of variance components. In Case 1, each PC axis's variance contributions were presented as percentages of the total dataset variance, directly reflecting the geometric significance of each component's effects. In Case 2, residual variance was excluded, and non-residual components ('Treatment', 'Population', 'Interaction') were rescaled to 100%, emphasizing systematic effects while minimizing the influence of noise. This analysis reveals how 'Treatment' effects, population differences, and their interactions are distributed across the PC spectrum, clarifying which axes carry the most biologically meaningful systematic variation.



<div style="margin-left: 0px; margin-right: 0px; font-size: 14px;">
<b>Table A5. </b> (A) presents the variance component decomposition for the first six principal components, which collectively account for 95% of the total shape variance in the dataset. The PC_total variance for each component is additively decomposed into four sources: 'Treatment' effects, 'Population' effects, 'Treatment×Population interaction' effects, and Residual variance. Values in parentheses represent the percentage contribution of each variance component within that specific PC. (B) presents the variance component decomposition based on non-residual variance, where systematic effects ('Treatment', 'Population', and 'Interaction') total 100%. This rescaling excludes residual variance to focus specifically on biologically interpretable variation patterns.

**(A)**: Total variance of all PCs = 100%

| PC   | PC_total     | Treatment      | Population     | Interaction    | Residual       |
| ---- | ------------ | -------------- | -------------- | -------------- | -------------- |
| PC1  | 64.59  (100) | 12.11  (18.75) | 23.36  (36.17) | 10.70  (16.57) | 18.41  (28.51) |
| PC2  | 10.43  (100) | 2.80  (26.79)  | 1.26  (12.11)  | 1.62  (15.51)  | 4.76  (45.58)  |
| PC3  | 7.92  (100)  | 1.41  (17.79)  | 0.43  (  5.44) | 0.90  (11.38)  | 5.19  (65.39)  |
| PC4  | 6.52  (100)  | 0.60  (  9.16) | 0.26  (  4.04) | 0.96  (14.70)  | 4.70  (72.10)  |
| PC5  | 2.50  (100)  | 0.13  (  5.16) | 0.33  (13.28)  | 0.22  (  8.06) | 1.83  (72.96)  |
| PC6  | 2.10  (100)  | 0.06  (  2.98) | 0.10  (  4.69) | 0.10  (  4.69) | 1.84  (87.63)  |

**(B)**: Total non-residual variance = 100%

| PC   | PC_total | Treatment | Population | Interaction |
| ---- | -------- | --------- | ---------- | ----------- |
| PC1  | 80.51    | 21.11     | 40.74      | 18.66       |
| PC2  | 9.90     | 4.87      | 2.20       | 2.82        |
| PC3  | 4.78     | 2.46      | 0.75       | 1.57        |
| PC4  | 3.17     | 1.04      | 0.46       | 1.67        |
| PC5  | 1.18     | 0.23      | 0.58       | 0.38        |
| PC6  | 0.45     | 0.11      | 0.17       | 0.17        |



<div style="display: flex; justify-content: center;">
  <img src="./images/image-20251101184424730.png" alt="image-20251101184424730" style="zoom:34%;" />
</div>
<div style="margin-left: 150px; margin-right: 160px; font-size: 14px;">
<b>Fig. A5.</b>  Variance component decomposition of principal components (PCs) for generalized Procrustes coordinate shape data. (A) Total variance decomposition showing the contribution of each PC to 100% of total shape variance, with each bar representing the proportional breakdown of variance components ('Treatment', 'Population', 'Interaction', and 'Residual') within each PC. (B) Non-residual variance decomposition where systematic effects ('Treatment', 'Population', and 'Interaction') are rescaled to 100%, excluding residual variance. 
</div>



The variance distribution among the principal components reveals a highly skewed distribution. PC1 dominates the morphometric variation, accounting for approximately 64.59% of total variance, while even the variance explained by the second principal component was a small value of about 10% (Table A5(A) and Figure A5(A)). This pattern indicates that shape variation is primarily captured along a single principal axis of variation. Furthermore, PC1 showed a relatively small proportion of residual variance (28.51%). PC1 shows a relatively balanced distribution of systematic effects, with ‘Population’ effect contributing the largest proportion (36.17%), followed by ‘Treatment’ (18.75%) and 'Interaction' (16.57%) effects (Table A5(A) and Figure A5(A)). The substantial systematic components in PC1 indicate that this primary morphometric axis captures biologically meaningful variation related to both genetic ('Population') and environmental ('Treatment') factors. In contrast, higher-order components (PC2-PC6) show progressively larger residual variance proportions, ranging from 45.58% in PC2 to 87.63% in PC6 (Table A5(A) and Figure A5(A)). This pattern suggests that while these secondary axes contribute to total morphometric variation, they primarily reflect individual-level noise rather than systematic experimental effects. 

Analysis of variance, excluding residual variance, reveals further characteristics between PCs and the systematic effects. PC1 dominates the systematic variance structure, accounting for 80.51% of all non-residual variation (Table A5(B) and Figure A5(B)). Within PC1, 'population' effects contribute the largest proportion (40.74%), followed by 'Treatment' effects (21.11%) and Interaction effects (18.66%)(Table A5(B) and Figure A5(B)). This distribution indicates that population-level differences represent the primary source of systematic shape variation, while treatment responses and population-specific treatment interactions also contribute substantially. This decomposition demonstrates that when individual-level noise is excluded, shape responses to experimental treatments and population differences are highly structured along a primary gradient of variation, with secondary axes contributing progressively less to systematic biological effects. This concentration indicates that the primary axis of shape variation captures the majority of biologically meaningful changes, while the remaining components contain little shape change information and gradually decrease in their contributions. 

This pattern suggests that while morphometric shape variation is multidimensional, the biologically meaningful responses to experimental treatment and population differences are concentrated along a primary axis of variation, with secondary axes contributing less to systematic effects. Based on these overall considerations, in order to understand and evaluate the differences in shape development caused by the two different types of food, only the PC1 score will be meaningfully employed, and information on the remaining PC axes will be excluded.

---



## A6. Thin-plate spline presentation

To interpret the biological significance of variation captured by principal component analysis, we visualized shape deformation patterns using thin-plate spline (TPS) deformation grids. This method reconstructs the projections of data samples onto principal component vectors in the generalized Procrustes shape space as landmark coordinates, facilitating the interpretation of biologically meaningful shape patterns.

We primarily focused on PC1, which captured the largest proportion of systematic morphological variation (80.5% of non-residual variance), and examined PC2 and PC3 as secondary components (9.9% and 4.8%, respectively)(Figure A6). For each principal component, we generated deformation grids at positive and negative extremes along each component axis. Deformation visualization scales were proportionally adjusted to each component's non-residual variance contribution: PC1 (scale ±0.1), PC2 (scale ±0.012), and PC3 (scale ±0.006). This proportional scaling ensures that visual emphasis reflects the relative magnitude of morphological variation captured by each component.



<div style="display: flex; justify-content: center;">
  <img src="./images/image-20250824121437904.png" alt="image-20250824121437904" style="zoom:60%;" />
</div>
<div style="margin-left: 150px; margin-right: 150px; font-size: 14px;">  
<b>Fig. A6.</b> Thin-plate spline deformation grids illustrating shape variation along principal components PC1-PC3 based on non-residual variance contributions. Each panel displays deformation patterns at positive (upper row) and negative (lower row) extremes of each component axis. Black dots represent the consensus (mean) landmark configuration, while red dots indicate deformed landmark positions. Blue vectors show the magnitude and direction of landmark displacement from the consensus shape. Grid distortions visualize the relative magnitude and spatial patterns of morphological variation along each principal component axis.</div>



These deformation grids reveal a clear hierarchical organization of morphological effects. PC1 (80.5% variance) dominates, exhibiting large-scale deformations affecting overall head shape. PC2 (9.9%) and PC3 (4.8%) show progressively smaller and more localized deformations. This pattern demonstrates that head shape development in response to feeding treatment and population differences is fundamentally governed by a single primary axis of variation, with secondary axes contributing only minor, refined modifications to morphological structure.













