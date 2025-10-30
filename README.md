# A situated perspective on CONIC: Evidence of compensatory effects of conscientiousness and situational interest on the task level

Conscientiousness and interest each have a substantial impact on learning. The conscientiousness × interest compensation (CONIC) model conceptualizes their interplay as compensatory in predicting academic effort. Previous research has primarily applied the CONIC model to aggregated measures of effort, leaving a gap in understanding the nuanced mechanisms underlying the activation of interest in specific learning situations. To address this gap, we investigated whether the compensatory relationship between conscientiousness and interest does not only pertain to the level of a course or subject but does also exists when working on individual tasks. Specifically, we examined how students’ conscientiousness and situational interest interact in predicting task effort. Our sample consisted of 1,839 secondary school students in Germany (*M*<sub>age</sub> = 16.4, *SD* = 1.5, 42.7% female). Using latent moderated structural equation modeling, we observed positive main effects of conscientiousness and situational interest on task effort, alongside negative interaction effects between these variables. The findings support the compensatory effects proposed by the CONIC model at the task level. This research contributes to a deeper understanding of how conscientiousness and situational interest interact dynamically to influence academic effort, offering insights into how these effects persist over time to influence broader learning outcomes.

## Citation

Höft, L., Meyer, J., Bernholt, S. & Jansen, T. (2025). A situated perspective on CONIC: Evidence of compensatory effects of conscientiousness and situational interest on the task level. Contemporary Educational Psychology, 102375. <https://doi.org/10.1016/j.cedpsych.2025.102375>

## Overview

This repository contains the materials accompanying the article  
> **A situated perspective on CONIC: Evidence of compensatory effects of conscientiousness and situational interest on the task level**  
published in *Contemporary Educational Psychology* (2025).

It provides the Quarto source file used to generate an APA-compliant manuscript using the [**apaquarto**](https://github.com/crsh/apaquarto) template, as well as scripts for reproducing the statistical analyses performed in **Mplus** via the **MplusAutomation** R package.

## Repository Contents

-   `article.qmd` — Quarto source file for generating the APA-formatted manuscript.
-   `ref.bib` — Bibliographic reference file for Quarto rendering.
-   `dat_mplus.csv` — Dataset used for the statistical analyses.
-   `codebook.csv` — Codebook describing all variables and their labels.
-   `img/` — Directory containing figures and visual materials.

## Reproduction Instructions

### 1. Requirements

-   [R](https://www.r-project.org/)
-   [Mplus](https://www.statmodel.com/)
-   [quarto](https://quarto.org/)
-   [apaquarto](https://github.com/wjschne/apaquarto)

### 2. Run the Analysis and Render the APA-compliant manuscript

``` bash
quarto render article.qmd
```

This will produce a `.docx` document formatted according to APA guidelines.

