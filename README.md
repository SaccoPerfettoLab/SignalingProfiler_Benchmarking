# **SignalingProfiler** 2.0: a network-based approach to bridge multi-omic data to phenotypic hallmarks

## Abstract

Unraveling the cellular signaling remodeling upon a perturbation is a fundamental challenge to understand disease mechanisms and to identify potential drug targets. In this pursuit, computational tools that generate mechanistic hypotheses from multi-omic data have an invaluable potential. Here we present SignalingProfiler 2.0, a multi-step pipeline to systematically derive context-specific signaling models by integrating proteogenomic data with prior knowledge-causal networks. This is a freely accessible and flexible tool that incorporates statistical, footprint-based and graph algorithms to accelerate the integration and interpretation of multi-omic data. Through benchmarking and rigorous parameter selection on a proof-of-concept study, performed in metformin treated breast cancer cells, we demonstrate the tool's ability to generate a hierarchical mechanistic network that recapitulates novel and known drug-perturbed signaling and phenotypic outcomes. In summary, SignalingProfiler 2.0 addresses the emergent need to derive biologically-relevant information from complex multi-omic data by extracting interpretable networks.

Read the full paper on [bioxRiv]().

## Repository

The repository showcases the code used for the benchmarking of *SignalingProfiler* in the paper. We present five key notebooks reporting the code and some plots of the analyses:

1.  [Data cleaning and preparation for *SignaligProfiler*](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/1.Data_preparation_for_SP.html)

2.  [Training dataset gold standard creation for the parameters tuning process (Sacco et al., 2016)](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/2.Gold_standard_creation.html)

3.  [Step 1 parameters tuning process](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/3.Step1_benchmark.html)

4.  [Step 2 parameters tuning process](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/4.Step2_benchmark.html)

5.  [Step 3 parameters tuning process](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/5.Step3_benchmark.html)

6.  [Visualization of best network of the parameters tuning process on training dataset](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/6.Visualization.html)

7.  [Validation of optimal parameters on two independent datasets](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/10.Benchmarking_validation.html)

8.  [EGF analysis](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/8.HeLa_EGF_SP_analysis.html)

9.  [AML analysis](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/9.JMD_quizartinib_analysis.html)

10. [Metformin analysis with custom regulons](https://raw.githack.com/SaccoPerfettoLab/SignalingProfiler_Benchmarking/main/scripts/7.Custom_PKN_analysis.html)
