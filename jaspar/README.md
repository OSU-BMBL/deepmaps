# JASPAR TFBS tracks

DeepMAPS used TFBS prediction data from [JASPAR database](https://jaspar.genereg.net/genome-tracks/), and has converted to serilized R obejct for quick loading. We filtered out TFBS using the [cutoff of 500 (p-value < 10^-5)](https://github.com/wassermanlab/JASPAR-UCSC-tracks) The JASPAR TFBS data was used for calculating regulons ([_Calregulon function_](https://github.com/OSU-BMBL/deepmaps/blob/master/scRNA_scATAC1.r#L817)).

To read `.qsave` file, you need to install `qs` R package, which is faster to load compared to reading `rds` file. For example:

```
tfbs_df <- qs::qread("jaspar_hg38_500.qsave")
```

## Download

- hg38 (221 MB): https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/jaspar_hg38_500.qsave
- mm10 (10 MB): https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/jaspar_mm10_500.qsave
