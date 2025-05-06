Xtail
=========

Genome-wide assessment of differential translations with ribosome profiling data

# VERSION

1.2.0

# Authors

Zhengtao Xiao

# DESCRIPTION

The xtail package is designed to identify genes exhibiting differential translation 
between two experimental conditions by simultaneously analyzing changes in RNA-seq 
and ribo-seq data. The estimation of changes in read counts are performed 
using the DESeq2 package.

# REQUIREMENTS

* R >= 3.2
* DESeq2
* Rcpp >= 0.10.1
* parallel 

# INSTALLATION

Details of the installation process can be found in the INSTALL file located in the root directory of the package.
Additionally, a Docker image of xtail is available; please refer to this [page](https://hub.docker.com/r/yanglab/xtail) for more information.

# CONTENTS

The R install packages are located in the "releases" tag.

# DOCUMENTATION

To use Xtail, please refer to the instructions in xtail.pdf in this subdirectory inst/doc,
or type: vignette("xtail") at the R-prompt.

# LICENSE

Xtail is licensed under the GPL version 3 or any later version


For more information please contact  xzt13[at]xjtu.edu.cn
