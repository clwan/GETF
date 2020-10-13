# GETF
Introduction
------------

Boolean tensor has been broadly utilized in representing high
dimensional logical data collected on spatial, temporal and/or other
relational domains. Boolean Tensor Decomposition (BTD) factorizes a
binary tensor into the Boolean sum of multiple rank-1 tensors, which is
an NP-hard problem. Existing BTD methods have been limited by their high
computational cost, in applications to large scale or higher order
tensors. In this work, we presented a computationally efficient BTD
algorithm, namely *Geometric Expansion for all-order Tensor
Factorization* (GETF), that sequentially identifies the rank-1 basis
components for a tensor from a geometric perspective. We conducted
rigorous theoretical analysis on the validity as well as algorithemic
efficiency of GETF in decomposing all-order tensor. Experiments on both
synthetic and real-world data demonstrated that GETF has significantly
improved performance in reconstruction accuracy, extraction of latent
structures and it is an order of magnitude faster than other
state-of-the-art methods.

Currently, GETF has been accepetd by NeurIPS 2020. For more detail,
please refer and cite our paper [Wan, C., Chang, W., Zhao, T., Cao, S.,
& Zhang, C. (2020). Geometric All-Way Boolean Tensor Decomposition.
arXiv preprint arXiv:2007.15821.](https://arxiv.org/abs/2007.15821) And
if you run into some issue, please contact wan82@purdue.edu or
wancl0422@gmail.com for faster reply. I did not get notified by Github.

Process
-------

Currently, GETF takes the form of CP decomposition and all the functions
related to GETF are stored at GETF\_CP.R.

We provide a Tensor\_Simulate function for the ease of testing. It has
four inputs: Dims for the dimension of simulated tensor, pattern is the
latent pattern size, default 5, density is the density of latent
pattern, default is 0.2, noise is added noise ratio of data, default is
0.01

GETF\_CP is the function to conduct the sequential Boolean tensor
decomposition. The decomposition is based on the geometric property of
the tensor structure. We omitted the filtering process as we find it
does not contribute to the final outcome. GETF\_CP has 5 inputs: TENS is
the to-be-decomposed tensor data, Thres is the noise endurance level,
defualt is 0.6, meaning, current decomposition can take as much as 40%
noise. B\_num is the maximum number of decomposized patterns, default
20. COVER is limit of coverage ratio of decomposed patterns to cover the
original tensor, default 0.9. Exhausive is the indicator of whether use
exhuasive search for best patterns, default is F. Sometimes increase the
noise endurance level would generate better decomposition without
introducing much noise. As it could capture better geometric structure
of the tensor. In general GETF is very fast, Exhuasive search can be
conducted without took much additional time.

We also provide a function to test the reconstructed error of decomposed
pattern with original tensor.

    source("GETF_CP.R")

    TENS<-Tensor_Simulate(Dims = c(100,100,100),pattern = 5,density = 0.2,Noise = 0.01)

    Patterns<-GETF_CP(TENS = TENS,Thres = 0.6,B_num = 20,COVER = 0.9,Exhausive = F)

    Reconstruct_error(TENS,Patterns)

    ## [1] 0.01792

