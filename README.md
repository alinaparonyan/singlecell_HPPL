# Parallelized Single Cell Processing

__Navigation:__
* _/scripts_ folder contains our implementation of QC calculation & normalization.
* _speedup_benchmark.ipynb_ notebook contains comparison of three approaches: plain Python (Numpy), Numba and Cuda for two tasks. These tasks are mithocondrial QC calculation and normalization.

<br>

__Workflow:__
* The data is taken from NeurIPS 2021 Bioinformatics Competition and can be loaded from [this GoogleDrive link](https://drive.google.com/file/d/1q7gn2VLmgJRH50gOlKPlj8lwJiLDl2Eb/view?usp=sharing).
* To reproduce comparison notebook, do the following:
    ```python3
    conda create env --name sc_speedup python=3.9
    conda activate sc_speedup
    pip install -r requirements.txt
    ``` 
