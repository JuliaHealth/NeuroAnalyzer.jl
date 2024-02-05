# NeuroAnalyzer benchmarking

Set `JULIA_NUM_THREADS` and/or `--threads` according to your system:
```sh
julia -O3 -g0 --threads 24 --cpu-target=native benchmarks/na_benchmarking.jl > benchmarks/results.txt
```

## Results

Here are some results:

- [Lenovo T14](benchmarks/t14.txt) (AMD Ryzen 5 PRO 4650U, 12 threads, CUDA disabled)
- [AMD Threadripper workstation](benchmarks/3960x_nocuda.txt) (AMD Threadripper 3960X, 24 threads, CUDA disabled)
- [Intel i5 workstation](benchmarks/i5-4570.txt) (i5-4570, 2 threads)