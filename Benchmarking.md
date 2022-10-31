# NeuroAnalyzer benchmarking

Set JULIA_NUM_THREADS and --threads according to your system:
```sh
julia -O3 -g0 --threads 24 --cpu-target=native benchmarks/na_benchmarking.jl > benchmarks/results.txt
```

## Results

Here are some results:
- [Lenovo T14](benchmarks/t14.md)
- [Threadripper 3960X, NVIDIA Quadro P2200, CUDA enabled](benchmarks/3960x_cuda.md) (48 threads)
- [Threadripper 3960X, CUDA disabled](benchmarks/3960x_nocuda.md) (24 threads)