# imagePartitionining

#### Authors: Mark Beers, Klim Rayskiy, Brian Xy, Tianhe Wang

In this repository we implement a basic image segmentation algorithm based on the Geman and Geman loss function for CS168/CS268. This function makes a determination regarding whether or not an edge exists between two pixels by comparing the intensity values of those pixels. For a given image "im", the Geman and Geman cost function attempts to generate a smoother corresponding image "f", a matrix of estimated horizontal edges "h", and a matrix of estimated vertical edges "v". In practice,there is an edge, or there is not. However, during the optimization, this is relaxed so that the presence of an edge can be represented on a continuous scale from zero to one. As the optimization progresses, these edge values approach their boundaries and eventually we classify them as horizontal or vertical edges by rounding to the nearest integer. 

An interior point algorithm is used to optimize the cost function. Gradient descent with backwards line search is used at each step in this algorithm. We write code to ensure that updates to h and v never cause any value in h or v to leave the [0,1] interval. For this project, the noisy rectangular images can be found in the images2 folder. Optimization code can be found in "objective3.jl". Functions for computing F-values, as measures of performance, can be found in "numerical_statistics.jl". The Geman and Geman cost function contains two parameters whose values have to be selected by hand for the specific image dataset. "Sensitiviy Analysis.ipynb" contains some code for searching over a grid to find good values for these parameters. Below we have an implementation of how to use the code to generate an outline of an image. This code was taken from "Sensitiviy Analysis.ipynb" and should produce the attached image. 

```julia
include("objective3.jl")
include("numerical_statistics.jl")
using Images, FileIO, Statistics, LinearAlgebra

im_path = "mjolsness.png"
img = load(im_path);
imgg = Gray.(img);
im = convert(Array{Float64}, imgg);
colorview(Gray, im)

imheight, imwidth = size(im) 
f = .5 .+ randn(imheight, imwidth)/10;
h = .5 .+ randn(imheight-1, imwidth)/10;
v = .5 .+ randn(imheight, imwidth-1)/10;

niters = 7
ρs = [1.5625*.2^k for k in 4:niters]
αws = repeat([.005], 4) # .005 seems to work with p, .002 seems to work well with q. .05 maximum
λs = repeat([1.1], 4)
f,h,v = interior_point(gg, ∇gg, q, ∇q, f, h, v, λs, ρs, αws, im)

edge_image = 1 .- max.(v[1:imheight-1, 1:imwidth-1],h[1:imheight-1, 1:imwidth-1])
edge_image_binary = edge_image_ant .> .5
mosaic(colorview(Gray, im), colorview(Gray, edge_image_binary), nrow= 1)
```


![alt text](https://github.com/krayskiy/imagePartitionining/blob/main/anteater_edge.png)
