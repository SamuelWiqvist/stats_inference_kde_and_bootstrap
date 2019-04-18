# load packages
using PyPlot
using RDatasets
using Statistics

# load data
data = dataset("MASS", "geyser")
waiting = zeros(size(data,1))
for i = 1:size(data,1); waiting[i] = data[i,1]; end # load the waiting variable

# plot data
PyPlot.figure()
h = PyPlot.plt[:hist](waiting,10,density= true)


# function to compute the kde
function kde(data::Vector, grid::Vector, kernel::Function, bandwidth::Real)

    n = length(data)
    nbr_grid_points = length(grid)
    kde_est = zeros(length(grid))

    for i = 1:nbr_grid_points
        for j = 1:n
            kde_est[i] = kde_est[i] + kernel((grid[i]-data[j])/bandwidth)
        end
    end


    return 1/(n*bandwidth)*kde_est

end

# kernels
k_gaussian(u) = 1/(2*pi)*exp(-1/2*u^2)
k_uniform(u) = abs(u) <= 1 ? 1/2 : 0

# run kde est
kde_grid = collect(LinRange(40, 110, 50))
bandwidth = 3
kde_est = kde(waiting, kde_grid, k_uniform, bandwidth)

PyPlot.figure()
h = PyPlot.plt[:hist](waiting,20, density= true)
PyPlot.plot(kde_grid, kde_est)

# task a) Uniform kernel with bandwidth = 3 work ok, we also tried using the Gaussian kernel,
# but the unifrom kernal seemed to work better

# rejection sampling algortihm
function rejectionsampling(p::Function, g::Function, sample::Function, M::Real, N::Int)

    samples = zeros(N)

    for i = 1:N

        generate_proposal = true

        while generate_proposal
            proposal = sample()
            r = p(proposal)/(M*g(proposal))
            if rand() <= r
                samples[i] = proposal
                generate_proposal = false
            end
        end

    end

    return samples

end


# sample from kde est
p(x) = kde(waiting, [x], k_uniform, bandwidth)[1]
g(x) = 1
sample() = 40 + (110-40)*rand()
M = 0.05

samples = @time rejectionsampling(p, g, sample, M, 10000)


PyPlot.figure()
h = PyPlot.plt[:hist](samples, 100, density= true)
PyPlot.plot(kde_grid, kde_est)

# task b) We use unifrom proposal dist on [40,110], we get 10000 samples in 0.05 sec. Hard to
# find better prop dist since the kde is bimodal, could use a mixture of Gaussians but that
# is too complicated for this simple problem


# compute quantile bootstrap confidence band
B = 100 # bootstrap sampels
nbr_samples = length(waiting) # nbr of samples in each bootstrap sample

# pre-allocated vectors
samples = zeros(nbr_samples, B)
kde_ests = zeros(length(kde_grid), B)
kde_quantiles = zeros(2, length(kde_grid))

# generate bootstrap samples
for i = 1:B
    samples[:,i] = rejectionsampling(p, g, sample, M, nbr_samples) # generate bootstrap sample
    kde_ests[:,i] = kde(samples[:,i], kde_grid, k_uniform, bandwidth) # compute kde for bootstrap sample
end

# compute bootstrap quantiles for each point where we evalute the kde
for i = 1:length(kde_grid)
    kde_quantiles[:,i] = quantile(kde_ests[i,:], [0.025, 0.975])
end

# plot quantile bootstrap confidence band
PyPlot.figure()
PyPlot.plot(kde_grid, kde_est, "r")
PyPlot.plot(kde_grid, kde_ests[:,4], "k")
PyPlot.fill_between(kde_grid, kde_quantiles[1,:], kde_quantiles[2,:],alpha=0.5)

# task c) the quantile bootstrap confidence band seem to be ok
