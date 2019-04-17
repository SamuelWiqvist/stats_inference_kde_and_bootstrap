# load packages
using PyPlot
using RDatasets

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
            kde_est[i] = kde_est[i] + kernel((gird[i]-data[j])/bandwidth)
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


# rejection sampling algortihm
function rejectionsampling(p::Function, g::Function, sample::Function, M::Real, N::Int)


    samples = zeros(N)

    for i = 1:N

        generate_proposal = true

        while true
            proposal = sample()
            r = p(proposal)/(M*g(proposal))
            if rand() < r
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
M = 0.005

samples = rejectionsampling(p, g, sample, M, 10)

PyPlot.figure()
h = PyPlot.plt[:hist](samples, 100, density= true)
PyPlot.plot(kde_grid, kde_est)
