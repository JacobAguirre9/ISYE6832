module PlotsModule

using Plots

export plot_queue_length

function plot_queue_length(time_series::Vector{Float64}, queue_lengths::Vector{Float64})
    step = plot(time_series, queue_lengths, seriestype = :step, label = "Queue Length", xlabel = "Time", ylabel = "Number of Customers", title = "Queue Length Over Time")
    display(step)
end

end # module
