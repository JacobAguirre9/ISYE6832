
module GSMPModelModule

using ..SimulationEngine.SchedulerModule
using ..Utilities.StatisticsModule

export GSMPModel

struct GSMPModel
    states::Vector{Any}
    current_state::Any
    transition_matrix::Dict{Any, Dict{Any, Float64}} # Example structure
    scheduler::SchedulerModule.Scheduler
    stats::StatisticsModule.Statistics
end

function GSMPModel(states::Vector{Any}, transition_matrix::Dict{Any, Dict{Any, Float64}}, scheduler::SchedulerModule.Scheduler, stats::StatisticsModule.Statistics)
    return GSMPModel(states, states[1], transition_matrix, scheduler, stats)
end

function transition!(model::GSMPModel)
  # add later, this will be a lil difficult
end

function run!(model::GSMPModel, until::Float64)
    run!(model.scheduler, until)
end

end # module
