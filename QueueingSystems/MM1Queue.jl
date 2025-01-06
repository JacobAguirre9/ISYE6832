module MM1QueueModule

using ..SimulationEngine.SchedulerModule
using ..SimulationEngine.EventModule
using ..Utilities.RandomGenerators
using ..Utilities.StatisticsModule

export MM1Queue

struct MM1Queue
    λ::Float64
    μ::Float64
    scheduler::SchedulerModule.Scheduler
    stats::StatisticsModule.Statistics
    queue::Vector{Float64} # Stores arrival times
    server_busy::Bool
    customer_id::Int
end

function MM1Queue(λ::Float64, μ::Float64, scheduler::SchedulerModule.Scheduler, stats::StatisticsModule.Statistics)
    return MM1Queue(λ, μ, scheduler, stats, Float64[], false, 0)
end

function arrival_action!(queue::MM1Queue, entity::EntityModule.Entity)
    queue.stats.customers += 1
    if !queue.server_busy
        queue.server_busy = true
        service_time = exponential(queue.μ)
        departure_event = EventModule.Event(
            time = queue.scheduler.current_time + service_time,
            priority = 1,
            action = () -> departure_action!(queue),
            description = "Departure"
        )
        schedule_event!(queue.scheduler, departure_event)
    else
        push!(queue.queue, queue.scheduler.current_time)
    end

    # Schedule next arrival
    inter_arrival = exponential(queue.λ)
    next_arrival_time = queue.scheduler.current_time + inter_arrival
    new_entity = EntityModule.Entity(queue.customer_id += 1, next_arrival_time)
    arrival_event = EventModule.Event(
        time = next_arrival_time,
        priority = 0,
        action = () -> arrival_action!(queue, new_entity),
        description = "Arrival"
    )
    schedule_event!(queue.scheduler, arrival_event)
end

function departure_action!(queue::MM1Queue)
    if !isempty(queue.queue)
        arrival_time = popfirst!(queue.queue)
        wait_time = queue.scheduler.current_time - arrival_time
        queue.stats.total_wait_time += wait_time
        service_time = exponential(queue.μ)
        departure_event = EventModule.Event(
            time = queue.scheduler.current_time + service_time,
            priority = 1,
            action = () -> departure_action!(queue),
            description = "Departure"
        )
        schedule_event!(queue.scheduler, departure_event)
    else
        queue.server_busy = false
    end
end

function initialize!(queue::MM1Queue)
    # Schedule first arrival
    first_arrival_time = exponential(queue.λ) + queue.scheduler.current_time
    first_entity = EntityModule.Entity(1, first_arrival_time)
    arrival_event = EventModule.Event(
        time = first_arrival_time,
        priority = 0,
        action = () -> arrival_action!(queue, first_entity),
        description = "Arrival"
    )
    schedule_event!(queue.scheduler, arrival_event)
end

end # module
