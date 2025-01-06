module SchedulerModule

using Base.Dates
using DataStructures
using .EventModule

export Scheduler

struct Scheduler
    event_queue::PriorityQueue{Int, EventModule.Event}
    current_time::Float64
end

function Scheduler()
    pq = PriorityQueue{Int, EventModule.Event}(lt = (x, y) -> x < y)
    return Scheduler(pq, 0.0)
end

function schedule_event!(scheduler::Scheduler, event::EventModule.Event)
    # PriorityQueue in Julia pops the largest key first, so we use negative time for min-heap behavior
    enqueue!(scheduler.event_queue, -event.time, event)
end

function run!(scheduler::Scheduler, until::Float64)
    while !isempty(scheduler.event_queue)
        _, event = dequeue_pair!(scheduler.event_queue)
        event_time = event.time
        if event_time > until
            break
        end
        scheduler.current_time = event_time
        event.action()
    end
    scheduler.current_time = until
end

end # module
