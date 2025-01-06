# src/SimulationEngine/Event.jl

module EventModule

export Event

struct Event
    time::Float64
    priority::Int
    action::Function
    description::String
end

# To make Event sortable by time and priority
import Base: <, ==
function Base.<(e1::Event, e2::Event)
    return (e1.time < e2.time) || (e1.time == e2.time && e1.priority < e2.priority)
end

function Base.==(e1::Event, e2::Event)
    return e1.time == e2.time && e1.priority == e2.priority
end

end # module
