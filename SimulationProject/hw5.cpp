#include <iostream>
#include <queue>
#include <random>
#include <vector>
#include <functional>
#include <memory>
#include <cassert>

// Define simulation parameters.
const double LAMBDA = 10.0; 
const double MU1 = 15.0;
const double MU2 = 20.0;
const double ROUTE_PROB = 0.1;
const double SIM_TIME = 100.0;
const int NUM_REPLICATIONS = 1000;

enum class EventType { ARRIVAL, DEPARTURE1, DEPARTURE2 };

struct Customer {
    double arrivalTime;
    Customer(double t) : arrivalTime(t) {}
};

struct Event {
    double time; 
    EventType type;
    std::shared_ptr<Customer> cust;
    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

std::pair<double, double> runSimulation(std::mt19937 &rng) {
    // Exponential distributions for interarrival and service times.
    std::exponential_distribution<double> expArrival(LAMBDA);
    std::exponential_distribution<double> expService1(MU1);
    std::exponential_distribution<double> expService2(MU2);
    std::uniform_real_distribution<double> uniformDist(0.0, 1.0);
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> eventList;
    std::queue<std::shared_ptr<Customer>> queue1;
    std::queue<std::shared_ptr<Customer>> queue2;
    bool busy1 = false;
    bool busy2 = false;
    double timeLast = 0.0;
    double areaNumInSystem = 0.0;
    int numInSystem = 0;
    
  // Sojourn time statistics.
    double totalSojournTime = 0.0;
    int departedCustomers = 0;

    // Schedule the first external arrival at time 0.
    {
        auto cust = std::make_shared<Customer>(0.0);
        Event firstArrival {0.0, EventType::ARRIVAL, cust};
        eventList.push(firstArrival);
    }

    // Also schedule the next external arrival.
    double nextArrivalTime = expArrival(rng);
    {
        auto cust = std::make_shared<Customer>(nextArrivalTime);
        Event arrivalEvent {nextArrivalTime, EventType::ARRIVAL, cust};
        eventList.push(arrivalEvent);
    }
    numInSystem++;  // first arrival counted in system

    // Main simulation loop.
    while (!eventList.empty()) {
        Event ev = eventList.top();
        eventList.pop();

        if (ev.time > SIM_TIME) {
            areaNumInSystem += numInSystem * (SIM_TIME - timeLast);
            break;
        }

        areaNumInSystem += numInSystem * (ev.time - timeLast);
        timeLast = ev.time;
        if (ev.type == EventType::ARRIVAL) {
            queue1.push(ev.cust);

            if (!busy1) {
                busy1 = true;
                auto cust = queue1.front();
                queue1.pop();
                double serviceTime = expService1(rng);
                Event dep1 {ev.time + serviceTime, EventType::DEPARTURE1, cust};
                eventList.push(dep1);
            }

            if (std::abs(ev.cust->arrivalTime - ev.time) < 1e-8) {
                double newArrivalTime = ev.time + expArrival(rng);
                if (newArrivalTime <= SIM_TIME) {
                    auto newCust = std::make_shared<Customer>(newArrivalTime);
                    Event newArr {newArrivalTime, EventType::ARRIVAL, newCust};
                    eventList.push(newArr);
                    numInSystem++; // new arrival added
                }
            }
        } else if (ev.type == EventType::DEPARTURE1) {
            if (!queue1.empty()) {
                auto cust = queue1.front();
                queue1.pop();
                double serviceTime = expService1(rng);
                Event dep1 {ev.time + serviceTime, EventType::DEPARTURE1, cust};
                eventList.push(dep1);
            } else {
                busy1 = false;
            }
            queue2.push(ev.cust);
            if (!busy2) {
                busy2 = true;
                auto cust = queue2.front();
                queue2.pop();
                double serviceTime = expService2(rng);
                Event dep2 {ev.time + serviceTime, EventType::DEPARTURE2, cust};
                eventList.push(dep2);
            }
        } else if (ev.type == EventType::DEPARTURE2) {
            double r = uniformDist(rng);
            if (r <= 0.9) {
                double sojournTime = ev.time - ev.cust->arrivalTime;
                totalSojournTime += sojournTime;
                departedCustomers++;
                numInSystem--;  // customer leaves system
            } else {
                queue1.push(ev.cust);
                if (!busy1) {
                    busy1 = true;
                    auto cust = queue1.front();
                    queue1.pop();
                    double serviceTime = expService1(rng);
                    Event dep1 {ev.time + serviceTime, EventType::DEPARTURE1, cust};
                    eventList.push(dep1);
                }
            }
            if (!queue2.empty()) {
                auto cust = queue2.front();
                queue2.pop();
                double serviceTime = expService2(rng);
                Event dep2 {ev.time + serviceTime, EventType::DEPARTURE2, cust};
                eventList.push(dep2);
            } else {
                busy2 = false;
            }
        }
    } // end simulation loop

    if (timeLast < SIM_TIME) {
        areaNumInSystem += numInSystem * (SIM_TIME - timeLast);
    }

    double avgNumInSystem = areaNumInSystem / SIM_TIME;
    double avgSojournTime = (departedCustomers > 0) ? totalSojournTime / departedCustomers : 0.0;

    return {avgNumInSystem, avgSojournTime};
}

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());

    double sumN = 0.0;
    double sumT = 0.0;
    for (int rep = 0; rep < NUM_REPLICATIONS; rep++) {
        auto result = runSimulation(rng);
        sumN += result.first;
        sumT += result.second;
    }

    double overallAvgN = sumN / NUM_REPLICATIONS;
    double overallAvgT = sumT / NUM_REPLICATIONS;

    std::cout << "Simulation results over " << NUM_REPLICATIONS << " replications (100 hours each):" << std::endl;
    std::cout << "Average number of customers in the system (n): " << overallAvgN << std::endl;
    std::cout << "Average time in the system per customer (t): " << overallAvgT << " hours" << std::endl;

    return 0;
}
