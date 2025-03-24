/*
QueueingNetworkSimulation_OutputAnalysis.cpp
*/

#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <limits>
#include <cmath>
#include <numeric>
#include <cassert>
#include <algorithm>

using namespace std;

const double INF = numeric_limits<double>::infinity();

// ------------------------------
// Structures for simulation results
// ------------------------------
struct SimulationResult {
    double totalTime;           // Total simulation time
    double integratedCustomers; // ∫ X(t) dt over simulation period
    double integratedIndicator; // ∫ I{X(t) ≤ 2} dt over simulation period
    int eventCount;             // Number of events processed
};

// Structure to hold the output analysis result from one complete run
struct MethodResult {
    double L_est;       // Point estimate for L
    double L_lower;     // Lower bound of 95% CI for L
    double L_upper;     // Upper bound of 95% CI for L
    double p_est;       // Point estimate for p
    double p_lower;     // Lower bound of 95% CI for p
    double p_upper;     // Upper bound of 95% CI for p
};

// ------------------------------
// Exponential RV generator.
double expRand(double rate, mt19937 &rng) {
    exponential_distribution<double> dist(rate);
    return dist(rng);
}

// ------------------------------
// Simulation procedures (inner methods)
// ------------------------------

// (A) Event-driven simulation until time T, starting from state (0,0).
SimulationResult simulateUntilTime(double T, double lambda, double mu1, double mu2, double q, mt19937 &rng) {
    int x1 = 0, x2 = 0;
    double t = 0.0;
    SimulationResult res {0.0, 0.0, 0.0, 0};
    double t_next_A = t + expRand(lambda, rng);
    double t_next_S1 = INF, t_next_S2 = INF;
    if (x1 > 0) t_next_S1 = t + expRand(mu1, rng);
    if (x2 > 0) t_next_S2 = t + expRand(mu2, rng);
    while (t < T) {
        double t_next = min({t_next_A, t_next_S1, t_next_S2});
        if (t + (t_next - t) > T) {
            double dt = T - t;
            int tot = x1 + x2;
            res.integratedCustomers += tot * dt;
            res.integratedIndicator += (tot <= 2 ? dt : 0.0);
            t = T;
            break;
        }
        double dt = t_next - t;
        int tot = x1 + x2;
        res.integratedCustomers += tot * dt;
        res.integratedIndicator += (tot <= 2 ? dt : 0.0);
        t = t_next;
        res.eventCount++;
        if (fabs(t - t_next_A) < 1e-12) {
            x1++;
            t_next_A = t + expRand(lambda, rng);
            if (x1 == 1) t_next_S1 = t + expRand(mu1, rng);
        } else if (fabs(t - t_next_S1) < 1e-12) {
            x1--;
            x2++;
            if (x1 > 0) t_next_S1 = t + expRand(mu1, rng);
            else t_next_S1 = INF;
            if (x2 == 1) t_next_S2 = t + expRand(mu2, rng);
        } else if (fabs(t - t_next_S2) < 1e-12) {
            x2--;
            double u = uniform_real_distribution<double>(0.0, 1.0)(rng);
            if (u > q) {
                x1++;
                if (x1 == 1) t_next_S1 = t + expRand(mu1, rng);
            }
            if (x2 > 0) t_next_S2 = t + expRand(mu2, rng);
            else t_next_S2 = INF;
        }
    }
    res.totalTime = t;
    return res;
}

// (B) Regenerative cycle simulation.
// A regenerative cycle is defined as the time interval between successive visits
// to state (0,0) after having left (0,0).
SimulationResult simulateRegenerativeCycle(double lambda, double mu1, double mu2, double q, mt19937 &rng) {
    int x1 = 0, x2 = 0;
    double t = 0.0;
    SimulationResult cycle {0.0, 0.0, 0.0, 0};
    double t_next_A = t + expRand(lambda, rng);
    double t_next_S1 = INF, t_next_S2 = INF;
    bool cycleStarted = false;
    while (true) {
        double t_next = min({t_next_A, t_next_S1, t_next_S2});
        double dt = t_next - t;
        int tot = x1 + x2;
        cycle.integratedCustomers += tot * dt;
        cycle.integratedIndicator += (tot <= 2 ? dt : 0.0);
        t = t_next;
        cycle.eventCount++;
        if (fabs(t - t_next_A) < 1e-12) {
            x1++;
            if (!cycleStarted && (x1 + x2) > 0)
                cycleStarted = true;
            t_next_A = t + expRand(lambda, rng);
            if (x1 == 1) t_next_S1 = t + expRand(mu1, rng);
        } else if (fabs(t - t_next_S1) < 1e-12) {
            x1--;
            x2++;
            if (x1 > 0) t_next_S1 = t + expRand(mu1, rng);
            else t_next_S1 = INF;
            if (x2 == 1) t_next_S2 = t + expRand(mu2, rng);
        } else if (fabs(t - t_next_S2) < 1e-12) {
            x2--;
            double u = uniform_real_distribution<double>(0.0, 1.0)(rng);
            if (u > q) {
                x1++;
                if (x1 == 1) t_next_S1 = t + expRand(mu1, rng);
            }
            if (x2 > 0) t_next_S2 = t + expRand(mu2, rng);
            else t_next_S2 = INF;
        }
        if (cycleStarted && x1 == 0 && x2 == 0) break;
    }
    cycle.totalTime = t;
    return cycle;
}

// ------------------------------
// Utility functions for computing sample mean and variance.
double sampleMean(const vector<double>& data) {
    return accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double sampleVariance(const vector<double>& data, double mean) {
    double sumSq = 0.0;
    for (double x : data) sumSq += (x - mean) * (x - mean);
    return sumSq / (data.size() - 1);
}

// t-critical value for 95% CI (using normal approximation for large df)
double tCritical(int df) {
    if (df >= 30) return 2.0; // approximate
    // For small df, one might use a table; here we use a simple approximation.
    return 2.045;
}

// ------------------------------
// Methods to perform a complete output analysis procedure run.
// Each function returns a MethodResult structure with point estimates and 95% CIs.

// Method I: Multiple Independent Replications.
MethodResult runMultipleReplications(double T, int R, double lambda, double mu1, double mu2, double q, mt19937 &rng) {
    vector<double> Lvals, pvals;
    for (int r = 0; r < R; r++) {
        SimulationResult res = simulateUntilTime(T, lambda, mu1, mu2, q, rng);
        double L_rep = res.integratedCustomers / res.totalTime;
        double p_rep = res.integratedIndicator / res.totalTime;
        Lvals.push_back(L_rep);
        pvals.push_back(p_rep);
    }
    double meanL = sampleMean(Lvals);
    double varL = sampleVariance(Lvals, meanL);
    double meanP = sampleMean(pvals);
    double varP = sampleVariance(pvals, meanP);
    int df = R - 1;
    double tCrit = tCritical(df);
    double seL = sqrt(varL / R);
    double seP = sqrt(varP / R);
    MethodResult result;
    result.L_est = meanL;
    result.L_lower = meanL - tCrit * seL;
    result.L_upper = meanL + tCrit * seL;
    result.p_est = meanP;
    result.p_lower = meanP - tCrit * seP;
    result.p_upper = meanP + tCrit * seP;
    return result;
}

// Method II: Batch Means.
MethodResult runBatchMeans(double T_total, double T_init, int B, double lambda, double mu1, double mu2, double q, mt19937 &rng) {
    double batchLength = (T_total - T_init) / B;
    vector<double> Lvals(B, 0.0), Pvals(B, 0.0), batchTimes(B, 0.0);
    double t = 0.0;
    int x1 = 0, x2 = 0;
    double t_next_A = t + expRand(lambda, rng);
    double t_next_S1 = INF, t_next_S2 = INF;
    while (t < T_total) {
        double t_next = min({t_next_A, t_next_S1, t_next_S2});
        if (t_next > T_total) t_next = T_total;
        double dt = t_next - t;
        int tot = x1 + x2;
        if (t >= T_init) {
            int batchIndex = min(static_cast<int>((t - T_init) / batchLength), B - 1);
            Lvals[batchIndex] += tot * dt;
            Pvals[batchIndex] += (tot <= 2 ? dt : 0.0);
            batchTimes[batchIndex] += dt;
        }
        t = t_next;
        if (fabs(t - t_next_A) < 1e-12) {
            x1++;
            t_next_A = t + expRand(lambda, rng);
            if (x1 == 1) t_next_S1 = t + expRand(mu1, rng);
        } else if (fabs(t - t_next_S1) < 1e-12) {
            x1--;
            x2++;
            if (x1 > 0) t_next_S1 = t + expRand(mu1, rng);
            else t_next_S1 = INF;
            if (x2 == 1) t_next_S2 = t + expRand(mu2, rng);
        } else if (fabs(t - t_next_S2) < 1e-12) {
            x2--;
            double u = uniform_real_distribution<double>(0.0, 1.0)(rng);
            if (u > q) {
                x1++;
                if (x1 == 1) t_next_S1 = t + expRand(mu1, rng);
            }
            if (x2 > 0) t_next_S2 = t + expRand(mu2, rng);
            else t_next_S2 = INF;
        }
    }
    vector<double> L_batch(B), p_batch(B);
    for (int i = 0; i < B; i++) {
        if (batchTimes[i] > 0) {
            L_batch[i] = Lvals[i] / batchTimes[i];
            p_batch[i] = Pvals[i] / batchTimes[i];
        }
    }
    double meanL = sampleMean(L_batch);
    double varL = sampleVariance(L_batch, meanL);
    double meanP = sampleMean(p_batch);
    double varP = sampleVariance(p_batch, meanP);
    int df = B - 1;
    double tCrit = tCritical(df);
    double seL = sqrt(varL / B);
    double seP = sqrt(varP / B);
    MethodResult result;
    result.L_est = meanL;
    result.L_lower = meanL - tCrit * seL;
    result.L_upper = meanL + tCrit * seL;
    result.p_est = meanP;
    result.p_lower = meanP - tCrit * seP;
    result.p_upper = meanP + tCrit * seP;
    return result;
}

// Method III: Regenerative Method.
MethodResult runRegenerativeMethod(int N, double lambda, double mu1, double mu2, double q, mt19937 &rng) {
    vector<double> cycleT, cycleR, cycleU;
    int totalEvents = 0;
    for (int i = 0; i < N; i++) {
        SimulationResult cycle = simulateRegenerativeCycle(lambda, mu1, mu2, q, rng);
        cycleT.push_back(cycle.totalTime);
        cycleR.push_back(cycle.integratedCustomers);
        cycleU.push_back(cycle.integratedIndicator);
        totalEvents += cycle.eventCount;
    }
    double sumT = accumulate(cycleT.begin(), cycleT.end(), 0.0);
    double sumR = accumulate(cycleR.begin(), cycleR.end(), 0.0);
    double sumU = accumulate(cycleU.begin(), cycleU.end(), 0.0);
    double L_est = sumR / sumT;
    double p_est = sumU / sumT;
    vector<double> d(N), d_p(N);
    for (int i = 0; i < N; i++) {
        d[i] = cycleR[i] - L_est * cycleT[i];
        d_p[i] = cycleU[i] - p_est * cycleT[i];
    }
    double s2 = sampleVariance(d, 0.0);
    double s2_p = sampleVariance(d_p, 0.0);
    double meanT = sumT / N;
    double seL = sqrt(s2 / N) / meanT;
    double seP = sqrt(s2_p / N) / meanT;
    int df = N - 1;
    double tCrit = tCritical(df);
    MethodResult result;
    result.L_est = L_est;
    result.L_lower = L_est - tCrit * seL;
    result.L_upper = L_est + tCrit * seL;
    result.p_est = p_est;
    result.p_lower = p_est - tCrit * seP;
    result.p_upper = p_est + tCrit * seP;
    return result;
}

// ------------------------------
// Functions for outer replications of output analysis procedures.
// For each method, we run outerRep outer replications of the method and compute the sample mean, variance, and coverage probability for L and p.
struct OuterStats {
    double meanL;
    double varL;
    double meanp;
    double varp;
    double coverageL; // fraction of outer replications where true L lies in CI
    double coveragep; // fraction for p
};

OuterStats outerReplication_MR(double T, int R, int outerRep, double lambda, double mu1, double mu2, double q,
                                double L_true, double p_true, mt19937 &rng) {
    vector<double> L_outer, p_outer;
    int coverL = 0, coverp = 0;
    for (int i = 0; i < outerRep; i++) {
        MethodResult res = runMultipleReplications(T, R, lambda, mu1, mu2, q, rng);
        L_outer.push_back(res.L_est);
        p_outer.push_back(res.p_est);
        if (L_true >= res.L_lower && L_true <= res.L_upper) coverL++;
        if (p_true >= res.p_lower && p_true <= res.p_upper) coverp++;
    }
    OuterStats stats;
    stats.meanL = sampleMean(L_outer);
    stats.varL = sampleVariance(L_outer, stats.meanL);
    stats.meanp = sampleMean(p_outer);
    stats.varp = sampleVariance(p_outer, stats.meanp);
    stats.coverageL = static_cast<double>(coverL) / outerRep;
    stats.coveragep = static_cast<double>(coverp) / outerRep;
    return stats;
}

OuterStats outerReplication_BM(double T_total, double T_init, int B, int outerRep, double lambda, double mu1, double mu2, double q,
                                double L_true, double p_true, mt19937 &rng) {
    vector<double> L_outer, p_outer;
    int coverL = 0, coverp = 0;
    for (int i = 0; i < outerRep; i++) {
        MethodResult res = runBatchMeans(T_total, T_init, B, lambda, mu1, mu2, q, rng);
        L_outer.push_back(res.L_est);
        p_outer.push_back(res.p_est);
        if (L_true >= res.L_lower && L_true <= res.L_upper) coverL++;
        if (p_true >= res.p_lower && p_true <= res.p_upper) coverp++;
    }
    OuterStats stats;
    stats.meanL = sampleMean(L_outer);
    stats.varL = sampleVariance(L_outer, stats.meanL);
    stats.meanp = sampleMean(p_outer);
    stats.varp = sampleVariance(p_outer, stats.meanp);
    stats.coverageL = static_cast<double>(coverL) / outerRep;
    stats.coveragep = static_cast<double>(coverp) / outerRep;
    return stats;
}

OuterStats outerReplication_Reg(int N, int outerRep, double lambda, double mu1, double mu2, double q,
                                 double L_true, double p_true, mt19937 &rng) {
    vector<double> L_outer, p_outer;
    int coverL = 0, coverp = 0;
    for (int i = 0; i < outerRep; i++) {
        MethodResult res = runRegenerativeMethod(N, lambda, mu1, mu2, q, rng);
        L_outer.push_back(res.L_est);
        p_outer.push_back(res.p_est);
        if (L_true >= res.L_lower && L_true <= res.L_upper) coverL++;
        if (p_true >= res.p_lower && p_true <= res.p_upper) coverp++;
    }
    OuterStats stats;
    stats.meanL = sampleMean(L_outer);
    stats.varL = sampleVariance(L_outer, stats.meanL);
    stats.meanp = sampleMean(p_outer);
    stats.varp = sampleVariance(p_outer, stats.meanp);
    stats.coverageL = static_cast<double>(coverL) / outerRep;
    stats.coveragep = static_cast<double>(coverp) / outerRep;
    return stats;
}

// ------------------------------
// True values for L and p based on open Jackson network formulas.
void computeTrueValues(double rho1, double rho2, double &L_true, double &p_true) {
    L_true = rho1 / (1 - rho1) + rho2 / (1 - rho2);
    double pi00 = (1 - rho1) * (1 - rho2);
    double pi10 = (1 - rho1) * (1 - rho2) * rho1;
    double pi01 = (1 - rho1) * (1 - rho2) * rho2;
    double pi11 = (1 - rho1) * (1 - rho2) * rho1 * rho2;
    double pi20 = (1 - rho1) * (1 - rho2) * rho1 * rho1;
    double pi02 = (1 - rho1) * (1 - rho2) * rho2 * rho2;
    p_true = pi00 + pi10 + pi01 + pi11 + pi20 + pi02;
}

// ------------------------------
// Main: Run outer replications for each method and parameter set.
int main() {
    mt19937 rng(static_cast<unsigned>(chrono::high_resolution_clock::now().time_since_epoch().count()));
    double mu1 = 1.0, mu2 = 1.0, q = 0.5;
    vector<double> targetRho = {0.1, 0.5, 0.9};
    int outerRep = 30;  // Number of outer replications for output analysis
    cout << "Running outer replications for output analysis...\n\n";
    for (double rho : targetRho) {
        double lambda = q * mu1 * rho; // λ = q * μ * ρ, with μ = 1.
        double rho1 = lambda / (q * mu1);
        double rho2 = lambda / (q * mu2);
        double L_true, p_true;
        computeTrueValues(rho1, rho2, L_true, p_true);
        cout << "========================================\n";
        cout << "Parameters: λ = " << lambda << ", μ₁ = " << mu1
             << ", μ₂ = " << mu2 << ", q = " << q << "\n";
        cout << "Traffic intensities: ρ₁ = " << rho1 << ", ρ₂ = " << rho2 << "\n";
        cout << "True L = " << L_true << ", True p = " << p_true << "\n\n";
        
        // Settings for inner simulation procedures.
        int R_inner = 30;     // for Multiple Replications method inner run
        double T = 1e5;       // simulation horizon for MR method
        double T_total = 1e6; // total simulation time for Batch Means
        double T_init = 1e4;  // truncation period for Batch Means
        int B = 30;           // number of batches for Batch Means
        int N = 100;          // number of regenerative cycles for Regenerative Method
        
        // Outer replications for each method.
        OuterStats statsMR = outerReplication_MR(T, R_inner, outerRep, lambda, mu1, mu2, q, L_true, p_true, rng);
        OuterStats statsBM = outerReplication_BM(T_total, T_init, B, outerRep, lambda, mu1, mu2, q, L_true, p_true, rng);
        OuterStats statsReg = outerReplication_Reg(N, outerRep, lambda, mu1, mu2, q, L_true, p_true, rng);
        
        cout << "Method I: Multiple Independent Replications:\n";
        cout << "  Outer Replications: " << outerRep << "\n";
        cout << "  Mean L estimate = " << statsMR.meanL << ", Variance = " << statsMR.varL
             << ", Coverage = " << statsMR.coverageL * 100 << "%\n";
        cout << "  Mean p estimate = " << statsMR.meanp << ", Variance = " << statsMR.varp
             << ", Coverage = " << statsMR.coveragep * 100 << "%\n\n";
        
        cout << "Method II: Batch Means:\n";
        cout << "  Outer Replications: " << outerRep << "\n";
        cout << "  Mean L estimate = " << statsBM.meanL << ", Variance = " << statsBM.varL
             << ", Coverage = " << statsBM.coverageL * 100 << "%\n";
        cout << "  Mean p estimate = " << statsBM.meanp << ", Variance = " << statsBM.varp
             << ", Coverage = " << statsBM.coveragep * 100 << "%\n\n";
        
        cout << "Method III: Regenerative Method:\n";
        cout << "  Outer Replications: " << outerRep << "\n";
        cout << "  Mean L estimate = " << statsReg.meanL << ", Variance = " << statsReg.varL
             << ", Coverage = " << statsReg.coverageL * 100 << "%\n";
        cout << "  Mean p estimate = " << statsReg.meanp << ", Variance = " << statsReg.varp
             << ", Coverage = " << statsReg.coveragep * 100 << "%\n\n";
        
        cout << "-----------------------------------------------------\n";
    }
    
    
    return 0;
}

